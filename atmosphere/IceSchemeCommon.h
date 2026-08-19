#pragma once

#include "MixtureAtm.h"
#include "SaturationH2O.h"
#include "cAtmosphereModel.h"

// ============================================================================
// Shared machinery for the ice / precipitation microphysics schemes
// (Zero/One/Two/ThreeCatIceScheme). The four schemes share a large amount of
// physics-independent boilerplate — boundary extrapolation, topography fill,
// saturation, the finiteness flux cap — which was previously duplicated (and
// independently mis-fixed) in each file. It is factored out here so a fix lives
// in ONE place. The per-scheme MICROPHYSICS stays in each scheme's
// computeColumns(); only the common scaffolding is shared.
// ============================================================================
namespace IceSchemeCommon {

    // Hard cap on a precipitation flux [kg/(m2*s)] (~260 mm/d, well above any
    // physical precip). Backstop against the ∝Snow riming/deposition runaways.
    constexpr double P_max_flux = 3.0e-3;

    // c43/c13 quadratic edge extrapolation of ONE precipitation-flux field over
    // the top (radial), latitude (theta) and longitude (phi) boundaries. The phi
    // seam optionally gets periodicity averaging: rain uses it in every scheme;
    // snow/graupel use it in the >=2-category schemes but NOT in ZeroCat, so it
    // is a per-field flag to preserve each scheme's exact behaviour.
    inline void extrapolateBC(cAtmosphereModel& m, Array& P, bool periodic_avg) {
        // Top boundary
        #pragma omp parallel for collapse(2)
        for(int j = 0; j < m.jm; j++)
            for(int k = 0; k < m.km; k++)
                P.x[m.im-1][j][k] = m.c43 * P.x[m.im-2][j][k] - m.c13 * P.x[m.im-3][j][k];

        // Latitude (theta) boundaries
        #pragma omp parallel for collapse(2)
        for(int k = 0; k < m.km; k++)
            for(int i = 0; i < m.im; i++){
                P.x[i][0][k]      = m.c43 * P.x[i][1][k]      - m.c13 * P.x[i][2][k];
                P.x[i][m.jm-1][k] = m.c43 * P.x[i][m.jm-2][k] - m.c13 * P.x[i][m.jm-3][k];
            }

        // Longitude (phi) boundaries — von Neumann + optional periodicity average
        #pragma omp parallel for collapse(2)
        for(int i = 0; i < m.im; i++)
            for(int j = 0; j < m.jm; j++){
                P.x[i][j][0]      = m.c43 * P.x[i][j][1]      - m.c13 * P.x[i][j][2];
                P.x[i][j][m.km-1] = m.c43 * P.x[i][j][m.km-2] - m.c13 * P.x[i][j][m.km-3];
                if(periodic_avg)
                    P.x[i][j][0] = P.x[i][j][m.km-1] =
                        (P.x[i][j][0] + P.x[i][j][m.km-1]) * 0.5;
            }
    }

    // Fill the sub-surface land cells (i below i_topography) of ONE field with the
    // surface (i_mount) value. Each field independently copies its own ground value
    // downward, so calling this per field is equivalent to the old snapshot-all form.
    inline void fillTopography(cAtmosphereModel& m, Array& F) {
        #pragma omp parallel for collapse(2)
        for(int j = 0; j < m.jm; j++)
            for(int k = 0; k < m.km; k++){
                int i_mount = m.i_topography[j][k];
                double v = F.x[i_mount][j][k];
                for(int i = i_mount - 1; i >= 0; i--)
                    if(AtomUtils::is_land(m.h, i, j, k)) F.x[i][j][k] = v;
            }
    }

    // Saturation mass fraction over water / ice [kg/kg].
    //
    // Was the Magnus form with the dilute conversion ep*E/(p-E). Both parts fail here:
    // Magnus is calibrated to ~320 K against a column spanning 1500 K, and the dilute
    // conversion assumes water is a trace when it is 67 % of the mass. Note the old form
    // did not even guard p <= E, so it returned a NEGATIVE saturation humidity wherever
    // the extrapolated E exceeded the local pressure.
    inline double qSatWater(cAtmosphereModel& m, double t_u, int i, int j, int k) {
        const double M_other = AtmMixture::M_nonwater(m.c.x[i][j][k], m.co2.x[i][j][k],
                                                     m.m_comp.M_bg);
        return SaturationH2O::saturationMassFraction(
                   SaturationH2O::saturationPressure(t_u), m.p_stat.x[i][j][k], M_other);
    }
    inline double qSatIce(cAtmosphereModel& m, double t_u, int i, int j, int k) {
        const double M_other = AtmMixture::M_nonwater(m.c.x[i][j][k], m.co2.x[i][j][k],
                                                     m.m_comp.M_bg);
        return SaturationH2O::saturationMassFraction(
                   SaturationH2O::sublimationPressure(t_u), m.p_stat.x[i][j][k], M_other);
    }

    // ------------------------------------------------------------------------
    // ATHAD: can a condensed phase exist in this cell AT ALL?
    //
    // Two ways it cannot, and none of the four ice schemes checked for either — a grep
    // for the critical temperature across all of them, and across MoistConvection,
    // returned nothing:
    //
    //   SUPERCRITICAL (T >= 647.096 K). Liquid and vapour are one phase. There is no
    //     droplet, no surface tension, nothing to nucleate on, nothing to fall.
    //   SUPERHEATED (p_sat(T) > p). Subcritical, but the saturation pressure exceeds the
    //     total local pressure, so the vapour cannot reach saturation however much of it
    //     there is. qSatWater returns 1 in both cases, which is what "no saturation
    //     limit" means.
    //
    // Latent on Earth: every terrestrial cell is far below 647 K and far above its own
    // saturation pressure, so the question never arises. On ATHAD it is the state of the
    // column from the ground to ~240 km — everything below the cloud deck. The schemes
    // were running their full warm and cold microphysics there, and holding condensate
    // that had sedimented down from the one level that can genuinely condense.
    inline bool canCondense(cAtmosphereModel& m, double t_u, int i, int j, int k) {
        if (t_u >= AtmMixture::T_CRIT_H2O) return false;
        return qSatWater(m, t_u, i, j, k) < 1.0;
    }

    // Enforce it. Condensate that finds itself in such a cell — advected in, or fallen in
    // from the deck above — evaporates completely and immediately, and the microphysical
    // sources and precipitation fluxes have nothing to act on.
    //
    // Water is conserved: the condensate goes back into the vapour. Energy is conserved
    // where the phase change is real, i.e. below the critical point, where evaporation
    // absorbs L(T) and cools the cell. ABOVE the critical point there is no phase change
    // and no latent heat — the "condensate" was never a separate phase, so the mass moves
    // to c and the temperature is untouched. Watson's L(T) is not defined there anyway.
    inline void evaporateWhereImpossible(cAtmosphereModel& m, double t_u, int i, int j, int k) {
        const double q_cond = std::max(0.0, m.cloud.x[i][j][k])
                            + std::max(0.0, m.ice.x[i][j][k])
                            + std::max(0.0, m.gr.x[i][j][k]);

        if (q_cond > 0.0) {
            // The vapour has a physical ceiling: the mass fractions sum to 1 and the
            // background is carried as the remainder 1 - c - co2, so c > 1 - co2 is a
            // negative background mass. Nothing in the inherited code enforced it, and
            // AtmMixture::split() renormalises defensively, so the violation showed only
            // as a gas constant that had quietly saturated. Returning the condensate to
            // the vapour must not create one: take it up to the ceiling and no further.
            // Item 57: with the carrier renormalised, q_v + q_CO2 + q_bg == 1 - q_cond for ANY
            // q_v — the carrier shrinks as water grows — so the old c <= 1 - co2 ceiling has no
            // meaning under that convention and 1.0 is the real bound on a mass fraction. It
            // mattered more here than in ATHAD: with co2 pinned at 0.6233 the ceiling sat at
            // 0.3767, only 9 % above the sea-surface water content.
            const double c_max  = AtmMixture::co2_dilute()
                                ? 1.0 : std::max(0.0, 1.0 - m.co2.x[i][j][k]);
            const double c_new  = std::min(c_max, m.c.x[i][j][k] + q_cond);
            const double q_used = std::max(0.0, c_new - m.c.x[i][j][k]);

            m.c.x[i][j][k]     = c_new;
            m.cloud.x[i][j][k] = 0.0;
            m.ice.x[i][j][k]   = 0.0;
            m.gr.x[i][j][k]    = 0.0;

            if (t_u < AtmMixture::T_CRIT_H2O) {
                const double cp_loc = AtmMixture::cp_of(m.c.x[i][j][k], m.co2.x[i][j][k],
                                                        t_u, m.m_comp.M_bg);
                if (cp_loc > 0.0) {
                    // q_used, not q_cond: only the mass that actually became vapour
                    // absorbed its latent heat.
                    double T_new = t_u - SaturationH2O::latentHeat(t_u) * q_used / cp_loc;
                    if (T_new < 1.0) T_new = 1.0;                  // defensive only
                    m.t.x[i][j][k] = T_new / m.t_0;
                }
            }
        }

        // Precipitation falling into such a cell evaporates too, and its mass must go
        // somewhere. The first version simply zeroed the fluxes, which DESTROYED the
        // water: the deck at 243 km condensed, the ice scheme autoconverted it to rain,
        // the rain fell one level into the superheated band, and it vanished — leaving a
        // column with vapour everywhere and no cloud, ice or rain anywhere, which is
        // exactly what the output showed.
        //
        // A downward flux P [kg/(m2 s)] falling at speed v corresponds to a mass
        // concentration P/v [kg/m3] in the air it is passing through, hence a mass
        // fraction P/(v*rho). Convert, then clear. The fall speeds are the ones the
        // schemes' own residence times are built from (dt_rain_dim = step/1.6,
        // dt_snow_dim = step/0.96).
        {
            constexpr double v_rain = 1.6, v_snow = 0.96, v_graupel = 2.0;   // [m/s]
            const double rho = std::max(1.0e-6, m.r_humid.x[i][j][k]);
            const double q_p = std::max(0.0, m.P_rain.x[i][j][k])    / (v_rain    * rho)
                             + std::max(0.0, m.P_snow.x[i][j][k])    / (v_snow    * rho)
                             + std::max(0.0, m.P_graupel.x[i][j][k]) / (v_graupel * rho);
            if (q_p > 0.0) {
                const double c_max = AtmMixture::co2_dilute()          // item 57, see above
                                   ? 1.0 : std::max(0.0, 1.0 - m.co2.x[i][j][k]);
                m.c.x[i][j][k] = std::min(c_max, m.c.x[i][j][k] + q_p);
            }
        }

        m.P_rain.x[i][j][k]        = 0.0;
        m.P_snow.x[i][j][k]        = 0.0;
        m.P_graupel.x[i][j][k]     = 0.0;
        m.Precipitation.x[i][j][k] = 0.0;
        m.S_v.x[i][j][k]   = 0.0;
        m.S_c.x[i][j][k]   = 0.0;
        m.S_i.x[i][j][k]   = 0.0;
        m.S_r.x[i][j][k]   = 0.0;
        m.S_s.x[i][j][k]   = 0.0;
        m.S_g.x[i][j][k]   = 0.0;
        m.S_c_c.x[i][j][k] = 0.0;
    }

    // ---- Vapour -> ice -> snow throttle (Seifert-Beheng / COSMO; TwoCat's) ----
    // The stable way to make snow: grow it through the CLOUD-ICE reservoir, not
    // directly from vapour. Deposition onto ice is SUPERSATURATION-LIMITED
    // (S_i_dep proportional to c - q_Ice, and to the ice particle population
    // N_i*m_i^(1/3)), and the ice->snow autoconversions (S_i_au aggregation,
    // S_d_au deposition) are bounded — so snow cannot run away. This is exactly
    // what OneCat lacked (it deposited vapour straight onto snow proportional to
    // Snow^0.58, an unbounded feedback). TwoCat carries the same physics inline;
    // sharing it here lets OneCat inherit the throttle.
    // Also returns the ice particle number density N_i and mean mass m_i, which the
    // callers need for their own nucleation / ice-rain-collection terms.
    struct IceSnowRates { double S_i_dep = 0.0, S_i_au = 0.0, S_d_au = 0.0,
                                 N_i = 0.0, m_i = 1.0e-9; };

    inline IceSnowRates depositionThrottle(cAtmosphereModel& m, int i, int j, int k,
                                           double t_u, double q_Ice, double dt_snow_dim) {
        constexpr double N_i_0   = 1.0e2;    // 1/m3
        constexpr double m_i_0   = 1.0e-12;  // kg
        constexpr double m_i_max = 1.0e-9;   // kg
        constexpr double m_s_0   = 3.0e-9;   // kg (snow-size threshold)
        constexpr double c_i_dep = 1.3e-5;   // m3/(s*kg^1/3)
        constexpr double c_i_au  = 1.0e-3;   // 1/s
        constexpr double t_hn    = 236.15;   // K (-37C) homogeneous freezing floor

        IceSnowRates r;
        const double c   = m.c.x[i][j][k];
        const double ice = m.ice.x[i][j][k];

        // ice number density + mean particle mass (only in the mixed-phase band)
        double N_i = 0.0, m_i = m_i_max;
        if(t_u <= m.t_0 && t_u > t_hn){
            N_i = N_i_0 * std::exp(0.2 * (m.t_0 - t_u));
            m_i = std::min(m.r_humid.x[i][j][k] * ice / N_i, m_i_max);
            m_i = std::max(m_i_0, std::min(m_i, m_i_max));
        }

        // deposition (supersaturation-limited) / sublimation (ice-mass-limited)
        if(c > q_Ice)      r.S_i_dep = c_i_dep * N_i * std::pow(m_i, 1.0/3.0) * (c - q_Ice);
        else if(c < q_Ice) r.S_i_dep = std::max(-ice / dt_snow_dim, (c - q_Ice) / dt_snow_dim);

        // ice -> snow: aggregation (bounded) + depositional autoconversion
        if(ice > 0.0)       r.S_i_au = std::max(c_i_au * ice, 0.0);
        if(r.S_i_dep > 0.0) r.S_d_au = r.S_i_dep / (1.5 * (std::pow(m_s_0 / m_i, 2.0/3.0) - 1.0));
        r.N_i = N_i;  r.m_i = m_i;
        return r;
    }
}

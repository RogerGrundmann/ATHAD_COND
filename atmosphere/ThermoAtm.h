#pragma once

#include "MixtureAtm.h"
#include "SaturationH2O.h"
#include "cAtmosphereModel.h"
#include "Utils.h"

#include <vector>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AtomUtils;

// Physical constants for precipitable water column integration
namespace PrecipWaterConstants {
    constexpr double HPA_TO_PA        = 100.0;
    constexpr double MIN_SAFE_TEMP    = 100.0;   // [K]
    constexpr double MIN_SAFE_PRESSURE = 1e-6;   // [hPa]
}

class ThermoAtm {
public:
    explicit ThermoAtm(cAtmosphereModel& model)
        : m(model)
    {}

    // ------------------------------------------------------------------
    void latentSensibleHeat()
    {
        using namespace std;
        cout << endl << endl << endl << "      LatentSensibleHeat" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        const double inv_2dr    = 1.0 / (2.0 * m.dr);
        const double inv_2dthe  = 1.0 / (2.0 * m.dthe);
        const double inv_2dphi  = 1.0 / (2.0 * m.dphi);
        const double inv_sqrt3  = 1.0 / sqrt(3.0);

        std::vector<double> sinthe_table(m.jm);
        for (int j = 0; j < m.jm; j++) {
            sinthe_table[j] = sin(m.the.z[j]);
            if (std::abs(sinthe_table[j]) < 1e-10)
                sinthe_table[j] = 1e-10;
        }

        #pragma omp parallel for collapse(2)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                m.Q_Latent.x[0][j][k]   = 0.0;
                m.Q_Sensible.x[0][j][k] = 0.0;
            }
        }

        #pragma omp parallel for collapse(2) schedule(dynamic, 4)
        for (int j = 1; j < m.jm-1; j++) {
            for (int k = 1; k < m.km-1; k++) {

                double sinthe = sinthe_table[j];

                for (int i = 1; i < m.im-1; i++) {

                    if (is_land(m.h, i, j, k)) {
                        m.Q_Latent.x[i][j][k]   = 0.0;
                        m.Q_Sensible.x[i][j][k] = 0.0;
                        continue;
                    }

                    double rm           = m.rad.z[i];
                    double exp_rm       = 1.0 / (rm + 1.0);
                    double inv_rm       = 1.0 / rm;
                    double inv_rmsinthe = 1.0 / (rm * sinthe);

                    // Q_Latent is owned by RHS_Atm_Turb (the signed advective
                    // tendency × coeff_L written at every saturated cell). The
                    // earlier ThermoAtm path here wrote an unsigned magnitude
                    // diagnostic (lv · |∇q_sat|) every "moist" iter, racing the
                    // RHS writer on alternate iters and producing a 2Δt sawtooth
                    // in the field that confused diagnostics and (combined with
                    // the now-removed RHS_Atm_Turb.cpp:758 source term) drove the
                    // iter-358 NaN at 62°N upper-tropo.  See
                    // [[project-iter358-62n-dry-nan]] and
                    // [[project-precip-chain-fixes]].

                    // --- sensible heat ---
                    double dtdr   = (m.t.x[i+1][j][k] - m.t.x[i-1][j][k])
                                    * inv_2dr * exp_rm;
                    double dtdthe = (m.t.x[i][j+1][k] - m.t.x[i][j-1][k])
                                    * inv_2dthe * inv_rm;
                    double dtdphi = (m.t.x[i][j][k+1] - m.t.x[i][j][k-1])
                                    * inv_2dphi * inv_rmsinthe;

                    m.Q_Sensible.x[i][j][k] = m.cp_l * m.t_0 * inv_sqrt3
                        * sqrt(dtdr*dtdr + dtdthe*dtdthe + dtdphi*dtdphi);

                }  // i
            }  // k
        }  // j

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for LatentSensibleHeat\n", elapsed.count() * 1e-9);
        cout << "      LatentSensibleHeat ended" << endl;
    }

    // ------------------------------------------------------------------
    void waterVapourEvaporation()
    {
        using namespace std;
        cout << "\n\n\n      WaterVapourEvaporation  (model: " << m.evap_model << ")" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        const double conv_factor  = 8.64e4;                             // [mm/s] -> [mm/d]
        const double hPa_to_mmHg  = 0.750062;                           // [mmHg/hPa]
        const double ms_to_kmh    = 3.6;                                // [(km/h)/(m/s)]
        // Meyer (1915): E[mm/month] = C·(1 + W/16)·(e_s − e_a), with the deficit in
        // mmHg and C the open-water coefficient ≈ 11 (deep/large water, e.g. ocean) or
        // 15 (shallow ponds).  C is unit-identical for E in mm/month & deficit in mmHg
        // because inch→mm and inHg→mmHg share the 25.4 factor.  The previous value 0.36
        // was ~30× too small, giving ~0.18 mm/d (≈66 mm/yr) instead of the realistic
        // ~5 mm/d (≈1800 mm/yr) and far below the Dalton/Rohwer diagnostics.
        const double K_Meyer      = 11.0;                               // Meyer (1915) deep open-water coefficient [mm/month/mmHg]

        // Vertical spread: distribute c_eq over n_spread+1 levels with exp(-i) weights,
        // normalised so that sum_{i=0}^{n_spread} exp(-i) = 1 (moisture conserved).
//        const int    n_spread = 5;
        const int    n_spread = 3;
//        const int    n_spread = 2;
        const double r        = std::exp(-1.0);
        const double w_norm   = (1.0 - r) / (1.0 - std::pow(r, n_spread + 1)); // 1/sum

        // Determine which formula drives the surface-humidity (c) update.
        // All three are always computed for diagnostic output.
        enum class EvapModel { Dalton, Meyer, Rohwer } active;

        if      (m.evap_model == "Meyer")  active = EvapModel::Meyer;
        else if (m.evap_model == "Rohwer") active = EvapModel::Rohwer;
        else                               active = EvapModel::Dalton;

        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                if (is_land(m.h, 0, j, k)) {
                    m.Evaporation_Dalton.y[j][k] = 0.0;
                    m.Evaporation_Meyer.y[j][k]  = 0.0;
                    m.Evaporation_Rohwer.y[j][k] = 0.0;
                    m.Evaporation.y[j][k]        = 0.0;
                    continue;
                }

                // ATHAD: no evaporation from a supercritical surface.
                //
                // Dalton, Meyer and Rohwer all drive evaporation from the saturation
                // DEFICIT (e_s - e_a) across a liquid-vapour interface. Above the critical
                // point there is no interface and no e_s: the surface and the atmosphere
                // are one fluid. Evaluated anyway at 1500 K, the Magnus e_s returns
                // ~1.2e7 hPa against a 250 bar column, and the resulting c_eq drove the
                // surface water mass fraction to 21 — twenty times the mass present.
                //
                // Guarding here rather than clamping downstream keeps the surface cell's
                // mixture properties (R, cp, density) meaningful, which matters because
                // the whole column is anchored on them.
                //
                // ATHAD_COND KEEPS THIS GUARD — the surface here is 513 K, well below the
                // critical point, so it never fires — and adds the branch below, which is
                // the one that runs. Removing the guard rather than passing it would have
                // been the obvious move and the wrong one: it is still correct, and a
                // 27 bar corner of the stated input range boils at 501 K.
                if (m.t.x[0][j][k] * m.t_0 >= AtmMixture::T_CRIT_H2O) {
                    m.Evaporation_Dalton.y[j][k] = 0.0;
                    m.Evaporation_Meyer.y[j][k]  = 0.0;
                    m.Evaporation_Rohwer.y[j][k] = 0.0;
                    m.Evaporation.y[j][k]        = 0.0;
                    continue;
                }

                // ==============================================================
                // ATHAD_COND: THE SEA IS A SATURATION BOUNDARY, NOT A DALTON FLUX.
                //
                // The three formulas below this block are Earth fits: E = coeff * (e_s -
                // e_a) with coefficients calibrated on lakes and oceans at ~1 bar and
                // ~288 K, where the saturation deficit is a few tens of hPa. Here
                // p_sat(513 K) = 33.5 bar, so the deficit at the initial state is of order
                // 33 bar = 24 700 mmHg, and Meyer's K*(1+u/16)*deficit returns of order
                // 9 000 mm/day against Earth's ~5. Enabling them by simply dropping the
                // guard above would reproduce, from the other side, exactly the failure the
                // guard was written to stop: c driven to a mass fraction larger than the
                // mass present.
                //
                // The physics does not need them. A liquid ocean in contact with the air
                // above it holds that air at its own vapour pressure; the flux is whatever
                // it takes to maintain that, and the empirical coefficient is a statement
                // about how fast, not about where it ends up. So the surface layer is
                // RELAXED TOWARD SATURATION with the same exponential vertical spread the
                // Earth branch uses (n_spread = 3, weights summing to 1), and w_norm sets
                // the rate.
                //
                // The three empirical numbers are still computed and still printed, marked
                // as Earth-calibrated and unused, so the gap between them and the boundary
                // condition stays visible instead of being deleted. If one of them ever
                // becomes plausible here, that is worth seeing.
                {
                    const double T_s   = std::max(180.0, m.t.x[0][j][k] * m.t_0);   // [K]
                    const double p_s   = m.p_stat.x[0][j][k];                       // [hPa]
                    const double M_nw  = AtmMixture::M_nonwater(m.c.x[0][j][k],
                                             m.co2.x[0][j][k], m.m_comp.M_bg);
                    const double E_s   = SaturationH2O::saturationPressureAuto(T_s);
                    const double q_sea = SaturationH2O::saturationMassFraction(E_s, p_s, M_nw);

                    // Wind speed only enters the diagnostics; the boundary condition is
                    // thermodynamic.
                    const double vel = sqrt((m.u.x[0][j][k] * m.u.x[0][j][k]
                                           + m.v.x[0][j][k] * m.v.x[0][j][k]
                                           + m.w.x[0][j][k] * m.w.x[0][j][k]) / 3.0) * m.u_0;
                    const double u_kmh_d = vel * 3.6;
                    // Exact q -> e, the inverse of saturationMassFraction. Written inline
                    // here first; SaturationH2O::vapourPressureFromMassFraction is that same
                    // expression, so the two cannot drift apart any more.
                    const double e_air   = (q_sea > 0.0 && q_sea < 1.0)
                                         ? SaturationH2O::vapourPressureFromMassFraction(
                                               m.c.x[0][j][k], p_s, M_nw)
                                         : 0.0;                                     // [hPa], exact
                    const double sd_d    = std::max(0.0, E_s - e_air);              // [hPa]

                    m.Evaporation_Dalton.y[j][k] =
                        AtomUtils::C_Dalton(0, j, k, m.coeff_Dalton, m.u_0, m.u, m.v, m.w) * 24.0 * sd_d;
                    m.Evaporation_Meyer.y[j][k]  =
                        K_Meyer * hPa_to_mmHg * (1.0 + u_kmh_d / 16.0) / 30.0 * sd_d;
                    m.Evaporation_Rohwer.y[j][k] =
                        0.771 * (1.465 - 0.000732 * p_s * hPa_to_mmHg)
                              * (0.44 + 0.0733 * u_kmh_d) * hPa_to_mmHg * sd_d;
                    m.Evaporation.y[j][k] = m.Evaporation_Dalton.y[j][k];

                    // Surface layer toward saturation, then the exponential spread above it.
                    m.c_fix.y[j][k] = m.c.x[0][j][k];
                    m.c.x[0][j][k]  = m.c_fix.y[j][k] + (q_sea - m.c_fix.y[j][k]) * w_norm;

                    for (int i = 1; i <= n_spread && i < m.im; i++) {
                        const double weight = std::pow(r, i) * w_norm;
                        const double T_i    = std::max(180.0, m.t.x[i][j][k] * m.t_0);
                        const double p_i    = m.p_stat.x[i][j][k];
                        const double M_nw_i = AtmMixture::M_nonwater(m.c.x[i][j][k],
                                                  m.co2.x[i][j][k], m.m_comp.M_bg);
                        const double q_s_i  = SaturationH2O::saturationMassFractionAt(
                                                  T_i, p_i, M_nw_i);
                        m.c.x[i][j][k] = std::min(
                            m.c.x[i][j][k] + (q_sea - m.c_fix.y[j][k]) * weight, q_s_i);
                    }
                    continue;
                }

                double c_Dalton    = AtomUtils::C_Dalton(0, j, k, m.coeff_Dalton, m.u_0, m.u, m.v, m.w);
                                                                        // [mm/(h*hPa)]
                double p_stat_0jk  = m.p_stat.x[0][j][k];               // [hPa]
                // Floor the surface temperature at 180 K (below Earth's coldest-ever ~184 K, so
                // it never touches real physics) to backstop the marginal Gulf-of-Alaska crash:
                // an undamped i=0 coastal-surface T oscillation drove t.x[0] to ~160 K, making
                // E_sat/c_eq in the evaporation formula singular → c NaN. See
                // [[project_upper_velocity_secular_growth]].
                double t_u_base    = std::max(180.0, m.t.x[0][j][k] * m.t_0);  // [K], floored

                // EVERYTHING FROM HERE DOWN IS THE INHERITED EARTH PATH AND IS UNREACHABLE
                // IN THIS FORK. The ATHAD_COND branch above ends in `continue` for every
                // subcritical cell, and the supercritical guard above that catches the rest,
                // so no cell arrives here. It is kept because the guard is a temperature test
                // and the 27 bar corner of the stated input range has no liquid ocean at all.
                //
                // It is made consistent rather than left as found: it used the scalar config
                // ep = R_Air/R_v, and R_Air here is the background EXCLUDING CO2 (N2,
                // 28.014 g/mol). The "other" gas in a saturation formula is everything that is
                // not water — CO2 and the background, M_nonwater = 42.888 g/mol at this sea
                // surface. The scalar therefore stands in 0.6431 for a true 0.4201, which is
                // the 29 % error test/cond_column_selftest.cpp:158-165 describes. That test
                // asserts against a value it computes itself, so it never covered this code;
                // what protects the model is the branch above, not the test.
                //
                // Converting it is a trap removal, NOT a bug fix: it changes no output, and it
                // must not be quoted as one. If the guard above is ever loosened, this path
                // becomes live with the right molar ratio instead of the wrong one.
                const double M_nw_0 = AtmMixture::M_nonwater(m.c.x[0][j][k],
                                          m.co2.x[0][j][k], m.m_comp.M_bg);
                double precip_term = conv_factor * m.Precipitation.x[0][j][k]; // [mm/d]

                double vel_ms = sqrt((m.u.x[0][j][k] * m.u.x[0][j][k]
                                    + m.v.x[0][j][k] * m.v.x[0][j][k]
                                    + m.w.x[0][j][k] * m.w.x[0][j][k]) / 3.0) * m.u_0; // [m/s]
                double u_kmh  = vel_ms * ms_to_kmh;                     // [km/h]
                double p_mmHg = p_stat_0jk * hPa_to_mmHg;               // [mmHg]

                double E_sat  = (t_u_base >= m.t_0)                     // [hPa]
                    ? SaturationH2O::saturationPressure(t_u_base)
                    : SaturationH2O::sublimationPressure(t_u_base);

                // Per-formula coefficients [mm/(d*hPa)]: E = coeff * sat_deficit_hPa
                //   Dalton: from wind-dependent mass-transfer coefficient
                //   Meyer (1915): K_M*(1+u/16)/30 [mm/month/mmHg -> mm/d/hPa]
                //   Rohwer (1931): 0.771*(1.465-0.000732*p)*(0.44+0.0733*u) [mm/d/mmHg -> mm/d/hPa]
                double coeff_D = c_Dalton * 24.0;
                double coeff_M = K_Meyer * hPa_to_mmHg * (1.0 + u_kmh / 16.0) / 30.0;
                double coeff_R = 0.771 * (1.465 - 0.000732 * p_mmHg)
                               * (0.44  + 0.0733   * u_kmh) * hPa_to_mmHg;

                // Calm conditions (zero wind): skip c update, Dalton = 0.
                // Meyer and Rohwer retain their still-air (u=0) terms.
                if (c_Dalton <= 0.0) {
                    double e  = SaturationH2O::vapourPressureFromMassFraction(
                                    m.c.x[0][j][k], p_stat_0jk, M_nw_0);   // [hPa]  exact, at c_Dalton = 0
                    double sd = std::max(0.0, E_sat - e);               // [hPa]  saturation deficit at c_Dalton = 0
                    m.Evaporation_Dalton.y[j][k] = 0.0;
                    m.Evaporation_Meyer.y[j][k]  = coeff_M * sd;        // [mm/d]
                    m.Evaporation_Rohwer.y[j][k] = coeff_R * sd;        // [mm/d]
                    m.Evaporation.y[j][k] = (active == EvapModel::Meyer)  ? m.Evaporation_Meyer.y[j][k]
                                          : (active == EvapModel::Rohwer) ? m.Evaporation_Rohwer.y[j][k]
                                          :                                 m.Evaporation_Dalton.y[j][k];
                    // Calm ocean cell: the wind-driven evaporation FLUX is zero (above), but the
                    // surface-layer vapour CONCENTRATION still equilibrates toward saturation (a
                    // windy no-precip cell converges to ~c_sat too). Relax c.x[0] toward c_sat with
                    // the same w_norm as the windy branch below, instead of leaving the
                    // RK4-unintegrated i=0 layer at a stale near-zero value — which showed as
                    // near-zero water-vapour patches over calm ocean (with normal "spots" at windy
                    // cells) in zonal cross-sections. Concentration only; the flux stays wind-limited.
                    double c_sat_calm = SaturationH2O::saturationMassFraction(
                                            E_sat, p_stat_0jk, M_nw_0);     // [kg/kg] exact
                    m.c_fix.y[j][k] = m.c.x[0][j][k];
                    m.c.x[0][j][k]  = m.c_fix.y[j][k] + (c_sat_calm - m.c_fix.y[j][k]) * w_norm;
                    continue;
                }

                // Saturation specific humidity — exact mole-to-mass conversion on M_nonwater.
                // The form that stood here, ep*E/(p-(1-ep)E), is exact only in the molar
                // ratio it is given; it was given the CO2-free one. Both the ratio and the
                // conversion are now the model's shared ones.
                double c_sat = SaturationH2O::saturationMassFraction(
                                   E_sat, p_stat_0jk, M_nw_0);          // [kg/kg]

                // Active formula solves for the c_eq that balances precipitation:
                //   E_active = coeff_active * (E_sat - e(c_eq)) = precip_term
                // so the equilibrium VAPOUR PRESSURE is e_eq = E_sat - precip_term/coeff,
                // and c_eq is that pressure converted to a mass fraction by the same exact
                // route as c_sat. Solving for c directly (the old ep/p * (...) line) is the
                // dilute inverse and disagrees with the c_sat two lines above by 29 % here.
                // If e_eq <= 0, heavy rain saturates the air -> clamp to c_sat.
                double coeff_active = (active == EvapModel::Meyer)  ?   coeff_M
                                    : (active == EvapModel::Rohwer) ?   coeff_R
                                    :                                   coeff_D;

                double e_eq = (coeff_active != 0.0)
                            ? E_sat - precip_term / coeff_active        // [hPa]
                            : E_sat;
                double c_eq = (e_eq > 0.0)
                            ? SaturationH2O::saturationMassFraction(e_eq, p_stat_0jk, M_nw_0)
                            : c_sat;
                c_eq = std::min(c_eq, c_sat);

                m.c_fix.y[j][k] = m.c.x[0][j][k];

                // Evaporations from current (pre-update) saturation deficit
                double e_cur  = SaturationH2O::vapourPressureFromMassFraction(
                                    m.c_fix.y[j][k], p_stat_0jk, M_nw_0); // [hPa] exact current vapour pressure
                double sd_cur = std::max(0.0, E_sat - e_cur);           // [hPa]  current saturation deficit
 
                m.Evaporation_Dalton.y[j][k] = coeff_D * sd_cur;        // [mm/d]
                m.Evaporation_Meyer.y[j][k]  = coeff_M * sd_cur;        // [mm/d]
                m.Evaporation_Rohwer.y[j][k] = coeff_R * sd_cur;        // [mm/d]
                m.Evaporation.y[j][k] = (active == EvapModel::Meyer)  ? m.Evaporation_Meyer.y[j][k]
                                      : (active == EvapModel::Rohwer) ? m.Evaporation_Rohwer.y[j][k]
                                      :                                 m.Evaporation_Dalton.y[j][k];

                // c update: distribute c_eq across surface + n_spread levels above,
                // conserving total moisture (weights sum to 1).
                m.c.x[0][j][k] = m.c_fix.y[j][k] + (c_eq - m.c_fix.y[j][k]) * w_norm;  // i=0: weight = exp(0)*w_norm

                if (c_eq > 0.0) {
                    for (int i = 1; i <= n_spread; i++) {
                        double weight  = std::pow(r, i) * w_norm;
                        double t_i     = m.t.x[i][j][k] * m.t_0;
                        double p_i     = m.p_stat.x[i][j][k];
                        double E_i     = (t_i >= m.t_0)
                            ? SaturationH2O::saturationPressure(t_i)
                            : SaturationH2O::sublimationPressure(t_i);
                        double M_nw_i  = AtmMixture::M_nonwater(m.c.x[i][j][k],
                                             m.co2.x[i][j][k], m.m_comp.M_bg);
                        double c_sat_i = SaturationH2O::saturationMassFraction(E_i, p_i, M_nw_i);
                        m.c.x[i][j][k] = std::min(m.c.x[i][j][k] + c_eq * weight, c_sat_i);
                    }
                }


/*
                cout.precision(8);
                cout.setf(ios::fixed);
                if ((j == 90) && (k == 180)) cout << endl
                    << "  WaterVapourEvaporation" << endl
                    << "  j = " << j << "  k = " << k << endl
                    << "  p_stat_0jk = "     << p_stat_0jk
                    << "  p_mmHg = "         << p_mmHg << endl
                    << "  c_fix = "          << m.c_fix.y[j][k] * 1e3
                    << "  c_eq = "           << c_eq * 1e3
                    << "  c_sat = "          << c_sat * 1e3
                    << "  c = "              << m.c.x[0][j][k] * 1e3 << endl
                    << "  E_sat = "          << E_sat
                    << "  e_cur = "          << e_cur
                    << "  sd_cur = "         << sd_cur << endl
                    << "  c_Dalton = "       << c_Dalton
                    << "  coeff_D = "        << coeff_D
                    << "  coeff_M = "        << coeff_M
                    << "  coeff_R = "        << coeff_R << endl
                    << "  Evap_Dalton = "    << m.Evaporation_Dalton.y[j][k]
                    << "  Evap_Meyer = "     << m.Evaporation_Meyer.y[j][k]
                    << "  Evap_Rohwer = "    << m.Evaporation_Rohwer.y[j][k]
                    << "  Evap = "           << m.Evaporation.y[j][k]
                    << "  Prec = "           << precip_term << endl;
*/
            }  // k
        }  // j

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for WaterVapourEvaporation\n", elapsed.count() * 1e-9);
        cout << "      WaterVapourEvaporation ended" << endl;
    }

    // ------------------------------------------------------------------
    void standAtm_DewPoint_HumidRel()
    {
        using namespace std;
        cout << endl << endl << endl << "      StandAtm_DewPoint_HumidRel" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        std::vector<double> height_table(m.im);
        for (int i = 0; i < m.im; i++)
            height_table[i] = m.get_layer_height(i);

        const double lapse_rate = 6.5e-3;
        const double inv_ep     = 1.0 / m.ep;

        std::vector<int> i_trop_table(m.jm);
        for (int j = 0; j < m.jm; j++) {
            i_trop_table[j] = (j == 90)
                ? (int)m.tropopause_layers[j-1]
                : (int)m.tropopause_layers[j];
        }

        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                int i_trop = i_trop_table[j];

                double TempStand_surface   = m.t.x[0][j][k] * m.t_0 - m.t_0;
                m.TempStand.x[0][j][k]    = TempStand_surface;

                // ATHAD: compute over the WHOLE column, not just up to the tropopause.
                //
                // The loop ran to i_trop and the block below copied the tropopause value to
                // every level above it. On Earth that is reasonable — the stratosphere is
                // dry and dull. Here the tropopause index sits near 200 km and the ONLY
                // levels where water can condense are above it, so the humidity and dew
                // point in the one region of interest were not computed at all: they were
                // the copied value from below, which is why the condensation level read a
                // few per cent RH while being genuinely supersaturated.
                for (int i = 0; i < m.im; i++) {

                    double t_u = m.t.x[i][j][k] * m.t_0;

                    double E = (t_u > m.t_0)
                        ? SaturationH2O::saturationPressure(t_u)
                        : SaturationH2O::sublimationPressure(t_u);

                    // ATHAD: the EXACT water partial pressure, e = x_H2O * p.
                    //
                    // Was e = c * p / ep, the dilute approximation. At 67 % water by mass
                    // it returns 243 600 hPa against a true 200 000 — and against a TOTAL
                    // pressure of 250 000, so the "partial" pressure was approaching the
                    // whole column. x_H2O_of does the mass-to-mole conversion properly.
                    const double x_v = AtmMixture::x_H2O_of(m.c.x[i][j][k],
                                                            m.co2.x[i][j][k], m.m_comp.M_bg);
                    double e           = x_v * m.p_stat.x[i][j][k];
                    bool   zero_vapour = (e <= 0.0);
                    if (zero_vapour) e = 1e-3;

                    m.TempStand.x[i][j][k]    = TempStand_surface - lapse_rate * height_table[i];

                    // Dew point by bisecting the IAPWS curve. The closed-form Magnus
                    // inverse it replaces also returned a value in DEGREES CELSIUS while
                    // everything around it is in kelvin — preserved here as the same
                    // convention (subtract t_0) so the ParaView output is unchanged.
                    m.TempDewPoint.x[i][j][k] = SaturationH2O::dewPoint(e) - m.t_0;

                    // ATHAD: report the RELATIVE HUMIDITY, not a clamp.
                    //
                    // Two defects. E is SaturationH2O::NO_SATURATION = -1 above the
                    // critical point, so e/E*100 came out around -2.4e7 through the whole
                    // supercritical column — the enormous values visible in ParaView. And
                    // the std::min(..., 100.0) hid genuine supersaturation, which on this
                    // column reaches ~900 % at the condensation level: the one number that
                    // would have shown the saturation adjustment was failing (README
                    // item 16) was capped at 100 the entire time.
                    //
                    // Now: the true ratio, uncapped, and NO_SATURATION where no saturation
                    // state exists at all — the same -1 sentinel the saturation module
                    // uses, so a supercritical cell is visibly distinct in ParaView from a
                    // dry one (0) rather than being a huge negative number.
                    m.HumidityRel.x[i][j][k]  = (E <= 0.0)
                        ? SaturationH2O::NO_SATURATION
                        : (zero_vapour ? 0.0 : e / E * 100.0);

                }  // i

                // TempStand is the International Standard Atmosphere profile — a
                // troposphere-only construction — so it alone keeps the copy above the
                // tropopause. The dew point and the relative humidity are now computed at
                // every level and are no longer overwritten here.
                double ts_trop = m.TempStand.x[i_trop][j][k];
                for (int i = i_trop + 1; i < m.im; i++)
                    m.TempStand.x[i][j][k] = ts_trop;

                // Below-terrain fill: copy from topography surface downward
                int    i_mount   = m.i_topography[j][k];
                double ts_mount  = m.TempStand.x[i_mount][j][k];
                double td_mount  = m.TempDewPoint.x[i_mount][j][k];
                double hr_mount  = m.HumidityRel.x[i_mount][j][k];

                for (int i = i_mount - 1; i >= 0; i--) {
                    if (is_land(m.h, i, j, k)) {
                        m.TempStand.x[i][j][k]    = ts_mount;
                        m.TempDewPoint.x[i][j][k] = td_mount;
                        m.HumidityRel.x[i][j][k]  = hr_mount;
                    }
                }
            }  // k
        }  // j

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for StandAtm_DewPoint_HumidRel\n", elapsed.count() * 1e-9);
        cout << "      StandAtm_DewPoint_HumidRel ended" << endl;
    }

    // ------------------------------------------------------------------
    void forces()
    {
        using namespace std;
        cout << endl << endl << endl << "      ATOM: Forces" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        const double inv_2dr    = 1.0 / (2.0 * m.dr);
        const double inv_2dthe  = 1.0 / (2.0 * m.dthe);
        const double inv_2dphi  = 1.0 / (2.0 * m.dphi);
        const double two_omega  = 2.0 * m.omega;
        const double omega2     = m.omega * m.omega;
        const double pres_coeff = 1e2 * m.p_0 / m.L_atm;
        const double r_Earth_m  = m.r_Earth * 1e3;

        std::vector<double> sinthe_table(m.jm);
        std::vector<double> costhe_table(m.jm);
        for (int j = 0; j < m.jm; j++) {
            sinthe_table[j] = sin(m.the.z[j]);
            costhe_table[j] = cos(m.the.z[j]);
        }

        #pragma omp parallel for collapse(2)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {

                double rm           = m.rad.z[i];
                double exp_rm       = 1.0 / (rm + 1.0);
                double inv_rm       = 1.0 / rm;
                double sinthe       = sinthe_table[j];
                double costhe       = costhe_table[j];
                double inv_rmsinthe = 1.0 / (rm * sinthe);
                double abs_sinthe   = fabs(sinthe);

                double rad_dist  = (double)i * m.L_atm * exp_rm;
                double rad_Earth = rad_dist + r_Earth_m;

                for (int k = 1; k < m.km-1; k++) {

                    double w_ijk   = m.w.x[i][j][k] * m.u_0;
                    // ATHAD: the DIAGNOSTIC must show what the dynamics actually applies.
                    //
                    // RHS_Atm_Turb uses the traditional approximation — it keeps only the
                    // Omega_r = Omega*cos(theta) component and drops the non-traditional
                    // ("Eotvos") sin(theta) pair unless ATOM_CORIOLIS_NONTRAD is set:
                    //     a_theta = +2*Omega*cos(theta)*w
                    //     a_phi   = -2*Omega*cos(theta)*v
                    // This diagnostic carried a +2*Omega*sin(theta)*u term the momentum
                    // equations do not have, and a radial term they drop, so the ParaView
                    // "Coriolis force" field showed a force that was never applied. Follow
                    // the same switch so the two cannot diverge.
                    //
                    // These signs agree with ATURAN 8b284cb / ATNEPT 024c37f once their
                    // opposite convention is accounted for — they store -a because their
                    // rhs subtracts, this file stores +a because its rhs adds.
                    const double nontrad_d = AtomUtils::coriolis_nontraditional() ? 1.0 : 0.0;
                    double Cor_r   = -nontrad_d * two_omega * sinthe * w_ijk;
                    double Cor_the =  two_omega * costhe * w_ijk;
                    double Cor_phi = -two_omega * (costhe * m.v.x[i][j][k]
                                        + nontrad_d * sinthe * m.u.x[i][j][k]) * m.u_0;

                    m.CoriolisForce.x[i][j][k] = m.Coriolis * m.r_air
                        * sqrt(Cor_r*Cor_r + Cor_the*Cor_the + Cor_phi*Cor_phi);

                    // Centrifugal acceleration points away from the ROTATION AXIS, whose
                    // distance is r*sin(theta) — so it vanishes at the poles and is largest
                    // at the equator. The previous factor (1 + |sin(theta)|) was maximal at
                    // the pole, where the true value is zero. Components:
                    //     a_r     = omega^2 * r * sin^2(theta)
                    //     a_theta = omega^2 * r * sin(theta)*cos(theta)
                    // ported from ATURAN 4201957 / ATJUP 8649675. Diagnostic only: this
                    // model's RHS carries no centrifugal term, the force being curl-free
                    // and absorbed into the pressure projection — which is what ATURAN
                    // found when it corrected the dynamical version.
                    const double cen_r   = omega2 * rad_Earth * sinthe * sinthe;
                    const double cen_the = omega2 * rad_Earth * sinthe * costhe;
                    m.CentrifugalForce.x[i][j][k] = m.centrifugal * m.r_air
                        * sqrt(cen_r * cen_r + cen_the * cen_the);

                    // Diagnostic buoyancy force for ParaView/Results — must match the
                    // perturbation-form body force applied in RHS_Atm.cpp (rhs_u), i.e.
                    // +g·ρ·(t − t̄(i)) referenced to the per-level mean t_ref_level[i].
                    const double t_ref_b = ((int)m.t_ref_level.size() == m.im)
                                           ? m.t_ref_level[i] : 1.0;
                    m.BuoyancyForce.x[i][j][k] = 1.0e-3 * m.buoyancy * m.r_humid.x[i][j][k] * m.g
                                                 * (m.t.x[i][j][k] - t_ref_b);

                    double dpdr   = (m.p_dyn.x[i+1][j][k] - m.p_dyn.x[i-1][j][k])
                                    * inv_2dr * exp_rm;
                    double dpdthe = (m.p_dyn.x[i][j+1][k] - m.p_dyn.x[i][j-1][k])
                                    * inv_2dthe * inv_rm;
                    double dpdphi = (m.p_dyn.x[i][j][k+1] - m.p_dyn.x[i][j][k-1])
                                    * inv_2dphi * inv_rmsinthe;

                    m.PresGradForce.x[i][j][k] =
                        - 1.0e-3 * sqrt(dpdr*dpdr + dpdthe*dpdthe + dpdphi*dpdphi) * pres_coeff;
                }  // k
            }  // j
        }  // i

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for Forces\n", elapsed.count() * 1e-9);
        cout << "      ATOM: Forces ended" << endl;
    }

    // ------------------------------------------------------------------
    void adjustTemperatureIC(double** t, int jm, int km)
    {
        const double inv_t0 = 1.0 / m.t_0;
        const int    k_half = (km - 1) / 2;

        #pragma omp parallel for
        for (int j = 0; j < jm; j++) {
            for (int k = 0; k < km; k++)
                t[j][k] = (t[j][k] + m.t_0) * inv_t0;
            t[j][k_half] = (t[j][k_half + 1] + t[j][k_half - 1]) * 0.5;
            m.temperature_NASA.y[j][k_half] =
                (m.temperature_NASA.y[j][k_half + 1] +
                 m.temperature_NASA.y[j][k_half - 1]) * 0.5;
        }
    }

    // ------------------------------------------------------------------
    void precipitableWater()
    {
        using namespace std;
        using namespace PrecipWaterConstants;

        cout << "\n\n\n      PrecipitableWater" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        #pragma omp parallel for collapse(2)
        for (int j = 0; j < m.jm; j++)
            for (int k = 0; k < m.km; k++)
                m.precipitable_water.y[j][k] = 0.0;

        #pragma omp parallel for collapse(2)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                double column_sum = 0.0;

                for (int i = 0; i < m.im - 1; i++) {

                    const double t_actual  = m.t.x[i][j][k] * m.t_0;
                    const double p_actual  = m.p_stat.x[i][j][k];
                    const double q_mixing  = m.c.x[i][j][k];

                    if (t_actual < MIN_SAFE_TEMP || p_actual < MIN_SAFE_PRESSURE) {
                        m.PrecipitableWaterLocal.x[i][j][k] = 0.0;
                        continue;
                    }

                    // Exact inverse, not the dilute q*p/ep: at this fork's water loading the
                    // two differ by ~3.4 %, and this integral is the precipitable-water
                    // diagnostic. Same repair as waterVapourEvaporation() above.
                    const double M_nw_pw     = AtmMixture::M_nonwater(m.c.x[i][j][k],
                                                   m.co2.x[i][j][k], m.m_comp.M_bg);
                    const double e          = HPA_TO_PA * SaturationH2O::vapourPressureFromMassFraction(
                                                  q_mixing, p_actual, M_nw_pw);
                    const double a          = e / (m.R_WaterVapour * t_actual);
                    const double step       = m.get_layer_height(i + 1) - m.get_layer_height(i);
                    const double local_mass = a * step;

                    m.PrecipitableWaterLocal.x[i][j][k] = local_mass;
                    column_sum += local_mass / m.r_0_water;
                }

                m.precipitable_water.y[j][k] = column_sum * 1e3;   // m -> mm
            }  // k
        }  // j

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" Time measured: %.3f seconds for PrecipitableWater\n", elapsed.count() * 1e-9);
        cout << "      PrecipitableWater ended" << endl;
    }

    // ------------------------------------------------------------------
    void vegetationLand()
    {
        using namespace std;
        cout << "\n\n\n      VegetationLand" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        constexpr double max_vegetation_height_m = 4400.0;
        constexpr double min_vegetation_temp_c   = -40.0;

        std::vector<double> height_table(m.im);
        for (int i = 0; i < m.im; i++)
            height_table[i] = m.get_layer_height(i);

        double max_Precipitation = 0.0;

        #pragma omp parallel for collapse(2) reduction(max:max_Precipitation)
        for (int j = 0; j < m.jm; j++)
            for (int k = 0; k < m.km; k++)
                max_Precipitation = std::max(max_Precipitation, m.Precipitation.x[0][j][k]);

        double inv_max_Precipitation = 1.0 / std::max(max_Precipitation, 1e-5);

        #pragma omp parallel for collapse(2)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                int i_mount = m.i_topography[j][k];
                if (is_land(m.h, 0, j, k)
                    && height_table[i_mount] < max_vegetation_height_m
                    && m.t.x[i_mount][j][k] * m.t_0 - m.t_0 > min_vegetation_temp_c) {
                    m.Vegetation.y[j][k] = m.Precipitation.x[0][j][k] * inv_max_Precipitation;
                } else {
                    m.Vegetation.y[j][k] = 0.0;
                }
            }
        }

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        cout << " time measured: " << fixed << setprecision(3)
             << elapsed.count() * 1e-9 << " seconds for VegetationLand" << endl;
        cout << "      VegetationLand ended" << endl;
    }

    // ------------------------------------------------------------------
    void printDataAtm()
    {
        using namespace std;
        cout << endl << endl << endl << "      print_data_atm" << endl;
        cout << endl << " properties of the atmosphere at the surface: " << endl;

        double temperature_average,
               precipitablewater_average, precipitation_average,
               co2_average, Evaporation_average, Evaporation_Dalton_average,
               Evaporation_Meyer_average, Evaporation_Rohwer_average;

        #pragma omp parallel sections
        {
            #pragma omp section
                { temperature_average         = (AtomUtils::GetMean_3D(m.jm, m.km, m.t) - 1.0) * m.t_0; }
            #pragma omp section
                { precipitablewater_average   = AtomUtils::GetMean_2D(m.jm, m.km, m.precipitable_water); }
            #pragma omp section
                { precipitation_average       = 365.0 * AtomUtils::GetMean_3D(m.jm, m.km, m.Precipitation); }
            #pragma omp section
                { co2_average                 = AtomUtils::GetMean_2D(m.jm, m.km, m.co2_total); }
            #pragma omp section
                { Evaporation_Dalton_average  = 365.0 * AtomUtils::GetMean_2D(m.jm, m.km, m.Evaporation_Dalton); }
            #pragma omp section
                { Evaporation_Meyer_average   = 365.0 * AtomUtils::GetMean_2D(m.jm, m.km, m.Evaporation_Meyer); }
            #pragma omp section
                { Evaporation_Rohwer_average  = 365.0 * AtomUtils::GetMean_2D(m.jm, m.km, m.Evaporation_Rohwer); }
            #pragma omp section
                { Evaporation_average         = 365.0 * AtomUtils::GetMean_2D(m.jm, m.km, m.Evaporation); }
        }
        // ATHAD: no "expected" mean temperature. That line looked the value up in the
        // Scotese global paleo-temperature curve, which is empty here (no curve reaches
        // 4.4 Ga) — and get_temperatures_from_curve() dereferences begin() and
        // decrements end() without checking emptiness, so the lookup was UB on an empty
        // map. Removed together with the curve machinery rather than guarded.

        cout.precision(2);
        cout << endl << endl;

        auto row = [](const char* n1, double v1, const char* u1,
                      const char* n2, double v2, const char* u2,
                      const char* n3, double v3, const char* u3) {
            cout << setiosflags(ios::left) << setw(40) << setfill('.')
                 << n1 << " = " << resetiosflags(ios::left) << setw(7)
                 << fixed << setfill(' ') << v1 << setw(6) << u1 << "   "
                 << setiosflags(ios::left) << setw(40) << setfill('.')
                 << n2 << " = " << resetiosflags(ios::left) << setw(7)
                 << fixed << setfill(' ') << v2 << setw(6) << u2 << "   "
                 << setiosflags(ios::left) << setw(40) << setfill('.')
                 << n3 << " = " << resetiosflags(ios::left) << setw(7)
                 << fixed << setfill(' ') << v3 << setw(6) << u3 << endl;
        };

        double precip_mm_a = precipitation_average * 8.64e4;

        row(" precipitable water average", precipitablewater_average, " mm",
            " precipitation average per year", precip_mm_a, " mm/a",
            " precipitation average per day", precip_mm_a / 365.0, " mm/d");

        // Per-component breakdown: which of P_rain / P_snow / P_graupel / P_conv
        // is climbing? Same units as the total: mean × 365·86400 gives mm/a.
        // Max scan reports the worst single cell with its (i,j,k) so a runaway
        // can be localised geographically.
        {
            const double s_per_year = 365.0 * 8.64e4;
            auto component = [&](const char* name, const Array& F){
                double mean_v = AtomUtils::GetMean_3D(m.jm, m.km, const_cast<Array&>(F));
                double mx = 0.0; int mi = 0, mj = 0, mk = 0;
                for(int i = 0; i < m.im; i++)
                    for(int j = 0; j < m.jm; j++)
                        for(int k = 0; k < m.km; k++){
                            double v = F.x[i][j][k];
                            if(v > mx){ mx = v; mi = i; mj = j; mk = k; }
                        }
                cout << "    " << setw(8) << setfill(' ') << name
                     << " mean = " << scientific << setprecision(3) << mean_v * s_per_year
                     << " mm/a   max = " << mx * s_per_year
                     << " mm/a  @(i=" << mi << ",j=" << mj
                     << ",k=" << mk << ")" << fixed << endl;
            };
            cout << " precipitation by component (mm/a, surface-equivalent):" << endl;
            component("Precip",    m.Precipitation);   // sanity check vs model total above
            component("P_rain",    m.P_rain);
            component("P_snow",    m.P_snow);
            component("P_graupel", m.P_graupel);
            component("P_conv",    m.P_conv);
        }

        row(" precipitable water average", precipitablewater_average, " mm",
            " Evaporation_Dalton_average per year", Evaporation_Dalton_average, " mm/a",
            " Evaporation_Dalton_average per day", Evaporation_Dalton_average / 365.0, " mm/d");

        row(" co2 column average", co2_average, " kg/m2",
            " Evaporation_average per year", Evaporation_Dalton_average, " mm/a",
            " Evaporation_average per day", Evaporation_Dalton_average / 365.0, " mm/d");

        row(" co2 column average", co2_average, " kg/m2",
            " Evaporation_Meyer_average per year", Evaporation_Meyer_average, " mm/a",
            " Evaporation_Meyer_average per day", Evaporation_Meyer_average / 365.0, " mm/d");

        row(" co2 column average", co2_average, " kg/m2",
            " Evaporation_Rohwer_average per year", Evaporation_Rohwer_average, " mm/a",
            " Evaporation_Rohwer_average per day", Evaporation_Rohwer_average / 365.0, " mm/d");

        cout << endl;

        cout << setiosflags(ios::left) << setw(40) << setfill('.')
             << " temperature_average" << " = " << resetiosflags(ios::left)
             << setw(7) << fixed << setfill(' ') << temperature_average
             << setw(6) << " deg" << "   "
             << setiosflags(ios::left) << setw(40) << setfill('.')
             << " temperature_average" << " = " << resetiosflags(ios::left)
             << setw(7) << fixed << setfill(' ') << temperature_average + m.t_0
             << setw(6) << " K" << endl << endl << endl;

        cout << "      print_data_atm ended" << endl;
    }

    // ------------------------------------------------------------------
    // ATHAD: CO2 is a well-mixed MASS FRACTION, uniform in the vertical.
    //
    // The Earth version built a ppm profile — a surface value scaled by local
    // temperature plus a paleo increment, decaying parabolically to a fixed tropopause
    // concentration, with separate vegetation / ocean / land ppm budgets and an Earth
    // regression in t_equat_modern for the "mean CO2 at modern times".
    //
    // None of that has a subject in the Hadean: there is no biosphere to draw CO2 down,
    // no carbonate-silicate ocean sink to absorb it, and no land. CO2 is simply 10 % of
    // the atmosphere by mole (20.5 % by mass) and stays where it is put. The vertical
    // gradient the parabola imposed encoded Earth's surface sources and stratospheric
    // depletion, so imposing it here would be inventing structure.
    //
    // Units: the field is now kg/kg, not ppm. At 20.5 % by mass ppm is meaningless, and
    // every consumer (the mixture properties, the radiative optical depth) wants a
    // fraction. co2_scale still multiplies the field for sensitivity experiments.
    void co2Atmosphere()
    {
        using namespace std;
        cout << endl << endl << endl << "      AGCM: co2_atmosphere" << endl;

        const double co2_ref = m.co2_0 * m.co2_scale;                   // [kg/kg] at the sea

        // ATHAD_COND: "well mixed" is a statement about the DRY air, not about the
        // mass fraction.
        //
        // ATHAD set co2 to one uniform mass fraction everywhere, and that was right there:
        // its water field is uniform too, so a uniform CO2 mass fraction and a uniform
        // CO2:background ratio are the same statement. Here the water mass fraction runs
        // from 0.346 at the sea to 0.004 above the cold trap, and the two statements come
        // apart. Holding the MASS fraction uniform would mean the non-water air changes
        // composition with height — the background would have to make up all 0.34 of the
        // mass the water vacates — and the column's gas constant then settles at 230 J/(kg K)
        // aloft instead of the 195 the dry composition actually has. An 18 % error in R
        // through the entire upper atmosphere, from a field that is not supposed to have
        // any structure at all.
        //
        // What is physically well mixed is the CO2:background MOLE ratio, and since that
        // ratio is fixed, so is their mass ratio within the dry air:
        //
        //     q_CO2(i) = (1 - c(i)) * f_CO2,     f_CO2 = q_CO2 / (q_CO2 + q_bg) at the sea
        //
        // This keeps CO2 a passive, source-free tracer — the conservation test in
        // co2Column() is unchanged — while making the dry mixture it belongs to uniform,
        // which is the thing that was meant. The uniformity diagnostic has to change with
        // it: min == max on co2 was ATHAD's test, and here the invariant is that
        // co2/(1-c) is constant instead.
        const double f_CO2 = (m.m_comp.q_CO2 + m.m_comp.q_bg > 0.0)
                           ? m.m_comp.q_CO2 / (m.m_comp.q_CO2 + m.m_comp.q_bg)
                           : co2_ref;

        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++)
            for (int k = 0; k < m.km; k++)
                for (int i = 0; i < m.im; i++)
                    m.co2.x[i][j][k] = (1.0 - m.c.x[i][j][k]) * f_CO2 * m.co2_scale;

        cout.precision(6);
        cout << "      AGCM: co2 well mixed in the DRY air at f_CO2 = " << f_CO2
             << " kg/kg of non-water (co2_0 = " << m.co2_0
             << " at the sea, co2_scale = " << m.co2_scale << ")" << endl;
        cout << "      AGCM: co2_atmosphere ended" << endl;
    }

    // ------------------------------------------------------------------
    // Global TOTAL-water conservation check.
    //
    // ATHAD has no water source and no water sink: the surface is supercritical so
    // waterVapourEvaporation() returns at once, and nothing else adds or removes water.
    // The mass-weighted mean of (vapour + cloud + ice + graupel) is therefore conserved
    // exactly, and any drift is scheme error or a limiter deleting mass. Nothing measured
    // it, which is why the sub-cloud band could sit at an impossible composition for as
    // long as it did.
    //
    // Reported alongside the mass deleted by the c <= 1 - co2 ceiling, so the two can be
    // compared: if the drift is accounted for by the clipping, the transport is conserving
    // and the problem is upstream of it; if it is not, something is creating water.
    void waterBudget(bool report, const char* stage = "")
    {
        using namespace std;

        double w_num = 0.0, w_den = 0.0;

        #pragma omp parallel for collapse(2) reduction(+:w_num,w_den) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                const double coslat = cos((j / (double)(m.jm - 1) - 0.5) * M_PI);

                for (int i = 0; i < m.im; i++) {
                    const double dp_Pa = (i < m.im - 1)
                        ? (m.p_stat.x[i][j][k] - m.p_stat.x[i+1][j][k]) * 100.0
                        :  m.p_stat.x[i][j][k] * 100.0;
                    if (!(dp_Pa > 0.0)) continue;

                    const double u_air = dp_Pa / m.g;                       // [kg/m2]
                    const double q_w   = std::max(0.0, m.c.x[i][j][k])
                                       + std::max(0.0, m.cloud.x[i][j][k])
                                       + std::max(0.0, m.ice.x[i][j][k])
                                       + std::max(0.0, m.gr.x[i][j][k]);

                    w_den += coslat * u_air;
                    w_num += coslat * q_w * u_air;
                }
            }
        }

        const double q_mean = (w_den > 0.0) ? w_num / w_den : 0.0;
        if (m.m_q_h2o_ref <= 0.0) m.m_q_h2o_ref = q_mean;

        // ------------------------------------------------------------------
        // The same mean against FROZEN weights, and the air mass those weights carry.
        //
        // q_mean above is water mass over air mass, and BOTH are read from p_stat, which
        // densities() re-integrates hydrostatically every iteration on the local R and cp —
        // which depend on the composition. So the weights move when the water moves, and a
        // drift in q_mean is not by itself evidence that any water was created: it can be
        // the column being re-weighed. Nothing separated the two, and three separate
        // repairs to the transport (the metric fix, the anelastic projection, the diffusive
        // flux term) each left the drift identical to three significant figures, which is
        // not how a transport error behaves.
        //
        // q_fixed applies the reference weights to the current field: the horizontal mean
        // of q at each level, weighted by the air mass that level had at the reference
        // time. If q_fixed is flat while q_mean drifts, the water did not move — the
        // weighing did. air_mean is the column air mass itself, the same test at one
        // remove: it is a fixed 250 bar column and must not drift at all.
        std::vector<double> qbar(m.im, 0.0), abar(m.im, 0.0);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < m.im; i++) {
            double qn = 0.0, qd = 0.0, an = 0.0;
            for (int j = 0; j < m.jm; j++) {
                const double coslat = cos((j / (double)(m.jm - 1) - 0.5) * M_PI);
                for (int k = 0; k < m.km; k++) {
                    const double dp_Pa = (i < m.im - 1)
                        ? (m.p_stat.x[i][j][k] - m.p_stat.x[i+1][j][k]) * 100.0
                        :  m.p_stat.x[i][j][k] * 100.0;
                    const double q_w = std::max(0.0, m.c.x[i][j][k])
                                     + std::max(0.0, m.cloud.x[i][j][k])
                                     + std::max(0.0, m.ice.x[i][j][k])
                                     + std::max(0.0, m.gr.x[i][j][k]);
                    qn += coslat * q_w;
                    qd += coslat;
                    if (dp_Pa > 0.0) an += coslat * dp_Pa / m.g;
                }
            }
            qbar[i] = (qd > 0.0) ? qn / qd : 0.0;                  // plain horizontal mean of q
            abar[i] = an;                                          // air mass carried by the level
        }
        if (m.m_air_levels.empty()) m.m_air_levels = abar;

        double qf_num = 0.0, qf_den = 0.0, air_now = 0.0, air_ref = 0.0;
        for (int i = 0; i < m.im; i++) {
            qf_num  += m.m_air_levels[i] * qbar[i];
            qf_den  += m.m_air_levels[i];
            air_now += abar[i];
            air_ref += m.m_air_levels[i];
        }
        const double q_fixed = (qf_den > 0.0) ? qf_num / qf_den : 0.0;
        if (m.m_q_h2o_fixed_ref <= 0.0) m.m_q_h2o_fixed_ref = q_fixed;

        // Per-level attribution: which levels the drift is actually appearing in. The
        // global number says water is created; only this says where, and "where" decides
        // whether the cause is the transport, a limiter, or the moist physics.
        std::vector<double> lev(m.im, 0.0);
        for (int i = 0; i < m.im; i++) {
            double num = 0.0;
            for (int j = 0; j < m.jm; j++) {
                const double coslat = cos((j / (double)(m.jm - 1) - 0.5) * M_PI);
                for (int k = 0; k < m.km; k++) {
                    const double dp_Pa = (i < m.im - 1)
                        ? (m.p_stat.x[i][j][k] - m.p_stat.x[i+1][j][k]) * 100.0
                        :  m.p_stat.x[i][j][k] * 100.0;
                    if (!(dp_Pa > 0.0)) continue;
                    const double q_w = std::max(0.0, m.c.x[i][j][k])
                                     + std::max(0.0, m.cloud.x[i][j][k])
                                     + std::max(0.0, m.ice.x[i][j][k])
                                     + std::max(0.0, m.gr.x[i][j][k]);
                    num += coslat * q_w * dp_Pa / m.g;
                }
            }
            lev[i] = (w_den > 0.0) ? num / w_den : 0.0;
        }
        if (m.m_q_h2o_levels.empty()) m.m_q_h2o_levels = lev;

        if (report) {
            const double drift = (m.m_q_h2o_ref > 0.0)
                               ? 100.0 * (q_mean / m.m_q_h2o_ref - 1.0) : 0.0;
            cout << "      AGCM: H2O" << (stage[0] ? " [" : "") << stage << (stage[0] ? "]" : "")
                 << " mass-weighted mean q = " << fixed << setprecision(6)
                 << q_mean << " kg/kg   (initial " << m.m_q_h2o_ref
                 << ", drift " << showpos << setprecision(4) << drift << " %"
                 << noshowpos << ";  deleted by the c ceiling so far "
                 << setprecision(6) << m.m_q_h2o_clipped << " kg/kg)" << endl;

            const double drift_fixed = (m.m_q_h2o_fixed_ref > 0.0)
                                     ? 100.0 * (q_fixed / m.m_q_h2o_fixed_ref - 1.0) : 0.0;
            const double drift_air   = (air_ref > 0.0)
                                     ? 100.0 * (air_now / air_ref - 1.0) : 0.0;
            cout << "            against FROZEN weights q = " << fixed << setprecision(6)
                 << q_fixed << " kg/kg (drift " << showpos << setprecision(4) << drift_fixed
                 << " %)" << noshowpos << ";   column air mass drift "
                 << showpos << setprecision(4) << drift_air << " %" << noshowpos << endl;

            // Rank the levels by how much of the drift they carry.
            std::vector<std::pair<double,int> > d;
            for (int i = 0; i < m.im; i++)
                d.push_back(std::make_pair(lev[i] - m.m_q_h2o_levels[i], i));
            std::sort(d.begin(), d.end(),
                      [](const std::pair<double,int>& a, const std::pair<double,int>& b){
                          return std::fabs(a.first) > std::fabs(b.first); });
            cout << "            drift by level (top 6):";
            for (int n = 0; n < 6 && n < (int)d.size(); n++)
                cout << "  i=" << d[n].second << " ("
                     << setprecision(0) << m.get_layer_height(d[n].second) * 1e-3 << " km) "
                     << showpos << scientific << setprecision(2) << d[n].first << noshowpos
                     << fixed;
            cout << endl;
        }
    }

    // ------------------------------------------------------------------
    // Column CO2 mass path [kg/m2] and the global CO2 mass-conservation check.
    //
    // This replaces the in-loop call to co2Atmosphere(). CO2 has a full transport
    // equation — RHS_Atm_Turb builds rhs_co2 and RungeKutta_Atm_Turb carries it through
    // all four stages — and co2Atmosphere() was overwriting the result with a uniform
    // field every iteration, so the CO2 in the output was exactly constant (min = max =
    // 0.205300 in every cell) and the "prognostic CO2" of CLAUDE.md was diagnostic.
    // co2Atmosphere() is now the INITIAL CONDITION only: the field starts well mixed, as
    // a non-condensable gas below the homopause should be, and is then transported.
    //
    // Because there is no CO2 source or sink anywhere in the model, the global
    // mass-weighted mean mass fraction is a conserved quantity, and any drift in it is
    // the transport scheme's error rather than physics. That is what the drift line
    // reports; watch it, since nothing else would notice.
    //
    // co2_total was declared ("areas of higher co2 concentration") and never filled, so
    // its min/max and the co2_average derived from it both read 0.000. It now carries the
    // column CO2 mass path, which is the CO2 analogue of precipitable water.
    void co2Column(bool report)
    {
        using namespace std;

        double w_num = 0.0, w_den = 0.0;

        #pragma omp parallel for collapse(2) reduction(+:w_num,w_den) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                const double coslat = cos((j / (double)(m.jm - 1) - 0.5) * M_PI);
                double col = 0.0;

                for (int i = 0; i < m.im; i++) {
                    // Layer air mass from the pressure drop across it; the top layer
                    // carries everything above it. Same construction as the optical depth.
                    const double dp_Pa = (i < m.im - 1)
                        ? (m.p_stat.x[i][j][k] - m.p_stat.x[i+1][j][k]) * 100.0
                        :  m.p_stat.x[i][j][k] * 100.0;
                    if (!(dp_Pa > 0.0)) continue;

                    const double u_air = dp_Pa / m.g;                       // [kg/m2]
                    const double q_c   = std::max(0.0, m.co2.x[i][j][k]);

                    col   += q_c * u_air;
                    w_den += coslat * u_air;
                    w_num += coslat * q_c * u_air;
                }
                m.co2_total.y[j][k] = col;                                  // [kg/m2]
            }
        }

        const double q_mean = (w_den > 0.0) ? w_num / w_den : 0.0;
        if (m.m_q_co2_ref <= 0.0) m.m_q_co2_ref = q_mean;                   // first call sets the reference

        if (report) {
            const double drift = (m.m_q_co2_ref > 0.0)
                               ? 100.0 * (q_mean / m.m_q_co2_ref - 1.0) : 0.0;
            cout << "      AGCM: CO2 mass-weighted mean q = " << fixed << setprecision(6)
                 << q_mean << " kg/kg   (initial " << m.m_q_co2_ref
                 << ", drift " << showpos << setprecision(4) << drift << " %"
                 << noshowpos << ")" << endl;
        }
    }

    // ------------------------------------------------------------------
    // ------------------------------------------------------------------------
    // Column profile diagnostic, ported from ATURAN 73406af.
    //
    // Prints height, temperature, pressure and density down one meridian column. This is
    // the single most useful check that a re-based atmosphere is hydrostatically sane —
    // it is what caught ATNEPT's pressure defect (ffee191). The top row answers the
    // question the domain sizing turns on: does the shell actually reach the radiating
    // level (~0.1 bar)?
    void printColumnProfile(int j_lat, const char* label)
    {
        using namespace std;
        const int k = 0;
        if (j_lat < 0 || j_lat >= m.jm) return;

        cout << endl << "      Column profile — " << label
             << "  (j = " << j_lat << ", k = " << k << ")" << endl;
        cout << "        "
             << setw(4)  << "i"      << setw(11) << "height[km]"
             << setw(10) << "T[K]"   << setw(13) << "p[bar]"
             << setw(11) << "rho"    << setw(10) << "R"
             << setw(9)  << "q_H2O"  << setw(10) << "q_sat"
             << setw(10) << "q_cld"  << setw(10) << "q_ice"
             << setw(14) << "phase" << endl;

        int i_cond_top = m.im;                 // lowest level where saturation can bite

        cout.precision(4);
        for (int i = 0; i < m.im; i++) {
            const double T   = m.t.x[i][j_lat][k] * m.t_0;
            const double p   = m.p_stat.x[i][j_lat][k];
            const double rho = m.r_humid.x[i][j_lat][k];
            const double R   = AtmMixture::R_of(m.c.x[i][j_lat][k],
                                                m.co2.x[i][j_lat][k], m.m_comp.R_bg);
            // Saturation state. q_sat = 1 means "no limit": either the cell is
            // supercritical (T >= 647.096 K, no liquid phase exists) or the vapour is
            // superheated (p_sat(T) exceeds the local pressure, so it cannot saturate).
            // Only where q_sat < q_H2O can water actually condense.
            const double M_other = AtmMixture::M_nonwater(m.c.x[i][j_lat][k],
                                                          m.co2.x[i][j_lat][k], m.m_comp.M_bg);
            const double q_sat   = SaturationH2O::saturationMassFractionAt(T, p, M_other);
            const double q_v     = m.c.x[i][j_lat][k];

            const double q_cld = m.cloud.x[i][j_lat][k];
            const double q_ice = m.ice.x[i][j_lat][k];

            const char* phase;
            if      (T >= AtmMixture::T_CRIT_H2O) phase = "supercrit";
            else if (q_sat >= 1.0)                phase = "superheat";
            else if (q_v > q_sat)                 phase = "CONDENSING";
            else                                  phase = "subsat";

            // Condensate where none can exist is the signature of a saturation formula
            // being evaluated outside its regime, which is how the dilute q_sat kept
            // manufacturing cloud here. Flag it in the phase column rather than leaving it
            // to be spotted in the density.
            if ((q_cld + q_ice) > 1.0e-6 && (T >= AtmMixture::T_CRIT_H2O || q_sat >= 1.0))
                phase = "CLOUD?!";

            if (q_sat < 1.0 && i < i_cond_top) i_cond_top = i;

            cout << "        " << setw(4) << i
                 << setw(11) << fixed << m.get_layer_height(i) * 1.0e-3
                 << setw(10) << setprecision(1) << T
                 << setw(13) << setprecision(5) << p * 1.0e-3
                 << setw(11) << setprecision(3) << rho
                 << setw(10) << setprecision(1) << R
                 << setw(9)  << setprecision(4) << q_v
                 << setw(10) << setprecision(4) << q_sat
                 << setw(10) << setprecision(4) << q_cld
                 << setw(10) << setprecision(4) << q_ice
                 << setw(14) << phase << endl;
        }

        // Where, if anywhere, water can condense in this column. For a runaway steam
        // atmosphere the honest answer may be "nowhere": below the critical temperature
        // the column must also be cool enough that p_sat(T) drops below the local
        // pressure, and a 250 bar column that is still ~700 K at its top never gets there.
        if (i_cond_top < m.im)
            cout << "        condensation possible from i = " << i_cond_top
                 << " (" << setprecision(1) << m.get_layer_height(i_cond_top) * 1.0e-3
                 << " km) upward" << endl;
        else
            cout << "        NO condensation anywhere in this column"
                 << " — supercritical or superheated at every level" << endl;

        // Outgoing longwave flux at the top of the atmosphere.
        //
        // THIS IS THE SHARPEST TEST OF THE RADIATION SCHEME. A runaway steam atmosphere
        // has a well-known property: its OLR saturates near the Nakajima /
        // Komabayashi-Ingersoll limit of ~280-310 W/m2 and becomes almost INDEPENDENT of
        // the surface temperature, because the emission level sits in the optically thick
        // water column rather than at the ground. If instead the OLR comes out near
        // sigma*T_surf^4 (287 kW/m2 at 1500 K) the atmosphere is radiatively transparent
        // and the opacity is wrong by orders of magnitude.
        //
        // Summed the same way as L_down but upward: each layer's emission attenuated by
        // every layer above it, plus the surface contribution attenuated by all of them.
        {
            double olr = 0.0, trans = 1.0, eps_sum = 0.0;
            for (int i = m.im - 1; i >= 0; i--) {
                const double eps = m.epsilon.x[i][j_lat][k];
                const double T_i = m.t.x[i][j_lat][k] * m.t_0;
                eps_sum += eps;
                olr   += eps * m.sigma * std::pow(T_i, 4.0) * trans;
                trans *= (1.0 - eps);
            }
            const double T_s = m.t.x[0][j_lat][k] * m.t_0;
            olr += m.sigma * std::pow(T_s, 4.0) * trans;          // surface, seen through the column

            const double sigT4_surf = m.sigma * std::pow(T_s, 4.0);

            if (eps_sum <= 0.0) {
                // Called before RadiationMultiLayer has ever run, so every layer emissivity
                // is still zero and the column is trivially transparent. Say so rather than
                // reporting sigma*T_surf^4 as if it were a computed OLR.
                cout << "        OLR: not yet meaningful — radiation has not run"
                     << " (all layer emissivities are zero)" << endl;
            } else {
                cout << "        OLR = " << setprecision(1) << olr << " W/m2"
                     << "   (sigma*T_surf^4 = " << sigT4_surf << " W/m2,"
                     << " suppression x" << setprecision(0) << sigT4_surf / std::max(olr, 1e-30)
                     << ", transmitted fraction " << setprecision(3) << trans << ")" << endl;
                if (olr > 1000.0)
                    cout << "        WARNING: OLR far above the ~280-310 W/m2 runaway limit"
                         << " — the column is too transparent" << endl;
                else if (olr >= 250.0 && olr <= 350.0)
                    cout << "        OLR is in the runaway-greenhouse range"
                         << " (Nakajima / Komabayashi-Ingersoll limit ~280-310 W/m2)" << endl;
            }
        }

        // The two facts the sizing depends on, stated rather than left to be read off.
        const double p_top = m.p_stat.x[m.im-1][j_lat][k] * 1.0e-3;
        cout << "        top of domain: " << setprecision(5) << p_top << " bar"
             << (p_top <= 0.1 ? "  (reaches the radiating level)"
                              : "  <-- ABOVE 0.1 bar: shell too shallow") << endl;
        if (p_top <= 0.0)
            cout << "        ERROR: non-positive pressure at the domain top" << endl;
        cout << endl;
    }

    // ------------------------------------------------------------------------
    // Per-level global summary: min / mean / max of temperature and pressure on every
    // radial level, plus an explicit non-positive-pressure count with its location.
    //
    // The column profile shows one meridian; this shows whether ANY column is
    // misbehaving, which is what a single profile cannot tell you.
    void printLevelSummary(const char* label)
    {
        using namespace std;
        cout << endl << "      Level summary — " << label << endl;
        cout << "        " << setw(4) << "i" << setw(11) << "height[km]"
             << setw(10) << "T_min" << setw(10) << "T_mean" << setw(10) << "T_max"
             << setw(13) << "p_min[bar]" << setw(13) << "p_mean[bar]"
             << setw(13) << "p_max[bar]" << setw(8) << "p<=0" << endl;

        int    bad_total = 0;
        int    bad_i = -1, bad_j = -1, bad_k = -1;
        double bad_p = 0.0;

        cout.precision(2);
        for (int i = 0; i < m.im; i++) {
            double t_min = 1e30, t_max = -1e30, t_sum = 0.0;
            double p_min = 1e30, p_max = -1e30, p_sum = 0.0;
            int    bad = 0;

            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    const double T = m.t.x[i][j][k] * m.t_0;
                    const double P = m.p_stat.x[i][j][k];
                    t_min = min(t_min, T); t_max = max(t_max, T); t_sum += T;
                    p_min = min(p_min, P); p_max = max(p_max, P); p_sum += P;
                    if (P <= 0.0) {
                        bad++;
                        if (P < bad_p) { bad_p = P; bad_i = i; bad_j = j; bad_k = k; }
                    }
                }
            }
            const double n = (double)(m.jm * m.km);
            bad_total += bad;

            cout << "        " << setw(4) << i
                 << setw(11) << fixed << setprecision(2) << m.get_layer_height(i) * 1.0e-3
                 << setw(10) << setprecision(1) << t_min
                 << setw(10) << t_sum / n
                 << setw(10) << t_max
                 << setw(13) << setprecision(5) << p_min * 1.0e-3
                 << setw(13) << p_sum / n * 1.0e-3
                 << setw(13) << p_max * 1.0e-3
                 << setw(8)  << bad << endl;
        }

        if (bad_total > 0) {
            cout << "        WARNING: " << bad_total
                 << " cells with non-positive static pressure; worst " << bad_p
                 << " hPa at (i=" << bad_i << ", j=" << bad_j << ", k=" << bad_k << ")" << endl;
        } else {
            cout << "        static pressure positive everywhere" << endl;
        }
        cout << endl;
    }

    // ------------------------------------------------------------------------
    // Planetary energy balance: the global means the column profile cannot show.
    //
    // The column diagnostic prints ONE meridian's OLR; the budget that has to close is
    // the global one, and the albedo that closes it is the model's OWN — built by
    // MultiLayerRadiation from the condensate the model actually made, not the clear-sky
    // value assumed in the config. Both are computed here, area-weighted by cos(latitude)
    // exactly as AtomUtils::GetMean_2D weights them.
    //
    // The implied skin temperature printed at the end is the fixed-point target for
    // t_skin: sigma*T_skin^4 = absorbed SW + geothermal, with THIS albedo. See
    // cAtmosphereModel::solveSkinTemperature().
    void printPlanetaryBalance(const char* label)
    {
        using namespace std;

        double alb_mean = 0.0, sw_mean = 0.0, abs_mean = 0.0;
        const bool measured = m.planetaryShortWave(alb_mean, sw_mean, abs_mean);

        double w_sum = 0.0, olr_w = 0.0, lid_w = 0.0, eps_top_w = 0.0;
        for (int j = 0; j < m.jm; j++) {
            const double w = cos((j / (double)(m.jm - 1) - 0.5) * M_PI);

            double olr_k = 0.0, lid_k = 0.0, eps_top_k = 0.0;
            for (int k = 0; k < m.km; k++) {
                // Upward flux at the top: every layer's emission attenuated by all the
                // layers above it, plus the surface seen through the whole column.
                double olr = 0.0, trans = 1.0;
                for (int i = m.im - 1; i >= 0; i--) {
                    const double eps = m.epsilon.x[i][j][k];
                    olr   += eps * m.sigma * pow(m.t.x[i][j][k] * m.t_0, 4.0) * trans;
                    trans *= (1.0 - eps);
                }
                olr += m.sigma * pow(m.t.x[0][j][k] * m.t_0, 4.0) * trans;
                olr_k     += olr;
                lid_k     += m.t.x[m.im-1][j][k] * m.t_0;
                eps_top_k += m.epsilon.x[m.im-1][j][k];
            }
            olr_k     /= (double)m.km;
            lid_k     /= (double)m.km;
            eps_top_k /= (double)m.km;

            w_sum     += w;
            olr_w     += w * olr_k;
            lid_w     += w * lid_k;
            eps_top_w += w * eps_top_k;
        }

        const double olr_mean = olr_w / w_sum;
        const double lid_mean = lid_w / w_sum;
        const double eps_top  = eps_top_w / w_sum;
        const double in_mean  = abs_mean + m.geothermal_flux;

        cout << endl << "      Planetary balance — " << label << endl;
        cout.precision(2);
        if (!measured) {
            cout << "        not yet meaningful — radiation has not run"
                 << " (albedo and layer emissivities are still zero)" << endl << endl;
            return;
        }
        cout << "        mean planetary albedo ............ = " << fixed << setprecision(4)
             << alb_mean << "   (clear-sky surface " << m.albedo_surface
             << ", thick cloud " << m.albedo_cloud << ")" << endl;
        cout << "        mean incident short wave ......... = " << setprecision(2)
             << sw_mean << " W/m2" << endl;
        cout << "        absorbed SW + geothermal ......... = " << abs_mean << " + "
             << m.geothermal_flux << " = " << in_mean << " W/m2" << endl;
        cout << "        outgoing long wave (OLR) ......... = " << olr_mean << " W/m2" << endl;
        cout << "        imbalance (in - out) ............. = " << (in_mean - olr_mean)
             << " W/m2" << endl;
        cout << "        implied skin temperature ......... = "
             << pow(in_mean / m.sigma, 0.25) << " K   (t_skin = " << m.t_skin << " K)" << endl;
        // Where the OLR actually comes from. The lid emissivity decides whether the OLR
        // is an integral over the column or simply sigma*T_lid^4 — and while the lid is
        // opaque it is the latter, so the OLR is a PRESCRIBED temperature reported back,
        // not a computed flux. Print the identity rather than let two
        // independently-derived-looking numbers read as an agreement.
        // See cAtmosphereModel::updateSkinTemperature and param.py (t_skin, L_atm).
        cout << "        lid temperature / emissivity ..... = " << setprecision(2) << lid_mean
             << " K  /  eps = " << setprecision(4) << eps_top << endl;
        cout << "        sigma*T_lid^4 .................... = " << setprecision(2)
             << m.sigma * pow(lid_mean, 4.0) << " W/m2";
        if (eps_top > 0.9)
            cout << "   <-- the OLR above IS this: the lid is opaque, so the model"
                 << " reports a prescribed temperature, not a computed flux";
        cout << endl;
        cout << endl;
    }

    // ------------------------------------------------------------------------
    // set_water_profile: diagnose c in the same sweep as T and p. INITIALISATION ONLY —
    // during the time loop c is prognostic and this must stay false, or the transport is
    // overwritten by the profile every iteration.
    void densities(bool set_water_profile = false)
    {
        using namespace std;
        cout << "\n\n\n      PressureDensity" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        // ATHAD: the gas constant is LOCAL, not a constant.
        //
        // The Earth version used R_Air everywhere because water vapour was a ~1 % trace on a
        // fixed N2/O2 carrier, so R varied by well under a percent. Here H2O is 67 % of the
        // mass, so R of a parcel is set by its own composition and moves as the water field
        // does: pure background is 317 J/(kg K), pure steam 462, the reference mixture 388.
        //
        // Two distinct gas constants appear below and must not be confused:
        //   R_mix  — the COLUMN reference, which sets the surface pressure and scale height.
        //            p_sl uses it because r_air was calibrated as p/(R_mix*T); using R_Air
        //            here instead yields 204 bar rather than the intended 250.
        //   R_of() — the LOCAL value at each cell, used for the densities.
        const double R_mix            = m.m_comp.R_mix;
        const double p_sl_factor      = 1e-2 * m.r_air * R_mix;
        const double M_bg             = m.m_comp.M_bg;
        const double R_bg             = m.m_comp.R_bg;

        // ATHAD: the column is integrated, not fitted.
        //
        // The inherited profile was the COSMO barometric form
        // T(h) = T0*sqrt(1 - 2*beta*g*h/(R*T0^2)), an empirical shape for Earth's
        // troposphere with beta tuned to reproduce ~5 K/km near the ground. Two things
        // break it here. It is a SQRT in height, so matching its near-surface slope to the
        // dry adiabat does not make it an adiabat — it plunges to zero at
        // h = R*T0^2/(2*beta*g), which for the dry-adiabatic beta is 156 km, well inside a
        // 300 km domain. And beta is a fitted Earth constant with no meaning for this
        // mixture.
        //
        // Replaced by the physics the profile is meant to represent, integrated layer by
        // layer:
        //     dry adiabat     dT/dz = -g/cp        (cp local: follows composition and T)
        //     hydrostatic     dp/dz = -p*g/(R*T)   (R local, on the layer-mean T)
        //     isothermal top  T = t_skin           where the adiabat falls below it
        //
        // Exact for a constant-cp adiabat, correct to O(dz^2) as cp and R vary, needs no
        // tuned constant, and cannot produce the zero temperature the sqrt form did.

        std::vector<double> height_table(m.im);
        for (int i = 0; i < m.im; i++)
            height_table[i] = m.get_layer_height(i);

        #pragma omp parallel for collapse(2) schedule(dynamic, 4)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                // Surface anchor. The floor guards against a transient cold surface making
                // the integration meaningless; 180 K is far below anything physical here.
                double T_prev = std::max(180.0, m.t.x[0][j][k] * m.t_0);
                double p_prev = p_sl_factor * T_prev;                   // [hPa]
                double q_sat_min = 1.0;                                 // cold trap, set_water only

                for (int i = 0; i < m.im; i++) {
                    const double q_v = m.c.x[i][j][k];
                    const double q_c = m.co2.x[i][j][k];
                    const double R_loc = AtmMixture::R_of(q_v, q_c, R_bg);

                    double T_i, p_i;
                    if (i == 0) {
                        T_i = T_prev;
                        p_i = p_prev;
                    } else {
                        const double dz = height_table[i] - height_table[i-1];

                        // ATHAD_COND: MOIST adiabat where the layer is saturated, dry where
                        // it is not.
                        //
                        // ATHAD integrates dT/dz = -g/cp throughout, and that is correct
                        // there for a stated reason: water is supercritical from the ground
                        // to ~177 km, so nothing condenses and there is no latent heat to
                        // release. Here the air over the sea IS saturated and stays so up to
                        // the cold trap, which is a third of the mass of the column. Using
                        // the dry adiabat through it overstates the lapse by 44 % and puts
                        // the cold trap kilometres too low.
                        //
                        // The saturation test uses the water actually present against the
                        // local q_sat, so the switch is a property of the state rather than
                        // a prescribed height — the same choice as the cold trap in
                        // initWaterWapour. moistLapse returns g/cp itself wherever no
                        // condensation is possible, so the branch below is about which
                        // physics APPLIES, not about avoiding a bad value.
                        // The step from i-1 to i is a property of the parcel AT i-1, so the
                        // composition, the heat capacity and the saturation test all come
                        // from level i-1 — the level T_prev and p_prev belong to. Using
                        // level i's water against level i-1's saturation is an index
                        // mismatch that looks harmless and is not: q_sat falls with height,
                        // so c(i) is always below q_sat(i-1) by exactly one level's worth,
                        // the test fails at EVERY level, and the moist branch never runs.
                        // It cost an hour here, silently, with a column that looked right.
                        const double q_v_p  = m.c.x[i-1][j][k];
                        const double q_c_p  = m.co2.x[i-1][j][k];
                        const double cp_loc = AtmMixture::cp_of(q_v_p, q_c_p, T_prev, M_bg);
                        const double M_nw   = AtmMixture::M_nonwater(q_v_p, q_c_p, M_bg);
                        const double q_s    = SaturationH2O::saturationMassFractionAt(
                                                  T_prev, p_prev, M_nw);

                        // 0.99 rather than 1.0: the profile is built from q_sat in the first
                        // place, so exact equality is a rounding coin-flip and the lapse rate
                        // would flicker between the two branches from level to level.
                        const double gamma = (q_v_p >= 0.99 * q_s)
                            ? SaturationH2O::moistLapse(T_prev, p_prev, M_nw,
                                                        cp_loc, R_loc, m.g)
                            : m.g / cp_loc;

                        const double T_ad = T_prev - gamma * dz;

                        // Isothermal once the adiabat drops below the radiative skin value.
                        T_i = std::max(m.t_skin, T_ad);

                        // Hydrostatic, integrated on the layer-mean temperature.
                        const double T_mean = 0.5 * (T_prev + T_i);
                        p_i = p_prev * exp(-m.g * dz / (R_loc * T_mean));
                    }

                    m.t.x[i][j][k]      = T_i / m.t_0;
                    m.p_stat.x[i][j][k] = p_i;

                    // ATHAD_COND: on the initial call the water profile is DIAGNOSED in this
                    // same sweep instead of being read from c.
                    //
                    // A saturated column is one ODE — dT/dz from the moist lapse, dp/dz
                    // hydrostatic, q = q_sat(T, p) — and its three fields have to be
                    // integrated together. Building c from a profile and then rebuilding the
                    // profile from that c is a fixed-point iteration over a POSITIVE feedback
                    // (a warmer column holds more water, which releases more latent heat,
                    // which warms it further), and it does not settle: measured here, 3 passes
                    // gave 446.7 K at 8.9 km, 8 passes gave 463.8, and it was still climbing.
                    // That is the runaway-greenhouse feedback arriving through the numerics
                    // instead of the physics, and no number of passes fixes it.
                    //
                    // One consistent sweep has no such problem. q at level i is set from the
                    // T and p just computed FOR level i, and the lapse rate that carries the
                    // sweep from i to i+1 reads it back at i — so the water profile, the
                    // temperature and the pressure are all the same integration.
                    if (set_water_profile && i > 0) {
                        const double M_nw_i = AtmMixture::M_nonwater(q_v, q_c, M_bg);
                        const double q_s_i  = SaturationH2O::saturationMassFractionAt(
                                                  T_i, p_i, M_nw_i);
                        if (q_s_i < q_sat_min) q_sat_min = q_s_i;
                        m.c.x[i][j][k] = std::max(q_sat_min, m.c_h2o_dry_top);
                    }

                    const double water_factor = std::max(0.5, 1.0
                                        - m.cloud.x[i][j][k] - m.ice.x[i][j][k]);
                    m.r_humid.x[i][j][k] = 1e2 * p_i / (R_loc * T_i * water_factor);

                    const double R_dry_loc = AtmMixture::R_of(0.0, q_c, R_bg);
                    m.r_dry.x[i][j][k]     = 1e2 * p_i / (R_dry_loc * T_i);

                    T_prev = T_i;
                    p_prev = p_i;
                }

                // Refresh the lid snapshot bcRadius pins t to. It used to be taken once
                // from the initial condition, which was right while t_skin was a constant
                // and wrong the moment it became a fixed-point iterate: densities() would
                // rebuild the column on the new t_skin and bcRadius would then pin the lid
                // back to the ORIGINAL value, so the top — which is the emission level, and
                // therefore the whole OLR — never moved. The pin still holds the lid against
                // the cubic extrapolation it replaced; it now holds it to the current
                // prescribed profile rather than to a stale one.
                if ((int)m.t_top_init.size() == m.jm)
                    m.t_top_init[j][k] = m.t.x[m.im-1][j][k];

                const int    i_m = m.i_topography[j][k];
                m.p_stat_landscape.y[j][k] = m.p_stat.x[i_m][j][k];

                if (i_m >= 0 && i_m < m.im) {
                    m.p_stat.x[0][j][k]  = m.p_stat.x[i_m][j][k];
                    m.r_dry.x[0][j][k]   = m.r_dry.x[i_m][j][k];
                    m.r_humid.x[0][j][k] = m.r_humid.x[i_m][j][k];
                }
            }
        }

        // ------------------------------------------------------------------
        // Anelastic base state rho_bar(z) — step 2 of the anelastic scope.
        //
        // Rebuilt here because this is where the profile is rebuilt: the base state the
        // pressure projection uses must be the same density field the rest of the model
        // sees, or the projection removes the divergence of a column that does not exist.
        // cos(latitude)-weighted, matching waterBudget() and planetaryShortWave().
        //
        // Only the horizontal MEAN enters. r_humid varies horizontally by a few per cent
        // (the surface temperature is a shallow parabola and there is no topography),
        // against five orders of magnitude vertically, so a one-dimensional base state
        // loses almost nothing and keeps the elliptic operator constant in time.
        m.m_rho_base.assign(m.im, 0.0);
        m.m_dlnrho_dr.assign(m.im, 0.0);
        {
            std::vector<double> lnrho(m.im, 0.0);
            for (int i = 0; i < m.im; i++) {
                double num = 0.0, den = 0.0;
                for (int j = 0; j < m.jm; j++) {
                    const double coslat = cos((j / (double)(m.jm - 1) - 0.5) * M_PI);
                    for (int k = 0; k < m.km; k++) {
                        const double rho = m.r_humid.x[i][j][k];
                        if (!(rho > 0.0)) continue;
                        num += coslat * rho;
                        den += coslat;
                    }
                }
                m.m_rho_base[i] = (den > 0.0) ? num / den : 0.0;
                lnrho[i] = (m.m_rho_base[i] > 0.0) ? log(m.m_rho_base[i]) : 0.0;
            }

            // d ln(rho_bar)/d(rad.z). rad.z is uniform (r0 + i*dr, dr = 1/(im-1)), so the
            // centred difference is the same stencil the solver uses on every other field;
            // the ends get the one-sided form so the first and last interior cells, which
            // are exactly where the projection's wall BC acts, are not fed a half-step.
            const double inv_2dr = 1.0 / (2.0 * m.dr);
            const double inv_dr  = 1.0 / m.dr;
            for (int i = 1; i < m.im - 1; i++)
                m.m_dlnrho_dr[i] = (lnrho[i+1] - lnrho[i-1]) * inv_2dr;
            if (m.im > 1) {
                m.m_dlnrho_dr[0]        = (lnrho[1] - lnrho[0]) * inv_dr;
                m.m_dlnrho_dr[m.im - 1] = (lnrho[m.im-1] - lnrho[m.im-2]) * inv_dr;
            }
        }

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for PressureDensity\n", elapsed.count() * 1e-9);
        cout << "      PressureDensity ended" << endl;
    }

private:
    cAtmosphereModel& m;
};

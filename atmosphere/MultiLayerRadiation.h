#pragma once

#include "MixtureAtm.h"
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

// ============================================================================
// MultiLayerRadiation — friend class of cAtmosphereModel
//
// Multi-layer long/short-wave radiation model for the surface-temperature
// distribution (computation of the local temperature from short- and long-wave
// radiation). For each (j,k) column it builds:
//   - the grey-body emission of every layer (Stefan-Boltzmann, sigma T^4),
//   - the layer emissivity epsilon from the water-vapour law (Bignami 1995,
//     valid -40..45 C) plus an optional CO2 term (Atwater & Ball) — the CO2
//     contribution is presently disabled (eps_co2 = 0; the CO2 coupling is
//     rough and to be updated),
// then solves the layer radiative-balance as a tridiagonal system with the
// Thomas algorithm and inverts the resulting radiation back to temperature.
// The surface layer receives the albedo-reduced incoming short-wave flux.
//
// Refactored from the former free member function
// cAtmosphereModel::RadiationMultiLayer() into the project's friend-class idiom
// (mirrors ThermoAtm / TurbulenceAtm / PressureSolverAtm). The (j,k) columns are
// independent, so the outer j loop is OpenMP-parallel with ALL Thomas-solver
// scratch (step/alfa/beta/AA/CA/CC and the per-column radiation_original)
// thread-local — the original reused a single shared set of scratch vectors,
// which would race under OpenMP.
//
// Reads : t, c, co2, p_stat, h, layer heights, albedo/short-wave/emissivity
//         constants (albedo_pole/equator, rad_pole/equator_short, sigma, ep,
//         co2_0, t_0, p_0).
// Writes: t (updated), radiation, epsilon, epsilon_2D, albedo, short_wave_radiation.
//
// Usage:  MultiLayerRadiation(*this).run();
// ============================================================================


class MultiLayerRadiation {
public:
    explicit MultiLayerRadiation(cAtmosphereModel& model)
        : m(model)
    {}

    void run()
    {
        using namespace std;
        cout << endl << "      RadiationMultiLayer" << endl;
        auto begin = std::chrono::high_resolution_clock::now();

        const int j_max  = m.jm - 1;
        const int j_half = j_max / 2;

        // ---- latitude profiles (computed once; read-only in the parallel loop) ----

        // Surface albedo with a temperature-dependent ICE/SNOW FEEDBACK (replaces the former
        // fixed pole->equator parabola, which froze the feedback out). Base value by surface
        // type (ocean vs land); where the surface is cold, blend toward a bright snow/sea-ice
        // albedo over a smooth ramp around freezing. Recomputed every MLR call on the CURRENT
        // surface T, so it is a LIVE feedback: warming -> less ice -> lower albedo -> more
        // absorbed shortwave -> extra (polar-amplified) warming that the dynamics cannot mix
        // away (unlike a forcing/nudge). Constants are tunable. NOTE the cloud SW bump below
        // overwrites this per column where cloud is present, so its net reach is cloud-limited.
        // ATHAD: the surface is a magma ocean, and there is no ice-albedo feedback.
        //
        // What was here: an ocean/land base albedo with a ramp toward snow/sea-ice as the
        // surface cooled through 275 -> 265 K. None of it applies. There is no land, no
        // snow, and at 1500 K nothing within 1200 K of the ice thresholds — the ramp was
        // dead code that always returned the ice-free ocean value.
        //
        // A quenching silicate melt is dark: measured basaltic-melt albedos are ~0.05-0.10.
        // This is the CLEAR-SKY value only; the cloud bump below raises it wherever the
        // model actually produces condensate. That separation is the point of this fix —
        // the reflective cloud deck a runaway greenhouse is supposed to have must be EARNED
        // by condensate the model generated, not asserted as a constant albedo. Asserting
        // 0.4 while the column condenses nothing was the inconsistency being removed.
        //
        // Now a PARAMETER (albedo_surface), not a literal. It was written here as a
        // constexpr while the config carried an inert albedo_pole/albedo_equator pair that
        // nothing read — so the file said one thing and the configuration said another.
        const double alb_surface_molten = m.albedo_surface;   // dark silicate melt, clear sky
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < m.jm; j++)
            for (int k = 0; k < m.km; k++)
                m.albedo.y[j][k] = alb_surface_molten;

        // Incoming short-wave radiation: pole -> equator parabola, hemispherically symmetric.
        m.short_wave_radiation = std::vector<double>(m.jm, m.rad_pole_short);
        const double rad_short_eff = m.rad_pole_short - m.rad_equator_short;
        for (int j = j_half; j >= 0; j--)
            m.short_wave_radiation[j] =
                rad_short_eff * parabola((double)j / (double)j_half) + m.rad_pole_short;
        for (int j = j_max; j > j_half; j--)
            m.short_wave_radiation[j] = m.short_wave_radiation[j_max - j];

        // ---- per-column radiative balance (columns independent -> OpenMP over j) ----
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < m.jm; j++) {

            // Thread-local sweep scratch; reused across k within this j. Every column
            // fully overwrites the entries it later reads, so reuse is race-free.
            //   up[i]  upward   long-wave flux leaving the TOP    of layer i  [W/m2]
            //   dn[i]  downward long-wave flux leaving the BOTTOM of layer i  [W/m2]
            //   T[i]   layer temperature [K], updated in place by the sweeps
            std::vector<double> up(m.im, 0.0), dn(m.im + 1, 0.0), T(m.im, 0.0);

            const int i_trop  = m.im - 1;   // top layer
            const int i_mount = 0;          // surface / bottom layer

            for (int k = 0; k < m.km; k++) {

                for (int i = i_mount; i <= i_trop; i++)
                    T[i] = m.t.x[i][j][k] * m.t_0;

                // ATHAD: layer optical depth from COLUMN MASS with pressure broadening.
                //
                // What was here: the Bignami (1995) clear-sky column emissivity
                // eps = 0.684 + 0.0056*e_surf, split into a "dry baseline" and a
                // water-vapour part and distributed over the layers, plus an Atwater & Ball
                // CO2 band with a scale factor tuned to 0.17 to reproduce Earth's ~30 W/m2
                // of CO2 greenhouse. Every one of those numbers is a regression on
                // present-day terrestrial columns — Bignami is a fit to MEDITERRANEAN SEA
                // SURFACE measurements. At 250 bar with a 67 %-by-mass water column they are
                // extrapolated some four orders of magnitude beyond their calibration, and
                // eps = 0.684 + 0.0056*e_surf saturates to 0.999 the instant e_surf exceeds
                // ~56 hPa (here it is ~2e5 hPa), so the whole scheme degenerates to "every
                // layer is a blackbody" and carries no information about the composition.
                //
                // What replaces it: a grey optical depth built from the absorber mass each
                // layer actually contains,
                //
                //     tau_i = SUM_s kappa_s * u_s,i * (p_i / p_ref)
                //     u_s,i = q_s * dp_i / g            [kg/m2]  column mass of species s
                //
                // The (p_i/p_ref) factor is pressure broadening: collisional line widths grow
                // in proportion to pressure, so the absorption per unit mass does too. It is
                // the term that matters most here and the one no Earth-calibrated emissivity
                // fit contains — at 250 bar it is a factor of 250 over the 1 bar reference,
                // and it is what makes the deep atmosphere opaque and the thin top
                // transparent, rather than a fit saturating everywhere.
                //
                // kappa values are grey-band mass absorption coefficients [m2/kg]. The water
                // value is the one conventionally used for grey runaway-greenhouse models;
                // CO2 is weaker per unit mass; the N2/CO/CH4 background is nearly
                // transparent in the thermal infrared. These are ASSUMPTIONS with roughly a
                // factor-of-two uncertainty and are the single biggest lever on the answer —
                // see README, and the OLR check below, which is what tests them.
                const double kappa_H2O = m.kappa_H2O;                 // [m2/kg]
                const double kappa_CO2 = m.kappa_CO2;                 // [m2/kg]
                // COMPOSITION-WEIGHTED, not the lumped constant. Equal to kappa_N2 here,
                // because this fork's background IS nitrogen; the sibling's is not, and the
                // same lumped 1e-6 was 1867x too small there. See param.py.
                const double kappa_bg  = m.kappaBackground();         // [m2/kg]
                constexpr double p_ref = 1.0e5;                       // [Pa] 1 bar broadening reference
                const double inv_g     = 1.0 / m.g;

                double lwp_col = 0.0, iwp_col = 0.0;                  // condensate paths [g/m2] (SW albedo bump)
                for (int i = i_mount; i <= i_trop; i++) {
                    // Layer mass, from the pressure drop across it. dp in hPa -> Pa.
                    double dp_Pa = (i < i_trop)
                                 ? (m.p_stat.x[i][j][k] - m.p_stat.x[i+1][j][k]) * 100.0
                                 :  m.p_stat.x[i][j][k] * 100.0;      // top layer carries all mass above
                    if (dp_Pa < 0.0) dp_Pa = 0.0;

                    const double q_v = std::max(0.0, m.c.x[i][j][k]);
                    // LOCAL CO2 mass fraction, not the field read raw (items 57/59). Here the
                    // difference is large: pinning q_c understated the CO2 optical depth aloft
                    // by ~50 %, because kappa_CO2 = 1e-3 against kappa_bg = 1e-6 and the old
                    // split was handing 0.31 of the mass to the background at the top.
                    const double q_c = std::max(0.0, AtmMixture::q_CO2_of(m.c.x[i][j][k],
                                                                          m.co2.x[i][j][k]));
                    const double q_b = std::max(0.0, 1.0 - q_v - q_c);

                    const double u_col = dp_Pa * inv_g;               // [kg/m2] total layer mass
                    const double broad = m.p_stat.x[i][j][k] * 100.0 / p_ref;   // pressure broadening

                    double tau_gas = (kappa_H2O * q_v + kappa_CO2 * q_c + kappa_bg * q_b)
                                   * u_col * broad;

                    // Cloud liquid + ice longwave greenhouse, unchanged in form (Stephens
                    // 1978 mass absorption), but the density now comes from the LOCAL mixture
                    // gas constant instead of the hard-coded 287.0 J/(kg K) of dry Earth air —
                    // which is 26 % off here and was applied to every layer.
                    constexpr double k_liq = 0.12, k_ice = 0.055;     // LW mass absorption [m2/g]
                    const double dz  = (i < i_trop) ? (m.get_layer_height(i+1) - m.get_layer_height(i))
                                                    : (m.get_layer_height(i) - m.get_layer_height(i-1));
                    const double T_i   = m.t.x[i][j][k] * m.t_0;
                    const double R_i   = AtmMixture::R_of(q_v, q_c, m.m_comp.R_bg);
                    const double rho_i = (T_i > 0.0) ? (m.p_stat.x[i][j][k] * 100.0) / (R_i * T_i) : 0.0;
                    const double cw_l  = std::max(0.0, m.cloud.x[i][j][k]);
                    const double cw_i  = std::max(0.0, m.ice.x[i][j][k]);
                    const double LWP_i = cw_l * rho_i * dz * 1000.0;  // [g/m2]
                    const double IWP_i = cw_i * rho_i * dz * 1000.0;  // [g/m2]
                    lwp_col += LWP_i;  iwp_col += IWP_i;

                    const double tau = tau_gas + k_liq * LWP_i + k_ice * IWP_i;

                    // exp(-tau) underflows to 0 for the huge optical depths the deep column
                    // carries, which is the correct answer (eps = 1, a perfect blackbody
                    // layer) — but guard it so no NaN can come out of the exponential.
                    m.epsilon.x[i][j][k] = (tau > 700.0) ? 1.0 : (1.0 - exp(-tau));

                    // Stash the LAYER optical depth; converted to the cumulative-from-the-top
                    // value in the downward pass below. Cannot be recovered from epsilon later
                    // because epsilon saturates at tau ~ 37 — see cAtmosphereModel.h.
                    m.tau_above.x[i][j][k] = tau;
                    m.tau_layer.x[i][j][k] = tau;      // kept: the resolution diagnostic
                }

                // tau_above: walk DOWNWARD from the lid, accumulating. After this,
                // tau_above[i] is the optical depth of everything ABOVE level i, so it is 0 at
                // the lid and monotonically increasing downward. One thread owns this whole
                // column (the parallel for is over j, k is inner), so this is race-free.
                {
                    double acc = 0.0;
                    for (int i = i_trop; i >= i_mount; i--) {
                        const double layer = m.tau_above.x[i][j][k];
                        m.tau_above.x[i][j][k] = acc;
                        acc += layer;
                    }
                }
                m.epsilon_2D.y[j][k] = m.epsilon.x[i_mount][j][k];

                // Cloud/ice SHORTWAVE albedo bump (stage 2): reflective clouds raise the column
                // albedo toward a cloud value, cutting absorbed SW (cooling) to compete with the LW
                // greenhouse above. alpha_eff = alpha_surf + (alpha_cloud - alpha_surf)*
                // (1 - exp(-k_sw*CWP_sw)), CWP_sw = LWP + f_ice*IWP (ice is optically thinner / less
                // reflective per unit path, so weighted down -> thin cirrus stays net-warming while
                // thick low liquid cloud net-cools). Overwrites the clear-sky latitude albedo for
                // this column; it is read below by the SW terms (tridiagonal source dd + surface
                // energy balance) and re-derived fresh each MLR call. See project_multilayer_radiation.
                {
                    // Two-stream-like cloud SW reflectivity on a CAPPED condensate path. The old
                    // exp(-k_sw*CWP) with k_sw=0.030 saturated at CWP~100 g/m2, so — and especially
                    // with the model's excessive cloud water (LWP ~1500 g/m2 vs observed ~100) —
                    // EVERY cloudy column was pinned at the asymptotic 0.60, wiping out the latitude
                    // gradient and MASKING the surface ice-albedo feedback (poles read the same 0.60
                    // as tropical ocean). Instead cap the path the radiation sees at a physical
                    // thick-cloud value and let reflectivity rise GENTLY as refl = tau/(tau+2)
                    // (0.5 at tau=2), composited over the surface (ice-feedback) albedo. Cloudy
                    // tropics now land ~0.3, thin/clear cells relax toward the surface value, and
                    // the polar ice albedo shows through — a physical gradient. (The same cloud-water
                    // excess still inflates the LW tau_cloud above; capping that is the next step.)
                    // alpha_cloud is a PARAMETER (albedo_cloud), and on ATHAD it is very
                    // nearly the entire planetary albedo: the deck's condensate path runs to
                    // 1e5 g/m2 against a cwp_tau of 100, so refl = tau/(tau+2) saturates and
                    // every cloudy column returns alpha_cloud to four decimals. Whatever
                    // this number is set to IS the model's albedo. It was a bare literal.
                    const double alpha_cloud     = m.albedo_cloud;   // thick cloud-top SW albedo
                    constexpr double f_ice_sw    = 0.50;   // ice SW reflectivity weight vs liquid
                    constexpr double cwp_tau     = 100.0;  // g/m2 per unit effective optical thickness
                    const double cwp_sw = lwp_col + f_ice_sw * iwp_col;          // [g/m2]
                    const double tau    = cwp_sw / cwp_tau;
                    const double refl   = tau / (tau + 2.0);                     // gentle saturation (0.5 at tau=2)
                    const double a0     = m.albedo.y[j][k];                      // surface (ice-feedback) albedo
                    if (alpha_cloud > a0)
                        m.albedo.y[j][k] = a0 + (alpha_cloud - a0) * refl;
                }

                // ============================================================
                // ATHAD: two-stream (Schwarzschild) flux sweeps, replacing the
                // tridiagonal Thomas solve.
                //
                // WHAT WAS HERE. The inherited scheme assembled a tridiagonal system whose
                // rows were products of the layer emissivity — sub-diagonal
                // eps_{i-1}*sigma*T_{i-1}^4, diagonal -2*eps_i*sigma*T_i^4, super-diagonal
                // eps_{i+1}*sigma*T_{i+1}^4 — with a right-hand side built from differences
                // of cumulative transmitted sums (AA[i]-AA[i-1], CA[i]-CA[i-1]), and solved
                // it for a correction to sigma*T^4. Every one of those entries vanishes with
                // eps, and eps goes as p^2 through the pressure-broadened optical depth, so
                // in the upper atmosphere adjacent rows differ by orders of magnitude while
                // the right-hand side is a difference of two nearly equal cumulative sums.
                //
                // It did not survive a domain that reaches the radiating level. Measured:
                // shells of 260 km (top 0.017 bar) and 300 km (top 3.2e-4 bar) built a
                // finite, sane initial state and then produced NaN across the ENTIRE field
                // in the first call here — the field is finite in the diagnostic immediately
                // before and NaN in the one immediately after, in both runs. That put a
                // ceiling on the shell at 230 km, where the top layer is still optically
                // thick (tau ~ 6), which in turn meant the model emitted sigma*T_lid^4 out
                // of its own boundary and reported a PRESCRIBED temperature as its OLR.
                //
                // WHAT REPLACES IT. The plane-parallel grey transfer the tridiagonal system
                // was standing in for, integrated directly:
                //
                //     up[i] = up[i-1]*(1 - eps_i) + eps_i*sigma*T_i^4        (surface -> top)
                //     dn[i] = dn[i+1]*(1 - eps_i) + eps_i*sigma*T_i^4        (top -> surface)
                //
                // with up[i_mount] = sigma*T_surf^4 and dn above the top layer = 0 (no
                // downward flux from space). Radiative equilibrium of layer i is then the
                // statement that it absorbs what it emits, eps_i*(up[i-1] + dn[i+1]) =
                // 2*eps_i*sigma*T_i^4, and the emissivity CANCELS:
                //
                //     sigma*T_i^4 = (up[i-1] + dn[i+1]) / 2
                //
                // No division by eps anywhere, and the eps -> 0 limit is exactly right
                // rather than merely survivable: a transparent layer passes both streams
                // unchanged and takes the mean of what goes by. At the top, where dn -> 0,
                // it reduces to sigma*T^4 = up/2 — the classical skin temperature, which
                // this model has until now been PRESCRIBING as t_skin.
                //
                // The temperature and the fluxes depend on each other, so the pair is
                // iterated (Lambda iteration). Convergence is slow in optically thick
                // layers, which is the known weakness of the method — but there it is also
                // nearly a no-op, since up and dn both approach the local sigma*T^4 and the
                // update becomes a three-point average. The fluxes themselves, and hence
                // the OLR, are exact for the current temperature field at every iteration.
                constexpr int n_lambda = 4;

                for (int it = 0; it < n_lambda; it++) {

                    // Upward sweep. The surface is layer i_mount and emits as a black body;
                    // its temperature is set by the energy balance further down.
                    up[i_mount] = m.sigma * pow(T[i_mount], 4.0);
                    for (int i = i_mount + 1; i <= i_trop; i++) {
                        const double eps = m.epsilon.x[i][j][k];
                        up[i] = up[i-1] * (1.0 - eps) + eps * m.sigma * pow(T[i], 4.0);
                    }

                    // Downward sweep. Nothing comes down from space.
                    dn[i_trop + 1] = 0.0;
                    for (int i = i_trop; i >= i_mount + 1; i--) {
                        const double eps = m.epsilon.x[i][j][k];
                        dn[i] = dn[i+1] * (1.0 - eps) + eps * m.sigma * pow(T[i], 4.0);
                    }

                    // Radiative-equilibrium temperature of every atmospheric layer. The
                    // surface is excluded: it has its own energy balance below, which
                    // carries the shortwave, the geothermal flux and the turbulent flux.
                    for (int i = i_mount + 1; i <= i_trop; i++) {
                        const double emit = 0.5 * (up[i-1] + dn[i+1]);
                        // emit is a sum of non-negative fluxes, so this is defensive only —
                        // but a negative or NaN value here would propagate into every layer
                        // above through the next sweep.
                        // is_finite_safe, not std::isfinite: this file is compiled with
                        // -ffast-math, under which the standard predicate may be folded away.
                        T[i] = (emit > 0.0 && AtomUtils::is_finite_safe(emit))
                             ? pow(emit / m.sigma, 0.25)
                             : T[i];
                    }
                }

                // Final sweeps on the converged temperatures, so the fluxes written out
                // below are the ones this profile actually produces.
                up[i_mount] = m.sigma * pow(T[i_mount], 4.0);
                for (int i = i_mount + 1; i <= i_trop; i++) {
                    const double eps = m.epsilon.x[i][j][k];
                    up[i] = up[i-1] * (1.0 - eps) + eps * m.sigma * pow(T[i], 4.0);
                }
                dn[i_trop + 1] = 0.0;
                for (int i = i_trop; i >= i_mount + 1; i--) {
                    const double eps = m.epsilon.x[i][j][k];
                    dn[i] = dn[i+1] * (1.0 - eps) + eps * m.sigma * pow(T[i], 4.0);
                }

                for (int i = i_mount + 1; i <= i_trop; i++)
                    m.t.x[i][j][k] = T[i] / m.t_0;

                // ---- Surface energy balance (radiative-CONVECTIVE) ----
                // Unchanged in form. The downwelling long-wave is now dn[i_mount + 1] — the
                // flux the sweep above actually delivers to the ground — instead of a
                // separate hand-rolled loop that re-derived the same quantity.
                //  - c_H*(T_s - T_air1): bulk turbulent (sensible+latent) flux to the lowest
                //    air layer, which keeps the surface off pure-radiative overheating.
                //  - sigma*T_s^4 is linearised about the current T_s0 (one Newton step).
                const double L_down = dn[i_mount + 1];
                const double SW_abs = (1.0 - m.albedo.y[j][k]) * m.short_wave_radiation[j];
                const double T_air1 = T[i_mount + 1];                      // lowest air-layer T [K]
                const double T_s0   = T[i_mount];                          // linearisation point [K]
                const double c_H    = 15.0;                                // bulk turbulent transfer [W/m2/K]
                const double dsigT4 = 4.0 * m.sigma * T_s0 * T_s0 * T_s0;
                // ATHAD: the surface is molten, so it supplies heat from below as well as
                // absorbing it from above. A quenching magma ocean radiates far more than the
                // modern Earth's 0.09 W/m2, and at 1500 K this term is plausibly comparable to
                // the absorbed solar — omitting it would let the surface cool as if it were
                // rock. It enters the balance exactly as absorbed shortwave does.
                const double T_s    = (SW_abs + m.geothermal_flux + L_down
                                       - m.sigma * pow(T_s0, 4.0)
                                       + dsigT4 * T_s0 + c_H * T_air1) / (dsigT4 + c_H);
                m.t.x[i_mount][j][k] = T_s / m.t_0;

                // radiation.x is now the UPWARD LONG-WAVE FLUX at the top of each layer,
                // not sigma*T^4 of that layer. It is a diagnostic field (nothing feeds it
                // back into the dynamics), and as a flux it is the more useful one: it is
                // continuous across the surface by construction, so the 1-2-1 smoothing pass
                // that used to hide the surface kink is gone, and radiation.x[im-1] is the
                // OLR — which is what the mode-5 cloud diagnostic already assumed it was.
                m.radiation.x[i_mount][j][k] = m.sigma * pow(T_s, 4.0);
                for (int i = i_mount + 1; i <= i_trop; i++)
                    m.radiation.x[i][j][k] = up[i];

            }  // k
        }  // j

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for MultiLayerRadiation\n", elapsed.count() * 1e-9);
        cout << "      RadiationMultiLayer ended" << endl;
    }

private:
    cAtmosphereModel& m;
};

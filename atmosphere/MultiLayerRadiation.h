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

    // ATM_SFC_COUPLE (default 0 = the per-call Newton step below, byte-identical; 2026-11-27, user: "the
    // coupled ground-air flux", with "Physical + speed-up knob" and "L_atm" as the time unit).
    //
    // WHY. The ground had no working energy balance in the time loop. The Newton step below solves the
    // balance against the RADIATIVE-EQUILIBRIUM air (which densities() then overwrites), has no heat
    // capacity, and apply_radiative_heating keeps 5 % of it per call, so the ground relaxed on an iteration
    // count, not on a time. The only exchange with the REAL air was the radial damp_wiggles(t) at the ground
    // row (removed by ATM_T_DAMP_GROUND): on ATHAD_PERID it cooled the ground (air above colder), on ATHAD
    // it warmed it (level 1 heated above the ground after the diffusion ramp), in neither case conserving
    // energy (ATM_SFC_DIAG + stage probes, ATHAD_PERID/CLAUDE.md, 2026-11-27).
    //
    // WHAT. With 1, inside apply_radiative_heating (radiation_mode 2, the time loop) only:
    //   ground   C dT_s/dt = SW + geo + Ldown - sigma T_s^4 - H,   H = c_H (T_s - T_1),  T_1 = the ACTUAL
    //            level-1 temperature as it entered; linearised backward Euler over tau = dt L_atm/u_0 * ACCEL
    //            (stable for any ACCEL; ACCEL -> infinity is the Newton step with the actual air);
    //   level 1  + H dt L_atm/u_0 / (rho_1 cp_1 dz_1)   (real time, not accelerated);
    //   both written AFTER the omega_rad blend, so the blend no longer touches the ground.
    // ATM_SFC_HEATCAP [J/m2/K] default 2.1e8 (a 50 m water layer), ATM_SFC_ACCEL default 1. At ACCEL 1 the
    // ground is practically constant over any affordable run: one iteration is 0.03 s (COND, PERID) to
    // 0.2 s (ATHAD) on L_atm (CLAUDE.md item 47), so +180 W/m2 moves it ~1e-8 K per iteration.
    // Other callers of run() (diagnostic refreshes that restore t) keep the Newton step.
    static bool sfc_couple_on() {
        static const bool v = [](){ const char* e = getenv("ATM_SFC_COUPLE");
                                    return e && atoi(e) != 0; }();
        return v;
    }
    static bool& sfc_coupled_call() { static bool b = false; return b; }   // set by apply_radiative_heating
    static std::vector<double>& sfc_Ts()  { static std::vector<double> v; return v; }   // [K], j*km + k
    static std::vector<double>& sfc_dT1() { static std::vector<double> v; return v; }   // [K], j*km + k
    static double sfc_heatcap() {
        static const double v = [](){ const char* e = getenv("ATM_SFC_HEATCAP");
                                      return (e && atof(e) > 0.0) ? atof(e) : 2.1e8; }();
        return v;
    }
    static double sfc_accel() {
        static const double v = [](){ const char* e = getenv("ATM_SFC_ACCEL");
                                      return (e && atof(e) > 0.0) ? atof(e) : 1.0; }();
        return v;
    }

    void run()
    {
        using namespace std;
        cout << endl << "      RadiationMultiLayer" << endl;
        auto begin = std::chrono::high_resolution_clock::now();

        const int j_max  = m.jm - 1;
        const int j_half = j_max / 2;

        // ATM_SFC_DIAG=<n> (print-only, default 0 = off; 2026-11-27, user: "add the surface-term
        // print and confirm it"). Every n-th call of this routine: the terms of the surface energy
        // balance below at k = 0 on four rows (90 N, 60 N, 30 N, equator), the net they give at the
        // linearisation point T_s0, and the same net with the air the model actually KEEPS — the
        // adiabat densities() hung from T_s0, plus the RK step, as it entered this routine — instead of the
        // radiative-equilibrium air this routine computes and densities() then overwrites (with
        // ATM_PROGNOSTIC_T=0). Question it answers: why the 12-bar ATHAD_PERID ground cools by a
        // fixed amount per call. Reads only; the state is untouched.
        static const int sfc_diag_n = [](){ const char* e = getenv("ATM_SFC_DIAG");
                                            return e ? std::max(0, atoi(e)) : 0; }();
        static long sfc_calls = 0;
        const bool sfc_diag_on = sfc_diag_n > 0 && (sfc_calls++ % sfc_diag_n == 0);
        const int sfc_rows[4] = { 0, j_half / 3, (2 * j_half) / 3, j_half };
        struct SfcRec { double Ts0, Ts, Ta_ad, Ta_re, SW, geo, Ld_re, Ld_ad, sT4, H_re, H_ad,
                        net_re, net_ad; bool set = false; };
        SfcRec sfc_rec[4];

        // ATM_SFC_COUPLE: active only on the call apply_radiative_heating flags.
        const bool   sfc_couple = sfc_couple_on() && sfc_coupled_call();
        const double dt_phys    = m.dt * m.L_atm / m.u_0;                  // [s] per iteration (L_atm)
        if (sfc_couple) {
            sfc_Ts().assign((size_t)m.jm * m.km, 0.0);
            sfc_dT1().assign((size_t)m.jm * m.km, 0.0);
        }

        // ---- latitude profiles (computed once; read-only in the parallel loop) ----

        // Surface albedo with a temperature-dependent ICE/SNOW FEEDBACK (replaces the former
        // fixed pole->equator parabola, which froze the feedback out). Base value by surface
        // type (ocean vs land); where the surface is cold, blend toward a bright snow/sea-ice
        // albedo over a smooth ramp around freezing. Recomputed every MLR call on the CURRENT
        // surface T, so it is a LIVE feedback: warming -> less ice -> lower albedo -> more
        // absorbed shortwave -> extra (polar-amplified) warming that the dynamics cannot mix
        // away (unlike a forcing/nudge). Constants are tunable. NOTE the cloud SW bump below
        // overwrites this per column where cloud is present, so its net reach is cloud-limited.
        // ATHAD_COND: the surface is a LIQUID WATER SEA at 503-513 K, and there is still no
        // ice-albedo feedback — the sea is 240 K above freezing, so the inherited ramp
        // toward snow/sea-ice through 275 -> 265 K remains dead code, exactly as it was in
        // ATHAD for the opposite reason (1200 K above the thresholds rather than 240).
        //
        // THE VALUE IS WATER'S SINCE 2026-08-20, AND IT WAS THE MAGMA OCEAN'S BEFORE THAT.
        // `albedo_surface` carried 0.08 with the justification "a quenching silicate melt is
        // dark, measured basaltic-melt albedos are 0.05-0.10" — ATHAD's surface, inherited
        // unchanged into a model whose surface is a sea. It is now 0.06, a water value. See
        // param.py for the derivation and for why the low-latitude weighting is the right
        // one on a zero-obliquity planet.
        //
        // This is the CLEAR-SKY value only; the cloud bump below raises it wherever the
        // model actually produces condensate — which in THIS fork is essentially everywhere,
        // so the number's reach is small and measured rather than assumed. That separation
        // is the point: a reflective deck must be EARNED by condensate the model generated,
        // not asserted as a constant albedo.
        //
        // It is a PARAMETER (albedo_surface), not a literal. It was written here as a
        // constexpr while the config carried an inert albedo_pole/albedo_equator pair that
        // nothing read — so the file said one thing and the configuration said another.
        const double alb_surface_clear = m.albedo_surface;    // clear-sky sea surface
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < m.jm; j++)
            for (int k = 0; k < m.km; k++)
                m.albedo.y[j][k] = alb_surface_clear;

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
            std::vector<double> up(m.im, 0.0), dn(m.im + 1, 0.0), T(m.im, 0.0),
                                Gf(m.im, 0.0),   // dB/dF for the direct solve
                                T_in, dn_ad;     // ATM_SFC_DIAG only: entering T, its down sweep
            int sfc_slot = -1;                   // ATM_SFC_DIAG: which of sfc_rows this j is
            if (sfc_diag_on)
                for (int s = 0; s < 4; s++) if (sfc_rows[s] == j) sfc_slot = s;

            const int i_trop  = m.im - 1;   // top layer
            const int i_mount = 0;          // surface / bottom layer

            for (int k = 0; k < m.km; k++) {

                for (int i = i_mount; i <= i_trop; i++)
                    T[i] = m.t.x[i][j][k] * m.t_0;
                if (sfc_slot >= 0 && k == 0) T_in = T;   // ATM_SFC_DIAG: the air as it entered
                const double T_air_act = T[i_mount + 1]; // ATM_SFC_COUPLE: the actual level-1 air

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
                // ATM_N_LAMBDA — the sweep count, made a knob so the Lambda iteration can be
                // run to convergence and compared against the closed form below. 4 is the
                // inherited value and the default, so this is bit-identical when unset.
                static const int n_lambda = [](){
                    const char* e = getenv("ATM_N_LAMBDA"); return e ? atoi(e) : 4; }();

                // ---- DIRECT SOLUTION, no iteration at all (ATM_RAD_DIRECT) ----
                // Ported from ATHAD, README items 30 and 71.
                //
                // The Lambda iteration below is a Jacobi relaxation on an im-link chain:
                // information moves one layer per sweep, so it needs O(N^2) sweeps, and
                // n_lambda = 4 is an Earth constant — a loop bound nobody reads as a physical
                // assumption, adequate on a 1 bar column and under-converged on a thick one.
                //
                // It does not need iterating. In radiative equilibrium the NET flux is
                // constant with height, and that closes the system in one pass. Writing
                // a_i = 1 - eps_i/2, b_i = eps_i/2, a layer in equilibrium
                // (B_i = (U_{i-1} + D_i)/2, which is the eps-cancelled condition above)
                // transfers
                //     U_i     = a_i U_{i-1} + b_i D_i
                //     D_{i-1} = b_i U_{i-1} + a_i D_i
                // and a_i + b_i = 1 exactly, so U - D takes the same value at both faces:
                // F is constant, which is the definition of radiative equilibrium and here
                // falls out of the discretisation rather than being imposed. Then
                //     U_{i-1} = B_i + F/(2 a_i),   D_i = B_i - F/(2 a_i)
                // and eliminating the fluxes leaves an explicit march for the source
                // function, LINEAR in F:
                //     B_{i+1} - B_i = (F/2) [ (1 - eps_i)/a_i - 1/a_{i+1} ]
                // March G = dB/dF from the surface boundary U_m = sigma*T_s^4 (so
                // B_{m+1} = sigma*T_s^4 - F/(2 a_{m+1}), i.e. G_{m+1} = -1/(2 a_{m+1})),
                // then close on the top boundary D_n = 0 (so B_n = F/(2 a_n)):
                //     F = sigma*T_s^4 / ( 1/(2 a_n) - G_n )
                // and B_i = sigma*T_s^4 + F G_i. Two O(N) passes, exact.
                //
                // Validated in ATHAD against the Lambda iteration run to convergence on
                // synthetic columns — identical to all printed digits transparent, optically
                // thick and 250-bar-graded — and it reproduces the two exact limits: with
                // eps -> 0 it gives F = sigma*T_s^4 and T_i = T_s/2^(1/4), the classical skin
                // temperature of a freely radiating surface.
                //
                // **DEFAULT ON, MATCHING ATHAD SINCE 2026-08-22 (its item 71).
                // `ATM_RAD_DIRECT=0` restores the 4-sweep Lambda iteration; every OLR this
                // tree recorded before the port used it.**
                //
                // WHY THE DEFAULT IS ON HERE, AND WHAT IS NOT MEASURED. The argument that
                // flipped it in ATHAD is not a tuning result: the Lambda iteration converges
                // MONOTONICALLY onto this closed form, which makes the closed form the answer
                // the sweeps are trying to reach rather than an alternative to them. Measured
                // there at 40 iterations, 24 threads:
                //
                //     n_lambda = 4 (was shipped) .... 243.43 W/m2    +15.06 %
                //     n_lambda = 64 ................. 227.55 W/m2     +7.55 %
                //     n_lambda = 512 ................ 212.54 W/m2     +0.46 %
                //     ATM_RAD_DIRECT=1 .............. 211.57 W/m2      0.00 %
                //
                // The exact answer is also free: 280.97 s against 277.34 s for the 4-sweep arm
                // and 467 s for n_lambda = 512. **THE SIZE OF THE SHIFT IS ATHAD'S 250 BAR
                // COLUMN AND IS NOT MEASURED IN THIS FORK** — this atmosphere is a quarter of
                // that pressure, so its optical depth per layer, and hence how far 4 sweeps
                // fall short, is a different number. Only the direction is inherited.
                static const bool rad_direct = [](){
                    const char* e = getenv("ATM_RAD_DIRECT");
                    return !(e && atoi(e) == 0); }();

                if (rad_direct) {
                    const int m1 = i_mount + 1, ntop = i_trop;
                    if (ntop > m1) {
                        auto aof = [&](int i){ return 1.0 - 0.5 * m.epsilon.x[i][j][k]; };
                        const double B_s = m.sigma * pow(T[i_mount], 4.0);
                        Gf[m1] = -0.5 / aof(m1);
                        for (int i = m1; i < ntop; i++) {
                            const double e = m.epsilon.x[i][j][k];
                            Gf[i+1] = Gf[i] + 0.5 * ((1.0 - e) / aof(i) - 1.0 / aof(i+1));
                        }
                        const double den = 0.5 / aof(ntop) - Gf[ntop];
                        const double F   = (std::fabs(den) > 1.0e-30) ? B_s / den : 0.0;
                        for (int i = m1; i <= ntop; i++) {
                            const double B = B_s + F * Gf[i];
                            if (B > 0.0 && AtomUtils::is_finite_safe(B))
                                T[i] = pow(B / m.sigma, 0.25);
                        }
                    }
                }
                else
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
                double T_s;
                if (!sfc_couple) {
                    T_s = (SW_abs + m.geothermal_flux + L_down
                           - m.sigma * pow(T_s0, 4.0)
                           + dsigT4 * T_s0 + c_H * T_air1) / (dsigT4 + c_H);
                } else {
                    // ATM_SFC_COUPLE (see sfc_couple_on): heat capacity, the actual air, a time.
                    const double H   = c_H * (T_s0 - T_air_act);
                    const double F   = SW_abs + m.geothermal_flux + L_down
                                     - m.sigma * pow(T_s0, 4.0) - H;
                    const double tau = dt_phys * sfc_accel();
                    T_s = T_s0 + tau * F / (sfc_heatcap() + tau * (dsigT4 + c_H));
                    const int    m1   = i_mount + 1;
                    const double rho1 = m.r_humid.x[m1][j][k];
                    const double cp1  = AtmMixture::cp_of(m.c.x[m1][j][k], m.co2.x[m1][j][k],
                                                          T_air_act, m.m_comp.M_bg);
                    const double dz1  = 0.5 * (m.get_layer_height(m1 + 1) - m.get_layer_height(i_mount));
                    const double cap1 = rho1 * cp1 * dz1;                  // [J/m2/K]
                    sfc_Ts()[(size_t)j * m.km + k]  = T_s;
                    sfc_dT1()[(size_t)j * m.km + k] = (cap1 > 0.0 && AtomUtils::is_finite_safe(cap1))
                                                    ? H * dt_phys / cap1 : 0.0;
                }
                m.t.x[i_mount][j][k] = T_s / m.t_0;

                // ATM_SFC_DIAG (print-only): record the balance; the adiabat's down-welling
                // long-wave is the same sweep on the entering temperatures.
                if (sfc_slot >= 0 && k == 0) {
                    dn_ad.assign(m.im + 1, 0.0);
                    for (int i = i_trop; i >= i_mount + 1; i--) {
                        const double eps = m.epsilon.x[i][j][k];
                        dn_ad[i] = dn_ad[i+1] * (1.0 - eps) + eps * m.sigma * pow(T_in[i], 4.0);
                    }
                    SfcRec& r = sfc_rec[sfc_slot];
                    r.Ts0 = T_s0;  r.Ts = T_s;
                    r.Ta_re = T_air1;  r.Ta_ad = T_in[i_mount + 1];
                    r.SW = SW_abs;  r.geo = m.geothermal_flux;
                    r.Ld_re = L_down;  r.Ld_ad = dn_ad[i_mount + 1];
                    r.sT4 = m.sigma * pow(T_s0, 4.0);
                    r.H_re = c_H * (T_s0 - r.Ta_re);  r.H_ad = c_H * (T_s0 - r.Ta_ad);
                    r.net_re = r.SW + r.geo + r.Ld_re - r.sT4 - r.H_re;
                    r.net_ad = r.SW + r.geo + r.Ld_ad - r.sT4 - r.H_ad;
                    r.set = true;
                }

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

        if (sfc_diag_on) {
            printf("      [sfc] call %ld, k = 0; W/m2, K. net = SW + geo + Ldown - sigma Ts0^4 - 15 (Ts0 - T_air1);"
                   " _re = with the radiative-equilibrium air (what the balance uses), _ad = with the air as it"
                   " entered (densities() adiabat + the RK step), which densities() keeps\n", sfc_calls - 1);
            for (int s = 0; s < 4; s++) {
                const SfcRec& r = sfc_rec[s];
                if (!r.set) continue;
                printf("      [sfc] lat %5.1f  Ts0 %7.2f -> %7.2f (%+6.3f)  T_air1 ad %7.2f re %7.2f"
                       "  SW %6.1f geo %5.1f  Ldown re %6.1f ad %6.1f  sT4 %6.1f  H re %6.1f ad %6.1f"
                       "  net re %+7.2f ad %+7.2f\n",
                       90.0 - 180.0 * sfc_rows[s] / (double)j_max, r.Ts0, r.Ts, r.Ts - r.Ts0,
                       r.Ta_ad, r.Ta_re, r.SW, r.geo, r.Ld_re, r.Ld_ad, r.sT4, r.H_re, r.H_ad,
                       r.net_re, r.net_ad);
            }
        }

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for MultiLayerRadiation\n", elapsed.count() * 1e-9);
        cout << "      RadiationMultiLayer ended" << endl;
    }

private:
    cAtmosphereModel& m;
};

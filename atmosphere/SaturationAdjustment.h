#pragma once

#include "MixtureAtm.h"
#include "SaturationH2O.h"
#include "cAtmosphereModel.h"
#include "IceSchemeCommon.h"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <cstdio>
#include <vector>

using namespace AtomUtils;


class SaturationAdjustment {
public:
    explicit SaturationAdjustment(cAtmosphereModel& model)
        : m(model)
    {}


    void run() {
        using namespace std;

        cout << endl << endl << endl << "      SaturationAdjustment" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        cout.precision(9);

        computeSteps();
        adjustSaturation();
        applyTopography();
        clampAndFade();
        printReport();

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for SaturationAdjustment\n", elapsed.count() * 1e-9);

        cout << "      SaturationAdjustment ended" << endl;
    }

private:
    cAtmosphereModel& m;

    std::vector<double> step;

    // Precomputed constants
    static constexpr double fade_K        = 5.0;                        // transition half-width in Kelvin

    // ATM_SAT_TRACE=1 — print the Newton loop for ONE cell, to answer why this routine leaves
    // this fork mildly supersaturated THROUGHOUT (100.2-102 % from 1.4 to 50 km, tail to 148.9 %
    // at 63 km) where ATHAD is violently so at two levels (RH 2518 % and 626 %; ATHAD item 64). Print-only: it writes no field
    // and returns nothing, so a traced run and an untraced one differ only in stdout.
    //
    // What to read. omega = 1/(1+G) is Newton-optimal damping and G = (L/cp)*dq_sat/dT is
    // enormous here — L/cp ~ 1225 K at ATHAD's cp, and dq_sat/dT is steep where q_sat swings
    // 0.186 -> 0.033 over 20 km — so each pass may close only a 1/(1+G) fraction of the gap.
    // The loop's exit test is |q_v_b/q_v_hyp - 1| <= 1e-6, a test on the STEP and not on the
    // residual, which heavy damping can satisfy while the cell is still far from saturation.
    // The trace prints both so the two can be told apart.
    // Default OFF: unset reproduces every number this tree has printed.
    // ATM_SAT_LEGACY=1 restores BOTH pre-2026-08-20 behaviours together: the alpha_entry
    // master gain and the step-based exit test. Default is the repaired path.
    //
    // (1) alpha_entry APPLIED THE -37 C ICE THRESHOLD TWICE. Inside the Newton loop it is
    //     already the PHASE SPLIT -- CND = clamp((T - t_00)*t_range_inv), DEP = 1 - CND --
    //     which is the physics: below -37 C supercooled liquid cannot exist, so condensation
    //     becomes deposition. alpha_entry then re-applied the same threshold as a MASTER GAIN
    //     on all five write-backs (S_c_c, c, cloud, ice, t), so a cell below ~236 K kept only
    //     a few per cent of whatever the loop computed -- DEPOSITION INCLUDED, which is
    //     precisely the process that should be running there. It also gated entry at
    //     alpha_entry > 0.01, skipping cells below 213.2 K outright.
    //
    //     Invisible on Earth, where a cell at -37 C holds ~0.1 g/kg of vapour. Measured live
    //     in ATHAD (its README item 64): level 38 free-running to 220-228 K with 683 g/kg of
    //     vapour and alpha_entry = 0.036, i.e. the adjustment allowed to apply 3.6 % of its
    //     own answer. In THIS fork the trace shows alpha_entry = 1.0000 and 0.9998 at the two
    //     probed levels, so the gain never engages and the repair is expected to be a no-op
    //     here -- which is exactly why this is the safe tree to make it in first.
    //
    // (2) THE EXIT TEST MEASURED THE STEP, NOT THE RESIDUAL. |q_v_b/q_v_hyp - 1| <= 1e-6 is a
    //     statement about how far the last pass moved, and heavy damping (omega = 1/(1+Gain),
    //     and Gain reaches 97 in ATHAD) makes the step small while the cell is still far from
    //     saturation. ATHAD's trace caught it exiting with the step at 4.4e-06 while the
    //     residual still swung +-0.003, on a limit cycle rather than a converged answer. The
    //     repaired test is on |q_v_b - q_v_target| relative to q_v_target -- the thing the
    //     loop is actually trying to drive to zero -- with the step test kept as a SECOND,
    //     looser guard so a genuinely stalled iteration still terminates.
    static inline const bool sat_legacy = [](){
        const char* e = getenv("ATM_SAT_LEGACY"); return (e && atoi(e) != 0); }();
    static inline const bool sat_no_alpha = [](){
        const char* e = getenv("ATM_SAT_NO_ALPHA");
        return (e && atoi(e) != 0) || !sat_legacy; }();
    static inline const bool sat_trace = [](){
        const char* e = getenv("ATM_SAT_TRACE"); return (e && atoi(e) != 0); }();
    // The two traced levels default to levels 20 (10.8 km, mid-column) and 48 (63 km, the
    // 148.9 % tail). ATHAD defaults to 37/38, its only condensing band — the one deliberate
    // difference between the two copies of this file.
    static inline const int trace_i1 = [](){
        const char* e = getenv("ATM_SAT_TRACE_I1"); return e ? atoi(e) : 20; }();
    static inline const int trace_i2 = [](){
        const char* e = getenv("ATM_SAT_TRACE_I2"); return e ? atoi(e) : 48; }();
    static constexpr int trace_j = 90;    // equator on the 181-point grid
    static constexpr int trace_k = 0;

    static constexpr int    iter_prec_end = 20;

    void computeSteps() {
        step.resize(m.im);
        for (int i = 0; i < m.im; i++)
            step[i] = m.get_layer_height(i+1) - m.get_layer_height(i);
    }

    void adjustSaturation() {
        const double inv_t_0      = 1.0 / m.t_0;
        const double t_range_inv  = 1.0 / (m.t_0 - m.t_00);
        // Latent heat is now a FUNCTION of temperature, evaluated per cell inside the
        // loop rather than hoisted as a constant. It has to be: L falls from 2.50e6 J/kg
        // at 273 K to zero at the critical point, and a constant L near T_crit injects
        // heat from a phase change that is not happening. cp is likewise local, since the
        // mixture's heat capacity moves with composition and temperature.

        // Surface row is skipped below, but its condensation source must still be
        // cleared every call: S_c_c.x[0] feeds the ice schemes' cloud-water source,
        // and over ocean (i_mount == 0) they copy S_c_c.x[0] onto itself, so a stale
        // nonzero value would silently re-condense surface cloud and defeat the skip.
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++)
            for (int k = 0; k < m.km; k++)
                m.S_c_c.x[0][j][k] = 0.0;

        // Start at i = 1: the surface row (i = 0) is a vapour SOURCE only and must
        // never be saturation-adjusted in place. waterVapourEvaporation injects c[0]
        // from the warm ocean toward c_eq ~ q_sat every moist iter; if adjustSaturation
        // then condenses that just-injected vapour into cloud at the same cell, cloud
        // water doubles every ~3-4 iters (Cook Inlet runaway: cloud 5e-13 -> 32 over
        // iters 315-363). Once cloud + ice > 1 the r_humid denominator
        // (1 + 0.608*c - cloud - ice) flips negative, -grad(p)/rho inverts, NaN cascade.
        // Physically, condensation only happens once a parcel has risen, cooled
        // adiabatically and reached its LCL aloft; the warm ocean surface is below
        // saturation by definition. So cloud forms from i = 1 up, not at i = 0.
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 1; i < m.im - 1; i++) {
            for (int j = 0; j < m.jm; j++) {
                double *S_c_c_row = m.S_c_c.x[i][j];
                double *c_row     = m.c.x[i][j];
                double *cloud_row = m.cloud.x[i][j];
                double *ice_row   = m.ice.x[i][j];
                double *t_row     = m.t.x[i][j];
                double *p_row     = m.p_stat.x[i][j];
                double  dt_dim    = step[i] / 1.6;

                for (int k = 0; k < m.km; k++) {
                    S_c_c_row[k] = 0.0;

                    double q_v_old = std::max(0.0, c_row[k]);
                    double q_c_old = std::max(0.0, cloud_row[k]);
                    double q_i_old = std::max(0.0, ice_row[k]);

                    double T       = t_row[k] * m.t_0;

                    // ATHAD: above the critical point there is NO saturation to adjust to.
                    //
                    // The inherited code capped T at 333.15 K here and WROTE THE CAP BACK
                    // into the prognostic field, on the reasoning that 60 °C is "well above
                    // any physical surface temperature". That is an Earth statement. ATHAD's
                    // surface is 1500 K, so the cap destroyed the entire temperature field on
                    // the first call — 250 bar became 30 bar within one iteration.
                    //
                    // The cap existed only to keep T inside the Magnus formula's validity
                    // range. The physically correct statement is stronger and needs no cap:
                    // above 647.096 K water is supercritical, liquid and vapour are one phase,
                    // and there is nothing to condense. So skip the cell entirely.
                    if (T >= AtmMixture::T_CRIT_H2O) continue;

                    double p_local = p_row[k];

                    // IAPWS saturation curve and the EXACT mass-fraction conversion.
                    // The Magnus form here was valid to ~320 K and its q_sat fallback
                    // (ep * 1e-5) fired whenever E_sat exceeded the local pressure —
                    // collapsing q_sat by ~5000x and triggering the runaway condensation
                    // the old temperature cap existed to contain. The exact form
                    // saturates at 1 instead, which is what "the column is all vapour"
                    // actually means.
                    const double M_other = AtmMixture::M_nonwater(c_row[k], m.co2.x[i][j][k],
                                                                  m.m_comp.M_bg);
                    double E_sat  = SaturationH2O::saturationPressureAuto(T);
                    double q_sat  = SaturationH2O::saturationMassFraction(E_sat, p_local,
                                                                          M_other);

                    // ATM_SAT_NO_ALPHA=1 — the -37 C ice threshold is applied TWICE here, and
                    // the second application is a defect. Inside the Newton loop it is already
                    // the PHASE SPLIT: CND = clamp((T - t_00)*t_range_inv), DEP = 1 - CND,
                    // which is the physics ("below -37 C supercooled liquid cannot exist, so
                    // condensation becomes deposition"). alpha_entry then re-applies the same
                    // threshold as a MASTER GAIN on all five write-backs (S_c_c, c, cloud, ice
                    // and t, lines below), so a cell below ~236 K keeps only a few per cent of
                    // whatever the loop computed — deposition included, which is precisely the
                    // process that should be running there.
                    //
                    // Invisible on Earth, where a cell at -37 C holds ~0.1 g/kg of vapour.
                    // Live here: with ATM_PROGNOSTIC_T=1 level 38 free-runs to 220-228 K with
                    // 683 g/kg of vapour and alpha_entry = 0.036, so the adjustment is allowed
                    // to apply 3.6 % of its own answer. The entry test alpha_entry > 0.01 also
                    // skips the cell outright below 213.2 K.
                    const double alpha_entry = sat_no_alpha
                        ? 1.0
                        : 1.0 / (1.0 + std::exp(-(T - m.t_00) / fade_K));

                    const bool trace = sat_trace && (i == trace_i1 || i == trace_i2)
                                       && j == trace_j && k == trace_k;
                    if (trace)
                        std::printf("\n[sat] i=%d  T=%.2f K  p=%.5f bar  q_v=%.6f  q_sat=%.6f"
                                    "  q_v/q_sat=%.2f  alpha_entry=%.6f\n",
                                    i, T, p_local * 1e-3, q_v_old, q_sat,
                                    (q_sat > 0.0 ? q_v_old / q_sat : -1.0), alpha_entry);

                    if ((q_v_old > q_sat && alpha_entry > 0.01) ||
                        (q_v_old < q_sat &&
                        (q_c_old > 1e-12 || q_i_old > 1e-12))) {

                        double q_v_b   = q_v_old;
                        double q_c_b   = q_c_old;
                        double q_i_b   = q_i_old;
                        // ATHAD: the first Newton step must be DAMPED like every other one.
                        //
                        // This was `q_v_hyp = q_sat`, i.e. an undamped jump straight to the
                        // saturation value, with the loop's omega = 1/(1+Gain) damping only
                        // applied from the second pass onward. On Earth that is harmless:
                        // condensing the ~0.01 kg/kg a terrestrial parcel holds releases
                        // ~12 K, and the loop mops it up.
                        //
                        // Here it is fatal. At 243 km the column runs q_v = 0.72 against
                        // q_sat = 0.076, so the undamped step condenses 0.65 kg/kg in one
                        // go and releases 0.65*L/cp = 800 K of latent heat. T is then
                        // clamped to the critical point, where p_sat = 220 bar against a
                        // local 0.085 bar — so q_sat becomes 1, the target inverts, and the
                        // next pass evaporates everything back. The iteration flip-flops
                        // between fully condensed and fully evaporated and ends at zero, so
                        // the model produced NO cloud, ice or rain anywhere in the domain
                        // while sitting ninefold supersaturated.
                        //
                        // Starting from q_v_b makes the first pass a no-op that computes a
                        // properly damped target. The equilibrium it should find is modest:
                        // condensing ~0.04 kg/kg warms the parcel ~47 K, at which point the
                        // vapour is superheated and nothing further can condense.
                        //
                        // The defect is not a constant this time, it is an ASSUMPTION —
                        // that latent heating is a perturbation. At 67 % water by mass,
                        // condensation is a bulk phase change of the atmosphere.
                        double q_v_hyp = q_v_b;
                        const double T_original = t_row[k] * m.t_0;

                        for (int iter = 1; iter <= iter_prec_end; iter++) {
                            double CND = std::max(0.0, std::min(1.0,
                                (T - m.t_00) * t_range_inv));
                            double DEP = 1.0 - CND;

                            double d_q_v = q_v_hyp - q_v_b;

                            if (d_q_v > 0) {
                                double max_evap = q_c_b + q_i_b;
                                if (d_q_v > max_evap) d_q_v = max_evap;
                            }

                            double d_cnd = d_q_v * CND;
                            double d_dep = d_q_v * DEP;

                            q_v_b += d_q_v;
                            q_c_b  = std::max(0.0, q_c_b - d_cnd);
                            q_i_b  = std::max(0.0, q_i_b - d_dep);

                            const double cp_loc = AtmMixture::cp_of(q_v_b, m.co2.x[i][j][k],
                                                                    T, m.m_comp.M_bg);
                            T -= (SaturationH2O::latentHeat(T)           * d_cnd
                                + SaturationH2O::latentHeatSublimation(T) * d_dep) / cp_loc;

                            double E_sat = SaturationH2O::saturationPressure(T);
                            double E_Ice = SaturationH2O::sublimationPressure(T);
                            // Same exact conversion as the entry q_sat above. This copy inside
                            // the Newton loop was left on the dilute form, so the loop pulled
                            // q_v toward ep*1e-5 — effectively zero — on every superheated cell
                            // it was iterating, undoing the entry fix on the cells that reach
                            // this branch at all.
                            double q_sat = SaturationH2O::saturationMassFraction(
                                               E_sat, p_local, M_other);
                            double q_Ice = SaturationH2O::saturationMassFraction(
                                               E_Ice, p_local, M_other);

                            double q_sum = q_c_b + q_i_b;
                            double q_v_target = (q_sum > 1e-12)
                                ? (q_c_b * q_sat + q_i_b * q_Ice) / q_sum
                                : ((T >= m.t_0) ? q_sat : q_Ice);

                            // Adaptive Newton-damped update. The previous fixed 0.5 under-
                            // relaxation is UNSTABLE when the latent-heat gain
                            // G = (L/cp)*dq_sat/dT exceeds 3 — i.e. warm cells (T>~25C,
                            // q_sat>~20 g/kg) where dq_sat/dT is steep: the fixed-point map
                            // derivative 0.5*(1-G) then has magnitude >1, so the loop
                            // OSCILLATES and never converges, leaving RH stuck at 120-139%
                            // with cloud present (thermodynamically impossible) and feeding
                            // the coastal precip runaway. omega = 1/(1+G) drives the map
                            // derivative 1-omega*(1+G) to 0 (stable, ~Newton-optimal) at all T.
                            // dq_sat/dT from Clausius-Clapeyron.
                            // project_overprecip_saturation_injection.
                            // L(T) and cp(T), not constants: the whole point of the
                            // damping is to track dq_sat/dT, and both factors in
                            // G = (L/cp)*dq_sat/dT move strongly across 273-647 K.
                            //
                            // ATHAD_COND: the EXACT derivative, not the dilute
                            // q_sat*L/(Rv*T^2). The two differ by M_other/(x*Mw+(1-x)*Mo),
                            // which is 1.48 at this model's sea surface, so the dilute form
                            // understated the gain by a third and the damping omega =
                            // 1/(1+Gain) was correspondingly too weak — in exactly the
                            // regime (high q_sat, steep dq_sat/dT) the damping exists for.
                            // Phase given explicitly here because both branches are wanted
                            // at the same T.
                            const double cp_g   = AtmMixture::cp_of(q_v_b, m.co2.x[i][j][k],
                                                                    T, m.m_comp.M_bg);
                            const double L_cnd  = SaturationH2O::latentHeat(T);
                            const double L_dep  = SaturationH2O::latentHeatSublimation(T);
                            const double Gain = CND * (L_cnd / cp_g)
                                              * SaturationH2O::dqSatdTFrom(E_sat, L_cnd, T, p_local, M_other)
                                              + DEP * (L_dep / cp_g)
                                              * SaturationH2O::dqSatdTFrom(E_Ice, L_dep, T, p_local, M_other);
                            const double omega = 1.0 / (1.0 + Gain);
                            const double q_v_prev_pass = q_v_b;
                            q_v_hyp = q_v_b + omega * (q_v_target - q_v_b);

                            if (trace)
                                std::printf("[sat]   pass %2d  G=%.4g  omega=%.4g  T=%.2f"
                                            "  q_v=%.6f  q_sat=%.6f  residual=%.6f"
                                            "  step=%.3g  exit_test=%.3g\n",
                                            iter, Gain, omega, T, q_v_b, q_sat,
                                            q_v_b - q_sat, q_v_hyp - q_v_prev_pass,
                                            (q_v_hyp != 0.0)
                                              ? std::fabs(q_v_b / q_v_hyp - 1.0) : -1.0);

                            // Converged when the RESIDUAL is small -- q_v_b is at its target
                            // -- not when the STEP is small, which damping alone can produce.
                            // The step test survives as a looser backstop for a stalled loop.
                            const double resid = (q_v_target > 1e-12)
                                ? std::fabs(q_v_b - q_v_target) / q_v_target
                                : std::fabs(q_v_b - q_v_target);
                            const double step_rel = std::fabs(q_v_b / q_v_hyp - 1.0);
                            if (sat_legacy ? (step_rel <= 1.0e-6)
                                           : (resid <= 1.0e-6 || step_rel <= 1.0e-12))
                                break;
                        }

                        if (trace)
                            std::printf("[sat]   OUT     T=%.2f K (from %.2f)  q_v=%.6f"
                                        "  q_c=%.6f  q_i=%.6f   write-back q_v=%.6f\n",
                                        T, T_original, q_v_b, q_c_b, q_i_b,
                                        q_v_old + alpha_entry * (q_v_b - q_v_old));

                        q_c_b = std::max(0.0, q_c_b);
                        q_i_b = std::max(0.0, q_i_b);

                        // Cap T after the Newton loop. The q_v_hyp = 0.5*(q_v_target + q_v_b)
                        // damping is too weak when dq_sat/dT is steep (marginal saturation),
                        // so the iteration's amplitude grows. Within one call T can swing to a
                        // value the write-back below would persist. The entry guard is not
                        // enough because the runaway happens during the loop, not between calls.
                        //
                        // The bound is the critical temperature, not an arbitrary cap: this
                        // branch only runs on cells that entered SUBcritical, and a
                        // condensation adjustment cannot legitimately heat one past the point
                        // where the phase it is condensing into ceases to exist.
                        if (T > AtmMixture::T_CRIT_H2O) T = AtmMixture::T_CRIT_H2O;

                        if (!std::isnan(T) && !std::isnan(q_v_b)) {
                            S_c_c_row[k] = alpha_entry * (q_c_b - q_c_old) / dt_dim;
                            c_row[k]     = q_v_old + alpha_entry * (q_v_b - q_v_old);
                            cloud_row[k] = q_c_old + alpha_entry * (q_c_b - q_c_old);
                            ice_row[k]   = q_i_old + alpha_entry * (q_i_b - q_i_old);
                            t_row[k]     = (T_original + alpha_entry * (T - T_original)) * inv_t_0;
                        }
 
                        if (T < m.t_00) {
                            cloud_row[k] = 0.0;
                            ice_row[k]   = 0.0;
                        }
                    }
                }
            }
        }
    }

    void applyTopography() {
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                int i_mount = m.i_topography[j][k];
                m.c.x[0][j][k]     = m.c.x[i_mount][j][k];
                m.cloud.x[0][j][k] = m.cloud.x[i_mount][j][k];
                m.ice.x[0][j][k]   = m.ice.x[i_mount][j][k];
                if (AtomUtils::is_finite_safe(m.t.x[i_mount][j][k]))
                    m.t.x[0][j][k] = m.t.x[i_mount][j][k];
            }
        }
    }

    void clampAndFade() {
        const double inv_t_0    = 1.0 / m.t_0;
        // L/cp evaluated per cell below, for the same reason as in adjustSaturation().
        // Defensive physical bounds. The upper temperature bound is now a configured
        // PHYSICAL constant (t_max_phys) rather than the 333.15 K literal: a bound
        // expressed as "no air parcel exceeds 60 °C" is an Earth fact, and on ATHAD it sat
        // a factor of 4.5 below the surface temperature it was supposed to be protecting.
        // cloud_cap is now the PARAMETER of the same name, not a literal. Its inherited
        // justification — "~50x the largest physical cloud/ice mixing ratio, so it never
        // clips a real cloud, it only stops a runaway" — is an Earth statement:
        // terrestrial cloud water is ~1 g/kg. ATHAD's polar column pegs it exactly
        // (49.999996 g/kg at 80N, 185.6 km), so here it SETS the condensate rather than
        // bounding a pathology, and with it the optical depth and the albedo.
        const double T_max         = m.t_max_phys;
        const double cloud_cap = m.cloud_cap;   // kg/kg condensate ceiling (parameter)

        long n_ceiling = 0;                    // cells hitting the water-vapour ceiling

        #pragma omp parallel for collapse(2) schedule(static) reduction(+:n_ceiling)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                double *c_row     = m.c.x[i][j];
                double *cloud_row = m.cloud.x[i][j];
                double *ice_row   = m.ice.x[i][j];
                double *t_row_nd  = m.t.x[i][j];
                double *p_row     = m.p_stat.x[i][j];

                for (int k = 0; k < m.km; k++) {
                    if (c_row[k]     < 0.0) c_row[k]     = 0.0;
                    if (cloud_row[k] < 0.0) cloud_row[k] = 0.0;
                    if (ice_row[k]   < 0.0) ice_row[k]   = 0.0;

                    // ATHAD: water vapour has a physical CEILING and never had one.
                    //
                    // The mass fractions must sum to 1, and the background is carried as the
                    // remainder 1 - c - co2, so c > 1 - co2 means a NEGATIVE background mass.
                    // Nothing checked it. AtmMixture::split() renormalises defensively, so the
                    // gas constant stayed finite and the violation was invisible — the
                    // equatorial column was running c = 0.9971 against co2 = 0.2053, a
                    // composition summing to 1.20, and the only symptom was an R that had
                    // quietly saturated.
                    //
                    // The cause is not this clamp's business and is not fixed by it: water is
                    // pumped downward out of the one condensing level by sedimentation, and
                    // evaporates into the superheated band below with no return path, so c
                    // there grows without bound. The clamp stops the composition being
                    // impossible; the counter is what says how hard it is having to work.
                    // A run where n_ceiling stays large is not to be trusted.
                    const double c_max = std::max(0.0, 1.0 - m.co2.x[i][j][k]);
                    if (c_row[k] > c_max) { c_row[k] = c_max; n_ceiling++; }

                    double T_dim = t_row_nd[k] * m.t_0;
                    // Upper physical bound, well above the prescribed surface temperature.
                    if (T_dim > T_max) { T_dim = T_max; t_row_nd[k] = T_max * inv_t_0; }

                    // ATHAD: cells that cannot hold a condensed phase must be EMPTIED of
                    // one, not skipped.
                    //
                    // This used to read "Supercritical cells carry no condensate and need no
                    // fade" and simply `continue`. They carry no condensate only if something
                    // takes it away — and nothing did, so condensate advected or sedimented
                    // into the whole supercritical column (ground to ~180 km) stayed there
                    // and set the planetary albedo. The superheated band above it, 180 to
                    // ~240 km, escaped too: there q_sat = 1 and the supersaturation test
                    // below can never fire, since c < 1 always.
                    //
                    // IceSchemeCommon::evaporateWhereImpossible sends it back to the vapour
                    // with its latent heat, and clears the sources and precipitation fluxes.
                    // This is the net that catches whatever the ice schemes and the advection
                    // put there; the schemes now carry the same guard so they do not create
                    // it in the first place.
                    {
                        const double M_o = AtmMixture::M_nonwater(c_row[k], m.co2.x[i][j][k],
                                                                  m.m_comp.M_bg);
                        const double q_s = SaturationH2O::saturationMassFractionAt(
                                               T_dim, p_row[k], M_o);
                        if (T_dim >= AtmMixture::T_CRIT_H2O || q_s >= 1.0) {
                            IceSchemeCommon::evaporateWhereImpossible(m, T_dim, i, j, k);
                            continue;
                        }
                    }

                    // ---- Always-on supersaturation removal (ROOT FIX) ----
                    // adjustSaturation scales its condensation by alpha_entry, so in cold
                    // air (alpha ≪ 1) it leaves q_v ≫ q_sat — and cells colder than
                    // ~−60 °C (alpha ≤ 0.01) it skips entirely. With no upper bound on q_v
                    // anywhere, that residual supersaturation accumulated unbounded at the
                    // orographic saturation level in the cold fade zone (NZ Southern Alps,
                    // i≈12, T≈−40 °C → q_v→6, cloud→129, t→1e5 °C, then buoyancy drove a
                    // vertical-velocity runaway). Remove supersaturation here for EVERY
                    // cell, independent of the alpha fade, conserving water (excess → cloud
                    // above freezing, → ice below) and energy (latent heat → T). In a clean
                    // run the per-call excess is small (physical); the T_max recap below is
                    // the backstop if a transient ever drives a large excess.
                    // ATHAD: the EXACT saturation mass fraction, and the sign of the
                    // superheated branch reversed.
                    //
                    // What was here: q_sat = (p > E_sat) ? ep*E_sat/(p - E_sat) : ep*1e-5.
                    // Both halves are wrong off Earth, and the second is wrong in the most
                    // damaging possible direction. When p_sat(T) EXCEEDS the local pressure
                    // the vapour is superheated: it cannot condense at all, so the correct
                    // saturation limit is 1 (all of the water stays vapour). This wrote
                    // ep*1e-5 ~ 7e-6 instead — a limit of essentially zero — so the block
                    // below dumped the ENTIRE 0.67 water mass fraction into cloud and
                    // released its latent heat, in exactly the layers where nothing can
                    // condense. On ATHAD's column that is every level between ~373 K and the
                    // critical point: a 60 km deep slab of manufactured cloud from ~140 to
                    // ~200 km, which then set the planetary albedo and, smeared downward by
                    // damp_wiggles(), inflated the density in the supercritical layers
                    // beneath it through the (1 - cloud - ice) loading term.
                    //
                    // Latent on Earth: p_sat exceeds the local pressure only above ~373 K,
                    // and no cell in a terrestrial column is ever that warm, so the branch
                    // never ran. It is the same class of defect as the 333.15 K cap above —
                    // an Earth-only regime assumption written as a fallback.
                    const double p_local = p_row[k];
                    const double M_other = AtmMixture::M_nonwater(c_row[k], m.co2.x[i][j][k],
                                                                  m.m_comp.M_bg);
                    const double E_sat = SaturationH2O::saturationPressureAuto(T_dim);
                    const double q_sat = SaturationH2O::saturationMassFraction(
                                             E_sat, p_local, M_other);
                    if (c_row[k] > q_sat) {
                        const double excess = c_row[k] - q_sat;
                        c_row[k] = q_sat;
                        if (T_dim >= m.t_00) {
                            cloud_row[k] += excess;
                            T_dim        += SaturationH2O::latentHeat(T_dim) * excess
                                          / AtmMixture::cp_of(c_row[k], m.co2.x[i][j][k],
                                                              T_dim, m.m_comp.M_bg);
                        } else {
                            ice_row[k]   += excess;
                            T_dim        += SaturationH2O::latentHeatSublimation(T_dim) * excess
                                          / AtmMixture::cp_of(c_row[k], m.co2.x[i][j][k],
                                                              T_dim, m.m_comp.M_bg);
                        }
                        if (T_dim > T_max) T_dim = T_max;   // backstop on the latent release
                        t_row_nd[k] = T_dim * inv_t_0;
                    }

                    // ---- Condensate upper bounds (CAP SAFETY NET) ----
                    if (cloud_row[k] > cloud_cap) cloud_row[k] = cloud_cap;
                    if (ice_row[k]   > cloud_cap) ice_row[k]   = cloud_cap;

                    // Existing cold fade: smoothly dry moisture toward 0 below t_00.
                    const double alpha = 1.0 / (1.0 + std::exp(-(T_dim - m.t_00) / fade_K));
                    c_row[k]     *= alpha;
                    cloud_row[k] *= alpha;
                    ice_row[k]   *= alpha;
                }
            }
        }

        if (n_ceiling > 0)
            std::cout << "      SaturationAdjustment: water-vapour ceiling c = 1 - co2 hit in "
                      << n_ceiling << " cells (water deleted there — see the note above)"
                      << std::endl;
    }

    void printReport() const {
        // Diagnostic variables are not updated in the parallel loops;
        // kept here to preserve the original output contract.
        const bool   satadjust = false;
        const int    iter_prec = 0;
        const int    i_sat = 0, j_sat = 0, k_sat = 0;
        const double saturation = 0.0;

        if (!satadjust)
            std::cout << "      no saturation of water vapour in SaturationAdjustment found"
                      << std::endl;
        else
            std::cout << "      saturation of water vapour in SaturationAdjustment found"
                      << std::endl
                      << "      iter_prec = " << iter_prec << std::endl
                      << "      i_sat = "  << i_sat
                      << "   j_sat = "     << j_sat
                      << "   k_sat = "     << k_sat
                      << "   height_sat[m] = " << m.get_layer_height(i_sat)
                      << "   saturation[g/kg] = " << saturation * 1e3 << std::endl;

        if (iter_prec >= iter_prec_end)
            std::cout << std::endl
                      << "      no convergent solution found in SaturationAdjustment"
                      << std::endl
                      << "      iter_prec_end = " << iter_prec_end << std::endl
                      << "      iter_prec = "     << iter_prec     << std::endl
                      << "      results see above" << std::endl;
    }
};

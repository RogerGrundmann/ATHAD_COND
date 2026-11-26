/*
 * Dry convective adjustment — Manabe & Strickler (1964).
 *
 * PORTED 2026-11-26 to ATHAD_PERID / ATHAD_COND (user: "port ATM_PROGNOSTIC_T to ATHAD_PERID and COND",
 * "With the moist option"), copied from ATHAD's working tree including ATM_CONV_ADJ_MOIST. Here it runs
 * ONLY under ATM_PROGNOSTIC_T=1 (cAtmosphereModel.cpp, before densities() in the time loop), so a run
 * without the knob never calls it and stays bit-identical. Invariant 1 (no topography) holds here as in
 * ATHAD, so i0 = 0 is right; SaturationH2O::moistLapse is this tree's own (ATHAD_COND is its origin).
 *
 * Ported from the family's shared planet/ConvectiveAdjustment.h (ATURAN/ATJUP/ATSAT), whose
 * header comment is the reference for why a model needs this at all and is worth reading. The
 * two properties that matter: it removes the instability COMPLETELY rather than damping it, and
 * it does not create or destroy energy.
 *
 * WHY ATHAD NEEDS IT, SPECIFICALLY. Until now ThermoAtm::densities() overwrote t with its own
 * adiabat every iteration, so the column could not go unstable — the prescription was doing the
 * convective adjustment's job, invisibly. README item 7 measured what that costs: with
 * ATM_PROGNOSTIC_T=1 the energy imbalance falls from -433 to -19 W/m2, but the lower column then
 * cools at 5.97 K/km against a 4.19 K/km dry adiabat (superadiabatic, i.e. convectively
 * unstable) and the top 50 km develops a grid-scale sawtooth. Radiative equilibrium in an
 * optically thick atmosphere IS superadiabatic; nothing else in this model restores a column to
 * its adiabat. This file is the prerequisite for finishing item 7.
 *
 * WHAT WAS NOT COPIED, AND WHY. The shared file is deliberately byte-identical across the
 * planets that use it, and two of its assumptions do not hold here:
 *
 *   1. IT ASSUMES A UNIFORM GRID. `dz_m = L_atm*1e3/(im-1)` is a layer thickness only when the
 *      layers are equal. ATHAD's radial coordinate is exponentially stretched — dz runs from
 *      0.81 km at the surface to ~15 km at the lid, a factor of 19 — so a single critical drop
 *      would be ~19x too strict at the bottom and far too lax at the top. This is the same
 *      defect class as init_tropopause_layers' round(h/L_atm) (CLAUDE.md): a height-to-index
 *      conversion that is only valid on a uniform mesh. The critical drop is therefore computed
 *      PER LAYER from get_layer_height(), and the segment adiabat becomes a cumulative sum
 *      rather than a linear ramp.
 *
 *   2. IT ASSUMES A CONSTANT cp_mix. ATHAD's cp varies by a factor of two across 300-1500 K and
 *      follows the composition, so cp is taken locally from AtmMixture::cp_of() at each layer —
 *      which is strictly better than the shared file's constant and is the reason its own
 *      comment warns not to substitute a textbook value.
 *
 * WHAT IT DELIBERATELY DOES NOT DO, carried over unchanged: it mixes temperature only, not
 * composition, and it uses the DRY adiabat even where cloud is present. A moist adjustment would
 * use the saturated lapse — ATHAD_COND has SaturationH2O::moistLapse for exactly that — but that
 * is a modelling decision with consequences for the microphysics already running in the
 * saturation adjustment, and it is not made here.
 *
 * ATM_CONV_ADJ_MOIST (default 0, 2026-10-10) is that moist version, as an option: the critical
 * drop of a layer is the SATURATED one (SaturationH2O::moistLapse, ported from ATHAD_COND) where
 * the lower cell of the layer is saturated (c >= 0.99 q_sat, the test of ATHAD_COND's
 * densities()), and the dry one elsewhere. It still mixes temperature only. The saturated lapse
 * is the SMALLER one, so the option makes the scheme stricter in saturated layers and changes
 * nothing in the supercritical column, where moistLapse returns g/cp. The critical drops are
 * computed once per column from the state before the sweep and are not re-evaluated while it
 * mixes. With the prescribed profile (ATM_PROGNOSTIC_T=0) densities() overwrites the result,
 * as it does the dry one's.
 */

#pragma once

#include "cAtmosphereModel.h"
#include "MixtureAtm.h"
#include "SaturationH2O.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

class ConvectiveAdjustment {
public:
    explicit ConvectiveAdjustment(cAtmosphereModel& model) : m(model) {}

    void run()
    {
        using namespace std;
        // ATM_CONV_ADJ_MOIST (default 0): the saturated lapse where the cell is saturated.
        static const bool moist = [](){ const char* e = getenv("ATM_CONV_ADJ_MOIST");
                                        return e && atoi(e) != 0; }();
        cout << endl << "      AGCM: ConvectiveAdjustment (" << (moist ? "moist where saturated" : "dry")
             << ", Manabe-Strickler)" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        // ATM_CONV_ADJ_LAPSE scales the critical lapse rate: 0 gives an isothermal criterion,
        // values above 1 make the scheme stricter than the dry adiabat, which is one crude way
        // to stand in for a moist adiabat in a condensing region.
        static const double lapse_fac = [](){ const char* e = getenv("ATM_CONV_ADJ_LAPSE");
                                              return e ? atof(e) : 1.0; }();
        // Sweeps repeat until the column is stable; the cap only exists so a pathological
        // column cannot spin here.
        static const int max_pass = [](){ const char* e = getenv("ATM_CONV_ADJ_PASSES");
                                          const int v = e ? atoi(e) : 64;
                                          return v > 0 ? v : 64; }();

        // A column is left alone unless it is superadiabatic by more than this, so round-off
        // does not make the scheme fire on a column that is already neutral.
        constexpr double tol_nd = 1.0e-12;

        // Layer thicknesses of the stretched grid, once.
        std::vector<double> dz(m.im, 0.0);
        for (int i = 0; i + 1 < m.im; i++)
            dz[i] = m.get_layer_height(i+1) - m.get_layer_height(i);

        long long n_columns_adjusted = 0, n_layers_mixed = 0, n_moist_layers = 0;
        int    max_passes_used = 0;
        double max_dT_K = 0.0, max_rel_drift = 0.0;

        std::vector<double> w(m.im), t_col(m.im), dT_ad(m.im);

        #pragma omp parallel for collapse(2) schedule(static) firstprivate(w, t_col, dT_ad) \
            reduction(+:n_columns_adjusted, n_layers_mixed, n_moist_layers) \
            reduction(max:max_passes_used, max_dT_K, max_rel_drift)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                // No topography in ATHAD (invariant 1), so the fluid column is the whole
                // column. Kept as i0/i1 rather than hard-coded so the routine still reads
                // correctly if a sibling with topography picks it up.
                const int i0 = 0, i1 = m.im - 1;
                if (i1 - i0 < 1) continue;

                // Mass weights and the per-layer critical drop.
                //
                // The enthalpy per unit area of a layer is (cp/g)*dp, so the conserved quantity
                // is the dp-weighted temperature. dp comes from the hydrostatic p_stat with
                // faces at the midpoints; the end layers get their half-cell. Using dp rather
                // than a density avoids a trap: rho = p/(R*T) makes rho*T identically p/R, so a
                // density weight would conserve nothing at all.
                bool ok = true;
                for (int i = i0; i <= i1; i++) {
                    const double p_lo = (i > i0) ? 0.5 * (m.p_stat.x[i-1][j][k] + m.p_stat.x[i][j][k])
                                                 : m.p_stat.x[i0][j][k];
                    const double p_hi = (i < i1) ? 0.5 * (m.p_stat.x[i][j][k] + m.p_stat.x[i+1][j][k])
                                                 : m.p_stat.x[i1][j][k];
                    w[i]     = p_lo - p_hi;
                    t_col[i] = m.t.x[i][j][k];
                    if (!(w[i] > 0.0) || !std::isfinite(w[i]))     ok = false;
                    if (!std::isfinite(t_col[i]) || t_col[i] <= 0.0) ok = false;

                    // Critical drop across layer i -> i+1, non-dimensional, on the LOCAL cp.
                    if (i < i1) {
                        const double cp_loc = AtmMixture::cp_of(m.c.x[i][j][k], m.co2.x[i][j][k],
                                                                t_col[i] * m.t_0, m.m_comp.M_bg);
                        double gamma = (cp_loc > 0.0) ? m.g / cp_loc : 0.0;
                        if (moist && cp_loc > 0.0) {
                            const double T_dim = t_col[i] * m.t_0;
                            const double M_nw  = AtmMixture::M_nonwater(m.c.x[i][j][k],
                                                     m.co2.x[i][j][k], m.m_comp.M_bg);
                            const double q_s   = SaturationH2O::saturationMassFractionAt(
                                                     T_dim, m.p_stat.x[i][j][k], M_nw);
                            if (q_s < 1.0 && m.c.x[i][j][k] >= 0.99 * q_s) {
                                const double q_l   = m.cloud.x[i][j][k] + m.ice.x[i][j][k]
                                                   + m.gr.x[i][j][k];
                                const double R_loc = AtmMixture::R_of(m.c.x[i][j][k],
                                                         m.co2.x[i][j][k], m.m_comp.R_bg, q_l);
                                gamma = SaturationH2O::moistLapse(T_dim, m.p_stat.x[i][j][k],
                                                                  M_nw, cp_loc, R_loc, m.g);
                                n_moist_layers++;
                            }
                        }
                        dT_ad[i] = lapse_fac * gamma * dz[i] / m.t_0;
                    }
                }
                // A column with a non-monotonic p_stat or a bad temperature is left untouched:
                // this routine is not the place to repair either, and mixing across a NaN would
                // spread it through the whole column.
                if (!ok) continue;

                double sum_before = 0.0;
                for (int i = i0; i <= i1; i++) sum_before += w[i] * t_col[i];

                // Whole unstable SEGMENTS are mixed at once, not adjacent pairs: pairwise mixing
                // moves heat one layer per sweep, so a deep unstable block needs as many sweeps
                // as it has layers. A segment [a..b] is put on the adiabat T_q = C - D_q, where
                // D_q is the CUMULATIVE adiabatic drop from a to q (the shared file's linear
                // dT_ad*(q-a), generalised to a stretched grid), with C fixed by conserving the
                // dp-weighted temperature:  C = [sum w_q T_q + sum w_q D_q] / sum w_q.
                int passes = 0;
                long long layers_mixed = 0;
                bool changed = true;
                std::vector<double> D;
                while (changed && passes < max_pass) {
                    changed = false;
                    passes++;
                    int i = i0;
                    while (i < i1) {
                        if (t_col[i] - t_col[i+1] <= dT_ad[i] + tol_nd) { i++; continue; }

                        int a = i, b = i + 1;
                        for (;;) {
                            D.assign(b - a + 1, 0.0);
                            for (int q = a + 1; q <= b; q++)
                                D[q - a] = D[q - a - 1] + dT_ad[q - 1];

                            double sw = 0.0, swt = 0.0, swd = 0.0;
                            for (int q = a; q <= b; q++) {
                                sw  += w[q];
                                swt += w[q] * t_col[q];
                                swd += w[q] * D[q - a];
                            }
                            const double C = (swt + swd) / sw;
                            for (int q = a; q <= b; q++) {
                                const double t_new = C - D[q - a];
                                const double dK = std::fabs(t_new - t_col[q]) * m.t_0;
                                if (dK > max_dT_K) max_dT_K = dK;
                                t_col[q] = t_new;
                            }

                            // Mixing can destabilise the joint with the layer below or above, so
                            // the segment grows in whichever direction is still too steep and is
                            // re-mixed. It can only grow, and only within the column, so this
                            // terminates.
                            bool extended = false;
                            if (a > i0 && t_col[a-1] - t_col[a] > dT_ad[a-1] + tol_nd) { a--; extended = true; }
                            if (b < i1 && t_col[b] - t_col[b+1] > dT_ad[b]   + tol_nd) { b++; extended = true; }
                            if (!extended) break;
                        }

                        layers_mixed += (b - a + 1);
                        changed = true;
                        i = b;                          // carry on above the block just mixed
                    }
                }

                if (layers_mixed == 0) continue;

                double sum_after = 0.0;
                for (int i = i0; i <= i1; i++) sum_after += w[i] * t_col[i];
                const double drift = (sum_before != 0.0)
                                   ? std::fabs(sum_after - sum_before) / std::fabs(sum_before) : 0.0;
                if (drift > max_rel_drift) max_rel_drift = drift;

                for (int i = i0; i <= i1; i++) m.t.x[i][j][k] = t_col[i];

                n_columns_adjusted++;
                n_layers_mixed += layers_mixed;
                if (passes > max_passes_used) max_passes_used = passes;
            }
        }

        const double frac = 100.0 * double(n_columns_adjusted) / double(m.jm * m.km);
        printf("      AGCM: convective adjustment — %lld of %d columns (%.2f %%), %lld layers,"
               " worst column %d sweeps of %d, max dT %.3f K, enthalpy drift %.2e\n",
               n_columns_adjusted, m.jm * m.km, frac, n_layers_mixed,
               max_passes_used, max_pass, max_dT_K, max_rel_drift);
        if (moist)
            printf("      AGCM: convective adjustment — %lld layers took the saturated lapse\n",
                   n_moist_layers);
        if (max_passes_used >= max_pass)
            cout << "      AGCM: WARNING - the sweep cap was reached, a column may still be "
                    "superadiabatic (raise ATM_CONV_ADJ_PASSES)" << endl;

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for ConvectiveAdjustment\n", elapsed.count() * 1e-9);
        cout << "      AGCM: ConvectiveAdjustment ended" << endl;
    }

private:
    cAtmosphereModel& m;
};

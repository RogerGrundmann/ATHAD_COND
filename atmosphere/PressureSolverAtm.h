#pragma once

#include "cAtmosphereModel.h"
#include "Utils.h"

#include <vector>
#include <cstdint>
#include <cstdlib>   // getenv/atof for the ATM_POISSON_METRIC_FIX A/B knob
#include <cmath>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AtomUtils;

class PressureSolverAtm {
public:
    explicit PressureSolverAtm(cAtmosphereModel& model)
        : m(model)
    {}

    void run(bool verbose = true)
    {
        using namespace std;
        if (verbose) cout << endl << endl << endl << "      ATOM: PressureSolverAtm" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        // precompute sin(the) table — only depends on j, avoids redundant sin() calls
        std::vector<double> sinthe_table(m.jm);
        for (int j = 0; j < m.jm; j++) {
            sinthe_table[j] = sin(m.the.z[j]);
            if (sinthe_table[j] < 0.55) sinthe_table[j] = 0.55;   // metric floor ~57° — keep in sync with the RHS geometry floor (RungeKutta_Atm*.cpp)
        }

        // Precompute land mask — eliminates repeated function call overhead
        // Allocate flat mask: 1 = land, 0 = air. Also enforce no-penetration on the
        // provisional velocity: bcSolidGround masks u/v/w at land but NOT aux_u/v/w, and
        // RHS_Atm writes aux unconditionally — so land cells carry stale aux (advection of
        // the adjacent air jet + the blown-up p_dyn gradient). The fluid-cell divergence
        // source uses central differences at cliff edges, differencing against that stale
        // land aux. Zero it here so every sweep sees a clean solid wall.
        std::vector<int8_t> land(m.im * m.jm * m.km);
        #pragma omp parallel for collapse(3)
        for (int i = 0; i < m.im; i++)
            for (int j = 0; j < m.jm; j++)
                for (int k = 0; k < m.km; k++) {
                    const bool isl = is_land(m.h, i, j, k);
                    land[i*m.jm*m.km + j*m.km + k] = isl ? 1 : 0;
                    if (isl) {
                        m.aux_u.x[i][j][k] = 0.0;
                        m.aux_v.x[i][j][k] = 0.0;
                        m.aux_w.x[i][j][k] = 0.0;
                    }
                }

        #define LAND(i,j,k) land[(i)*m.jm*m.km + (j)*m.km + (k)]

        // ==================================================================
        // ANELASTIC CONTINUITY (ATM_ANELASTIC, default 0 = bit-identical).
        //
        // Steps 2-7 of the anelastic scope in the README. The projection below enforces
        // div(u) = 0; a column whose density spans five orders of magnitude needs
        // div(rho_bar u) = 0. Written as one knob because the pieces are not separable:
        // a divergence source without the matching Poisson operator, or either without the
        // wall BC, is a projection onto the wrong subspace and measures nothing.
        //
        //   source     D = (1/rho_bar) div(rho_bar u*) = div(u*) + u*_r dln(rho_bar)/dr
        //   operator   (1/rho_bar) div(rho_bar grad p) = lap(p) + dln(rho_bar)/dr dp/dr
        //   BC         rho_bar u_r = 0 at surface and lid, i.e. aux_u = 0 there
        //
        // Both extra terms carry the SAME dln(rho_bar)/dr array, which is what keeps the
        // discrete divergence and gradient adjoint to each other — the lesson of step 1
        // (ATM_POISSON_METRIC_FIX), where an inconsistent pair meant the projection could
        // not remove what it had measured. The operator gains a first-derivative term and
        // no new solver: the 7-point stencil keeps its shape, the radial off-diagonals
        // become num1 +/- A. Diagonal dominance survives because A/num1 = |dln rho|*dr/2
        // and rho_bar changes by ~0.19 in the log per level here, so the ratio is ~0.1.
        //
        // rho_bar is 1-D and rebuilt by ThermoAtm::densities(); before the first call the
        // vectors are empty and this falls back to Boussinesq automatically.
        //
        // Step 7 (the buoyancy reference) needed no change: computeLevelMeanTemperature()
        // already builds t_ref_level[i] as the sin(colatitude)-weighted horizontal mean of
        // t at each level, which is the same average, of the same prescribed profile, that
        // rho_bar is. The two base states agree by construction.
        static const bool anelastic_knob = [](){ const char* e = getenv("ATM_ANELASTIC");
                                                 return e ? (atof(e) != 0.0) : false; }();
        const bool anelastic = anelastic_knob
                            && (int)m.m_dlnrho_dr.size() == m.im;
        const double* const dlnrho = anelastic ? m.m_dlnrho_dr.data() : nullptr;

        // Fuse the three boundary loops into one pass
        //
        // aux_u at the two radial walls: c43/c13 extrapolation is zero-GRADIENT, which
        // permits a through-wall mass flux and leaves the Poisson problem's compatibility
        // condition to be satisfied by whatever the extrapolation happens to produce. Under
        // ATM_ANELASTIC it becomes zero normal mass flux (step 5) — rho_bar > 0 everywhere,
        // so rho_bar u_r = 0 is u_r = 0. Note this is the PROVISIONAL velocity only; the
        // prognostic u keeps bcRadius's zero-gradient wall, which is a separate question.
        #pragma omp parallel for collapse(2)
        for (int j = 1; j < m.jm-1; j++) {
            for (int k = 1; k < m.km-1; k++) {
                if (anelastic) {
                    m.aux_u.x[0][j][k]      = 0.0;
                    m.aux_u.x[m.im-1][j][k] = 0.0;
                } else {
                    m.aux_u.x[0][j][k]      = m.c43 * m.aux_u.x[1][j][k]      - m.c13 * m.aux_u.x[2][j][k];
                    m.aux_u.x[m.im-1][j][k] = m.c43 * m.aux_u.x[m.im-2][j][k] - m.c13 * m.aux_u.x[m.im-3][j][k];
                }
                m.aux_v.x[0][j][k]      = m.c43 * m.aux_v.x[1][j][k]      - m.c13 * m.aux_v.x[2][j][k];
                m.aux_v.x[m.im-1][j][k] = m.c43 * m.aux_v.x[m.im-2][j][k] - m.c13 * m.aux_v.x[m.im-3][j][k];
                m.aux_w.x[0][j][k]      = m.c43 * m.aux_w.x[1][j][k]      - m.c13 * m.aux_w.x[2][j][k];
                m.aux_w.x[m.im-1][j][k] = m.c43 * m.aux_w.x[m.im-2][j][k] - m.c13 * m.aux_w.x[m.im-3][j][k];
            }
        }

        // Grid-spacing reciprocals — constant for the entire grid
        const double inv_2dr   = 1.0 / (2.0 * m.dr);
        const double inv_2dthe = 1.0 / (2.0 * m.dthe);
        const double inv_2dphi = 1.0 / (2.0 * m.dphi);
        const double inv_dr2   = 1.0 / (m.dr   * m.dr);
        const double inv_dthe2 = 1.0 / (m.dthe * m.dthe);
        const double inv_dphi2 = 1.0 / (m.dphi * m.dphi);
        const double inv_dthe  = 1.0 / m.dthe;
        const double inv_dphi  = 1.0 / m.dphi;

        // Cap on the divergence source of the pressure Poisson update. Over steep
        // orography (cliff faces: Patagonian Andes 49°S/74°W, Atlas 36°N/2°E) the
        // central-difference divergence of the provisional velocity can spike, inject a
        // huge Poisson source → runaway p_dyn → pressure-gradient force → velocity that
        // pegs the ±100 m/s clamp → larger divergence (a closed dry instability,
        // independent of moist physics). Bounding the per-cell source to denom·p_dyn_cap
        // caps its contribution to p_dyn at p_dyn_cap (discrete max principle). p_dyn is
        // non-dimensionalised by p_0=1013.25 hPa, so healthy dynamic pressure is ≪1
        // non-dim (tens of hPa); the runaway reached 18–43 non-dim (≈18000–44000 hPa).
        // p_dyn_cap is a GENEROUS backstop only. A tight value (1.0) under-removed the
        // velocity divergence (the source IS the divergence the projection must cancel),
        // breaking incompressibility and triggering a worse dry CFL blow-up at the NZ Alps
        // (46°S/170°E, iter 288). 2.0 leaves the projection intact while still catching the
        // gross Andes injection; the orographic Shapiro filter (post-solve, on u/v/w) is now
        // the primary stabiliser for the steep-orography CFL mode, not this clamp.
        constexpr double p_dyn_cap = 2.0;

        // ⚠️ A/B KNOB 2026-07-21 (ATM_POISSON_METRIC_FIX, default 0 = bit-identical). The θ/φ
        // Poisson Laplacian coefficients (num2/num3/denom) use SINGLE-power inv_rm /
        // inv_rmsinthe, which is inconsistent with this solver's own divergence source and
        // gradient correction (both single-power). The consistent discrete div·grad needs
        // inv_rm² (θ) and inv_rm²/sin²θ (φ) — exactly the pattern the RADIAL term already
        // uses (exp_2_rm in the Laplacian vs single exp_rm in the divergence). The ocean
        // diagnosed this single-power form as the non-idempotent collocated projection and
        // moved its per-iter solve to project_velocity with inv_rm2 / inv_rm2sinthe2
        // (PressureSolverHyd.h:428). This knob brings the atm θ/φ Laplacian to the same
        // consistent metric so it can be A/B tested against the tuned steep-orography
        // stabilisers (p_dyn_cap, p_dyn_ceiling, topo Dirichlet pins) that were calibrated on
        // the old operator. NB: this repairs the metric POWER only; the collocated checkerboard
        // (Rhie-Chow face reconstruction) is a separate, larger port not done here.
        // ATHAD: ON by default since the A/B measured it. With the consistent metric the
        // residual divergence of the projected velocity falls from rms 2.625e-02 to
        // 2.153e-02 — an 18 % improvement — while the water-mass drift is bit-identical
        // (+0.1698 % over 20 iterations either way), so this is a genuine repair of the
        // discrete div/grad adjointness and NOT the cause of the tracer mass error. The
        // env var still forces it off for A/B (ATM_POISSON_METRIC_FIX=0).
        //
        // Caveat carried from the note above: p_dyn_cap, the p_dyn_ceiling and the
        // topography Dirichlet pins were calibrated against the OLD operator. Watch a long
        // run for instability; re-deriving them is step 6 of the anelastic scope.
        const bool poisson_metric_fix = [](){ const char* e = getenv("ATM_POISSON_METRIC_FIX");
                                              return e ? (atof(e) != 0.0) : true; }();

        // ATOM_METRIC_DIVERGENCE — hoisted out of the cell loop; see lib/Utils.h.
        const bool metric_div = AtomUtils::metric_divergence();

        // How often the p_dyn_cap source clamp actually binds. Step 6 of the anelastic
        // scope says the stabilisers were calibrated against the old operator and may clip
        // the anelastic projection before it acts; this is the number that says whether
        // they do, rather than leaving it to be assumed either way.
        long n_src_clamped = 0, n_src_cells = 0;

        // ==================================================================
        // THE RELAXATION IS RED-BLACK, AND WAS NOT ALWAYS.
        //
        // Each solve is two passes over a checkerboard colouring of (i+j+k): every cell of
        // one colour has all six of its stencil neighbours in the other, so within a pass
        // nothing is read while it is being written.
        //
        // This replaced `#pragma omp parallel for collapse(2) schedule(dynamic, 4)` over
        // (i,j) writing p_dyn IN PLACE while reading p_dyn[i±1][j±1] — over the very two
        // indices the stencil reads across. Cell (i,j,k) was read by the thread owning
        // (i+1,j) or (i,j+1) while its owner was writing it, and schedule(dynamic) made it
        // worse than a thread-count dependence: which thread got which chunk varied with
        // timing, so the SAME binary at the SAME thread count gave different answers RUN TO
        // RUN. Measured here before the fix, two runs of one binary at 24 threads, one
        // iteration: max u-component 0.080795 against 0.080751, the extremum wandering in
        // longitude and hemisphere. The varying cell sat at 37905 m, i ~ 39 of 61 —
        // interior, which is where an in-place stencil race puts it and not where a
        // boundary defect would.
        //
        // Red-black rather than Jacobi because of what the old loop was reaching for: k runs
        // serially inside a thread, so k-1 is current and k+1 one sweep old, i.e.
        // lexicographic Gauss-Seidel — correct in serial, broken only by the (i,j)
        // parallelism. Jacobi would have been the easier fix and would have cost the
        // convergence rate. The colour is selected with a `continue` rather than by striding
        // k, because the k loop carries a sliding window over the land mask that assumes
        // consecutive k.
        //
        // Ported from ATHAD, which took it from the family's shared PressureSolver.h
        // (ATURAN `ffd0e0e` diagnosed the same defect and cured it by serialising, which was
        // right there — its computePressure is 0.003 s of a 5.3 s step; here the projection
        // is far too expensive for that).
        //
        // NOTE this does NOT reproduce the old 1-thread answer: red-black is a different
        // sweep order from lexicographic Gauss-Seidel, so it converges to the same solution
        // by a different path. Every measurement in this repo taken before this commit moves
        // in its last digits.
        //
        // The two colours together visit each cell exactly once, so the clamp counters below
        // are not double-counted.
        // ==================================================================
        for (int colour = 0; colour < 2; colour++) {

        // Main compute loop — land mask lookups + hoisted j-invariants + k sliding window
        #pragma omp parallel for collapse(2) schedule(static) \
                reduction(+:n_src_clamped,n_src_cells)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {

                // Build geometry struct ONCE per (i,j)
                cAtmosphereModel::CellGeometry geo;

                geo.rm           = m.rad.z[i];
                geo.rm2          = geo.rm * geo.rm;
                geo.exp_rm       = 1.0 / (geo.rm + 1.0);
                geo.exp_2_rm     = geo.exp_rm * geo.exp_rm;
                geo.sinthe       = sinthe_table[j];
                geo.sinthe2      = geo.sinthe * geo.sinthe;
                geo.costhe       = cos(m.the.z[j]);
                // ATM_METRIC_RADIUS — must match RungeKutta_Atm_Turb exactly; identity when off.
                const double rmet  = m.metricRadius(geo.rm);
                const double rmet2 = rmet * rmet;
                geo.inv_rm       = 1.0 / rmet;
                geo.inv_rm2      = 1.0 / rmet2;
                geo.inv_rmsinthe         = 1.0 / (rmet * geo.sinthe);
                geo.inv_rm2sinthe        = geo.inv_rm2 / geo.sinthe;
                geo.inv_rm2sinthe2       = geo.inv_rm2 / geo.sinthe2;
                geo.costhe_inv_rm2sinthe = geo.costhe * geo.inv_rm2sinthe;
                geo.inv_2dr   = inv_2dr;
                geo.inv_2dthe = inv_2dthe;
                geo.inv_2dphi = inv_2dphi;
                geo.inv_dr2   = inv_dr2;
                geo.inv_dthe2 = inv_dthe2;
                geo.inv_dphi2 = inv_dphi2;

                // θ/φ Laplacian metric: single-power (legacy) or the consistent double-power
                // (inv_rm2, inv_rm2sinthe2) when ATM_POISSON_METRIC_FIX is set. See knob note above.
                const double m_the = poisson_metric_fix ? geo.inv_rm2        : geo.inv_rm;
                const double m_phi = poisson_metric_fix ? geo.inv_rm2sinthe2 : geo.inv_rmsinthe;
                const double denom = 2.0 * geo.exp_2_rm * inv_dr2
                                   + 2.0 * m_the * inv_dthe2
                                   + 2.0 * m_phi * inv_dphi2;
                const double inv_denom = 1.0 / denom;
                const double num1 = geo.exp_2_rm * inv_dr2;
                const double num2 = m_the * inv_dthe2;
                const double num3 = m_phi * inv_dphi2;

                // Anelastic first-derivative term of the operator (step 4). Both radial
                // factors are exp_rm — dln(rho_bar)/dr_physical = exp_rm * dln(rho_bar)/d(rad.z)
                // and dp/dr_physical = exp_rm * dp/d(rad.z) — hence exp_2_rm, matching the
                // Laplacian's radial term. Zero when the knob is off, so denom, num1..num3
                // and the update below are untouched.
                const double num_a = anelastic
                                   ? geo.exp_2_rm * dlnrho[i] * inv_2dr : 0.0;

                const bool i_in_range = (i < m.im-2);
                const bool j_inner    = (j > 2) && (j < m.jm-2);

                // sliding window for k-direction land/air status
                int8_t lnd_k0 = LAND(i,j,0), lnd_k1 = LAND(i,j,1);

                for (int k = 1; k < m.km-1; k++) {
                    const int8_t lnd_ijk = lnd_k1;
                    const int8_t lnd_kp1 = LAND(i,j,k+1);
                    const int8_t lnd_km1 = lnd_k0;
                    const int8_t lnd_kp2 = (k < m.km-2) ? LAND(i,j,k+2) : 0;
                    const int8_t lnd_km2 = (k > 2)      ? LAND(i,j,k-2) : 0;

                    lnd_k0 = lnd_k1;
                    lnd_k1 = lnd_kp1;

                    // Cells of the other colour are skipped AFTER the window bookkeeping
                    // above — lnd_k0/lnd_k1 slide with k and assume every k is visited.
                    if (((i + j + k) & 1) != colour) continue;

                    double du_dr, dv_dthe, dw_dphi;
                    bool r_flag   = false;
                    bool the_flag = false;
                    bool phi_flag = false;

                    // EXPERIMENT 2 (2026-06-01): fluid-cell one-sided differences AWAY from
                    // a land neighbour. Without this, the central-difference fallback below
                    // samples aux_u/v/w at the zeroed land cell (line 54), producing a
                    // ~3× overestimate of div_src at coastal fluid cells → spuriously large
                    // negative p_dyn → −∂p/∂r drives u upward → divergence grows → loop.
                    // Use one-sided 2nd-order forward/backward stencils that don't touch
                    // the wall value at all.
                    //
                    // r direction (fluid above land or below overhang)
                    if (!lnd_ijk && i >= 1 && i+2 < m.im) {
                        if (LAND(i-1,j,k)) {
                            du_dr = (-3.0 * m.aux_u.x[i][j][k] + 4.0 * m.aux_u.x[i+1][j][k]
                                     - m.aux_u.x[i+2][j][k]) * inv_2dr;
                            r_flag = true;
                        } else if (i-2 >= 0 && LAND(i+1,j,k)) {
                            du_dr = (3.0 * m.aux_u.x[i][j][k] - 4.0 * m.aux_u.x[i-1][j][k]
                                     + m.aux_u.x[i-2][j][k]) * inv_2dr;
                            r_flag = true;
                        }
                    }

                    // r direction (land cell with air above — existing land-stencil case;
                    // value gets discarded by the Neumann mean update at line 241)
                    if (!r_flag && i_in_range && lnd_ijk && !LAND(i+1,j,k)) {
                        du_dr = (-3.0 * m.aux_u.x[i][j][k] + 4.0 * m.aux_u.x[i+1][j][k]
                                 - m.aux_u.x[i+2][j][k]) * inv_2dr;
                        r_flag = true;
                    }

                    // theta direction
                    if (j_inner) {
                        const int8_t air_jp1 = !LAND(i,j+1,k);
                        const int8_t air_jm1 = !LAND(i,j-1,k);

                        // Fluid cell with land neighbour in theta
                        if (!the_flag && !lnd_ijk) {
                            if (!air_jm1 && air_jp1 && !LAND(i,j+2,k)) {
                                dv_dthe = (-3.0 * m.aux_v.x[i][j][k] + 4.0 * m.aux_v.x[i][j+1][k]
                                           - m.aux_v.x[i][j+2][k]) * inv_2dthe;
                                the_flag = true;
                            } else if (!air_jp1 && air_jm1 && !LAND(i,j-2,k)) {
                                dv_dthe = (3.0 * m.aux_v.x[i][j][k] - 4.0 * m.aux_v.x[i][j-1][k]
                                           + m.aux_v.x[i][j-2][k]) * inv_2dthe;
                                the_flag = true;
                            }
                        }

                        if (lnd_ijk && air_jp1 && !LAND(i,j+2,k)) {
                            dv_dthe = (-3.0 * m.aux_v.x[i][j][k] + 4.0 * m.aux_v.x[i][j+1][k]
                                       - m.aux_v.x[i][j+2][k]) * inv_2dthe;
                            the_flag = true;
                        }
                        if (lnd_ijk && air_jm1 && !LAND(i,j-2,k)) {
                            dv_dthe = -(-3.0 * m.aux_v.x[i][j][k] + 4.0 * m.aux_v.x[i][j-1][k]
                                        - m.aux_v.x[i][j-2][k]) * inv_2dthe;
                            the_flag = true;
                        }
                        if ((lnd_ijk && air_jp1 && LAND(i,j+2,k))
                            || (j == m.jm-2 && !lnd_ijk && LAND(i,j+1,k))) {
                            dv_dthe = (m.aux_v.x[i][j+1][k] - m.aux_v.x[i][j][k]) * inv_dthe;
                            the_flag = true;
                        }
                        if ((lnd_ijk && air_jm1 && LAND(i,j-2,k))
                            || (j == 1 && lnd_ijk && air_jm1)) {
                            dv_dthe = (m.aux_v.x[i][j-1][k] - m.aux_v.x[i][j][k]) * inv_dthe;
                            the_flag = true;
                        }
                    }

                    // phi direction
                    const bool k_inner = (k > 2) && (k < m.km-2);
                    if (k_inner) {
                        const int8_t air_kp1 = !lnd_kp1;
                        const int8_t air_km1 = !lnd_km1;

                        // Fluid cell with land neighbour in phi
                        if (!phi_flag && !lnd_ijk) {
                            if (!air_km1 && air_kp1 && !lnd_kp2) {
                                dw_dphi = (-3.0 * m.aux_w.x[i][j][k] + 4.0 * m.aux_w.x[i][j][k+1]
                                           - m.aux_w.x[i][j][k+2]) * inv_2dphi;
                                phi_flag = true;
                            } else if (!air_kp1 && air_km1 && !lnd_km2) {
                                dw_dphi = (3.0 * m.aux_w.x[i][j][k] - 4.0 * m.aux_w.x[i][j][k-1]
                                           + m.aux_w.x[i][j][k-2]) * inv_2dphi;
                                phi_flag = true;
                            }
                        }

                        if (lnd_ijk && air_kp1 && !lnd_kp2) {
                            dw_dphi = (-3.0*m.aux_w.x[i][j][k] + 4.0*m.aux_w.x[i][j][k+1]
                                       - m.aux_w.x[i][j][k+2]) * inv_2dphi;
                            phi_flag = true;
                        }
                        if (lnd_ijk && air_km1 && !lnd_km2) {
                            dw_dphi = -(-3.0*m.aux_w.x[i][j][k] + 4.0*m.aux_w.x[i][j][k-1]
                                        - m.aux_w.x[i][j][k-2]) * inv_2dphi;
                            phi_flag = true;
                        }
                        if ((lnd_ijk && air_kp1 && lnd_kp2)
                            || (k == m.km-2 && !lnd_ijk && lnd_kp1)) {
                            dw_dphi = (m.aux_w.x[i][j][k+1] - m.aux_w.x[i][j][k]) * inv_dphi;
                            phi_flag = true;
                        }
                        if ((lnd_ijk && air_km1 && lnd_km2)
                            || (k == 1 && lnd_ijk && air_km1)) {
                            dw_dphi = (m.aux_w.x[i][j][k-1] - m.aux_w.x[i][j][k]) * inv_dphi;
                            phi_flag = true;
                        }
                    }

                    // central-difference fallbacks (only for fluid cells fully surrounded by fluid)
                    if (!r_flag)
                        du_dr   = (m.aux_u.x[i+1][j][k] - m.aux_u.x[i-1][j][k]) * inv_2dr;
                    if (!the_flag)
                        dv_dthe = (m.aux_v.x[i][j+1][k] - m.aux_v.x[i][j-1][k]) * inv_2dthe;
                    if (!phi_flag)
                        dw_dphi = (m.aux_w.x[i][j][k+1] - m.aux_w.x[i][j][k-1]) * inv_2dphi;

                    // Solid-wall treatment. A land cell must NOT act as a Poisson pressure
                    // source: the one-sided divergence stencils above sample the adjacent
                    // ocean velocity (e.g. the Drake jet), which injected a spurious high
                    // +p_dyn onto coastal land cells, bled into the ocean pressure and
                    // blocked the through-flow. The first attempt (zero the source ⇒ harmonic
                    // mean of ALL neighbours) over-corrected: it pulled the gap pressure toward
                    // the land average and flattened the venturi gradient, so the Drake flow
                    // stagnated (−u/−w piled up upstream). Instead impose a TRUE Neumann
                    // ∂p/∂n=0 wall — set the land cell to the mean of its FLUID neighbours
                    // only. The fluid then develops its pressure (and the through-gap gradient)
                    // freely; the wall mirrors it, with no spurious source and no cross-coast
                    // force, so the passage flow is preserved. Buried cell (no fluid neighbour) → 0.
                    if (lnd_ijk) {
                        double psum = 0.0; int pn = 0;
                        if (!LAND(i+1,j,k)) { psum += m.p_dyn.x[i+1][j][k]; pn++; }
                        if (!LAND(i-1,j,k)) { psum += m.p_dyn.x[i-1][j][k]; pn++; }
                        if (!LAND(i,j+1,k)) { psum += m.p_dyn.x[i][j+1][k]; pn++; }
                        if (!LAND(i,j-1,k)) { psum += m.p_dyn.x[i][j-1][k]; pn++; }
                        if (!lnd_kp1)       { psum += m.p_dyn.x[i][j][k+1]; pn++; }
                        if (!lnd_km1)       { psum += m.p_dyn.x[i][j][k-1]; pn++; }
                        m.p_dyn.x[i][j][k] = (pn > 0) ? psum / pn : 0.0;
                    } else {
                        // pressure update (fluid interior) — clamp the divergence source
                        // so a steep-orography velocity spike cannot drive p_dyn unbounded.
                        double div_src = du_dr   * geo.exp_rm
                                       + dv_dthe * geo.inv_rm
                                       + dw_dphi * geo.inv_rmsinthe;

                        // Anelastic divergence source (step 3). Because rho_bar depends on
                        // r only, (1/rho_bar) div(rho_bar u*) - div(u*) is this one extra
                        // term on the radial component. Added to du_dr's own metric factor,
                        // not to the whole source, and formed on aux_u for the same reason
                        // the rest of the source is: it is the divergence of the PROVISIONAL
                        // velocity the projection has to remove.
                        if (anelastic)
                            div_src += m.aux_u.x[i][j][k] * dlnrho[i] * geo.exp_rm;

                        // ATOM_METRIC_DIVERGENCE (lib/Utils.h) — the spherical divergence also
                        // carries +2u/r and +v*cot(theta)/r. Formed from aux_*, because the source
                        // is the divergence of the PROVISIONAL velocity that the projection has to
                        // remove. Read the warning in lib/Utils.h before switching this on.
                        if (metric_div) {
                            div_src += (2.0 * m.aux_u.x[i][j][k]
                                        + m.aux_v.x[i][j][k] * geo.costhe / geo.sinthe)
                                       * geo.inv_rm;
                        }
                        const double src_max = denom * p_dyn_cap;

                        n_src_cells++;
                        if (!is_finite_safe(div_src))   div_src = 0.0;
                        else if (div_src >  src_max)  { div_src =  src_max; n_src_clamped++; }
                        else if (div_src < -src_max)  { div_src = -src_max; n_src_clamped++; }

                        m.p_dyn.x[i][j][k] =
                            ((m.p_dyn.x[i+1][j][k] + m.p_dyn.x[i-1][j][k]) * num1
                           + (m.p_dyn.x[i][j+1][k] + m.p_dyn.x[i][j-1][k]) * num2
                           + (m.p_dyn.x[i][j][k+1] + m.p_dyn.x[i][j][k-1]) * num3
                           + (m.p_dyn.x[i+1][j][k] - m.p_dyn.x[i-1][j][k]) * num_a
                           - div_src) * inv_denom;
                    }
                } // k
            } // j
        } // i

        } // colour

        #undef LAND

        // Radial boundary extrapolation.
        //
        // At i=0 the old cubic p[0] = p[3] - 3*p[2] + 3*p[1] has condition number ~7
        // and amplifies any cliff-cell pressure into the reference layer by up to 3×.
        // Over many pressure-solve calls this compounds into huge ± spikes at steep
        // topography (e.g. Himalaya). Replaced by:
        //   - i_topography[j][k] >= 1 (i=0 inside the mountain): hold p_dyn = 0, matching
        //     bcSolidGround's treatment of fully buried interior cells; acts as a Dirichlet
        //     pin for the otherwise all-Neumann pressure Poisson.
        //   - i_topography[j][k] == 0 (i=0 is a real ocean surface): von Neumann zero-gradient
        //     ∂p/∂n = 0 — physically correct at a free-slip wall.
        // Tried fully zero-gradient (no Dirichlet anywhere) — Poisson becomes singular and
        // diverges to NaN within 100 iters at the pole.
        // At i=im-1 the column is always air, so the cubic is fine and is kept.
        #pragma omp parallel for collapse(2)
        for (int k = 0; k < m.km; k++) {
            for (int j = 0; j < m.jm; j++) {
                if (m.i_topography[j][k] >= 1) {
                    m.p_dyn.x[0][j][k] = 0.0;
                } else {
                    m.p_dyn.x[0][j][k] = m.c43 * m.p_dyn.x[1][j][k]
                                       - m.c13 * m.p_dyn.x[2][j][k];
                }
                // Top (rigid lid): zero-gradient plain copy. The old cubic
                // p[im-1] = p[im-4] − 3·p[im-3] + 3·p[im-2] is the SAME high-condition (~7,
                // up to 3× amplification) stencil removed at i=0 above, and it has the same
                // failure mode here: over high orography (user-observed Himalaya k=87 / j=62)
                // the strong near-top vertical p_dyn gradient is amplified 3× into a spurious
                // pressure MAXIMUM pegged at the |p_dyn| ceiling (+3 non-dim ≈ 3000 hPa) right
                // at the lid, which then imprints on the orography-oriented velocities. The old
                // "column is air so the cubic is fine" rationale was wrong — the amplification
                // is independent of air/land. At a rigid lid ∂p_dyn/∂r ≈ 0 (no through-lid
                // flow), so a non-amplifying zero-gradient copy is both physical and stable
                // (amplification 1, matching the θ-pole plain-copy BC below). 2026-06-10.
                m.p_dyn.x[m.im-1][j][k] = m.p_dyn.x[m.im-2][j][k];
            }
        }

        // Theta-pole BC: plain copy (zero-gradient, no extrapolation).
        // The 2nd-order Neumann form p[0] = (4/3)·p[1] − (1/3)·p[2] amplifies polar
        // grid noise by 4/3 per call. The Poisson iteration is rerun every moist_stride,
        // and the metric 1/sin²θ already blows up near the pole — compounding the 4/3
        // factor drove p_dyn to NaN at 90°N within ~150 iters, then dpdr/dpdthe/dpdphi
        // propagated NaN into u,v,w in the next RHS call. Plain copy = amplification 1
        // and matches the axisymmetric-pole assumption used for the same fields in bcTheta.
        #pragma omp parallel for collapse(2)
        for (int k = 0; k < m.km; k++) {
            for (int i = 0; i < m.im; i++) {
                m.p_dyn.x[i][0][k]      = m.p_dyn.x[i][1][k];
                m.p_dyn.x[i][m.jm-1][k] = m.p_dyn.x[i][m.jm-2][k];
            }
        }

        // Phi boundary average
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                m.p_dyn.x[i][j][0]      = m.c43 * m.p_dyn.x[i][j][1]      - m.c13 * m.p_dyn.x[i][j][2];
                m.p_dyn.x[i][j][m.km-1] = m.c43 * m.p_dyn.x[i][j][m.km-2] - m.c13 * m.p_dyn.x[i][j][m.km-3];
                m.p_dyn.x[i][j][0] = m.p_dyn.x[i][j][m.km-1]
                    = (m.p_dyn.x[i][j][0] + m.p_dyn.x[i][j][m.km-1]) / 2.0;
            }
        }

        // Hard ceiling on |p_dyn| (applied after every sweep, incl. boundaries). The
        // gradient of p_dyn is the dominant RHS force (presgrad ≫ buoyancy/Coriolis), so
        // bounding p_dyn directly bounds the velocity forcing and breaks the dry steep-
        // orography loop p_dyn → ∇p → w pegs ±100 m/s → divergence → p_dyn. The source
        // clamp above limits injection, but p_dyn still accumulates across the cliff via
        // the Laplacian terms (reached ~7.7 non-dim ≈ 7800 hPa with the source clamp
        // alone); this caps the accumulated result. p_dyn is non-dim'd by p_0=1013.25 hPa
        // p_dyn_ceiling is a NaN/extreme backstop only — NOT the primary control. A tight
        // value (1.0) clipped p_dyn below the level the projection legitimately needs during
        // the violent dry spin-up, leaving residual divergence that blew up at the NZ Alps
        // (iter 288). 10.0 (≈10000 hPa) sits above the accumulated steep-orography value
        // (~7.7 with p_dyn_cap=2) so it never clips normal operation, only the runaway and
        // NaN/Inf (→0).
        //
        // PHASE-DEPENDENT TEST (2026-06-09): the iter-483 runaway is in the VISCOUS/moist
        // phase, where p_dyn pegs −10 over the Tian Shan/Pamir and the −10→0 vertical gradient
        // becomes a −16.8 pgr that drives the velocity (see [[project_upper_velocity_secular_growth]]).
        // Bound the FORCE there by lowering the ceiling to 3.0 once past moist onset (iter 300),
        // while KEEPING 10.0 through the dry spin-up so we don't re-trigger the NZ-Alps failure
        // the tight uniform value caused. Risk: the divergence may legitimately need a larger
        // p_dyn (converging the Poisson gave an even bigger gradient), so a too-tight cap can
        // leave residual velocity divergence that blows up elsewhere — this run tests that.
        // NOTE (2026-06-11): a SMOOTH tanh saturation p=C·tanh(p/C) was tried here in
        // place of this hard clamp (hypothesis: the flat ±ceiling block's discontinuous
        // edge feeds a spurious vertical pgr). FALSIFIED — from atm_restart_500.bin it
        // still NaN'd at the SAME iter 532, merely RELOCATING the seed from i=10/44°N/128°E
        // (near-surface) to i=26/42°N/68°E (Pamir, ~4 km, the documented secular-growth
        // band). Symptom-only / whack-a-mole; the crash time is p_dyn-clip-independent.
        // Reverted to the hard clamp. See [[project_upper_velocity_secular_growth]].
        const double p_dyn_ceiling = (m.total_iter_count > 300) ? 3.0 : 10.0;
        #pragma omp parallel for collapse(3)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    double p = m.p_dyn.x[i][j][k];
                    if (!is_finite_safe(p))      m.p_dyn.x[i][j][k] = 0.0;
                    else if (p >  p_dyn_ceiling) m.p_dyn.x[i][j][k] =  p_dyn_ceiling;
                    else if (p < -p_dyn_ceiling) m.p_dyn.x[i][j][k] = -p_dyn_ceiling;
                }
            }
        }

        // ------------------------------------------------------------------
        // ATHAD: measure the divergence the projection is supposed to be removing.
        //
        // Nothing measured it, so "the projection enforces div u = 0" was an assumption
        // about the code rather than an observation of it — and the file's own note says
        // div and grad are NOT discretely adjoint here (the Laplacian uses squared metric
        // factors where the source and the gradient correction are single-power), which
        // means the projection cannot remove exactly what it measured. This is the number
        // that says how far off it is, and the A/B control for ATM_POISSON_METRIC_FIX.
        //
        // Formed on the ACTUAL velocity u/v/w with the same expression as div_src above,
        // so the two are directly comparable.
        // Both forms are reported, always, and that is the point: div(u) is the Boussinesq
        // residual and div(rho_bar u)/rho_bar the anelastic one, and a projection that
        // drives one to zero drives the other to the size of the term between them. Which
        // number is "the" residual depends on which continuity the model claims, so
        // printing only the enforced one would make either scheme look good.
        if (verbose) {
            double d2 = 0.0, dmax = 0.0, a2 = 0.0, amax = 0.0;
            long   n  = 0;
            const bool have_rho = ((int)m.m_dlnrho_dr.size() == m.im);
            #pragma omp parallel for collapse(2) schedule(static) reduction(+:d2,a2,n) \
                    reduction(max:dmax,amax)
            for (int i = 1; i < m.im-1; i++) {
                for (int j = 1; j < m.jm-1; j++) {
                    const double rm      = m.rad.z[i];
                    const double exp_rm  = 1.0 / (rm + 1.0);
                    const double rmet    = m.metricRadius(rm);
                    const double sinthe  = sinthe_table[j];
                    const double inv_rm  = 1.0 / rmet;
                    const double inv_rms = 1.0 / (rmet * sinthe);
                    const double cotanthe = cos(m.the.z[j]) / sinthe;
                    const double dlr     = have_rho ? m.m_dlnrho_dr[i] : 0.0;

                    for (int k = 1; k < m.km-1; k++) {
                        if (land[i*m.jm*m.km + j*m.km + k]) continue;
                        double d = (m.u.x[i+1][j][k] - m.u.x[i-1][j][k]) * inv_2dr   * exp_rm
                                 + (m.v.x[i][j+1][k] - m.v.x[i][j-1][k]) * inv_2dthe * inv_rm
                                 + (m.w.x[i][j][k+1] - m.w.x[i][j][k-1]) * inv_2dphi * inv_rms;
                        if (metric_div)
                            d += (2.0 * m.u.x[i][j][k] + m.v.x[i][j][k] * cotanthe) * inv_rm;
                        if (!is_finite_safe(d)) continue;
                        const double a = d + m.u.x[i][j][k] * dlr * exp_rm;
                        d2 += d * d;
                        a2 += a * a;
                        if (std::fabs(d) > dmax) dmax = std::fabs(d);
                        if (std::fabs(a) > amax) amax = std::fabs(a);
                        n++;
                    }
                }
            }
            cout << "      ATOM: div(u) after projection  rms = "
                 << std::scientific << std::setprecision(3) << ((n > 0) ? sqrt(d2 / n) : 0.0)
                 << "   max = " << dmax << std::fixed
                 << "   (metric_fix = " << (poisson_metric_fix ? "on" : "off")
                 << ", anelastic = " << (anelastic ? "on" : "off") << ")" << endl;
            if (have_rho)
                cout << "      ATOM: div(rho u)/rho           rms = "
                     << std::scientific << std::setprecision(3) << ((n > 0) ? sqrt(a2 / n) : 0.0)
                     << "   max = " << amax << std::fixed << endl;
            cout << "      ATOM: divergence source clamped at p_dyn_cap in "
                 << n_src_clamped << " of " << n_src_cells << " fluid cells ("
                 << std::fixed << std::setprecision(2)
                 << ((n_src_cells > 0) ? 100.0 * (double)n_src_clamped / (double)n_src_cells : 0.0)
                 << " %)" << endl;
        }

        auto end = std::chrono::high_resolution_clock::now();
        if (verbose) {
            auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            printf(" time measured: %.3f seconds for PressureSolverAtm\n", elapsed.count() * 1e-9);
            cout << "      ATOM: PressureSolverAtm ended" << endl;
        }
    }

    // One-shot Helmholtz projection of the prescribed initial velocity onto its
    // divergence-free subspace.  Without this, the incremental projection that runs
    // every moist_stride iters spends the first ~100 iters silently destroying the
    // dilatational component of the Hadley/Ferrel profile, taking the prescribed
    // global circulation with it (max u decays from 30 m/s to <0.5 m/s by iter 150).
    // Performing the projection once at t=0 lets the simulation start from a clean
    // divergence-free state; the surviving solenoidal part is preserved and only the
    // unphysical dilatational artefact of the analytical profile is removed.
    //
    // The standard solver loop assembles ∇²p = ∇·v* where v* is the post-RHS
    // intermediate velocity; here we feed v itself as the source by copying v into
    // aux_u/aux_v/aux_w. After n_sweeps Jacobi passes, p_dyn approximates the
    // projection pressure; we then apply v ← v − ∇p in the same metric form used
    // by the time-stepping RHS, and reset p_dyn to 0 so the next RK4 call does not
    // double-correct via its own −∂p/∂r term.
    void project_initial_velocity(int n_sweeps = 200)
    {
        using namespace std;
        cout << endl << endl << "      ATOM: project_initial_velocity ("
             << n_sweeps << " Jacobi sweeps)" << endl;
        auto t0 = std::chrono::high_resolution_clock::now();

        // Step 1 — seed Poisson source from current velocity, zero p_dyn.
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    m.aux_u.x[i][j][k] = m.u.x[i][j][k];
                    m.aux_v.x[i][j][k] = m.v.x[i][j][k];
                    m.aux_w.x[i][j][k] = m.w.x[i][j][k];
                    m.p_dyn.x[i][j][k] = 0.0;
                }
            }
        }

        // Step 2 — quiet Jacobi iteration until p_dyn approximates the projection
        // pressure.  Each call to run() also re-applies the i, theta, and phi BCs on
        // p_dyn, so polar/topographic anchors stay consistent with the time loop.
        for (int s = 0; s < n_sweeps; s++) {
            run(false);
        }

        // Step 3 — gradient correction v ← v − ∇p_dyn in the interior.
        // Metric factors match the rhs_u/v/w pressure-gradient term so the magnitudes
        // are consistent with the rest of the code.  Boundary cells (i, theta, phi
        // outer faces) are left untouched; bcRadius / bcTheta / bcPhi will re-impose
        // their patterns at the next call.
        const double inv_2dr   = 1.0 / (2.0 * m.dr);
        const double inv_2dthe = 1.0 / (2.0 * m.dthe);
        const double inv_2dphi = 1.0 / (2.0 * m.dphi);

        std::vector<double> sinthe_tab(m.jm);
        for (int j = 0; j < m.jm; j++) {
            sinthe_tab[j] = sin(m.the.z[j]);
            if (sinthe_tab[j] < 0.55) sinthe_tab[j] = 0.55;   // metric floor ~57° — keep in sync with the RHS geometry floor (RungeKutta_Atm*.cpp)
        }

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {
                const double rm           = m.rad.z[i];
                const double exp_rm       = 1.0 / (rm + 1.0);   // grid coordinate, not the radius
                // The gradient correction must use the SAME metric the source and the RHS use,
                // or the projection stops being a projection. ATM_METRIC_RADIUS; identity when off.
                const double rmet         = m.metricRadius(rm);
                const double inv_rm       = 1.0 / rmet;
                const double inv_rmsinthe = 1.0 / (rmet * sinthe_tab[j]);

                for (int k = 1; k < m.km-1; k++) {
                    const double dpdr   = (m.p_dyn.x[i+1][j][k] - m.p_dyn.x[i-1][j][k]) * inv_2dr;
                    const double dpdthe = (m.p_dyn.x[i][j+1][k] - m.p_dyn.x[i][j-1][k]) * inv_2dthe;
                    const double dpdphi = (m.p_dyn.x[i][j][k+1] - m.p_dyn.x[i][j][k-1]) * inv_2dphi;

                    m.u.x[i][j][k] -= dpdr   * exp_rm;
                    m.v.x[i][j][k] -= dpdthe * inv_rm;
                    m.w.x[i][j][k] -= dpdphi * inv_rmsinthe;
                }
            }
        }

        // Step 4 — clear p_dyn and aux so the time loop starts fresh.
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    m.p_dyn.x[i][j][k] = 0.0;
                    m.aux_u.x[i][j][k] = 0.0;
                    m.aux_v.x[i][j][k] = 0.0;
                    m.aux_w.x[i][j][k] = 0.0;
                }
            }
        }

        auto t1 = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0);
        printf("      ATOM: project_initial_velocity ended (%.3fs)\n", elapsed.count() * 1e-9);
    }

private:
    cAtmosphereModel& m;
};

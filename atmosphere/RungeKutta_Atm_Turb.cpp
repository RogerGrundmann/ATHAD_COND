/*
 * Atmosphere General Circulation Modell (AGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * class to produce results by the Runge-Kutta solution scheme
*/

#include "cAtmosphereModel.h"

#include <cstdint>
#include <cstring>

using namespace std;


// Boussinesq buoyancy base state: the area-weighted (sin θ, θ = colatitude)
// horizontal mean of the non-dim temperature at each radial level i, taken over
// air cells only. Buoyancy is then driven by (t − t_ref_level[i]) so the body
// force has zero mean at every height and only horizontal temperature contrasts
// accelerate the flow. Refilled once per RK4 step. (Moved here from the former
// RungeKutta_Atm.cpp when the separate laminar solver was dropped.)
void cAtmosphereModel::computeLevelMeanTemperature(){
    if((int)t_ref_level.size() != im) t_ref_level.assign(im, 1.0);

    #pragma omp parallel for schedule(static)
    for(int i = 0; i < im; i++){
        double sum = 0.0, wsum = 0.0;
        for(int j = 0; j < jm; j++){
            const double w = sin(the.z[j]);                 // spherical area weight
            for(int k = 0; k < km; k++){
                if(AtomUtils::is_air(h, i, j, k)){
                    sum  += w * t.x[i][j][k];
                    wsum += w;
                }
            }
        }
        t_ref_level[i] = (wsum > 0.0) ? sum / wsum : 1.0;   // all-land level → neutral
    }
}


void cAtmosphereModel::solveRungeKutta_Atmosphere_Turb(){
    cout << endl << "      solveRungeKutta_Atmosphere_Turb begin" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    computeLevelMeanTemperature();   // refresh buoyancy base state t_ref_level[i]

    const double half_dt  = 0.5 * dt;
    const double dt_sixth = dt / 6.0;

    // Velocity caps removed (see RungeKutta_Atm.cpp note).  Multi-sweep
    // PressureSolverAtm is the structural fix; the caps were breaking mass
    // conservation by clipping the velocity at the same time the projection was
    // trying to enforce ∇·v = 0.  TKE cap retained — it's mass-neutral.

    // TKE cap.  After the buoyancy/projection/pole fixes unlocked the global circulation,
    // tke was reaching ~1e98 m²/s² at a single high-altitude polar cell by iter 100.
    // Removing the conservative-form divergence correction on transport_tke knocked the
    // magnitude back but the runaway persisted via the −0.667·tke·∇·v term inside prod
    // (positive contribution in convergent regions, linear in tke).  Pending a structural
    // fix to the source assembly, hard-clamp |k| to a physically generous bound at every
    // RK4 stage so the simulation cannot diverge.  1000 m²/s² = 50% above the strongest
    // hurricane-core TKE on record; any legitimate physics fits well inside this cap.
    const double tke_max_phys = 1000.0;                                  // [m²/s²]
    const double tke_max_nd   = tke_max_phys / (u_0 * u_0);

    // NaN-safe clamp.  std::clamp uses unguarded `<` comparisons; under -ffast-math
    // (-ffinite-math-only) the compiler assumes operands are finite and the result
    // on NaN input is undefined.  Observed: tke reaching 7×10²⁰⁰ at one cell despite
    // a 1000 m²/s² cap, because a NaN intermediate leaked through std::clamp.  Detect
    // non-finite values by IEEE-754 exponent bits and replace with lo (no ffast-math
    // hazard).  Same trick we used in Paraview_Atm.cpp::safe_val.
    auto safe_clamp = [](double v, double lo, double hi) -> double {
        std::uint64_t bits;
        std::memcpy(&bits, &v, sizeof(bits));
        if ((bits & 0x7FF0000000000000ULL) == 0x7FF0000000000000ULL) return lo;
        return (v < lo) ? lo : (v > hi) ? hi : v;
    };

    // Grid-spacing reciprocals — constant for the entire grid
    const double inv_2dr   = 1.0 / (2.0 * dr);
    const double inv_2dthe = 1.0 / (2.0 * dthe);
    const double inv_2dphi = 1.0 / (2.0 * dphi);
    const double inv_dr2   = 1.0 / (dr * dr);
    const double inv_dthe2 = 1.0 / (dthe * dthe);
    const double inv_dphi2 = 1.0 / (dphi * dphi);

    // Precompute sin/cos tables — only depend on j
    std::vector<double> sinthe_tbl(jm), costhe_tbl(jm);
    for(int j = 0; j < jm; j++){
        sinthe_tbl[j] = sin(the.z[j]);
        if(sinthe_tbl[j] < 0.55) sinthe_tbl[j] = 0.55;   // metric floor ~57° (was 0.4/~66°): caps the 1/sinθ amplification of the high-lat coastal pressure-gradient force that blows up at i=1
        costhe_tbl[j] = cos(the.z[j]);
    }

    #pragma omp parallel for collapse(2) schedule(static)
    for(int i = 1; i < im-1; i++){
        for(int j = 1; j < jm-1; j++){

            // ---- Build geometry struct ONCE per (i,j) ----
            CellGeometry geo;

            geo.rm       = rad.z[i];
            geo.rm2      = geo.rm * geo.rm;
            geo.exp_rm   = metricExpRm(geo.rm);
            geo.exp_2_rm = geo.exp_rm * geo.exp_rm;
            geo.curv     = metricCurv(geo.rm);

            geo.sinthe  = sinthe_tbl[j];
            geo.sinthe2 = geo.sinthe * geo.sinthe;
            geo.costhe  = costhe_tbl[j];

            // The 1/r factors carry the PLANETARY radius; exp_rm above carries the stretched grid
            // coordinate. Those are two different things and ATM_METRIC_RADIUS separates them.
            // Identity when the knob is off (metricRadius returns rm), so this is bit-identical
            // by default. Keep in sync with PressureSolverAtm and TurbulenceAtm: a solver that
            // disagrees with the RHS about the metric is the defect class this repo keeps finding.
            const double rmet  = metricRadius(geo.rm);
            const double rmet2 = rmet * rmet;

            geo.inv_rm               = 1.0 / rmet;
            geo.inv_rm2              = 1.0 / rmet2;
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

            for(int k = 1; k < km-1; k++){

                // Cache time-level-n values (read once, used 4×)
                double tn_ijk   = tn.x[i][j][k];
                double un_ijk   = un.x[i][j][k];
                double vn_ijk   = vn.x[i][j][k];
                double wn_ijk   = wn.x[i][j][k];
                double cn_ijk   = cn.x[i][j][k];
                double cldn_ijk = cloudn.x[i][j][k];
                double icen_ijk = icen.x[i][j][k];
                double grn_ijk  = grn.x[i][j][k];
                double co2n_ijk = co2n.x[i][j][k];
                double tken_ijk = tken.x[i][j][k];
                double disn_ijk = disn.x[i][j][k];

                // --- Stage 1 ---
                cAtmosphereModel::RHS_Atmosphere_Turb(i, j, k, geo);

                double kt1   = rhs_t.x[i][j][k];
                double ku1   = rhs_u.x[i][j][k];
                double kv1   = rhs_v.x[i][j][k];
                double kw1   = rhs_w.x[i][j][k];
                double kc1   = rhs_c.x[i][j][k];
                double kcl1  = rhs_cloud.x[i][j][k];
                double ki1   = rhs_ice.x[i][j][k];
                double kg1   = rhs_g.x[i][j][k];
                double kco1  = rhs_co2.x[i][j][k];
                double ktke1 = rhs_tke.x[i][j][k];
                double kdis1 = rhs_dis.x[i][j][k];

                t.x[i][j][k]     = tn_ijk   + kt1   * half_dt;
                u.x[i][j][k]     = un_ijk   + ku1   * half_dt;
                v.x[i][j][k]     = vn_ijk   + kv1   * half_dt;
                w.x[i][j][k]     = wn_ijk   + kw1   * half_dt;
                // Positivity clamp on the moisture scalars at every RK4 stage (as tke/dis
                // below): mass mixing ratios are physically >= 0. The spherical-operator
                // update can drive a near-zero cloud/c/ice cell slightly negative at sharp
                // coastal gradients; clamping here guarantees the prognostic field is never
                // negative (was previously only fixed downstream by storeIntermediateData3D).
                c.x[i][j][k]     = std::max(0.0, cn_ijk   + kc1   * half_dt);
                cloud.x[i][j][k] = std::max(0.0, cldn_ijk + kcl1  * half_dt);
                ice.x[i][j][k]   = std::max(0.0, icen_ijk + ki1   * half_dt);
                gr.x[i][j][k]    = std::max(0.0, grn_ijk  + kg1   * half_dt);
                co2.x[i][j][k]   = co2n_ijk + kco1  * half_dt;
                tke.x[i][j][k]   = safe_clamp(tken_ijk + ktke1 * half_dt, 0.0,     tke_max_nd);
                dis.x[i][j][k]   = std::max(1.0e-10, disn_ijk + kdis1 * half_dt);

                // --- Stage 2 ---
                cAtmosphereModel::RHS_Atmosphere_Turb(i, j, k, geo);

                double kt2   = rhs_t.x[i][j][k];
                double ku2   = rhs_u.x[i][j][k];
                double kv2   = rhs_v.x[i][j][k];
                double kw2   = rhs_w.x[i][j][k];
                double kc2   = rhs_c.x[i][j][k];
                double kcl2  = rhs_cloud.x[i][j][k];
                double ki2   = rhs_ice.x[i][j][k];
                double kg2   = rhs_g.x[i][j][k];
                double kco2  = rhs_co2.x[i][j][k];
                double ktke2 = rhs_tke.x[i][j][k];
                double kdis2 = rhs_dis.x[i][j][k];

                t.x[i][j][k]     = tn_ijk   + kt2   * half_dt;
                u.x[i][j][k]     = un_ijk   + ku2   * half_dt;
                v.x[i][j][k]     = vn_ijk   + kv2   * half_dt;
                w.x[i][j][k]     = wn_ijk   + kw2   * half_dt;
                c.x[i][j][k]     = std::max(0.0, cn_ijk   + kc2   * half_dt);
                cloud.x[i][j][k] = std::max(0.0, cldn_ijk + kcl2  * half_dt);
                ice.x[i][j][k]   = std::max(0.0, icen_ijk + ki2   * half_dt);
                gr.x[i][j][k]    = std::max(0.0, grn_ijk  + kg2   * half_dt);
                co2.x[i][j][k]   = co2n_ijk + kco2  * half_dt;
                tke.x[i][j][k]   = safe_clamp(tken_ijk + ktke2 * half_dt, 0.0,     tke_max_nd);
                dis.x[i][j][k]   = std::max(1.0e-10,  disn_ijk + kdis2 * half_dt);

                // --- Stage 3 ---
                cAtmosphereModel::RHS_Atmosphere_Turb(i, j, k, geo);

                double kt3   = rhs_t.x[i][j][k];
                double ku3   = rhs_u.x[i][j][k];
                double kv3   = rhs_v.x[i][j][k];
                double kw3   = rhs_w.x[i][j][k];
                double kc3   = rhs_c.x[i][j][k];
                double kcl3  = rhs_cloud.x[i][j][k];
                double ki3   = rhs_ice.x[i][j][k];
                double kg3   = rhs_g.x[i][j][k];
                double kco3  = rhs_co2.x[i][j][k];
                double ktke3 = rhs_tke.x[i][j][k];
                double kdis3 = rhs_dis.x[i][j][k];

                t.x[i][j][k]     = tn_ijk   + kt3   * dt;
                u.x[i][j][k]     = un_ijk   + ku3   * dt;
                v.x[i][j][k]     = vn_ijk   + kv3   * dt;
                w.x[i][j][k]     = wn_ijk   + kw3   * dt;
                c.x[i][j][k]     = std::max(0.0, cn_ijk   + kc3   * dt);
                cloud.x[i][j][k] = std::max(0.0, cldn_ijk + kcl3  * dt);
                ice.x[i][j][k]   = std::max(0.0, icen_ijk + ki3   * dt);
                gr.x[i][j][k]    = std::max(0.0, grn_ijk  + kg3   * dt);
                co2.x[i][j][k]   = co2n_ijk + kco3  * dt;
                tke.x[i][j][k]   = safe_clamp(tken_ijk + ktke3 * dt, 0.0,     tke_max_nd);
                dis.x[i][j][k]   = std::max(1.0e-10,  disn_ijk + kdis3 * dt);

                // --- Stage 4 + Final combination ---
                cAtmosphereModel::RHS_Atmosphere_Turb(i, j, k, geo);

                double kt4   = rhs_t.x[i][j][k];
                double ku4   = rhs_u.x[i][j][k];
                double kv4   = rhs_v.x[i][j][k];
                double kw4   = rhs_w.x[i][j][k];
                double kc4   = rhs_c.x[i][j][k];
                double kcl4  = rhs_cloud.x[i][j][k];
                double ki4   = rhs_ice.x[i][j][k];
                double kg4   = rhs_g.x[i][j][k];
                double kco4  = rhs_co2.x[i][j][k];
                double ktke4 = rhs_tke.x[i][j][k];
                double kdis4 = rhs_dis.x[i][j][k];

                t.x[i][j][k]     = tn_ijk   + (kt1   + 2.0*kt2   + 2.0*kt3   + kt4)   * dt_sixth;
                u.x[i][j][k]     = un_ijk   + (ku1   + 2.0*ku2   + 2.0*ku3   + ku4)   * dt_sixth;
                v.x[i][j][k]     = vn_ijk   + (kv1   + 2.0*kv2   + 2.0*kv3   + kv4)   * dt_sixth;
                w.x[i][j][k]     = wn_ijk   + (kw1   + 2.0*kw2   + 2.0*kw3   + kw4)   * dt_sixth;
                c.x[i][j][k]     = std::max(0.0, cn_ijk   + (kc1   + 2.0*kc2   + 2.0*kc3   + kc4)   * dt_sixth);
                cloud.x[i][j][k] = std::max(0.0, cldn_ijk + (kcl1  + 2.0*kcl2  + 2.0*kcl3  + kcl4)  * dt_sixth);
                ice.x[i][j][k]   = std::max(0.0, icen_ijk + (ki1   + 2.0*ki2   + 2.0*ki3   + ki4)   * dt_sixth);
                gr.x[i][j][k]    = std::max(0.0, grn_ijk  + (kg1   + 2.0*kg2   + 2.0*kg3   + kg4)   * dt_sixth);
                co2.x[i][j][k]   = co2n_ijk + (kco1  + 2.0*kco2  + 2.0*kco3  + kco4)  * dt_sixth;
                tke.x[i][j][k]   = safe_clamp(tken_ijk + (ktke1 + 2.0*ktke2 + 2.0*ktke3 + ktke4) * dt_sixth, 0.0, tke_max_nd);
                dis.x[i][j][k]   = std::max(1.0e-10, disn_ijk + (kdis1 + 2.0*kdis2 + 2.0*kdis3 + kdis4) * dt_sixth);
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for solveRungeKutta_Atmosphere_Turb\n", elapsed.count() * 1e-9);

    cout << "      solveRungeKutta_Atmosphere_Turb end" << endl;
}

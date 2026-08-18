/*
 * Atmosphere General Circulation Modell (AGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 *
 * class to combine the right hand sides of the differential equations for the Runge-Kutta scheme
*/

#include "SaturationH2O.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>   // getenv/atof for the remaining A/B knobs
#include "cAtmosphereModel.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;

// minmod slope limiter on two undivided differences. Returns 0 when a,b have
// opposite signs (a local extremum) — the property that makes the limited
// advective gradient positivity-preserving: at a local min/max the advection
// term vanishes, so it cannot manufacture a new under/overshoot.
static inline double minmod(double a, double b){
    return (a * b <= 0.0) ? 0.0 : (std::fabs(a) < std::fabs(b) ? a : b);
}


void cAtmosphereModel::RHS_Atmosphere_Turb(int i, int j, int k, const CellGeometry& geo){

    // All geometric quantities come from the precomputed struct —
    // NO sin(), cos(), division, or reciprocal computation here.
    const double exp_rm             = geo.exp_rm;
    const double exp_2_rm           = geo.exp_2_rm;
    const double sinthe             = geo.sinthe;
    const double sinthe2            = geo.sinthe2;
    const double costhe             = geo.costhe;
    const double inv_rm             = geo.inv_rm;
    const double inv_rm2            = geo.inv_rm2;
    const double inv_rmsinthe       = geo.inv_rmsinthe;
    const double inv_rm2sinthe      = geo.inv_rm2sinthe;
    const double inv_rm2sinthe2     = geo.inv_rm2sinthe2;
    const double costhe_inv_rm2sinthe = geo.costhe_inv_rm2sinthe;

    const double inv_2dr   = geo.inv_2dr;
    const double inv_2dthe = geo.inv_2dthe;
    const double inv_2dphi = geo.inv_2dphi;
    const double inv_dr2   = geo.inv_dr2;
    const double inv_dthe2 = geo.inv_dthe2;
    const double inv_dphi2 = geo.inv_dphi2;

    // Cache local cell velocities
    double u_ijk = u.x[i][j][k];
    double v_ijk = v.x[i][j][k];
    double w_ijk = w.x[i][j][k];

    bool land_ijk = is_land(h, i, j, k);

    bool r_flag = false, the_flag = false, phi_flag = false;

    // ---- First-order derivatives storage ----
    double dudr, dvdr, dwdr, dtdr, dpdr, dcdr, dclouddr, dicedr, dgdr, dcodr, dtkedr, ddisdr, dnuedr;
    double dudthe, dvdthe, dwdthe, dtdthe, dpdthe, dcdthe, dclouddthe, dicedthe, dgdthe, dcodthe, dtkedthe, ddisdthe, dnuedthe;
    double dudphi, dvdphi, dwdphi, dtdphi, dpdphi, dcdphi, dclouddphi, dicedphi, dgdphi, dcodphi, dtkedphi, ddisdphi, dnuedphi;

    // ---- Second-order derivatives storage ----
    double d2udr2, d2vdr2, d2wdr2, d2tdr2, d2cdr2, d2clouddr2, d2icedr2, d2gdr2, d2codr2, d2tkedr2, d2disdr2;
    double d2udthe2, d2vdthe2, d2wdthe2, d2tdthe2, d2cdthe2, d2clouddthe2, d2icedthe2, d2gdthe2, d2codthe2, d2tkedthe2, d2disdthe2;
    double d2udphi2, d2vdphi2, d2wdphi2, d2tdphi2, d2cdphi2, d2clouddphi2, d2icedphi2, d2gdphi2, d2codphi2, d2tkedphi2, d2disdphi2;


    // ===== R-direction derivatives =====
    if(i < im-2 && land_ijk && is_air(h, i+1, j, k)){
        // Forward-biased stencil near land-air boundary
        #define COMPUTE_DR_B(FIELD, d1, d2) \
            d1 = (-3.0*FIELD.x[i][j][k] + 4.0*FIELD.x[i+1][j][k] - FIELD.x[i+2][j][k]) * inv_2dr; \
            d2 = (FIELD.x[i][j][k] - 2.0*FIELD.x[i+1][j][k] + FIELD.x[i+2][j][k]) * inv_dr2;

        COMPUTE_DR_B(u, dudr, d2udr2)
        COMPUTE_DR_B(v, dvdr, d2vdr2)
        COMPUTE_DR_B(w, dwdr, d2wdr2)
        COMPUTE_DR_B(t, dtdr, d2tdr2)
        COMPUTE_DR_B(c, dcdr, d2cdr2)
        COMPUTE_DR_B(cloud, dclouddr, d2clouddr2)
        COMPUTE_DR_B(ice, dicedr, d2icedr2)
        COMPUTE_DR_B(gr, dgdr, d2gdr2)
        COMPUTE_DR_B(co2, dcodr, d2codr2)
        COMPUTE_DR_B(tke, dtkedr, d2tkedr2)
        COMPUTE_DR_B(dis, ddisdr, d2disdr2)
        dpdr   = (-3.0*p_dyn.x[i][j][k] + 4.0*p_dyn.x[i+1][j][k] - p_dyn.x[i+2][j][k]) * inv_2dr;
        dnuedr = (-3.0*nue.x[i][j][k]   + 4.0*nue.x[i+1][j][k]   - nue.x[i+2][j][k])   * inv_2dr;
        r_flag = true;
        #undef COMPUTE_DR_B
    }
    else if(i == 0){
        #define COMPUTE_DR_F(FIELD, d1, d2) \
            d1 = (-3.0*FIELD.x[0][j][k] + 4.0*FIELD.x[1][j][k] - FIELD.x[2][j][k]) * inv_2dr; \
            d2 = (FIELD.x[0][j][k] - 2.0*FIELD.x[1][j][k] + FIELD.x[2][j][k]) * inv_dr2;

        COMPUTE_DR_F(u, dudr, d2udr2)
        COMPUTE_DR_F(v, dvdr, d2vdr2)
        COMPUTE_DR_F(w, dwdr, d2wdr2)
        COMPUTE_DR_F(t, dtdr, d2tdr2)
        COMPUTE_DR_F(c, dcdr, d2cdr2)
        COMPUTE_DR_F(cloud, dclouddr, d2clouddr2)
        COMPUTE_DR_F(ice, dicedr, d2icedr2)
        COMPUTE_DR_F(gr, dgdr, d2gdr2)
        COMPUTE_DR_F(co2, dcodr, d2codr2)
        COMPUTE_DR_F(tke, dtkedr, d2tkedr2)
        COMPUTE_DR_F(dis, ddisdr, d2disdr2)
        dpdr   = (-3.0*p_dyn.x[0][j][k] + 4.0*p_dyn.x[1][j][k] - p_dyn.x[2][j][k]) * inv_2dr;
        dnuedr = (-3.0*nue.x[0][j][k]   + 4.0*nue.x[1][j][k]   - nue.x[2][j][k])   * inv_2dr;
        r_flag = true;
        #undef COMPUTE_DR_F
    }

    if(!r_flag){
        // Central differences (default path — most common)
        #define COMPUTE_DR_C(FIELD, d1, d2) \
            d1 = (FIELD.x[i+1][j][k] - FIELD.x[i-1][j][k]) * inv_2dr; \
            d2 = (FIELD.x[i+1][j][k] - 2.0*FIELD.x[i][j][k] + FIELD.x[i-1][j][k]) * inv_dr2;

        COMPUTE_DR_C(u, dudr, d2udr2)
        COMPUTE_DR_C(v, dvdr, d2vdr2)
        COMPUTE_DR_C(w, dwdr, d2wdr2)
        COMPUTE_DR_C(t, dtdr, d2tdr2)
        COMPUTE_DR_C(c, dcdr, d2cdr2)
        COMPUTE_DR_C(cloud, dclouddr, d2clouddr2)
        COMPUTE_DR_C(ice, dicedr, d2icedr2)
        COMPUTE_DR_C(gr, dgdr, d2gdr2)
        COMPUTE_DR_C(co2, dcodr, d2codr2)
        COMPUTE_DR_C(tke, dtkedr, d2tkedr2)
        COMPUTE_DR_C(dis, ddisdr, d2disdr2)
        dpdr   = (p_dyn.x[i+1][j][k] - p_dyn.x[i-1][j][k]) * inv_2dr;
        dnuedr = (nue.x[i+1][j][k]   - nue.x[i-1][j][k])   * inv_2dr;
        #undef COMPUTE_DR_C

        // At the first fluid cell directly above terrain (i-1 is a land cell), the
        // centered tke/dis stencil reads tke[i-1]=0 (land), producing a large
        // negative Laplacian that pulls tke → 0 and hence nue → 0 at surface heights.
        // Apply a Neumann (zero-gradient) BC at the terrain surface: replace
        // tke[i-1] with tke[i] so the diffusion flux is one-sided into the atmosphere.
        if (!land_ijk && i > 0 && is_land(h, i-1, j, k)) {
            dtkedr   = (tke.x[i+1][j][k] - tke.x[i][j][k]) * inv_2dr;
            d2tkedr2 = (tke.x[i+1][j][k] - tke.x[i][j][k]) * inv_dr2;
            ddisdr   = (dis.x[i+1][j][k] - dis.x[i][j][k]) * inv_2dr;
            d2disdr2 = (dis.x[i+1][j][k] - dis.x[i][j][k]) * inv_dr2;
        }
    }


    // ===== THETA-direction derivatives =====
    if(j > 2 && j < jm-3){
        if(land_ijk && is_air(h, i, j+1, k) && is_air(h, i, j+2, k)){
            #define COMPUTE_DTHE_B(FIELD, d1, d2) \
                d1 = (-3.0*FIELD.x[i][j][k] + 4.0*FIELD.x[i][j+1][k] - FIELD.x[i][j+2][k]) * inv_2dthe; \
                d2 = (FIELD.x[i][j][k] - 2.0*FIELD.x[i][j+1][k] + FIELD.x[i][j+2][k]) * inv_dthe2;

            COMPUTE_DTHE_B(u, dudthe, d2udthe2) COMPUTE_DTHE_B(v, dvdthe, d2vdthe2)
            COMPUTE_DTHE_B(w, dwdthe, d2wdthe2) COMPUTE_DTHE_B(t, dtdthe, d2tdthe2)
            COMPUTE_DTHE_B(c, dcdthe, d2cdthe2) COMPUTE_DTHE_B(cloud, dclouddthe, d2clouddthe2)
            COMPUTE_DTHE_B(ice, dicedthe, d2icedthe2) COMPUTE_DTHE_B(gr, dgdthe, d2gdthe2)
            COMPUTE_DTHE_B(co2, dcodthe, d2codthe2)
            COMPUTE_DTHE_B(tke, dtkedthe, d2tkedthe2) COMPUTE_DTHE_B(dis, ddisdthe, d2disdthe2)
            dpdthe   = (-3.0*p_dyn.x[i][j][k] + 4.0*p_dyn.x[i][j+1][k] - p_dyn.x[i][j+2][k]) * inv_2dthe;
            dnuedthe = (-3.0*nue.x[i][j][k]   + 4.0*nue.x[i][j+1][k]   - nue.x[i][j+2][k])   * inv_2dthe;
            the_flag = true;
            #undef COMPUTE_DTHE_B
        }
        if(!the_flag && land_ijk && is_air(h, i, j-1, k) && is_air(h, i, j-2, k)){
            #define COMPUTE_DTHE_C(FIELD, d1, d2) \
                d1 = -(-3.0*FIELD.x[i][j][k] + 4.0*FIELD.x[i][j-1][k] - FIELD.x[i][j-2][k]) * inv_2dthe; \
                d2 = -(FIELD.x[i][j][k] - 2.0*FIELD.x[i][j-1][k] + FIELD.x[i][j-2][k]) * inv_dthe2;

            COMPUTE_DTHE_C(u, dudthe, d2udthe2) COMPUTE_DTHE_C(v, dvdthe, d2vdthe2)
            COMPUTE_DTHE_C(w, dwdthe, d2wdthe2) COMPUTE_DTHE_C(t, dtdthe, d2tdthe2)
            COMPUTE_DTHE_C(c, dcdthe, d2cdthe2) COMPUTE_DTHE_C(cloud, dclouddthe, d2clouddthe2)
            COMPUTE_DTHE_C(ice, dicedthe, d2icedthe2) COMPUTE_DTHE_C(gr, dgdthe, d2gdthe2)
            COMPUTE_DTHE_C(co2, dcodthe, d2codthe2)
            COMPUTE_DTHE_C(tke, dtkedthe, d2tkedthe2) COMPUTE_DTHE_C(dis, ddisdthe, d2disdthe2)
            dpdthe   = -(-3.0*p_dyn.x[i][j][k] + 4.0*p_dyn.x[i][j-1][k] - p_dyn.x[i][j-2][k]) * inv_2dthe;
            dnuedthe = -(-3.0*nue.x[i][j][k]   + 4.0*nue.x[i][j-1][k]   - nue.x[i][j-2][k])   * inv_2dthe;
            the_flag = true;
            #undef COMPUTE_DTHE_C
        }
    }
    else if(j == 0){
        #define COMPUTE_DTHE_NP(FIELD, d1, d2) \
            d1 = (-3.0*FIELD.x[i][0][k] + 4.0*FIELD.x[i][1][k] - FIELD.x[i][2][k]) * inv_2dthe; \
            d2 = (FIELD.x[i][0][k] - 2.0*FIELD.x[i][1][k] + FIELD.x[i][2][k]) * inv_dthe2;

        COMPUTE_DTHE_NP(u, dudthe, d2udthe2) COMPUTE_DTHE_NP(v, dvdthe, d2vdthe2)
        COMPUTE_DTHE_NP(w, dwdthe, d2wdthe2) COMPUTE_DTHE_NP(t, dtdthe, d2tdthe2)
        COMPUTE_DTHE_NP(c, dcdthe, d2cdthe2) COMPUTE_DTHE_NP(cloud, dclouddthe, d2clouddthe2)
        COMPUTE_DTHE_NP(ice, dicedthe, d2icedthe2) COMPUTE_DTHE_NP(gr, dgdthe, d2gdthe2)
        COMPUTE_DTHE_NP(co2, dcodthe, d2codthe2)
        COMPUTE_DTHE_NP(tke, dtkedthe, d2tkedthe2) COMPUTE_DTHE_NP(dis, ddisdthe, d2disdthe2)
        dpdthe   = (-3.0*p_dyn.x[i][0][k] + 4.0*p_dyn.x[i][1][k] - p_dyn.x[i][2][k]) * inv_2dthe;
        dnuedthe = (-3.0*nue.x[i][0][k]   + 4.0*nue.x[i][1][k]   - nue.x[i][2][k])   * inv_2dthe;
        the_flag = true;
        #undef COMPUTE_DTHE_NP
    }
    else if(j == jm-1){
        #define COMPUTE_DTHE_SP(FIELD, d1, d2) \
            d1 = -(-3.0*FIELD.x[i][jm-1][k] + 4.0*FIELD.x[i][jm-2][k] - FIELD.x[i][jm-3][k]) * inv_2dthe; \
            d2 = -(FIELD.x[i][jm-1][k] - 2.0*FIELD.x[i][jm-2][k] + FIELD.x[i][jm-3][k]) * inv_dthe2;

        COMPUTE_DTHE_SP(u, dudthe, d2udthe2) COMPUTE_DTHE_SP(v, dvdthe, d2vdthe2)
        COMPUTE_DTHE_SP(w, dwdthe, d2wdthe2) COMPUTE_DTHE_SP(t, dtdthe, d2tdthe2)
        COMPUTE_DTHE_SP(c, dcdthe, d2cdthe2) COMPUTE_DTHE_SP(cloud, dclouddthe, d2clouddthe2)
        COMPUTE_DTHE_SP(ice, dicedthe, d2icedthe2) COMPUTE_DTHE_SP(gr, dgdthe, d2gdthe2)
        COMPUTE_DTHE_SP(co2, dcodthe, d2codthe2)
        COMPUTE_DTHE_SP(tke, dtkedthe, d2tkedthe2) COMPUTE_DTHE_SP(dis, ddisdthe, d2disdthe2)
        dpdthe   = -(-3.0*p_dyn.x[i][jm-1][k] + 4.0*p_dyn.x[i][jm-2][k] - p_dyn.x[i][jm-3][k]) * inv_2dthe;
        dnuedthe = -(-3.0*nue.x[i][jm-1][k]   + 4.0*nue.x[i][jm-2][k]   - nue.x[i][jm-3][k])   * inv_2dthe;
        the_flag = true;
        #undef COMPUTE_DTHE_SP
    }

    if(!the_flag){
        #define COMPUTE_DTHE_A(FIELD, d1, d2) \
            d1 = (FIELD.x[i][j+1][k] - FIELD.x[i][j-1][k]) * inv_2dthe; \
            d2 = (FIELD.x[i][j+1][k] - 2.0*FIELD.x[i][j][k] + FIELD.x[i][j-1][k]) * inv_dthe2;

        COMPUTE_DTHE_A(u, dudthe, d2udthe2) COMPUTE_DTHE_A(v, dvdthe, d2vdthe2)
        COMPUTE_DTHE_A(w, dwdthe, d2wdthe2) COMPUTE_DTHE_A(t, dtdthe, d2tdthe2)
        COMPUTE_DTHE_A(c, dcdthe, d2cdthe2) COMPUTE_DTHE_A(cloud, dclouddthe, d2clouddthe2)
        COMPUTE_DTHE_A(ice, dicedthe, d2icedthe2) COMPUTE_DTHE_A(gr, dgdthe, d2gdthe2)
        COMPUTE_DTHE_A(co2, dcodthe, d2codthe2)
        COMPUTE_DTHE_A(tke, dtkedthe, d2tkedthe2) COMPUTE_DTHE_A(dis, ddisdthe, d2disdthe2)
        dpdthe   = (p_dyn.x[i][j+1][k] - p_dyn.x[i][j-1][k]) * inv_2dthe;
        dnuedthe = (nue.x[i][j+1][k]   - nue.x[i][j-1][k])   * inv_2dthe;
        #undef COMPUTE_DTHE_A

        // Neumann BC for tke/dis at theta-direction land boundaries.
        // Air cells adjacent to land would read tke/dis=0 from the zeroed land
        // cell, creating a false large gradient that amplifies D_w cross-diffusion.
        // Replace the land neighbour's value with the current cell's (zero-gradient).
        if (!land_ijk) {
            const bool neumann_jm1 = (j > 0    && is_land(h, i, j-1, k));
            const bool neumann_jp1 = (j < jm-1 && is_land(h, i, j+1, k));
            if (neumann_jm1 || neumann_jp1) {
                const double tke_jm1_n = neumann_jm1 ? tke.x[i][j][k] : tke.x[i][j-1][k];
                const double tke_jp1_n = neumann_jp1 ? tke.x[i][j][k] : tke.x[i][j+1][k];
                const double dis_jm1_n = neumann_jm1 ? dis.x[i][j][k] : dis.x[i][j-1][k];
                const double dis_jp1_n = neumann_jp1 ? dis.x[i][j][k] : dis.x[i][j+1][k];
                dtkedthe   = (tke_jp1_n - tke_jm1_n) * inv_2dthe;
                d2tkedthe2 = (tke_jp1_n - 2.0*tke.x[i][j][k] + tke_jm1_n) * inv_dthe2;
                ddisdthe   = (dis_jp1_n - dis_jm1_n) * inv_2dthe;
                d2disdthe2 = (dis_jp1_n - 2.0*dis.x[i][j][k] + dis_jm1_n) * inv_dthe2;
            }
        }
    }


    // ===== PHI-direction derivatives =====
    if(k > 2 && k < km-3){
        if(land_ijk && is_air(h, i, j, k+1) && is_air(h, i, j, k+2)){
            #define COMPUTE_DPHI_B(FIELD, d1, d2) \
                d1 = (-3.0*FIELD.x[i][j][k] + 4.0*FIELD.x[i][j][k+1] - FIELD.x[i][j][k+2]) * inv_2dphi; \
                d2 = (FIELD.x[i][j][k] - 2.0*FIELD.x[i][j][k+1] + FIELD.x[i][j][k+2]) * inv_dphi2;

            COMPUTE_DPHI_B(u, dudphi, d2udphi2) COMPUTE_DPHI_B(v, dvdphi, d2vdphi2)
            COMPUTE_DPHI_B(w, dwdphi, d2wdphi2) COMPUTE_DPHI_B(t, dtdphi, d2tdphi2)
            COMPUTE_DPHI_B(c, dcdphi, d2cdphi2) COMPUTE_DPHI_B(cloud, dclouddphi, d2clouddphi2)
            COMPUTE_DPHI_B(ice, dicedphi, d2icedphi2) COMPUTE_DPHI_B(gr, dgdphi, d2gdphi2)
            COMPUTE_DPHI_B(co2, dcodphi, d2codphi2)
            COMPUTE_DPHI_B(tke, dtkedphi, d2tkedphi2) COMPUTE_DPHI_B(dis, ddisdphi, d2disdphi2)
            dpdphi   = (-3.0*p_dyn.x[i][j][k] + 4.0*p_dyn.x[i][j][k+1] - p_dyn.x[i][j][k+2]) * inv_2dphi;
            dnuedphi = (-3.0*nue.x[i][j][k]   + 4.0*nue.x[i][j][k+1]   - nue.x[i][j][k+2])   * inv_2dphi;
            phi_flag = true;
            #undef COMPUTE_DPHI_B
        }
        if(!phi_flag && land_ijk && is_air(h, i, j, k-1) && is_air(h, i, j, k-2)){
            #define COMPUTE_DPHI_C(FIELD, d1, d2) \
                d1 = -(-3.0*FIELD.x[i][j][k] + 4.0*FIELD.x[i][j][k-1] - FIELD.x[i][j][k-2]) * inv_2dphi; \
                d2 = -(FIELD.x[i][j][k] - 2.0*FIELD.x[i][j][k-1] + FIELD.x[i][j][k-2]) * inv_dphi2;

            COMPUTE_DPHI_C(u, dudphi, d2udphi2) COMPUTE_DPHI_C(v, dvdphi, d2vdphi2)
            COMPUTE_DPHI_C(w, dwdphi, d2wdphi2) COMPUTE_DPHI_C(t, dtdphi, d2tdphi2)
            COMPUTE_DPHI_C(c, dcdphi, d2cdphi2) COMPUTE_DPHI_C(cloud, dclouddphi, d2clouddphi2)
            COMPUTE_DPHI_C(ice, dicedphi, d2icedphi2) COMPUTE_DPHI_C(gr, dgdphi, d2gdphi2)
            COMPUTE_DPHI_C(co2, dcodphi, d2codphi2)
            COMPUTE_DPHI_C(tke, dtkedphi, d2tkedphi2) COMPUTE_DPHI_C(dis, ddisdphi, d2disdphi2)
            dpdphi   = -(-3.0*p_dyn.x[i][j][k] + 4.0*p_dyn.x[i][j][k-1] - p_dyn.x[i][j][k-2]) * inv_2dphi;
            dnuedphi = -(-3.0*nue.x[i][j][k]   + 4.0*nue.x[i][j][k-1]   - nue.x[i][j][k-2])   * inv_2dphi;
            phi_flag = true;
            #undef COMPUTE_DPHI_C
        }
    }
    else if(k == 0){
        #define COMPUTE_DPHI_W(FIELD, d1, d2) \
            d1 = (-3.0*FIELD.x[i][j][0] + 4.0*FIELD.x[i][j][1] - FIELD.x[i][j][2]) * inv_2dphi; \
            d2 = (FIELD.x[i][j][km-1] - 2.0*FIELD.x[i][j][km+2] + FIELD.x[i][j][km+3]) * inv_dphi2;

        COMPUTE_DPHI_W(u, dudphi, d2udphi2) COMPUTE_DPHI_W(v, dvdphi, d2vdphi2)
        COMPUTE_DPHI_W(w, dwdphi, d2wdphi2) COMPUTE_DPHI_W(t, dtdphi, d2tdphi2)
        COMPUTE_DPHI_W(c, dcdphi, d2cdphi2) COMPUTE_DPHI_W(cloud, dclouddphi, d2clouddphi2)
        COMPUTE_DPHI_W(ice, dicedphi, d2icedphi2) COMPUTE_DPHI_W(gr, dgdphi, d2gdphi2)
        COMPUTE_DPHI_W(co2, dcodphi, d2codphi2)
        COMPUTE_DPHI_W(tke, dtkedphi, d2tkedphi2) COMPUTE_DPHI_W(dis, ddisdphi, d2disdphi2)
        dpdphi   = (-3.0*p_dyn.x[i][j][0] + 4.0*p_dyn.x[i][j][1] - p_dyn.x[i][j][2]) * inv_2dphi;
        dnuedphi = (-3.0*nue.x[i][j][0]   + 4.0*nue.x[i][j][1]   - nue.x[i][j][2])   * inv_2dphi;
        phi_flag = true;
        #undef COMPUTE_DPHI_W
    }
    else if(k == km-1){
        #define COMPUTE_DPHI_E(FIELD, d1, d2) \
            d1 = -(-3.0*FIELD.x[i][j][km-1] + 4.0*FIELD.x[i][j][km-2] - FIELD.x[i][j][km-3]) * inv_2dphi; \
            d2 = -(FIELD.x[i][j][km-1] - 2.0*FIELD.x[i][j][km-2] + FIELD.x[i][j][km-3]) * inv_dphi2;

        COMPUTE_DPHI_E(u, dudphi, d2udphi2) COMPUTE_DPHI_E(v, dvdphi, d2vdphi2)
        COMPUTE_DPHI_E(w, dwdphi, d2wdphi2) COMPUTE_DPHI_E(t, dtdphi, d2tdphi2)
        COMPUTE_DPHI_E(c, dcdphi, d2cdphi2) COMPUTE_DPHI_E(cloud, dclouddphi, d2clouddphi2)
        COMPUTE_DPHI_E(ice, dicedphi, d2icedphi2) COMPUTE_DPHI_E(gr, dgdphi, d2gdphi2)
        COMPUTE_DPHI_E(co2, dcodphi, d2codphi2)
        COMPUTE_DPHI_E(tke, dtkedphi, d2tkedphi2) COMPUTE_DPHI_E(dis, ddisdphi, d2disdphi2)
        dpdphi   = -(-3.0*p_dyn.x[i][j][km-1] + 4.0*p_dyn.x[i][j][km-2] - p_dyn.x[i][j][km-3]) * inv_2dphi;
        dnuedphi = -(-3.0*nue.x[i][j][km-1]   + 4.0*nue.x[i][j][km-2]   - nue.x[i][j][km-3])   * inv_2dphi;
        phi_flag = true;
        #undef COMPUTE_DPHI_E
    }

    if(!phi_flag){
        #define COMPUTE_DPHI_A(FIELD, d1, d2) \
            d1 = (FIELD.x[i][j][k+1] - FIELD.x[i][j][k-1]) * inv_2dphi; \
            d2 = (FIELD.x[i][j][k+1] - 2.0*FIELD.x[i][j][k] + FIELD.x[i][j][k-1]) * inv_dphi2;

        COMPUTE_DPHI_A(u, dudphi, d2udphi2) COMPUTE_DPHI_A(v, dvdphi, d2vdphi2)
        COMPUTE_DPHI_A(w, dwdphi, d2wdphi2) COMPUTE_DPHI_A(t, dtdphi, d2tdphi2)
        COMPUTE_DPHI_A(c, dcdphi, d2cdphi2) COMPUTE_DPHI_A(cloud, dclouddphi, d2clouddphi2)
        COMPUTE_DPHI_A(ice, dicedphi, d2icedphi2) COMPUTE_DPHI_A(gr, dgdphi, d2gdphi2)
        COMPUTE_DPHI_A(co2, dcodphi, d2codphi2)
        COMPUTE_DPHI_A(tke, dtkedphi, d2tkedphi2) COMPUTE_DPHI_A(dis, ddisdphi, d2disdphi2)
        dpdphi   = (p_dyn.x[i][j][k+1] - p_dyn.x[i][j][k-1]) * inv_2dphi;
        dnuedphi = (nue.x[i][j][k+1]   - nue.x[i][j][k-1])   * inv_2dphi;
        #undef COMPUTE_DPHI_A

        // Neumann BC for tke/dis at phi-direction land boundaries.
        // Without this, an ocean cell adjacent to land reads tke/dis=0 from the
        // zeroed land cell — the 1/(r sinθ)² amplification then makes D_w
        // enormous near the poles, causing blow-up within a few iterations.
        if (!land_ijk) {
            const bool neumann_km1 = (k > 0    && is_land(h, i, j, k-1));
            const bool neumann_kp1 = (k < km-1 && is_land(h, i, j, k+1));
            if (neumann_km1 || neumann_kp1) {
                const double tke_km1_n = neumann_km1 ? tke.x[i][j][k] : tke.x[i][j][k-1];
                const double tke_kp1_n = neumann_kp1 ? tke.x[i][j][k] : tke.x[i][j][k+1];
                const double dis_km1_n = neumann_km1 ? dis.x[i][j][k] : dis.x[i][j][k-1];
                const double dis_kp1_n = neumann_kp1 ? dis.x[i][j][k] : dis.x[i][j][k+1];
                dtkedphi   = (tke_kp1_n - tke_km1_n) * inv_2dphi;
                d2tkedphi2 = (tke_kp1_n - 2.0*tke.x[i][j][k] + tke_km1_n) * inv_dphi2;
                ddisdphi   = (dis_kp1_n - dis_km1_n) * inv_2dphi;
                d2disdphi2 = (dis_kp1_n - 2.0*dis.x[i][j][k] + dis_km1_n) * inv_dphi2;
            }
        }
    }


    // ===== Physics: coefficients =====
    double buoyancy = 1.0;

    double coeff_buoy     = r_air * g;
    double coeff_energy_p = u_0 * u_0 / (cp_l * t_0);
    double coeff_energy   = 1.0 / (c_0 * cp_l * t_0);
    double coeff_u_p      = 0.5 * r_air * u_0 * u_0 / L_atm;
    double coeff_trans    = 1.0;
    double coeff_MC_vel   = L_atm / (u_0 * u_0);
    double coeff_MC_q     = L_atm / (u_0 * c_0);
    double coeff_MC_t     = L_atm / (u_0 * t_0);

    // Coriolis
    double coriolis = 1.0;
    double coeff_Coriolis = r_air * u_0 * omega;
    // Coriolis acceleration -2 Ω × v in (r, θ, φ) with θ = colatitude (the0 = 0
    // at N pole), v south-positive (= +ê_θ), w east-positive. TRADITIONAL
    // APPROXIMATION: keep only Ω_r = Ω cosθ; drop the non-traditional ("Eötvös")
    // terms from Ω_θ = -Ω sinθ:
    //   F_r = 0                          (was +2Ω sinθ w — dropped)
    //   F_θ = +2Ω cos(θ) w
    //   F_φ = -2Ω cos(θ) v               (the -2Ω sinθ u part dropped)
    // (coefficient 2 ≡ 2Ω in this file's Ω-based non-dim time.) The earlier form
    // had sin↔cos swapped + wrong signs (the active k_omega_SST bug); after the
    // sign fix the kept +2Ω sinθ w drove the eastward jets steadily upward under
    // the lagging pressure projection, so it is dropped here.
    // ATOM_CORIOLIS_NONTRAD (AtomUtils::coriolis_nontraditional, lib/Utils.h) — the SAME switch
    // the ocean reads. Default 0 keeps exactly what stood here; =1 restores the non-traditional
    // pair, which is what the hydrosphere used to carry unilaterally. The pair goes together:
    // dropping only one of the two injects energy.
    const double nontrad = AtomUtils::coriolis_nontraditional() ? 1.0 : 0.0;
    double coriolis_rad =  nontrad * 2.0 * sinthe * w_ijk;
    double coriolis_the =  2.0 * costhe * w_ijk;
    double coriolis_phi = -2.0 * (costhe * v_ijk + nontrad * sinthe * u_ijk);

    // NONDIM-CONVENTION FIX (2026-06-19). This turbulent path had written its forcing
    // terms (Coriolis "2"=2Ω, surf_drag/HS as k/omega, buoyancy g/(omega*L_atm)) in an
    // "Ω-based time" convention that is INCONSISTENT with its pressure/advection terms,
    // which carry NO dt (advective-time, identical to the laminar RHS_Atm.cpp). That made
    // the forcing ~5e5× too "hot" per step → polar vertical-velocity runaway once the
    // buoyancy was strong. The laminar path (turb_model="none") is the calibrated, STABLE
    // reference (proper closed cells, w~23 m/s, u settling). force_nd converts the clean
    // Ω-time coefficients to the laminar advective-time scaling so the WHOLE turbulent RHS
    // is consistent: turbulent path = laminar dynamics + eddy viscosity (in diffusion_vel_re).
    // TRADE-WIND FIX 2026-07-10: dropped the extra "* dt". force_nd scales ONLY
    // the Coriolis (coriolis_rad/the/phi below); RK4 already applies its own dt
    // (v += rhs*dt), so the old "* dt" made Coriolis enter as dt^2 -> ~1/dt (1e4
    // at dt=1e-4) too weak vs advection (which enters as dt). Result: Rossby
    // number ~1/dt, quasi-non-rotating -> Hadley overturning works but the
    // equatorward surface flow is never deflected westward -> NO easterly trades
    // -> westerly (super-rotating) surface wind everywhere. This is the SAME bug
    // the ocean RHS already fixed (project_hydro_coriolis_dt_scaling: dropped its
    // extra *dt). Coriolis is rotational (energy-conserving) with a tiny CFL
    // coefficient (2*force_nd*dt ~ 7e-7), so this does NOT re-trigger the thermal
    // polar-runaway the *dt was added to tame (that is buoyancy, scaled separately
    // by g*dt/u_0). Advective-time nondim: Coriolis coeff = 2 Omega L/u0.
    // ATM_CORIOLIS_SCALE IS GONE (2026-07-28). It was a multiplier for sweeping force_nd toward
    // the value that would make Coriolis balance horizontal advection, because force_nd was
    // ~16000x too weak and nobody knew where the factor lived. It is now known, and it factorises:
    //
    //     16000  =  397.5   x   40
    //               ^metric     ^this
    //
    // The 397.5 was the horizontal metric using rad.z ~ 1.5 as the planetary radius; that is fixed
    // and is now the default (ATM_METRIC_RADIUS). The 40 is here, and it is exactly 1/dr: force_nd
    // divided by L_atm = 400 m, while the gradients it must balance are expressed per rad.z unit,
    // and one rad.z unit is metricShellLength() ~ 16023 m, not 400 m (L_atm is the amplitude of
    // the exponential stretch, not a grid step — see metricShellLength in cAtmosphereModel.h).
    //
    // So the length-consistent coefficient is omega * metricShellLength() / u_0 = 1.46e-01 against
    // the old 3.65e-03, and the knob has nothing left to sweep. Tuning a multiplier was the right
    // move while the factor was unexplained; keeping it now would only hide the fix.
    //
    // This puts Coriolis on the SAME length as the two terms it has to balance — horizontal
    // advection and the pressure gradient — both of which take their length from the metric
    // (inv_rm) rather than from a coefficient. Geostrophy is that three-way balance, so making
    // this one consistent is not the piecemeal rescaling that ATJUP measured going wrong; it is
    // the term that was out of step with the other two.
    //
    // STILL ON L_atm, and therefore still 40x too weak: the Held-Suarez relaxation (k_T*L_atm/u_0),
    // the Rayleigh surface drag (rayleigh_kf*L_atm/u_0), coeff_MC_*, coeff_S, coeff_L and nue_max.
    // Those are damping and forcing terms rather than members of the geostrophic balance, and
    // fixing them is the remaining piece of one coherent nondimensionalisation.
    const double force_nd = omega * metricShellLength() / u_0;

    // Acting forces
    CoriolisForce.x[i][j][k] = coriolis * coeff_Coriolis
        * sqrt((pow(coriolis_rad, 2)
        + pow(coriolis_the, 2)
        + pow(coriolis_phi, 2)) / 3.0);
    BuoyancyForce.x[i][j][k] = buoyancy * coeff_buoy * (t.x[i][j][k] - t_ref_level[i]);
    PressureGradientForce.x[i][j][k] = -coeff_u_p
        * sqrt((pow(dpdr, 2)
        + pow(dpdthe * inv_rm, 2)
        + pow(dpdphi * inv_rmsinthe, 2)) / 3.0);

    double t_u    = t.x[i][j][k] * t_0;
    double E_Rain = SaturationH2O::saturationPressure(t_u);
    double E_Ice  = SaturationH2O::sublimationPressure(t_u);
    // Exact saturation mass fractions. The dilute ep*E/(p - E) is negative wherever
    // p_sat exceeds the local pressure — the whole column above ~140 km here — and these
    // two feed the 0.85*q thresholds that switch the latent-heat term on. A negative
    // threshold switches it on unconditionally.
    const double M_other_ijk = AtmMixture::M_nonwater(c.x[i][j][k], co2.x[i][j][k],
                                                      m_comp.M_bg);
    double q_Rain = SaturationH2O::saturationMassFraction(E_Rain, p_hydro.x[i][j][k],
                                                          M_other_ijk);
    double q_Ice  = SaturationH2O::saturationMassFraction(E_Ice,  p_hydro.x[i][j][k],
                                                          M_other_ijk);
    double Q_Latent_Ice = 0.0;

    double coeff_S = lamda * t_0 / L_atm;
    double coeff_L = r_air * c_0 * lv * u_0 / L_atm;

    if(c.x[i][j][k] >= 0.85 * q_Rain){
        Q_Latent.x[i][j][k] = u_ijk * dclouddr
            + v_ijk * inv_rm     * dclouddthe
            + w_ijk * inv_rmsinthe * dclouddphi;
    } else Q_Latent.x[i][j][k] = 0.0;

    if(c.x[i][j][k] >= 0.85 * q_Ice){
        Q_Latent_Ice = u_ijk * dicedr
            + v_ijk * inv_rm     * dicedthe
            + w_ijk * inv_rmsinthe * dicedphi;
    } else Q_Latent_Ice = 0.0;

    Q_Latent.x[i][j][k] = coeff_L * (Q_Latent.x[i][j][k] + Q_Latent_Ice);

    Q_Sensible.x[i][j][k] = coeff_S
        * (d2tdr2 + dtdr * 2.0 * inv_rm + d2tdthe2 * inv_rm2
        + dtdthe * costhe_inv_rm2sinthe + d2tdphi2 * inv_rm2sinthe2);

    if(land_ijk){
        BuoyancyForce.x[i][j][k]         = 0.0;
        PressureGradientForce.x[i][j][k] = 0.0;
        CoriolisForce.x[i][j][k]         = 0.0;
        Q_Latent.x[i][j][k]              = 0.0;
        Q_Sensible.x[i][j][k]            = 0.0;
    }


    // ===== Turbulence model: diffusion coefficients =====
    double diffusion_t_re   = 0.0;
    double diffusion_vel_re = 0.0;
    double diffusion_tke_re = 0.0;
    double diffusion_dis_re = 0.0;
    const double diff_prec_re_inv = 1.0 / (sc_WaterVapour * re_turb);
    const double diff_co2_re_inv  = 1.0 / (sc_CO2 * re_turb);

    double d_i       = 0.0;
    double d_i_mount = 0.0;
    int i_mount = i_topography[j][k];
    if(use_stretched_coordinate_system){
        d_i       = log(zeta * (rad.z[i]       + 1.0));
        d_i_mount = log(zeta * (rad.z[i_mount] + 1.0));
    } else {
        d_i       = (double)i;
        d_i_mount = (double)i_mount;
    }

    // ---- Refresh prod.x at every RK4 sub-stage ----
    // compute_sources() updates prod.x only at even iterations, so the stored
    // value is stale for intermediate RK4 stages.  Recompute it here from the
    // velocity derivatives already available above, scaled to match the
    // spherical-coord convention in compute_sources() (exp_rm / inv_rm /
    // inv_rmsinthe factors).  Constants mirror TurbulenceAtm private members.
    if(use_turbulence_model && !land_ijk){
        constexpr double C_nue_loc   = 0.028;    // mirrors TurbulenceAtm::C_nue
        constexpr double dis_min_loc = 1.0e-10;  // mirrors TurbulenceAtm::dis_min
        const double dudr_s   = dudr   * exp_rm;
        const double dvdr_s   = dvdr   * exp_rm;
        const double dwdr_s   = dwdr   * exp_rm;
        const double dudthe_s = dudthe * inv_rm;
        const double dvdthe_s = dvdthe * inv_rm;
        const double dwdthe_s = dwdthe * inv_rm;
        const double dudphi_s = dudphi * inv_rmsinthe;
        const double dvdphi_s = dvdphi * inv_rmsinthe;
        const double dwdphi_s = dwdphi * inv_rmsinthe;

        const double dis_here = std::max(dis.x[i][j][k], dis_min_loc);
        const double tke_here = std::max(tke.x[i][j][k], 0.0);
        const double nue_max  = 1000.0 / (u_0 * L_atm);
        double cnue = use_k_epsilon_turbulence_model
                    ? C_nue_loc * tke_here * tke_here / dis_here
                    : tke_here / dis_here;
        cnue = std::min(cnue, nue_max);

        const double der = 0.66667 * (dudr_s + dvdthe_s + dwdphi_s);

        prod.x[i][j][k] = std::max(0.0,
              (cnue * (2.0 * dudr_s   - der) - 0.66667 * tke.x[i][j][k]) * dudr_s
            + (cnue * (dudthe_s + dvdr_s))                               * dudthe_s
            + (cnue * (dudphi_s + dwdr_s))                               * dudphi_s
            + (cnue * (2.0 * dvdthe_s - der) - 0.66667 * tke.x[i][j][k]) * dvdthe_s
            + (cnue * (dvdr_s + dudthe_s))                               * dvdr_s
            + (cnue * (dvdphi_s + dwdthe_s))                             * dvdphi_s
            + (cnue * (2.0 * dwdphi_s - der) - 0.66667 * tke.x[i][j][k]) * dwdphi_s
            + (cnue * (dwdr_s + dudphi_s))                               * dwdr_s
            + (cnue * (dwdthe_s + dvdphi_s))                             * dwdthe_s);
    }

    if(use_k_epsilon_turbulence_model){                                 // KE_chien 1982
        const double sig_k   = 1.0;
        const double sig_w   = 1.3;
        const double C_eps_1 = 1.35;
        const double C_eps_2 = 1.80;
        // Diffusion coefficients: molecular background 1/re_turb plus turbulent nue*.
        // nue* = ν_T/(u_0·L_atm) enters directly — NOT divided by re_turb — because
        // the non-dimensionalisation already bakes the (u_0·L_atm) scale into nue*.
        // pr_turb (turbulent Prandtl number, configurable) replaces the laminar pr for
        // both the molecular background and the turbulent contribution to heat transport.
        {
            // Diffusion: molecular background 1/re_turb plus turbulent nue* directly.
            // nue* = ν_T/(u_0·L_atm) is already in the same non-dim units as the
            // diffusion coefficient — it must NOT be divided by re_turb again.
            diffusion_t_re   = (1.0/re_turb + nue.x[i][j][k]) / pr_turb;
            diffusion_vel_re = 1.0/re_turb + nue.x[i][j][k];
            diffusion_tke_re = 1.0/re_turb + nue.x[i][j][k] / sig_k;
            diffusion_dis_re = 1.0/re_turb + nue.x[i][j][k] / sig_w;
        }
        // Source terms — computed fresh from current tke/dis so every RK4 sub-stage
        // is consistent (compute_k_epsilon runs only at even iterations).
        // L_k / L_w wall-damping corrections (Chien 1982) are kept from the stored
        // tke_source / dis_source set by compute_k_epsilon; only the P–Y balance is
        // updated here so the damping is not lost.
        {
            const double tke_s = std::max(tke.x[i][j][k], 1.0e-10);
            const double dis_s = std::max(dis.x[i][j][k], 1.0e-10);
            const double P_k   = prod.x[i][j][k];
            const double Y_k   = dis_s;
            const double P_w   = C_eps_1 * (dis_s / tke_s) * P_k;
            const double Y_w   = C_eps_2 * dis_s * dis_s / tke_s;
            const double dis_src_max = 20.0 * dis_s * dis_s;
            const double tke_src_max = 20.0 * tke_s * dis_s;
            tke_source.x[i][j][k] = std::clamp(P_k - Y_k, -tke_src_max, tke_src_max);
            dis_source.x[i][j][k] = std::clamp(P_w - Y_w, -dis_src_max, dis_src_max);
        }
    }

    if(use_k_omega_turbulence_model){                                   // Wilcox 1988
        const double sig_k    = 0.6;
        const double sig_w    = 0.5;
        const double bet_star = 0.09;
        const double gam      = 0.52;
        const double bet_0    = 0.0708;
        {
            diffusion_t_re   = (1.0/re_turb + nue.x[i][j][k]) / pr_turb;
            diffusion_vel_re = 1.0/re_turb + nue.x[i][j][k];
            diffusion_tke_re = 1.0/re_turb + sig_k * nue.x[i][j][k];
            diffusion_dis_re = 1.0/re_turb + sig_w * nue.x[i][j][k];
        }
        // Source terms — fresh tke/dis, stored prod
        {
            const double tke_s  = std::max(tke.x[i][j][k], 1.0e-10);
            const double dis_s  = std::max(dis.x[i][j][k], 1.0e-10);
            const double P_k    = std::min(prod.x[i][j][k], 20.0 * bet_star * tke_s * dis_s);
            const double Y_k    = bet_star * tke_s * dis_s;
            const double P_w    = gam * (dis_s / tke_s) * P_k;
            const double Y_w    = bet_0 * dis_s * dis_s;
            const double cross  = dtkedr   * ddisdr
                                + dtkedthe * inv_rm       * ddisdthe  * inv_rm
                                + dtkedphi * inv_rmsinthe * ddisdphi  * inv_rmsinthe;
            const double sig_d  = (cross <= 0.0) ? 0.0 : 0.125;
            const double D_w    = sig_d * cross / dis_s;
            const double dis_src_max = 20.0 * dis_s * dis_s;
            const double tke_src_max = 20.0 * tke_s * dis_s;
            tke_source.x[i][j][k] = std::clamp(P_k - Y_k, -tke_src_max, tke_src_max);
            dis_source.x[i][j][k] = std::clamp(P_w - Y_w + D_w, -dis_src_max, dis_src_max);
        }
    }

    if(use_k_omega_SST_turbulence_model){                               // Menter 1994
        double bet_star = 0.09;    // β* destruction coefficient — NOT C_μ=0.028
        double sig_k1   = 1.176;
        double sig_k2   = 1.0;
        double sig_w1   = 2.0;
        double sig_w2   = 1.168;

        int j_beg = 0, k_beg = 0;
        double d_j_mount = 0.0, d_k_mount = 0.0;
        double d_j_south = 0.0, d_j_north = 0.0;
        double d_k_east  = 0.0, d_k_west  = 0.0;
        double d_j = (double)j;
        double d_k = (double)k;

        if(land_ijk && is_water(h, i, j+1, k))   j_beg = j;
        for(int jj = j_beg; jj < jm-1; jj++)
            if(is_land(h, i, jj, k) || jj == jm-1) d_j_mount = jj;
        d_j_south = fabs(d_j - d_j_mount);

        if(land_ijk && is_water(h, i, j-1, k))   j_beg = j;
        for(int jj = jm-1; jj > j_beg; jj--)
            if(is_land(h, i, jj, k) || jj == 0) d_j_mount = jj;
        d_j_north = fabs(d_j - d_j_mount);

        if(land_ijk && is_water(h, i, j, k+1))   k_beg = k;
        for(int kk = k_beg; kk < km-1; kk++)
            if(is_land(h, i, j, kk) || kk == km-1) d_k_mount = kk;
        d_k_east = fabs(d_k - d_k_mount);

        if(land_ijk && is_water(h, i, j, k-1))   k_beg = k;
        for(int kk = km-1; kk > k_beg; kk--)
            if(is_land(h, i, j, kk) || kk == 0) d_k_mount = kk;
        d_k_west = fabs(d_k - d_k_mount);

        double d_close = d_i - d_i_mount;
        if(d_j_south <= d_close) d_close = d_j_south;
        if(d_j_north <= d_close) d_close = d_j_north;
        if(d_k_east  <= d_close) d_close = d_k_east;
        if(d_k_west  <= d_close) d_close = d_k_west;

        double CD_kw = std::max(
            (2.0 * sig_w2) / dis.x[i][j][k]
            * (dtkedr   * ddisdr
            + dtkedthe * inv_rm      * ddisdthe  * inv_rm
            + dtkedphi * inv_rmsinthe * ddisdphi * inv_rmsinthe),
            1.0e-20);

        // Menter (1994): arg1 uses molecular viscosity ν_air, not turbulent ν_T.
        // nue_air_nd = ν_air / (u_0·L_atm) is the non-dimensional molecular viscosity
        // (same constant as TurbulenceAtm uses for nue_air_nd).
        constexpr double nue_air_phys = 1.5e-5;   // [m²/s]
        const double nue_air_nd = nue_air_phys / (u_0 * L_atm);
        double arg1 = std::min(std::max(
            sqrt(tke.x[i][j][k]) / (bet_star * dis.x[i][j][k] * d_close),
            (500.0 * nue_air_nd) / (d_close * d_close * dis.x[i][j][k])),
            (4.0 * sig_w2 * tke.x[i][j][k]) / (CD_kw * d_close * d_close));
        double F1 = tanh(pow(arg1, 4));

        {
            diffusion_t_re   = (1.0/re_turb + nue.x[i][j][k]) / pr_turb;
            diffusion_vel_re =  1.0/re_turb + nue.x[i][j][k];
            diffusion_tke_re =  1.0/re_turb + nue.x[i][j][k] / blend(sig_k1, sig_k2, F1);
            diffusion_dis_re =  1.0/re_turb + nue.x[i][j][k] / blend(sig_w1, sig_w2, F1);
        }

        // Source terms — computed fresh from current tke/dis at every RK4 sub-stage.
        // compute_k_omega_SST runs only at even iterations and would otherwise leave
        // tke_source/dis_source one outer iteration behind.
        {
            const double bet1   = 0.0333;
            const double bet2   = 0.0368;
            const double gam1   = 0.413;
            const double gam2   = 0.2;
            const double tke_s  = std::max(tke.x[i][j][k], 1.0e-10);
            const double dis_s  = std::max(dis.x[i][j][k], 1.0e-10);
            const double P_k    = std::min(prod.x[i][j][k], 20.0 * bet_star * tke_s * dis_s);
            const double Y_k    = bet_star * tke_s * dis_s;
            const double P_w    = blend(gam1, gam2, F1) * P_k * dis_s / tke_s;
            const double Y_w    = blend(bet1,  bet2, F1) * dis_s * dis_s;
            const double D_w    = 2.0 * (1.0 - F1) * sig_w2 / dis_s
                                * (dtkedr   * ddisdr
                                + dtkedthe * inv_rm       * ddisdthe  * inv_rm
                                + dtkedphi * inv_rmsinthe * ddisdphi  * inv_rmsinthe);
            const double dis_src_max = 20.0 * dis_s * dis_s;
            const double tke_src_max = 20.0 * tke_s * dis_s;
            tke_source.x[i][j][k] = std::clamp(P_k - Y_k, -tke_src_max, tke_src_max);
            dis_source.x[i][j][k] = std::clamp(P_w - Y_w + D_w, -dis_src_max, dis_src_max);
        }
    }

    if(!use_turbulence_model){
        diffusion_t_re        = 1.0 / (re_turb * pr_turb);
        diffusion_vel_re      = 1.0 / re_turb;
        diffusion_tke_re      = 0.0;
        diffusion_dis_re      = 0.0;
        tke.x[i][j][k]        = 0.0;
        dis.x[i][j][k]        = 0.0;
        nue.x[i][j][k]        = 0.0;
        prod.x[i][j][k]       = 0.0;
        tke_source.x[i][j][k] = 0.0;
        dis_source.x[i][j][k] = 0.0;
        tken.x[i][j][k]       = 0.0;
        disn.x[i][j][k]       = 0.0;
    }

    // Inviscid spin-up: ramp the diffusion (molecular background + eddy nue) from 0, mirroring
    // the former laminar RHS_Atm.cpp. This lets the turbulent path also serve the inviscid phase
    // (diffusion_ramp=0 => Euler) and the laminar case (use_turbulence_model=false => molecular
    // only, above), so it is now the single dynamical core. No-op post-spinup (diffusion_ramp=1).
    diffusion_t_re   *= diffusion_ramp;
    diffusion_vel_re *= diffusion_ramp;
    diffusion_tke_re *= diffusion_ramp;
    diffusion_dis_re *= diffusion_ramp;

    // ===== Transport terms (advection) =====
    double u_exp   = u_ijk * exp_rm;
    double v_invrm = v_ijk * inv_rm;
    double w_invrs = w_ijk * inv_rmsinthe;

    double pressure_t = coeff_energy_p * (u_exp * dpdr + v_invrm * dpdthe + w_invrs * dpdphi);

    // Positivity-preserving (minmod-limited) advective gradients for the moisture
    // scalars c/cloud/ice. Used ONLY in transport_{c,cloud,ice}; diffusion and
    // Q_Latent keep the centered stencil. Centered advection overshoots at the sharp
    // coastal cloud/vapour gradient under the strong coastal jet, driving cloud/c
    // negative (then silently clamped by storeIntermediateData3D — a mass leak + 2Δt
    // sawtooth). minmod→0 at a local extremum, so advection can no longer create the
    // undershoot. Only the interior central stencils are limited (!flag); the one-sided
    // boundary/pole/seam stencils fall back to the centered value already computed.
    const double inv_dr   = 2.0 * inv_2dr;
    const double inv_dthe = 2.0 * inv_2dthe;
    const double inv_dphi = 2.0 * inv_2dphi;

    double dcdr_adv = dcdr, dclouddr_adv = dclouddr, dicedr_adv = dicedr;
    double dcdthe_adv = dcdthe, dclouddthe_adv = dclouddthe, dicedthe_adv = dicedthe;
    double dcdphi_adv = dcdphi, dclouddphi_adv = dclouddphi, dicedphi_adv = dicedphi;
    // Temperature gets the same minmod-limited advective gradient as the moisture
    // scalars. A term-by-term rhs_t budget at the steep BC/Alaska coast (57°N, the low
    // ocean-air column at the foot of the cliff) showed the growing cold bubble is driven
    // entirely by centered -transport_t (~ -11) running as a 2Δt odd/even checkerboard,
    // with diffusion/latent/MC all negligible — the temperature analogue of the c/cloud
    // overshoot. minmod→0 at the vertical local extremum (i=6 cold between warmer i=5/i=7)
    // so advection can no longer deepen the cold pool. diffusion_t keeps its centered stencil.
    double dtdr_adv = dtdr, dtdthe_adv = dtdthe, dtdphi_adv = dtdphi;
    // Velocity gets the same minmod-limited advective gradient (transport_{u,v,w} only;
    // diffusion/pressure/Coriolis/metric keep centered stencils). A momentum budget at the
    // Greenwich-seam blow-up cell (29°N/1°E, k=1, upper troposphere) showed an explosive
    // 2Δ grid-scale velocity checkerboard driven ENTIRELY by the centered zonal self-
    // advection -w·∂w/∂φ (transport_w's phi term), pressure/Coriolis/diffusion negligible —
    // the velocity analogue of the coastal-t and cloud overshoots. minmod→0 at a 2Δ extremum
    // kills the odd-even mode; smooth multi-cell jets/vortices keep ~centered values.
    double dudr_adv = dudr, dudthe_adv = dudthe, dudphi_adv = dudphi;
    double dvdr_adv = dvdr, dvdthe_adv = dvdthe, dvdphi_adv = dvdphi;
    double dwdr_adv = dwdr, dwdthe_adv = dwdthe, dwdphi_adv = dwdphi;

    if(!r_flag){
        dtdr_adv     = minmod(t.x[i][j][k]     - t.x[i-1][j][k],     t.x[i+1][j][k]     - t.x[i][j][k])     * inv_dr;
        dudr_adv     = minmod(u.x[i][j][k]     - u.x[i-1][j][k],     u.x[i+1][j][k]     - u.x[i][j][k])     * inv_dr;
        dvdr_adv     = minmod(v.x[i][j][k]     - v.x[i-1][j][k],     v.x[i+1][j][k]     - v.x[i][j][k])     * inv_dr;
        dwdr_adv     = minmod(w.x[i][j][k]     - w.x[i-1][j][k],     w.x[i+1][j][k]     - w.x[i][j][k])     * inv_dr;
        dcdr_adv     = minmod(c.x[i][j][k]     - c.x[i-1][j][k],     c.x[i+1][j][k]     - c.x[i][j][k])     * inv_dr;
        dclouddr_adv = minmod(cloud.x[i][j][k] - cloud.x[i-1][j][k], cloud.x[i+1][j][k] - cloud.x[i][j][k]) * inv_dr;
        dicedr_adv   = minmod(ice.x[i][j][k]   - ice.x[i-1][j][k],   ice.x[i+1][j][k]   - ice.x[i][j][k])   * inv_dr;
    }
    if(!the_flag){
        dtdthe_adv     = minmod(t.x[i][j][k]     - t.x[i][j-1][k],     t.x[i][j+1][k]     - t.x[i][j][k])     * inv_dthe;
        dudthe_adv     = minmod(u.x[i][j][k]     - u.x[i][j-1][k],     u.x[i][j+1][k]     - u.x[i][j][k])     * inv_dthe;
        dvdthe_adv     = minmod(v.x[i][j][k]     - v.x[i][j-1][k],     v.x[i][j+1][k]     - v.x[i][j][k])     * inv_dthe;
        dwdthe_adv     = minmod(w.x[i][j][k]     - w.x[i][j-1][k],     w.x[i][j+1][k]     - w.x[i][j][k])     * inv_dthe;
        dcdthe_adv     = minmod(c.x[i][j][k]     - c.x[i][j-1][k],     c.x[i][j+1][k]     - c.x[i][j][k])     * inv_dthe;
        dclouddthe_adv = minmod(cloud.x[i][j][k] - cloud.x[i][j-1][k], cloud.x[i][j+1][k] - cloud.x[i][j][k]) * inv_dthe;
        dicedthe_adv   = minmod(ice.x[i][j][k]   - ice.x[i][j-1][k],   ice.x[i][j+1][k]   - ice.x[i][j][k])   * inv_dthe;
    }
    if(!phi_flag){
        dtdphi_adv     = minmod(t.x[i][j][k]     - t.x[i][j][k-1],     t.x[i][j][k+1]     - t.x[i][j][k])     * inv_dphi;
        dudphi_adv     = minmod(u.x[i][j][k]     - u.x[i][j][k-1],     u.x[i][j][k+1]     - u.x[i][j][k])     * inv_dphi;
        dvdphi_adv     = minmod(v.x[i][j][k]     - v.x[i][j][k-1],     v.x[i][j][k+1]     - v.x[i][j][k])     * inv_dphi;
        dwdphi_adv     = minmod(w.x[i][j][k]     - w.x[i][j][k-1],     w.x[i][j][k+1]     - w.x[i][j][k])     * inv_dphi;
        dcdphi_adv     = minmod(c.x[i][j][k]     - c.x[i][j][k-1],     c.x[i][j][k+1]     - c.x[i][j][k])     * inv_dphi;
        dclouddphi_adv = minmod(cloud.x[i][j][k] - cloud.x[i][j][k-1], cloud.x[i][j][k+1] - cloud.x[i][j][k]) * inv_dphi;
        dicedphi_adv   = minmod(ice.x[i][j][k]   - ice.x[i][j][k-1],   ice.x[i][j][k+1]   - ice.x[i][j][k])   * inv_dphi;
    }

    double transport_u     = u_exp * dudr_adv     + v_invrm * dudthe_adv     + w_invrs * dudphi_adv;
    double transport_v     = u_exp * dvdr_adv     + v_invrm * dvdthe_adv     + w_invrs * dvdphi_adv;
    double transport_w     = u_exp * dwdr_adv     + v_invrm * dwdthe_adv     + w_invrs * dwdphi_adv;

    // ATOM_METRIC_CURVATURE (AtomUtils, lib/Utils.h) — see there for the form and for why these
    // must not be switched on without ATM_METRIC_RADIUS. Signs follow from rhs_x = ... - transport_x,
    // and the expressions are ATJUP's (RHS_Jup_Turb.cpp:405-411) term for term.
    if(AtomUtils::metric_curvature()){
        const double cotanthe = costhe / sinthe;
        transport_u += -(v_ijk * v_ijk + w_ijk * w_ijk) * inv_rm;
        transport_v +=  (u_ijk * v_ijk - w_ijk * w_ijk * cotanthe) * inv_rm;
        transport_w +=  (w_ijk * u_ijk + v_ijk * w_ijk * cotanthe) * inv_rm;
    }
    double transport_t     = u_exp * dtdr_adv     + v_invrm * dtdthe_adv     + w_invrs * dtdphi_adv;
    double transport_c     = u_exp * dcdr_adv     + v_invrm * dcdthe_adv     + w_invrs * dcdphi_adv;
    double transport_cloud = u_exp * dclouddr_adv + v_invrm * dclouddthe_adv + w_invrs * dclouddphi_adv;
    double transport_ice   = u_exp * dicedr_adv   + v_invrm * dicedthe_adv   + w_invrs * dicedphi_adv;
    double transport_g     = u_exp * dgdr     + v_invrm * dgdthe     + w_invrs * dgdphi;
    double transport_co2   = u_exp * dcodr    + v_invrm * dcodthe    + w_invrs * dcodphi;
    double transport_tke   = u_exp * dtkedr   + v_invrm * dtkedthe   + w_invrs * dtkedphi;
    double transport_dis   = u_exp * ddisdr   + v_invrm * ddisdthe   + w_invrs * ddisdphi;

    // Advective form only (v·∇k, v·∇ω) — matches the standard k-ε / k-ω derivations
    // (Pope, Wilcox, Menter).  A previous version added the conservative-form
    // correction `+tke·∇·v` and `+dis·∇·v`; in an incompressible flow ∇·v=0 and the
    // two forms agree, but ATOM is compressible and the term produced exponential
    // growth of tke in any region with persistent negative divergence (top of an
    // updraft, tropopause spreading).  Observed: max tke = 6.5e98 at 64°S 63°W
    // 14566 m by iter 100 once the global circulation activated div_vel.  tke's sink
    // Y_k = β*·tke·dis is linear in tke so it can never catch the runaway; dis was
    // bounded only because Y_w = β₀·dis² is quadratic and self-limiting.


    // NOTE (2026-06-09): a background velocity-diffusion floor (nue_floor=1e-3) was tested
    // here and REVERTED. It gave negligible benefit against the iter-491 NaN (dif only rose
    // 0.001→0.009, crash unchanged — the cause is the pressure/divergence loop, not a
    // diffusion deficit) AND it adds 1e-3 isotropic diffusion wherever turbulent nue<1e-3
    // (i.e. most of the domain, including i=0), which smooths the emergent i=0 vortices that
    // are a healthy-run validation signal ([[project_emergence_validation]]). Net negative.
    // See [[project_upper_velocity_secular_growth]].

    // ===== Diffusion terms =====
    double two_over_rm_exp  = 2.0 * inv_rm * exp_rm;
    double cos_rm2sin       = costhe_inv_rm2sinthe;
    // v_metric = (1/sin²θ)/r² = (1 + cot²θ)/r² = (1 + cos²θ/sin²θ)/r². 2026-06-19: the
    // turbulent path had costhe (not costhe²), the pre-fix form of the
    // project_vmetric_antidiffusion bug that was corrected in RHS_Atm.cpp/RHS_Hyd.cpp but
    // missed here — it over-damped v,w in the NH and anti-diffused poleward of ~38°S.
    double v_metric         = (1.0 + costhe * costhe / sinthe2) * inv_rm2;

    double diffusion_t = (d2tdr2 * exp_2_rm + dtdr * two_over_rm_exp
        + d2tdthe2 * inv_rm2 + dtdthe * cos_rm2sin
        + d2tdphi2 * inv_rm2sinthe2) * diffusion_t_re;

    double diffusion_u = (d2udr2 * exp_2_rm + 2.0 * u_ijk * inv_rm2
        + d2udthe2 * inv_rm2 + 4.0 * dudr * inv_rm * exp_rm
        + dudthe * cos_rm2sin + d2udphi2 * inv_rm2sinthe2) * diffusion_vel_re;

    double diffusion_v = (d2vdr2 * exp_2_rm + dvdr * two_over_rm_exp
        + d2vdthe2 * inv_rm2 + dvdthe * cos_rm2sin
        - v_metric * v_ijk + d2vdphi2 * inv_rm2sinthe2
        + 2.0 * dudthe * inv_rm2
        - dwdphi * 2.0 * costhe * inv_rm2sinthe2) * diffusion_vel_re;

    double diffusion_w = (d2wdr2 * exp_2_rm + dwdr * two_over_rm_exp
        + d2wdthe2 * inv_rm2 + dwdthe * cos_rm2sin
        - v_metric * w_ijk + d2wdphi2 * inv_rm2sinthe2
        + 2.0 * dudphi * inv_rm2sinthe
        + dvdphi * 2.0 * costhe * inv_rm2sinthe2) * diffusion_vel_re;

    // ATM_TRACER_DIFF_FLUX (default 0 = bit-identical) — the DIFFUSIVE half of the tracer
    // mass budget, which is a different question from the advective half and is not
    // answered by it.
    //
    // For a mass fraction the conservative diffusive flux is d(rho q)/dt = div(rho K grad q),
    // i.e. dq/dt = K lap(q) + K grad(ln rho) . grad(q). The second term is missing here: the
    // operators below are the plain Laplacian, which conserves the integral of q dV and NOT
    // of rho q dV. That is the same failure the advective form has — but with a crucial
    // difference. The advective one is repaired by continuity, so it belongs in the flow
    // (README item 15). This one is not repaired by anything: it is a term the tracer
    // equation is simply missing, and it survives however good the velocity field is.
    //
    // It is also safe in the way the advective flux-form correction was not. It is
    // proportional to grad(q), so a uniform tracer stays uniform and CO2 stays well mixed —
    // the test that ruled the other correction out.
    //
    // rho_bar depends on r only, so this is one term on the radial derivative, and both
    // radial factors are exp_rm (dln(rho)/dr_phys and dq/dr_phys) — hence exp_2_rm, the
    // same factor the radial Laplacian term carries. The centred dcdr is used, not the
    // minmod-limited dcdr_adv: this is a diffusive flux, not advection.
    static const bool tracer_diff_flux = [](){ const char* e = getenv("ATM_TRACER_DIFF_FLUX");
                                               return e ? (atof(e) != 0.0) : false; }();
    const double dlnrho_dr_i = (tracer_diff_flux && (int)m_dlnrho_dr.size() == im)
                             ? m_dlnrho_dr[i] * exp_2_rm : 0.0;

    double diffusion_c = (d2cdr2 * exp_2_rm + dcdr * two_over_rm_exp
        + d2cdthe2 * inv_rm2 + dcdthe * cos_rm2sin
        + d2cdphi2 * inv_rm2sinthe2
        + dcdr * dlnrho_dr_i) * diff_prec_re_inv;

    double diffusion_cloud = (d2clouddr2 * exp_2_rm + dclouddr * two_over_rm_exp
        + d2clouddthe2 * inv_rm2 + dclouddthe * cos_rm2sin
        + d2clouddphi2 * inv_rm2sinthe2
        + dclouddr * dlnrho_dr_i) * diff_prec_re_inv;

    double diffusion_ice = (d2icedr2 * exp_2_rm + dicedr * two_over_rm_exp
        + d2icedthe2 * inv_rm2 + dicedthe * cos_rm2sin
        + d2icedphi2 * inv_rm2sinthe2
        + dicedr * dlnrho_dr_i) * diff_prec_re_inv;

    double diffusion_g = (d2gdr2 * exp_2_rm + dgdr * two_over_rm_exp
        + d2gdthe2 * inv_rm2 + dgdthe * cos_rm2sin
        + d2gdphi2 * inv_rm2sinthe2
        + dgdr * dlnrho_dr_i) * diff_prec_re_inv;

    double diffusion_co2 = (d2codr2 * exp_2_rm + dcodr * two_over_rm_exp
        + d2codthe2 * inv_rm2 + dcodthe * cos_rm2sin
        + d2codphi2 * inv_rm2sinthe2
        + dcodr * dlnrho_dr_i) * diff_co2_re_inv;

    double diffusion_tke = (d2tkedr2 * exp_2_rm + dtkedr * two_over_rm_exp
        + d2tkedthe2 * inv_rm2 + dtkedthe * cos_rm2sin
        + d2tkedphi2 * inv_rm2sinthe2) * diffusion_tke_re;

    double diffusion_dis = (d2disdr2 * exp_2_rm + ddisdr * two_over_rm_exp
        + d2disdthe2 * inv_rm2 + ddisdthe * cos_rm2sin
        + d2disdphi2 * inv_rm2sinthe2) * diffusion_dis_re;


    // ===== RHS assembly =====
    double dpdr_exp      = dpdr   * exp_rm;
    double dpdthe_invrm  = dpdthe * inv_rm;
    double dpdphi_invrs  = dpdphi * inv_rmsinthe;

    // Latent heating: signed microphysical mass-exchange × specific latent heat ×
    // density.  S_c+S_r (vapour↔liquid, factor lv), S_i+S_s+S_g (vapour↔ice/solid,
    // factor ls).  The laminar path (RHS_Atm.cpp) uses the same form.
    //
    // A previous version added `+ Q_Latent.x[i][j][k] / coeff_L` here; that was
    // an unsigned *diagnostic* magnitude written by ThermoAtm.latentSensibleHeat()
    // for ParaView output (= lv·|∇q_sat|, always positive).  Because ThermoAtm
    // runs only every other iter (moist_stride=2) while RHS reads Q_Latent every
    // iter, it produced a 2Δt sawtooth heating source — huge on iters where
    // ThermoAtm had just overwritten the field, near-zero on alternate iters —
    // that drove the iter-358 NaN at 62°N upper-troposphere
    // (see [[project-iter358-62n-dry-nan]]).  The signed S_*·lv/ls terms above
    // already account for the latent heat; the diagnostic term was a duplicate.
    rhs_t.x[i][j][k] = pressure_t - transport_t + diffusion_t
        + coeff_MC_t * MC_t.x[i][j][k]
        + coeff_energy * (S_c.x[i][j][k] + S_r.x[i][j][k]) * lv * r_humid.x[i][j][k]
        + coeff_energy * (S_i.x[i][j][k] + S_s.x[i][j][k] + S_g.x[i][j][k]) * ls * r_humid.x[i][j][k];

    // ----- Held-Suarez Newtonian thermal relaxation (PROTOTYPE, turbulent path) -----
    // Supplies the maintained APE the general circulation otherwise lacks (rhs_t had no
    // radiative forcing, so a forced-equilibrium IC like the Grotjahn mean cells spins down).
    // Relax T toward the radiative-equilibrium target t_eq (snapshot of the initial Scotese+
    // topography+land-sea field). Held-Suarez (1994) rate k_T = k_a + (k_s-k_a)*
    // max(0,(sigma-sigma_b)/(1-sigma_b))*cos^4(lat); sigma = p/p_surf, cos(lat)=sinthe.
    // Non-dimensionalised in advective time (k*L_atm/u_0*dt) to match the laminar RHS_Atm.cpp
    // and the rest of this path's forcing — see the force_nd note above (2026-06-19 fix).
    constexpr bool HELD_SUAREZ_RELAX = true;          // prototype toggle (keep in sync with RHS_Atm.cpp)
    if (HELD_SUAREZ_RELAX) {
        constexpr double k_a     = 1.0 / ( 4.0 * 86400.0);   // free-atmosphere relax rate [1/s] (10x H-S: matches ATOM's 10x friction scaling)
        constexpr double k_s     = 1.0 / ( 0.4 * 86400.0);   // near-surface relax rate   [1/s] (10x H-S)
        constexpr double sigma_b = 0.7;
        const double p_surf = p_stat.x[i_topography[j][k]][j][k];
        const double sigma  = (p_surf > 0.0) ? p_stat.x[i][j][k] / p_surf : 0.0;
        double bl = (sigma - sigma_b) / (1.0 - sigma_b);
        if (bl < 0.0) bl = 0.0;
        const double cos4lat = sinthe2 * sinthe2;            // cos^4(lat) = sin^4(colatitude)
        const double k_T = k_a + (k_s - k_a) * bl * cos4lat; // [1/s]
        // HS-dt TEST (2026-07-04): removed the extra *dt. Advective-time nondim of a Newtonian
        // relaxation is k_T*L_atm/u_0 (dimensionless), integrated by RK4's own *dt like every
        // other rhs_t term (advection/diffusion/latent, all dt^1). The extra *dt made HS enter
        // as dt^2 = ~1e-4 too weak -> inert (t_eq shift never reached T; CO2 forcing invisible).
        rhs_t.x[i][j][k] -= (k_T * L_atm / u_0) * (t.x[i][j][k] - t_eq.x[i][j][k]);
    }


    // Boussinesq buoyancy body force, in the laminar RHS_Atm.cpp advective-time form
    // g*dt/u_0 (the calibrated, STABLE reference). The turbulent path originally had the
    // bare anomaly (coeff 1.0). An intermediate "fix" used g/(omega*L_atm)≈336 to match the
    // buoyancy/Coriolis RATIO, but that combined with this path's over-hot Ω-time Coriolis to
    // drive a polar vertical runaway — see the force_nd note above. Using g*dt/u_0 (and
    // force_nd on Coriolis) makes the whole RHS share the laminar scaling.
    rhs_u.x[i][j][k] = -dpdr_exp - transport_u + diffusion_u
        + buoyancy_ramp * buoyancy * g * dt / u_0 * (t.x[i][j][k] - t_ref_level[i])
        + coriolis * force_nd * coriolis_rad;

    // ---- RADIAL momentum-budget term capture (checkpoint iters only) ----
    // Mirror of the vbud_*/wbud_* blocks below, for the VERTICAL wind. Ported from ATHAD
    // 2026-08-18 (its item 42). This is the component that had no instrument: ATHAD's
    // item 28 spurious 293-rms radial acceleration was invisible because nothing
    // decomposed rhs_u, and Psi is built from v alone so it cannot see a radial failure.
    // Note what the split makes checkable at a glance — with the default switches
    // coriolis_nontraditional() is false, so ubud_cor should be identically zero, and
    // metric_curvature() is false too; if either is nonzero, a term this fork's CLAUDE.md
    // assumes inert is live. ubud_buoy isolates ATHAD item 34's remaining extra-*dt
    // suspect, the half NOT repaired by the surf_drag fix above.
    if(ubudget_capture){
        ubud_pgf.x[i][j][k]  = -dpdr_exp;
        ubud_cor.x[i][j][k]  =  coriolis * force_nd * coriolis_rad;
        ubud_advv.x[i][j][k] = -(u_exp * dudr_adv);
        ubud_advh.x[i][j][k] = -(v_invrm * dudthe_adv + w_invrs * dudphi_adv);
        ubud_diff.x[i][j][k] =  diffusion_u;
        ubud_buoy.x[i][j][k] =  buoyancy_ramp * buoyancy * g * dt / u_0
                             * (t.x[i][j][k] - t_ref_level[i]);
    }

    // ----- Near-surface Rayleigh (boundary-layer) drag on the horizontal wind -----
    // See RHS_Atm.cpp for the rationale: the free-slip wall (bcSolidGround) + init-only
    // bcVelSurfSur leave the near-surface tangential wind with NO momentum sink, so the
    // coastal-jet v,w accelerate unopposed → the recurring i=1 CFL blow-up at steep
    // coasts (Antarctic Peninsula / Baffin / Norway / Sea of Okhotsk). Held-Suarez-style
    // linear drag on v,w only (NOT radial u), full strength at the first air cell above
    // the LOCAL surface (i_topography[j][k]) and ramping to zero over the boundary layer.
    // rayleigh_kf and drag_n_layers are CONFIG PARAMETERS now (ATHAD README item 44, ported
    // 2026-08-18). The Earth calibration comments that used to live here are in param.py,
    // where the values are — including the reason a CELL COUNT is not a depth: 5 cells is
    // 236 m on ATOM_Precipitation's grid and 1786 m on this one, and it moves with im.
    //
    // DRAG-dt FIX (2026-08-18, ported from ATHAD). Removed the extra *dt, exactly as the
    // Held-Suarez block ~30 lines above did on 2026-07-04 and for the identical reason —
    // the two had been inconsistent with each other in this file ever since. The
    // advective-time non-dimensionalisation of a linear (Rayleigh) damping rate is
    // k_f*L_atm/u_0, already dimensionless; RK4 then integrates it with its OWN *dt
    // (RungeKutta_Atm_Turb.cpp:169,201,228,255) like every other term in rhs_v/rhs_w. The
    // extra *dt made the drag enter as dt^2.
    //
    // IT IS WORSE HERE THAN IN ATHAD, because the error is quadratic in dt and this fork
    // runs dt_visc = 4e-5 against ATHAD's 1e-4. With L_atm = 6287.5, u_0 = 8:
    //     surf_drag                  = 3.6386e-07   (already carried *dt)
    //     RK4's own *dt -> per step   = 1.4554e-11
    //     over a 200-iteration run    = 2.9e-09     <- BELOW the ~1e-8 reproducibility floor
    // ATHAD's equivalent is 4.5e-08, so the shipped drag here was a further 15.6x weaker.
    // In ATHAD a four-arm scan spanning 100x in rayleigh_kf returned Psi_max agreeing to 8
    // significant figures at every checkpoint: that scan measured the noise floor, not the
    // drag. Any such scan here would have been even more purely noise.
    //
    // A SIDE EFFECT WORTH NAMING: under the dt^2 form the effective drag scaled as dt^2, so
    // ATHAD and ATHAD_COND differed in actual drag strength by (1e-4/4e-5)^2 = 6.25x purely
    // from a timestep neither chose as a physics setting. After this fix the coefficient is
    // independent of dt, as a physical rate should be, and the two forks differ only through
    // L_atm as intended.
    //
    // This is a correctness fix, not a regime change: the fixed baseline damps v by
    // 3.6e-07 per step, ~7.3e-05 over 200 iterations. A 1/day drag cannot do more, because
    // 200 iterations is roughly 6 s (on L_atm) to 2 min (on the 120 km shell) of physical
    // time — see ATHAD README item 47 on the ambiguous time unit. It is also strictly safer
    // than the same repair on the buoyancy term that ATHAD item 34 pairs it with: drag is a
    // DAMPING term, so a larger coefficient is stabilising, and item 34's "336x drove a
    // polar vertical runaway" caveat is about buoyancy, a body force.
    //
    // THIS CHANGES RESULTS. Every Psi, KE and wind number recorded here predates it.

    double drag_profile = 1.0 - (double)(i - i_topography[j][k]) / drag_n_layers;   // both from config
    if(drag_profile < 0.0) drag_profile = 0.0;
    if(drag_profile > 1.0) drag_profile = 1.0;
    double surf_drag = (rayleigh_kf * L_atm / u_0) * drag_profile;

    rhs_v.x[i][j][k] = -dpdthe_invrm - transport_v + diffusion_v
        + coriolis * force_nd * coriolis_the
        + coeff_MC_vel * MC_v.x[i][j][k]
        - surf_drag * v_ijk;

    // ---- zonal-mean v momentum-budget term capture (checkpoint iters only) ----
    // Store each rhs_v contribution (nondim dv*/dt*) so write_v_momentum_budget can
    // attribute the meridional-cell spin-down. transport_v splits into vertical
    // (-u·∂v/∂r) and horizontal advection. The four RK4 stages overwrite the same
    // cell; the last (stage 4) wins, which differs from stage 1 only by O(dt)~0.05%.
    if(vbudget_capture){
        vbud_pgf.x[i][j][k]   = -dpdthe_invrm;
        vbud_cor.x[i][j][k]   =  coriolis * force_nd * coriolis_the;
        vbud_advv.x[i][j][k]  = -(u_exp * dvdr_adv);
        vbud_advh.x[i][j][k]  = -(v_invrm * dvdthe_adv + w_invrs * dvdphi_adv);
        vbud_diff.x[i][j][k]  =  diffusion_v;
        vbud_other.x[i][j][k] =  coeff_MC_vel * MC_v.x[i][j][k] - surf_drag * v_ijk;
    }

    rhs_w.x[i][j][k] = -dpdphi_invrs - transport_w + diffusion_w
        + coriolis * force_nd * coriolis_phi
        + coeff_MC_vel * MC_w.x[i][j][k]
        - surf_drag * w_ijk;

    // ---- zonal-mean w (zonal-wind / trade) momentum-budget term capture ----
    // Mirror of the vbud_* block above, for the ZONAL wind: attributes the trade-easterly
    // spin-down. wbud_cor (Coriolis) is the trade SOURCE (turns the meridional Hadley inflow
    // into easterlies); a growing wbud_diff / wbud_other (diffusion / drag) or the post-RK4
    // filters (captured separately in run_3D_loop) would be the SINK.
    if(wbudget_capture){
        wbud_pgf.x[i][j][k]   = -dpdphi_invrs;
        wbud_cor.x[i][j][k]   =  coriolis * force_nd * coriolis_phi;
        wbud_advv.x[i][j][k]  = -(u_exp * dwdr_adv);
        wbud_advh.x[i][j][k]  = -(v_invrm * dwdthe_adv + w_invrs * dwdphi_adv);
        wbud_diff.x[i][j][k]  =  diffusion_w;
        wbud_other.x[i][j][k] =  coeff_MC_vel * MC_w.x[i][j][k] - surf_drag * w_ijk;
    }

    rhs_c.x[i][j][k] = -transport_c + diffusion_c
        + coeff_trans * S_v.x[i][j][k] * r_humid.x[i][j][k]
        + coeff_MC_q * MC_q.x[i][j][k];

    rhs_cloud.x[i][j][k] = -transport_cloud + diffusion_cloud
        + coeff_trans * S_c.x[i][j][k] * r_humid.x[i][j][k];

    rhs_ice.x[i][j][k] = -transport_ice + diffusion_ice
        + coeff_trans * S_i.x[i][j][k] * r_humid.x[i][j][k];

    rhs_g.x[i][j][k] = -transport_g + diffusion_g
        + coeff_trans * S_g.x[i][j][k] * r_humid.x[i][j][k];

    rhs_co2.x[i][j][k] = -transport_co2 + diffusion_co2;

    rhs_tke.x[i][j][k] = -transport_tke + diffusion_tke + tke_source.x[i][j][k];
    rhs_dis.x[i][j][k] = -transport_dis + diffusion_dis + dis_source.x[i][j][k];

    aux_u.x[i][j][k] = rhs_u.x[i][j][k] + dpdr_exp;
    aux_v.x[i][j][k] = rhs_v.x[i][j][k] + dpdthe_invrm;
    aux_w.x[i][j][k] = rhs_w.x[i][j][k] + dpdphi_invrs;
}

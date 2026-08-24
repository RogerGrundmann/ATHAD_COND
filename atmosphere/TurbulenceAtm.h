/*
 * Atmosphere General Circulation Model (AGCM) applied to turbulent flow
 * Standalone turbulence class: k-epsilon (Chien 1982),
 *                              k-omega   (Wilcox 1988),
 *                              k-omega SST (Menter 1994)
 * Turbulence arrays (tke, dis, nue, prod, tke_source, dis_source) are owned
 * by cAtmosphereModel and accessed here through the model reference m.
*/

#pragma once

#include "cAtmosphereModel.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>

using namespace AtomUtils;

class TurbulenceAtm {
public:
    enum Model { k_epsilon, k_omega, k_omega_SST };

    explicit TurbulenceAtm(cAtmosphereModel& model,
                            double           vel_star_ref = 0.4,        // reference friction velocity for re_turb [m/s]
                            double           z_0          = 0.1,        // roughness length    [m]
                            double           nue_air      = 1.5e-5)     // kin. viscosity air  [m²/s]
        : m(model),
          turb_model(parse_model(model.turb_model)),
          vel_star_ref(vel_star_ref),
          z_0(z_0),
          nue_air(nue_air),
          re_turb(vel_star_ref * z_0 / nue_air)
    {}

    // -----------------------------------------------------------------------
    void init() {
        using namespace std;
        cout << endl << "      AGCM: TurbulenceAtm::init" << endl;
        print_model_name();

        auto begin = chrono::high_resolution_clock::now();

        m.re_turb = re_turb;
        compute_vel_star();
        init_fields();
        apply_wall_bc();
        compute_sources();   // prime nue, tke_source, dis_source from the ABL profile
        zero_land_cells();
        clamp_nue();

        auto end     = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for TurbulenceAtm::init\n",
               elapsed.count() * 1e-9);
        cout << "      AGCM: TurbulenceAtm::init ended" << endl;
    }

    // -----------------------------------------------------------------------
    void run() {
        using namespace std;
        cout << endl << "      AGCM: TurbulenceAtm::run" << endl;
        print_model_name();

        auto begin = chrono::high_resolution_clock::now();

        compute_vel_star();
        compute_sources();
        zero_land_cells();
        clamp_nue();

        auto end     = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for TurbulenceAtm::run\n",
               elapsed.count() * 1e-9);
        cout << "      AGCM: TurbulenceAtm::run ended" << endl;
    }

    // -----------------------------------------------------------------------
    // Compute the friction velocity u_τ per (j,k) column from the horizontal
    // wind speed at the first air cell above the topographic surface.
    //   U_horiz = u_0 · sqrt(v² + w²)  at i = i_mount + 1
    //   u_τ     = max(vel_star_min, sqrt(C_D) · U_horiz)
    // C_D = 0.002 (neutral drag coefficient); vel_star_min = 0.05 m/s so
    // that calm-wind columns still receive a small but non-zero ABL seed.
    void compute_vel_star() {
        constexpr double C_D          = 0.002;
        constexpr double vel_star_min = 0.05;   // [m/s]

        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                const int i_surf = std::min(m.i_topography[j][k] + 1, m.im - 1);
                const double v_s = m.v.x[i_surf][j][k];
                const double w_s = m.w.x[i_surf][j][k];
                const double U_horiz = m.u_0 * std::sqrt(v_s * v_s + w_s * w_s);
                m.vel_star.y[j][k] = std::max(vel_star_min, std::sqrt(C_D) * U_horiz);
            }
        }
    }

    // -----------------------------------------------------------------------
    // Wall boundary condition for ω at i = 0 (ground surface).
    // Physical formulas produce dimensional ω [s⁻¹]; normalise to ω* = ω · L_atm/u_0.
    // Public so run_3D_loop can reassert the wall values after bcRadius
    // overwrites tke[0] / dis[0] with cubic Neumann extrapolation.
    void apply_wall_bc() {
        const double y1      = m.get_layer_height(1);                   // first-cell height above surface [m]
        const double y1_safe = std::max(y1, 1.0e-6);
        const double nd_omega = m.L_atm / m.u_0;                        // ω_phys [s⁻¹] → ω*
        const double inv_u02  = 1.0 / (m.u_0 * m.u_0);

        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                const double vs = m.vel_star.y[j][k];

                // tke wall BC — log-layer equilibrium k = u_τ²/√Cμ (all models).
                // Reasserts the ghost-cell value after bcRadius cubic extrapolation
                // so the diffusion flux at i=1 always has the correct lower face value,
                // preventing near-surface tke from being drained as the profile flattens.
                m.tke.x[0][j][k] = vs * vs / (std::sqrt(C_nue) * m.u_0 * m.u_0);

                // nue wall BC — turbulent viscosity vanishes at the no-slip wall
                // (viscous sublayer: ν_T → 0).  bcRadius cubic extrapolation can
                // push nue[0] negative; clamp it to 0 here unconditionally.
                m.nue.x[0][j][k] = 0.0;

                if (turb_model == k_epsilon) {
                    // Chien (1982) low-Re k-ε: ε_wall = 0 (Dirichlet).
                    // bcRadius cubic extrapolation can overshoot to negative ε;
                    // the no-slip condition requires ε→0 as y→0.
                    m.dis.x[0][j][k] = 0.0;
                } else if (turb_model == k_omega) {
                    // Wilcox 1988: ω_wall = 6ν / (β₁ y₁²),  β₁ = 0.075
                    const double om_phys = 6.0 * nue_air / (0.075 * y1_safe * y1_safe);
                    m.dis.x[0][j][k] = om_phys * nd_omega;
                } else if (turb_model == k_omega_SST) {
                    // Menter 1994: ω_wall = 60ν / (β₁ y₁²),  β₁ = 0.0333
                    const double om_phys = 60.0 * nue_air / (0.0333 * y1_safe * y1_safe);
                    m.dis.x[0][j][k] = om_phys * nd_omega;
                }
            }
        }
    }

private:
    cAtmosphereModel& m;
    Model  turb_model;
    double vel_star_ref;   // reference value for re_turb; per-column u_τ lives in m.vel_star
    double z_0;
    double nue_air;
    double re_turb;    // = vel_star · z_0 / nue_air  (computed, not a free parameter)

    // ABL model constants
    static constexpr double C_nue   = 0.028;
    static constexpr double Karman  = 0.42;
    // Safety floors / caps (dimensionless)
    static constexpr double dis_min = 1.0e-10;                          // minimum dis* to avoid nue → ∞

    // Map the string parameter m.turb_model to the internal enum
    static Model parse_model(const std::string& s) {
        if (s == "k_epsilon")   return k_epsilon;
        if (s == "k_omega")     return k_omega;
        if (s == "k_omega_SST") return k_omega_SST;
        return k_omega_SST;                                             // default / "none" fallback
    }

    // SST blending: inner (zone 1, near wall) and outer (zone 2, free stream)
    static double blend(double inner, double outer, double F1) {
        return F1 * inner + (1.0 - F1) * outer;
    }

    // -----------------------------------------------------------------------
    void print_model_name() const {
        using namespace std;
        switch(turb_model) {
            case k_epsilon:   cout << "      k-epsilon turbulence model (Chien 1982)" << endl; break;
            case k_omega:     cout << "      k-omega turbulence model (Wilcox 1988)" << endl;  break;
            case k_omega_SST: cout << "      k-omega SST turbulence model (Menter 1994)" << endl; break;
        }
    }

    // -----------------------------------------------------------------------
    // Initialise tke and dis throughout the domain (ABL profile, Narjisse).
    // All source quantities (vel_star, get_layer_height) are dimensional;
    // results are normalised to the model's dimensionless convention:
    //   k*  = k_phys / u_0²
    //   ω*  = ω_phys · L_atm / u_0        (k-ω, k-ω SST)
    //   ε*  = ε_phys · L_atm / u_0³       (k-ε)
    void init_fields() {
        const double inv_u02  = 1.0 / (m.u_0 * m.u_0);
        const double nd_omega = m.L_atm / m.u_0;                        // ω_phys → ω*
        const double nd_eps   = m.L_atm / (m.u_0 * m.u_0 * m.u_0);      // ε_phys → ε*

        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 1; j < m.jm - 1; j++) {
            for (int k = 1; k < m.km - 1; k++) {
                // Per-column friction velocity and derived background floors.
                // k_bg: 1% of peak surface TKE — prevents k→0 and ω/ε→0 at altitude.
                const double vs    = m.vel_star.y[j][k];
                const double k_bg  = 0.01 * vs * vs / sqrt(C_nue);      // [m²/s²]
                const double om_bg = pow(C_nue, -0.25) * sqrt(k_bg) / (Karman * m.L_atm); // [s⁻¹]
                const double eps_bg = pow(C_nue, 0.75) * pow(k_bg, 1.5) / (Karman * m.L_atm); // [m²/s³]

                // ABL is anchored to the LOCAL surface, not sea level: heights are
                // measured above the terrain top, so the boundary layer rides on top
                // of mountains (ABL top = mountain_height + abl_height) instead of
                // being clipped away wherever the ground itself exceeds abl_height.
                // Consistent with the terrain-relative wall distance used for ω/ε.
                const double y_mount = m.get_layer_height(m.i_topography[j][k]); // surface height [m]

                for (int i = 1; i < m.im - 1; i++) {
                    const double z_i    = m.get_layer_height(i) - y_mount; // height above local surface [m]
                    const double z_safe = std::max(z_i, 1.0e-6);

                    // Clamp z_frac to [0,1] so the parabolic profile never
                    // exceeds 1 above the physical ABL depth abl_height (prevents
                    // huge TKE aloft). abl_height (~1.5 km, config-tunable) is the
                    // meteorological boundary-layer top and is DELIBERATELY decoupled
                    // from the grid length scale L_atm (=400 m): pinning it to L_atm
                    // collapsed the surface-driven TKE profile to its background floor
                    // by i≈8 (~400 m), so turbulence was seeded only in the lowest
                    // ~400 m instead of the full boundary layer.
                    const double z_frac = std::min(z_i, m.abl_height) / m.abl_height;

                    // Physical TKE [m²/s²]: ABL parabolic profile with background floor
                    const double k_abl  = vs * vs / sqrt(C_nue)
                        * pow(1.0 - z_frac, 2);
                    const double k_phys = std::max(k_abl, k_bg);

                    double dis_nd = 0.0;
                    if (turb_model == k_epsilon) {
                        // Physical ε [m²/s³]: ε = C_μ^0.75 · k^1.5 / (κ · z)
                        const double eps_phys = pow(C_nue, 0.75)
                            * pow(k_phys, 1.5) / (Karman * z_safe);
                        dis_nd = std::max(eps_phys, eps_bg) * nd_eps;
                    } else {
                        // Physical ω [s⁻¹]: ω = C_μ^-0.25 · k^0.5 / (κ · z)
                        // Floor at om_bg to prevent ω→0 aloft driving nue=k/ω → ∞.
                        const double om_phys = pow(C_nue, -0.25)
                            * sqrt(k_phys) / (Karman * z_safe);
                        dis_nd = std::max(om_phys, om_bg) * nd_omega;
                    }

                    m.tke.x[i][j][k]  = std::max(0.0, k_phys * inv_u02);
                    m.dis.x[i][j][k]  = std::max(dis_min, dis_nd);
                    m.tken.x[i][j][k] = m.tke.x[i][j][k];
                    m.disn.x[i][j][k] = m.dis.x[i][j][k];
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Global nue clamp: enforce 0 ≤ nue* ≤ nue_max everywhere, including pole
    // rows (j=0, j=jm-1) that are not touched by compute_sources and may receive
    // negative values from cubic BC extrapolation in bcTheta.
    //
    // Eddy viscosity is ALSO confined to the boundary layer: above the ABL ceiling
    // (abl_height) the free troposphere is treated as non-turbulent, so nue* is
    // tapered to zero there and only the molecular 1/re_turb diffusion remains in
    // the RHS. The background TKE/ω floors (k_bg, om_bg) keep nue=k/ω finite aloft,
    // but that residual eddy viscosity homogenises the meridional momentum / vertical
    // shear the Hadley–Ferrel cells live on — the same spin-down failure mode the
    // radial-Shapiro fix addressed ([[project_radial_shapiro_spindown]]). The taper
    // runs over ~one local grid layer above abl_height to avoid an abrupt
    // diffusion-coefficient step at the interface.
    void clamp_nue() {
        const double nue_max = 1000.0 / (m.u_0 * m.L_atm);

        // Taper depth = thickness of the (stretched) grid layer straddling the ABL top.
        int i_abl = 0;
        while (i_abl < m.im - 1 && m.get_layer_height(i_abl) <= m.abl_height) i_abl++;
        const double taper_dz = std::max(1.0,
            (double)m.get_layer_height(std::min(i_abl, m.im - 1))
          - (double)m.get_layer_height(std::max(i_abl - 1, 0)));

        #pragma omp parallel for collapse(3) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    // Height above the LOCAL surface: the ABL ceiling rides on top of
                    // terrain (cutoff at mountain_height + abl_height), matching the
                    // terrain-relative seeding in init_fields.
                    const double z_agl = (double)m.get_layer_height(i)
                                       - (double)m.get_layer_height(m.i_topography[j][k]);
                    // f = 1 at/below abl_height (above ground), ramps to 0 over taper_dz above.
                    const double f = std::clamp(
                        (m.abl_height + taper_dz - z_agl) / taper_dz, 0.0, 1.0);
                    m.nue.x[i][j][k] = std::clamp(m.nue.x[i][j][k], 0.0, nue_max) * f;
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Zero all turbulence quantities inside land cells
    void zero_land_cells() {
        #pragma omp parallel for collapse(3) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    if (is_land(m.h, i, j, k)) {
                        m.tke.x[i][j][k]        = 0.0;
                        m.tken.x[i][j][k]       = 0.0;
                        m.dis.x[i][j][k]        = 0.0;
                        m.disn.x[i][j][k]       = 0.0;
                        m.nue.x[i][j][k]        = 0.0;
                        m.prod.x[i][j][k]       = 0.0;
                        m.tke_source.x[i][j][k] = 0.0;
                        m.dis_source.x[i][j][k] = 0.0;
                    }
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Main turbulence source computation
    void compute_sources() {
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 1; j < m.jm - 1; j++) {
            for (int k = 1; k < m.km - 1; k++) {
                const int    i_mount    = m.i_topography[j][k];
                const double y_mount    = m.get_layer_height(i_mount);  // surface physical height [m]

                // Background tke floor for this column [dimensionless k* = k/u_0²].
                // Prevents tke from decaying below the ABL seed level in low-shear
                // regions where P_k is insufficient to offset destruction.
                const double vs      = m.vel_star.y[j][k];
                const double k_bg_nd = 0.01 * vs * vs / (std::sqrt(C_nue) * m.u_0 * m.u_0);

                for (int i = 1; i < m.im - 1; i++) {
                    if (is_land(m.h, i, j, k)) continue;

                    // ---- geometry ----
                    const double rm           = m.rad.z[i];
                    const double exp_rm       = m.metricExpRm(rm);
                    double sinthe             = sin(m.the.z[j]);
                    if (sinthe == 0.0) sinthe = 1.0e-5;
                    const double rmsinthe     = rm * sinthe;
                    const double inv_2dr      = 1.0 / (2.0 * m.dr);
                    const double inv_2dthe    = 1.0 / (2.0 * m.dthe);
                    const double inv_2dphi    = 1.0 / (2.0 * m.dphi);
                    // ATM_METRIC_RADIUS — the 1/r in a horizontal derivative is the planetary
                    // radius, not the grid coordinate. Identity when the knob is off.
                    const double rmet             = m.metricRadius(rm);
                    const double inv_rm2dthe      = inv_2dthe / rmet;
                    const double inv_rmsinthe2dphi = inv_2dphi / (rmet * sinthe);

                    // Enforce floors on this cell before any reads — the RK4 can push
                    // tke or dis below their safe limits between compute_sources calls.
                    m.tke.x[i][j][k] = std::max(0.0,     m.tke.x[i][j][k]);
                    m.dis.x[i][j][k] = std::max(dis_min, m.dis.x[i][j][k]);

                    // ---- velocity gradients ----
                    const double dudr   = (m.u.x[i+1][j][k] - m.u.x[i-1][j][k]) * inv_2dr * exp_rm;
                    const double dvdr   = (m.v.x[i+1][j][k] - m.v.x[i-1][j][k]) * inv_2dr * exp_rm;
                    const double dwdr   = (m.w.x[i+1][j][k] - m.w.x[i-1][j][k]) * inv_2dr * exp_rm;

                    const double dudthe = (m.u.x[i][j+1][k] - m.u.x[i][j-1][k]) * inv_rm2dthe;
                    const double dvdthe = (m.v.x[i][j+1][k] - m.v.x[i][j-1][k]) * inv_rm2dthe;
                    const double dwdthe = (m.w.x[i][j+1][k] - m.w.x[i][j-1][k]) * inv_rm2dthe;

                    const double dudphi = (m.u.x[i][j][k+1] - m.u.x[i][j][k-1]) * inv_rmsinthe2dphi;
                    const double dvdphi = (m.v.x[i][j][k+1] - m.v.x[i][j][k-1]) * inv_rmsinthe2dphi;
                    const double dwdphi = (m.w.x[i][j][k+1] - m.w.x[i][j][k-1]) * inv_rmsinthe2dphi;

                    // ---- tke / dis gradients with Neumann BC at land faces ----
                    // Replace each land neighbour's value with the current cell's value
                    // (zero-gradient) so stencils do not read the zeroed land-cell values.
                    const bool neumann_bot = (i > 0      && is_land(m.h, i-1, j, k));
                    const bool neumann_jm1 = (j > 0      && is_land(m.h, i, j-1, k));
                    const bool neumann_jp1 = (j < m.jm-1 && is_land(m.h, i, j+1, k));
                    const bool neumann_km1 = (k > 0      && is_land(m.h, i, j, k-1));
                    const bool neumann_kp1 = (k < m.km-1 && is_land(m.h, i, j, k+1));

                    const double tke_im1 = neumann_bot ? m.tke.x[i][j][k] : m.tke.x[i-1][j][k];
                    const double dis_im1 = neumann_bot ? m.dis.x[i][j][k] : m.dis.x[i-1][j][k];
                    const double tke_jm1 = neumann_jm1 ? m.tke.x[i][j][k] : m.tke.x[i][j-1][k];
                    const double tke_jp1 = neumann_jp1 ? m.tke.x[i][j][k] : m.tke.x[i][j+1][k];
                    const double dis_jm1 = neumann_jm1 ? m.dis.x[i][j][k] : m.dis.x[i][j-1][k];
                    const double dis_jp1 = neumann_jp1 ? m.dis.x[i][j][k] : m.dis.x[i][j+1][k];
                    const double tke_km1 = neumann_km1 ? m.tke.x[i][j][k] : m.tke.x[i][j][k-1];
                    const double tke_kp1 = neumann_kp1 ? m.tke.x[i][j][k] : m.tke.x[i][j][k+1];
                    const double dis_km1 = neumann_km1 ? m.dis.x[i][j][k] : m.dis.x[i][j][k-1];
                    const double dis_kp1 = neumann_kp1 ? m.dis.x[i][j][k] : m.dis.x[i][j][k+1];

                    const double dtkedr   = (m.tke.x[i+1][j][k] - tke_im1) * inv_2dr * exp_rm;
                    const double ddisdr   = (m.dis.x[i+1][j][k] - dis_im1) * inv_2dr * exp_rm;
                    const double dtkedthe = (tke_jp1 - tke_jm1) * inv_rm2dthe;
                    const double ddisdthe = (dis_jp1 - dis_jm1) * inv_rm2dthe;
                    const double dtkedphi = (tke_kp1 - tke_km1) * inv_rmsinthe2dphi;
                    const double ddisdphi = (dis_kp1 - dis_km1) * inv_rmsinthe2dphi;

                    // ---- production tensor P_k ----
                    // Compute νt fresh from current tke/dis rather than reading m.nue,
                    // which lags one call behind (nue is updated inside compute_k_*
                    // AFTER prod is needed here).  On the first call from init() m.nue
                    // is still 0, giving prod≈0 and breaking the turbulence bootstrap.
                    // k-ε: νt = Cμ·k²/ε   k-ω / k-ω SST: νt ≈ k/ω  (F2 limiter is
                    // a minor correction; the full SST value is recomputed below).
                    double cnue;
                    {
                        // Use dis_min (not 1e-20) so that RK4-driven dips below dis_min
                        // do not produce cnue=k/1e-20 → inf in prod.
                        const double dis_here = std::max(m.dis.x[i][j][k], dis_min);
                        const double tke_here = std::max(m.tke.x[i][j][k], 0.0);
                        const double nue_max  = 1000.0 / (m.u_0 * m.L_atm);
                        cnue = (turb_model == k_epsilon)
                             ? C_nue * tke_here * tke_here / dis_here
                             : tke_here / dis_here;
                        cnue = std::min(cnue, nue_max);                 // guard against large tke spikes
                    }
                    const double der  = 0.66667 * (dudr + dvdthe + dwdphi);

                    m.prod.x[i][j][k] = std::max(0.0,
                          (cnue * (2.0 * dudr   - der) - 0.66667 * m.tke.x[i][j][k]) * dudr
                        + (cnue * (dudthe + dvdr))                                   * dudthe
                        + (cnue * (dudphi + dwdr))                                   * dudphi
                        + (cnue * (2.0 * dvdthe - der) - 0.66667 * m.tke.x[i][j][k]) * dvdthe
                        + (cnue * (dvdr + dudthe))                                   * dvdr
                        + (cnue * (dvdphi + dwdthe))                                 * dvdphi
                        + (cnue * (2.0 * dwdphi - der) - 0.66667 * m.tke.x[i][j][k]) * dwdphi
                        + (cnue * (dwdr + dudphi))                                   * dwdr
                        + (cnue * (dwdthe + dvdphi))                                 * dwdthe);

                    // ---- vorticity magnitude Ω = |∇×u| ----
                    // Each cross-derivative pair must be squared and summed — not
                    // added with its negated partner (which would cancel to zero).
                    const double W12 = dudthe - dvdr;
                    const double W13 = dudphi - dwdr;
                    const double W23 = dvdphi - dwdthe;
                    const double Omega = sqrt(W12*W12 + W13*W13 + W23*W23);

                    // ================================================================
                    if (turb_model == k_epsilon) {
                        compute_k_epsilon(i, j, k, y_mount,
                            dudr, dvdr, dwdr, dudthe, dvdthe, dwdthe, dudphi, dvdphi, dwdphi,
                            Omega);
                    }
                    // ================================================================
                    else if (turb_model == k_omega) {
                        compute_k_omega(i, j, k,
                            dtkedr, ddisdr, dtkedthe, ddisdthe, dtkedphi, ddisdphi,
                            Omega, Omega);
                    }
                    // ================================================================
                    else if (turb_model == k_omega_SST) {
                        compute_k_omega_SST(i, j, k, y_mount, rm, sinthe,
                            dtkedr, ddisdr, dtkedthe, ddisdthe, dtkedphi, ddisdphi,
                            Omega, Omega);
                    }

                    // nue_max: 1000 m²/s physical cap converted to dimensionless ν* = ν/(u_0·L_atm)
                    const double nue_max = 1000.0 / (m.u_0 * m.L_atm);
                    m.tke.x[i][j][k] = std::max(k_bg_nd, m.tke.x[i][j][k]);
                    m.dis.x[i][j][k] = std::max(dis_min, m.dis.x[i][j][k]);
                    m.nue.x[i][j][k] = std::clamp(m.nue.x[i][j][k], 0.0, nue_max);
                    // dis_source_max: 20·ω*² cap — prevents polar blow-up from the 1/sinθ²
                    // amplification in the cross-diffusion term D_w near j→0/jm-1.
                    // Factor 20 gives ~20×β_max ≈ 1.8× max destruction, ample headroom.
                    const double dis_src_max = 20.0 * m.dis.x[i][j][k] * m.dis.x[i][j][k];
                    m.dis_source.x[i][j][k] = std::clamp(m.dis_source.x[i][j][k],
                                                          -dis_src_max, dis_src_max);
                    // tke_source_max: 20·k*·ω* cap — mirrors dis_source logic;
                    // bounds the P_k-Y_k balance without suppressing legitimate physics.
                    const double tke_src_max = 20.0 * m.tke.x[i][j][k] * m.dis.x[i][j][k];
                    m.tke_source.x[i][j][k] = std::clamp(m.tke_source.x[i][j][k],
                                                          -tke_src_max, tke_src_max);
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // y_mount: physical height of the topographic surface [m] (get_layer_height(i_mount))
    void compute_k_epsilon(int i, int j, int k, double y_mount,
        double dudr, double dvdr, double dwdr,
        double dudthe, double dvdthe, double dwdthe,
        double dudphi, double dvdphi, double dwdphi,
        double /*Omega*/)
    {
        const double C_eps_1 = 1.35;
        const double C_eps_2 = 1.80;

        // Physical wall distance [m] and its dimensionless form y* = y/L_atm.
        // Floor at L_atm (one grid cell) so that y_star ≥ 1 and the 1/y_star² in
        // the Chien D term stays bounded for cells straddling mountain tops.
        const double y_phys = std::max(m.get_layer_height(i) - y_mount, m.L_atm);  // [m]
        const double y_star = y_phys / m.L_atm;  // dimensionless, always ≥ 1

        // ν_T* = C_μ · k*² / ε*   (already dimensionless — see header)
        m.nue.x[i][j][k] = C_nue * m.tke.x[i][j][k] * m.tke.x[i][j][k]
                          / std::max(m.dis.x[i][j][k], dis_min);

        // y⁺ = y_phys · u_τ / ν_air  (viscous wall units, dimensionless)
        const double d_plus = y_phys * m.vel_star.y[j][k] / nue_air;
        const double f_nue  = 1.0 - exp(-0.0115 * d_plus);

        // Turbulent Reynolds number Re_T = k*² · u_0 · L_atm / (ε* · ν_air)
        // (uses molecular viscosity, not eddy viscosity)
        const double Re_T = m.tke.x[i][j][k] * m.tke.x[i][j][k]
                            * m.u_0 * m.L_atm
                            / std::max(m.dis.x[i][j][k] * nue_air, dis_min * nue_air);
        const double f_2  = 1.0 - 0.4 / 1.8 * exp(-Re_T * Re_T / 36.0);

        m.nue.x[i][j][k] *= f_nue;

        // Chien (1982) viscous near-wall damping terms D and E.
        // y_star ≥ 1 (enforced above), so y_star² ≥ 1 and L_k is bounded.
        // For ABL cells (d_plus ≈ 10⁷), exp(-d_plus/2) ≈ 0 so L_w ≈ 0.
        const double y_star2 = y_star * y_star;
        const double L_k = -2.0 * m.tke.x[i][j][k] / y_star2;
        const double L_w = -2.0 * m.dis.x[i][j][k] / y_star2 * exp(-d_plus / 2.0);

        const double P_k = m.prod.x[i][j][k];
        const double Y_k = m.dis.x[i][j][k];
        // Standard k-ε: P_ε = C_ε1 · (ε/k) · P_k  — NOT k·ε·P_k (wrong by k²).
        const double tke_safe = std::max(m.tke.x[i][j][k], dis_min);
        const double P_w = C_eps_1 * m.dis.x[i][j][k] / tke_safe * P_k;
        const double Y_w = C_eps_2 * f_2 * m.dis.x[i][j][k] * m.dis.x[i][j][k]
                           / tke_safe;

        m.tke_source.x[i][j][k] = P_k - Y_k + L_k / re_turb;
        m.dis_source.x[i][j][k] = P_w - Y_w + L_w / re_turb;
    }

    // -----------------------------------------------------------------------
    void compute_k_omega(int i, int j, int k,
        double dtkedr, double ddisdr,
        double dtkedthe, double ddisdthe,
        double dtkedphi, double ddisdphi,
        double Omega_mag, double /*W_unused*/)
    {
        const double bet_star = 0.09;
        const double gam      = 0.52;
        const double C_lim    = 0.875;
        const double bet_0    = 0.0708;
        const double sig_k    = 1.0;
        const double sig_w    = 1.3;

        const double rm        = m.rad.z[i];
        double sinthe          = sin(m.the.z[j]);
        if (sinthe == 0.0) sinthe = 1.0e-5;
        const double rmsinthe  = rm * sinthe;
        const double rm2       = rm * rm;
        const double sinthe2   = sinthe * sinthe;
        const double rm2sinthe = rm2 * sinthe;

        // Geometry for exp-stretching — placed here so Neumann corrections below can use exp_rm.
        // exp_rm = 1/(rm+1) is the Jacobian of the radial coordinate transformation;
        // exp_2_rm scales the raw grid second derivative to physical space.
        // dtkedr is already physical (multiplied by exp_rm in compute_sources),
        // so its first-derivative coefficient is simply 2/rm (not 2*exp_rm/rm).
        // dtkedthe already carries the 1/rm factor, so the cosθ coefficient is
        // cosθ/(rm sinθ) rather than cosθ/(rm² sinθ).
        const double exp_rm   = m.metricExpRm(rm);
        const double exp_2_rm = exp_rm * exp_rm;
        // ATM_METRIC_RADIUS — identity when off. See RungeKutta_Atm_Turb for the rationale.
        const double inv_rm   = 1.0 / m.metricRadius(rm);
        const double inv_rm2  = inv_rm * inv_rm;

        // Neumann BC at all six faces: replace a land neighbour's tke/dis with the
        // current cell's value (zero-gradient) so that sig_d, D_w, and the diffusion
        // Laplacian are not corrupted by the zeroed land-cell values.
        const bool neumann_bot = (i > 0      && is_land(m.h, i-1, j, k));
        const bool neumann_jm1 = (j > 0      && is_land(m.h, i, j-1, k));
        const bool neumann_jp1 = (j < m.jm-1 && is_land(m.h, i, j+1, k));
        const bool neumann_km1 = (k > 0      && is_land(m.h, i, j, k-1));
        const bool neumann_kp1 = (k < m.km-1 && is_land(m.h, i, j, k+1));

        const double tke_im1 = neumann_bot ? m.tke.x[i][j][k] : m.tke.x[i-1][j][k];
        const double dis_im1 = neumann_bot ? m.dis.x[i][j][k] : m.dis.x[i-1][j][k];
        const double tke_jm1 = neumann_jm1 ? m.tke.x[i][j][k] : m.tke.x[i][j-1][k];
        const double tke_jp1 = neumann_jp1 ? m.tke.x[i][j][k] : m.tke.x[i][j+1][k];
        const double dis_jm1 = neumann_jm1 ? m.dis.x[i][j][k] : m.dis.x[i][j-1][k];
        const double dis_jp1 = neumann_jp1 ? m.dis.x[i][j][k] : m.dis.x[i][j+1][k];
        const double tke_km1 = neumann_km1 ? m.tke.x[i][j][k] : m.tke.x[i][j][k-1];
        const double tke_kp1 = neumann_kp1 ? m.tke.x[i][j][k] : m.tke.x[i][j][k+1];
        const double dis_km1 = neumann_km1 ? m.dis.x[i][j][k] : m.dis.x[i][j][k-1];
        const double dis_kp1 = neumann_kp1 ? m.dis.x[i][j][k] : m.dis.x[i][j][k+1];

        const double dtkedr_neu  = neumann_bot
            ? (m.tke.x[i+1][j][k] - m.tke.x[i][j][k]) * (1.0 / (2.0 * m.dr)) * exp_rm
            : dtkedr;
        const double ddisdr_neu  = neumann_bot
            ? (m.dis.x[i+1][j][k] - m.dis.x[i][j][k]) * (1.0 / (2.0 * m.dr)) * exp_rm
            : ddisdr;
        const double dtkedthe_neu = (tke_jp1 - tke_jm1) / (rm  * 2.0 * m.dthe);
        const double ddisdthe_neu = (dis_jp1 - dis_jm1) / (rm  * 2.0 * m.dthe);
        const double dtkedphi_neu = (tke_kp1 - tke_km1) / (rmsinthe * 2.0 * m.dphi);
        const double ddisdphi_neu = (dis_kp1 - dis_km1) / (rmsinthe * 2.0 * m.dphi);

        const double dudr   = (m.u.x[i+1][j][k] - m.u.x[i-1][j][k]) / (2.0 * m.dr) / (rm + 1.0);
        const double dvdthe = (m.v.x[i][j+1][k] - m.v.x[i][j-1][k]) / (rm * 2.0 * m.dthe);
        const double dwdphi = (m.w.x[i][j][k+1] - m.w.x[i][j][k-1]) / (rmsinthe * 2.0 * m.dphi);

        // Strain-rate magnitude S = sqrt(2 S_ij S_ij).
        // Off-diagonal pairs (S_12, S_13, S_23) each contribute twice (ij + ji),
        // giving factor 4 per pair; diagonal terms (S_11, S_22, S_33) contribute
        // 2 each.  The divergence-correction (−2/3 div) is omitted here because
        // the Wilcox realizability limiter uses the full strain rate, not the
        // deviatoric part.
        const double dudthe_u  = (m.u.x[i][j+1][k] - m.u.x[i][j-1][k]) / (rm       * 2.0 * m.dthe);
        const double dudphi_u  = (m.u.x[i][j][k+1] - m.u.x[i][j][k-1]) / (rmsinthe * 2.0 * m.dphi);
        const double dvdr_u    = (m.v.x[i+1][j][k] - m.v.x[i-1][j][k]) / (2.0 * m.dr) / (rm + 1.0);
        const double dwdr_u    = (m.w.x[i+1][j][k] - m.w.x[i-1][j][k]) / (2.0 * m.dr) / (rm + 1.0);
        const double dvdphi_u  = (m.v.x[i][j][k+1] - m.v.x[i][j][k-1]) / (rmsinthe * 2.0 * m.dphi);
        const double dwdthe_u  = (m.w.x[i][j+1][k] - m.w.x[i][j-1][k]) / (rm       * 2.0 * m.dthe);
        // S_ij = (∂u_i/∂x_j + ∂u_j/∂x_i)/2;  2 S_ij S_ij sums all 9 components
        const double S11 = dudr,   S22 = dvdthe, S33 = dwdphi;
        const double S12 = 0.5 * (dudthe_u + dvdr_u);
        const double S13 = 0.5 * (dudphi_u + dwdr_u);
        const double S23 = 0.5 * (dvdphi_u + dwdthe_u);
        const double S_mag = sqrt(2.0 * (S11*S11 + S22*S22 + S33*S33
                                       + 2.0*(S12*S12 + S13*S13 + S23*S23)));

        const double w_hat = std::max(m.dis.x[i][j][k],
            C_lim * S_mag / sqrt(bet_star));

        double sig_d = dtkedr_neu * ddisdr_neu + dtkedthe_neu * ddisdthe_neu + dtkedphi_neu * ddisdphi_neu;
        sig_d = (sig_d <= 0.0) ? 0.0 : 0.125;

        // Wilcox (2006) vorticity–strain correction: χ_ω = |Ω²·S| / (β*·ω)³
        // Using scalar magnitudes as a standard rank-3 tensor approximation.
        const double chi_w  = fabs(Omega_mag * Omega_mag * S_mag
            / pow(bet_star * std::max(m.dis.x[i][j][k], 1.0e-20), 3));
        const double f_bet  = (1.0 + 85.0 * chi_w) / (1.0 + 100.0 * chi_w);
        const double bet_wc = bet_0 * f_bet;

        m.nue.x[i][j][k] = m.tke.x[i][j][k] / std::max(w_hat, 1.0e-20);

        const double P_k = std::min(m.prod.x[i][j][k],
            20.0 * bet_star * m.tke.x[i][j][k] * m.dis.x[i][j][k]);
        const double Y_k = bet_star * m.tke.x[i][j][k] * m.dis.x[i][j][k];
        const double P_w = gam * m.dis.x[i][j][k]
                           / std::max(m.tke.x[i][j][k], 1.0e-20) * P_k;
        const double Y_w = bet_wc * m.dis.x[i][j][k] * m.dis.x[i][j][k];
        const double D_w = sig_d / std::max(m.dis.x[i][j][k], 1.0e-20)
            * (dtkedr_neu * ddisdr_neu + dtkedthe_neu * ddisdthe_neu + dtkedphi_neu * ddisdphi_neu);

        // tke_source / dis_source: store only the production–destruction balance.
        // Turbulent diffusion of tke is handled by the RHS (diffusion_tke Laplacian)
        // and is recomputed fresh at every RK4 sub-stage, so it must NOT be added
        // here to avoid double-counting.
        m.tke_source.x[i][j][k] = P_k - Y_k;
        m.dis_source.x[i][j][k] = P_w - Y_w + D_w;
    }

    // -----------------------------------------------------------------------
    // y_mount: physical height of the topographic surface [m] (get_layer_height(i_mount))
    void compute_k_omega_SST(int i, int j, int k,
        double y_mount, double rm, double sinthe,
        double dtkedr, double ddisdr,
        double dtkedthe, double ddisdthe,
        double dtkedphi, double ddisdphi,
        double Omega, double /*W*/)
    {
        // Menter 1994 ABL constants
        const double a1       = 0.31;
        const double bet_star = 0.09;    // β* destruction coefficient — NOT C_μ=0.028
        const double bet1     = 0.0333;
        const double bet2     = 0.0368;
        const double gam1     = 0.413;
        const double gam2     = 0.2;
        const double sig_k1   = 1.176;
        const double sig_k2   = 1.0;
        const double sig_w1   = 2.0;
        const double sig_w2   = 1.168;

        // Dimensionless molecular viscosity ν* = ν_air / (u_0 · L_atm)
        // Used in blending-function arguments that require molecular viscosity.
        const double nue_air_nd = nue_air / (m.u_0 * m.L_atm);

        // Physical wall distance normalised to dimensionless y* = y/L_atm.
        // Menter's blending functions require the wall-normal distance, not a
        // stretched coordinate difference; use get_layer_height for consistency.
        const double y_phys = std::max(m.get_layer_height(i) - y_mount, 1.0e-6);
        const double y_star = y_phys / m.L_atm;  // dimensionless

        const double rm2      = rm * rm;
        const double sinthe2  = sinthe * sinthe;
        const double rmsinthe = rm * sinthe;

        // Geometry factors matching compute_k_omega convention
        const double exp_rm   = m.metricExpRm(rm);
        const double exp_2_rm = exp_rm * exp_rm;
        // ATM_METRIC_RADIUS — identity when off. See RungeKutta_Atm_Turb for the rationale.
        const double inv_rm   = 1.0 / m.metricRadius(rm);
        const double inv_rm2  = inv_rm * inv_rm;

        // Neumann BC at all six faces: replace a land neighbour's tke/dis with the
        // current cell's value (zero-gradient) so that CD_kw, D_w, and the diffusion
        // Laplacian are not corrupted by the zeroed land-cell values.
        const bool neumann_bot = (i > 0      && is_land(m.h, i-1, j, k));
        const bool neumann_jm1 = (j > 0      && is_land(m.h, i, j-1, k));
        const bool neumann_jp1 = (j < m.jm-1 && is_land(m.h, i, j+1, k));
        const bool neumann_km1 = (k > 0      && is_land(m.h, i, j, k-1));
        const bool neumann_kp1 = (k < m.km-1 && is_land(m.h, i, j, k+1));

        const double tke_im1 = neumann_bot ? m.tke.x[i][j][k] : m.tke.x[i-1][j][k];
        const double dis_im1 = neumann_bot ? m.dis.x[i][j][k] : m.dis.x[i-1][j][k];
        const double tke_jm1 = neumann_jm1 ? m.tke.x[i][j][k] : m.tke.x[i][j-1][k];
        const double tke_jp1 = neumann_jp1 ? m.tke.x[i][j][k] : m.tke.x[i][j+1][k];
        const double dis_jm1 = neumann_jm1 ? m.dis.x[i][j][k] : m.dis.x[i][j-1][k];
        const double dis_jp1 = neumann_jp1 ? m.dis.x[i][j][k] : m.dis.x[i][j+1][k];
        const double tke_km1 = neumann_km1 ? m.tke.x[i][j][k] : m.tke.x[i][j][k-1];
        const double tke_kp1 = neumann_kp1 ? m.tke.x[i][j][k] : m.tke.x[i][j][k+1];
        const double dis_km1 = neumann_km1 ? m.dis.x[i][j][k] : m.dis.x[i][j][k-1];
        const double dis_kp1 = neumann_kp1 ? m.dis.x[i][j][k] : m.dis.x[i][j][k+1];

        const double dtkedr_neu  = neumann_bot
            ? (m.tke.x[i+1][j][k] - m.tke.x[i][j][k]) * (1.0 / (2.0 * m.dr)) * exp_rm
            : dtkedr;
        const double ddisdr_neu  = neumann_bot
            ? (m.dis.x[i+1][j][k] - m.dis.x[i][j][k]) * (1.0 / (2.0 * m.dr)) * exp_rm
            : ddisdr;
        const double dtkedthe_neu = (tke_jp1 - tke_jm1) / (rm       * 2.0 * m.dthe);
        const double ddisdthe_neu = (dis_jp1 - dis_jm1) / (rm       * 2.0 * m.dthe);
        const double dtkedphi_neu = (tke_kp1 - tke_km1) / (rmsinthe * 2.0 * m.dphi);
        const double ddisdphi_neu = (dis_kp1 - dis_km1) / (rmsinthe * 2.0 * m.dphi);

        // Local second-order derivatives for diffusion Laplacian (raw, exp_rm applied later)
        const double d2tkedr2   = (m.tke.x[i+1][j][k] - 2.0*m.tke.x[i][j][k] + tke_im1) / (m.dr   * m.dr);
        const double d2tkedthe2 = (tke_jp1 - 2.0*m.tke.x[i][j][k] + tke_jm1) / (m.dthe * m.dthe);
        const double d2tkedphi2 = (tke_kp1 - 2.0*m.tke.x[i][j][k] + tke_km1) / (m.dphi * m.dphi);

        // Cross-diffusion term CD_kω using Neumann-corrected gradients
        const double CD_kw = std::max(
            2.0 * sig_w2 / std::max(m.dis.x[i][j][k], 1.0e-20)
            * (dtkedr_neu * ddisdr_neu + dtkedthe_neu * ddisdthe_neu + dtkedphi_neu * ddisdphi_neu),
            1.0e-20);

        // Blending arguments (all terms dimensionless after substitution):
        //   term1 = sqrt(k*) / (β* ω* y*)
        //   term2 = 500 ν* / (ω* y*²)     where ν* = nue_air_nd
        //   term3 = 4 σ_w2 k* / (CD_kω y*²)
        const double dis_safe = std::max(m.dis.x[i][j][k], 1.0e-20);
        const double arg1 = std::min(
            std::max(sqrt(m.tke.x[i][j][k]) / (bet_star * dis_safe * y_star),
                     500.0 * nue_air_nd / (y_star * y_star * dis_safe)),
            4.0 * sig_w2 * m.tke.x[i][j][k] / (CD_kw * y_star * y_star));
        const double arg2 = std::max(
            2.0 * sqrt(m.tke.x[i][j][k]) / (bet_star * dis_safe * y_star),
            500.0 * nue_air_nd / (y_star * y_star * dis_safe));

        const double F1 = tanh(pow(arg1, 4));
        const double F2 = tanh(pow(arg2, 2));

        m.nue.x[i][j][k] = a1 * m.tke.x[i][j][k]
            / std::max(a1 * dis_safe, Omega * F2);

        const double tke_safe = std::max(m.tke.x[i][j][k], 1.0e-20);
        const double P_k = std::min(m.prod.x[i][j][k],
            20.0 * bet_star * tke_safe * dis_safe);
        const double Y_k = bet_star * tke_safe * dis_safe;
        // P_w uses ω/k directly (Menter 1994 §2.2) — avoids blow-up when ν_T is Ω-limited
        const double P_w = blend(gam1, gam2, F1) * P_k * dis_safe / tke_safe;
        const double Y_w = blend(bet1, bet2, F1) * dis_safe * dis_safe;
        const double D_w = 2.0 * (1.0 - F1) * sig_w2 / dis_safe
            * (dtkedr_neu * ddisdr_neu + dtkedthe_neu * ddisdthe_neu + dtkedphi_neu * ddisdphi_neu);

        // tke_source / dis_source: store only the production–destruction balance.
        // Turbulent diffusion of tke is handled by the RHS (diffusion_tke Laplacian)
        // and is recomputed fresh at every RK4 sub-stage, so it must NOT be added
        // here to avoid double-counting.
        m.tke_source.x[i][j][k] = P_k - Y_k;
        m.dis_source.x[i][j][k] = P_w - Y_w + D_w;
    }
};

#pragma once

#include "MixtureAtm.h"
#include "SaturationH2O.h"
#include "cAtmosphereModel.h"

#include <vector>
#include <cmath>
#include <chrono>
#include <iostream>
#include <cstdio>
#include <algorithm>
#include <cstring>
#include <cstdint>

using namespace AtomUtils;



namespace AtomMoistConvection {
    // PORTED FROM ATHAD 2026-08-19 (its README item 52, its commits of the same day). Every
    // defect below was found there and verified present here by reading THIS file, not by
    // assuming the fork shares it — the family's own rule that a cross-reference is not a
    // check. The magnitudes differ because cp_l does.
    //
    // ATM_MC_S_LEGACY=1 restores the three s-scaling defects together:
    //
    //  (1) UPDRAFT RECURRENCE divided its whole mass-flux bracket by s_0 while the four
    //      sibling recurrences beside it (q_v_u, q_c_u, v_u, w_u) divided by nothing. Every
    //      term in it is already in normalised s, so s_u came out ~2.7e5x too small.
    //  (2) DOWNDRAFT RECURRENCE divided the whole bracket when only L_latent*r_h*e_d needs
    //      it — that term is J/(m^3 s), its neighbours are kg/(m^2 s) times a normalised s.
    //  (3) MC_t CONVERTED s BACK TO KELVIN WITH t_0, which inverts s = cp_l*T/s_0 only if
    //      s_0 == cp_l*t_0. s_0 = 274515.75 = 1005*273.15 is EARTH'S dry-air cp times t_0 and
    //      is correct in ATOM_Precipitation; ATHAD_COND sets cp_l = 1349 and kept it, so the
    //      true factor is s_0/cp_l = 203.50 K, not 273.15, and the convective transport
    //      heating was cp_l*t_0/s_0 = 1.342x too large. (ATHAD, at cp_l = 2040: 2.03x.)
    //
    // MC_t's SECOND term, (L/cp)*conv_src*t_0, is fixed unconditionally below: it was already
    // K/s before the *t_0. Inherited from ATOM_Precipitation, where it is equally wrong.
    static const bool s_consistent = [](){
        const char* e = getenv("ATM_MC_S_LEGACY"); return !(e && atoi(e) != 0); }();

    // ATM_MC_UNBOUNDED_UPDRAFT=1 restores the unbounded updraft scalar recurrence.
    //
    // THE REPAIR IS NOT A CLAMP. d(M phi_u)/dz = E*phi - D*phi_u with dM/dz = E - D is
    // algebraically a convex mixing — expanding gives M dphi_u/dz = E*(phi - phi_u), so
    // DETRAINMENT CANCELS and phi_u cannot leave the interval spanned by its previous value
    // and the environment. Discretely that identity needs the denominator to be exactly the
    // mass the numerator was built from, and three things break it: clamp_M and the
    // is_land / t_00 kills rewrite M(i) after the fact; E_u is a moisture CONVERGENCE and
    // goes negative where the flow diverges moisture, which an entrainment rate cannot; and
    // step*D can exceed M(i-1) on the deep layers. Form the two weights, floor both at zero,
    // divide by their sum. Bounded by construction, identical to the old expression whenever
    // nothing was clamped, and it needs no ceiling constant.
    static const bool updraft_bounded = [](){
        const char* e = getenv("ATM_MC_UNBOUNDED_UPDRAFT"); return !(e && atoi(e) != 0); }();

    constexpr double a_ev = 1.0e-3;
    constexpr double b_ev = 5.9;
    constexpr double t_00 = 236.15;
    constexpr double t_00_minus2 = t_00 - 2.0;
    constexpr double bet = 42.0;
    constexpr double b_u = 0.33;
    constexpr double alf_1 = 0.05;                                      // evaporation rate coefficient [1/s]
    constexpr double alf_2 = 0.05;                                      // evaporation rate coefficient [1/s]
    constexpr double bet_p = 2.0e-3;                                    // in [1/s]
    constexpr double R_cloud = 100.0;                                   // cloud radius in [m] by ECMWF
    double eps_u = 0.2 / R_cloud;                                       // in [1/m] by ECMWF
    double del_u = eps_u;                                               // in [1/m] by ECMWF
//    double eps_u = 1.0e-4;                                              // in [1/m] by COSMO for penetrative (deep) and midlevel clouds
//    double del_u = eps_u;                                               // in [1/m] by COSMO for penetrative (deep) and midlevel clouds
//    constexpr double eps_u = 1.0e-4;                                    // in [1/m] by COSMO for shallow clouds
//    constexpr double del_u = eps_u;                                     // in [1/m] by COSMO for shallow clouds
    constexpr double eps_d = 1.0e-4;                                    // turbulent mixing by COSMO only applied in lower levels for deep clouds
    constexpr double del_d = eps_d;                                     // turbulent mixing by COSMO only applied in lower levels for deep clouds
    constexpr double gam_d = -0.3;                                      // by COSMO
    constexpr double a_u = 3.0e-2;                                      // [·/·]
    constexpr double a_d = 3.0e-2;                                      // [·/·]
    constexpr double vel = 1.1;                                         // [·/·]
    // ATHAD_COND: the convective thresholds are FRACTIONS OF THE LOCAL SURFACE PRESSURE,
    // not absolute hPa.
    //
    // These were 1000/970/900/800/200/700 hPa — Earth surface pressures, and the exact
    // pattern CLAUDE.md warns about: a number that is true at 1013.25 hPa written as a bare
    // literal inside a physics kernel. At ATHAD's 250 bar the whole column sits above
    // 1000 hPa, every trigger fails, and deep convection is simply inactive — wrong, but
    // loudly wrong, and ATHAD documented it as such.
    //
    // At 60 bar it is worse than wrong. p_stat <= 1000 hPa now selects everything above
    // ~53 km, i.e. the top of the saturated column rather than the sea surface, so the
    // triggers do not fail: they fire, in the wrong place, and produce a plausible-looking
    // deep convection scheme keyed to the stratosphere. That is the danger this fork was
    // expected to bring, and it is why these are converted rather than left documented.
    //
    // Fractions are the Earth values divided by 1013.25, so the scheme behaves identically
    // on Earth and travels correctly to any surface pressure.
    constexpr double f_stat_beg        = 1000.0 / 1013.25;              // 0.9869 of p_surf
    constexpr double f_stat_end        =  970.0 / 1013.25;              // 0.9573
    constexpr double f_stat_Cloud_Base =  900.0 / 1013.25;              // 0.8882
    constexpr double f_stat_base       =  800.0 / 1013.25;              // 0.7895
    constexpr double f_stat_diff       =  200.0 / 1013.25;              // 0.1974 (a DEPTH)
    constexpr double f_stat_midlevel   =  700.0 / 1013.25;              // 0.6908, Bechtold (2001)
    constexpr double t_add_u = 0.2;                                     // in [K]
    constexpr double q_v_u_add = 1.0e-4;                                // in [kg/kg]
    constexpr double coeff_recurr = 0.1; 
    constexpr double scale = 0.98;                                      // this value creates deep cloud initial values q_c_u because of t_add_u and t_add_u
    constexpr double q_c_u_max  = 1.0e-2;                               // [kg/kg] physical ceiling on updraft cloud water (~10 g/kg); stops the q_c_u recurrence runaway
    constexpr double w_conv_max = 50.0;                                 // [m/s] physical ceiling on |u_u| convective updraft velocity; stops the safe_r_humid 1/ρ amplification
    // In-loop cap on |M_u|, |M_d| applied at every recurrence write inside the
    // iter_prec loop (E_u contains eps_u·M_u, growth factor ≈ 1+eps_u·step ≈ 1.8
    // per level → ×10⁸ over the column-deep updraft column). rhsForcing alone
    // capped only the final value used in the RHS; the internal c_u, q_v_u,
    // q_v_d, P_conv chain still saw the runaway. Same magnitude as the safe_cap
    // M_max in rhsForcing (~10× any realistic value).
    constexpr double M_max = 3.0;                                       // [kg/(m²s)] (also used in rhsForcing safe_cap)

    // Cloud-base mass-flux coefficient (fix #3, 2026-06-23). The old seed M_u,base = ρ·u
    // used the RESOLVED vertical velocity (~cm/s at the LCL) → M_u≈0.03 kg/m²s, below
    // coeff_recurr(0.1) in 100% of columns → the updraft branch never activated → P_conv≡0
    // ([[project_convective_precip_zero]]). Convective updrafts are SUB-GRID: seed from the
    // parcel buoyancy instead — M_u,base = ρ·c_mb·√(2·CAPE) (convective velocity scale
    // w_b = c_mb·√(2·CAPE)). c_mb≈0.03 gives w_b~0.4–1.3 m/s and M_u~0.4–1.3 kg/m²s for
    // CAPE~100–1000 J/kg (healthy ~0.3, clamp_M ceiling 3.0), and activates the branch for
    // any CAPE≳6 J/kg. Cf. the CAPE-seed prototype in [[project_arabian_coast_precip_spike]].
    constexpr double c_mb = 0.003;                                      // [·] cloud-base mass-flux coefficient (0.03 over-drove: M_u pegged clamp, maxT→57C, stratiform→16mm/d; 10x down → M_u~0.3 healthy)
}
/*
*
*/
// ===================== CALCULATION OF DEEP MOIST CONVECTION =======================================
class MoistConvection {
public:
    explicit MoistConvection(cAtmosphereModel& model)
        : m(model)
        , cloud(model.cloud)
        , R_W_R_A(model.R_WaterVapour / model.R_Air)
        , eps_thermo(1.0 / R_W_R_A)
        , alf((1.0 - eps_thermo) / eps_thermo)
    {}

    int iter;

    void run(int iter) {
        using namespace std;
        using namespace AtomMoistConvection;

        cout << endl << endl << "      MoistConvection  class starts" << endl;
        auto begin = std::chrono::high_resolution_clock::now();

        iter = m.iter_n;

        precompute();
        findCloudBaseLFS();
        surfaceFluxPerturbations();
        initUpdraft();
        CloudBaseUpdate();


        for(int iter_prec = 1; iter_prec <= m.iter_prec; iter_prec++){
            updraftEntrainment(iter_prec);
            CloudBaseUpdate();
            updraftRecurrence();

            downdraftEntrainment(iter_prec);
            downdraftRecurrence();
        }

        rhsForcing();
        subTerrainFill();



        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for MoistConvection  class ended\n", elapsed.count() * 1e-9);
    }

private:
    cAtmosphereModel& m;
    const Array& cloud;
    const double R_W_R_A;
    const double eps_thermo;
    const double alf;


    // Precomputed data
    std::vector<double> height_table;
    std::vector<double> step;
    std::vector<int8_t> land_surf;
    std::vector<int8_t> air_surf;

    // Column indices (local to this run)
    std::vector<std::vector<int>> i_Base_local;
    std::vector<std::vector<int>> i_LFS_local;
    std::vector<std::vector<int>> i_deep_beg_local;
    std::vector<std::vector<int>> i_deep_end_local;
    std::vector<std::vector<double>> M_d_LFS_local;
    std::vector<std::vector<double>> cape_col;        // parcel CAPE per column [m²/s²]  (fix #3, seeds M_u)
    std::vector<std::vector<double>> delta_T_sfp;    // surface-flux T perturbation per column [K]      (Bechtold 2008)
    std::vector<std::vector<double>> delta_q_sfp;    // surface-flux q perturbation per column [kg/kg]  (Bechtold 2008)


// ==================== NaN-SAFE HELPERS ====================
    // Saturation MASS FRACTION, exact conversion.
    //
    // Was q_sat = ep*E/(p - E), the dilute approximation, with the denominator floored at
    // 1e-10 to survive p <= E. That floor is the tell: when E exceeds p the dilute form
    // has no answer, and the floor turned it into an enormous positive number rather than
    // the correct "the column is entirely vapour, q_sat = 1". Water is 67 % of ATHAD's
    // mass, so there is no small parameter to expand in and the exact form is used.
    //
    // M_other is built from the reference CO2 mass fraction but the LOCAL water mass
    // fraction q_v. ATHAD's CO2 is well mixed by construction (co2Atmosphere fills it
    // uniformly), so the reference and the cell agree; water is not, and it is 67 % of
    // the mass, so it decides how much of the non-water carrier is CO2 (heavy) versus
    // background (light) — 28.6 g/mol at q_v = 0 against 35.1 at q_v = 0.67.
    double safe_q_sat(double E_sat, double p_u, double q_v) const noexcept {
        return SaturationH2O::saturationMassFraction(
                   E_sat, p_u, AtmMixture::M_nonwater(q_v, m.co2_0, m.m_comp.M_bg));
    }

    // Clausius-Clapeyron scaling factor for the convective moisture-perturbation
    // ceiling (Part B of the moisture-coupling fix). Returns q_sat(T)/q_sat(T_ref)
    // evaluated at the same pressure, so warmer columns are allowed a proportionally
    // larger updraft-moisture seed and convective precipitation follows CC (~7%/K)
    // instead of falling as q_sat outruns a fixed cap. >1 warm, <1 cold, 1 at T_ref.
    double cc_factor(double T_K, double p_hPa) const {
        // T_ref_cc WAS 288.15 — EARTH'S mean surface temperature as a bare literal, and the
        // ratio this helper returns is taken against it. On Earth that makes the factor O(1)
        // by construction, which is the design ("warmer columns get proportionally more seed
        // moisture, ~7 %/K"). Off Earth it is a ratio between two different planets: it
        // inflated the moisture-seed cap q_v_u_add*cc_factor far above any physical mass
        // fraction (in ATHAD, to 3.47 kg/kg, which is where its impossible q_v_u = 3565 g/kg
        // came from — injected by the seed, not grown by the recurrence). Referencing the
        // model's own surface temperature restores the intended meaning on any planet.
        const double T_ref_cc = 0.5 * (m.t_surf_equator + m.t_surf_pole);
        auto E_sat_of = [&](double T) {
            return (T >= m.t_0)
                ? SaturationH2O::saturationPressure(T)            // over water
                : SaturationH2O::sublimationPressure(T);          // over ice
        };
        // The reference water mass fraction is used for both, so M_other cancels in the
        // ratio this helper returns.
        double qs_T   = safe_q_sat(E_sat_of(T_K),     p_hPa, m.c_0);
        double qs_ref = safe_q_sat(E_sat_of(T_ref_cc), p_hPa, m.c_0);
        return (qs_ref > 0.0) ? qs_T / qs_ref : 1.0;
    }

    // Moist-air density floor: prevents 1/r_humid blow-up when r_humid is
    // clamped to 0 by the limiter in UtilsAtm.
    static double safe_r_humid(double r_h) noexcept {
        return std::max(r_h, 1e-6);
    }

    // In-loop saturation cap on convective mass flux. Applied at every M_u/M_d
    // recurrence write so the eps_u·M_u feedback in E_u (≈ ×1.8 per level)
    // cannot drive c_u/q_v_u/q_v_d/P_conv into runaway before rhsForcing's
    // safe_cap runs at the end of the iter_prec loop. Bit-level NaN/Inf reset
    // (under -ffast-math std::min/max are unreliable on NaN).
    static inline double clamp_M(double v) noexcept {
        std::uint64_t bits;
        std::memcpy(&bits, &v, sizeof(bits));
        if ((bits & 0x7FF0000000000000ULL) == 0x7FF0000000000000ULL) return 0.0;
        const double M_max = AtomMoistConvection::M_max;
        return (v < -M_max) ? -M_max : (v > M_max) ? M_max : v;
    }

// ==================== ENTRAINMENT GRADIENTS (shared helper) ====================
    // Returns the moisture gradients in r, theta and phi directions plus c_ijk = c.x[i][j][k] (floored at 1e-4).
    struct EntrainmentGradients { double dcdr, dcdthe, dcdphi, c_ijk; };  // c_ijk = c.x[i][j][k], floored at 1e-4

    // Neumann (zero-normal-gradient) extrapolation of c into a terrain cell,
    // using the same second-order formula as bcScalarSurfSur:
    //   c_land = (4/3)*c_near_air - (1/3)*c_far_air
    // Called when one neighbor of the centered-difference stencil is terrain.
    // 'near' is the current (air) cell value; 'far' is the cell on the opposite
    // side of the stencil (which may itself be terrain — fall back to 'near' then).
    static inline double neumannC(double c_near, double c_far) noexcept {
        return (4.0/3.0) * c_near - (1.0/3.0) * c_far;
    }

    inline EntrainmentGradients entrainmentGradients(int i, int j, int k) {
        const double costhe       = cos(m.the.z[90] - m.the.z[j]);
        const double height       = m.get_layer_height(i);              // [m]
        const double r_the_height = (m.r_Earth * 1.0e3 + height) * costhe;
        const double U_Earth      = r_the_height * 2.0 * M_PI;          // [m]
        const double step_dr      = step[i];
        const double inv_step_dr  = 0.5 / step_dr;
        const double inv_step_the = 0.5 / (U_Earth / (m.jm-1));         // [1/m]
        const double inv_step_phi = 1.0 / (U_Earth / (m.km-1));
        const int    i_p          = std::min(i+1, m.im-1);
        const int    i_m          = std::max(i-1, 0);

        // Surface-gradient BC (bcScalarSurfSur approach): when a stencil neighbor
        // is terrain, replace its c value with a second-order Neumann extrapolation
        // from the current (air) cell and the cell on the opposite side, so the
        // normal gradient at the land face evaluates to the air-side gradient
        // rather than a spurious c_ocean − c_land=0 difference.
        // If both neighbors are terrain the component is set to zero.
        const double c0 = m.c.x[i][j][k];

        // ---- radial (i) ----
        const bool land_ip = is_land(m.h, i_p, j, k);
        const bool land_im = is_land(m.h, i_m, j, k);
        const double c_ip = land_ip
            ? (land_im ? c0 : neumannC(c0, m.c.x[i_m][j][k]))
            : m.c.x[i_p][j][k];
        const double c_im = land_im
            ? (land_ip ? c0 : neumannC(c0, m.c.x[i_p][j][k]))
            : m.c.x[i_m][j][k];
        const double dcdr = (c_ip - c_im) * inv_step_dr;

        // ---- meridional (j / theta) ----
        const bool land_jp = is_land(m.h, i, j+1, k);
        const bool land_jm = is_land(m.h, i, j-1, k);
        const double c_jp = land_jp
            ? (land_jm ? c0 : neumannC(c0, m.c.x[i][j-1][k]))
            : m.c.x[i][j+1][k];
        const double c_jm = land_jm
            ? (land_jp ? c0 : neumannC(c0, m.c.x[i][j+1][k]))
            : m.c.x[i][j-1][k];
        const double dcdthe = (c_jp - c_jm) * inv_step_the;

        // ---- zonal (k / phi) ----
        const bool land_kp = is_land(m.h, i, j, k+1);
        const bool land_km = is_land(m.h, i, j, k-1);
        const double c_kp = land_kp
            ? (land_km ? c0 : neumannC(c0, m.c.x[i][j][k-1]))
            : m.c.x[i][j][k+1];
        const double c_km = land_km
            ? (land_kp ? c0 : neumannC(c0, m.c.x[i][j][k+1]))
            : m.c.x[i][j][k-1];
        const double dcdphi = (c_kp - c_km) * inv_step_phi;

        return {
            dcdr,
            dcdthe,
            dcdphi,
            std::max(c0, 1e-4)   // floor: was exact ==0.0 check, misses small positives
        };
    }

// ==================== PRECOMPUTE ====================
    void precompute() {
        height_table.resize(m.im);
        step.assign(m.im, 0.0);
        for(int i = 0; i < m.im; i++){
            height_table[i] = m.get_layer_height(i);
        }
        for(int i = 0; i < m.im - 1; i++){
            step[i] = height_table[i+1] - height_table[i];
        }
        if(m.im >= 2) step[m.im-1] = step[m.im-2];

        land_surf.resize(m.jm * m.km);
        air_surf.resize(m.jm * m.km);

        #pragma omp parallel for collapse(2)
        for(int j = 0; j < m.jm; j++){
            for(int k = 0; k < m.km; k++){
                int idx = j * m.km + k;
                land_surf[idx] = is_land(m.h, 0, j, k) ? 1 : 0;
                air_surf[idx]  = is_air(m.h, 0, j, k)  ? 1 : 0;
            }
        }

        i_Base_local.assign(m.jm, std::vector<int>(m.km, 0));
        i_LFS_local.assign(m.jm, std::vector<int>(m.km, 0));
        i_deep_beg_local.assign(m.jm, std::vector<int>(m.km, -1));
        i_deep_end_local.assign(m.jm, std::vector<int>(m.km, -1));
        M_d_LFS_local.assign(m.jm, std::vector<double>(m.km, 0));
        cape_col.assign(m.jm, std::vector<double>(m.km, 0.0));
        delta_T_sfp.assign(m.jm, std::vector<double>(m.km, 0.0));
        delta_q_sfp.assign(m.jm, std::vector<double>(m.km, 0.0));

        m.K_u = std::vector<double>(m.im, 0.0);
        m.K_d = std::vector<double>(m.im, 0.0);

        // Reset 3-D convection fields so stale values from a previous call
        // cannot corrupt M_u/M_d when the saturation
        // condition is no longer met in cells that were active before.
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                for (int i = 0; i < m.im; i++) {
                    // Scalar fields driven by s = cp_l*T/s_0.
                    // s must be computed from current t before any other function reads it;
                    // if left stale, downdraftRecurrence propagates the error downward.
                    // THE REMAINING cp_l USES IN THIS FILE ARE THE s <-> T NORMALISATION AND
                    // MUST NOT BE MADE LOCAL. s = cp_l*T/s_0 and T = s*s_0/cp_l are an exact
                    // inverse pair, so cp_l CANCELS and never reaches a result; making one side
                    // local would break the round trip and turn s into an implicit function of
                    // T. The genuine cp uses here — the ones converting ENERGY to TEMPERATURE,
                    // which do not cancel — are cp_sfc, cp_b, cp_e, cp_u and cp_mc, now local.
                    //
                    // INHERITED, recorded not fixed: s_0 = 274515.75 = 1005 * t_0 carries
                    // EARTH'S dry-air cp, not this fork's cp_l = 1349 (which would give
                    // 368472), while param.py calls s_0 "cp_l * t_0" so it reads as derived.
                    // It changes no result precisely BECAUSE cp_l cancels. The one place the
                    // scale would matter, s = 1.0 as a boundary value (BC_Atm.h), sits behind
                    // `if(!is_land) continue` and is dead under invariant 1.
                    m.s.x[i][j][k]      = m.cp_l * m.t.x[i][j][k] * m.t_0 / m.s_0;

                    // Mass-flux fields — reset so cells outside the active [i_base..i_lfs]
                    // range cannot contribute to rhsForcing flux divergences.
                    m.M_u.x[i][j][k]    = 0.0;
                    m.M_d.x[i][j][k]    = 0.0;

                    // Updraft scalars — stale values corrupt MC_t / MC_q (causes #5, #2, #3).
                    m.s_u.x[i][j][k]    = 0.0;
                    m.q_v_u.x[i][j][k]  = 0.0;
                    m.q_c_u.x[i][j][k]  = 0.0;

                    // Source terms used in rhsForcing conv_src:
                    // (L/cp)*c_u*t_0 ≈ 68 K/s per 1e-4 kg/kg/s; stale c_u drives 1e7 K spikes.
                    m.c_u.x[i][j][k]    = 0.0;
                    m.e_l.x[i][j][k]    = 0.0;
                    m.e_d.x[i][j][k]    = 0.0;
                    m.e_p.x[i][j][k]    = 0.0;
                    m.g_p.x[i][j][k]    = 0.0;

                    // Downdraft scalars — stale s_d / q_v_d corrupt the downdraft recurrence
                    // and the flux divergence in rhsForcing.
                    m.s_d.x[i][j][k]    = 0.0;
                    m.q_v_d.x[i][j][k]  = 0.0;

                    // Convective precipitation — stale P_conv propagates through sigma_p
                    // and e_p into MC_q/PrecipitationConv.
                    m.P_conv.x[i][j][k] = 0.0;
                }
            }
        }

        // Was a second, independent copy of the radial coordinate definition, re-imposing r0/dr on
        // every precompute() call — which would have silently undone any change made at the one
        // real definition site. Routed through that definition instead; the defensive intent is
        // kept, the duplication is not.
        m.initGridCoordinates();
    }
/*
*
*/
// ==================== SURFACE-FLUX PERTURBATIONS (BECHTOLD 2008) ====================
    // Computes per-column T and q perturbations for shallow convection based on surface fluxes.
    // Populates delta_T_sfp / delta_q_sfp for every column, controlled by
    // m.convection_perturbation (set in param.py / config XML):
    //
    //   0 – fixed offsets: δT = t_add_u, δq = q_v_u_add  for all columns
    //   1 – Bechtold (2008) for shallow columns (cloud depth < f_stat_diff * m.p_stat.x[0][j][k]):
    //           δT = α · H_s / (ρ · c_p · w*)
    //           δq = α · E  / (ρ · w*)
    //       where w* = (g/T · z_BL · H_s/(ρ·c_p))^(1/3)  [Deardorff 1970]
    //       and z_BL ≈ cloud-base height as proxy for BL top.
    //       Deep columns and shallow columns with H_s ≤ 0 receive zero perturbation,
    //       suppressing deep convection triggering.
    //
    // Reference: Bechtold et al. (2008) J. Atmos. Sci. 65:3077, Sect. 2b.
    void surfaceFluxPerturbations() {
        using namespace AtomMoistConvection;
        constexpr double alpha_sfp = 0.2;   // perturbation amplitude (Bechtold 2008, ECMWF)

        #pragma omp parallel for collapse(2) schedule(static)
        for(int j = 1; j < m.jm-1; j++){
            for(int k = 1; k < m.km-1; k++){

                if(m.convection_perturbation == 0){
                    // switch 0: fixed offsets for all columns (original behaviour)
                    delta_T_sfp[j][k] = t_add_u;
                    delta_q_sfp[j][k] = q_v_u_add;
                    continue;
                }

                // switch 1: zero by default → deep columns get no perturbation,
                // suppressing deep-convection triggering
                delta_T_sfp[j][k] = 0.0;
                delta_q_sfp[j][k] = 0.0;

                int i_base = i_Base_local[j][k];
                int i_lfs  = i_LFS_local[j][k];
                if(i_base < 0 || i_lfs <= i_base) continue;

                // Surface state (i = 0 is the lowest model level), shared by all branches.
                double T_sfc   = m.t.x[0][j][k] * m.t_0;             // [K]
                double rho_sfc = safe_r_humid(m.r_humid.x[0][j][k]); // [kg/m³]
                // LOCAL mixture cp, not the constant cp_l. Ported from ATHAD e7bd455, and the
                // error is LARGER here than there, because this fork's subject is that water
                // condenses OUT: cp_l = 1349 is exact at the saturated 513 K sea and 1.574x too
                // high above the cold trap, where q_H2O has fallen 0.3464 -> 0.0042 and cp_of
                // reads 857. param.py already calls cp_l "only the fallback"; this makes that
                // true at the call sites. (H_s is Q_sensible_2D, which nothing writes, so these
                // two surface uses multiply zero — converted so the file has one rule.)
                const double cp_sfc = AtmMixture::cp_of(m.c.x[0][j][k], m.co2.x[0][j][k],
                                                        T_sfc, m.m_comp.M_bg);

                // Part B — Clausius-Clapeyron ceiling for the moisture perturbation.
                // The old fixed cap q_v_u_add clipped away the warming signal (q_sat
                // rises ~7%/K but the seed stayed constant, so condensation
                // q_c_u = q_v_u − scale·q_sat FELL with warming → wrong-sign precip).
                // Scaling the cap by q_sat(T_sfc)/q_sat(T_ref) lets warmer columns
                // inject proportionally more updraft moisture so convective precip
                // follows CC. Applied to both the midlevel and shallow moisture seeds.
                double q_cap = q_v_u_add * cc_factor(T_sfc, m.p_stat.x[0][j][k]);

                // Midlevel convection (Bechtold 2001): triggered by free tropospheric
                // instability, not surface fluxes → fixed T offset, CC-scaled moisture.
                // convection_mode controls whether midlevel clouds are non-precipitating
                // (updraftEntrainment), not how they are triggered.
                if(m.p_stat.x[i_base][j][k] < f_stat_midlevel * m.p_stat.x[0][j][k]){
                    delta_T_sfp[j][k] = t_add_u;
                    delta_q_sfp[j][k] = q_cap;
                    continue;
                }

                // Deep convection: leave at zero (Bechtold 2008: α = 0 for deep)
                if((m.p_stat.x[i_base][j][k] - m.p_stat.x[i_lfs][j][k]) >= f_stat_diff * m.p_stat.x[0][j][k])
                    continue;

                // Part A — surface MOISTURE flux from the T-responsive Dalton evaporation
                // (ThermoAtm: Evaporation ∝ E_sat(T_sfc), i.e. Clausius-Clapeyron). The
                // Bechtold branch previously read Q_latent_2D/Q_sensible_2D, which are
                // never populated (always 0) → δq/δT were identically zero for every
                // shallow column, so warming never reached the convective moisture seed.
                // Sensible flux is still unavailable (Q_sensible_2D ≈ 0); it only affects
                // the small δT and the buoyancy flux, so convection is moisture-driven.
                // 1 mm/d over 1 m² = 1 kg/m² per 86400 s.
                double H_s = m.Q_sensible_2D.y[j][k];                         // [W/m²] (≈0)
                double E   = std::max(0.0, m.Evaporation.y[j][k]) / 86400.0;  // [kg/(m²s)]
                if(E <= 0.0 && H_s <= 0.0) continue;                          // no surface forcing

                // BL depth: cloud-base height as proxy for BL top
                double z_BL = height_table[i_base];     // [m]
                if(z_BL <= 0.0) continue;

                // Deardorff (1970) convective velocity scale from the surface VIRTUAL
                // heat flux B = g/T·[H_s/(ρ·c_p) + 0.61·T·E/ρ], so moisture-driven
                // convection is captured even when the sensible flux is unavailable.
                double B_sfc  = (m.g / T_sfc)
                    * (H_s / (rho_sfc * cp_sfc) + 0.61 * T_sfc * E / rho_sfc);  // [m²/s³]
                double w_star = std::cbrt(std::max(0.0, B_sfc * z_BL));          // [m/s]
                if(w_star <= 0.0) continue;

                double inv_rho_w  = 1.0 / (rho_sfc * w_star);
                // Cap at the CC-scaled ceiling: Bechtold gives smaller values for strong
                // BL (w* ~ 2 m/s → δq ~ 5e-6 kg/kg); the cap prevents blow-up when w* is
                // near zero but positive (weak-flux edge case).
                delta_T_sfp[j][k] = std::min(alpha_sfp * H_s * inv_rho_w / cp_sfc, t_add_u);
                delta_q_sfp[j][k] = std::min(alpha_sfp * E  * inv_rho_w, q_cap);
            }
        }
    }
/*
*
*/
// ==================== UPDRAFT INITIAL VALUES ====================
    void initUpdraft() {
        using namespace AtomMoistConvection;

        #pragma omp parallel for collapse(2)
        for(int j = 1; j < m.jm-1; j++){
            for(int k = 1; k < m.km-1; k++){

                int local_i_Base = i_Base_local[j][k];
                int local_i_beg = i_deep_beg_local[j][k];
                int local_i_end = i_deep_end_local[j][k];
                if(local_i_beg == -1 || local_i_end == -1) continue;

                int local_i_LFS = i_LFS_local[j][k];

                // Perturbations set by surfaceFluxPerturbations() for all columns
                const double t_pert = delta_T_sfp[j][k];
                const double q_pert = delta_q_sfp[j][k];

                // Cloud-base mass-flux seed from parcel buoyancy (fix #3): M_u = ρ·c_mb·√(2·CAPE).
                // Replaces the old ρ·u (resolved w) seed that left M_u below coeff_recurr everywhere.
                m.M_u.x[local_i_beg][j][k] = clamp_M(
                    m.r_humid.x[local_i_beg][j][k] * c_mb * std::sqrt(2.0 * cape_col[j][k]));

                for(int i = local_i_beg+1; i <= local_i_end; i++){

                    double t_u     = m.t.x[i][j][k] * m.t_0;
                    double t_u_add = t_u + t_pert;
                    double r_h_i   = m.r_humid.x[i][j][k];
                    double p_u     = m.p_stat.x[i][j][k];

                    double E_sat_add = SaturationH2O::saturationPressureAuto(t_u_add);
                    double q_sat_add = safe_q_sat(E_sat_add, p_u, m.c.x[i][j][k]);


                    m.q_v_u.x[i][j][k]     = m.c.x[i][j][k] + q_pert;
                    if(m.q_v_u.x[i][j][k] >= scale * q_sat_add && t_u >= t_00)
                        m.q_c_u.x[i][j][k] = m.q_v_u.x[i][j][k] - scale * q_sat_add;
                    if(m.q_v_u.x[i][j][k] <= 0.0) m.q_v_u.x[i][j][k] = 0.0;
                    if(m.q_c_u.x[i][j][k] <= 0.0) m.q_c_u.x[i][j][k] = 0.0;


                    // Entrainment
                    auto [dcdr, dcdthe, dcdphi, c_ijk] = entrainmentGradients(i, j, k);

                    if((m.p_stat.x[local_i_Base][j][k]
                        - m.p_stat.x[local_i_LFS][j][k]) >= f_stat_diff * m.p_stat.x[0][j][k]) {
                        m.E_u.x[i][j][k] =                              // local moisture convergence by entrainment
                            - r_h_i / c_ijk
                            * (m.u.x[i][j][k] * dcdr
                            +  m.v.x[i][j][k] * dcdthe
                            +  m.w.x[i][j][k] * dcdphi);
                    } else  m.E_u.x[i][j][k] = 0.0;

                    if(t_u <= t_00_minus2) m.E_u.x[i][j][k] = 0.0;


                    // Detrainment
                    m.D_u.x[i][j][k] = 0.0;                             // no detrainment


                    double d_Mu      = m.E_u.x[i][j][k] - m.D_u.x[i][j][k];// in [kg/(m³s)]   updraft mass flux based on the local moisture convergence
                    m.M_u.x[i][j][k] = clamp_M(m.M_u.x[i-1][j][k] + d_Mu * step[i]);// in [kg/(m³s)]
                    if(is_land(m.h, i, j, k)) m.M_u.x[i][j][k] = 0.0;

                    m.s.x[i][j][k]   = m.cp_l * t_u / m.s_0;
                    m.s_u.x[i][j][k] = m.cp_l * t_u_add / m.s_0;

                    m.u_u.x[i][j][k] = vel * m.u.x[i][j][k];
                    m.v_u.x[i][j][k] = vel * m.v.x[i][j][k];
                    m.w_u.x[i][j][k] = vel * m.w.x[i][j][k];

                }

                M_d_LFS_local[j][k] = clamp_M(gam_d * m.M_u.x[local_i_beg][j][k]);// downdraft mass flux based on the mass flux at cloud base
                m.M_d.x[local_i_LFS][j][k] = M_d_LFS_local[j][k];


                for(int i = local_i_end+1; i <= local_i_LFS-1; i++){

                    double inv_a_u = 1.0 / (a_u * m.u_0);

                    double t_u     = m.t.x[i][j][k] * m.t_0;
                    double r_h_i   = m.r_humid.x[i][j][k];
                    double t_u_add = t_u + t_pert;
                    double p_u     = m.p_stat.x[i][j][k];

                    double E_sat_add = SaturationH2O::saturationPressureAuto(t_u_add);
                    double q_sat_add = safe_q_sat(E_sat_add, p_u, m.c.x[i][j][k]);

                    if(i == local_i_end+1)  m.M_u.x[i-1][j][k] = m.M_u.x[local_i_beg][j][k];

                    m.s.x[i][j][k]   = m.cp_l * t_u / m.s_0;
                    m.s_u.x[i][j][k] = m.cp_l * t_u / m.s_0;

                    m.q_v_u.x[i][j][k]     = m.c.x[i][j][k] + q_pert;
                    if(m.q_v_u.x[i][j][k] >= scale * q_sat_add && t_u_add >= t_00)
                        m.q_c_u.x[i][j][k] = m.q_v_u.x[i][j][k] - scale * q_sat_add;
//                    if(m.q_v_u.x[i][j][k] <= 0.0) m.q_v_u.x[i][j][k] = 0.0;
//                    if(m.q_c_u.x[i][j][k] <= 0.0) m.q_c_u.x[i][j][k] = 0.0;


                    // Entrainment
                    auto [dcdr, dcdthe, dcdphi, c_ijk] = entrainmentGradients(i, j, k);

                    if((m.p_stat.x[local_i_Base][j][k]
                        - m.p_stat.x[local_i_LFS][j][k]) >= f_stat_diff * m.p_stat.x[0][j][k]) {
                        m.E_u.x[i][j][k] =                              // local moisture convergence by entrainment
                            - r_h_i / c_ijk
                            * (m.u.x[i][j][k] * dcdr
                            +  m.v.x[i][j][k] * dcdthe
                            +  m.w.x[i][j][k] * dcdphi);
                    } else  m.E_u.x[i][j][k] = 0.0;

                    if(t_u <= t_00_minus2) m.E_u.x[i][j][k] = 0.0;


                    // Detrainment
                    m.D_u.x[i][j][k] = 0.0;                             // no detrainment


                    double d_Mu = m.E_u.x[i][j][k] - m.D_u.x[i][j][k];  // in [kg/(m³s)]   updraft mass flux based on the local moisture convergence
                    m.M_u.x[i][j][k] = clamp_M(m.M_u.x[i-1][j][k] + d_Mu * step[i]);  // in [kg/(m²s)]
                    if(is_land(m.h, i, j, k)) m.M_u.x[i][j][k] = 0.0;

                    m.u_u.x[i][j][k] = m.u.x[i][j][k] + m.M_u.x[i][j][k] * inv_a_u / safe_r_humid(r_h_i);
                    const double u_u_cap_e = w_conv_max / m.u_0;
                    m.u_u.x[i][j][k] = std::max(-u_u_cap_e, std::min(m.u_u.x[i][j][k], u_u_cap_e));
                    m.v_u.x[i][j][k] = m.v.x[i][j][k];
                    m.w_u.x[i][j][k] = m.w.x[i][j][k];
                }
            }
        }
    }
/*
*
*/
// ==================== FIND CLOUD BASE & LFS, DEEP BEG/END ====================
void findCloudBaseLFS() {                                                                                                                                                      
    using namespace AtomMoistConvection;                                                                                                                                                      
                                                                                                                                                                                              
    // ---- Phase 1: compute raw indices (parallel) ----                                                                                                                                    
    #pragma omp parallel
    {
        std::vector<double> q_sat_col(m.im);

        #pragma omp for collapse(2) schedule(dynamic, 4)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                for (int i = 0; i < m.im; i++) {
                    if (m.p_stat.x[i][j][k] <= f_stat_beg * m.p_stat.x[0][j][k]) {
                        i_deep_beg_local[j][k] = i;
                        m.Deep_beg.x[i][j][k]  = height_table[i];
                        break;
                    }
                }

                for (int i = m.im - 1; i >= 0; i--) {
                    if (m.p_stat.x[i][j][k] >= f_stat_end * m.p_stat.x[0][j][k]) {
                        i_deep_end_local[j][k] = i;
                        m.Deep_end.x[i][j][k]  = height_table[i];
                        break;
                    }
                }



                // Both scans above may find nothing, leaving the index at its -1 sentinel:
                // the first fires only where p_stat <= f_stat_beg * m.p_stat.x[0][j][k] (1000 hPa), which on an
                // atmosphere whose whole column sits above that threshold never happens.
                // The two other consumers of this index (initUpdraft, CloudBaseUpdate)
                // already guard the sentinel; these two loops did not, and indexed
                // t.x[-1][j][k]. On Earth the scan always fires by level 1 (p_surf ~1013 hPa)
                // so the defect is latent there; ATHAD's 250 bar column makes it live.
                // NOTE: skipping the column is only a bounds fix. The convective thresholds
                // f_stat_beg/end/Cloud_Base/base are now FRACTIONS of the local surface
                // pressure (ATHAD_COND), so the scan fires at the sea surface as intended.
                // The bounds guard below is kept: a column can still find nothing.
                if (i_deep_beg_local[j][k] < 0 || i_deep_end_local[j][k] < 0) continue;

                for (int i = i_deep_beg_local[j][k]; i < m.im; i++) {
                    const double t_u     = m.t.x[i][j][k] * m.t_0;
                    const double t_u_add = t_u + t_add_u;
                    const double p_u     = m.p_stat.x[i][j][k];

                    double E_sat_add = SaturationH2O::saturationPressureAuto(t_u_add);

                    q_sat_col[i] = safe_q_sat(E_sat_add, p_u, m.c.x[i][j][k]);
                }

                for (int i = i_deep_beg_local[j][k]; i < m.im; i++) {
                    const double t_u            = m.t.x[i][j][k] * m.t_0;
                    const double cloud_t_weight = 1.0 / (1.0 + std::exp(-(t_u - t_00) / 2.0));

                    m.q_v_u.x[i][j][k] = m.c.x[i][j][k] + q_v_u_add;
                    if(m.q_v_u.x[i][j][k] >= scale * q_sat_col[i] && t_u >= t_00)
                        m.q_c_u.x[i][j][k] = m.q_v_u.x[i][j][k] - scale * q_sat_col[i];
                    if(m.q_v_u.x[i][j][k] <= 0.0) m.q_v_u.x[i][j][k] = 0.0;
                    if(m.q_c_u.x[i][j][k] <= 0.0) m.q_c_u.x[i][j][k] = 0.0;

                    if (m.q_v_u.x[i][j][k] >= scale * q_sat_col[i] && cloud_t_weight > 0.01
                            && m.p_stat.x[i][j][k] <= f_stat_Cloud_Base * m.p_stat.x[0][j][k]) {  // cloud base must be above 900 hPa
                        i_Base_local[j][k]     = i;
                        m.CloudBase.x[i][j][k] = height_table[i] * cloud_t_weight;
                        break;
                    }
                }

                // LFS / cloud top = level of neutral buoyancy (LNB) of a parcel lifted
                // from cloud base conserving equivalent potential temperature θe.
                // History: the old top-down 50/50-mix scan pinned i_LFS at the model top
                // (im-2) → ~30-level entrainment column → eps_u·M_u blow-up; the replacement
                // upward 50/50-mix scan over-corrected, tripping at i_base+1 in 100% of
                // columns (measured 2026-06-23: depth ~15 hPa, zero deep clouds) → E_u and
                // the precip coefficient K_p forced to 0 → P_conv identically 0. Both compared
                // a water-CONTENT mixture, never the lifted parcel's buoyancy. Here: the parcel
                // θe (conserved from base) vs the environment saturation θes(i); the parcel is
                // buoyant where θe_parcel > θes_env. i_LFS = highest contiguous buoyant level
                // = the LNB (physical cloud top): deep where conditionally unstable, shallow
                // otherwise. M_max-clamped M_u (clamp_M) guards the deep entrainment column,
                // so the LNB ceiling does not reintroduce the old top-of-model blow-up.
                {
                    const int i_base_col = i_Base_local[j][k];
                    if (i_base_col > 0) {
                        // ATHAD: the potential-temperature reference is the SURFACE pressure,
                        // not Earth's 1000 hPa — at 250 bar that constant would put the
                        // reference level ~250 scale-pressures above the ground and make
                        // theta_e meaningless. The Poisson exponent uses the local mixture
                        // R and cp, both of which move with composition and temperature here.
                        const double p0    = m.p_0;                           // reference pressure [hPa]
                        const double T_b   = m.t.x[i_base_col][j][k] * m.t_0 + t_add_u;
                        const double kappa = AtmMixture::R_of(m.c.x[i_base_col][j][k],
                                                              m.co2.x[i_base_col][j][k],
                                                              m.m_comp.R_bg)
                                           / AtmMixture::cp_of(m.c.x[i_base_col][j][k],
                                                               m.co2.x[i_base_col][j][k],
                                                               T_b, m.m_comp.M_bg);
                        const double p_b   = m.p_stat.x[i_base_col][j][k];
                        const double qsb   = q_sat_col[i_base_col];
                        const double cp_b  = AtmMixture::cp_of(m.c.x[i_base_col][j][k],
                                                 m.co2.x[i_base_col][j][k], T_b, m.m_comp.M_bg);
                        const double L_b   = (T_b >= m.t_0) ? m.lv : m.ls;
                        const double thetae_parcel =
                            T_b * std::pow(p0 / p_b, kappa)
                                * std::exp(L_b * qsb / (cp_b * T_b));

                        bool became_buoyant = false;
                        int  i_top = i_base_col;                              // never buoyant → shallow
                        double cape = 0.0;                                    // ∫ buoyancy dz over the buoyant ascent [m²/s²]
                        for (int i = i_base_col + 1; i < m.im; i++) {
                            const double T_e = m.t.x[i][j][k] * m.t_0;
                            const double p_e = m.p_stat.x[i][j][k];
                            const double qse = q_sat_col[i];
                            const double cp_e = AtmMixture::cp_of(m.c.x[i][j][k],
                                                    m.co2.x[i][j][k], T_e, m.m_comp.M_bg);
                            const double L_e = (T_e >= m.t_0) ? m.lv : m.ls;
                            const double thetaes_env =
                                T_e * std::pow(p0 / p_e, kappa)
                                    * std::exp(L_e * qse / (cp_e * T_e));
                            if (thetae_parcel - thetaes_env > 0.0) {
                                became_buoyant = true;
                                i_top = i;                                    // raise the LNB
                                // buoyancy b = g·(θe_parcel − θes_env)/θes_env [m/s²]; CAPE = ∫ b dz
                                const double dz = height_table[i] - height_table[i-1];
                                cape += m.g * (thetae_parcel - thetaes_env) / thetaes_env * dz;
                            } else if (became_buoyant) {
                                break;                                        // lost buoyancy after LFC → LNB reached
                            }
                        }
                        const int i_lnb = std::min(std::max(i_top, i_base_col + 1), m.im - 2);
                        i_LFS_local[j][k]                 = i_lnb;
                        cape_col[j][k]                    = std::max(0.0, cape);
                        m.LevelFreeSinking.x[i_lnb][j][k] = height_table[i_lnb];
                    }
                }
            }
        }
    }

    for(int j = 0; j < m.jm; j++){
        for(int k = 0; k < m.km; k++){
            for(int i = 0; i < m.im; i++){
                if(is_land(m.h, i, j, k)){
                    m.Deep_beg.x[i][j][k]         = 0.0;
                    m.Deep_end.x[i][j][k]         = 0.0;
                    m.CloudBase.x[i][j][k]        = 0.0;
                    m.LevelFreeSinking.x[i][j][k] = 0.0;
                }
            }
        }
    }
}

// ==================== UPDRAFT ENTRAINMENT ====================
    void updraftEntrainment(int iter_prec) {
        using namespace AtomMoistConvection;

        #pragma omp parallel for collapse(2) schedule(static)
        for(int j = 1; j < m.jm-1; j++){
            for(int k = 1; k < m.km-1; k++){

                int i_base = i_Base_local[j][k];
                int i_lfs  = i_LFS_local[j][k];

                int idx_surf = j * m.km + k;
                bool is_land_surf_jk = land_surf[idx_surf];

                if(i_base < 0 || i_lfs <= i_base) continue;

                // Shallow and midlevel convection are non-precipitating:
                // set delta_i_c beyond any reachable cloud depth so K_p stays zero.
                // Midlevel: cloud base above f_stat_midlevel * p_surf (0.69), Bechtold (2001).
                bool is_shallow   = m.convection_mode >= 1
                                    && (m.p_stat.x[i_base][j][k]
                                        - m.p_stat.x[i_lfs][j][k]) < f_stat_diff * m.p_stat.x[0][j][k];
                bool is_midlevel  = m.convection_mode >= 2
                                    && m.p_stat.x[i_base][j][k] < f_stat_midlevel * m.p_stat.x[0][j][k];
                double delta_i_c  = (is_shallow || is_midlevel) ? 1.0e9
//                                  : (is_land_surf_jk ? 3000.0 : 1500.0); // Tiedtke (1989)
                                  : (is_land_surf_jk ? 1000.0 : 500.0); // Bechtold (2001)
                double height_base = height_table[i_base];

                // Cause #1: initUpdraft left non-zero M_u in [0..i_base-1].  Those cells
                // are never touched by this function, so their values create a flux
                // discontinuity at i_base in rhsForcing.  Zero them here on every
                // iter_prec call so the chain below i_base is clean.
                for (int ii = 0; ii < i_base; ii++)
                    m.M_u.x[ii][j][k] = 0.0;

                m.M_u.x[i_base][j][k] = clamp_M(                        // fix #3: CAPE-based seed (was ρ·u)
                    m.r_humid.x[i_base][j][k] * c_mb * std::sqrt(2.0 * cape_col[j][k]));  // [kg/(m²s)]
                if(is_land(m.h, i_base, j, k)) m.M_u.x[i_base][j][k] = 0.0;

                // Causes #2 and #3: updraftRecurrence starts at i_base+1 and reads
                // s_u[i_base] and q_v_u[i_base] as i-1 seed values, but no function
                // ever assigns them at exactly i_base.  Set them here.
                {
                    const double t_base = m.t.x[i_base][j][k] * m.t_0;
                    const double t_pert = delta_T_sfp[j][k];
                    const double q_pert = delta_q_sfp[j][k];
                    m.s_u.x[i_base][j][k]   = m.cp_l * (t_base + t_pert) / m.s_0;
                    m.q_v_u.x[i_base][j][k] = std::max(m.c.x[i_base][j][k] + q_pert, 0.0);
                }

                M_d_LFS_local[j][k] = clamp_M(gam_d * m.M_u.x[i_base][j][k]);

                for(int i = i_base+1; i <= i_lfs; i++){

                    double height_i = height_table[i];
                        double t_u      = m.t.x[i][j][k] * m.t_0;
                        double r_h_i    = m.r_humid.x[i][j][k];

                        // Entrainment
                        auto [dcdr, dcdthe, dcdphi, c_ijk] = entrainmentGradients(i, j, k);

                        if((m.p_stat.x[i_base][j][k]
                            - m.p_stat.x[i_lfs][j][k]) >= f_stat_diff * m.p_stat.x[0][j][k]) {
                            m.E_u.x[i][j][k] =                          // local moisture convergence by entrainment
                                - r_h_i / c_ijk
                                * (m.u.x[i][j][k] * dcdr
                                +  m.v.x[i][j][k] * dcdthe
                                +  m.w.x[i][j][k] * dcdphi)
                                + eps_u * m.M_u.x[i][j][k];
                        } else  m.E_u.x[i][j][k] = 0.0;

                        if(t_u <= t_00_minus2) m.E_u.x[i][j][k] = 0.0;



                        // Detrainment
                        // Turbulent detrainment: applies on all iterations below LFS
                        m.D_u.x[i][j][k] = (i <= i_lfs)
                            ? eps_u * m.M_u.x[i][j][k] : 0.0;           // turbulent detrainment [kg/(m³s)]
                            

                        if(iter_prec > 1 && i == i_lfs)
                            m.D_u.x[i][j][k] = del_u * m.M_u.x[i][j][k] // turbulent detrainment
                                + b_u * m.M_u.x[i][j][k]/step[i];       // dynamic detrainment

                        if(iter_prec > 1 && i == i_lfs-1)
                            m.D_u.x[i][j][k] = del_u * m.M_u.x[i][j][k] // turbulent detrainment
                                + (1.0 - b_u) * m.M_u.x[i][j][k]/step[i];// dynamic detrainment

                        if(t_u <= t_00_minus2) m.D_u.x[i][j][k] = 0.0;


                        double d_Mu = m.E_u.x[i][j][k] - m.D_u.x[i][j][k];// in [(kg/kg)/s]
                        m.M_u.x[i][j][k] = clamp_M(m.M_u.x[i-1][j][k] + d_Mu * step[i]);// in [(kg/kg)/s]
                        if(is_land(m.h, i, j, k)) m.M_u.x[i][j][k] = 0.0;


                        // Condensation rate c_u is now set in updraftRecurrence (fix #2) as the
                        // ACTUAL parcel condensation increment, replacing the old circular
                        // c_u = q_c_u·M_u/(r_h·step) that double-counted latent heating. Zero here.
                        m.c_u.x[i][j][k] = 0.0;

                        // Evaporation of detrained cloud water in updraft
                        if(iter_prec == 1)
                            m.e_l.x[i][j][k] = 0.0;
                        else{
                            const double r_h_floor = std::max(r_h_i, 0.01);
                            m.e_l.x[i][j][k] = m.D_u.x[i][j][k] * m.q_c_u.x[i][j][k] / r_h_floor;// in [(kg/kg)/s]
                            if(m.e_l.x[i][j][k] <= 0.0) m.e_l.x[i][j][k] = 0.0;
                        }

                        // Precipitation formation
                        double K_p = (height_i > height_base + delta_i_c) ? bet_p : 0.0;

                        if(t_u >= m.t_0){
                            m.g_p.x[i][j][k] = (iter_prec == 1) ? 0.0 : K_p * m.q_c_u.x[i][j][k];// in [(kg/kg)/s]
                        }else{
                            m.g_p.x[i][j][k] = 0.0;
                        }
                } // end i

                // Zero M_u at the LFS so rhsForcing sees no updraft flux jump at the
                // cloud top: the updraft fully detrains at this level.
                if (i_base >= 0 && i_lfs > i_base)
                    m.M_u.x[i_lfs][j][k] = 0.0;

           } // end k
        } // end j
    }
/*
*
*/
// ==================== DOWNDRAFT ENTRAINMENT ====================
    void downdraftEntrainment(int iter_prec) {
        using namespace AtomMoistConvection;

        #pragma omp parallel for collapse(2) schedule(static)
        for(int j = 0; j < m.jm; j++){
            for(int k = 0; k < m.km; k++){

                int i_base = i_Base_local[j][k];
                int i_lfs  = i_LFS_local[j][k];
                // Cause #7: without this guard, columns with no detected cloud
                // (i_lfs still at its initialised value 0) still set M_d[0],
                // which rhsForcing then applies to every row/column.
                if(i_lfs <= 0 || i_lfs <= i_base) continue;

                double t_u = m.t.x[i_lfs][j][k] * m.t_0;
                double p_u = m.p_stat.x[i_lfs][j][k];

                m.M_d.x[i_lfs][j][k] = clamp_M(M_d_LFS_local[j][k]);

                m.E_d.x[i_lfs][j][k] = eps_d * fabs(m.M_d.x[i_lfs][j][k]);
                m.D_d.x[i_lfs][j][k] = del_d * fabs(m.M_d.x[i_lfs][j][k]);

                double E_sat = SaturationH2O::saturationPressureAuto(t_u);
                double q_sat = safe_q_sat(E_sat, p_u, m.c.x[i_lfs][j][k]);

                m.q_v_d.x[i_lfs][j][k] = 0.5 * (cloud.x[i_lfs][j][k] + scale * q_sat);

                // ATHAD: local mixture gas constant, not R_Air with a virtual-temperature
                // expansion that assumes water is a trace. Matches ThermoAtm::densities().
                double r_humid_lfs = 1e2 * p_u
                    / (AtmMixture::R_of(m.c.x[i_lfs][j][k], m.co2.x[i_lfs][j][k],
                                        m.m_comp.R_bg) * t_u);
                double inv_a_d = 1.0 / (a_d * m.u_0);

                m.u_d.x[i_lfs][j][k] = m.u.x[i_lfs][j][k]
                    + m.M_d.x[i_lfs][j][k] * inv_a_d / safe_r_humid(r_humid_lfs);

                m.v_d.x[i_lfs][j][k] = m.v.x[i_lfs][j][k];
                m.w_d.x[i_lfs][j][k] = m.w.x[i_lfs][j][k];
                m.s_d.x[i_lfs][j][k] = m.s.x[i_lfs][j][k];

                m.P_conv.x[i_lfs][j][k] = 0.0;

                // Nordeng (1994): convective precipitation area fraction σ = M_u / (ρ·w_u)
                // evaluated at cloud base; w_u is in model units → multiply by u_0 for [m/s]
                double sigma_p = 0.0;
                {
                    double w_u_cb = fabs(m.w_u.x[i_base][j][k]) * m.u_0;
                    if(w_u_cb > 0.0)
                        sigma_p = fabs(m.M_u.x[i_base][j][k])
                                  / (safe_r_humid(m.r_humid.x[i_base][j][k]) * w_u_cb);
                    sigma_p = std::min(sigma_p, 1.0);
                }

                for(int i = i_lfs-1; i >= 0; i--){

                    p_u = m.p_stat.x[i][j][k];
                    t_u = m.t.x[i][j][k] * m.t_0;

                    double r_h_i = m.r_humid.x[i][j][k];

                    m.s_d.x[i][j][k] = m.cp_l * t_u / m.s_0;

                    double E_sat = SaturationH2O::saturationPressureAuto(t_u);
                    double q_sat = safe_q_sat(E_sat, p_u, m.c.x[i][j][k]);

                    m.q_v_d.x[i][j][k] = 0.5 * (cloud.x[i][j][k] + scale * q_sat);

                    m.E_d.x[i][j][k] = eps_d * fabs(m.M_d.x[i][j][k]);
                    m.D_d.x[i][j][k] = del_d * fabs(m.M_d.x[i][j][k]);

                    double d_Md      = m.E_d.x[i][j][k] - m.D_d.x[i][j][k];  // in [kg/(m³s)]
                    m.M_d.x[i][j][k] = clamp_M(m.M_d.x[i+1][j][k] - d_Md * step[i]);  // in [kg/(m²s)]
                    if(is_land(m.h, i, j, k)) m.M_d.x[i][j][k] = 0.0;

                    double inv_a_d   = 1.0 / (a_d * m.u_0);
                    m.u_d.x[i][j][k] = m.u.x[i][j][k]
                        + m.M_d.x[i][j][k] * inv_a_d / safe_r_humid(r_h_i);
                    m.v_d.x[i][j][k] = m.v.x[i][j][k];
                    m.w_d.x[i][j][k] = m.w.x[i][j][k];

                    // Evaporation below cloud base
                    if(i <= i_base){
                        if(iter_prec == 1)
                          m.e_p.x[i][j][k] = 0.0;
                        else{
                            double precip_i = m.P_conv.x[i][j][k];
                            if(precip_i > 0.0 && sigma_p > 0.0)
                                m.e_p.x[i][j][k] = sigma_p * alf_1      // in [(kg/kg)/s]
                                    * (q_sat - m.q_v_d.x[i][j][k])
                                    * sqrt(1e1 / alf_2 * precip_i / sigma_p);
                            if(m.e_p.x[i][j][k] <= 0.0) m.e_p.x[i][j][k] = 0.0;
                        }
                    }

                    // Evaporation within downdraft
                    if(iter_prec == 1)
                        m.e_d.x[i][j][k] = 0.0;
                    else{
                        double precip_i = m.P_conv.x[i][j][k];
                        if(precip_i > 0.0){
                          double P_16 = pow(precip_i, 1.0/6.0);
                          double P_49 = pow(precip_i, 4.0/9.0);
                          m.e_d.x[i][j][k] = a_ev * (1.0 + b_ev * P_16) // in [(kg/kg)/s]
                              * (q_sat - m.q_v_d.x[i][j][k]) * P_49;
                        }else{
                          m.e_d.x[i][j][k] = 0.0;
                        }
                        if(m.e_d.x[i][j][k] <= 0.0) m.e_d.x[i][j][k] = 0.0;
                    }

                } // end i downdraft
            }
        }
    }
/*
*
*/
// ==================== UPDRAFT MAIN RECURRENCE ====================
    void updraftRecurrence() {
        using namespace AtomMoistConvection;

        #pragma omp parallel for collapse(2) schedule(static)
        for(int j = 1; j < m.jm-1; j++){
            for(int k = 1; k < m.km-1; k++){

                int i_base = i_Base_local[j][k];
                int i_LFS  = i_LFS_local[j][k];
                if(i_base < 0 || i_LFS <= i_base) continue;

                for(int i = i_base+1; i <= i_LFS-1; i++){

                    double t_u       = m.t.x[i-1][j][k] * m.t_0;
                    double r_h_prev  = m.r_humid.x[i-1][j][k];
                    double step_prev = step[i-1];
                    double r_h_i     = m.r_humid.x[i][j][k];

                    double d_Mu = m.E_u.x[i-1][j][k] - m.D_u.x[i-1][j][k];

                    m.M_u.x[i][j][k] = clamp_M(m.M_u.x[i-1][j][k] + d_Mu * step_prev);
                    if(is_land(m.h, i, j, k)) m.M_u.x[i][j][k] = 0.0;
                    if(t_u <= t_00_minus2) m.M_u.x[i][j][k] = 0.0;

                    double M_u_prev = m.M_u.x[i-1][j][k];

                    // Condensation thermodynamics (vapour removal AND latent heating) are now
                    // handled entirely by the parcel saturation adjustment (fix #2) below, using
                    // the ACTUAL condensation increment. The old c_u-based terms here used the
                    // circular c_u = q_c_u·M_u/(r_h·step) and double-counted the heating, driving
                    // a slow low-level thermal runaway (maxT→48-57°C). Transport only here:
                    double dummy_q_v_u = M_u_prev * m.q_v_u.x[i-1][j][k]
                        + step_prev * (m.E_u.x[i-1][j][k] * m.c.x[i-1][j][k]
                        - m.D_u.x[i-1][j][k] * m.q_v_u.x[i-1][j][k]);

                    double dummy_q_c_u = M_u_prev * m.q_c_u.x[i-1][j][k]
                        + step_prev * (-m.D_u.x[i-1][j][k] * cloud.x[i-1][j][k]
                        - r_h_prev * m.g_p.x[i-1][j][k]);

                    double dummy_vel_v_u = M_u_prev * m.v_u.x[i-1][j][k]
                        + step_prev * m.E_u.x[i-1][j][k] * m.v.x[i-1][j][k];
                    double dummy_vel_w_u = M_u_prev * m.w_u.x[i-1][j][k]
                        + step_prev * m.E_u.x[i-1][j][k] * m.w.x[i-1][j][k];

                    // The /s_0 is spurious: every term here is already in normalised s. The
                    // four recurrences above carry no such division, which is the tell.
                    double dummy_s_u = (M_u_prev * m.s_u.x[i-1][j][k]
                        + step_prev * (m.E_u.x[i-1][j][k] * m.s.x[i-1][j][k]
                        - m.D_u.x[i-1][j][k] * m.s_u.x[i-1][j][k]));
                    if(!s_consistent) dummy_s_u /= m.s_0;

                    // The consistent denominator: the SAME mass the numerator was built from,
                    // as two non-negative weights (parcel air retained, environment air
                    // entrained). See ATM_MC_UNBOUNDED_UPDRAFT at the top of this file.
                    const double w_env = std::max(step_prev * m.E_u.x[i-1][j][k], 0.0);
                    const double w_par = std::max(M_u_prev - step_prev * m.D_u.x[i-1][j][k], 0.0);
                    const double M_u_mix = w_par + w_env;

                    double M_u_i = m.M_u.x[i][j][k];
                    // BOTH tests, deliberately: the repair changes the DENOMINATOR, not which
                    // cells are active. |M_u_i| > coeff_recurr is the shipped activation test
                    // and it stays, so the environment fallback fires exactly where it fired
                    // before. (ATHAD measured what dropping it costs: max MC_t 8.3e-5 -> its
                    // 0.01 cap at iteration 10. Widening a scheme's active set is a change of
                    // physics and does not belong in a units repair.)
                    if(updraft_bounded && M_u_mix > coeff_recurr
                       && fabs(M_u_i) > coeff_recurr){
                        const double inv_mix = 1.0 / M_u_mix;
                        // q_c_u keeps the flux form: its bracket carries genuine SOURCES
                        // (detrained cloud water, precipitation conversion), so it is not a
                        // convex combination of anything. It is still divided by the
                        // consistent mass, which is what removes the M_max/|M| amplification.
                        m.q_v_u.x[i][j][k] = (w_par * m.q_v_u.x[i-1][j][k]
                                            + w_env * m.c.x[i-1][j][k]) * inv_mix;
                        m.q_c_u.x[i][j][k] = dummy_q_c_u * inv_mix;
                        m.v_u.x[i][j][k]   = (w_par * m.v_u.x[i-1][j][k]
                                            + w_env * m.v.x[i-1][j][k]) * inv_mix;
                        m.w_u.x[i][j][k]   = (w_par * m.w_u.x[i-1][j][k]
                                            + w_env * m.w.x[i-1][j][k]) * inv_mix;
                        m.s_u.x[i][j][k]   = s_consistent
                            ? (w_par * m.s_u.x[i-1][j][k] + w_env * m.s.x[i-1][j][k]) * inv_mix
                            : dummy_s_u / M_u_i;
                    } else if(!updraft_bounded && fabs(M_u_i) > coeff_recurr){
                        double inv_M_u = 1.0 / M_u_i;
                        m.q_v_u.x[i][j][k] = dummy_q_v_u * inv_M_u;
                        m.q_c_u.x[i][j][k] = dummy_q_c_u * inv_M_u;
                        m.v_u.x[i][j][k]   = dummy_vel_v_u * inv_M_u;
                        m.w_u.x[i][j][k]   = dummy_vel_w_u * inv_M_u;
                        m.s_u.x[i][j][k]   = dummy_s_u * inv_M_u;
                    } else {
                        // When M_u ≤ coeff_recurr the updraft is negligible; set
                        // updraft scalars to their environment counterparts so every
                        // flux term M_u*(scalar_u - env) evaluates to ~0 in rhsForcing.
                        // Setting s_u=0 with s≈0.8 gave flux=M_u*(0-0.8)→13M °C blow-up.
                        m.s_u.x[i][j][k]   = m.s.x[i][j][k];
                        m.q_v_u.x[i][j][k] = m.c.x[i][j][k];
                        m.q_c_u.x[i][j][k] = 0.0;
                        m.v_u.x[i][j][k]   = m.v.x[i][j][k];
                        m.w_u.x[i][j][k]   = m.w.x[i][j][k];
                    }

                    // Updraft condensation source (fix #2, 2026-06-23). The q_c_u recurrence
                    // above only CARRIES and LOSES condensate (detrainment, precip) — it has
                    // NO in-ascent source, because the c_u term that used to feed it is
                    // circularly defined (c_u = q_c_u·M_u/(r_h·step)), so removing it for the
                    // "double-count" left q_c_u seeded only at cloud base and depleting upward
                    // → q_c_u≈0 at the K_p precip levels → P_conv≡0 ([[project_convective_precip_zero]]).
                    // Genuine source = the rising parcel condensing as it cools: saturation-
                    // adjust the parcel here, moving q_v_u in excess of q_sat(T_u,p) into q_c_u
                    // and releasing latent heat into s_u. Adaptive 1/(1+G) damping (G = latent
                    // gain L/cp·dq_sat/dT) stops the warm-cell overshoot (cf. SaturationAdjustment).
                    m.c_u.x[i][j][k] = 0.0;
                    if(fabs(m.M_u.x[i][j][k]) > coeff_recurr){
                        double T_u       = m.s_u.x[i][j][k] * m.s_0 / m.cp_l;   // parcel temp [K]
                        const double p_u = m.p_stat.x[i][j][k];
                        double dcond_tot = 0.0;
                        for(int it = 0; it < 2; ++it){
                            const double E_s     = SaturationH2O::saturationPressureAuto(T_u);
                            const double q_sat_u = safe_q_sat(E_s, p_u, m.c.x[i][j][k]);
                            const double dq      = m.q_v_u.x[i][j][k] - q_sat_u;
                            if(dq <= 0.0) break;
                            const double cp_u  = AtmMixture::cp_of(m.c.x[i][j][k],
                                                    m.co2.x[i][j][k], T_u, m.m_comp.M_bg);
                            const double L_u   = (T_u >= m.t_0) ? m.lv : m.ls;
                            // ATHAD_COND: exact dq_sat/dT, which needs the pressure and the
                            // non-water molar mass rather than q_sat alone. See
                            // SaturationH2O::satDerivFactor — the dilute form this replaces
                            // is the x -> 0 limit and is a third too small at 56 % water.
                            const double M_o_u = AtmMixture::M_nonwater(m.c.x[i][j][k],
                                                     m.co2.x[i][j][k], m.m_comp.M_bg);
                            const double dqsdT = SaturationH2O::dqSatdT(T_u, p_u, M_o_u);
                            const double G     = (L_u / cp_u) * dqsdT;          // latent gain
                            const double dcond = dq / (1.0 + G);               // damped condensation
                            m.q_v_u.x[i][j][k] -= dcond;
                            m.q_c_u.x[i][j][k] += dcond;
                            T_u                += (L_u / cp_u) * dcond;         // latent heating (sole source)
                            dcond_tot          += dcond;
                        }
                        m.s_u.x[i][j][k] = m.cp_l * T_u / m.s_0;
                        // c_u = the ACTUAL condensation rate [(kg/kg)/s] (was the circular
                        // q_c_u·M_u/(r_h·step)); feeds only the environment forcing conv_src
                        // (MC_t heating / MC_q drying) in rhsForcing — no longer the s_u/q_v_u budgets.
                        m.c_u.x[i][j][k] = dcond_tot * m.M_u.x[i][j][k]
                                           / (std::max(r_h_i, 0.01) * step[i]);
                    }

                    double inv_a_u = 1.0 / (a_u * m.u_0);
                    m.u_u.x[i][j][k] = m.u.x[i][j][k] + m.M_u.x[i][j][k] * inv_a_u / safe_r_humid(r_h_i);
                    const double u_u_cap = w_conv_max / m.u_0;
                    m.u_u.x[i][j][k] = std::max(-u_u_cap, std::min(m.u_u.x[i][j][k], u_u_cap));

                    if(m.q_v_u.x[i][j][k] <= 0.0) m.q_v_u.x[i][j][k] = 0.0;
                    if(m.q_c_u.x[i][j][k] <= 0.0) m.q_c_u.x[i][j][k] = 0.0;
                    if(m.q_c_u.x[i][j][k] > q_c_u_max) m.q_c_u.x[i][j][k] = q_c_u_max;
                    if(m.s_u.x[i][j][k] <= 0.0)   m.s_u.x[i][j][k]   = 0.0;

                    if(is_land(m.h, i, j, k) || t_u <= t_00){
                        m.E_u.x[i][j][k]   = 0.0; m.D_u.x[i][j][k] = 0.0; m.M_u.x[i][j][k] = 0.0;
                        m.q_v_u.x[i][j][k] = 0.0; m.q_c_u.x[i][j][k] = 0.0;
                        m.u_u.x[i][j][k]   = 0.0; m.v_u.x[i][j][k] = 0.0; m.w_u.x[i][j][k] = 0.0;
                        m.s_u.x[i][j][k]   = 1.0;
                        m.g_p.x[i][j][k]   = 0.0; m.c_u.x[i][j][k] = 0.0; m.e_l.x[i][j][k] = 0.0;
                    }
                } // end i
            }
        }
    }
/*
*
*/
// ==================== DOWNDRAFT MAIN RECURRENCE ====================
    void downdraftRecurrence() {
        using namespace AtomMoistConvection;

        #pragma omp parallel for collapse(2) schedule(static)
        for(int j = 0; j < m.jm; j++){
            for(int k = 0; k < m.km; k++){

                int i_lfs_loc = i_LFS_local[j][k];

                for(int i = i_lfs_loc-1; i >= 0; i--){

                    double t_u      = m.t.x[i][j][k] * m.t_0;
                    double step_ip1 = step[i+1];
                    double r_h_ip1  = m.r_humid.x[i+1][j][k];
                    double r_h_i    = m.r_humid.x[i][j][k];

                    double d_Md      = m.E_d.x[i][j][k] - m.D_d.x[i][j][k];
                    m.M_d.x[i][j][k] = clamp_M(m.M_d.x[i+1][j][k] - d_Md * step_ip1);
                    if(t_u <= t_00_minus2) m.M_d.x[i][j][k] = 0.0;

                    double L_latent = (t_u >= m.t_0) ? m.lv : m.ls;
                    double M_d_ip1  = m.M_d.x[i+1][j][k];

                    double dummy_q_v_d = M_d_ip1 * m.q_v_d.x[i+1][j][k]
                        - step_ip1 * (m.E_d.x[i+1][j][k] * m.c.x[i+1][j][k]
                        - m.D_d.x[i+1][j][k] * m.q_v_d.x[i+1][j][k]
                        + r_h_ip1 * m.e_d.x[i+1][j][k]);

                    double dummy_v_d = M_d_ip1 * m.v_d.x[i+1][j][k]
                        - step_ip1 * (m.E_d.x[i+1][j][k] * m.v.x[i+1][j][k]
                        - m.D_d.x[i+1][j][k] * m.v_d.x[i+1][j][k]);

                    double dummy_w_d = M_d_ip1 * m.w_d.x[i+1][j][k]
                        - step_ip1 * (m.E_d.x[i+1][j][k] * m.w.x[i+1][j][k]
                        - m.D_d.x[i+1][j][k] * m.w_d.x[i+1][j][k]);

                    // Only the L_latent term needs the /s_0 — it is J/(m^3 s) where its two
                    // neighbours are kg/(m^2 s) times an already-normalised s. The shipped
                    // form divided all three. See ATM_MC_S_LEGACY at the top of this file.
                    double dummy_s_d;
                    if(s_consistent){
                        dummy_s_d = M_d_ip1 * m.s_d.x[i+1][j][k]
                            - step_ip1 * (m.E_d.x[i+1][j][k] * m.s.x[i+1][j][k]
                            - m.D_d.x[i+1][j][k] * m.s_d.x[i+1][j][k]
                            - L_latent * r_h_ip1 * m.e_d.x[i+1][j][k] / m.s_0);
                    } else {
                        dummy_s_d = (M_d_ip1 * m.s_d.x[i+1][j][k]
                            - step_ip1 * (m.E_d.x[i+1][j][k] * m.s.x[i+1][j][k]
                            - m.D_d.x[i+1][j][k] * m.s_d.x[i+1][j][k]
                            - L_latent * r_h_ip1 * m.e_d.x[i+1][j][k])) / m.s_0;
                    }

                    double M_d_i = m.M_d.x[i][j][k];
                    if(fabs(M_d_i) > coeff_recurr){
                        double inv_M_d = 1.0 / M_d_i;
                        m.q_v_d.x[i][j][k] = dummy_q_v_d * inv_M_d;
                        m.v_d.x[i][j][k]   = dummy_v_d   * inv_M_d;
                        m.w_d.x[i][j][k]   = dummy_w_d   * inv_M_d;
                        m.s_d.x[i][j][k]   = dummy_s_d   * inv_M_d;
                    } else {
                        // Mirrors the updraft else-branch above. When |M_d| ≤ coeff_recurr
                        // the downdraft is negligible, but skipping the assignment leaves
                        // stale q_v_d / s_d from the previous MoistConvection call. The
                        // next iteration then feeds that stale value through the D_d*q_v_d
                        // term in dummy_q_v_d and through the M_d*(q_v_d-c) flux in
                        // rhsForcing — amplifying it by the M_d ratio across marginal cells
                        // (observed blow-up at 1°N 119°E, 0–38 m once spin-up was off).
                        m.q_v_d.x[i][j][k] = m.c.x[i][j][k];
                        m.v_d.x[i][j][k]   = m.v.x[i][j][k];
                        m.w_d.x[i][j][k]   = m.w.x[i][j][k];
                        m.s_d.x[i][j][k]   = m.s.x[i][j][k];
                    }

                    double inv_a_d = 1.0 / (a_d * m.u_0);
                    m.u_d.x[i][j][k] = m.u.x[i][j][k] + m.M_d.x[i][j][k] * inv_a_d / safe_r_humid(r_h_i);

                    if(m.q_v_d.x[i][j][k] <= 0.0) m.q_v_d.x[i][j][k] = 0.0;
                    if(m.s_d.x[i][j][k] <= 0.0)   m.s_d.x[i][j][k]   = 0.0;

                    if(is_land(m.h, i, j, k) || t_u <= t_00){
                        m.E_d.x[i][j][k]   = 0.0; m.D_d.x[i][j][k] = 0.0; m.M_d.x[i][j][k] = 0.0;
                        m.q_v_d.x[i][j][k] = 0.0;
                        m.u_d.x[i][j][k]   = 0.0; m.v_d.x[i][j][k] = 0.0; m.w_d.x[i][j][k] = 0.0;
                        m.s_d.x[i][j][k]   = 1.0; m.e_d.x[i][j][k] = 0.0; m.e_p.x[i][j][k] = 0.0;
                    }

                    // Precipitation from convection (skip i=0: flux arrives unchanged)
                    if(i > 0)
                        m.P_conv.x[i][j][k] = std::max(0.0, m.P_conv.x[i+1][j][k]
                            + step[i] * m.r_humid.x[i][j][k]
                            * (m.g_p.x[i][j][k] - m.e_d.x[i][j][k] - m.e_p.x[i][j][k]));

                } // end i

                m.P_conv.x[0][j][k] = m.P_conv.x[1][j][k];
            }
        }
    }
/*
*
*/
// ==================== CAPE ====================
    // computeCAPE() DELETED 2026-08-19, ported from ATHAD. It was a SECOND, wrong CAPE that
    // nothing read: m.CAPE was written here and never read anywhere in the tree, and the
    // calculation was wrong three ways — it divided a PHYSICAL thickness step[i] by exp_rm
    // (the core's quadratic-stretch Jacobian, not a metric this module uses and not the right
    // one for the grid anyway); its "parcel" was the environment temperature plus a fixed
    // t_add_u = 0.2 K at every level, so it never lifted anything and never consulted s_u;
    // and m.CAPE was a 1-D array indexed by LEVEL written inside a (j,k) loop, so the last
    // column overwrote every other one. cape_col[j][k] in findCloudBaseLFS is the real CAPE —
    // a theta_e-conserving ascent with local cp and true thicknesses — and it is what seeds M_u.

// ==================== RHS FORCING ====================
    void rhsForcing() {
        // Physical caps that break the convective feedback loop
        //   v,w → E_u (∝ v·∇c) → M_u (recurrence) → MC_v,MC_w (∝ M_u·(v_u−v))
        //        → coeff_MC_vel·MC in rhs_v/rhs_w → v,w …
        // which runs away unbounded when moist physics first activates over steep
        // orography (observed 49–51°S / 72–76°W, the Patagonian Andes, ~iter 309).
        // M_u/M_d are rebuilt every MoistConvection call, so bounding them here each
        // call — together with the tendencies that feed the prognostic fields — keeps
        // the loop bounded across iterations.
        //
        // Healthy M_u/M_d are ~0.3 kg/m²s and the tendencies far below these caps, so
        // realistic convection is untouched; the caps are ~3–10× above any physical
        // value and only clip the runaway.  Tune down if convection still over-drives.
        // Bit-level non-finite reset (NaN/Inf → 0): under -ffast-math std::min/max are
        // unreliable on NaN, so use the same trick as BC_Atm.h Pass 5 / Paraview safe_val.
        // Tightened 2026-05-31 (was 10/0.05/5e-4/0.05): the Gulf-of-Alaska 60°N/151°W
        // near-surface cell pegged MC_t/MC_w at the previous caps every step (saturating
        // u/v/w at ±100 m/s), so the old caps were the velocity forcing in disguise. New
        // values are 3–5× physical, so realistic convection is still untouched.
        constexpr double M_max   = 3.0;        // [kg/(m²s)]  up/downdraft mass flux  (was 10.0; healthy ~0.3)
        constexpr double MCt_max = 0.01;       // [K/s]       convective heating      (was 0.05; ~36 K/hr, still 3× realistic)
        // 2026-06-25: 2.0e-4 (=0.2 g/kg/s) let the convective moisture pump flood the
        // near-surface layer over high tropical orography (Ethiopian highlands 9°N/37°E,
        // cloud base ≈ surface at ~2500 m) — SaturationAdjustment then condensed it
        // explosively each moist step (+8→+32 K latent heat, ratcheting against the
        // Held-Suarez recool into a condense↔recool limit cycle, surface T→48 °C). A
        // physical convective environmental-moistening rate is ~0.02 g/kg/s, so 10× down
        // throttles the SatAdjust fuel while leaving realistic convection (P_conv,
        // detrainment aloft) untouched.  See project_convective_precip_zero.
        constexpr double MCq_max = 2.0e-5;     // [(kg/kg)/s] convective moistening    (was 2.0e-4; flooded surface → SatAdjust runaway)
        // Healthy convective momentum transport is M_u·(v_u−v)/(step·ρ) ≈ 0.3·3/(400·1)
        // ≈ 2e-3 m/s². The old 0.5 cap (×coeff_MC_vel=6.25 → 3.125 non-dim RHS forcing)
        // was pinned at the Patagonian Andes (49°S/74°W) and drove the pressure solver to
        // p_dyn≈44000 hPa / t≈−1877°C at moist onset. 0.01 (~5× physical, 0.063 non-dim)
        // clips the runaway while leaving realistic convection untouched.
        constexpr double MCv_max = 0.01;       // [m/s²]      convective momentum transport  (was 0.05)

        auto safe_cap = [](double v, double mag) -> double {
            std::uint64_t bits;
            std::memcpy(&bits, &v, sizeof(bits));
            if ((bits & 0x7FF0000000000000ULL) == 0x7FF0000000000000ULL) return 0.0;
            return (v < -mag) ? -mag : (v > mag) ? mag : v;
        };

        // Bound the mass-flux arrays in place first, so the fluxes below — and u_u/v_u/
        // w_u and the next call which also read M_u/M_d — all see physical values.
        #pragma omp parallel for collapse(2) schedule(static)
        for(int j = 0; j < m.jm; j++){
            for(int k = 0; k < m.km; k++){
                for(int i = 0; i < m.im; i++){
                    m.M_u.x[i][j][k] = safe_cap(m.M_u.x[i][j][k], M_max);
                    m.M_d.x[i][j][k] = safe_cap(m.M_d.x[i][j][k], M_max);
                }
            }
        }

        #pragma omp parallel for collapse(2)
        for(int j = 0; j < m.jm; j++){
            for(int k = 0; k < m.km; k++){
                for(int i = 0; i < m.im - 1; i++){

                    double L_latent = (m.t.x[i][j][k] * m.t_0 >= m.t_0) ? m.lv : m.ls;                  // J/kg
                    // Use the same physical floor as c_u/e_l (0.01 kg/m³ ≈ air at 50 hPa).
                    // safe_r_humid's 1e-6 floor allowed 1/r → 1e6 amplification of flux
                    // divergence at near-vacuum altitudes if M_u/M_d leaked upward.
                    double inv_step_rh = 1.0 / (step[i] * std::max(m.r_humid.x[i][j][k], 0.01));        // m²/kg

                    double flux_s_ip1 = m.M_u.x[i+1][j][k] * (m.s_u.x[i+1][j][k] - m.s.x[i+1][j][k])    // kg/(m²s)
                                    + m.M_d.x[i+1][j][k] * (m.s_d.x[i+1][j][k] - m.s.x[i+1][j][k]);
                    double flux_s_i   = m.M_u.x[i][j][k] * (m.s_u.x[i][j][k] - m.s.x[i][j][k])          // kg/(m²s)    s non-dimensional
                                    + m.M_d.x[i][j][k] * (m.s_d.x[i][j][k] - m.s.x[i][j][k]);

                    double conv_src = m.c_u.x[i][j][k] - m.e_d.x[i][j][k]                               // (kg/kg)/s
                                  - m.e_l.x[i][j][k] - m.e_p.x[i][j][k];

                    const double cp_mc = AtmMixture::cp_of(m.c.x[i][j][k], m.co2.x[i][j][k],
                                             m.t.x[i][j][k] * m.t_0, m.m_comp.M_bg);
                    // s -> K uses t_0, which inverts s = cp_l*T/s_0 only if s_0 == cp_l*t_0.
                    // Here cp_l = 1349 and s_0 was not re-derived, so the true factor is
                    // s_0/cp_l = 203.50 K and the shipped one over-heats by 1.342x.
                    // SECOND TERM: (L/cp)*conv_src is ALREADY K/s, so the *t_0 that used to
                    // stand there made it 273x its own scale, and it is what held MC_t at its
                    // NEGATIVE cap wherever e_d was appreciable. Removed unconditionally.
                    const double s_to_K = AtomMoistConvection::s_consistent
                                        ? (m.s_0 / m.cp_l) : m.t_0;
                    m.MC_t.x[i][j][k] = safe_cap(
                        -(flux_s_ip1 - flux_s_i) * inv_step_rh * s_to_K                                 // K/s
                        + (L_latent / cp_mc) * conv_src, MCt_max);                                      // K/s



                    double flux_q_ip1 = m.M_u.x[i+1][j][k] * (m.q_v_u.x[i+1][j][k] - m.c.x[i+1][j][k])  // kg/(m²s)(kg/kg)
                                    + m.M_d.x[i+1][j][k] * (m.q_v_d.x[i+1][j][k] - m.c.x[i+1][j][k]);
                    double flux_q_i   = m.M_u.x[i][j][k] * (m.q_v_u.x[i][j][k] - m.c.x[i][j][k])
                                    + m.M_d.x[i][j][k] * (m.q_v_d.x[i][j][k] - m.c.x[i][j][k]);

                    m.MC_q.x[i][j][k] = safe_cap(
                                        - (flux_q_ip1 - flux_q_i) * inv_step_rh                         // (kg/kg)/s
                                        - conv_src, MCq_max);                                           // (kg/kg)/s




                    double flux_v_ip1 = m.M_u.x[i+1][j][k] * (m.v_u.x[i+1][j][k] - m.v.x[i+1][j][k])    // kg/(m²s)
                                    + m.M_d.x[i+1][j][k] * (m.v_d.x[i+1][j][k] - m.v.x[i+1][j][k]);
                    double flux_v_i   = m.M_u.x[i][j][k] * (m.v_u.x[i][j][k] - m.v.x[i][j][k])
                                    + m.M_d.x[i][j][k] * (m.v_d.x[i][j][k] - m.v.x[i][j][k]);

                    m.MC_v.x[i][j][k] = safe_cap(
                        -(flux_v_ip1 - flux_v_i) * inv_step_rh * m.u_0, MCv_max);                       // (m/s)(1/s)



                    double flux_w_ip1 = m.M_u.x[i+1][j][k] * (m.w_u.x[i+1][j][k] - m.w.x[i+1][j][k])
                                    + m.M_d.x[i+1][j][k] * (m.w_d.x[i+1][j][k] - m.w.x[i+1][j][k]);
                    double flux_w_i   = m.M_u.x[i][j][k] * (m.w_u.x[i][j][k] - m.w.x[i][j][k])
                                    + m.M_d.x[i][j][k] * (m.w_d.x[i][j][k] - m.w.x[i][j][k]);

                    m.MC_w.x[i][j][k] = safe_cap(
                        -(flux_w_ip1 - flux_w_i) * inv_step_rh * m.u_0, MCv_max);                       // (m/s)(1/s)
                }
            }
        }
    }
/*
*
*/
// ==================== SIMPSON INTEGRATION OF -c * E_u OVER i ====================
    // Accumulates ∫ -c(i,j,k)·E_u(i,j,k) dz from i_beg upward in a single (j,k) column.
    // Uses composite non-uniform Simpson's 1/3 rule on the variable vertical steps step[i].
    // Falls back to the trapezoidal rule for one remaining interval when the count is odd.
    // Returns the first i-index where the cumulative sum becomes positive, or -1 if never.
    int simpsonCxEu(int j, int k, int i_beg, int i_end) const {
        if(i_end <= i_beg) return -1;
        double cumsum = 0.0;
        int i = i_beg;
        // --- composite non-uniform Simpson pairs ---
        for(; i + 2 <= i_end; i += 2){
            double h1 = step[i];
            double h2 = step[i+1];
            double f0 = -(m.c.x[i  ][j][k] * m.E_u.x[i  ][j][k]);
            double f1 = -(m.c.x[i+1][j][k] * m.E_u.x[i+1][j][k]);
            double f2 = -(m.c.x[i+2][j][k] * m.E_u.x[i+2][j][k]);
            // non-uniform 3-point quadratic rule (reduces to h/3*(f0+4f1+f2) when h1==h2)
            double pair_contrib = (h1 + h2) / 6.0
                * ((2.0 - h2 / h1) * f0
                +  (h1 + h2) * (h1 + h2) / (h1 * h2) * f1
                +  (2.0 - h1 / h2) * f2);
            if(cumsum + pair_contrib > 0.0){
                // Refine: check trapezoidal partial [i, i+1] to find which sub-interval crosses
                return (cumsum + 0.5 * h1 * (f0 + f1) > 0.0) ? i + 1 : i + 2;
            }
            cumsum += pair_contrib;
        }
        // --- one remaining interval: trapezoidal ---
        if(i + 1 <= i_end){
            double f0 = -(m.c.x[i  ][j][k] * m.E_u.x[i  ][j][k]);
            double f1 = -(m.c.x[i+1][j][k] * m.E_u.x[i+1][j][k]);
            if(cumsum + 0.5 * step[i] * (f0 + f1) > 0.0) return i + 1;
        }
        return -1;
    }
/*
*
*/
    // Compute i_Base_local[j][k] = first i where ∫ -c·E_u dz becomes positive, per column.
    // Integrates from i_deep_beg_local[j][k] to i_LFS_local[j][k].
    // Rejects high-altitude results: cloud base below f_stat_base * p_surf (0.79).
    void CloudBaseUpdate() {
        using namespace AtomMoistConvection;
        #pragma omp parallel for collapse(2) schedule(static)
        for(int j = 1; j < m.jm-1; j++){
            for(int k = 1; k < m.km-1; k++){
                const int i_beg = i_deep_beg_local[j][k];
                const int i_end = i_LFS_local[j][k];
                if(i_beg < 0 || i_end <= i_beg) continue;
                const int i_simpson = simpsonCxEu(j, k, i_beg, i_end);
                if(i_simpson >= 0 && i_simpson < m.im
                        && m.p_stat.x[i_simpson][j][k] <= f_stat_Cloud_Base * m.p_stat.x[0][j][k]   // reject high-altitude outliers
                        // Cause #8: old condition (i_simpson > i_Base_local) only raised
                        // i_base, widening the stale-initUpdraft region below cloud base
                        // on every call.  Take the lowest (surface-nearest) valid base instead.
                        && (i_Base_local[j][k] <= 0 || i_simpson <= i_Base_local[j][k])){
                    i_Base_local[j][k] = i_simpson;
                    m.CloudBase.x[i_simpson][j][k] = height_table[i_simpson];
                }
            }
        }
    }
/*
*
*/
// ==================== SUB-TERRAIN FILL ====================
    void subTerrainFill() {
        #pragma omp parallel for collapse(2)
        for(int j = 0; j < m.jm; j++){
            for(int k = 0; k < m.km; k++){
              if(!land_surf[j * m.km + k]) continue;

                int i_mount = m.i_topography[j][k];

                // Cache mount-level values once
                double pc_m = m.P_conv.x[i_mount][j][k];
                double eu_m = m.E_u.x[i_mount][j][k], ed_m = m.E_d.x[i_mount][j][k];
                double du_m = m.D_u.x[i_mount][j][k], dd_m = m.D_d.x[i_mount][j][k];
                double mu_m = m.M_u.x[i_mount][j][k], md_m = m.M_d.x[i_mount][j][k];
                double mct_m = m.MC_t.x[i_mount][j][k], mcq_m = m.MC_q.x[i_mount][j][k];
                double mcv_m = m.MC_v.x[i_mount][j][k], mcw_m = m.MC_w.x[i_mount][j][k];
                double qvu_m = m.q_v_u.x[i_mount][j][k], qvd_m = m.q_v_d.x[i_mount][j][k];
                double qcu_m = m.q_c_u.x[i_mount][j][k];
                double cu_m = m.c_u.x[i_mount][j][k], gp_m = m.g_p.x[i_mount][j][k];
                double el_m = m.e_l.x[i_mount][j][k], ep_m = m.e_p.x[i_mount][j][k];
                double edm = m.e_d.x[i_mount][j][k];
                double uu_m = m.u_u.x[i_mount][j][k], ud_m = m.u_d.x[i_mount][j][k];
                double vu_m = m.v_u.x[i_mount][j][k], vd_m = m.v_d.x[i_mount][j][k];
                double wu_m = m.w_u.x[i_mount][j][k], wd_m = m.w_d.x[i_mount][j][k];
                double s_m = m.s.x[i_mount][j][k], su_m = m.s_u.x[i_mount][j][k];
                double sd_m = m.s_d.x[i_mount][j][k];

                for(int i = i_mount - 1; i >= 0; i--){
                    m.P_conv.x[i][j][k] = pc_m;
                    m.E_u.x[i][j][k] = eu_m; m.E_d.x[i][j][k] = ed_m;
                    m.D_u.x[i][j][k] = du_m; m.D_d.x[i][j][k] = dd_m;
                    m.M_u.x[i][j][k] = mu_m; m.M_d.x[i][j][k] = md_m;
                    m.MC_t.x[i][j][k] = mct_m; m.MC_q.x[i][j][k] = mcq_m;
                    m.MC_v.x[i][j][k] = mcv_m; m.MC_w.x[i][j][k] = mcw_m;
                    m.q_v_u.x[i][j][k] = qvu_m; m.q_v_d.x[i][j][k] = qvd_m;
                    m.q_c_u.x[i][j][k] = qcu_m;
                    m.c_u.x[i][j][k] = cu_m; m.g_p.x[i][j][k] = gp_m;
                    m.e_l.x[i][j][k] = el_m; m.e_p.x[i][j][k] = ep_m;
                    m.e_d.x[i][j][k] = edm;
                    m.u_u.x[i][j][k] = uu_m; m.u_d.x[i][j][k] = ud_m;
                    m.v_u.x[i][j][k] = vu_m; m.v_d.x[i][j][k] = vd_m;
                    m.w_u.x[i][j][k] = wu_m; m.w_d.x[i][j][k] = wd_m;
                    m.s.x[i][j][k] = s_m; m.s_u.x[i][j][k] = su_m;
                    m.s_d.x[i][j][k] = sd_m;
                }
            }
        }
    }
};
/*
*
*/

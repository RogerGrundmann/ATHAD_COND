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

// ----------------------------------------------------------------------------------
// THE PRECIPITATION CONVERGENCE PROBE WAS A BARE LEVEL INDEX.
//
// The iteration below sweeps the column top-down up to iter_prec_end (20) times and decides
// it has converged when the rain flux at ONE level stops changing:
//
//     Rain_check  = P_rain.x[23]
//     ... sweep ...
//     if(|P_rain.x[23] - Rain_check| * conv_mmd <= 1e-3) break;
//
// 23 is an index, not a height, so where it lands depends entirely on the grid:
//
//     ATHAD        im=41, L_atm 15719  ->  72.5 km
//     COND/PERID   im=61, L_atm  6287  ->  13.6 km
//     COND/PERID   im=41, L_atm  6287  ->  29.0 km
//
// COND's convective column tops at 20-28 km, so at im=61 the probe is inside the rain-bearing
// layer and the test measures something real; at im=41 it sits ABOVE it, where P_rain is ~0
// always. The difference is then ~0 on the FIRST sweep -- not because the column converged but
// because nothing was ever there to change -- so the 20-sweep iteration silently collapses to
// one and the precipitation column is used unconverged. ATHAD already runs it at 72.5 km, i.e.
// degenerate today, though only past moist_phys_start_iter = 300.
//
// NEITHER VALUE IS RIGHT. 13.6 km is not a physically motivated probe height either; it is an
// accident of im=61 that happens to land where the field varies. The correct probe is the
// SURFACE flux P_rain.x[0]: it is the accumulated result of the whole downward sweep, so it
// reflects convergence of the entire column rather than of one layer; it is what the scheme
// already normalises by (P_rain_0 below); it is what the model reports; and index 0 is the
// surface on every grid, so it is grid-independent by construction and drops out of the im and
// zeta decisions entirely.
//
// ATM_PRECIP_PROBE=1 selects the surface probe. DEFAULT OFF (legacy bare 23) so the four-arm
// im x probe measurement can separate the resolution effect from the probe effect instead of
// letting one hide inside the other. Off-branch is bit-identical.
//
// NOTE: the same bare 23 is still present in ZeroCatIceScheme.h (2 sites) and
// OneCatIceScheme.h (4 sites). Neither is the configured scheme, so both are latent rather
// than live -- checked in the control flow, not inferred from the grep.
static inline int precipProbeLevel(const cAtmosphereModel& m)
{
    static const bool surface_probe = [](){
        const char* e = getenv("ATM_PRECIP_PROBE"); return e && atoi(e) != 0; }();
    return surface_probe ? 0 : std::min(23, m.im - 1);
}


using namespace AtomUtils;


namespace ThreeCatIce {
    // ice particle parameters
    constexpr double N_i_0 = 1.0e2;                                    // 1/m3
    constexpr double N_g_0 = 4.0e6;                                    // 1/m3
    constexpr double m_i_0 = 1.0e-12;                                  // kg
    constexpr double m_i_max = 1.0e-9;                                  // kg
    constexpr double m_s_0 = 3.0e-9;                                    // kg
    constexpr double N_r_0 = 8.0e+6;                                   // 1/m3
    constexpr double N_s_0 = 8.0e+5;                                   // 1/m3

    // fall velocity prefactors
    constexpr double v_s_0 = 130.0;
    constexpr double v_r_0 = 4.9;
    constexpr double v_g_0 = 442.0;

    // microphysical rate coefficients
    constexpr double c_c_au = 4.0e-4;                                  // 1/s COSMO
    constexpr double c_i_au = 1.0e-3;                                  // 1/s COSMO
    constexpr double c_ac = 0.24;                                       // m2/kg
    constexpr double c_rim = 18.6;                                      // m2/kg
    constexpr double c_agg = 10.3;                                      // m2/kg
    constexpr double c_i_cri = 0.24;                                    // m2
    constexpr double c_r_cri = 3.2e-5;                                  // m2
    constexpr double b_ev = 5.9;                                        // m2*s/kg
    constexpr double c_r_frz = 3.75e-2;                                 // m2/(K*kg)
    constexpr double c_i_dep = 1.3e-5;                                  // m3/(s*kg)

    // graupel-specific coefficients
    constexpr double z_csg = 0.5;
    constexpr double c_g_rim = 4.43;
    constexpr double c_g_agg = 2.46;

    // mass-size relation prefactors
    constexpr double a_s_m = 0.038;
    constexpr double a_g_m = 169.6;

    // temperature thresholds
    constexpr double t_nuc = 267.15;                                    // K  (-6 C)
    constexpr double t_d = 248.15;                                      // K  (-25 C)
    constexpr double t_hn = 236.15;                                     // K  (-37 C)
    constexpr double t_r_frz = 271.15;                                  // K  (-2 C)

    // iteration / conversion
    constexpr double conv_mmd = 8.64e4;                                 // s/d conversion to mm/d
    constexpr int iter_prec_end = 2;                                    // COSMO iterations
}


// Three-Category-Ice-Scheme, COSMO-module from the German Weather Forecast,
// yielding the precipitation distribution from rain, snow, and graupel.
class ThreeCatIceScheme {
public:
    explicit ThreeCatIceScheme(cAtmosphereModel& model)
        : m(model)
    {}

    void run() {
        using namespace std;
        using namespace ThreeCatIce;

        cout << endl << endl << endl << "      ThreeCategoryIceScheme" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        initArrays();
        computeColumns();
        applyBoundaryConditions();
        applyTopography();

        printReport();

        cout << "      ThreeCategoryIceScheme ended" << endl;

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for ThreeCategoryIceScheme\n", elapsed.count() * 1e-9);
    }

private:
    cAtmosphereModel& m;

    // diagnostic output state (set during computeColumns)
    bool rain = false;
    bool snow = false;
    bool graupel = false;
    double global_max_rain = 0.0;
    double global_max_snow = 0.0;
    double global_max_graupel = 0.0;
    int i_rain = 0, j_rain = 0, k_rain = 0;
    int i_snow = 0, j_snow = 0, k_snow = 0;
    int i_graupel = 0, j_graupel = 0, k_graupel = 0;

    // layer geometry tables (filled in initArrays, read in computeColumns/printReport)
    std::vector<double> step_table;
    std::vector<double> height_table;


    // ==================== INIT ====================
    void initArrays() {
        // Precompute layer thickness table
        step_table.assign(m.im, 0.0);
        height_table.resize(m.im);
        for(int i = 0; i < m.im; i++){
            height_table[i] = m.get_layer_height(i);
            if(i < m.im - 1)
                step_table[i] = m.get_layer_height(i + 1) - m.get_layer_height(i);
        }

        // Clamping + copy pass
        #pragma omp parallel for collapse(2)
        for(int j = 0; j < m.jm; j++){
            for(int k = 0; k < m.km; k++){
                #pragma omp simd
                for(int i = 0; i < m.im; i++){
                    m.c.x[i][j][k]     = std::max(0.0, m.c.x[i][j][k]);
                    m.cloud.x[i][j][k] = std::max(0.0, m.cloud.x[i][j][k]);
                    m.ice.x[i][j][k]   = std::max(0.0, m.ice.x[i][j][k]);
                    m.gr.x[i][j][k]    = std::max(0.0, m.gr.x[i][j][k]);
                    m.P_rain.x[i][j][k] = m.P_rainn.x[i][j][k];
                    m.P_snow.x[i][j][k] = m.P_snown.x[i][j][k];
                }
            }
        }
    }


    // ==================== COLUMN COMPUTATION ====================
    void computeColumns() {
        using namespace ThreeCatIce;

        // Precompute gamma and composite constants
        const double gamma_4_5 = tgamma(4.5);
        const double gamma_4_0 = tgamma(4.0);
        const double gamma_3_25 = tgamma(3.25);

        const double A_r = m.r_0_water * M_PI * N_r_0;
        const double A_s = 2.0 * a_s_m * N_s_0;
        const double A_g = 2.0 * a_g_m * N_g_0;
        const double B_rad = m.r_0_water * M_PI * N_r_0 * v_r_0 * gamma_4_5 / gamma_4_0;
        const double B_s = N_s_0 * v_s_0 * a_s_m * gamma_3_25;
        const double B_g = N_g_0 * v_g_0 * a_g_m * gamma_3_25;

        // Precompute constant pow() bases — called im×jm×km times otherwise
        const double B_rad_neg89 = pow(B_rad, -8.0/9.0);
        const double B_s_neg1213 = pow(B_s, -12.0/13.0);
        const double B_g_neg1213 = pow(B_g, -12.0/13.0);

        // Precompute E_sat at t_0 (constant, used in every cell)
        const double E_sat_t_0_Pa = 1e2 * SaturationH2O::saturationPressure(m.t_0);

        // Thread-local accumulators for the OMP reduction
        double local_max_rain = 0.0, local_max_snow = 0.0, local_max_graupel = 0.0;

        #pragma omp parallel for collapse(2) schedule(static) \
            reduction(max:local_max_rain, local_max_snow, local_max_graupel)
        for(int j = 1; j < m.jm - 1; j++){
            for(int k = 1; k < m.km - 1; k++){

                m.P_rain.x[precipProbeLevel(m)][j][k] = 0.0;

                for(int iter_prec = 1; iter_prec <= iter_prec_end; iter_prec++){

                    double Rain_check = m.P_rain.x[precipProbeLevel(m)][j][k];

                    m.P_rain.x[m.im-1][j][k] = 0.0;
                    m.P_snow.x[m.im-1][j][k] = 0.0;
                    m.P_graupel.x[m.im-1][j][k] = 0.0;

                    for(int i = m.im - 2; i >= 0; i--){


                        // ATHAD: no condensed phase can exist here — see
                        // IceSchemeCommon::canCondense. Everything below the cloud deck
                        // (ground to ~240 km) is supercritical or superheated, and the
                        // microphysics below has no meaning there.
                        {
                            const double t_guard = m.t.x[i][j][k] * m.t_0;
                            if (!IceSchemeCommon::canCondense(m, t_guard, i, j, k)) {
                                IceSchemeCommon::evaporateWhereImpossible(m, t_guard, i, j, k);
                                continue;
                            }
                        }

                        // Normalize precipitation by the surface flux. FLOORED denominator: a
                        // tiny-but-nonzero surface flux made Snow = P_snow[i]/P_snow_0 explode,
                        // and S_s_rim ∝ Snow then amplified P_snow geometrically down the column
                        // to overflow (the ThreeCat NaN blow-up). Flooring at P_norm_floor bounds
                        // the normalized ratios; the flux cap below guarantees finiteness.
                        constexpr double P_norm_floor = 1.0e-6;      // kg/(m2*s) ~0.09 mm/d
                        double P_rain_0    = std::max(m.P_rain.x[0][j][k],    P_norm_floor);
                        double P_snow_0    = std::max(m.P_snow.x[0][j][k],    P_norm_floor);
                        double P_graupel_0 = std::max(m.P_graupel.x[0][j][k], P_norm_floor);

                        double Rain    = m.P_rain.x[i][j][k]    / P_rain_0;
                        double Snow    = m.P_snow.x[i][j][k]    / P_snow_0;
                        double Graupel = m.P_graupel.x[i][j][k] / P_graupel_0;

                        // pow() calls guarded — skip entirely when base is zero
                        double Rain_79  = (Rain > 0.0) ? pow(Rain, 7.0/9.0)  : 0.0;
                        double Rain_89  = (Rain > 0.0) ? pow(Rain, 8.0/9.0)  : 0.0;
                        double Rain_139 = (Rain > 0.0) ? pow(Rain, 13.0/9.0) : 0.0;
                        double Rain_16  = (Rain > 0.0) ? pow(Rain, 1.0/6.0)  : 0.0;
                        double Rain_49  = (Rain > 0.0) ? pow(Rain, 4.0/9.0)  : 0.0;

                        double r_q_r = A_r * B_rad_neg89 * Rain_89;
                        double r_q_s = A_s * B_s_neg1213 * ((Snow > 0.0)    ? pow(Snow,    12.0/13.0) : 0.0);
                        double r_q_g = A_g * B_g_neg1213 * ((Graupel > 0.0) ? pow(Graupel, 12.0/13.0) : 0.0);

                        double t_u   = m.t.x[i][j][k] * m.t_0;
                        double p_u   = m.p_stat.x[i][j][k];
                        double p_u_0 = 1e2 * p_u;
                        double r_h_i = m.r_humid.x[i][j][k];
                        double c_ijk = m.c.x[i][j][k];
                        double cl_i  = m.cloud.x[i][j][k];
                        double ice_i = m.ice.x[i][j][k];

                        double step_i = step_table[i];

                        double q_sat = IceSchemeCommon::qSatWater(m, t_u, i, j, k);
                        double q_Ice = IceSchemeCommon::qSatIce(m, t_u, i, j, k);

                        double dt_rain_dim = step_i / 1.6;
                        double dt_snow_dim = step_i / 0.96;

                        // Precompute shared pow() terms for snow/graupel deposition/melting
                        double rh_rqs_08 = (r_q_s > 0.0) ? pow(r_h_i * r_q_s, 0.8)     : 0.0;
                        double rh_rqg_06 = (r_q_g > 0.0) ? pow(r_h_i * r_q_g, 0.6)     : 0.0;
                        double rh_rqg_95 = (r_q_g > 0.0) ? pow(r_h_i * r_q_g, 0.94878) : 0.0;

                        // --- Ice particle properties + vapour->ice->snow throttle (shared) ---
                        IceSchemeCommon::IceSnowRates thr =
                            IceSchemeCommon::depositionThrottle(m, i, j, k, t_u, q_Ice, dt_snow_dim);
                        double m_i = thr.m_i;
                        double N_i = thr.N_i;

                        // --- Nucleation ---
                        double S_nuc = 0.0;
                        if(ice_i == 0.0){
                            if((t_u < t_d && c_ijk >= q_Ice)
                                || (t_d <= t_u && t_u <= t_nuc && c_ijk >= q_sat))
                                S_nuc = m_i_0 / (r_h_i * dt_snow_dim) * N_i;
                        }

                        double S_c_frz = (t_u < t_hn && cl_i > 0.0)
                            ? cl_i / dt_rain_dim : 0.0;

                        // --- Deposition growth of cloud ice (shared throttle) ---
                        double S_i_dep = thr.S_i_dep;

                        // --- Autoconversion ---
                        double S_c_au = (t_u >= m.t_0 && cl_i > 0.0)
                            ? std::max(c_c_au * (cl_i - 0.0002), 0.0) : 0.0;

                        double S_i_au = thr.S_i_au;   // shared throttle (ice aggregation -> snow)
                        double S_d_au = thr.S_d_au;   // shared throttle (ice deposition -> snow)

                        // --- Collection ---
                        double S_ac = (t_u >= m.t_0) ? c_ac * cl_i * Rain_79 : 0.0;

                        double S_s_rim, S_g_rim, S_s_shed, S_g_shed;
                        if(t_u < m.t_0){
                            S_s_rim  = c_rim * cl_i * Snow;
                            S_g_rim  = c_rim * cl_i * rh_rqg_95;
                            S_s_shed = 0.0;
                            S_g_shed = 0.0;
                        }else{
                            S_s_rim  = 0.0;
                            S_g_rim  = 0.0;
                            S_s_shed = c_rim * cl_i * Snow;
                            S_g_shed = c_g_rim * cl_i * rh_rqg_95;
                        }

                        double S_s_agg, S_g_agg, S_i_cri, S_r_cri;
                        if(t_u <= m.t_0){
                            S_s_agg = c_agg * ice_i * Snow;
                            S_g_agg = c_g_agg * ice_i * rh_rqg_95;
                            S_i_cri = c_i_cri * ice_i * Rain_79;
                            S_r_cri = c_r_cri * ice_i / m_i * Rain_139;
                        }else{
                            S_s_agg = 0.0;
                            S_g_agg = 0.0;
                            S_i_cri = 0.0;
                            S_r_cri = 0.0;
                        }

                        // --- Evaporation ---
                        double S_ev = 0.0;
                        if(t_u >= m.t_0 && Rain > 0.0){
                            double a_ev = 2.76e-3 * exp(0.055 * (m.t_0 - t_u));
                            S_ev = a_ev * (1.0 + b_ev * Rain_16)
                                 * (q_sat - c_ijk) * Rain_49;
                        }

                        // --- Deposition/sublimation of snow and graupel ---
                        double S_s_dep = 0.0, S_g_dep = 0.0;
                        double inv_p_u_0 = 1.0 / p_u_0;
                        double c_minus_qIce = c_ijk - q_Ice;

                        // Compute d_v, l_h, t_crit for melting branch
                        double d_v = (t_u >= m.t_0)
                            ? 101325.0 * inv_p_u_0 * (2.22e-5 + 1.46e-7 * (t_u - m.t_0))
                            : 101325.0 * inv_p_u_0 * (2.22e-5 + 1.25e-7 * (t_u - m.t_0));
                        double l_h = 0.024 + 8.0e-5 * (t_u - m.t_0);
                        // Exact conversion; E and p are both in Pa, and only their ratio
                        // is used, so the units are consistent.
                        double q_sat_t_0 = SaturationH2O::saturationMassFraction(
                            E_sat_t_0_Pa, p_u_0,
                            AtmMixture::M_nonwater(c_ijk, m.co2.x[i][j][k], m.m_comp.M_bg));

                        double t_crit = (c_ijk >= q_sat_t_0)
                            ? m.t_0 - (1.0 / l_h) * m.ls * d_v * r_h_i * (c_ijk - q_sat_t_0)
                            : m.t_0;

                        if(t_u < m.t_0){
                            S_s_dep = (2.91955 - 0.0109928*t_u + 15871.3*inv_p_u_0 + 1.74744e-6*p_u_0)
                                    * c_minus_qIce * rh_rqs_08;
                            S_g_dep = 0.0;  // as in original
                        }else if(t_u < t_crit){
                            S_s_dep = (0.28003 - 0.146293e-6*p_u_0)
                                    * c_minus_qIce * rh_rqs_08;
                            S_g_dep = (0.0418521 - 4.7524e-8*p_u_0)
                                    * c_minus_qIce * rh_rqg_06;
                        }else{
                            S_s_dep = (2.41897 + 31282.3*inv_p_u_0)
                                    * c_minus_qIce * rh_rqs_08;
                            S_g_dep = (0.153907 - 7.86703e-7*p_u_0)
                                    * c_minus_qIce * rh_rqg_06;
                        }

                        // --- Melting ---
                        double S_s_melt = 0.0, S_g_melt = 0.0;
                        if(t_u >= m.t_0){
                            constexpr double a_melt = 2.95e+3;
                            double melt_driver = (t_u - m.t_0) + a_melt * (c_ijk - q_sat_t_0);
                            S_s_melt = melt_driver * (0.612654e-3 + 79.6863 * inv_p_u_0) * rh_rqs_08;
                            S_g_melt = melt_driver * (7.39441e-5  + 12.31698 * inv_p_u_0) * rh_rqg_06;
                        }

                        double S_i_melt = (t_u > m.t_0 && ice_i > 0.0)
                            ? ice_i / dt_snow_dim : 0.0;

                        // --- Freezing ---
                        double S_r_frz;
                        if(t_u > t_hn){
                            double t_frz = std::max(t_r_frz - t_u, 0.0);
                            S_r_frz = c_r_frz * pow(t_frz, 1.5)
                                    * pow(r_h_i * r_q_r, 27.0/16.0);
                        }else{
                            S_r_frz = Rain / dt_rain_dim;
                        }

                        // --- Snow-to-graupel conversion ---
                        double S_csg = (t_u <= m.t_0 && cl_i > 0.0002)
                            ? z_csg * cl_i * ((r_q_s > 0.0) ? pow(r_h_i * r_q_s, 0.75) : 0.0)
                            : 0.0;

                        // --- Sinks and sources ---
                        double S_c_c_ijk = m.S_c_c.x[i][j][k];

                        m.S_v.x[i][j][k] = -S_c_c_ijk + S_ev - S_i_dep - S_s_dep - S_g_dep - S_nuc;
                        m.S_c.x[i][j][k] =  S_c_c_ijk - S_c_au - S_ac - S_c_frz + S_i_melt
                                           - S_s_rim - S_g_rim - S_s_shed - S_g_shed;
                        m.S_i.x[i][j][k] =  S_nuc + S_c_frz + S_i_dep - S_i_melt
                                           - S_i_au - S_d_au - S_s_agg - S_g_agg - S_i_cri;
                        m.S_r.x[i][j][k] =  S_c_au + S_ac - S_ev + S_s_shed + S_g_shed
                                           - S_r_cri - S_r_frz + S_s_melt + S_g_melt;
                        m.S_s.x[i][j][k] =  S_i_au + S_d_au + S_s_agg + S_s_rim + S_s_dep
                                           + S_i_cri + S_r_cri - S_s_melt - S_csg;
                        m.S_g.x[i][j][k] =  S_g_agg + S_g_rim + S_g_dep + S_i_cri + S_r_cri
                                           + S_r_frz - S_g_melt + S_csg;

                        // --- Column integration (top-down) ---
                        // Each flux is clamped to [0, P_max_flux]. The upper cap is a hard backstop
                        // against the P_snow riming runaway (see the floored normalization above):
                        // even if a source spikes, the accumulated flux can no longer overflow to
                        // inf, so the scheme stays finite (over-precipitates at worst, like the
                        // other untuned schemes — usable, not NaN). ~260 mm/d, well above any
                        // physical precip.
                        constexpr double P_max_flux = 3.0e-3;        // kg/(m2*s) ~260 mm/d hard cap
                        m.P_rain.x[i][j][k] = (t_u >= m.t_0)
                            ? std::min(P_max_flux, std::max(0.0, m.P_rain.x[i+1][j][k]
                                + m.r_humid.x[i+1][j][k] * m.S_r.x[i+1][j][k] * step_i))
                            : 0.0;

                        m.P_snow.x[i][j][k] = (t_u < m.t_0 && t_u >= m.t_000)
                            ? std::min(P_max_flux, std::max(0.0, m.P_snow.x[i+1][j][k]
                                + m.r_humid.x[i+1][j][k] * m.S_s.x[i+1][j][k] * step_i))
                            : 0.0;

                        m.P_graupel.x[i][j][k] = (t_u < m.t_0 && t_u >= m.t_00)
                            ? std::min(P_max_flux, std::max(0.0, m.P_graupel.x[i+1][j][k]
                                + m.r_humid.x[i+1][j][k] * m.S_g.x[i+1][j][k] * step_i))
                            : 0.0;

                        // Track column maxima (thread-local via reduction)
                        local_max_rain    = std::max(local_max_rain,    m.P_rain.x[i][j][k]);
                        local_max_snow    = std::max(local_max_snow,    m.P_snow.x[i][j][k]);
                        local_max_graupel = std::max(local_max_graupel, m.P_graupel.x[i][j][k]);

                        m.Precipitation.x[i][j][k] = m.P_rain.x[i][j][k]
                            + m.P_snow.x[i][j][k] + m.P_graupel.x[i][j][k];

                    } // end i

                    double P_rain_diff = fabs(m.P_rain.x[precipProbeLevel(m)][j][k] - Rain_check) * conv_mmd;
                    if(P_rain_diff <= 1.0e-3) break;

                } // end iter_prec
            } // end k
        } // end j

        global_max_rain    = local_max_rain;
        global_max_snow    = local_max_snow;
        global_max_graupel = local_max_graupel;

        // Locate global max positions (single cheap pass, outside hot loop)
        if(global_max_rain > 0.0 || global_max_snow > 0.0 || global_max_graupel > 0.0){
            #pragma omp parallel for collapse(2)
            for(int j = 1; j < m.jm - 1; j++){
                for(int k = 1; k < m.km - 1; k++){
                    for(int i = 0; i < m.im; i++){
                        if(m.P_rain.x[i][j][k] == global_max_rain && global_max_rain > 0.0){
                            #pragma omp critical(rain_loc)
                            { i_rain = i; j_rain = j; k_rain = k; rain = true; }
                        }
                        if(m.P_snow.x[i][j][k] == global_max_snow && global_max_snow > 0.0){
                            #pragma omp critical(snow_loc)
                            { i_snow = i; j_snow = j; k_snow = k; snow = true; }
                        }
                        if(m.P_graupel.x[i][j][k] == global_max_graupel && global_max_graupel > 0.0){
                            #pragma omp critical(graupel_loc)
                            { i_graupel = i; j_graupel = j; k_graupel = k; graupel = true; }
                        }
                    }
                }
            }
        }
    }


    // ==================== BOUNDARY CONDITIONS ====================
    void applyBoundaryConditions() {
        IceSchemeCommon::extrapolateBC(m, m.P_rain,    true);   // rain:    phi-seam periodicity averaged
        IceSchemeCommon::extrapolateBC(m, m.P_snow,    true);   // snow:    phi-seam periodicity averaged
        IceSchemeCommon::extrapolateBC(m, m.P_graupel, true);   // graupel: phi-seam periodicity averaged
    }


    // ==================== TOPOGRAPHY FILL ====================
    void applyTopography() {
        IceSchemeCommon::fillTopography(m, m.P_rain);
        IceSchemeCommon::fillTopography(m, m.P_snow);
        IceSchemeCommon::fillTopography(m, m.P_graupel);
    }


    // ==================== DIAGNOSTIC REPORT ====================
    void printReport() const {
        using namespace std;
        using namespace ThreeCatIce;

        if(!rain)
            cout << endl << "      no rain fall in ThreeCategoryIceScheme found" << endl << endl;
        else
            cout << endl << "      rain fall found  i=" << i_rain << " j=" << j_rain << " k=" << k_rain
                 << " h=" << height_table[i_rain] << "m  P=" << global_max_rain * conv_mmd << "mm/d" << endl;

        if(!snow)
            cout << "      no snow fall found" << endl;
        else
            cout << "      snow fall found  i=" << i_snow << " j=" << j_snow << " k=" << k_snow
                 << " h=" << height_table[i_snow] << "m  P=" << global_max_snow * conv_mmd << "mm/d" << endl;

        if(!graupel)
            cout << "      no graupel fall found" << endl;
        else
            cout << "      graupel fall found  i=" << i_graupel << " j=" << j_graupel << " k=" << k_graupel
                 << " h=" << height_table[i_graupel] << "m  P=" << global_max_graupel * conv_mmd << "mm/d" << endl;
    }
};

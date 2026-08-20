#pragma once

#include "cAtmosphereModel.h"
#include "Utils.h"

#include <iostream>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AtomUtils;

class UtilsAtm {
public:
    explicit UtilsAtm(cAtmosphereModel& model)
        : m(model)
    {}

    // ------------------------------------------------------------------
    void storeIntermediateData3D(float coeff = 1.0f)
    {
        const double cf = coeff;
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                #pragma omp simd
                for (int k = 0; k < m.km; k++) {
                    m.un.x[i][j][k]     = cf * m.u.x[i][j][k];
                    m.vn.x[i][j][k]     = cf * m.v.x[i][j][k];
                    m.wn.x[i][j][k]     = cf * m.w.x[i][j][k];
                    m.tn.x[i][j][k]     = cf * m.t.x[i][j][k];
                    m.cn.x[i][j][k]     = cf * m.c.x[i][j][k];
                    m.cloudn.x[i][j][k] = cf * m.cloud.x[i][j][k];
                    m.icen.x[i][j][k]   = cf * m.ice.x[i][j][k];
                    m.grn.x[i][j][k]    = cf * m.gr.x[i][j][k];
                    m.P_rainn.x[i][j][k]= cf * m.P_rain.x[i][j][k];
                    m.P_snown.x[i][j][k]= cf * m.P_snow.x[i][j][k];
                    m.co2n.x[i][j][k]   = cf * m.co2.x[i][j][k];
                    m.tken.x[i][j][k]   = cf * m.tke.x[i][j][k];
                    m.disn.x[i][j][k]   = cf * m.dis.x[i][j][k];
                }
            }
        }
    }

    // ------------------------------------------------------------------
    void precipitationSum()
    {
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                #pragma omp simd
                for (int k = 0; k < m.km; k++) {
                    m.Precipitation.x[i][j][k] =
                          m.P_rain.x[i][j][k]
                        + m.P_snow.x[i][j][k]
                        + m.P_graupel.x[i][j][k]
                        + m.P_conv.x[i][j][k];
                }
            }
        }
    }

    // ------------------------------------------------------------------
    void findResiduumAtm()
    {
        using namespace std;
        cout << endl << "      AGCM: find_residuum_atm" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        // The residuum of the PREVIOUS call, captured before anything can overwrite it.
        // It is what the "declining" test and the reported errors below are supposed to
        // compare against; see the race note at the maximum test.
        const double residuum_prev = m.residuum_old;

        struct ErrorState { double val; int i, j, k; };
        ErrorState global_max = {0.0, 0, 0, 0};
        double global_sum = 0.0;
        long   global_cnt = 0;
        ErrorState global_first_nan = {0.0, -1, -1, -1};
        long   global_nan_cnt = 0;

        #pragma omp parallel
        {
            ErrorState local_max = {0.0, 0, 0, 0};
            double local_sum = 0.0;
            long   local_cnt = 0;
            ErrorState local_first_nan = {0.0, -1, -1, -1};
            long   local_nan_cnt = 0;

            #pragma omp for collapse(2) schedule(static)
//            for (int j = 45; j <= 140; j++) {
//                for (int k = 150; k <= 240; k++) {
            for (int j = 1; j < m.jm - 1; j++) {
                for (int k = 1; k < m.km - 1; k++) {

                    double sinthe = std::max(0.55, sin(m.the.z[j]));   // metric floor ~57°, in sync with RHS/pressure-solver (residuum diagnostic only)

                    for (int i = 10; i <= 30; i++) {
                        double rm       = m.rad.z[i];
                        double exp_rm   = 1.0 / (rm + 1.0);
                        double rmsinthe = rm * sinthe;

                        double dudr   = (m.u.x[i+1][j][k] - m.u.x[i-1][j][k])
                                        / (2.0 * m.dr) * exp_rm;
                        double dvdthe = (m.v.x[i][j+1][k] - m.v.x[i][j-1][k])
                                        / (2.0 * rm * m.dthe);
                        double dwdphi = (m.w.x[i][j][k+1] - m.w.x[i][j][k-1])
                                        / (2.0 * rmsinthe * m.dphi);

                        double res = sqrt((dudr * dudr
                                        + dvdthe * dvdthe
                                        + dwdphi * dwdphi) / 3.0);

                        // Keep non-finite residuals out of the running sum so the
                        // grid-average error stays meaningful; record the first cell.
                        if(!is_finite_safe(res)){
                            if(local_nan_cnt == 0) local_first_nan = {res, i, j, k};
                            local_nan_cnt++;
                            continue;
                        }

                        local_sum += res;
                        local_cnt++;

                        // `m.residuum_old = res;` used to sit here. It is a SHARED model
                        // member written from inside every thread's loop with no
                        // synchronisation — a data race. Each thread wrote its own running
                        // maximum and whichever wrote last survived, so the value was not
                        // the previous iteration's residuum at all. It is read three times
                        // below — the "declining" vs "too high" message and both reported
                        // errors — so a raced value drove the line a human reads to decide
                        // whether the run is converging. Now set once, after the reduction,
                        // from the previous call's value. Ported from ATHAD (its item 18).
                        if(res > local_max.val) local_max = {res, i, j, k};
                   }
                }
            }

            #pragma omp critical
            {
                // Deterministic tie-break. A plain `>` lets the FIRST thread into the
                // critical section win an equal maximum, so the reported error LOCATION
                // depended on thread arrival order even though its value did not. Ties are
                // not rare here: the model is hemispherically symmetric by construction, so
                // equal maxima are the expected case, not an accident. Preferring the
                // lexicographically smallest (i,j,k) makes the choice a property of the
                // field rather than of the schedule.
                const bool better = (local_max.val > global_max.val)
                                 || (local_max.val == global_max.val
                                     && (local_max.i < global_max.i
                                         || (local_max.i == global_max.i
                                             && (local_max.j < global_max.j
                                                 || (local_max.j == global_max.j
                                                     && local_max.k < global_max.k)))));
                if (better) global_max = local_max;
                global_sum += local_sum;
                global_cnt += local_cnt;
                if (local_nan_cnt > 0){
                    if (global_nan_cnt == 0) global_first_nan = local_first_nan;
                    global_nan_cnt += local_nan_cnt;
                }
            }
        }

        const double avg_abs = (global_cnt > 0) ? global_sum / global_cnt : 0.0;
        const double avg_rel = (global_max.val > 0.0) ? avg_abs / global_max.val : 0.0;

        cout.precision(8);
        const bool declining = (residuum_prev - global_max.val) > 0.0;
        cout << endl
             << (declining
                 ? "      AGCM: find_residuum_atm, absolute error declining .......................\n"
                 : "      AGCM: find_residuum_atm, absolute error is too high .....................\n")
             << "      residuum_atm = " << global_max.val
             << "      residuum_old = " << residuum_prev
             << "      eps_residuum = " << m.eps_residuum << endl << endl
             << "      i_error = " << global_max.i
             << "   j_error = "    << global_max.j
             << "   k_error = "    << global_max.k << endl << endl
             << "      absolute error = " << fabs(residuum_prev - global_max.val) << endl
             << "      relative error = " << fabs(global_max.val / residuum_prev - 1.0) << endl << endl
             << "      error location: lat = " << (90 - global_max.j) << " deg N"
             << "   lon = " << global_max.k << " deg E"
             << "   height = " << global_max.i * 400 << " m" << endl
             << "      i_topography at error = " << m.i_topography[global_max.j][global_max.k]
             << "   surface height = " << m.get_layer_height(m.i_topography[global_max.j][global_max.k]) << " m" << endl
             << "      h.x[i][j][k] = " << m.h.x[global_max.i][global_max.j][global_max.k]
             << "   (0=air, 1=land)" << endl << endl
             << "      grid average absolute error = " << avg_abs
             << "   (" << global_cnt << " cells)" << endl
             << "      grid average relative error = " << avg_rel
             << "   (avg/max)" << endl << endl;

        if(global_nan_cnt > 0){
            cout << "      AGCM: find_residuum_atm WARNING — " << global_nan_cnt
                 << " non-finite residual cell(s) in scan band i=10..30; first at"
                 << " i=" << global_first_nan.i
                 << "  j=" << global_first_nan.j << " (lat " << (90 - global_first_nan.j) << " N)"
                 << "  k=" << global_first_nan.k << " (lon " << global_first_nan.k << " E)"
                 << "  height=" << global_first_nan.i * 400 << " m"
                 << "  (excluded from the average)"
                 << endl << endl;
        }

        // Set ONCE, after the reduction, so the next call compares against this call's
        // converged maximum rather than against whatever a thread wrote last.
        m.residuum_old = global_max.val;

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for find_residuum_atm\n", elapsed.count() * 1e-9);
        cout << "      AGCM: find_residuum_atm ended" << endl;
    }

    // ------------------------------------------------------------------
    void valueLimitationAtm()
    {
        using namespace std;
        cout << endl << endl << endl << "      ValueLimitationAtm" << endl;

        for (int k = 0; k < m.km; k++) {
            for (int j = 0; j < m.jm; j++) {
                for (int i = 0; i < m.im; i++) {
/*
                    if (m.u.x[i][j][k] >=  0.5)   m.u.x[i][j][k] =  0.5;
                    if (m.u.x[i][j][k] <= -0.5)   m.u.x[i][j][k] = -0.5;
                    if (m.v.x[i][j][k] >=  0.5)   m.v.x[i][j][k] =  0.5;
                    if (m.v.x[i][j][k] <= -0.5)   m.v.x[i][j][k] = -0.5;
                    if (m.w.x[i][j][k] >=  6.25)  m.w.x[i][j][k] =  6.25;
                    if (m.w.x[i][j][k] <= -2.0)   m.w.x[i][j][k] = -2.0;
*/
                    if (m.p_stat.x[i][j][k] <= 0.0)  m.p_stat.x[i][j][k] = 0.0;

                    if (m.r_dry.x[i][j][k]   <= 0.0)  m.r_dry.x[i][j][k]   = 1e-6;
                    if (m.r_humid.x[i][j][k] <= 0.0)  m.r_humid.x[i][j][k] = 1e-6;

                    // A mass fraction lives in [0,1]. The upper bound is not defensive
                    // padding: without it the Earth surface-evaporation scheme, driven by a
                    // Magnus saturation pressure evaluated far outside its validity range,
                    // wrote c = 20.8 — twenty times the total mass present — and every
                    // mixture property derived from it followed.
                    //
                    // ATHAD: the ceiling is 1 - co2, not 1. The three mass fractions sum to
                    // one and the background is carried as the remainder 1 - c - co2, so
                    // c > 1 - co2 is a NEGATIVE background mass. Clamping at 1 let the
                    // equatorial column sit at exactly c = 1.0000 through the sub-cloud band,
                    // a composition summing to 1.21; AtmMixture::split() renormalises
                    // defensively, so the only visible symptom was a gas constant pinned at
                    // 415.1 instead of moving with the composition. On Earth water is ~1 % of
                    // the mass and CO2 is ppm, so 1 and 1 - co2 are the same number.
                    if (m.c.x[i][j][k] < 0.0)  m.c.x[i][j][k] = 0.0;
                    {
                        const double c_ceiling = std::max(0.0, 1.0 - m.co2.x[i][j][k]);
                        if (m.c.x[i][j][k] > c_ceiling)  m.c.x[i][j][k] = c_ceiling;
                    }
                    if (m.co2.x[i][j][k] < 0.0)  m.co2.x[i][j][k] = 0.0;
                    if (m.co2.x[i][j][k] > 1.0)  m.co2.x[i][j][k] = 1.0;
                    if (m.cloud.x[i][j][k] >= 0.02)  m.cloud.x[i][j][k] = 0.02;
                    if (m.cloud.x[i][j][k] <  0.0)   m.cloud.x[i][j][k] = 0.0;
                    if (m.ice.x[i][j][k]   >= 0.01)  m.ice.x[i][j][k]   = 0.01;
                    if (m.ice.x[i][j][k]   <  0.0)   m.ice.x[i][j][k]   = 0.0;

                    double t_u = m.t.x[i][j][k] * m.t_0;
                    if (t_u <= m.t_00) {
                        m.c.x[i][j][k]     = 0.0;
                        m.cloud.x[i][j][k] = 0.0;
                        m.ice.x[i][j][k]   = 0.0;
                        m.gr.x[i][j][k]    = 0.0;

                        m.S_c_c.x[i][j][k] = 0.0;
                        m.S_c.x[i][j][k]   = 0.0;
                        m.S_v.x[i][j][k]   = 0.0;
                        m.S_r.x[i][j][k]   = 0.0;
                        m.S_s.x[i][j][k]   = 0.0;
                    }

                    if (m.P_rain.x[i][j][k]        < 0.0)  m.P_rain.x[i][j][k]        = 0.0;
                    if (m.P_snow.x[i][j][k]        < 0.0)  m.P_snow.x[i][j][k]        = 0.0;
                    if (m.P_graupel.x[i][j][k]     < 0.0)  m.P_graupel.x[i][j][k]     = 0.0;
                    if (m.Precipitation.x[i][j][k] < 0.0)  m.Precipitation.x[i][j][k] = 0.0;
                    if (m.P_conv.x[i][j][k]   < 0.0)  m.P_conv.x[i][j][k]   = 0.0;

                    if (m.Q_Latent.x[i][j][k]   <= 0.0)  m.Q_Latent.x[i][j][k]        = 0.0;
                    if (m.Q_Sensible.x[i][j][k] <= 0.0)  m.Q_Sensible.x[i][j][k]      = 0.0;

                }
            }
        }

        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                int i_mount = m.i_topography[j][k];

                if (is_land(m.h, i_mount, j, k)) {

                    m.t.x[0][j][k]     = m.t.x[i_mount][j][k];
                    m.c.x[0][j][k]     = m.c.x[i_mount][j][k];
                    m.cloud.x[0][j][k] = m.t.x[i_mount][j][k];
                    m.ice.x[0][j][k]   = m.t.x[i_mount][j][k];

                    m.S_c_c.x[0][j][k] = m.S_c_c.x[i_mount][j][k];
                    m.S_v.x[0][j][k]   = m.S_v.x[i_mount][j][k];
                    m.S_c.x[0][j][k]   = m.S_c.x[i_mount][j][k];
                    m.S_i.x[0][j][k]   = m.S_i.x[i_mount][j][k];
                    m.S_r.x[0][j][k]   = m.S_r.x[i_mount][j][k];
                    m.S_s.x[0][j][k]   = m.S_s.x[i_mount][j][k];

                    m.P_rain.x[0][j][k] = m.P_rain.x[i_mount][j][k];
                    m.P_snow.x[0][j][k] = m.P_snow.x[i_mount][j][k];
                }
            }
        }

        cout << "      ValueLimitationAtm ended" << endl;
    }

    // ------------------------------------------------------------------
    void writeFile(std::string& bathymetry_name, std::string& output_path, bool is_final_result)
    {
        using namespace std;
        cout << endl << endl << endl << "      AGCM: write_file" << endl;

        // Skip VTK slices on the final-result call when iter_n==0 (Ma==0 case):
        // both writeFile calls would write to the same filename; the second
        // truncates the file before writing, so a failure there leaves 0 bytes.
        // For Ma>0 the 3D loop advances iter_n, so the filenames differ and
        // neither call is skipped.
        if (!is_final_result || m.iter_n != 0) {
            int i_radial = 0;                                           // sea level
//            int i_radial = m.im-1;
            m.paraview_vtk_radial(bathymetry_name, i_radial, m.iter_n);

            int j_longal = 62;                                          // Mount Everest/Himalaya
            m.paraview_vtk_longal(bathymetry_name, j_longal, m.iter_n);

            int k_zonal = 87;                                           // Mount Everest/Himalaya
            m.paraview_vtk_zonal(bathymetry_name, k_zonal, m.iter_n);
        }

        // Panorama VTS fires whenever iter_n is a multiple of panorama_print.
        // (The previous panorama_cnt counter was incremented at the end of each iter and
        // could never reach panorama_print inside writeFile due to an off-by-one, and
        // it also failed to align with heavy_block_stride during inviscid spin-up.)
        if (m.paraview_panorama_vts_flag && m.panorama_print > 0
            && m.iter_n > 0 && m.iter_n % m.panorama_print == 0) {
            m.paraview_panorama_vts(bathymetry_name, m.iter_n);
//            m.paraview_sphere_vts(bathymetry_name, m.iter_n);
        }

        m.AtmosphereDataTransfer(bathymetry_name);
        m.AtmospherePlotData(bathymetry_name);

        cout << endl << "      AGCM: write_file ended " << endl;
    }

    // ------------------------------------------------------------------
    void resetArrays()
    {
        using namespace std;
        cout << endl << endl << endl << "      AGCM: reset_arrays" << endl;

        // 1D arrays.
        // rad/the/phi are NOT reset to the old 1.0/2.0/3.0 sentinels here. Those were overwritten
        // moments later by RunTimeSlice anyway, and rad's sentinel was exactly r0 = 1.0, so a
        // reordering would have left the metric at r = 1 everywhere while still looking sane.
        // They now come from the single definition instead.
        m.initGridCoordinates();
        m.aux_grad_v.initArray_1D(m.im, 4.0);

        // 2D arrays — all initialized to 0.0
        Array_2D* arrays_2d[] = {
            &m.Topography, &m.Vegetation, &m.precipitable_water, &m.precipitation_NASA,
            &m.temperature_NASA, &m.velocity_v_NASA, &m.velocity_w_NASA, &m.temp_pot,
            &m.temp_reconst, &m.temp_landscape, &m.p_stat_landscape, &m.relative_humidity,
            &m.Q_radiation_2D, &m.Q_latent_2D, &m.Q_sensible_2D,
            &m.Q_bottom_2D, &m.vapour_evaporation, &m.Evaporation_Dalton,
            &m.Evaporation_Meyer, &m.Evaporation_Rohwer, &m.Evaporation, &m.co2_total,
            &m.dew_point_temperature, &m.condensation_level, &m.c_fix, &m.Landscape,
            &m.Tropopause,
            &m.albedo, &m.epsilon_2D                                    // MultiLayerRadiation surface fields
        };
        // vel_star is initialised separately — must not be zero when TurbulenceAtm::init()
        // reads it before the first compute_vel_star() call.
        m.vel_star.initArray_2D(m.jm, m.km, 0.4);

        #pragma omp parallel for schedule(static)
        for (size_t n = 0; n < sizeof(arrays_2d) / sizeof(arrays_2d[0]); n++)
            arrays_2d[n]->initArray_2D(m.jm, m.km, 0.0);

        // 3D arrays initialized to 0.0
        Array* arrays_3d_zero[] = {
            &m.h, &m.u, &m.v, &m.w, &m.c, &m.cloud, &m.gr, &m.ice, &m.co2,
            &m.un, &m.vn, &m.wn, &m.cn, &m.cloudn, &m.icen, &m.grn,
            &m.p_dyn, &m.p_stat,
            &m.CloudBase, &m.LevelFreeSinking, &m.Deep_beg, &m.Deep_end,
            &m.rhs_t, &m.rhs_u, &m.rhs_v, &m.rhs_w, &m.rhs_c, &m.rhs_cloud, &m.rhs_ice, &m.rhs_g, &m.rhs_co2,
            &m.aux_u, &m.aux_v, &m.aux_w, &m.aux_t,
            &m.Q_Latent, &m.Q_Sensible,
            &m.BuoyancyForce, &m.CoriolisForce, &m.CentrifugalForce, &m.PresGradForce,
            &m.ubud_pgf, &m.ubud_cor, &m.ubud_advv, &m.ubud_advh, &m.ubud_diff, &m.ubud_buoy,
            &m.vbud_pgf, &m.vbud_cor, &m.vbud_advv, &m.vbud_advh, &m.vbud_diff, &m.vbud_other,
            &m.wbud_pgf, &m.wbud_cor, &m.wbud_advv, &m.wbud_advh, &m.wbud_diff, &m.wbud_other,
            &m.epsilon, &m.tau_above, &m.tau_layer, &m.brunt_N2, &m.Psi, &m.radiation,
            &m.P_rain, &m.P_snow, &m.P_rainn, &m.P_snown, &m.P_graupel,
            &m.P_conv,
            &m.Precipitation, &m.PrecipitableWaterLocal,
            &m.TempStand, &m.TempDewPoint,
            &m.S_v, &m.S_c, &m.S_i, &m.S_r, &m.S_s, &m.S_g, &m.S_c_c,
            &m.M_u, &m.M_d,
            &m.MC_t, &m.MC_q, &m.MC_v, &m.MC_w,
            &m.r_dry, &m.r_humid,
            &m.g_p, &m.c_u, &m.e_d, &m.e_l, &m.e_p,
            &m.u_u, &m.u_d, &m.v_u, &m.v_d, &m.w_u, &m.w_d,
            &m.q_v_u, &m.q_v_d, &m.q_c_u,
            &m.E_u, &m.D_u, &m.E_d, &m.D_d,
            // turbulence arrays
            &m.tke, &m.tken, &m.dis, &m.disn,
            &m.nue, &m.prod, &m.tke_source, &m.dis_source,
            &m.rhs_tke, &m.rhs_dis,
            &m.p_hydro, &m.PressureGradientForce
        };

        #pragma omp parallel for schedule(static)
        for (size_t n = 0; n < sizeof(arrays_3d_zero) / sizeof(arrays_3d_zero[0]); n++)
            arrays_3d_zero[n]->initArray(m.im, m.jm, m.km, 0.0);

        // 3D arrays initialized to 1.0
        Array* arrays_3d_one[] = { &m.t, &m.tn, &m.t_eq, &m.co2n, &m.s, &m.s_u, &m.s_d };

        #pragma omp parallel for schedule(static)
        for (size_t n = 0; n < sizeof(arrays_3d_one) / sizeof(arrays_3d_one[0]); n++)
            arrays_3d_one[n]->initArray(m.im, m.jm, m.km, 1.0);

        // 3D array initialized to 100.0
        m.HumidityRel.initArray(m.im, m.jm, m.km, 100.0);

        // integer 2D arrays
        for (auto& i : m.i_topography)
            std::fill(i.begin(), i.end(), 0);
        for (auto& i : m.i_tropopause)
            std::fill(i.begin(), i.end(), 0);

        cout << "      AGCM: reset_arrays ended" << endl;
    }

private:
    cAtmosphereModel& m;
};

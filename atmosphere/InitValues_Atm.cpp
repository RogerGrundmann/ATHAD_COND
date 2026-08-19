#include "MixtureAtm.h"
#include "SaturationH2O.h"
#include "cAtmosphereModel.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;


// ============================================================================
// Water Vapor and Cloud Initialization Module - Final Improved Version
// ============================================================================

#include <chrono>
#include <iostream>
#include <algorithm>
#include <cmath>

// ============================================================================
// Physical and Numerical Constants
// ============================================================================
namespace VaporCloudConstants {
    // Surface evaporation coefficients
    constexpr double COEFF_LAND = 0.74;                                 // Land surface evaporation coefficient
    constexpr double COEFF_OCEAN = 0.98;                                // Ocean surface evaporation coefficient
    
    // Moisture limits
    constexpr double Q_LIMIT = 3.0e-3;                                  // Maximum specific humidity [kg/kg]
    constexpr double Q_SCALING = 0.84;                                  // Empirical moisture scaling factor
    constexpr double FALLBACK_Q_FACTOR = 1e-5;                          // Fallback saturation factor

    // Relative humidity and cloud thresholds
    constexpr double RH_THRESHOLD = 85.0;                               // Cloud formation RH threshold [%]
    constexpr double CLOUD_SCALING = 0.0005;                            // Cloud density scaling factor
    
    // Stability parameters
    constexpr double LAPSE_RATE_REF = -0.0065;                          // Reference lapse rate [K/m]
    constexpr double STABILITY_SCALING = 500.0;                         // Stability weight scaling
    constexpr double STABILITY_IMPACT = 0.3;                            // Max stability influence
    
    // Dewpoint spread threshold
    constexpr double SPREAD_THRESHOLD = 2.0;                            // Temperature-dewpoint spread [K]

    // Magnus-Tetens formula coefficients
    // The Magnus coefficients that lived here are gone: every saturation calculation now
    // goes through SaturationH2O.h (IAPWS), which is valid across ATHAD's whole 1500 K
    // column rather than the ~320 K the Magnus fit was calibrated on.

    // Dewpoint calculation (inverse Magnus)
    constexpr double MAG_A = 17.27;
    constexpr double MAG_B = 237.3;
    constexpr double E0_HPA = 6.1078;                                   // Reference vapor pressure [hPa]
    
    // Safety limits
    constexpr double MIN_PRESSURE = 1e-10;
    constexpr double MIN_Q_SATUR = 1e-12;
}

using namespace VaporCloudConstants;

// ============================================================================
// Main Water Vapor and Cloud Initialization
// ============================================================================
void cAtmosphereModel::init_vapour_cloud() {                            // calculates initial water vapour and cloud distribution, no ice clouds are prepared
    std::cout << "\n\n\n      AGCM: init_vapour_cloud" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // Thread-local variables
    double p_u = 0.0;
    double t_u = 0.0;

    // Precompute inverse for efficiency
    const double inv_rh_range = 1.0 / (100.0 - RH_THRESHOLD);

    // ========================================================================
    // Main Computation Loop: Calculate Vapor and Cloud Fields
    // ========================================================================
    #pragma omp parallel for collapse(2) private(p_u, t_u)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {

            // Temperature thresholds (local constants for each thread)
            const double T_ice_end = t_000;
            const double T_freeze = t_0;

            const int i_mount = i_topography[j][k];

            // Process vertical column
            for (int i = 0; i < im-1; i++) {

                // Local variables for this grid point
                double E_sat_loc = 0.0;
                double q_sat_loc = 0.0;

                // Get local temperature and pressure
                t_u = t.x[i][j][k] * t_0;                               // [K]
                p_u = p_stat.x[i][j][k];                                // [hPa]

                // ------------------------------------------------------------
                // Calculate Saturation Vapor Pressure (Magnus-Tetens)
                // ------------------------------------------------------------
                const double E_wat = SaturationH2O::saturationPressure(t_u);
                const double E_ice = SaturationH2O::sublimationPressure(t_u);
                double E_sat;
                                                                                                                                                                                                        
                // Temperature-dependent phase transition
                if (t_u >= T_freeze) {
                    E_sat = E_wat;                                      // Pure liquid water
                } else if (t_u <= T_ice_end) {
                    E_sat = E_ice;                                      // Pure ice
                } else {
                    // Linear interpolation in mixed phase region
                    const double w = (t_u - T_ice_end) / (T_freeze - T_ice_end);
                    E_sat = w * E_wat + (1.0 - w) * E_ice;
                }                                                                                                                                                                                                     

                // Saturation MASS FRACTION, exact conversion. The dilute
                // ep*E/(p - E) with a fallback constant when E >= p is wrong in both
                // branches once water is the bulk gas — see SaturationH2O.h. This
                // routine is not called (see cAtmosphereModel::RunTimeSlice), but the
                // forbidden form left sitting in it is one uncomment away from being live.
                const double q_Satur = SaturationH2O::saturationMassFraction(
                    E_sat, p_u, AtmMixture::M_nonwater(c.x[i][j][k], co2.x[i][j][k],
                                                       m_comp.M_bg));

                // ------------------------------------------------------------
                // Surface Evaporation (only above topography)
                // ------------------------------------------------------------
                if (i >= i_mount) {
                    const double current_coeff = is_land(h, i, j, k) ? COEFF_LAND : COEFF_OCEAN;
                    E_sat_loc = E_sat * current_coeff;
                } else {
                    E_sat_loc = 0.0;
                }

                // Calculate specific humidity from evaporation
                if (E_sat_loc > 0.0) {
                    q_sat_loc = SaturationH2O::saturationMassFraction(
                        E_sat_loc, p_u, AtmMixture::M_nonwater(c.x[i][j][k],
                                                               co2.x[i][j][k], m_comp.M_bg));
                } else {
                    q_sat_loc = 0.0;
                }

                // ------------------------------------------------------------
                // Dewpoint Temperature and Spread
                // ------------------------------------------------------------
                const double e_actual   = (p_u * q_sat_loc) / (ep + q_sat_loc);
                const double log_val    = std::log(std::max(e_actual, MIN_PRESSURE) / E0_HPA);
//                const double t_dewpoint = (MAG_B * log_val) / (MAG_A - log_val) + 273.15;
//                const double spread     = t_u - t_dewpoint;

                // ------------------------------------------------------------
                // Atmospheric Stability Assessment
                // ------------------------------------------------------------
                // Bounds check: ensure i+1 is valid
                const double dT = (i < im-2) ? (t.x[i+1][j][k] - t.x[i][j][k]) : 0.0;
//                const double diff = dT - LAPSE_RATE_REF;

                // Stability weight: reduces moisture in stable conditions
//                const double stability_weight = 1.0 - STABILITY_IMPACT * std::tanh(diff * STABILITY_SCALING);

                // ------------------------------------------------------------
                // Cloud Formation Signal (based on relative humidity)
                // ------------------------------------------------------------
                const double rel_hum = (q_Satur > MIN_Q_SATUR) ? 
                                       (q_sat_loc / q_Satur * 100.0) : 0.0;
                double cloud_signal = (rel_hum - RH_THRESHOLD) * inv_rh_range;
                cloud_signal = std::clamp(cloud_signal, 0.0, 1.0);

                // ------------------------------------------------------------
                // Final Moisture Content (with stability and spread constraints)
                // ------------------------------------------------------------
//                const double q_final = Q_SCALING * q_sat_loc;

                // Apply moisture only when conditions are favorable:
                // - Small dewpoint spread (near saturation)
                // - Above topography
/*
                if (spread < SPREAD_THRESHOLD * stability_weight && i >= i_mount) {
                    const double q_weighted = q_final * stability_weight * cloud_signal;
                    c.x[i][j][k] = std::min(Q_LIMIT, q_weighted);
                } else {
                    c.x[i][j][k] = 0.0;
                }
*/

                const double RH_init = is_land(h, i, j, k) ? 0.60 : 0.75;
                c.x[i][j][k] = (i >= i_mount) ? RH_init * q_Satur : 0.0;                                                                                                                                         

                // Cloud density field
                cloud.x[i][j][k] = cloud_signal * CLOUD_SCALING;
                if (t_u < t_00)  cloud.x[i][j][k] = 0.0;
            }  // end i loop

            // Boundary conditions at top of domain
            c.x[im-1][j][k] = 0.0;
            cloud.x[im-1][j][k] = 0.0;

        }  // end k loop
    }  // end j loop

    // ========================================================================
    // Apply Surface Boundary Condition (copy topography values to surface)
    // ========================================================================
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            const int i_mount = i_topography[j][k];
            
            // Bounds check before accessing array
            if (i_mount >= 0 && i_mount < im && is_land(h, i_mount, j, k)) {
                c.x[0][j][k] = c.x[i_mount][j][k];
            }
        }
    }
/*
    // ========================================================================                                                                                                             
    // Post-processing: smooth c and cloud horizontally
    // ========================================================================                                                                                                             
    {           
        // Precompute Gaussian weights for the stencil (constant for all levels/cells)
        const int R = 2;
        const int D = 2 * R + 1;

        std::vector<double> gauss_w(D * D);

        for (int dj = -R; dj <= R; dj++)
            for (int dk = -R; dk <= R; dk++)
                gauss_w[(dj + R) * D + (dk + R)] =
                    std::exp(-0.5 * (dj*dj + dk*dk) / double(R*R));


        #pragma omp parallel
        {
            std::vector<double> tmp(jm * km);       // thread-local scratch buffer

            auto smoothLevel = [&](auto& field, int i) {
                for (int j = 0; j < jm; j++)
                    for (int k = 0; k < km; k++)
                        tmp[j * km + k] = field.x[i][j][k];

                for (int j = 0; j < jm; j++) {
                    for (int k = 0; k < km; k++) {
                        if (i < i_topography[j][k]) continue;

                        double sum = 0.0;
                        double cnt = 0.0;

                        for (int dj = -R; dj <= R; dj++) {
                            const int jj = std::clamp(j + dj, 0, jm - 1);

                            for (int dk = -R; dk <= R; dk++) {
                                const int kk = (k + dk + km) % km;
//                                const double w = std::exp(-0.5 * (dj*dj + dk*dk) / double(R*R));   // if not use Gaussian weights
                                const double w = gauss_w[(dj + R) * D + (dk + R)];

                                if (i >= i_topography[jj][kk]) {
                                    sum += w * tmp[jj * km + kk];
                                    cnt += w;
                                }
                            }
                        }
                        field.x[i][j][k] = cnt > 0.0 ? sum / cnt : 0.0;
                    }
                }
            };

            // Each i-level is independent: threads take different levels
            #pragma omp for schedule(dynamic, 4)
            for (int i = 0; i < im - 1; i++) {     // im-1 stays 0 (top BC)
                smoothLevel(c,     i);
                smoothLevel(cloud, i);
            }

        } // end omp parallel

        // Re-apply surface BC
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < jm; j++) {
            for (int k = 0; k < km; k++) {
                const int i_mount = i_topography[j][k];

                if (i_mount >= 0 && i_mount < im && is_land(h, i_mount, j, k))
                    c.x[0][j][k] = c.x[i_mount][j][k];
            }
        }
    }
*/

    // ========================================================================
    // Performance Timing and Output
    // ========================================================================
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" Time measured: %.3f seconds for init_vapour_cloud\n", elapsed.count() * 1e-9);

    std::cout << "      AGCM: init_vapour_cloud ended" << std::endl;
}
/*
*
*/
// ============================================================================
// Optional: Debug Output Function (compile with -DDEBUG_VAPOR_CLOUD)
// ============================================================================
#ifdef DEBUG_VAPOR_CLOUD

void cAtmosphereModel::debug_vapor_output(int i, int j, int k, 
                                          double t_u, double p_u, 
                                          double E_sat, double q_sat) {
    if ((j == 60) && (k == 87)) {
        std::cout.precision(5);
        std::cout.setf(std::ios::fixed);
        std::cout << "\n  WaterVapour Debug Output"
                  << "\n  Position: i=" << i << " j=" << j << " k=" << k
                  << "\n  p_static = " << p_u
                  << "\n  T = " << t_u << " K  (" << t_u - t_0 << " °C)"
                  << "\n  E_sat = " << E_sat
                  << "\n  0.84 * q_sat = " << q_sat
                  << "\n  Moisture c = " << c.x[i][j][k]
                  << "\n  Diff_c = " << c.x[i][j][k] - q_sat
                  << "\n" << std::endl;
    }
}
#endif
/*
*
*/
void cAtmosphereModel::init_tropopause_layers(){                                                                                                                                                         
    cout << endl << endl << endl << "      AGCM: init_tropopause_layers" << endl;                                                                                                                        
                                                                                                                                                                                                           

    int j_max = jm - 1;                                                                                                                                                                                  
    int j_half = j_max / 2;                                                                                                                                                                              
                  
    // Derive x_max so that Agnesi(tropopause_equator, x_max) == tropopause_pole exactly.                                                                                                                
    // Agnesi: a^3/(a^2+x^2) = b  =>  x = a * sqrt(a/b - 1)
    // Requires tropopause_equator > tropopause_pole (always true physically).                                                                                                                           
    double x_max = tropopause_equator                                                                                                                                                                    
                   * std::sqrt(tropopause_equator / tropopause_pole - 1.0);                                                                                                                              
                  
    // ATHAD: height -> level index must INVERT the stretched grid.
    //
    // What was here: round(h / L_atm). That is a level index only if the layers were
    // uniformly L_atm thick, and they are not — init_layer_heights builds
    //     height(i) = (exp(zeta * i / (im-1)) - 1) * L_atm
    // so L_atm is the AMPLITUDE of an exponential stretch, not a layer spacing (the same
    // confusion the metric-length comment in cAtmosphereModel.h warns about). The exact
    // inverse is
    //     i = (im-1) * ln(1 + h / L_atm) / zeta.
    //
    // The error is not small and it is not latent here. At the 300 km shell the pole's
    // 195 km convective top sits at level 52; the old formula returned round(195000/15719)
    // = 12, which is 13 km. VelocityInitializer builds the whole jet profile between the
    // surface and this level and applies a linear taper to zero from it to the domain top,
    // so the initial wind structure was compressed into the bottom 4 % of the atmosphere
    // and the remaining 96 % got the taper. On Earth (L_atm = 400 m, h = 11 km) it returned
    // 28 of 41 levels, which is wrong too — 4.9 km, not 11 — but wrong in a way that still
    // lands inside the troposphere, so it never showed.
    const double idx_scale = (double)(im - 1) / zeta;
    auto height_to_level = [&](double h) -> double {
        if(!(h > 0.0) || !(L_atm > 0.0) || !(zeta > 0.0)) return 0.0;
        const double i = round(idx_scale * std::log(1.0 + h / L_atm));
        return std::min(std::max(i, 0.0), (double)(im - 1));
    };
    // The forward map, written out rather than taken from get_layer_height(): this
    // routine runs in an omp section CONCURRENT with init_layer_heights(), so reading
    // m_layer_heights here would race the vector that fills it.
    auto level_to_height = [&](double i) {
        return (std::exp(zeta * i / (double)(im - 1)) - 1.0) * L_atm;
    };

    cout << "tropopause_pole=" << tropopause_pole
         << " x_max=" << x_max
         << " pole_index=" << height_to_level(tropopause_pole)
         << " (height " << level_to_height(height_to_level(tropopause_pole))
         << " m of " << level_to_height(im-1) << " m)" << endl;

    // Build symmetric cache of heights [m] and grid indices in one pass.
    std::vector<double> tropo_height_cache(jm);
    tropopause_layers = std::vector<double>(jm);

    for(int j = 0; j <= j_half; j++){
        double x = x_max * (double)(j_half - j) / (double)j_half;
        double h = AtomUtils::Agnesi(tropopause_equator, x);
        tropo_height_cache[j]       = h;
        tropo_height_cache[j_max-j] = h;
        tropopause_layers[j]        = height_to_level(h);
        tropopause_layers[j_max-j]  = tropopause_layers[j];
    }
                                                                                                                                                                                                           
    #pragma omp parallel for schedule(static)                                                                                                                                                            
   for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){                                                                                                                                                                     
            Tropopause.y[j][k] = tropo_height_cache[j];
        }                                                                                                                                                                                                
    }           
                                                                                                                                                                                                           
    cout << "      AGCM: init_tropopause_layers ended" << endl;                                                                                                                                          
}
/*
*
*/
// ============================================================================
// Atmosphere Model Initialization - Improved Version
// ============================================================================

#include <chrono>
#include <iostream>
#include <algorithm>
#include <cmath>

// ============================================================================
// Physical Constants (should ideally be in a separate constants header)
// ============================================================================
namespace AtmosphereConstants {
    // Temperature constants
    constexpr double BETA_COSMO = 44.0;              // [K] COSMO parameter
//    constexpr double BETA_COSMO = 42.0;              // [K] COSMO parameter
//    constexpr double BETA_COSMO = 38.0;              // [K] COSMO parameter
//    constexpr double BETA_COSMO = 35.0;              // [K] COSMO parameter
//    constexpr double BETA_COSMO = 30.0;              // [K] COSMO parameter
    constexpr double MIN_SAFE_TEMP = 100.0;          // [K] Minimum safe temperature
    
    // Humidity thresholds
    constexpr double MIN_HUMIDITY = 0.0;             // [%]
    constexpr double MAX_HUMIDITY = 100.0;           // [%]
    
    // Water/air ratio factors
    constexpr double MIN_WATER_FACTOR = 0.5;         // Minimum total water factor
}

using namespace AtmosphereConstants;

// ============================================================================
// Helper Function: Project temperature to sea level (algebraic COSMO inversion)
// ============================================================================
/**
 * Inverts the COSMO temperature profile T(h) = sqrt(T0^2 - 2*beta*g*h/R)
 * exactly, so that feeding T0 back into the vertical reconstruction reproduces
 * t_u_init at h_mt without any round-trip error.
 *
 * T0 = sqrt(t_u_init^2 + 2*beta*g*h_mt / R_Air)
 *
 * @param t_u_init  Surface temperature at mountain height h_mt [K]
 * @param h_mt      Mountain height [m]
 * @return pair<T0_sea_level [K], p0_sea_level [hPa]>
 */
std::pair<double, double> project_to_sea_level(
    double t_u_init,
    double h_mt,
    double beta,
    double R_Air,
    double r_air,
    double g)
{
    double T0 = sqrt(t_u_init * t_u_init + (2.0 * beta * g * h_mt) / R_Air);
    double p0 = 1e-2 * (r_air * R_Air * T0);
    return {T0, p0};
}

// ============================================================================
// Main Initialization Function
// ============================================================================
void cAtmosphereModel::initTemperatureData(int Ma) {
    std::cout << "\n\n\n      AGCM: initTemperatureData" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // ========================================================================
    // ATHAD: prescribed Hadean surface temperature
    // ========================================================================
    //
    // The inherited Earth code built t.x[0] from the NASA surface field (Ma == 0) or a
    // Scotese pole→equator parabola (Ma > 0), with EarthByte reconstruction corrections
    // between slices. None of that exists for 4.4 Ga, so ATHAD prescribes the surface
    // temperature directly from two parameters and keeps the parabola SHAPE only:
    //
    //     T_s(j) = (t_surf_pole - t_surf_equator) * parabola(ratio) + t_surf_pole
    //
    // with parabola(x) = x² - 2x, ratio = j / ((jm-1)/2), so ratio = 1 at the equator
    // (parabola = -1, giving t_surf_equator) and ratio = 0 or 2 at the poles (parabola = 0,
    // giving t_surf_pole). This is the identical algebra the paleo branch used, which keeps
    // the correspondence with ATOM_Precipitation readable.
    //
    // t.x[0] is stored non-dimensional (T / t_0), the convention Step 8 reads back.
    //
    // ASSUMPTION, not a result: a 250 bar steam atmosphere is optically thick enough that
    // the equator-pole surface contrast should be small, and 1500/1450 K is a guess at it.
    // See CLAUDE.md, "Stated assumptions".

    const double inv_t0    = 1.0 / t_0;
    const double d_j_half  = 0.5 * (jm - 1);

    double t_surf_sum = 0.0;

    #pragma omp parallel for collapse(2) reduction(+: t_surf_sum)
    for (int k = 0; k < km; k++) {
        for (int j = 0; j < jm; j++) {
            const double ratio = (double)j / d_j_half;
            const double T_s   = (t_surf_pole - t_surf_equator) * AtomUtils::parabola(ratio)
                               + t_surf_pole;                            // [K]

            t.x[0][j][k] = T_s * inv_t0;                                 // non-dimensional
            t_surf_sum  += T_s;
        }
    }

    // Diagnostic mean, in °C to match the units the reporting code expects.
    t_global_mean = t_surf_sum / (double)(jm * km) - t_0;

    std::cout.precision(4);
    std::cout << "\n       ATHAD: prescribed Hadean surface temperature"
              << "\n       equator ......................................... t_surf_equator   = "
              << t_surf_equator << " K"
              << "\n       pole ............................................ t_surf_pole      = "
              << t_surf_pole << " K"
              << "\n       area mean ....................................... t_global_mean    = "
              << t_global_mean + t_0 << " K\n\n";

    // ========================================================================
    // Step 8: Vertical Temperature Profile & Potential Temperature
    // ========================================================================
    const double R_mix = m_comp.R_mix;
    const double R_bg  = m_comp.R_bg;
    const double M_bg  = m_comp.M_bg;

    #pragma omp parallel for collapse(2)
    for (int k = 0; k < km; k++) {
        for (int j = 0; j < jm; j++) {
            int i_mount = i_topography[j][k];
            double t_u_init = t.x[0][j][k] * t_0;                       // [K]

            // ATHAD: no land, so no land-sea thermal contrast and no elevation to project
            // through. i_mount is 0 in every column, so the surface IS the reference level
            // and t_u_init stands as prescribed by initTemperatureData.
            //
            // Surface pressure uses R_mix, the COLUMN reference gas constant, because r_air
            // was calibrated as p/(R_mix*T) — using R_Air (the non-condensable background,
            // 317.3) here instead gives 204 bar rather than the intended 250.
            p_stat.x[0][j][k] = 1e-2 * (r_air * R_mix * t_u_init);
            t.x[0][j][k]      = t_u_init;

            temp_pot.y[j][k]     = t.x[0][j][k];
            temp_reconst.y[j][k] = t.x[0][j][k];

            // ================================================================
            // Surface Humidity Calculation (Magnus Formula)
            // ================================================================
            const double t_u_0 = t_u_init;
            const double p_u_0 = p_stat.x[0][j][k];
            
            // Magnus coefficients depend on phase (water vs ice)
            const double E_Satur = SaturationH2O::saturationPressureAuto(t_u_0);
            const double e_curr  = c.x[0][j][k] * p_u_0 / (c.x[0][j][k] + ep);

            relative_humidity.y[j][k] = std::clamp(
                (e_curr / E_Satur) * 100.0,
                MIN_HUMIDITY,
                MAX_HUMIDITY
            );

            // ================================================================
            // Vertical Profile: Temperature, Pressure, Density
            // ================================================================
            // ATHAD: dry adiabat + isothermal skin, integrated hydrostatically.
            // Identical scheme to ThermoAtm::densities() — see the long note there for why
            // the inherited COSMO sqrt profile was replaced rather than re-tuned.
            double T_prev = std::max(MIN_SAFE_TEMP, t.x[0][j][k]);
            double p_prev = p_stat.x[0][j][k];

            for (int i = 0; i < im; i++) {
                const double q_v   = c.x[i][j][k];
                const double q_c   = co2.x[i][j][k];
                const double R_loc = AtmMixture::R_of(q_v, q_c, R_bg);

                double T_curr, p_val;
                if (i == 0) {
                    T_curr = T_prev;
                    p_val  = p_prev;
                } else {
                    const double dz     = get_layer_height(i) - get_layer_height(i-1);
                    const double cp_loc = AtmMixture::cp_of(q_v, q_c, T_prev, M_bg);
                    const double T_ad   = T_prev - (g / cp_loc) * dz;
                    T_curr = std::max(t_skin, T_ad);
                    const double T_mean = 0.5 * (T_prev + T_curr);
                    p_val  = p_prev * exp(-g * dz / (R_loc * T_mean));
                }

                t.x[i][j][k]      = T_curr;
                p_stat.x[i][j][k] = p_val;

                const double total_water_factor =
                    std::max(MIN_WATER_FACTOR, 1.0 - (cloud.x[i][j][k] + ice.x[i][j][k]));
                const double R_dloc = AtmMixture::R_of(0.0, q_c, R_bg);

                r_dry.x[i][j][k]   = p_val * 100.0 / (R_dloc * T_curr);
                r_humid.x[i][j][k] = p_val * 100.0 / (R_loc * T_curr * total_water_factor);

                T_prev = T_curr;
                p_prev = p_val;
            }

            temp_landscape.y[j][k]   = t.x[i_mount][j][k] - t_0;        // [°C]
            p_stat_landscape.y[j][k] = p_stat.x[i_mount][j][k];
        }
    }

    // ========================================================================
    // Step 9: Convert to Non-Dimensional and Apply Boundary Conditions
    // ========================================================================
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            temp_reconst.y[j][k] = t.x[0][j][k] - t_0;  // [°C]

            // Convert all vertical levels to non-dimensional
            for (int i = 0; i < im; i++) {
                t.x[i][j][k] *= inv_t0;
            }

            // Apply topography boundary condition
            int i_mount = i_topography[j][k];
            if (i_mount >= 0 && i_mount < im) {
                t.x[0][j][k] = t.x[i_mount][j][k];
                p_stat.x[0][j][k] = p_stat.x[i_mount][j][k];
                r_dry.x[0][j][k] = r_dry.x[i_mount][j][k];
                r_humid.x[0][j][k] = r_humid.x[i_mount][j][k];
            }
        }
    }

    // ========================================================================
    // Step 9b (OPTION A): smooth the topography-draped surface temperature.
    // The DEM enters the lapse through the INTEGER level i_topography, so the draped
    // terrain-surface temperature has grid-scale steps (~one layer ~400 m ~2 K) at
    // cliffs/coasts -- exactly where this model's 2dx coastal/pressure modes historically
    // ignite. Apply a few light 1-2-1 passes (phi periodic, poles fixed) to the 2D surface
    // field and write it back into the PROGNOSTIC surface cell t.x[i_mount]; bcSolidGround
    // copies i_mount->0 every step, so smoothing i=0 alone would be undone on iteration 1.
    // Paleo only -- the modern NASA field and ocean cells (i_mount=0) are left untouched.
    if (*get_current_time() != 0) {
        const int n_passes = 4;
        std::vector<double> Tsurf(jm * km), Ttmp(jm * km);
        for (int j = 0; j < jm; j++)
            for (int k = 0; k < km; k++)
                Tsurf[j * km + k] = t.x[i_topography[j][k]][j][k];

        for (int p = 0; p < n_passes; p++) {
            for (int j = 0; j < jm; j++)                                // phi (k) pass, periodic
                for (int k = 0; k < km; k++) {
                    int km1 = (k - 1 + km) % km, kp1 = (k + 1) % km;
                    Ttmp[j * km + k] = 0.25 * Tsurf[j * km + km1]
                                     + 0.50 * Tsurf[j * km + k]
                                     + 0.25 * Tsurf[j * km + kp1];
                }
            for (int k = 0; k < km; k++) {                              // theta (j) pass, poles held
                Tsurf[k]              = Ttmp[k];
                Tsurf[(jm - 1) * km + k] = Ttmp[(jm - 1) * km + k];
                for (int j = 1; j < jm - 1; j++)
                    Tsurf[j * km + k] = 0.25 * Ttmp[(j - 1) * km + k]
                                      + 0.50 * Ttmp[j * km + k]
                                      + 0.25 * Ttmp[(j + 1) * km + k];
            }
        }

        for (int j = 0; j < jm; j++)
            for (int k = 0; k < km; k++) {
                int i_mount = i_topography[j][k];
                t.x[i_mount][j][k]     = Tsurf[j * km + k];             // prognostic surface cell
                t.x[0][j][k]           = Tsurf[j * km + k];             // sea-level reference layer
                temp_landscape.y[j][k] = Tsurf[j * km + k] * t_0 - t_0; // keep the diagnostic in sync [degC]
            }
    }

/*
    // ========================================================================
    // Post-processing: smooth t, p_stat, r_dry, r_humid horizontally                                        
    // ========================================================================
    {
        // Precompute Gaussian weights for the stencil (constant for all levels/cells)
        const int R = 2;
        const int D = 2 * R + 1;

        std::vector<double> gauss_w(D * D);

        for (int dj = -R; dj <= R; dj++)
            for (int dk = -R; dk <= R; dk++)

              gauss_w[(dj + R) * D + (dk + R)] =
                    std::exp(-0.5 * (dj*dj + dk*dk) / double(R*R));


        #pragma omp parallel
        {
            std::vector<double> tmp(jm * km);                           // thread-local scratch: one slice per thread

          // Lambda captures thread-local tmp — safe to call from any thread
            auto smoothLevel = [&](auto& field, int i, auto pred) {
                // Snapshot level i before writing (read-then-write on same array)
                for (int j = 0; j < jm; j++)
                    for (int k = 0; k < km; k++)
                        tmp[j * km + k] = field.x[i][j][k];

                for (int j = 0; j < jm; j++) {
                    for (int k = 0; k < km; k++) {
                        if (!pred(j, k)) continue;

                        double sum = 0.0;
                        double cnt = 0;

                        for (int dj = -R; dj <= R; dj++) {
                            const int jj = std::clamp(j + dj, 0, jm - 1);

                            for (int dk = -R; dk <= R; dk++) {
                                const int kk = (k + dk + km) % km;

                                if (i >= i_topography[jj][kk]) {
//                                    const double w = std::exp(-0.5 * (dj*dj + dk*dk) / double(R*R));   // if not use Gaussian weights
                                    const double w = gauss_w[(dj + R) * D + (dk + R)];
                                    sum += w * tmp[jj * km + kk];
                                    cnt += w;
                                }
                            }
                        }
                        field.x[i][j][k] = cnt > 0.0 ? sum / cnt : 0.0;
                    }
                }
            };

            // --- t and p_stat: above topography only ---
            // Each i-level is independent: no data dependency between levels
            #pragma omp for schedule(dynamic, 4)
            for (int i = 1; i < im - 1; i++) {
                auto above_topo = [&](int j, int k) {
                    return i >= i_topography[j][k];
                };

                smoothLevel(t,      i, above_topo);
                smoothLevel(p_stat, i, above_topo);
            }
            // implicit barrier: all t/p_stat levels done before r_dry/r_humid start

            // --- r_dry and r_humid: ocean cells only (land stays 0) ---
            #pragma omp for schedule(dynamic, 4)
            for (int i = 1; i < im - 1; i++) {
                auto ocean_above_topo = [&](int j, int k) {
                    return i >= i_topography[j][k] && !is_land(h, i, j, k);
                };

                smoothLevel(r_dry,   i, ocean_above_topo);
                smoothLevel(r_humid, i, ocean_above_topo);
            }
        } // end omp parallel — implicit barrier before BC loop

        // Re-apply i=0 boundary conditions from smoothed i_mount values
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < jm; j++) {
            for (int k = 0; k < km; k++) {
                const int im0 = i_topography[j][k];

                if (im0 >= 0 && im0 < im) {
                    t.x[0][j][k]       = t.x[im0][j][k];
                    p_stat.x[0][j][k]  = p_stat.x[im0][j][k];
                    r_dry.x[0][j][k]   = r_dry.x[im0][j][k];
                    r_humid.x[0][j][k] = r_humid.x[im0][j][k];
                }

                temp_reconst.y[j][k]     = t.x[0][j][k] * t_0 - t_0;
                temp_landscape.y[j][k]   = t.x[im0][j][k] * t_0 - t_0;
                p_stat_landscape.y[j][k] = p_stat.x[im0][j][k];
            }
        }
    }
*/


    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "      Initialization completed in " << duration.count() << " ms\n";
}
/*
*  ATHAD: the Scotese et al. (2021) paleo-temperature curve loaders and the
*  interpolating accessor get_temperatures_from_curve() lived here. They are gone:
*  the curves stop at ~540 Ma and nothing reaches 4.4 Ga, so ATHAD prescribes its
*  surface temperature from t_surf_equator / t_surf_pole instead (initTemperatureData).
*
*  The accessor was also unsafe as written — it dereferenced m.begin() and decremented
*  m.end() BEFORE its own m.size() < 2 guard, so calling it on an empty map was
*  undefined behaviour rather than the intended NAN return. With the curves removed
*  that path was live on every printDataAtm() call.
*/
/*
*
*/
// ATHAD: the surface is entirely water, so this only confirms the invariant.
//
// The Earth version also reported the per-point CO2 budget added by the ocean and land
// surfaces and removed by vegetation. The Hadean has no vegetation, no land and no
// carbonate ocean sink, so those terms are gone with their parameters. Note this
// function divided by h_land to form the ocean/land ratio — with no land that is a
// division by zero, which is the other reason it could not be left as it was.
void cAtmosphereModel::LandOceanFraction(){

    cout << endl << endl << endl << "      AGCM: LandOceanFraction" << endl;

    int h_land = 0;
    for(int j = 0; j < jm; j++)
        for(int k = 0; k < km; k++)
            if(is_land(h, 0, j, k)) h_land++;

    const int h_point_max = jm * km;

    cout.precision(3);
    cout << endl;
    cout << setiosflags(ios::left) << setw(50) << setfill('.')
        << "      total number of surface points " << " = "
        << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << h_point_max << endl << setiosflags(ios::left) << setw(50)
        << setfill('.') << "      number of points on the water surface " << " = "
        << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << h_point_max - h_land << endl << setiosflags(ios::left) << setw(50)
        << setfill('.') << "      number of points on the land surface " << " = "
        << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << h_land << endl << endl;

    // Invariant 1 in CLAUDE.md: there is no topography, so is_land() must be false
    // everywhere. If that ever stops holding, every land branch in the RHS, the BCs,
    // the turbulence and the ice schemes silently comes back to life.
    if(h_land != 0){
        cout << "      ERROR: " << h_land << " land points found on a surface that must be "
             << "entirely water. ATHAD prescribes no topography." << endl;
        throw std::logic_error("   ATHAD: land points found on the flat Hadean surface");
    }
    cout << "      flat Hadean surface confirmed: 100% water, no land points" << endl << endl;

    cout << "      AGCM: LandOceanFraction ended" << endl;
}
/*
*
*/
void cAtmosphereModel::initWaterWapour() {
    std::cout << "\n\n\n      AGCM: initWaterWapour" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();


    // ========================================================================
    // Main Computation Loop: water vapour field
    // ========================================================================
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            int i_mount = i_topography[j][k];

            // ATHAD_COND: water vapour is a SATURATION PROFILE over a dry-atmosphere floor.
            //
            // ATHAD initialised c uniform at the composition mass fraction, because there
            // water is supercritical from the ground to ~177 km and there is no saturation
            // curve to follow. Here there is: the air over a 240 C sea is saturated, and
            // q_sat falls by two orders of magnitude up the column. A uniform field would
            // put 35 % water by mass at 60 km, where the saturation adjustment would
            // condense essentially all of it on the first iteration — a manufactured global
            // cloud of the kind ATHAD spent item 9 removing, arrived at from the other side.
            //
            //     c(i) = max( min_{0..i} q_sat(T, p),  c_h2o_dry_top )
            //
            // The floor is the water that survives above the cold trap — the literature's
            // 0.4-2 % by mole of the dry CO2/N2 mixture. Below the crossing the column is
            // saturated and the profile is set by temperature; above it the column is dry
            // and the profile is set by composition. THE CROSSING IS NOT PRESCRIBED: it is
            // wherever the two curves meet, which makes the cold-trap height a computed
            // diagnostic of this model rather than one of its inputs.
            //
            // The RUNNING MINIMUM is not decoration, and the first version of this loop got
            // it wrong. saturationMassFraction returns 1.0 when E >= p — "no condensation
            // limit", the value ATHAD needs above the critical point — and above ~75 km the
            // pressure has collapsed far enough that p_sat(254 K) exceeds it, so q_sat comes
            // back as 1.0. A plain max() then read that as "the air is pure water" and filled
            // the top 200 km of the domain with q_H2O = 1.0 and R = 356.9. The running
            // minimum is also the physics: a parcel's water content is set by the driest
            // point it has passed through, which is what a cold trap is.
            //
            // Note the saturation call takes M_nonwater (CO2 + background, 42.88 g/mol),
            // NOT M_bg (N2 alone, 28.014). The two differ by 53 % here, and using M_bg puts
            // the sea surface at q_sat = 0.448 against the true 0.346 — a 29 % error in the
            // one quantity this model turns on. ATHAD's M_nonwater carried exactly this
            // confusion until it was fixed there; the fix matters here.
            (void)i_mount;
            const double M_nw = AtmMixture::M_nonwater(c_0, co2_0, m_comp.M_bg);
            double q_sat_min  = 1.0;                                    // cold-trap minimum so far

            for (int i = 0; i < im; i++) {
                const double T_i = t.x[i][j][k] * t_0;                  // [K]
                const double p_i = p_stat.x[i][j][k];                   // [hPa]

                const double q_s = SaturationH2O::saturationMassFractionAt(T_i, p_i, M_nw);
                if (q_s < q_sat_min) q_sat_min = q_s;

                c.x[i][j][k]     = std::max(q_sat_min, c_h2o_dry_top);
                cloud.x[i][j][k] = 0.0;
            }
        }
    }

    // ========================================================================
    // Surface Boundary Condition (copy topography values to surface)
    // ========================================================================
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            const int i_mount = i_topography[j][k];
            if (i_mount >= 0 && i_mount < im && is_land(h, 0, j, k)) {
                c.x[0][j][k] = c.x[i_mount][j][k];
                cloud.x[0][j][k] = cloud.x[i_mount][j][k];
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" Time measured: %.3f seconds for initWaterWapour\n", elapsed.count() * 1e-9);

    std::cout << "      AGCM: initWaterWapour ended" << std::endl;
}
/*
*
*/
// ===========================================================================
// initCO2 — the CO2 initial condition. Ported from ATHAD (its README items 56-57, 59), and
// the port is a SIMPLIFICATION here rather than a new idea, because this fork had already
// found half of it.
//
// WHAT WAS HERE. ThermoAtm::co2Atmosphere() stored the local mass fraction directly,
//
//     co2(i) = (1 - c(i)) * f_CO2,     f_CO2 = q_CO2/(q_CO2 + q_bg) at the sea = 0.9536
//
// with a comment explaining exactly why a uniform MASS fraction is wrong here: the water runs
// from 0.346 at the sea to 0.004 above the cold trap, so holding the mass fraction uniform
// makes the background absorb 0.34 of the mass and puts R aloft at 230 J/(kg K) instead of
// 195 — an 18 % error through the whole upper atmosphere. That reasoning was right.
//
// WHAT IT COULD NOT DO. It baked the distribution in ONCE, against the initial water field. As
// soon as the water evolves, the stored field is stale and the background silently resumes
// absorbing every change — the same defect, deferred by one initialisation.
//
// AtmMixture::q_CO2_of now applies q_CO2 = co2_stored*(1-q_v)/(1-c_0) CONTINUOUSLY, so the
// invariant is maintained instead of imprinted. The stored field therefore becomes UNIFORM at
// co2_0 and the height dependence is produced on demand. The two agree exactly at t = 0:
//
//     co2_0/(1 - c_0) = 0.6233/0.6536 = 0.9536 = f_CO2
//
// so q_CO2_of(c, co2_0) IS (1 - c)*f_CO2, the identical field, and it stays that field
// afterwards. One mechanism, applied continuously, replacing one applied once.
//
// ORDERING IS A CONTRACT (ATHAD item 22): this must precede initTemperatureData and
// densities(), which read the composition through R_of and cp_of. And note what the change
// buys upstream: co2Atmosphere() had to sit INSIDE the two-pass loop with densities() because
// the stored field followed c. It no longer does, so the CO2 half of that iteration is gone.
// ===========================================================================
void cAtmosphereModel::initCO2() {
    std::cout << "\n\n\n      AGCM: initCO2" << std::endl;
    auto begin = std::chrono::high_resolution_clock::now();

    const double co2_ref = co2_0 * co2_scale;                           // [kg/kg] at the sea

    // ATM_CO2_INIT_PERTURB — amplitude of a DELIBERATELY ARTIFICIAL vertical gradient on the
    // initial field, default 0.0 (bit-identical: the factor is exactly 1.0). A MEASUREMENT
    // TOOL, not a claim: a source-free tracer laid down uniformly has nothing to transport, so
    // the CO2 transport and its convective redistribution are untested rather than tested and
    // found working. Anything measured with this set describes the numerics.
    static const double perturb = [](){
        const char* e = getenv("ATM_CO2_INIT_PERTURB"); return e ? atof(e) : 0.0; }();

    const double z_top = std::max(1.0, (double)get_layer_height(im - 1));

    #pragma omp parallel for collapse(2)
    for (int j = 0; j < jm; j++)
        for (int k = 0; k < km; k++)
            for (int i = 0; i < im; i++) {
                // linear in TRUE height, not in the grid index — the stretch makes those very
                // different things.
                const double shape = 1.0 - 2.0 * (double)get_layer_height(i) / z_top;
                co2.x[i][j][k] = co2_ref * (1.0 + perturb * shape);
            }

    std::cout.precision(6);
    if (perturb != 0.0)
        std::cout << "      AGCM: co2 initialised with an ARTIFICIAL vertical gradient, amplitude "
                  << perturb << " — transport test only" << std::endl;
    else
        std::cout << "      AGCM: co2 uniform at " << co2_ref
                  << " kg/kg at the sea-surface water content; q_CO2_of makes the DRY-AIR ratio "
                  << "uniform with height (co2_0 = " << co2_0
                  << ", co2_scale = " << co2_scale << ")" << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" Time measured: %.3f seconds for initCO2\n", elapsed.count() * 1e-9);
    std::cout << "      AGCM: initCO2 ended" << std::endl;
}
/*
*
*/
void cAtmosphereModel::initCloudIce() {
    std::cout << "\n\n\n      AGCM: initCloudIce" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    const double alfa_s    = 1.5;
    const double Hu_cr_max = 1.0;
    const double Hu_cr_mid = 0.8;
    const double Hu_diff   = Hu_cr_max - Hu_cr_mid;
//    const double det_T_0   = t_0 - 3.0;
    const double det_T_0   = t_0;
    // Parabola H_crit(p): roots at p=0 and p=p_crit, minimum Hu_cr_mid at p=p_mid.
    // H_crit = Hu_cr_max - Hu_curv * x * (1 - x),  x = p / p_crit
    //
    // ATHAD_COND: p_crit and p_mid were 1000 and 550 hPa — Earth's surface pressure and
    // its mid-troposphere — so the parabola's shape was pinned to Earth's column depth.
    // At 60 bar the whole 0-to-p_crit range would be spent inside the top 53 km, x would
    // exceed 1 through the entire troposphere, and x*(1-x) would go NEGATIVE, turning the
    // critical humidity from a reduction into an amplification exactly where the cloud is.
    // Anchoring both to the surface pressure keeps the parabola the same SHAPE in
    // fractional depth, which is what it was always meant to be, and reproduces the Earth
    // values at an Earth surface pressure.
    const double p_crit     = p_0;                                      // was 1000 hPa
    const double p_mid      = 0.55 * p_0;                               // was 550 hPa
    const double x_mid      = p_mid / p_crit;                           // 0.55
    const double Hu_curv    = Hu_diff / (x_mid * (1.0 - x_mid));        // ~0.808
    const double inv_p_crit = 1.0 / p_crit;

    // ========================================================================
    // Pass 1: cloud_max[i] — parallel over i, sum thread-local
    // ========================================================================
    cloud_max = std::vector<double>(im, 0.0);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < im; i++) {
        double sum = 0.0;
        for (int j = 0; j < jm; j++) {
            for (int k = 0; k < km; k++) {
                const double t_u = t.x[i][j][k] * t_0;
                const double p_u = p_stat.x[i][j][k];

                // ATHAD: nothing condenses above the critical point, and nothing condenses
                // where the vapour is superheated either.
                //
                // The critical-point guard below was the first half of this fix. The second
                // half is the conversion itself: q_sat = ep*E_sat/(p_u - E_sat) goes NEGATIVE
                // whenever E_sat exceeds the local pressure, and (c - H_crit*q_sat) then
                // manufactures cloud out of the subtraction of a negative number. The guard
                // stopped that above 647 K but not in the SUBcritical band between ~373 K and
                // the critical point, where p_sat still exceeds the local pressure through the
                // whole upper column — 140 to 200 km here. The exact form returns 1 there,
                // which is what "no condensation is possible" actually means.
                if (t_u >= AtmMixture::T_CRIT_H2O) continue;

                const double E_sat = (t_u >= t_0)
                    ? SaturationH2O::saturationPressure(t_u)
                    : SaturationH2O::sublimationPressure(t_u);
                const double q_sat = SaturationH2O::saturationMassFraction(
                    E_sat, p_u, AtmMixture::M_nonwater(c.x[i][j][k], co2.x[i][j][k],
                                                       m_comp.M_bg));

//                sum += std::max(0.0, c.x[i][j][k] - q_sat); }
//                sum += std::max(0.0, c.x[i][j][k] - 0.84 * q_sat); }
                sum += std::max(0.0, c.x[i][j][k] - 0.74 * q_sat); }
        }
        cloud_max[i] = sum / ((jm-1) * (km-1));
        if (cloud_max[i] <= 1e-4)  cloud_max[i] = 1e-4;
    }

    // ========================================================================
    // Precompute per-level quantities (sequential: im is small)
    // ========================================================================
    std::vector<double> step(im, 0.0);
    std::vector<double> dt_dim(im, 0.0);
    std::vector<double> alfa_over_cmax(im, 0.0);
    std::vector<double> two_step(im, 0.0);

    for (int i = 0; i < im - 1; i++) {
        step[i]           = get_layer_height(i+1) - get_layer_height(i);
        dt_dim[i]         = step[i] / 1.6;
        alfa_over_cmax[i] = alfa_s / cloud_max[i];
        two_step[i]       = 2.0 * step[i];
    }
    step[im-1]           = step[im-2];
    dt_dim[im-1]         = dt_dim[im-2];
    alfa_over_cmax[im-1] = alfa_over_cmax[im-2];
    two_step[im-1]       = two_step[im-2];

    // ========================================================================
    // Pass 2: cloud, ice, graupel fields — collapse(j, k), i inner
    // ========================================================================
    #pragma omp parallel for collapse(2) schedule(static)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            for (int i = 0; i < im; i++) {
                const double t_u = t.x[i][j][k] * t_0;
                const double p_u = p_stat.x[i][j][k];

                // ATHAD: nothing condenses above the critical point, and nothing condenses
                // where the vapour is superheated either.
                //
                // The critical-point guard below was the first half of this fix. The second
                // half is the conversion itself: q_sat = ep*E_sat/(p_u - E_sat) goes NEGATIVE
                // whenever E_sat exceeds the local pressure, and (c - H_crit*q_sat) then
                // manufactures cloud out of the subtraction of a negative number. The guard
                // stopped that above 647 K but not in the SUBcritical band between ~373 K and
                // the critical point, where p_sat still exceeds the local pressure through the
                // whole upper column — 140 to 200 km here. The exact form returns 1 there,
                // which is what "no condensation is possible" actually means.
                if (t_u >= AtmMixture::T_CRIT_H2O) {
                    cloud.x[i][j][k] = 0.0;
                    ice.x[i][j][k]   = 0.0;
                    gr.x[i][j][k]    = 0.0;
                    continue;
                }

                const double E_sat = (t_u >= t_0)
                    ? SaturationH2O::saturationPressure(t_u)
                    : SaturationH2O::sublimationPressure(t_u);
                const double q_sat = SaturationH2O::saturationMassFraction(
                    E_sat, p_u, AtmMixture::M_nonwater(c.x[i][j][k], co2.x[i][j][k],
                                                       m_comp.M_bg));

                const double x_norm = p_u * inv_p_crit;
                double H_crit = Hu_cr_max - Hu_curv * x_norm * (1.0 - x_norm);
                if (H_crit > 1.0)  H_crit = 1.0;

                const double del_q_ls = std::max(0.0, c.x[i][j][k] - H_crit * q_sat);
                const double cloud_ls = cloud_max[i] * (1.0 - exp(-alfa_over_cmax[i] * del_q_ls));

                double cloud_conv = 0.0;
                if (P_rain.x[i][j][k] > 0.0 && i < im - 1) {
                    const double del_q_conv = std::max(0.0,
                        (P_rain.x[i+1][j][k] - P_rain.x[i][j][k])
                        / (two_step[i] * r_humid.x[i][j][k]) * dt_dim[i]);
                    cloud_conv = cloud_max[i] * (1.0 - exp(-alfa_over_cmax[i] * del_q_conv));
                }

                double cloud_val = cloud_ls + cloud_conv;
                if (is_land(h, i, j, k))  cloud_val = 0.0;
                cloud.x[i][j][k] = cloud_val;
 
                double h_T = 0.0;
                if (t_u < t_0) {
                    const double ratio = (t_u - t_0) / det_T_0;
                    h_T = 1.0 - exp(-0.5 * ratio * ratio);
                }

                ice.x[i][j][k] = cloud_val * h_T;
                gr.x[i][j][k]  = 0.1 * cloud_val * h_T;
 
                if (t_u <= t_00) {
                    cloud.x[i][j][k] = 0.0;
                    ice.x[i][j][k]   = 0.0;
                    gr.x[i][j][k]    = 0.0;
                }
            }  // i
        }  // k
    }  // j

    // ========================================================================
    // Surface Boundary Condition
    // ========================================================================
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            const int i_mount = i_topography[j][k];
//            if (i_mount >= 0 && i_mount < im && is_land(h, 0, j, k)) {
            if (i_mount >= 0 && i_mount < im && is_land(h, i_mount, j, k)) {
                cloud.x[0][j][k] = cloud.x[i_mount][j][k];
                ice.x[0][j][k]   = ice.x[i_mount][j][k];
                gr.x[0][j][k]    = gr.x[i_mount][j][k];
            }
            for (int i = i_mount-1; i >= 0; i--) {
                if (is_land(h, i, j, k)) {
                    cloud.x[i][j][k] = 0.0;
                    ice.x[i][j][k]   = 0.0;
                    gr.x[i][j][k]    = 0.0;
                }
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" Time measured: %.3f seconds for initCloudIce\n", elapsed.count() * 1e-9);

    std::cout << "      AGCM: initCloudIce ended" << std::endl;
}
/*
*
*/

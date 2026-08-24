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
        if (cap_probe) {
            printf("      S_r PROBE: max S_r %.6e at [%d][%d][%d]  t_u = %.2f K\n"
                   "                 S_c_au %.3e  S_ac %.3e  -S_ev %.3e  S_s_shed %.3e  S_g_shed %.3e\n"
                   "                 -S_r_cri %.3e  -S_r_frz %.3e  S_s_melt %.3e  S_g_melt %.3e\n",
                   g_max_Sr, g_Sr_i, g_Sr_j, g_Sr_k, g_Sr_T,
                   g_Sr_terms[0], g_Sr_terms[1], g_Sr_terms[2], g_Sr_terms[3], g_Sr_terms[4],
                   g_Sr_terms[5], g_Sr_terms[6], g_Sr_terms[7], g_Sr_terms[8]);
            g_max_Sr = 0.0;
        }
        if (cap_probe)
            printf("      CAP PROBE: P_rain cells at the cap %lld, largest value the recurrence"
                   " wanted %.6e kg/(m2 s)  (cap %.6e)\n",
                   g_cap_hits, g_max_want, 3.0e-3 * precip_cap_scale),
            g_cap_hits = 0, g_max_want = 0.0;
        cout << "      ThreeCategoryIceScheme ended" << endl;

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for ThreeCategoryIceScheme\n", elapsed.count() * 1e-9);
    }

private:
    cAtmosphereModel& m;

    // README item 76 — the precipitation-flux phase bands. Default ON;
    // ATM_PRECIP_BANDS=0 restores Earth's hard-zeroing bands. See the note at the
    // P_rain/P_snow/P_graupel assignment.
    static inline const bool precip_bands_fix = [](){
        const char* e = getenv("ATM_PRECIP_BANDS"); return e ? (atoi(e) != 0) : true; }();
    // Multiplier on P_max_flux (item 77). Default 1.0 = the shipped 3.0e-3 kg/(m2 s).
    static inline const double precip_cap_scale = [](){
        const char* e = getenv("ATM_PRECIP_CAP");
        const double v = e ? atof(e) : 1.0; return (v > 0.0) ? v : 1.0; }();
    static inline const bool cap_probe = [](){
        const char* e = getenv("ATM_ICE_CENSUS"); return e && atoi(e) != 0; }();
    static inline long long g_cap_hits = 0;
    static inline double g_max_want = 0.0;
    static inline double g_max_Sr = 0.0, g_Sr_T = 0.0, g_Sr_terms[9] = {0};
    static inline int g_Sr_i = -1, g_Sr_j = -1, g_Sr_k = -1;

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
                        // Hoisted: r_h_i and the fall-speed knob are needed by the
                        // Marshall-Palmer inversion below as well as by the residence times
                        // further down, so both are declared once here.
                        double r_h_i = m.r_humid.x[i][j][k];
                        static const bool fall_rho = [](){
                            const char* e = getenv("ATM_FALLSPEED_RHO"); return e && atoi(e) != 0; }();
                        constexpr double rho_ref_fall = 1.2041;      // kg/m3, Earth surface air

                        // ATM_PRECIP_DIMENSIONAL=1 — THE FLUX IS NORMALISED TWICE.
                        //
                        // These power laws invert the Marshall-Palmer flux/content relation. For
                        // an exponential size distribution the flux is P = B*lambda^-4.5 and the
                        // content is q = A*lambda^-4, so eliminating lambda gives
                        //
                        //     q_r = A_r * (P_rain / B_rad)^(8/9)
                        //
                        // and B_rad = rho_w*pi*N_r_0*v_r_0*Gamma(4.5)/Gamma(4) IS that physical
                        // normalisation — it is already applied, as B_rad_neg89, immediately
                        // beside the term below. Dividing by P_rain_0 as well is a SECOND
                        // normalisation, by a non-constant, non-physical quantity.
                        //
                        // The units say the same thing about accretion. c_ac is 0.24 m2/kg, so
                        // for S_ac = c_ac * q_c * X to come out in 1/s, X must carry kg/(m2 s) —
                        // the units of a mass flux. P_rain/P_rain_0 is dimensionless, so the
                        // shipped expression yields m2/kg.
                        //
                        // WHAT IT COSTS. P_rain rises monotonically downward and P_rain_0 is its
                        // surface value, so Rain runs 0 -> 1 and sits at ~1 near the ground. The
                        // flux dependence is therefore erased and S_ac collapses to c_ac*q_c. At
                        // COND's 49 g/kg that is 1.2e-2 /s, 12x over S_max; with the dimensional
                        // flux at the cap it is 1.3e-4 /s, comfortably under. ~92x.
                        //
                        // WHY IT IS LATENT UPSTREAM. ATOM_Precipitation has the identical code.
                        // On Earth q_c ~ 1e-3 kg/kg, so the inflated S_ac is ~2.4e-4 /s and never
                        // reaches the cap; here q_c is 49x larger and it does. The comment this
                        // replaces records that the normalisation was added to stop a P_snow
                        // riming overflow — i.e. a stability patch that changed the physics 92x.
                        // Both caps (S_max, P_max_flux) remain, so the dimensional branch cannot
                        // reintroduce the NaN it was guarding against.
                        //
                        // Default off; off-branch bit-identical.
                        static const bool precip_dim = [](){
                            const char* e = getenv("ATM_PRECIP_DIMENSIONAL"); return e && atoi(e) != 0; }();

                        constexpr double P_norm_floor = 1.0e-6;      // kg/(m2*s) ~0.09 mm/d
                        double P_rain_0    = precip_dim ? 1.0 : std::max(m.P_rain.x[0][j][k],    P_norm_floor);
                        double P_snow_0    = precip_dim ? 1.0 : std::max(m.P_snow.x[0][j][k],    P_norm_floor);
                        double P_graupel_0 = precip_dim ? 1.0 : std::max(m.P_graupel.x[0][j][k], P_norm_floor);

                        double Rain    = m.P_rain.x[i][j][k]    / P_rain_0;
                        double Snow    = m.P_snow.x[i][j][k]    / P_snow_0;
                        double Graupel = m.P_graupel.x[i][j][k] / P_graupel_0;

                        // pow() calls guarded — skip entirely when base is zero
                        double Rain_79  = (Rain > 0.0) ? pow(Rain, 7.0/9.0)  : 0.0;
                        double Rain_89  = (Rain > 0.0) ? pow(Rain, 8.0/9.0)  : 0.0;
                        double Rain_139 = (Rain > 0.0) ? pow(Rain, 13.0/9.0) : 0.0;
                        double Rain_16  = (Rain > 0.0) ? pow(Rain, 1.0/6.0)  : 0.0;
                        double Rain_49  = (Rain > 0.0) ? pow(Rain, 4.0/9.0)  : 0.0;

                        // ATM_FALLSPEED_RHO=1 — THE FALL SPEEDS THAT REACH THE FLUX.
                        //
                        // B_rad = rho_w*pi*N_r_0*v_r_0*Gamma(4.5)/Gamma(4), and B_s, B_g likewise,
                        // carry v_r_0 = 4.9, v_s_0 = 130.0, v_g_0 = 442.0 — the Marshall-Palmer
                        // fall-speed prefactors. THESE, not the 1.6/0.96 m/s in the residence
                        // times above, are the fall speeds with a path to the precipitation flux:
                        // they set the flux<->content inversion q = A*(P/B)^(8/9) that every
                        // collection and deposition term is built on.
                        //
                        // They are computed ONCE outside the cell loop, as const doubles at
                        // function scope, so they cannot vary with air density at all. That is an
                        // Earth assumption baked into the SHAPE of the code rather than into a
                        // literal — the sibling of item 39's exp_rm, where the defect was the
                        // structure and not the number.
                        //
                        // Terminal velocity balances gravity against drag, so v_t ~ 1/sqrt(rho_a)
                        // for a given particle. Writing f = sqrt(rho_ref/rho) for that factor,
                        // B ~ v_0 scales by f, and since B enters only as B^(-8/9) and B^(-12/13)
                        // the per-cell correction is a clean power of the density ratio:
                        //
                        //     B_rad^(-8/9)   -> B_rad^(-8/9)   * (rho/rho_ref)^(4/9)
                        //     B_s^(-12/13)   -> B_s^(-12/13)   * (rho/rho_ref)^(6/13)
                        //
                        // so the precomputed constants survive and cost one pow() per cell.
                        // COND's air is 39 kg/m3, so f = 0.176 and rain falls 5.7x slower than
                        // these constants assume; the correction RAISES the derived water content
                        // for a given flux, because slower particles need more of them.
                        //
                        // NO CLAMP HERE, unlike the residence-time branch: this is a smooth power
                        // of the density ratio with no divergence, and rho is bounded away from 0
                        // by the guard below.
                        double Brad89_loc = B_rad_neg89;
                        double Bs1213_loc = B_s_neg1213;
                        double Bg1213_loc = B_g_neg1213;
                        if(fall_rho && r_h_i > 0.0){
                            const double rr = r_h_i / rho_ref_fall;
                            Brad89_loc *= pow(rr, 4.0/9.0);
                            Bs1213_loc *= pow(rr, 6.0/13.0);
                            Bg1213_loc *= pow(rr, 6.0/13.0);
                        }

                        double r_q_r = A_r * Brad89_loc * Rain_89;
                        double r_q_s = A_s * Bs1213_loc * ((Snow > 0.0)    ? pow(Snow,    12.0/13.0) : 0.0);
                        double r_q_g = A_g * Bg1213_loc * ((Graupel > 0.0) ? pow(Graupel, 12.0/13.0) : 0.0);

                        double t_u   = m.t.x[i][j][k] * m.t_0;
                        double p_u   = m.p_stat.x[i][j][k];
                        double p_u_0 = 1e2 * p_u;
                        double c_ijk = m.c.x[i][j][k];
                        double cl_i  = m.cloud.x[i][j][k];
                        double ice_i = m.ice.x[i][j][k];

                        double step_i = step_table[i];

                        double q_sat = IceSchemeCommon::qSatWater(m, t_u, i, j, k);
                        double q_Ice = IceSchemeCommon::qSatIce(m, t_u, i, j, k);

                        // ATM_FALLSPEED_RHO=1 — DENSITY-SCALED TERMINAL FALL SPEEDS.
                        //
                        // 1.6 m/s (rain) and 0.96 m/s (snow) are Earth values, measured in air
                        // of ~1.2 kg/m3. Terminal velocity balances gravity against drag,
                        // v_t = sqrt(4 g D rho_particle / (3 C_d rho_air)), so v_t ~ 1/sqrt(rho_air)
                        // for a given particle. COND's surface air is 39 kg/m3 and PERID's 16.5,
                        // so rain there falls ~5.7x and ~3.7x slower than these constants assume.
                        //
                        // WHAT IT DOES AND DOES NOT REACH. These speeds enter ONLY through the
                        // residence times below, and every rate built on them has the form
                        // amount/dt (S_c_frz, S_nuc, S_i_melt, S_r_frz, and the deposition
                        // throttle). A slower fall is a LONGER residence, hence a SMALLER rate --
                        // the opposite of the intuition for a drop collecting as it falls, because
                        // this scheme's formulation is "the conversion completes in one residence
                        // time", not "the drop sweeps a volume". Autoconversion (c_c_au) and
                        // accretion (c_ac) do NOT use dt at all, and accretion is the term that
                        // breaches S_max, so this correction is not expected to move the
                        // precipitation much. It is made because it is right, and measured to find
                        // out, not because it is the rainfall fix.
                        //
                        // THE UPPER CLAMP IS REAL AND IS STATED RATHER THAN HIDDEN. 1/sqrt(rho)
                        // diverges as rho -> 0: at the lid rho ~ 1e-4 kg/m3 it would give 176 m/s,
                        // which is neither a terminal velocity nor a regime this drag law covers
                        // (drops break up well below it, and the continuum assumption fails).
                        // v_fac is clamped to [0.05, 1.0], i.e. the corrected speed is never
                        // FASTER than Earth's. The upper clamp binds wherever r_humid < 1.2041,
                        // which is the thin upper air where there is essentially no condensate to
                        // fall; the lower clamp corresponds to rho = 481 kg/m3 and is an inert
                        // guard. Default off; off-branch bit-identical.
                        double v_fac = 1.0;
                        if(fall_rho && r_h_i > 0.0)
                            v_fac = std::min(1.0, std::max(0.05, std::sqrt(rho_ref_fall / r_h_i)));

                        double dt_rain_dim = step_i / (1.6  * v_fac);
                        double dt_snow_dim = step_i / (0.96 * v_fac);

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
                        if (cap_probe && m.S_r.x[i][j][k] > 0.0) {
                            #pragma omp critical(srmax)
                            if (m.S_r.x[i][j][k] > g_max_Sr) {
                                g_max_Sr = m.S_r.x[i][j][k];
                                g_Sr_i = i; g_Sr_j = j; g_Sr_k = k; g_Sr_T = t_u;
                                g_Sr_terms[0]=S_c_au; g_Sr_terms[1]=S_ac; g_Sr_terms[2]=-S_ev;
                                g_Sr_terms[3]=S_s_shed; g_Sr_terms[4]=S_g_shed;
                                g_Sr_terms[5]=-S_r_cri; g_Sr_terms[6]=-S_r_frz;
                                g_Sr_terms[7]=S_s_melt; g_Sr_terms[8]=S_g_melt;
                            }
                        }
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
                        // ATHAD README item 77: the cap is a KNOB now, because it BINDS here.
                        // 3.0e-3 kg/(m2 s) ~ 260 mm/d is "well above any physical precip" on
                        // Earth; in ATHAD `max P_rain` sits exactly ON it (item 76). This file's
                        // own rule, learned from item 52's MC_t: when a cap binds, find out what
                        // it is holding back before trusting the field beneath it.
                        // ATM_PRECIP_CAP scales it; unset = 1.0 = the shipped 3.0e-3.
                        const double P_max_flux = 3.0e-3 * precip_cap_scale;   // kg/(m2*s)
                        // ATHAD README item 76: THE PHASE BANDS ARE EARTH'S, AND THEY DO TWO
                        // JOBS AT ONCE.
                        //
                        // As written, each category is `(band) ? (inherited + produced) : 0`,
                        // which conflates "can this phase be PRODUCED in this cell" with "can a
                        // flux PASS THROUGH this cell". The second has no temperature bound at
                        // all — falling ice does not cease to exist because the air it is
                        // passing through is cold — and the `: 0.0` does not merely withhold
                        // production, it DESTROYS a flux arriving from the level above.
                        //
                        // The bands bottom out at t_00 = 236.15 K (-37 C) and t_000 = 253.15 K
                        // (-20 C): Earth's homogeneous-freezing and mixed-phase thresholds. On
                        // Earth they bound a regime where precipitation is negligible anyway.
                        // ATHAD's ENTIRE condensing layer is at 221-230 K — colder than the
                        // coldest band — so all three fluxes were identically zero everywhere,
                        // in every run this project has produced, however much snow S_s made.
                        // Measured: S_s = 0.0119 g/kg/s at 277 km with P_snow = 0.000000.
                        //
                        // Fixed: the inherited flux always passes; only PRODUCTION is gated.
                        // Snow production runs wherever the air is below freezing (ice crystals
                        // form and fall at any temperature below t_0), so its t_000 floor goes.
                        // Graupel keeps its t_00 floor, which is physical rather than Earth-
                        // specific: riming needs supercooled LIQUID, and below -37 C there is
                        // none. Rain keeps t_0 for the same kind of reason.
                        //
                        // Where the flux then goes is already handled:
                        // IceSchemeCommon::evaporateWhereImpossible converts an incoming flux
                        // to vapour at P/(v*rho) when it falls into a superheated or
                        // supercritical cell, with the mass accounted for.
                        //
                        // ATM_PRECIP_BANDS=0 restores the old behaviour.
                        if (precip_bands_fix) {
                            const double prod_r = (t_u >= m.t_0)
                                ? m.r_humid.x[i+1][j][k] * m.S_r.x[i+1][j][k] * step_i : 0.0;
                            const double prod_s = (t_u <  m.t_0)
                                ? m.r_humid.x[i+1][j][k] * m.S_s.x[i+1][j][k] * step_i : 0.0;
                            const double prod_g = (t_u <  m.t_0 && t_u >= m.t_00)
                                ? m.r_humid.x[i+1][j][k] * m.S_g.x[i+1][j][k] * step_i : 0.0;
                            const double want_r = std::max(0.0, m.P_rain.x[i+1][j][k] + prod_r);
                            if (cap_probe && want_r >= P_max_flux) {
                                #pragma omp atomic
                                g_cap_hits++;
                                #pragma omp critical(capmax)
                                { if (want_r > g_max_want) g_max_want = want_r; }
                            }
                            m.P_rain.x[i][j][k]    = std::min(P_max_flux, want_r);
                            m.P_snow.x[i][j][k]    = std::min(P_max_flux, std::max(0.0,
                                m.P_snow.x[i+1][j][k]    + prod_s));
                            m.P_graupel.x[i][j][k] = std::min(P_max_flux, std::max(0.0,
                                m.P_graupel.x[i+1][j][k] + prod_g));
                        } else {
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
                        }

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

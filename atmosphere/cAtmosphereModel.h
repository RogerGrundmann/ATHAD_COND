#ifndef CATMOSPHEREMODEL_H
#define CATMOSPHEREMODEL_H

#include <fenv.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <ctime>    
#include <cmath>
#include <map>
#include <set>
#include <csignal>
#include <cstring>
#include <limits>
#include <functional>
#include <stdexcept>
#include <sys/stat.h>
#include <sys/types.h>
#include "Array.h"
#include "Array_1D.h"
#include "Array_2D.h"
#include "tinyxml2.h"
#include "Utils.h"
#include "Config.h"
#include "MixtureAtm.h"

#ifdef _OPENMP
#include <omp.h>
#endif


using namespace std;

class MoistConvection;
class ZeroCatIceScheme;
class OneCatIceScheme;
class TwoCatIceScheme;
class ThreeCatIceScheme;
class SaturationAdjustment;
class VelocityInitializer;
class PressureSolverAtm;
class ThermoAtm;
class UtilsAtm;
class BC_Atm;
class TurbulenceAtm;

class cAtmosphereModel{

    friend class MoistConvection;
    friend class ZeroCatIceScheme;
    friend class OneCatIceScheme;
    friend class TwoCatIceScheme;
    friend class ThreeCatIceScheme;
    friend class SaturationAdjustment;
    friend class VelocityInitializer;
    friend class PressureSolverAtm;
    friend class ThermoAtm;
    friend class UtilsAtm;
    friend class BC_Atm;
    friend class TurbulenceAtm;
    friend class MultiLayerRadiation;
    friend class RadiationSelfTest;   // offline standalone driver (test/ dir) — no run-loop call site

public:

    #include "AtmosphereParams.h.inc"

    double re_turb;    // set by TurbulenceAtm::init() as vel_star * z_0 / nue_air

    // Composition of the Hadean mixture, resolved once from the configured mole
    // fractions at the end of LoadConfig. Carries the mass fractions, the mean molar
    // mass and the background gas constant that R_of/cp_of need.
    AtmMixture::Composition m_comp;
    void initComposition();

    // Composition-weighted grey longwave opacity of the background, sum(f_bg[i]*kappa_i), set
    // by initComposition(). The six background gases are well mixed and source-free, so their
    // mass ratios are fixed and the per-species split collapses to this ONE number. At this
    // fork's composition the background is pure N2, so it equals kappa_N2 = 1e-6 — the same
    // value the lumped kappa_bg carried, which in ATHAD was 1867x too small. Right here, wrong
    // there, and only the split says which. ATM_BG_LUMPED=1 forces kappa_bg.
    double m_kappa_bg_eff = 0.0;
    double kappaBackground() const;

    // COSMO barometric lapse parameter, derived in initComposition() as
    // cosmo_lapse_fraction * R_mix * t_surf_equator / cp_l. Earth's hard-coded 42 K is
    // meaningless at ATHAD's R and T — see param.py, cosmo_lapse_fraction.
    double m_beta_cosmo = 42.0;

    // Initial global mass-weighted mean CO2 mass fraction, captured on the first
    // ThermoAtm::co2Column() call. CO2 has no source or sink, so this is conserved and any
    // later departure from it is transport error. Negative means "not yet captured".
    double m_q_co2_ref = -1.0;

    // Initial global mass-weighted mean TOTAL water (vapour + cloud + ice + graupel),
    // captured on the first ThermoAtm::waterBudget() call. There is no water source or
    // sink either — no evaporation from a supercritical surface, no influx — so this too
    // is conserved, and any drift is scheme error or a limiter deleting mass.
    double m_q_h2o_ref = -1.0;

    // Water mass fraction deleted by the c <= 1 - co2 ceiling, accumulated over the run.
    double m_q_h2o_clipped = 0.0;

    // Per-level water contribution at the reference time, for attributing the drift.
    std::vector<double> m_q_h2o_levels;

    // Per-level air mass [kg/m2, cos-lat weighted] at the reference time, and the water
    // mean computed against it. Together these separate "water was created" from "the
    // column was re-weighed": q_mean divides by an air mass that densities() rebuilds
    // every iteration, so a drift in it is not by itself a transport error. See
    // ThermoAtm::waterBudget().
    std::vector<double> m_air_levels;
    double m_q_h2o_fixed_ref = -1.0;

    // ---- Anelastic base state (step 2 of the anelastic scope in the README) ----
    //
    // rho_bar(z): the cos(latitude)-weighted horizontal mean of r_humid, rebuilt by
    // ThermoAtm::densities() every time it rebuilds the profile. It is a function of
    // HEIGHT ONLY on purpose — a fully three-dimensional rho makes the pressure Poisson
    // operator time-varying and costs solvability, so the anelastic projection is built
    // on the horizontal mean and the departures from it stay in the buoyancy term.
    //
    // m_dlnrho_dr is d ln(rho_bar) / d(rad.z), i.e. differentiated in the GRID radial
    // coordinate, not in height. That is what PressureSolverAtm needs: every radial
    // derivative there is an index difference times inv_2dr and then times exp_rm, so a
    // dln(rho)/dr in the same coordinate composes with the existing metric unchanged.
    // Both are empty until the first densities() call; PressureSolverAtm treats an empty
    // vector as "no base state" and falls back to the Boussinesq operator.
    std::vector<double> m_rho_base;
    std::vector<double> m_dlnrho_dr;

    // ATHAD shortwave budget, cos(latitude)-weighted over the whole sphere. albedo.y is
    // the model's OWN albedo — MultiLayerRadiation builds it from the condensate the model
    // made — so this is what the t_skin fixed point and the planetary-balance diagnostic
    // both have to be measured against. Returns false before the first radiation call,
    // when albedo.y is still zero and the numbers would be meaningless.
    bool planetaryShortWave(double& albedo_mean, double& sw_mean, double& absorbed_mean) const;

    // One relaxation step of the t_skin fixed point:
    //     sigma*t_skin^4 = (1 - albedo_mean)*sw_mean + geothermal_flux
    // No-op when t_skin_relax <= 0. See param.py (t_skin, t_skin_relax) for why the
    // reported OLR is identically sigma*t_skin^4 and what this does and does not buy.
    void updateSkinTemperature(bool report);

    cAtmosphereModel();
    ~cAtmosphereModel();

    cAtmosphereModel(const cAtmosphereModel&) = delete;
    cAtmosphereModel& operator=(const cAtmosphereModel&) = delete;

    static cAtmosphereModel* get_model(){
        if(!m_model){
            m_model = new cAtmosphereModel();
        }
        return m_model;
    }

    static const double pi180, the_degree, phi_degree, dthe, dphi, dr;
    double dt = 0.0;
    static const double the0, phi0, r0, residuum_ref_atm;

    const int c43 = 4.0/3.0, c13 = 1.0/3.0;

    // ATHAD_COND: 61 radial levels over a ~120 km shell (L_atm 6287.5, zeta 3.0 ->
    // (exp(zeta)-1)*L_atm = 120.0 km). The "~300 km" this comment used to say was ATHAD's
    // and never applied here. See param.py (L_atm, zeta) for the sizing and lib/Array.cpp
    // MAXI (= 81) for the matching assertion bound.
    //
    // im IS NOT A FREE RESOLUTION KNOB, despite dr = 1/(im-1) being derived
    // (cAtmosphereModel.cpp:57) so that the inherited dr = 0.025 coupling is already gone.
    // drag_n_layers is a count of CELLS, so changing im rescales a physical momentum sink:
    // the shipped 5.0 is 2861 m at im 41, 1786 m at 61 and 1297 m at 81. Change them
    // together, or reformulate the drag depth as a length first. ATHAD runs im = 41 for
    // wall clock (0ece6c7), not for physics, and its README item 39 later found the
    // photosphere reduced to a single grid cell there.
    static const int im = 61, jm = 181, km = 361;

    double residuum_old = 1.0e-5;

    int Ma;
    int n, n_print, n_paraview, panorama, panorama_step;


    struct CellGeometry {                                                                                                                                                                                            
        double rm, rm2, exp_rm, exp_2_rm;
        double curv;          // J'/J: the curvature term in d2f/dz2 = exp_2_rm*(f'' - curv*f'). See metricCurv().                                                                                                                                                                              
        double sinthe, sinthe2, costhe;                                                                                                                                                                                
        double inv_rm, inv_rm2;
        double inv_rmsinthe, inv_rm2sinthe, inv_rm2sinthe2;
        double costhe_inv_rm2sinthe;
        double inv_2dr, inv_2dthe, inv_2dphi;
        double inv_dr2, inv_dthe2, inv_dphi2;
    };


    bool is_final_result = false;
    bool is_print_mode = false;

    bool is_first_time_slice() const{
        if(m_time_list.empty()){
            throw("The time list is empty. It is likely the model has not started yet.");
        }
        return (m_current_time == m_time_list.begin());
    }

    std::vector<std::vector<int> > i_topography;
    std::vector<std::vector<int> > i_tropopause;
    std::vector<std::vector<int> > i_landscape;

    std::map<float,float> m;

    std::set<float>::const_iterator get_current_time() const{
        if(m_time_list.empty()){
            throw("The time list is empty. It is likely the model has not started yet.");
        }else{
            return m_current_time;
        }
    }

    std::set<float>::const_iterator get_previous_time() const{
        if(m_time_list.empty()){
            throw("The time list is empty. It is likely the model has not started yet.");
        }
        if(m_current_time != m_time_list.begin()){
            std::set<float>::const_iterator ret = m_current_time;
            ret--;
            return ret;
        }
        else{
            throw("The current time is the only time slice for now. There is no previous time yet.");
        }
    }

    /*
     * Given a latitude, return the layer index of tropopause
    */
    int get_tropopause_layer(int j){
        assert(j>=0);
        assert(j<jm);
        //refer to  BC_Thermo::TropopauseLocation and BC_Thermo::GetTropopauseHightAdd
        //tropopause height is proportional to the mean tropospheric temperature.
        //higher near the equator - warm troposphere
        //lower at the poles - cold troposphere
        return tropopause_layers[j];
    }
    /*
     *
    */
    int get_surface_layer(int j, int k){
        return i_topography[j][k];
    }
    /*
     * This function must be called after init_layer_heights()
     * Given a layer index i, return the height of this layer
    */
    float get_layer_height(int i){
        if(0>i || i>im-1){
            return -1;
        }
        return m_layer_heights[i];
    }
    std::vector<float> get_layer_heights(){
        return m_layer_heights;
    }   
    /*
    * Given a altitude, return the layer index
    */
    int get_layer_index(float height){
        std::size_t i = 0;
        for(; i<m_layer_heights.size(); i++){
            if(height<m_layer_heights[i])
                return i-1;
        }
        return i;
    }


    void LoadConfig(const char *filename);
    void Run();
    void RunTimeSlice(int time_slice);

    // THE definition of the grid coordinates. rad/the/phi are filled here and nowhere else;
    // every other site that used to write them (UtilsAtm::resetArrays, MoistConvection::precompute,
    // test/rad_selftest) calls this instead, so there is one place to change and no way for two
    // definitions to drift apart. See checkMetricConsistency() for why that matters.
    void initGridCoordinates();

    // Startup check that the radius the CORE metric uses and the radius the PHYSICS uses are the
    // same length. See the definition for the numbers and the reasoning.
    void checkMetricConsistency() const;

    // Startup check that exp_rm = 1/(rm+1) really is the Jacobian of the radial stretch it is
    // documented to be (TurbulenceAtm.h). It is not — see ATHAD README item 39; ported here
    // 2026-08-18, where the same zeta = 3.0 was live and unmeasured. Must be called AFTER
    // init_layer_heights(), since it compares against get_layer_height(). One summary line
    // always; ATM_METRIC_CHECK=1 adds the per-level table.
    void checkRadialMetric() const;

    // ATM_METRIC_RADIUS: reads the knob once and fills m_metric_r0. Call after the coordinates.
    void initMetricRadius();

    // The radius the HORIZONTAL metric should use, in rad.z units, given a grid coordinate rm.
    // Identity when the knob is off, so every call site is bit-identical by default.
    double metricRadius(double rm) const {
        return (m_metric_r0 > 0.0) ? (m_metric_r0 + (rm - rad.z[0])) : rm;
    }

    // Planetary radius in rad.z units when ATM_METRIC_RADIUS is on; 0.0 means off.
    double m_metric_r0 = 0.0;

    // How often the text diagnostics (column profile, level summary, OLR) are printed.
    //
    // Deliberately NOT tied to `checkpoint`, which controls the ParaView file writes:
    // those are expensive and wanted rarely, while the text diagnostics are cheap and
    // wanted often enough to watch a run develop. A short exploratory run should report
    // every 10 iterations; a long production run every 100, or the log becomes unreadable.
    // diagnostic_stride = 0 selects that automatically from the run length.
    int diagnosticStride() const {
        if (diagnostic_stride > 0) return diagnostic_stride;   // explicit override
        return (nm <= 100) ? 10 : 100;                         // short run : long run
    }

private:

    static cAtmosphereModel* m_model;

    string at = "AGCM";
    string bathymetry_name;

    bool has_printed_welcome_msg;


    int panorama_cnt, iter_n;

    // Inviscid spin-up state.
    // total_iter_count accumulates across time slices so the inviscid window is global, not per-slice.
    // diffusion_ramp ∈ [0,1] multiplies every diff_*_re coefficient in RHS_Atm and the no-slip flag.
    int total_iter_count = 0;
    double diffusion_ramp = 1.0;
    bool inviscid_phase = false;
    bool ubudget_capture = false;   // when true, rhs_u stores its per-term split into ubud_* (set on checkpoint iters)
    bool vbudget_capture = false;   // when true, rhs_v stores its per-term split into vbud_* (set on checkpoint iters)
    bool wbudget_capture = false;   // when true, rhs_w stores its per-term split into wbud_* (set on checkpoint iters)

    // Buoyancy ramp ∈ [0,1] — linearly increases the Boussinesq body force in rhs_u
    // from 0 at iter 0 to 1 at iter buoyancy_ramp_iters.  Set ramp_iters = 0 to
    // disable (buoyancy_ramp stays 1.0).  The earlier "500-iter ramp starved the
    // system" note applied to the OLD ~336× too-weak buoyancy; with the 2026-06-19
    // g/(omega*L_atm) scaling fix the body force is now ~336× stronger, so switching
    // it on cold CFL-blows up. Ramp it in over a few hundred iters so the pressure
    // field and circulation adjust gradually to the corrected forcing.
    double buoyancy_ramp = 1.0;
    static constexpr int buoyancy_ramp_iters = 300;

    double t_paleo_total = 0.0;
    double t_pole_total = 0.0;
    double t_global_mean = 0.0;

    // Turbulence model selection flags
    bool use_turbulence_model          = false;
    bool use_k_epsilon_turbulence_model    = false;
    bool use_k_omega_turbulence_model      = false;
    bool use_k_omega_SST_turbulence_model  = false;
    bool use_stretched_coordinate_system   = false;

    // NOTE: zeta (the radial coordinate-stretching factor) is a CONFIG PARAMETER in
    // ATHAD, declared via AtmosphereParams.h.inc — it had to become tunable to resize the
    // shell, so the hard-coded 3.715 that lived here is gone.

    // SST blending: inner (zone 1, near wall) and outer (zone 2, free stream)
    static double blend(double inner, double outer, double F1) {
        return F1 * inner + (1.0 - F1) * outer;
    }

    std::set<float> m_time_list;
    std::set<float>::const_iterator m_current_time;


    std::vector<double> alfa;
    std::vector<double> beta;

    std::vector<std::vector<double> > M_u_Base;
    std::vector<std::vector<double> > M_d_LFS;

    std::vector<std::vector<int> > i_Base;
    std::vector<std::vector<int> > i_LFS;

    std::vector<std::vector<int> > i_deep_beg;
    std::vector<std::vector<int> > i_deep_end;

    std::vector<double> short_wave_radiation;                           // lateral short wave radiation

    std::vector<double> radiation_original;                             // original radiation
    std::vector<double> tropopause_layers;                              // keep the tropopause layer index

    std::vector<double> cloud_loc;                                      // lateral cloudwater distribution
    std::vector<double> cloud_max;                                      // lateral cloud_max distribution
    std::vector<double> ice_loc;                                        // lateral ice distribution
    std::vector<double> ice_max;                                        // lateral ice_max distribution

    std::vector<double> c_land_red;
    std::vector<double> c_ocean_red;

    // CAPE removed 2026-08-19 with computeCAPE(): written, never read, wrong three ways.
    // cape_col[j][k] in MoistConvection is the real one.
    std::vector<double> K_u;
    std::vector<double> K_d;

    std::vector<double> lapse_rate;

    std::vector<float> m_layer_heights;
    std::vector<double> m_layer_J;   // dz/d(rad.z) per level; only filled for the pressure grid

    // Per-level horizontal-mean (non-dim) temperature, used as the Boussinesq
    // buoyancy base state so the body force has zero mean at every height and
    // does not fight the hydrostatic p_stat/p_dyn split. Refilled once per
    // RK4 step by computeLevelMeanTemperature().
    std::vector<double> t_ref_level;
    void computeLevelMeanTemperature();

    // Initial (non-dim) temperature at the model lid (i=im-1), snapshotted once
    // from the IC. bcRadius pins t at the lid to this so the isothermal-floor top
    // stays constant instead of drifting upward through the old cubic top
    // extrapolation. Sized [jm][km]; empty until RunStart populates it.
    std::vector<std::vector<double>> t_top_init;

    void SetDefaultConfig();
    void print_min_max_atm();
    void write_meridional_streamfunction(int iter);   // zonal-mean v + meridional mass streamfunction Ψ (Hadley/Ferrel cell diagnostic)
    // Zonal-mean meridional-wind momentum budget (Hadley/Ferrel spin-down attribution):
    // zonal_mean_v fills vbar[i][j] in m/s; write_v_momentum_budget writes the per-step
    // Δv̄ contributions (RK4 dynamics + each post-RK4 filter) differenced across one iteration.
    void zonal_mean_v(std::vector<std::vector<double> >& vbar);
    void write_v_momentum_budget(int iter,
        const std::vector<std::vector<double> >& dv_dyn,
        const std::vector<std::vector<double> >& dv_polar,
        const std::vector<std::vector<double> >& dv_orog,
        const std::vector<std::vector<double> >& dv_radial);
    // Zonal-mean ZONAL-wind (w) momentum budget — the trade / Walker component that the
    // meridional (v) budget cannot see. zonal_mean_w fills wbar[i][j] in m/s; the writer
    // attributes the per-iteration Δw̄ to RK4 dynamics + each post-RK4 filter, so the
    // trade-easterly spin-down can be pinned to a source (weakening Coriolis) or a sink
    // (a specific filter / diffusion). Same masking + scaling as the v budget.
    void zonal_mean_w(std::vector<std::vector<double> >& wbar);
    void write_w_momentum_budget(int iter,
        const std::vector<std::vector<double> >& dw_dyn,
        const std::vector<std::vector<double> >& dw_polar,
        const std::vector<std::vector<double> >& dw_orog,
        const std::vector<std::vector<double> >& dw_radial);
    void run_3D_loop(int Ma);
    void calculate_node_weights();
    void init_steps();
    void init_tropopause_layers();
    void RHS_Atmosphere_Turb(int i, int j, int k, const CellGeometry& geo);   // single dynamical core (laminar RHS_Atmosphere dropped 2026-07-08)
    void solveRungeKutta_Atmosphere_Turb();
    void fft(Array &);
    void LandOceanFraction();
    void initTemperatureData(int Ma);

    void initWaterWapour();
    void initCO2();                 // the CO2 initial condition; must precede initTemperatureData
    void initCloudIce();
    void init_vapour_cloud();
    void cloudiness_backup();
    void init_Maxwell();

    // ================= THE VERTICAL GRID =================
    //
    // ATM_GRID_PRESSURE=1 places the levels by MASS instead of by metres. Default OFF; the
    // legacy branch below is untouched and bit-identical.
    //
    // WHY. The shell is (exp(zeta)-1)*L_atm, a length, and `L_atm` has been inherited down the
    // family from ATOM, where its param.py entry still describes it as "total height is
    // 16000m*40 steps" -- i.e. it was conceived as a UNIFORM LAYER THICKNESS and only later
    // became the amplitude of an exponential stretch. A length is the wrong thing to hold fixed
    // across forks, because an atmosphere is organised by MASS and mass is exponential in
    // height with e-folding H = R*T/g. Measured across this family:
    //
    //     tree           H_surf     shell     shell/H
    //     ATOM (Earth)    8.43 km    16 km      1.90
    //     ATHAD          59.31 km   300 km      5.06
    //     ATHAD_COND     15.50 km   120 km      7.74
    //     ATHAD_PERID     7.72 km   120 km     15.54
    //
    // ATHAD_COND and ATHAD_PERID run the SAME grid (identical L_atm, zeta, im, so identical
    // 322 m bottom and 6.16 km top layers) at 7.7 and 15.5 scale heights. In ATHAD_PERID that
    // leaves 34 of 61 levels above the cold trap holding 5 % of the mass, and 13 of them above
    // 59.6 km holding 6 parts per million -- while the cloud water peaks at level 7 (2.6 km).
    //
    // WHAT THIS BRANCH DOES. Integrate a reference hydrostatic column with the model's own
    // adiabat (dT/dz = -g/cp(T), local cp, clamped at t_skin) and dp/dz = -p*g/(R*T), then
    // place level i where
    //
    //     ln(p_0/p_i) = x_i * Lambda,   Lambda = ln(p_0/p_top),
    //     x_i = (exp(beta*i/(im-1)) - 1) / (exp(beta) - 1)
    //
    // THE LADDER IS THE SAME EXPONENTIAL STRETCH THE LEGACY GRID USES -- only the coordinate
    // it stretches has changed, from metres to e-foldings of pressure. That is the whole idea
    // in one line, and it is why `beta` defaults to `zeta`: with beta = zeta the shape of the
    // grid is preserved and only its PLACEMENT moves onto the mass column.
    //
    // THE LAW MATTERS, AND A POWER LAW IS THE WRONG ONE. The first attempt here used
    // x_i = (i/(im-1))^alpha, which sets the bottom layer to (1/(im-1))^alpha and is
    // brutally sensitive to alpha. For ATHAD_PERID (im = 61, p_top/p_0 = 1e-6, H ~ 7.72 km):
    //
    //     legacy exponential-in-height ....  322 m      <- what the model runs today
    //     power law, alpha = 1.0 .......... 1778 m      equal mass per layer: far too coarse
    //     power law, alpha = 2.0 ..........   30 m      10x FINER than legacy: stiff and slow
    //     ln-p ladder, beta = zeta = 3 ....  286 m      <- matches, so this is the law used
    //
    // The 30 m case is not hypothetical: the 50 km-shell arm measured on 2026-08-24 ran 2.4x
    // finer at the bottom than its sibling and took ~7x the wall clock for the same 40
    // iterations. Refining the bottom of this grid is expensive, so the law has to hit the
    // existing spacing rather than approach it from either side.
    //
    //   ATM_GRID_PTOP   p_top/p_0, default 1e-6   (ATHAD_PERID: brings the lid 120 -> ~78 km)
    //   ATM_GRID_BETA   ln-p stretch, default = zeta
    static bool gridPressure(){
        static const bool v = [](){
            const char* e = getenv("ATM_GRID_PRESSURE"); return e && atoi(e) != 0; }();
        return v;
    }
    static double gridPTop(){
        static const double v = [](){
            const char* e = getenv("ATM_GRID_PTOP");
            const double d = e ? atof(e) : 1.0e-6;
            return (d > 0.0 && d < 1.0) ? d : 1.0e-6; }();
        return v;
    }
    // MEASURED PER TREE, 2026-08-24, and beta = zeta does NOT preserve the legacy
    // near-surface spacing outside ATHAD_PERID. The lid is invariant under beta (only
    // ATM_GRID_PTOP moves it); beta only redistributes:
    //
    //     dz_0 relative to the legacy grid
    //     beta        3.0     4.0     5.0     6.0     beta for dz_0 = legacy
    //     ATHAD      2.67x   1.29x   0.59x   0.26x          ~4.33
    //     ATHAD_COND 1.78x   0.85x   0.39x   0.17x          ~3.78
    //     ATHAD_PERID 0.88x    -       -       -            ~2.83
    //
    // WHY IT DIFFERS BY TREE. For an ISOTHERMAL column ln(p_0/p) = z/H, so a ladder uniform
    // in ln p IS a ladder uniform in z and the two grids coincide. These columns are not
    // isothermal: H = R*T/g falls with T up the column, so there are more e-foldings per
    // kilometre aloft and a ln-p ladder rides UPWARD relative to the legacy one, coarsening
    // the bottom. The size of the shift tracks the surface-to-top temperature contrast --
    // ATHAD 1500 -> 221 K (6.8x), ATHAD_COND 512 -> 221 K (2.3x), ATHAD_PERID 358 -> 221 K
    // (1.6x) -- which is the order the table shows.
    //
    // AND THAT IS WHY THIS BRANCH IS WORTH LESS HERE THAN IN ATHAD_PERID. At matched bottom
    // resolution the lid barely moves (300.0 -> 293.4 km in ATHAD, 120.0 -> 117.5 km in
    // ATHAD_COND, against 120 -> 79.6 km in ATHAD_PERID), so for these two trees the pressure
    // grid is a REDISTRIBUTION rather than a shell cut. Their shells were already sized about
    // right for their atmospheres (shell/H = 5.06 and 7.74); ATHAD_PERID's 15.54 was not.
    // Expect no payoff here and do not flip it expecting one.
    double gridBeta() const {
        const char* e = getenv("ATM_GRID_BETA");
        const double d = e ? atof(e) : zeta;
        return (d > 1.0e-6) ? d : zeta;
    }

    // Reference hydrostatic column, built once. Returns false if it cannot reach p_top.
    bool buildReferenceColumn(std::vector<double>& z_ref,
                              std::vector<double>& lnp_ref) const {
        const double T_s   = 0.5 * (t_surf_equator + t_surf_pole);   // representative surface
        const double R_loc = AtmMixture::R_of(c_0, co2_0, m_comp.R_bg);
        const double p_s   = p_0;                                    // [hPa], the anchor
        if(!(T_s > 0.0) || !(R_loc > 0.0) || !(p_s > 0.0) || !(g > 0.0)) return false;

        const double dz      = 25.0;          // m, fine enough that the ladder is smooth
        const double z_limit = 2.0e6;         // m, a guard and nothing more
        const double lnp_target = std::log(gridPTop());

        z_ref.clear(); lnp_ref.clear();
        z_ref.push_back(0.0); lnp_ref.push_back(0.0);            // ln(p/p_0) = 0 at the ground

        double T = T_s, lnp = 0.0, z = 0.0;
        while(z < z_limit && lnp > lnp_target){
            const double cp = AtmMixture::cp_of(c_0, co2_0, T, m_comp.M_bg);
            const double gamma = (cp > 0.0) ? (g / cp) : 0.0;     // dry adiabatic lapse
            double T_next = T - gamma * dz;
            if(T_next < t_skin) T_next = t_skin;                  // the same clamp densities() uses
            const double T_mid = 0.5 * (T + T_next);
            if(!(T_mid > 0.0)) break;
            lnp -= g * dz / (R_loc * T_mid);                      // d(ln p) = -g dz /(R T)
            z   += dz;
            T    = T_next;
            z_ref.push_back(z); lnp_ref.push_back(lnp);
        }
        return (lnp <= lnp_target) && (z_ref.size() > 2);
    }

    void init_layer_heights(){
        m_layer_heights.clear();
        m_layer_J.clear();

        if(gridPressure()){
            std::vector<double> z_ref, lnp_ref;
            if(buildReferenceColumn(z_ref, lnp_ref)){
                const double Lambda = -std::log(gridPTop());       // total e-foldings, > 0
                const double beta   = gridBeta();
                const double denom  = std::exp(beta) - 1.0;
                std::size_t k = 0;
                for(int i = 0; i < im; i++){
                    const double t = (double)i / (double)(im - 1);
                    const double x = (denom > 0.0)
                                   ? (std::exp(beta * t) - 1.0) / denom : t;
                    const double want = -x * Lambda;               // target ln(p/p_0), <= 0
                    while(k + 1 < lnp_ref.size() && lnp_ref[k + 1] > want) k++;
                    double zi = z_ref.back();
                    if(k + 1 < lnp_ref.size()){
                        const double d = lnp_ref[k] - lnp_ref[k + 1];
                        const double f = (d > 0.0) ? (lnp_ref[k] - want) / d : 0.0;
                        zi = z_ref[k] + f * (z_ref[k + 1] - z_ref[k]);
                    }
                    m_layer_heights.push_back((float)zi);
                }
                // Strictly increasing, or every radial derivative divides by zero.
                bool ok = true;
                for(int i = 1; i < im; i++)
                    if(!(m_layer_heights[i] > m_layer_heights[i-1])) ok = false;
                if(ok){ buildMetricTable(); return; }
                std::cout << "      ATOM: ATM_GRID_PRESSURE produced a non-monotonic grid"
                          << " - falling back to the legacy stretch" << std::endl;
                m_layer_heights.clear();
            } else {
                std::cout << "      ATOM: ATM_GRID_PRESSURE could not reach p_top ="
                          << gridPTop() << " p_0 - falling back to the legacy stretch"
                          << std::endl;
            }
        }

        float h = L_atm;
        for(int i=0; i<im; i++){
            // rad.z[0] instead of a hardcoded 1.0: the surface is wherever the radial coordinate
            // starts, which ATM_METRIC_RADIUS may move.
            m_layer_heights.push_back((exp(zeta
            * (rad.z[i] - rad.z[0])) - 1) * h);                         // in m      local atmospheric shell thickness
//            std::cout << m_layer_heights.back() << std::endl;
        }
        return;
    }

    // dz/d(rad.z) per level, central-differenced from the height table. Needed because a
    // pressure-placed grid has no closed-form Jacobian; the legacy branch keeps its analytic
    // one so it stays bit-identical.
    void buildMetricTable(){
        m_layer_J.assign(im, 0.0);
        const double dr_loc = (im > 1) ? (rad.z[im-1] - rad.z[0]) / (double)(im - 1) : 1.0;
        for(int i = 0; i < im; i++){
            const int lo = (i == 0) ? 0 : i - 1;
            const int hi = (i == im - 1) ? im - 1 : i + 1;
            const double dz = (double)m_layer_heights[hi] - (double)m_layer_heights[lo];
            const double dn = (double)(hi - lo) * dr_loc;
            m_layer_J[i] = (dn > 0.0) ? (dz / dn) : 1.0;
        }
    }

    // Physical length that ONE unit of rad.z represents, in metres.
    //
    // This is the number the metric needs and the one that is easy to get wrong. rad.z runs
    // 1.0 .. 2.0 and init_layer_heights above maps that span onto 0 .. 16.02 km, so one rad.z
    // unit is the SHELL THICKNESS, ~16 km — NOT L_atm. L_atm = 400 m is the amplitude of the
    // exponential stretch, not a grid step: the actual layer spacing runs from 39 m at the
    // surface to 1457 m at the top, and the "400 m for 1 radial step" in the dr comment
    // (cAtmosphereModel.cpp) is nominal only.
    //
    // The two readings differ by exactly 1/dr = 40, which is why the planetary radius expressed
    // in rad.z units is 6370/16.02 = 397.5 and not 6.37e6/400 = 15925.
    double metricShellLength() const {
        const double span = rad.z[im-1] - rad.z[0];
        if(!(span > 0.0)) return L_atm;
        // A pressure-placed grid has no analytic shell; read it off the table instead.
        if(gridPressure() && (int)m_layer_heights.size() == im)
            return ((double)m_layer_heights[im-1] - (double)m_layer_heights[0]) / span;
        return (exp(zeta * span) - 1.0) * L_atm / span;
    }

    // ================= THE RADIAL METRIC (README item 80) =================
    //
    // `exp_rm` is documented in TurbulenceAtm.h and PressureSolverAtm.h as the Jacobian of the
    // radial coordinate transformation. It is written `1/(rm+1)`, which is the Jacobian of a
    // QUADRATIC stretch z ~ (rm+1)^2/2, while init_layer_heights() above builds an EXPONENTIAL
    // one, z = (exp(zeta*(r-r0)) - 1)*L_atm. So every radial derivative in the core is
    // mis-scaled by a factor that varies across the column: 6.36 against 0.5 at the surface and
    // 0.317 against 0.333 at the top, i.e. 12.7x at the bottom and 0.95x at the top.
    // checkRadialMetric() has printed the spread at every startup since item 39.
    //
    // The true Jacobian is  J(r) = dz/d(rad.z) = zeta*L_atm*exp(zeta*(r-r0))  [m per rad.z unit],
    // and the core wants it as a DIMENSIONLESS factor against its own length unit, which is
    // metricShellLength() - the metres one rad.z unit represents on average.
    //
    // AND THE SECOND DERIVATIVE NEEDS A TERM THE CODE DOES NOT HAVE. With e = U/J,
    //     U^2 * d2f/dz2 = e^2 * ( d2f/dr2 - (J'/J) * df/dr )
    // and J'/J is `zeta` for the exponential stretch. The core computes `d2f/dr2 * exp_2_rm`
    // and stops, so the curvature term is missing in BOTH metrics - it is small under the
    // legacy one (J'/J = 1/(rm+1) <= 0.5) and the same order as the retained term under the
    // true one, since zeta = 3. metricCurv() returns 0 on the legacy branch so that branch
    // stays bit-identical; the missing legacy term is a separate, smaller defect, recorded
    // rather than silently fixed.
    //
    // ATM_METRIC_EXACT=1 switches to the true Jacobian. DEFAULT OFF: this moves every radial
    // derivative in the model and item 39 asks for the measurement before the flip.
    static bool metricExact(){
        static const bool v = [](){
            const char* e = getenv("ATM_METRIC_EXACT"); return e && atoi(e) != 0; }();
        return v;
    }

    // The dimensionless radial Jacobian factor a first derivative is multiplied by.
    double metricExpRm(double rm) const {
        if(!metricExact()) return 1.0 / (rm + 1.0);
        const double J = metricJ(rm);                                   // [m per rad.z unit]
        return (J > 0.0) ? (metricShellLength() / J) : (1.0 / (rm + 1.0));
    }

    // dz/d(rad.z) at rm. Analytic on the legacy exponential stretch; a table lookup on a
    // pressure-placed grid, where no closed form exists. rad.z is uniform in the index, so
    // the lookup is exact at the levels and nearest-level in between -- which is where every
    // caller evaluates it anyway.
    int metricLevelOf(double rm) const {
        const double dr_loc = (im > 1) ? (rad.z[im-1] - rad.z[0]) / (double)(im - 1) : 1.0;
        int i = (dr_loc > 0.0) ? (int)std::lround((rm - rad.z[0]) / dr_loc) : 0;
        if(i < 0) i = 0;
        if(i > im - 1) i = im - 1;
        return i;
    }
    double metricJ(double rm) const {
        if(gridPressure() && (int)m_layer_J.size() == im) return m_layer_J[metricLevelOf(rm)];
        return zeta * L_atm * exp(zeta * (rm - rad.z[0]));
    }

    // J'/J, the coefficient of the curvature term in d2f/dz2 = e^2*(f'' - curv*f').
    // Zero on the legacy branch, so the legacy operator is unchanged to the bit.
    double metricCurv(double rm) const {
        if(!metricExact()) return 0.0;
        if(gridPressure() && (int)m_layer_J.size() == im){
            // J'/J by central difference on the same table, in rad.z units.
            const int i = metricLevelOf(rm);
            const int lo = (i == 0) ? 0 : i - 1, hi = (i == im - 1) ? im - 1 : i + 1;
            const double dr_loc = (im > 1) ? (rad.z[im-1] - rad.z[0]) / (double)(im - 1) : 1.0;
            const double dn = (double)(hi - lo) * dr_loc;
            const double Ji = m_layer_J[i];
            return (dn > 0.0 && Ji > 0.0) ? ((m_layer_J[hi] - m_layer_J[lo]) / dn) / Ji : 0.0;
        }
        return zeta;
    }

    void init_topography();                                             // ATHAD: flat featureless surface, no file read
    void save_data();
    void save_array(const string& fn, const Array& a);
    std::vector<Array*> restart_arrays();   // the prognostic 3D fields a checkpoint serializes
    void save_state(int iter, int Ma);      // dump restart_arrays() + total_iter_count to a binary file (name carries Ma+iter)
    bool load_state(int iter, int Ma);      // restore them; returns false (and runs from scratch) if absent/mismatched
    void BC_pole();
    void print_welcome_msg();
    void print_final_remarks();
    void print_loop_3D_headings();
    void paraview_panorama_vts(string &Name_Bathymetry_File, int n);
    void paraview_sphere_vts(string &Name_Bathymetry_File, int n);
    void paraview_vtk_radial(string &Name_Bathymetry_File, int i_radial, int n);
    void paraview_vtk_zonal(string &Name_Bathymetry_File, int k_zonal, int n);
    void paraview_vtk_longal(string &Name_Bathymetry_File, int j_longal, int n); 
    void AtmospherePlotData(const string &Name_Bathymetry_File);
    void AtmosphereDataTransfer(const string &Name_Bathymetry_File);
    void read_Atmosphere_Surface_Data(int Ma);
    void searchMinMax_2D(const string &, const string &,
        const string &, Array_2D &, double coeff = 1.0);

    void searchMinMax_3D(const string &, const string &,
        const string &, Array &, double coeff = 1.0,
        std::function< double(double) > lambda = [](double i) -> double{return i;},
        bool print_heading = false);

public:
    Array_1D rad;                                                       // radial coordinate direction
    Array_1D the;                                                       // lateral coordinate direction
    Array_1D phi;                                                       // longitudinal coordinate direction
    Array_1D aux_grad_v;                                                // auxilliar array

    Array_2D Topography;                                                // topography
    Array_2D value_top;                                                 // auxiliar topography
    Array_2D Vegetation;                                                // vegetation via precipitation
    Array_2D precipitable_water;                                        // areas of precipitable water at the surface
    Array_2D precipitation_NASA;                                        // surface precipitation by NASA
    Array_2D temperature_NASA;                                          // surface temperature by NASA
    Array_2D velocity_v_NASA;                                           // surface v-velocity by NASA
    Array_2D velocity_w_NASA;                                           // surface w-velocity by NASA
    Array_2D temp_pot;                                                  // surface temperature by barotropic projection
    Array_2D temp_reconst;                                              // surface temperature by reconstuction tool
    Array_2D temp_landscape;                                            // landscape temperature
    Array_2D p_stat_landscape;                                          // landscape static pressure
    Array_2D relative_humidity;                                         // relative_humidity
    Array_2D Q_radiation_2D;                                            // heat from the radiation balance in [W/m2]
    Array_2D Q_latent_2D;                                               // latent heat from bottom values by the energy transport equation
    Array_2D Q_sensible_2D;                                             // sensible heat from bottom values by the energy transport equation
    Array_2D Q_bottom_2D;                                               // difference by Q_radiation - Q_latent - Q_sensible
    Array_2D albedo;                                                    // surface albedo (pole->equator parabola) — MultiLayerRadiation
    Array_2D epsilon_2D;                                                // surface emissivity (bottom layer of epsilon) — MultiLayerRadiation
    Array_2D vapour_evaporation;                                        // water vapour by evaporation in [mm/d]
    Array_2D Evaporation_Dalton;                                        // evaporation by Dalton in [mm/d]
    Array_2D Evaporation_Meyer;                                         // evaporation by Meyer (1915) in [mm/d]
    Array_2D Evaporation_Rohwer;                                        // evaporation by Rohwer (1931) in [mm/d]
    Array_2D Evaporation;                                               // evaporation by active model in [mm/d]
    Array_2D co2_total;                                                 // areas of higher co2 concentration
    Array_2D dew_point_temperature;                                     // dew point temperature
    Array_2D condensation_level;                                        // local condensation level
    Array_2D c_fix;                                                     // local surface water vapour fixed for iterations
    Array_2D Landscape;                                                 // local landscape
    Array_2D Tropopause;                                                // local tropopause
    Array_2D vel_star;                                                  // friction velocity u_tau per (j,k) column [m/s]

    Array h;                                                            // bathymetry, depth from sea level
    Array t;                                                            // temperature
    Array u;                                                            // u-component velocity component in r-direction
    Array v;                                                            // v-component velocity component in theta-direction
    Array w;                                                            // w-component velocity component in phi-direction
    Array c;                                                            // water vapour
    Array cloud;                                                        // cloud water
    Array ice;                                                          // cloud ice
    Array gr;                                                           // cloud graupel
    Array co2;                                                          // CO2
    Array tn;                                                           // temperature new
    Array t_eq;                                                         // Held-Suarez radiative-equilibrium temperature target (snapshot of the initial field)
    Array un;                                                           // u-velocity component in r-direction new
    Array vn;                                                           // v-velocity component in theta-direction new
    Array wn;                                                           // w-velocity component in phi-direction new
    Array cn;                                                           // water vapour new
    Array cloudn;                                                       // cloud water new
    Array icen;                                                         // cloud ice new
    Array grn;                                                          // cloud ice new
    Array co2n;                                                         // CO2 new
    Array p_dyn;                                                        // dynamic pressure
    Array p_stat;                                                       // static pressure
    Array TempStand;                                                    // US Standard Atmosphere Temperature
    Array TempDewPoint;                                                 // Dew Point Temperature
    Array HumidityRel;                                                  // relative humidity
    Array rhs_t;                                                        // auxilliar field RHS temperature
    Array rhs_u;                                                        // auxilliar field RHS u-velocity component
    Array rhs_v;                                                        // auxilliar field RHS v-velocity component
    Array rhs_w;                                                        // auxilliar field RHS w-velocity component
    Array rhs_c;                                                        // auxilliar field RHS water vapour
    Array rhs_cloud;                                                    // auxilliar field RHS cloud water
    Array rhs_ice;                                                      // auxilliar field RHS cloud ice
    Array rhs_g;                                                        // auxilliar field RHS cloud graupel
    Array rhs_co2;                                                      // auxilliar field RHS CO2
    Array aux_u;                                                        // auxilliar field u-velocity component
    Array aux_v;                                                        // auxilliar field v-velocity component
    Array aux_w;                                                        // auxilliar field w-velocity component
    Array aux_t;                                                        // auxilliar field t
    Array Q_Latent;                                                     // latent heat
    Array Q_Sensible;                                                   // sensible heat
    Array BuoyancyForce;                                                // buoyancy force, Boussinesque approximation
    Array CoriolisForce;                                                // Coriolis force terms
    Array CentrifugalForce;                                             // centrifugal force terms
    Array PresGradForce;                                                // Force caused by normal pressure gradient
    // Zonal-mean v momentum-budget term capture (diagnostic): per-cell rhs_v
    // RADIAL (vertical wind) momentum-budget term capture — ported from ATHAD 2026-08-18,
    // its item 42. This is the component that had NO instrument in either fork: ATHAD's
    // item 28 spurious 293-rms radial acceleration was invisible for exactly this reason,
    // and Psi is built from v alone so it cannot see a radial failure either. Added here
    // because attributing the cell decay needs all three components, not two.
    Array ubud_pgf;                                                     // -∂p/∂r ·exp_rm (radial pressure gradient)
    Array ubud_cor;                                                     // Coriolis (non-traditional; off by default)
    Array ubud_advv;                                                    // vertical advection  -u·∂u/∂r
    Array ubud_advh;                                                    // horizontal advection
    Array ubud_diff;                                                    // diffusion (molecular + turbulent, + metric)
    Array ubud_buoy;                                                    // buoyancy body force (ATHAD item 34: carries an extra *dt)
    // contributions, stored when vbudget_capture is set so write_v_momentum_budget
    // can attribute the Hadley/Ferrel spin-down to a specific dynamical term.
    Array vbud_pgf;                                                     // -∂p/∂θ /rm (meridional pressure gradient)
    Array vbud_cor;                                                     // Coriolis term (2·cosθ·w coupling to zonal wind)
    Array vbud_advv;                                                    // vertical advection  -u·∂v/∂r
    Array vbud_advh;                                                    // horizontal advection -(v/rm)∂v/∂θ -(w/rmsinθ)∂v/∂φ
    Array vbud_diff;                                                    // diffusion (molecular + turbulent, + metric terms)
    Array vbud_other;                                                   // surface drag + moist-convection momentum
    // Zonal-mean w (zonal-wind / trade) momentum-budget term capture — mirror of vbud_*,
    // stored when wbudget_capture is set so write_w_momentum_budget can attribute the
    // trade-easterly spin-down to a specific rhs_w term.
    Array wbud_pgf;                                                     // -∂p/∂φ /(rm·sinθ) (zonal pressure gradient)
    Array wbud_cor;                                                     // Coriolis term (coupling to radial/meridional wind)
    Array wbud_advv;                                                    // vertical advection  -u·∂w/∂r
    Array wbud_advh;                                                    // horizontal advection -(v/rm)∂w/∂θ -(w/rmsinθ)∂w/∂φ
    Array wbud_diff;                                                    // diffusion (molecular + turbulent, + metric terms)
    Array wbud_other;                                                   // surface drag + moist-convection momentum
    Array epsilon;                                                      // emissivity/ absorptivity
    // Ported from ATHAD's item 42 instrumentation, 2026-08-18. All four were absent here:
    // Psi existed only as a CSV of zonal means and the scalar Psi_max, and tau_above,
    // tau_layer and N2 were computed nowhere in this model at all.
    //
    // tau_above[i] = sum of the layer optical depths above level i, so it is 0 at the lid and
    // large at the surface. The level where it crosses 1 is the effective radiating level —
    // the photosphere — which is where the OLR is actually set.
    // WHY THIS IS A SEPARATE ARRAY AND NOT DERIVED FROM epsilon: epsilon = 1 - exp(-tau)
    // saturates to exactly 1.0 in double precision for tau > ~37, so tau = -ln(1 - eps) is
    // unrecoverable below the top few levels. The information is destroyed at write time; it
    // has to be accumulated where tau is still in scope.
    Array tau_above;                                                    // cumulative LW optical depth from the lid down
    // Per-LAYER long-wave optical depth, the quantity tau_above is the running sum of. It is
    // the model's own measure of how well it resolves its photosphere: dtau >~ 1 anywhere near
    // the crossing means the two-stream sweep, which is first order in dtau, is integrating
    // across an unresolved source function. ATHAD item 39 found a SINGLE cell of dtau = 55
    // spanning the crossing at im = 41, from an offline reconstruction because the model threw
    // this away. This fork runs im = 61 and has never measured it.
    Array tau_layer;                                                    // per-layer LW optical depth
    // Brunt-Vaisala frequency squared, N^2 = (g/theta) d(theta)/dz [1/s^2]. A direct test of
    // invariant 4: a column genuinely on its own integrated adiabat must have N^2 ~ 0 through
    // the convective part, so a departure there is an adiabat-integration defect, not weather.
    // NOTE this fork's invariant 4 is the MOIST adiabat, so N^2 ~ 0 is the sharper claim here.
    // RENAMED 2026-08-20 (was `N2`). Nitrogen is now written as N2 by the eight-species
    // ParaView block, and one VTK file carrying two different N2 fields is exactly the kind
    // of name that costs an afternoon. This is the frequency, not the gas.
    Array brunt_N2;                                                     // Brunt-Vaisala frequency squared [1/s2]
    // Zonal-mean meridional mass streamfunction [kg/s], replicated across k so the existing 3D
    // writers can emit it. Psi is genuinely 2D and the replication wastes memory, but it buys
    // the thing that was missing: the cells could not be looked at, only summarised by a
    // scalar. The streamfunction here was repaired in a85122c (density inside the integral);
    // that fixed what it measures, and this makes what it measures visible.
    Array Psi;                                                          // meridional mass streamfunction [kg/s]
    Array radiation;                                                    // radiation
    Array P_rain;                                                       // rain precipitation mass rate
    Array P_snow;                                                       // snow precipitation mass rate
    Array P_rainn;                                                      // rain precipitation mass rate   new
    Array P_snown;                                                      // snow precipitation mass rate   new
    Array P_graupel;                                                    // graupel precipitation mass rate
    Array P_conv;                                                       // rain formation by deep-level cloud convection
    Array Precipitation;                                                // areas of higher precipitation
    Array PrecipitableWaterLocal;                                       // precipitable water at each level
    Array S_v;                                                          // water vapour mass rate
    Array S_c;                                                          // cloud water mass rate
    Array S_i;                                                          // cloud ice mass rate
    Array S_r;                                                          // rain mass rate
    Array S_s;                                                          // snow mass rate
    Array S_g;                                                          // graupel mass rate
    Array S_c_c;                                                        // cloud water mass rate due to condensation and evaporation in the saturation adjustment technique
    Array M_u;                                                          // moist convection within the updraft
    Array M_d;                                                          // moist convection within the downdraft
    Array MC_t;                                                         // moist convection acting on dry static energy
    Array MC_q;                                                         // moist convection acting on water vapour development
    Array MC_v;                                                         // moist convection acting on v-velocity component
    Array MC_w;                                                         // moist convection acting on w-velocity component
    Array r_dry;                                                        // density of dry air
    Array r_humid;                                                      // density of humid air
    Array g_p;                                                          // conversion cloud droplets to raindrops
    Array c_u;                                                          // condensation in the updraft
    Array e_d;                                                          // evaporation of precipitation in the downdraft
    Array e_l;                                                          // evaporation of cloud water in the environment
    Array e_p;                                                          // evaporation of cloud water in the environment
    // NOT dry static energy and NOT entropy, despite both names being used for them in this
    // tree (param.py calls s_0 "entropy at 0 °C"; Results_Atm prints "entropies"). They hold
    // cp_l*T/s_0 — a NORMALISED TEMPERATURE, dimensionless. Dry static energy is cp*T + g*z
    // and the g*z is absent: on this 120 km shell it reaches 1.2e6 J/kg against cp*T = 3.4e5
    // at the cold top, so a parcel conserving s does not cool as it rises. ATHAD README item 52.
    Array s;                                                            // cp_l*T/s_0, environment
    Array s_u;                                                          // cp_l*T/s_0, updraft
    Array s_d;                                                          // cp_l*T/s_0, downdraft
    Array u_u;                                                          // u-velocity component in the updraft
    Array u_d;                                                          // u-velocity component in the downdraft
    Array v_u;                                                          // u-velocity component in the updraft
    Array v_d;                                                          // u-velocity component in the downdraft
    Array w_u;                                                          // u-velocity component in the updraft
    Array w_d;                                                          // u-velocity component in the downdraft
    Array q_v_u;                                                        // water vapour in the updraft
    Array q_v_d;                                                        // water vapour in the downdraft
    Array q_c_u;                                                        // cloud water in the updraft
    Array E_u;                                                          // moist entrainment in the updraft
    Array D_u;                                                          // moist detrainment in the updraft
    Array E_d;                                                          // moist entrainment in the downdraft
    Array D_d;                                                          // moist detrainment in the downdraft

    Array CloudBase;                                                    // cloud base
    Array LevelFreeSinking;                                             // level of free sinking

    Array Deep_beg;                                                     // cloud base
    Array Deep_end;                                                     // level of free sinking

    // Turbulence arrays
    Array tke;                                                          // turbulent kinetic energy
    Array tken;                                                         // turbulent kinetic energy (new time level)
    Array dis;                                                          // dissipation rate (ε or specific dissipation ω)
    Array disn;                                                         // dissipation rate (new time level)
    Array nue;                                                          // turbulent eddy viscosity
    Array prod;                                                         // turbulence production P_k
    Array tke_source;                                                   // RHS source term for k equation
    Array dis_source;                                                   // RHS source term for ε/ω equation
    Array rhs_tke;                                                      // auxiliary RHS field for tke
    Array rhs_dis;                                                      // auxiliary RHS field for dis
    Array p_hydro;                                                      // hydrostatic pressure
    Array PressureGradientForce;                                        // pressure gradient force magnitude
};
#endif

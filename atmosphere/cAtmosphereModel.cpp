/*
 * Atmosphere General Circulation Modell(AGCM) applied to laminar flow
 * program for the computation of geo-atmospherical circulating flows in a spherical shell
 * including the magnetic field of the Earth
 * finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 1 additional transport equations to describe the salinity
 * 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop
 * Poisson equation for the pressure solution in an outer iterational loop
 * temperature distribution given as a parabolic distribution from pole to pole, zonaly constant
 * code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz(roger.grundmann@web.de)
*/


#include "cAtmosphereModel.h"

#include <cstdlib>   // getenv/atof for the ATM_RADIAL_SHAPIRO_STRENGTH A/B knob

#include "MoistConvection.h"

#include "ZeroCatIceScheme.h"
#include "OneCatIceScheme.h"
#include "TwoCatIceScheme.h"
#include "ThreeCatIceScheme.h"
#include "SaturationAdjustment.h"
#include "VelocityInitializer.h"
#include "PressureSolverAtm.h"
#include "ThermoAtm.h"
#include "UtilsAtm.h"
#include "BC_Atm.h"
#include "TurbulenceAtm.h"
#include "MultiLayerRadiation.h"


using namespace std;
using namespace tinyxml2;
using namespace AtomUtils;


extern std::vector<std::vector<double> > m_node_weights;


cAtmosphereModel* cAtmosphereModel::m_model = NULL;

const double cAtmosphereModel::pi180 = 180.0/M_PI;                      // pi180 = 57.3

const double cAtmosphereModel::the_degree = 1.0;                        // compares to 1° step size laterally
const double cAtmosphereModel::phi_degree = 1.0;                        // compares to 1° step size longitudinally

// The radial step must span exactly 1.0 over the im levels, so that rad.z runs 1.0 .. 2.0
// and metricShellLength() reads one rad.z unit as the shell thickness.
//
// This was a hard-coded 0.025, silently tied to im = 41 (0.025 x 40 = 1.0). The two
// constants were independent and consistent only by coincidence: raising im to 61 left
// dr at 0.025, so rad.z ran 1.0 .. 2.5, one rad.z unit was read as 932 km instead of
// 300 km, and the domain top landed at 1399 km with the whole column clamped to t_00.
// Deriving it from im removes the coupling.
const double cAtmosphereModel::dr = 1.0/(double)(cAtmosphereModel::im - 1);

const double cAtmosphereModel::dthe = the_degree/pi180; 
const double cAtmosphereModel::dphi = phi_degree/pi180;
    
const double cAtmosphereModel::the0 = 0.0;                              // North Pole
const double cAtmosphereModel::phi0 = 0.0;                              // zero meridian in Greenwich

//earth's radius is r_earth = 6731 km, here it is assumed to be infinity, circumference of the earth 40074 km 
const double cAtmosphereModel::r0 = 1.0;                                // non-dimensional
//const double cAtmosphereModel::r_Earth = 6731.001;                      // in km

const double cAtmosphereModel::residuum_ref_atm = 1.0e-2;               // criterium to finish global iterations




cAtmosphereModel::cAtmosphereModel():
    i_topography(std::vector<std::vector<int> >(jm, std::vector<int>(km, 0))),
    i_tropopause(std::vector<std::vector<int> >(jm, std::vector<int>(km, 0))),
    i_landscape(std::vector<std::vector<int> >(jm, std::vector<int>(km, 0))),

    has_printed_welcome_msg(false){

    // If Ctrl-C is pressed, quit
    signal(SIGINT, exit);

    // set default configuration
    SetDefaultConfig();

    m_model = this;

    //  Coordinate system in form of a spherical shell
    //  rad for r-direction normal to the surface of the earth, the for lateral and phi for longitudinal direction
    rad.initArray_1D(im, 0);                                            // radial coordinate direction
    the.initArray_1D(jm, 0);                                            // lateral coordinate direction
    phi.initArray_1D(km, 0);                                            // longitudinal coordinate direction
}

cAtmosphereModel::~cAtmosphereModel(){
    m_model = NULL;
    logger().close();
}
 
#include "cAtmosphereDefaults.cpp.inc"
/*
*
*/
void cAtmosphereModel::LoadConfig(const char *filename){
    XMLDocument doc;
    XMLError err = doc.LoadFile(filename);
    if(err){
        doc.PrintError();
        throw std::invalid_argument("   couldn't load config file inside cAtmosphereModel");
    }
    XMLElement *atom = doc.FirstChildElement("atom");
    if(!atom){
        return;
    }
    XMLElement* elem_common = doc.FirstChildElement("atom")->FirstChildElement("common");
    if(!elem_common){
        return;
    }
    XMLElement* elem_atmosphere = doc.FirstChildElement("atom")->FirstChildElement("atmosphere");
    if(!elem_atmosphere){
        return;
    }
#include "AtmosphereLoadConfig.cpp.inc"

    initComposition();
}
/*
*
*/
// Resolve the configured mole fractions into the mass fractions, molar masses and gas
// constants the physics actually uses, and report them. Called once, at the end of
// LoadConfig, so every later consumer sees a populated m_comp.
void cAtmosphereModel::initComposition(){

    m_comp = AtmMixture::resolve(x_H2O, x_CO2, x_N2, x_CH4,
                                 x_NH3, x_H2,  x_CO, x_SO2);

    cout << endl << endl << "      AGCM: atmospheric composition (ATHAD)" << endl << endl;
    cout.precision(6);
    cout << "        species    mole frac.    mass frac." << endl;
    auto line = [&](const char* n, double x, double q){
        cout << "        " << setw(8) << left << n << resetiosflags(ios::left)
             << setw(12) << fixed << x << setw(14) << q << endl;
    };
    line("H2O", x_H2O, m_comp.q_H2O);
    line("CO2", x_CO2, m_comp.q_CO2);
    line("bkgnd", 1.0 - x_H2O - x_CO2, m_comp.q_bg);
    cout << endl;
    cout << "        mean molar mass ......... M_mean = " << m_comp.M_mean * 1.0e3 << " g/mol" << endl;
    cout << "        mixture gas constant .... R_mix  = " << m_comp.R_mix  << " J/(kg K)" << endl;
    cout << "        background molar mass ... M_bg   = " << m_comp.M_bg * 1.0e3 << " g/mol" << endl;
    cout << "        background gas constant . R_bg   = " << m_comp.R_bg   << " J/(kg K)" << endl;

    // The mole fractions are a hand-entered table; a typo that makes them sum to
    // something other than 1 would silently rescale every derived property.
    if(!m_comp.valid){
        cout << endl << "      ERROR: mole fractions sum to " << m_comp.x_sum
             << ", not 1. Fix x_* in param.py." << endl;
        throw std::invalid_argument("   ATHAD: atmospheric mole fractions do not sum to 1");
    }

    // R_Air is the CONFIGURED background gas constant; m_comp.R_bg is the one implied by
    // the mole fractions. They are independent inputs, so a composition edit that forgets
    // R_Air would leave the two disagreeing and the barometric profile subtly wrong.
    if(std::fabs(R_Air - m_comp.R_bg) > 1.0){
        cout << endl << "      WARNING: R_Air = " << R_Air
             << " J/(kg K) disagrees with the background gas constant implied by the mole "
             << "fractions, R_bg = " << m_comp.R_bg << " J/(kg K)." << endl;
    }

    // The reference density must reproduce the intended surface pressure, since
    // p_stat.x[0] is built as 1e-2*(r_air*R_mix*T) and NOT from p_0. R_mix, not R_Air:
    // r_air was calibrated as p/(R_mix*T), so the background constant would give 204 bar.
    const double p_surf_implied = 1.0e-2 * r_air * m_comp.R_mix * t_surf_equator;
    cout << "        implied equatorial surface pressure = " << p_surf_implied * 1.0e-3
         << " bar, p_0 = " << p_0 * 1.0e-3 << " bar" << endl;
    if(std::fabs(p_surf_implied - p_0) > 0.01 * p_0){
        cout << "      WARNING: r_air is inconsistent with p_0 by more than 1%. "
             << "Set r_air = p_0*1e2/(R_mix*t_surf_equator) = "
             << p_0 * 1.0e2 / (m_comp.R_mix * t_surf_equator) << " kg/m³." << endl;
    }

    // COSMO barometric lapse parameter. See param.py (cosmo_lapse_fraction) for why this
    // is derived and not the inherited constant 42 K.
    m_beta_cosmo = cosmo_lapse_fraction * m_comp.R_mix * t_surf_equator / cp_l;

    const double lapse_K_per_km = m_beta_cosmo * g / (m_comp.R_mix * t_surf_equator) * 1.0e3;
    const double adiabat_K_per_km = g / cp_l * 1.0e3;

    cout << "        COSMO lapse parameter ... beta   = " << m_beta_cosmo << " K" << endl;
    cout << "        near-surface lapse rate .......... = " << lapse_K_per_km
         << " K/km   (dry adiabat " << adiabat_K_per_km << " K/km)" << endl;

    // Skin temperature of the optically thin top, DERIVED from the energy budget rather
    // than prescribed.
    //
    // The top of the column is isothermal and radiative, and in equilibrium it must emit
    // what the planet absorbs: sigma*T_skin^4 = SW_absorbed + geothermal. Setting t_skin
    // by hand is circular — the first attempt took T_skin = (OLR/2sigma)^(1/4) from a
    // previously MEASURED OLR, and the model then dutifully reproduced an OLR of
    // sigma*T_skin^4, which says nothing except that the arithmetic is consistent.
    //
    // This is a one-shot estimate, not an iterated fixed point: the albedo it uses is the
    // clear-sky (molten-surface) value, whereas the cloud deck that forms above the
    // condensation level will raise it and lower the true balance. Closing that loop
    // properly means iterating t_skin against the model's own albedo. Stated in the README.
    {
        const double sw_mean = 0.5 * (rad_equator_short + rad_pole_short);   // [W/m2] TOA mean
        const double alb_clear = 0.08;                                       // molten surface, clear sky
        const double absorbed  = (1.0 - alb_clear) * sw_mean + geothermal_flux;
        const double t_skin_eq = std::pow(absorbed / sigma, 0.25);

        cout << "        absorbed SW + geothermal ......... = " << absorbed << " W/m2" << endl;
        cout << "        implied skin temperature ......... = " << t_skin_eq << " K"
             << "   (configured t_skin = " << t_skin << " K)" << endl;
        if(std::fabs(t_skin_eq - t_skin) > 5.0)
            cout << "      NOTE: configured t_skin differs from the energy-balance value by "
                 << std::fabs(t_skin_eq - t_skin) << " K. The top of the column will emit at the "
                 << "configured value, so the OLR will follow it rather than the budget." << endl;
    }

    cout << "        rotation rate ... omega = " << omega << " rad/s  ("
         << (2.0 * M_PI / omega) / 3600.0 << " h day, "
         << omega / 7.292e-5 << "x modern)" << endl;

    cout << "        text diagnostics every " << diagnosticStride() << " iterations"
         << (diagnostic_stride > 0 ? "  (set explicitly)" : "  (auto: nm = ")
         << (diagnostic_stride > 0 ? "" : std::to_string(nm) + ")") << endl;

    cout << endl << "      AGCM: atmospheric composition ended" << endl << endl;
}
/*
*
*/
// Cos(latitude)-weighted shortwave budget over the sphere, using the model's OWN albedo.
//
// The weight is cos(latitude) with latitude taken from j/(jm-1), which is the same weight
// AtomUtils::GetMean_2D uses — written here in a form that does not assume jm == 181, as
// the node-weight table does.
bool cAtmosphereModel::planetaryShortWave(double& albedo_mean, double& sw_mean,
                                          double& absorbed_mean) const
{
    albedo_mean = sw_mean = absorbed_mean = 0.0;
    if((int)short_wave_radiation.size() != jm || jm < 2 || km < 1) return false;

    double w_sum = 0.0, alb_w = 0.0, sw_w = 0.0, abs_w = 0.0, alb_max = 0.0;
    for(int j = 0; j < jm; j++){
        const double w = std::cos((j / (double)(jm - 1) - 0.5) * M_PI);
        double alb_k = 0.0;
        for(int k = 0; k < km; k++) alb_k += albedo.y[j][k];
        alb_k /= (double)km;
        alb_max = std::max(alb_max, alb_k);

        w_sum += w;
        alb_w += w * alb_k;
        sw_w  += w * short_wave_radiation[j];
        abs_w += w * (1.0 - alb_k) * short_wave_radiation[j];
    }
    if(!(w_sum > 0.0)) return false;

    albedo_mean   = alb_w / w_sum;
    sw_mean       = sw_w  / w_sum;
    absorbed_mean = abs_w / w_sum;

    // Before the first MultiLayerRadiation call albedo.y is all zeros, which would report a
    // perfectly black planet rather than "not measured yet".
    return alb_max > 0.0;
}
/*
*
*/
// One relaxation step of the t_skin fixed point.
//
// t_skin is the temperature of the isothermal top that ThermoAtm::densities() imposes
// wherever the adiabat falls below it. That top is optically thick, so the emission level
// sits there and the OLR the model reports is identically sigma*t_skin^4 — verified by
// running the model at 254 K and at 240 K and getting 236.01 and 188.13 W/m2, against
// sigma*T^4 = 236.01 and 188.13. Prescribing t_skin therefore prescribes the OLR.
//
// This does not repair that. What it repairs is the CONSISTENCY of the prescription: the
// configured 254 K came from a one-shot estimate using the clear-sky albedo (0.08), while
// the model's own albedo is ~0.50, so the planet was emitting 236 W/m2 while absorbing
// 204 — a 32 W/m2 imbalance that no amount of running would close, because nothing in the
// model could move t_skin. Iterating it against the model's own albedo closes the budget.
// The budget is then closed BY CONSTRUCTION and is not a test of anything; what the model
// still genuinely computes is the albedo, the optical depth and where water condenses.
//
// Making the OLR a real output means letting the top find its own temperature radiatively
// instead of having densities() re-impose it every iteration. See the README.
void cAtmosphereModel::updateSkinTemperature(bool report)
{
    if(t_skin_relax <= 0.0) return;

    double alb_mean = 0.0, sw_mean = 0.0, abs_mean = 0.0;
    if(!planetaryShortWave(alb_mean, sw_mean, abs_mean)) return;

    const double in_mean = abs_mean + geothermal_flux;
    if(!(in_mean > 0.0)) return;

    const double target = std::pow(in_mean / sigma, 0.25);
    const double t_old  = t_skin;
    t_skin += std::min(1.0, t_skin_relax) * (target - t_skin);

    if(report)
        cout << "      AGCM: t_skin " << std::fixed << std::setprecision(2) << t_old
             << " -> " << t_skin << " K   (fixed-point target " << target
             << " K from albedo " << std::setprecision(4) << alb_mean
             << ", absorbed + geothermal " << std::setprecision(2) << in_mean << " W/m2)"
             << endl;
}
/*
*
*/
// THE definition of the grid coordinates — see the declaration in cAtmosphereModel.h.
//
// This used to be written in three places with two different meanings, which is a trap rather
// than a redundancy:
//   cAtmosphereModel::RunTimeSlice   rad.Coordinates(im, r0, dr)      the real coordinate
//   UtilsAtm::resetArrays            rad.initArray_1D(im, 1.0)        a debug sentinel
//   MoistConvection::precompute      rad.Coordinates(im, r0, dr)      the real coordinate, again
// The order saved it: resetArrays runs before RunTimeSlice sets the coordinate, so the sentinel
// was always overwritten. But the sentinel value 1.0 is exactly r0, so had the order ever changed
// the metric would have collapsed to r = 1 everywhere and still looked plausible. And the
// MoistConvection copy re-imposes r0/dr on every call, so a change made at the RunTimeSlice site
// alone would have been silently undone once per iteration.
void cAtmosphereModel::initGridCoordinates(){
    rad.Coordinates(im, r0, dr);
    the.Coordinates(jm, the0, dthe);
    phi.Coordinates(km, phi0, dphi);
}
/*
*
*/
// Does the radius the CORE metric uses equal the radius the PHYSICS uses?
//
// Right now it does not, and this reports by how much rather than aborting, because the mismatch
// is the model's present state and not a fresh regression:
//
//   The core. RungeKutta_Atm_Turb sets geo.rm = rad.z[i] and inv_rm = 1/rm, and RHS_Atm_Turb uses
//   that as the PLANETARY radius in v/r d/dthe, 1/(r sinthe) d/dphi and in every Laplacian. With
//   r0 = 1.0 and dr = 0.025 over im = 41, rad.z runs 1.0 .. 2.0.
//
//   The physics. ThermoAtm (r_Earth_m, rad_Earth), MoistConvection ((r_Earth*1e3 + height)*costhe)
//   and MinMax_Atm (a = r_Earth*1000) all use the configured r_Earth = 6370.001 km.
//
//   The scaling. Every nondimensional coefficient in RHS_Atm_Turb — coeff_u_p, coeff_MC_*, coeff_S,
//   coeff_L, nue_max, force_nd — divides by L_atm, the 400 m RADIAL GRID STEP.
//
// In units of L_atm the Earth's radius is 6.370e6/400 = 15925, and the core puts 1.0 there. So the
// horizontal metric terms are ~16000x too large, and force_nd = omega*L_atm/u_0 is the same factor
// too small against the length-consistent omega*r_Earth/u_0 ~ 58. That is the arithmetic behind the
// note at RHS_Atm_Turb.cpp and behind ATM_CORIOLIS_SCALE having existed at all (now removed).
//
// Set ATM_METRIC_STRICT=1 to make the mismatch fatal instead of merely loud — for use AFTER the
// metric is made consistent, so it cannot drift back apart unnoticed.
// ATM_METRIC_RADIUS=<km> — opt-in, default off (bit-identical when unset).
//
// Gives the HORIZONTAL metric the planetary radius it is supposed to have. What it does NOT do is
// move rad.z, and that restraint is the whole design: rad.z is not only a radius, it is also the
// stretched radial coordinate that exp_rm = 1/(rm+1) and init_layer_heights are built on. Shifting
// it to 397.5 would divide every radial derivative by ~200 and send exp(zeta*(rad.z-1)) to
// overflow. Only the 1/r factors move; exp_rm stays a function of the grid coordinate.
//
// The value is the radius in km, so ATM_METRIC_RADIUS=6370 is the Earth. Converted to rad.z units
// with metricShellLength(), which is ~16 km per unit, so r0 lands at ~397.5 instead of 1.0.
// DEFAULT SINCE 2026-07-28: on, at the configured r_Earth. Unset means "use the real radius";
// ATM_METRIC_RADIUS=0 restores the old rad.z ~ 1.5 metric for A/B work.
//
// The evidence for flipping it, all measured at 20 iterations against the old default:
//   the radial wind stops being pumped   rms(u) 0.01816 -> 0.00124, p_dyn max 0.2452 -> 0.00108
//   the zonal wind is untouched          rms(w) 0.848 -> 0.858, temperature within 0.3 %
//   the eddy viscosity becomes physical  mean nu_T 1.45 -> 6.33 m2/s, max 22 -> 279 m2/s, and the
//                                        nue_max cap (1000 m2/s) goes from 2.2 % to 27.9 % used,
//                                        so it still does not bind. k is unchanged (-3 %); what
//                                        changed is omega, which fell 7.4x once the ~400x
//                                        inflated horizontal shear stopped feeding its production
//                                        and stopped holding the SST limiter down (median
//                                        nue/(k/omega) 0.026 -> 0.144).
// A mean eddy viscosity of 6 m2/s with boundary-layer peaks near 280 is the right order for the
// atmosphere; 1.45 was too small. The rise is a suppression being removed, not energy appearing.
void cAtmosphereModel::initMetricRadius(){
    static const double r_km = [this](){
        const char* e = getenv("ATM_METRIC_RADIUS");
        return e ? atof(e) : r_Earth; }();
    if(r_km <= 0.0){
        m_metric_r0 = 0.0;
        cout << "      AGCM: ATM_METRIC_RADIUS = 0 - horizontal metric left on the grid"
             << " coordinate (pre-2026-07-28 behaviour)" << endl;
        return;
    }

    const double L_unit = metricShellLength();
    m_metric_r0 = (L_unit > 0.0) ? (r_km * 1.0e3 / L_unit) : 0.0;

    cout << "      AGCM: ATM_METRIC_RADIUS = " << r_km << " km, one rad.z unit = "
         << L_unit << " m  ->  metric r0 = " << m_metric_r0
         << " (grid r0 stays " << rad.z[0] << ")" << endl;
}
/*
*
*/
void cAtmosphereModel::checkMetricConsistency() const {
    // ONE rad.z unit in metres — the shell thickness (~16 km), NOT L_atm. See metricShellLength().
    // An earlier version of this check used L_atm and therefore reported 15925 instead of 397.5;
    // the two readings differ by exactly 1/dr = 40.
    const double L_unit = metricShellLength();
    const double r_metric_m = metricRadius(rad.z[0]) * L_unit;   // radius the core metric implies [m]
    const double r_planet_m = r_Earth * 1.0e3;                   // radius the physics uses        [m]
    const double r0_consistent = r_planet_m / L_unit;            // what the metric r0 would have to be
    const double ratio = (r_metric_m > 0.0) ? r_planet_m / r_metric_m : 0.0;

    // 1 per mille is well inside anything a deliberate choice would produce.
    const bool consistent = std::fabs(ratio - 1.0) < 1.0e-3;

    cout << endl
         << "      AGCM: metric check - core radius = " << r_metric_m << " m"
         << " (metric r0 " << metricRadius(rad.z[0]) << " x " << L_unit << " m per rad.z unit)"
         << ", physics radius r_Earth = " << r_planet_m << " m" << endl;

    if(consistent){
        cout << "      AGCM: metric check - consistent." << endl;
        return;
    }

    cout << "      AGCM: metric check - MISMATCH by a factor of " << ratio << "." << endl
         << "            The core metric (inv_rm in RHS_Atm_Turb) and the physics modules"
         << " (ThermoAtm, MoistConvection, MinMax_Atm) are using different planetary radii." << endl
         << "            Horizontal metric terms are too large and force_nd too small by that"
         << " factor." << endl
         << "            The metric r0 would have to be " << r0_consistent
         << " instead of " << metricRadius(rad.z[0])
         << " - set ATM_METRIC_RADIUS=" << (r_planet_m / 1.0e3) << " to do that." << endl;

    static const bool strict = [](){
        const char* e = getenv("ATM_METRIC_STRICT"); return e && atoi(e) != 0; }();
    if(strict){
        cout << "      AGCM: metric check - ATM_METRIC_STRICT is set, stopping." << endl;
        std::exit(1);
    }
}
// Is exp_rm the Jacobian of the radial stretch? It is documented as one at TurbulenceAtm.h
// ("exp_rm = 1/(rm+1) is the Jacobian of the radial coordinate transformation") and it is not.
// ATHAD README item 39 has the measurement; this is the check that keeps it visible. Ported
// 2026-08-18 -- it was live here, unmeasured, at the same zeta = 3.0.
//
// THE TEST IS UNIT-FREE, which is what makes it conclusive. Every radial derivative in the core
// is (index difference) * inv_2dr * exp_rm. For that to be d/dz_physical under ONE constant
// length normalisation — whatever that normalisation happens to be — the quantity
//
//     ratio(i) = [inv_2dr * exp_rm(i)] / [1 / (z[i+1] - z[i-1])]
//
// must be CONSTANT in i. Its spread is therefore not a units question and cannot be argued away
// by picking a different L_unit: any variation is a wrong metric SHAPE.
//
// exp_rm = 1/(rm+1) implies dz/d(rm) = rm+1, i.e. a QUADRATIC stretch z ~ (rm+1)^2/2, while
// init_layer_heights builds an EXPONENTIAL one, z = (exp(zeta*(rm-rad.z[0])) - 1)*L_atm. Over
// rad.z in [1,2] the quadratic form spans 1.5x and the exponential spans exp(zeta), so the two
// agree only near zeta = ln(1.5) = 0.405. At the shipped zeta = 3.0 the spread is 11.8x.
//
// It is NOT a function of im: 11.8x at im = 41 and 12.3x at im = 61, so refining the grid does
// not help and marginally hurts. zeta is the only lever short of replacing exp_rm outright,
// which touches ~91 sites across 10 files including the Poisson operator.
void cAtmosphereModel::checkRadialMetric() const {
    if(m_layer_heights.size() < (std::size_t)im) return;      // called before init_layer_heights

    static const bool verbose = [](){
        const char* e = getenv("ATM_METRIC_CHECK"); return e && atoi(e) != 0; }();

    // setprecision/fixed are sticky on the stream, so save and restore rather than leaving
    // every later diagnostic in whatever format this one wanted.
    const std::ios::fmtflags saved_flags = cout.flags();
    const std::streamsize    saved_prec  = cout.precision();

    const double inv_2dr = 1.0 / (2.0 * dr);
    double ratio_min = 0.0, ratio_max = 0.0;
    bool first = true;

    if(verbose){
        cout << endl
             << "      AGCM: radial metric check - is exp_rm the Jacobian of the stretch?" << endl
             << "         i    rad.z    exp_rm     z[km]  dz_true[km]    J_true[m]  J_code[-]"
             << "     ratio[m]" << endl;
    }

    for(int i = 1; i < im-1; i++){
        const double rm     = rad.z[i];
        const double exp_rm = 1.0 / (rm + 1.0);
        const double dz2    = m_layer_heights[i+1] - m_layer_heights[i-1];       // [m]
        if(!(dz2 > 0.0)) continue;

        // ratio has units of metres: the length unit the core is implicitly working in at
        // this level. Constant iff exp_rm is the true Jacobian up to that one unit.
        const double ratio = (inv_2dr * exp_rm) * dz2;

        if(first){ ratio_min = ratio_max = ratio; first = false; }
        else { ratio_min = std::min(ratio_min, ratio); ratio_max = std::max(ratio_max, ratio); }

        if(verbose && (i % 4 == 1 || i == im-2)){
            const double J_true = zeta * L_atm * exp(zeta * (rm - rad.z[0]));    // dz/d(rad.z) [m]
            cout << "      " << setw(4) << i
                 << setw(9)  << setprecision(4) << fixed << rm
                 << setw(10) << setprecision(5) << exp_rm
                 << setw(10) << setprecision(1) << m_layer_heights[i] / 1000.0
                 << setw(13) << setprecision(3) << (m_layer_heights[i+1] - m_layer_heights[i]) / 1000.0
                 << setw(13) << setprecision(0) << J_true
                 << setw(11) << setprecision(3) << (1.0 / exp_rm)
                 << setw(13) << setprecision(1) << ratio << endl;
        }
    }

    if(first){ cout.flags(saved_flags); cout.precision(saved_prec); return; }
    const double spread = (ratio_min > 0.0) ? ratio_max / ratio_min : 0.0;

    // ln(1.5) is where the quadratic and exponential stretches span the same range over
    // rad.z in [1,2] — the zeta at which exp_rm would be right by construction.
    const double zeta_consistent = log(1.5);

    cout << "      AGCM: radial metric check - core length unit runs "
         << setprecision(0) << ratio_min << " m at the surface to "
         << ratio_max << " m at the top, spread " << setprecision(2) << spread
         << "x (1.00 = exp_rm is the Jacobian)." << endl;

    if(spread > 1.05){
        cout << "            NOT the Jacobian: exp_rm = 1/(rm+1) is a QUADRATIC stretch, the grid"
             << " is EXPONENTIAL (zeta = " << setprecision(3) << zeta << ")." << endl
             << "            Every radial derivative is mis-scaled by a factor varying "
             << setprecision(2) << spread << "x across the column, low at the surface."
             << endl
             << "            No length normalisation can absorb it. zeta = "
             << setprecision(3) << zeta_consistent
             << " would make exp_rm correct by construction; im does not move it." << endl
             << "            See ATHAD README item 39." << (verbose ? "" :
                " ATM_METRIC_CHECK=1 for the per-level table.") << endl;
    }

    cout.flags(saved_flags);
    cout.precision(saved_prec);
}
/*
*
*/
/*
*
*/
void cAtmosphereModel::RunTimeSlice(int Ma){

    #ifdef _OPENMP
        printf("\n\n   number of processors: %d\n\n", omp_get_num_procs());

    #pragma omp parallel
        {
            printf("   thread %d of %d in \"Desktop-Dell-XPS 8960\"\n", 
                omp_get_thread_num(), omp_get_num_threads());
        }
    #else
        printf("   OpenMP is not supported\n");
    #endif
        printf("   ended\n\n");




    cout << endl << endl << endl 
        << " ----------------------------- Ma = "  << Ma 
        << " million years back in time ------------------------------ " 
        << endl << endl;

    int start = AtomUtils::RunStart(at);

    m_current_time = m_time_list.insert(int(Ma)).first;

    use_presets(use_earthbyte_reconstruction, 
        use_NASA_temperature, use_NASA_velocity);

    UtilsAtm(*this).resetArrays();

    dt = dt_visc;                                                       //  no dimension; toggled to dt_inviscid during spin-up

    iter_n = 0;
    total_iter_count = 0;                                               // reset per time slice so the inviscid spin-up fires at the start of every Ma slice

    // "laminar" is the descriptive name for the no-turbulence path; accept the legacy
    // "none" and an empty string as synonyms so old configs keep working.
    if(turb_model.empty() || turb_model == "none") turb_model = "laminar";

    use_turbulence_model             = (turb_model != "laminar");
    use_k_epsilon_turbulence_model   = (turb_model == "k_epsilon");
    use_k_omega_turbulence_model     = (turb_model == "k_omega");
    use_k_omega_SST_turbulence_model = (turb_model == "k_omega_SST");

    initGridCoordinates();
    initMetricRadius();
    checkMetricConsistency();

    #pragma omp parallel sections
    {
        #pragma omp section
        { init_layer_heights(); }

        #pragma omp section
        { init_tropopause_layers(); }
    }

    // AFTER init_layer_heights: it compares exp_rm against the heights that call builds.
    checkRadialMetric();

    read_Atmosphere_Surface_Data(Ma);                                   // reading topography data and NASA measurements

    LandOceanFraction();
                                                // ratio of land to sea surface of the various time slices

    if(!use_NASA_velocity){
        VelocityInitializer(*this).compute();                           // construction of zonal initial velocities from measurements
    }

    // Skip the boundary-layer damping and the wiggle filter on the initial velocity
    // field when an inviscid spin-up is requested: those routines suppress exactly the
    // large-scale near-surface gradients that the Euler phase needs in order to organise
    // a global circulation. The viscous phase that follows will re-apply both.
//    if(inviscid_spinup_iters <= 0){
        BC_Atm(*this).bcVelSurfSur();                                   // velocities close to surfaces, resembling a boundary layer

        AtomUtils::damp_wiggles(u, &i_topography, true, true, true);
        AtomUtils::damp_wiggles(v, &i_topography, true, true, true);
        AtomUtils::damp_wiggles(w, &i_topography, true, true, true);
//    }

//    BC_Atm(*this).coastalCurrents();

    initTemperatureData(Ma);                                            // initialization of temperature, hydrostatic pressure and density of dry air, reconstruction of potential surface values
    AtomUtils::damp_wiggles(t, &i_topography, true, true, true);

    // Snapshot the lid temperature (i=im-1) from the IC. bcRadius pins t at the lid
    // to this fixed reference so the isothermal stratospheric top stays constant
    // instead of drifting up via the former cubic top extrapolation.
    t_top_init.assign(jm, std::vector<double>(km, 0.0));
    for(int j = 0; j < jm; j++)
        for(int k = 0; k < km; k++)
            t_top_init[j][k] = t.x[im-1][j][k];

    // ATHAD_COND: the initial state is a FIXED POINT of three fields, not a sequence.
    //
    // In ATHAD each of these could be built once, in order, because none of them fed back:
    // the water field was uniform, the CO2 field was uniform, and the profile depended on
    // them only through constants. Here every arrow is live —
    //
    //     c   depends on T and p   (it is the saturation profile, initWaterWapour)
    //     co2 depends on c         (well mixed in the DRY air, co2Atmosphere)
    //     T,p depend on c and co2  (through R_mix, cp and the moist adiabat, densities)
    //
    // — so building them once, on the provisional profile initTemperatureData leaves, gives
    // a column that is nowhere saturated: the stored c belongs to the old temperatures and
    // sits ~20 % below the q_sat of the new ones. The moist adiabat in densities() then
    // never fires, the profile stays dry-adiabatic, and the model silently runs the physics
    // it was forked to stop running. That is how this was found.
    //
    // Iterating c and (T,p) against each other does NOT converge — see the long note in
    // ThermoAtm::densities(): it is a fixed point over the runaway-greenhouse feedback, and
    // it climbed 17 K between three passes and eight without settling. densities(true)
    // integrates the water profile in the same upward sweep as the temperature and the
    // pressure instead, which is one pass and exact.
    //
    // initWaterWapour still runs first, to lay down a starting c for the surface anchor and
    // for the levels the sweep does not diagnose; co2Atmosphere then follows c. Only the CO2
    // dependence is left iterating, and that one is weak (co2 enters through R and cp, not
    // through the latent heat), so two passes settle it to the fourth digit.
    initWaterWapour();
    for (int pass = 0; pass < 2; pass++) {
        ThermoAtm(*this).co2Atmosphere();                               // co2 well mixed in the dry air
        ThermoAtm(*this).densities(true);                               // T, p AND the water profile
    }
    initCloudIce();
//    init_vapour_cloud();                                                // initialisation of water vapour and cloud/ice formation based on the temperature profile

//    goto Printout;

    SaturationAdjustment(*this).run();                                  // based on the initial distribution, recomputation of the cloud water and cloud ice formation in case of saturated water vapour detected
    AtomUtils::damp_wiggles(ice,   &i_topography, true, true, true);
    AtomUtils::damp_wiggles(c,     &i_topography, true, true, true);
    AtomUtils::damp_wiggles(cloud, &i_topography, true, true, true);

    ThermoAtm(*this).precipitableWater();

    switch(CategoryIceScheme){                                          // rain, snow graupel and precipitation production and reduction
        case -1: cout << endl << endl << endl                           // no CategoryIceScheme used
            << "  no CategoryIceScheme used" << endl;
                break;
        case 0: ZeroCatIceScheme(*this).run();                          // development of rain fall, water vapour and cloud water
                break;
        case 1: OneCatIceScheme(*this).run();                           // development of rain and snow fall, water vapour, cloud water and ice
                break;
        case 2: TwoCatIceScheme(*this).run();                           // development of rain and snow fall, water vapour, cloud water and ice
                break;
        case 3: ThreeCatIceScheme(*this).run();                         // development of rain and snow fall, water vapour, cloud water, ice and graupel
                break;
    }
    AtomUtils::damp_wiggles(P_rain, &i_topography, true, true, true);
    AtomUtils::damp_wiggles(P_snow, &i_topography, true, true, true);

//    goto Printout;

    MoistConvection(*this).run(iter_n);                                 // rainfall from convecting clouds

    AtomUtils::damp_wiggles(P_conv, &i_topography, true, true, true);
    AtomUtils::damp_wiggles(E_u,    &i_topography, true, true, true);
    AtomUtils::damp_wiggles(M_u,    &i_topography, true, true, true);
    AtomUtils::damp_wiggles(q_v_u,  &i_topography, true, true, true);
    AtomUtils::damp_wiggles(q_c_u,  &i_topography, true, true, true);
    AtomUtils::damp_wiggles(E_d,    &i_topography, true, true, true);
    AtomUtils::damp_wiggles(M_d,    &i_topography, true, true, true);
    AtomUtils::damp_wiggles(q_v_d,  &i_topography, true, true, true);

//    goto Printout;

    if(turb_model != "laminar") {
        TurbulenceAtm(*this).init();                                    // turbulence field initialisation (skipped in laminar mode)
        TurbulenceAtm(*this).run();                                     // call for turbulence models

        AtomUtils::damp_wiggles(tke,        &i_topography, true, true, true);
        AtomUtils::damp_wiggles(dis,        &i_topography, true, true, true);
        AtomUtils::damp_wiggles(nue,        &i_topography, true, true, true);
        AtomUtils::damp_wiggles(prod,       &i_topography, true, true, true);
        AtomUtils::damp_wiggles(tke_source, &i_topography, true, true, true);
        AtomUtils::damp_wiggles(dis_source, &i_topography, true, true, true);

    }

    UtilsAtm(*this) .precipitationSum();
    // ATHAD: co2Atmosphere() must run BEFORE densities(). The CO2 field is now a
    // mass fraction that enters the local mixture gas constant R_of(c, co2), so a
    // density built before it is set uses R_of(c, 0) = 414.2 instead of 387.9 —
    // a 7 % density error through the whole column. On Earth the ordering was
    // harmless because co2 was in ppm and never touched the density.
    ThermoAtm(*this).co2Atmosphere();                                   // INITIAL well-mixed CO2 mass fraction
    ThermoAtm(*this).co2Column(true);                                   // column CO2 path + the conservation reference
    ThermoAtm(*this).waterBudget(true);                                 // total-water conservation reference
    ThermoAtm(*this).densities();
    ThermoAtm(*this).forces();
    ThermoAtm(*this).standAtm_DewPoint_HumidRel();                      // International Standard Atmosphere temperature profile, dew point temperature, relative humidity profile
    ThermoAtm(*this).waterVapourEvaporation();                          // correction of surface water vapour by evaporation
    ThermoAtm(*this).latentSensibleHeat();                              // latent and sensible heat
    ThermoAtm(*this).vegetationLand();                                  // vegetation on land

    {   // ATHAD: hydrostatic sanity of the initial state, equator and pole. Placed after
        // co2Atmosphere() so every field the profile prints is populated.
        ThermoAtm probe(*this);
        probe.printColumnProfile((jm - 1) / 2, "equator");
        probe.printColumnProfile(4,            "near north pole");
        probe.printLevelSummary("initial state");
    }

//    UtilsAtm(*this) .valueLimitationAtm();                              // value limitation prevents local formation of NANs

    BC_Atm(*this).bcRadius();                                           // extrapolation in i-direction alomg grid boundaries
    BC_Atm(*this).bcTheta();                                            // extrapolation in j-direction alomg grid boundaries
    BC_Atm(*this).bcPhi();                                              // extrapolation in k-direction alomg grid boundaries

    BC_Atm(*this).bcScalarSurfSur();                                    // scalar variable at surfaces extrapolated by von Neumann
    BC_Atm(*this).bcSolidGround();                                      // values inside mountains

    // One-shot Helmholtz projection of the prescribed initial velocity.  The
    // analytical Hadley/Ferrel profile in VelocityInitializer is NOT divergence-
    // free; without this step the incremental projection in PressureSolverAtm
    // strips the dilatational component over the first ~100 iters and the
    // prescribed global circulation decays to noise (max u → <0.5 m/s by iter 150).
    // Doing the projection once here lets the time loop start from a clean
    // ∇·v = 0 state with the solenoidal part of the prescribed flow intact.
    PressureSolverAtm(*this).project_initial_velocity(200);

    UtilsAtm(*this).storeIntermediateData3D(1.0);
    UtilsAtm(*this).findResiduumAtm();


//    goto Printout;


    print_min_max_atm();
    ThermoAtm(*this).printDataAtm();

    UtilsAtm(*this).writeFile(bathymetry_name, output_path, false);     // printing files for ParaView, AtmosphereDataTransfer and AtmospherePlotData


    if(Ma == 0) run_3D_loop(Ma);                                        // iterational 3D loop to solve variables in 4-step Runge-Kutta time scheme

    cout << endl << endl;
    if(Ma > 0) run_3D_loop(Ma);                                         // iterational 3D loop to solve variables in 4-step Runge-Kutta time scheme
    cout << endl << endl;

/*
    Printout:
    print_min_max_atm();
    ThermoAtm(*this).printDataAtm();
    UtilsAtm(*this).writeFile(bathymetry_name, output_path, true);      // printing files for ParaView, AtmosphereDataTransfer and AtmospherePlotData
*/

    RunEnd(at, Ma, start);

    print_final_remarks();
}


/*
  Usage examples

  using namespace AtomUtils;

  // Horizontal wiggle damping (most common use case)
  AtomUtils::damp_wiggles(t, &i_topography, false, true, true);

  // Full 3D, two passes, reduced strength
  AtomUtils::damp_wiggles(cloud, &i_topography, true, true, true, 0.5, 2);

  // Vertical only, no topography mask
  AtomUtils::damp_wiggles(p_stat, nullptr, true, false, false);


    AtomUtils::damp_wiggles(u, &i_topography, true, true, true, 0.5, 2);
    AtomUtils::damp_wiggles(v, &i_topography, true, true, true, 0.5, 2);
    AtomUtils::damp_wiggles(w, &i_topography, true, true, true, 0.5, 2);

    AtomUtils::remove_peaks(BuoyancyForce, &i_bathymetry, true, true, true, 1.0, 2);


    temperature_NASA.printArray_2D(jm, km);
    temp_reconst.printArray_2D(jm, km);
    temp_landscape.printArray_2D(jm, km);

    t.printArray("AGCM", im, jm, km);
    c.printArray("AGCM", im, jm, km);
    cloud.printArray("AGCM", im, jm, km);
    P_rain.printArray("AGCM", im, jm, km);

    fft_gaussian_filter_3d(t,2);                                        // filter for smoothing temperature profiles, three grades of smoothing

    fft_gaussian_filter_3d(c,1);
    fft_gaussian_filter_3d(cloud,1);
    fft_gaussian_filter_3d(ice,1);

    fft_gaussian_filter_3d(P_rain,1);
    fft_gaussian_filter_3d(P_snow,1);
    fft_gaussian_filter_3d(P_graupel,1);
*/

/*
*
*/
void cAtmosphereModel::run_3D_loop(int Ma){ 

cout << endl << endl << endl << "      AGCM: run_3D_loop atm ..........................." << endl;

    // Held-Suarez relaxation target: snapshot the freshly-built initial temperature field
    // (used as-is by the LEGACY path; overwritten by MLR in the radiative path below).
    for (int i = 0; i < im; i++)
        for (int j = 0; j < jm; j++)
            for (int k = 0; k < km; k++)
                t_eq.x[i][j][k] = t.x[i][j][k];

    // ===== INTERACTIVE radiative target t_eq = MLR structure + CO2 surface forcing =====
    //
    // t_eq is the Held-Suarez relaxation target. MultiLayerRadiation computes the multi-layer
    // LW/SW radiative-equilibrium temperature (de-saturated water-vapour + tuned CO2 band
    // emissivity, albedo-reduced shortwave, tridiagonal solve) and we use it as t_eq — giving
    // t_eq a physical radiative/spatial structure (1). On TOP of it we add the CO2 climate
    // sensitivity (2): a Planck-equivalent shift dT_eq = lambda * 5.35*ln(C/C0) (Myhre 1998).
    // MLR's radiative-equilibrium leaves the SURFACE ~unchanged (its correction ~0 there), so
    // the CO2 band shapes OLR/the upper profile while the parametric shift supplies the surface
    // warming (3.7 W/m2/doubling) — complementary, not a surface double-count (option A).
    //
    // The radiation is INTERACTIVE: refresh_radiative_teq() re-solves MLR on the CURRENT field
    // both at start and every teq_refresh_stride iters, so t_eq tracks the evolving climate
    // (water-vapour feedback etc.) instead of being frozen at the initial state. MLR OVERWRITES
    // t with the RE field, so we back t up in t_eq, run MLR, then swap back — the dynamical t is
    // preserved and only the target t_eq is updated. MLR uses i_mount=0 (ignores topography);
    // BC_Atm handles sub-surface. See project_multilayer_radiation.
    // radiation_mode: 0 = legacy (no MLR; t_eq = Scotese snapshot + CO2 forcing),
    //                 1 = OPTION A (MLR radiative equilibrium -> t_eq target, refreshed every
    //                     teq_refresh_stride iters; Held-Suarez relaxes t toward it),
    //                 2 = OPTION B (MLR applied as a DIRECT radiative heating on T each iter:
    //                     operator-split relaxation of the dynamical t toward the RE field by
    //                     omega_rad. t_eq stays the Scotese snapshot + CO2 forcing so Held-
    //                     Suarez still supplies the mean-state APE and radiation heats T on top).
    //   3 = OPTION ii (DEFAULT): t_eq = Scotese snapshot + MLR CO2 PERTURBATION [MLR(C)-MLR(C0)].
    //       MLR now warms the surface with CO2 (surface energy balance), but its grey ABSOLUTE
    //       is ~7 K cold-biased; taking the DIFFERENCE MLR(C)-MLR(C0) cancels that bias and
    //       leaves the physical CO2 warming (surface ~+1.6 K/doubling + profile), added to the
    //       realistic observed baseline. Replaces the 5.35 parametric forcing (mode 0), no cold
    //       bias, no double-count. The parametric CO2 shift is now used ONLY in mode 0.
    //   4 = OPTION ii + SEED (DEFAULT): like mode 3, but the CO2 perturbation is ALSO added
    //       directly to the prognostic t (not only the t_eq target). The Held-Suarez relaxation
    //       is too gentle to transmit the t_eq shift to t within a single run (a co2x1/co2x2 pair
    //       otherwise differences to ~0 since both start t from the same Scotese snapshot); seeding
    //       t makes the two runs evolve to distinct, CO2-separated equilibria that the shifted t_eq
    //       then anchors. This is the working single-run CO2 sensitivity path.
    // ATHAD: mode 2 — direct sigma*T^4 radiative heating from the multi-layer solve.
    //
    // The inherited default was mode 5, Newtonian relaxation of t toward a target built
    // from the Scotese paleo-temperature curve. There IS no Scotese curve at 4.4 Ga: it
    // stops around 540 Ma. Modes 0, 1, 3 and 4 all lean on the same snapshot or on the
    // CO2-perturbation machinery calibrated to Earth's 280 ppm reference, so mode 2 is
    // the only one with a defensible meaning here — the radiation has to SET the profile
    // rather than nudge it toward a prescribed one.
    constexpr int    radiation_mode     = 2;               // mode 2: direct sigma*T^4 heating from MLR
    constexpr int    teq_refresh_stride = 20;              // option A: MLR re-solve cadence
    constexpr double omega_rad          = 0.05;            // option B: radiative heating (frac/iter toward RE)
    constexpr double omega_teq          = 0.20;            // mode 5: relaxation rate toward t_eq (frac/iter)
    constexpr double co2_ref_ppm        = 280.0;           // pre-industrial reference C0
    constexpr double co2_lambda         = 0.31;            // Planck (no-feedback) sensitivity [K/(W/m2)]
    // Deep (tropospheric-depth) CO2 perturbation shape (modes 3/4/5). The MLR CO2 response is
    // surface-trapped; imposing it only near the ground lets the strong t_eq relaxation aloft
    // (toward the UNperturbed Scotese profile) suppress the deep warming -> weak, shallow
    // sensitivity. Instead spread the SURFACE perturbation coherently through the troposphere
    // (uniform, optionally moist-adiabatically amplified toward the tropopause) then taper it to
    // zero across the tropopause, so t_eq carries the warming up the whole column.
    constexpr double co2_pert_ztrop     = 11000.0;         // uniform-deep up to here [m]
    constexpr double co2_pert_ztop      = 14000.0;         // perturbation tapered to 0 by here [m]
    constexpr double co2_pert_amp       = 1.0;             // amplification at tropopause (1 = uniform-deep)

    auto apply_co2_forcing = [&]() {                       // Planck-equivalent CO2 shift of t_eq
        if (co2_0 <= 0.0) return;
        const double dT_eq_nd = (co2_lambda * 5.35 * std::log(co2_0 / co2_ref_ppm)) / t_0;
        for (int i = 0; i < im; i++)
            for (int j = 0; j < jm; j++)
                for (int k = 0; k < km; k++)
                    t_eq.x[i][j][k] += dT_eq_nd;
    };

    // Option A helper: t_eq := MLR radiative equilibrium of the CURRENT field (+ CO2 forcing).
    // MLR OVERWRITES t, so t_eq doubles as the backup buffer (swap restores the dynamical t).
    auto refresh_radiative_teq = [&]() {
        for (int i = 0; i < im; i++)
            for (int j = 0; j < jm; j++)
                for (int k = 0; k < km; k++)
                    t_eq.x[i][j][k] = t.x[i][j][k];             // back up dynamical t
        MultiLayerRadiation(*this).run();                       // t := radiative equilibrium
        for (int i = 0; i < im; i++)
            for (int j = 0; j < jm; j++)
                for (int k = 0; k < km; k++) {                  // swap: restore t, set t_eq := RE
                    const double tmp = t.x[i][j][k];
                    t.x[i][j][k]     = t_eq.x[i][j][k];
                    t_eq.x[i][j][k]  = tmp;
                }
        // No separate CO2 forcing here: MLR now warms the surface with CO2 via its surface
        // energy balance, so adding the 5.35 shift would double-count. (Grey absolute stays
        // ~7 K cold-biased in this mode — mode 3 avoids that via the difference.)
    };

    // Mode 5 helper: strong Newtonian relaxation of the dynamical t toward t_eq (= Scotese/
    // observed baseline + the MLR CO2 perturbation, spread through the troposphere and refreshed
    // in-loop on the warming field by apply_co2_perturbation(false) for the water-vapour feedback).
    // REPLACES the near-inert in-RHS Held-Suarez term (k_a=1/4day enters as ~1e-8/iter,
    // never imposes t_eq); omega_teq per iter is dynamically competitive so the CO2 signal PERSISTS
    // to equilibrium AND the absolute climate stays anchored to the realistic baseline (no grey
    // cold bias, unlike mode 2). The realistic equator-pole gradient in t_eq preserves the mean
    // circulation driving; only transient eddies are damped (canonical Held-Suarez at proper
    // strength). t_new = (1-omega_teq)*t_old + omega_teq*t_eq.
    auto apply_teq_relaxation = [&]() {
        for (int i = 0; i < im; i++)
            for (int j = 0; j < jm; j++)
                for (int k = 0; k < km; k++)
                    t.x[i][j][k] = (1.0 - omega_teq) * t.x[i][j][k] + omega_teq * t_eq.x[i][j][k];
    };

    // Option B helper: nudge the dynamical t toward the MLR radiative-equilibrium field by
    // omega_rad (radiation as a direct heating/tendency). Own backup buffer — t_eq is busy
    // holding the Held-Suarez target here. t_new = (1-omega)*t_old + omega*RE.
    std::vector<double> t_backup;
    auto apply_radiative_heating = [&]() {
        if (t_backup.empty()) t_backup.assign((size_t)im * jm * km, 0.0);
        size_t idx = 0;
        for (int i = 0; i < im; i++)
            for (int j = 0; j < jm; j++)
                for (int k = 0; k < km; k++)
                    t_backup[idx++] = t.x[i][j][k];             // save dynamical t
        MultiLayerRadiation(*this).run();                       // t := radiative equilibrium
        idx = 0;
        for (int i = 0; i < im; i++)
            for (int j = 0; j < jm; j++)
                for (int k = 0; k < km; k++) {                  // blend t toward RE by omega_rad
                    const double re = t.x[i][j][k];
                    t.x[i][j][k] = t_backup[idx++] * (1.0 - omega_rad) + re * omega_rad;
                }
    };

    // Option (ii) helper: t_eq += MLR(current CO2) - MLR(reference 280 ppm). Two MLR solves on
    // the current field; MLR's grey absolute cold bias is identical in both and CANCELS, leaving
    // the physical CO2 warming (surface + profile) which is added to the observed t_eq baseline.
    // t and co2 are saved/restored around the two solves. Computed once at start (fixed CO2).
    // seed_prognostic (mode 4): in addition to shifting t_eq, also add the CO2 perturbation
    // directly to the prognostic t field, so the warming enters the DYNAMICAL temperature
    // immediately rather than having to be transmitted from t_eq by the (gentle) Held-Suarez
    // relaxation. Without this the CO2 signal never reaches t in a single run and a co2x1/co2x2
    // pair differences to ~0 (both runs start t from the identical Scotese snapshot).
    // Vertical weight for the deep perturbation: 1 through the troposphere (optionally amplified
    // toward the tropopause), linear taper to 0 across [ztrop, ztop]. Built once from the layer
    // heights (fixed grid). w_pert[0] == 1 so the surface diagnostic below is unaffected.
    std::vector<double> w_pert(im, 0.0);
    for (int i = 0; i < im; i++) {
        const double z = get_layer_height(i);
        w_pert[i] = (z <= co2_pert_ztrop)
            ? 1.0 + (co2_pert_amp - 1.0) * (z / co2_pert_ztrop)
            : (z <  co2_pert_ztop)
                ? co2_pert_amp * (co2_pert_ztop - z) / (co2_pert_ztop - co2_pert_ztrop)
                : 0.0;
    }
    // Pristine Scotese/observed baseline. apply_co2_perturbation rebuilds t_eq = t_eq_base + deep
    // perturbation from this each call (idempotent), so it can be re-run in-loop on the warming
    // field for the water-vapour feedback without accumulating.
    std::vector<double> t_eq_base((size_t)im * jm * km);
    { size_t idx = 0; for (int i = 0; i < im; i++) for (int j = 0; j < jm; j++) for (int k = 0; k < km; k++)
        t_eq_base[idx++] = t_eq.x[i][j][k]; }

    std::vector<double> t_save, mlr_hi, co2_save, dpert_surf, cloud_save, ice_save;
    auto apply_co2_perturbation = [&](bool seed_prognostic) {
        if (co2_0 <= 0.0) return;
        const size_t N = (size_t)im * jm * km;
        if (t_save.empty()) { t_save.assign(N, 0.0); mlr_hi.assign(N, 0.0); co2_save.assign(N, 0.0);
                              dpert_surf.assign((size_t)jm * km, 0.0);
                              cloud_save.assign(N, 0.0); ice_save.assign(N, 0.0); }
        // Compute the CO2 perturbation CLEAR-SKY: save then zero cloud/ice around the two MLR
        // solves so grey-opaque clouds don't LW-saturate and mask the CO2 marginal forcing (they
        // mask ~40-70% otherwise — a grey-radiation artifact, not a real cloud feedback). Cloud/ice
        // stay radiatively active for the absolute-radiation diagnostic and modes 1-2; only this
        // difference-based CO2 forcing is taken clear-sky. Restored in the final loop below.
        size_t idx = 0;
        for (int i = 0; i < im; i++) for (int j = 0; j < jm; j++) for (int k = 0; k < km; k++) {
            t_save[idx] = t.x[i][j][k]; co2_save[idx] = co2.x[i][j][k];
            cloud_save[idx] = cloud.x[i][j][k]; ice_save[idx] = ice.x[i][j][k];
            cloud.x[i][j][k] = 0.0; ice.x[i][j][k] = 0.0;         // clear-sky for both solves
            idx++;
        }
        MultiLayerRadiation(*this).run();                         // t := MLR(current CO2, clear-sky)
        idx = 0;
        for (int i = 0; i < im; i++) for (int j = 0; j < jm; j++) for (int k = 0; k < km; k++)
            mlr_hi[idx++] = t.x[i][j][k];
        idx = 0;                                                  // co2.x is ppm-valued (not a ratio),
        for (int i = 0; i < im; i++) for (int j = 0; j < jm; j++) for (int k = 0; k < km; k++) {
            t.x[i][j][k] = t_save[idx++]; co2.x[i][j][k] = co2_ref_ppm;   // so the reference field = 280 ppm directly
        }
        MultiLayerRadiation(*this).run();                         // t := MLR(reference CO2)
        idx = 0;
        double dsurf_sum = 0.0, dsurf_w = 0.0, dsurf_max = -1e30, dcol_sum = 0.0, dcol_w = 0.0;
        for (int i = 0; i < im; i++) for (int j = 0; j < jm; j++) for (int k = 0; k < km; k++) {
            const double dpert_raw = mlr_hi[idx] - t.x[i][j][k];   // raw MLR CO2 response (surface-trapped)
            if (i == 0) dpert_surf[(size_t)j * km + k] = dpert_raw;// i is the outer loop -> surface filled first
            const double dpert = dpert_surf[(size_t)j * km + k] * w_pert[i]; // deep tropospheric profile
            const double dK = dpert * t_0;                         // perturbation in Kelvin
            t_eq.x[i][j][k] = t_eq_base[idx] + dpert;             // rebuild from pristine baseline (idempotent)
            const double w = cos((j / (double)(jm - 1) - 0.5) * M_PI);
            dcol_sum += w * dK; dcol_w += w;
            if (i == 0) { dsurf_sum += w * dK; dsurf_w += w; if (dK > dsurf_max) dsurf_max = dK; }
            // mode 4 seeds the prognostic t with the perturbation; mode 3 restores the Scotese baseline.
            t.x[i][j][k] = seed_prognostic ? (t_save[idx] + dpert) : t_save[idx];
            co2.x[i][j][k] = co2_save[idx];
            cloud.x[i][j][k] = cloud_save[idx]; ice.x[i][j][k] = ice_save[idx];   // restore clouds
            idx++;
        }
        std::cout << "      AGCM: [co2-perturb DIAG] t_eq shift vs 280ppm (co2_scale=" << co2_scale
                  << "): surface cos-lat mean = " << (dsurf_w > 0 ? dsurf_sum / dsurf_w : 0.0)
                  << " K, surface max = " << dsurf_max
                  << " K, whole-domain mean = " << (dcol_w > 0 ? dcol_sum / dcol_w : 0.0) << " K" << std::endl;
    };

    // Mode 4 DIAGNOSTIC refresh: re-solve MLR on the CURRENT field so the radiation/epsilon
    // arrays track the evolving climate instead of staying frozen at the start-of-run solve
    // (mode 4 seeds the CO2 perturbation into t once and otherwise never calls MLR in-loop, so
    // radiation.x never changed — e.g. identical at iter 100 and 300). MLR overwrites t with its
    // radiative-equilibrium field, so back t up and restore it afterwards: the dynamical t and the
    // already-seeded CO2 perturbation are untouched; only the diagnostic radiation.x/epsilon (and
    // the σT^4-consistent RE profile behind them) are refreshed against the actual current CO2.
    std::vector<double> t_diag_backup;
    auto refresh_radiation_diag = [&]() {
        if (t_diag_backup.empty()) t_diag_backup.assign((size_t)im * jm * km, 0.0);
        size_t idx = 0;
        for (int i = 0; i < im; i++) for (int j = 0; j < jm; j++) for (int k = 0; k < km; k++)
            t_diag_backup[idx++] = t.x[i][j][k];                 // save dynamical t
        MultiLayerRadiation(*this).run();                        // fills radiation.x/epsilon; overwrites t
        idx = 0;
        for (int i = 0; i < im; i++) for (int j = 0; j < jm; j++) for (int k = 0; k < km; k++)
            t.x[i][j][k] = t_diag_backup[idx++];                 // restore dynamical t (radiation.x kept)
    };

    // Mode-5 CLOUD radiation diagnostic. Mode 5 otherwise calls MLR only clear-sky (inside the CO2
    // perturbation), so radiation.x/epsilon never see the cloud/ice fields. This runs MLR twice on
    // the CURRENT field — clear-sky then cloudy — leaving radiation.x/epsilon CLOUDY (visible in the
    // output) and printing the cloud longwave effect = cos-lat-mean top-layer radiation clear minus
    // cloudy. Diagnostic only: cloud/ice and t are saved and restored, dynamics untouched.
    std::vector<double> cld_diag_save, ice_diag_save;
    auto cloud_radiation_diag = [&]() {
        const size_t N = (size_t)im * jm * km;
        if (t_diag_backup.empty()) t_diag_backup.assign(N, 0.0);
        if (cld_diag_save.empty()) { cld_diag_save.assign(N, 0.0); ice_diag_save.assign(N, 0.0); }
        size_t idx = 0;
        for (int i = 0; i < im; i++) for (int j = 0; j < jm; j++) for (int k = 0; k < km; k++) {
            t_diag_backup[idx] = t.x[i][j][k];
            cld_diag_save[idx] = cloud.x[i][j][k]; ice_diag_save[idx] = ice.x[i][j][k];
            cloud.x[i][j][k] = 0.0; ice.x[i][j][k] = 0.0;        // clear-sky pass
            idx++;
        }
        MultiLayerRadiation(*this).run();
        double w_sum = 0.0, olr_clear = 0.0;
        for (int j = 0; j < jm; j++) { const double w = cos((j / (double)(jm - 1) - 0.5) * M_PI);
            w_sum += w * km;
            for (int k = 0; k < km; k++) olr_clear += w * radiation.x[im - 1][j][k]; }
        idx = 0;                                                 // restore t + clouds, cloudy pass
        for (int i = 0; i < im; i++) for (int j = 0; j < jm; j++) for (int k = 0; k < km; k++) {
            t.x[i][j][k] = t_diag_backup[idx];
            cloud.x[i][j][k] = cld_diag_save[idx]; ice.x[i][j][k] = ice_diag_save[idx];
            idx++;
        }
        MultiLayerRadiation(*this).run();                        // radiation.x/epsilon end CLOUDY (kept for output)
        double olr_cloudy = 0.0;
        for (int j = 0; j < jm; j++) { const double w = cos((j / (double)(jm - 1) - 0.5) * M_PI);
            for (int k = 0; k < km; k++) olr_cloudy += w * radiation.x[im - 1][j][k]; }
        idx = 0;                                                 // restore dynamical t (radiation.x kept cloudy)
        for (int i = 0; i < im; i++) for (int j = 0; j < jm; j++) for (int k = 0; k < km; k++)
            t.x[i][j][k] = t_diag_backup[idx++];
        std::cout << "      AGCM: [cloud-rad DIAG] cos-lat-mean top radiation: clear=" << olr_clear / w_sum
                  << " cloudy=" << olr_cloudy / w_sum << " W/m2  ->  cloud LW effect="
                  << (olr_clear - olr_cloudy) / w_sum << " W/m2" << std::endl;
    };

    if      (radiation_mode == 1) refresh_radiative_teq();       // A: MLR absolute -> t_eq
    else if (radiation_mode == 3) apply_co2_perturbation(false); // ii: Scotese t_eq + MLR CO2 perturbation
    else if (radiation_mode == 4) apply_co2_perturbation(true);  // ii+seed: also inject perturbation into t
    else if (radiation_mode == 5) apply_co2_perturbation(false); // 5: build t_eq = Scotese + CO2 perturbation
    else if (radiation_mode == 0) apply_co2_forcing();          // legacy: Scotese + 5.35 parametric forcing
    // mode 2: t_eq stays the Scotese snapshot; the in-loop MLR heating carries the radiation.
    std::cout << "      AGCM: radiation_mode=" << radiation_mode
              << " (0=legacy,1=MLR->t_eq,2=MLR->heat,3=Scotese+MLR-CO2-perturb,4=+seed-t,5=strong-relax-teq)  co2_0="
              << co2_0 << " kg/kg" << std::endl;

    // Restart shortcut: overwrite the just-built initial conditions with a saved 3D
    // state and resume from its total_iter_count (skips re-running the dry spin-up).
    const bool restarted = (restart_from_iter >= 0) && load_state(restart_from_iter, Ma);

    // On a successful restart, start the loop counter at the restart point so the
    // printed iter_n matches where the run actually resumed (not back at 1).
    const int iter_start = restarted ? restart_from_iter : 1;

    // ---- Convergence monitor (report-only) ----------------------------------------
    // Makes "is it equilibrated yet?" visible. Every conv_stride iters it samples two
    // slow, area+mass-weighted integral metrics — mean temperature (thermodynamics) and
    // mean kinetic energy (circulation) — and reports the relative drift of the LAST
    // conv_window samples' mean vs the PREVIOUS window's mean. Window-averaging cancels
    // the precipitation/limit-cycle wobble, so the drift is the SECULAR trend: the run is
    // equilibrated when BOTH metrics drift < conv_tol_pct per window. Circulation is
    // tracked alongside T because a run can look thermally settled while still spinning
    // down. Report-only (never auto-stops; nm is the hard ceiling); also logged to
    // output_path/convergence.csv. Cost: one grid reduction per conv_stride iters.
    // T convergence is gated on an ABSOLUTE K drift, not a percent of Kelvin: the mean sits
    // near 260 K, so the old 1 %-of-Kelvin test (~2.6 K/window) reported CONVERGED while the
    // circulation was still visibly evolving — i.e. the flag was effectively KE-only. KE
    // keeps the percent test (its zero IS physical, so a relative drift is meaningful there).
    // The percent T drift is still logged for continuity. Mirrors the hydro monitor.
    constexpr int    conv_stride  = 25;    // sample cadence [iters]
    constexpr int    conv_window  = 4;     // samples per window (-> 100-iter trailing window)
    constexpr double conv_tol_pct = 1.0;   // KE convergence threshold [% drift per window]
    constexpr double conv_tol_T_K = 0.01;  // T  convergence threshold [K drift per window]
    std::vector<double> conv_histT, conv_histKE;
    auto domain_mean_T = [&]() -> double {
        double sum = 0.0, wsum = 0.0;
        #pragma omp parallel for reduction(+:sum,wsum) schedule(static)
        for(int j = 0; j < jm; j++){
            const double coslat = cos((j/(double)(jm-1) - 0.5) * M_PI);
            for(int i = 0; i < im; i++){
                const double dz = (i < im-1) ? (get_layer_height(i+1) - get_layer_height(i)) : 0.0;
                const double wt = coslat * dz;
                for(int k = 0; k < km; k++){ sum += wt * t.x[i][j][k] * t_0; wsum += wt; }
            }
        }
        return (wsum > 0.0) ? sum / wsum : 0.0;                 // [K]
    };
    auto domain_mean_KE = [&]() -> double {
        double sum = 0.0, wsum = 0.0;
        #pragma omp parallel for reduction(+:sum,wsum) schedule(static)
        for(int j = 0; j < jm; j++){
            const double coslat = cos((j/(double)(jm-1) - 0.5) * M_PI);
            for(int i = 0; i < im; i++){
                const double dz = (i < im-1) ? (get_layer_height(i+1) - get_layer_height(i)) : 0.0;
                const double wt = coslat * dz;
                for(int k = 0; k < km; k++){
                    const double uu = u.x[i][j][k], vv = v.x[i][j][k], ww = w.x[i][j][k];
                    sum += wt * 0.5 * (uu*uu + vv*vv + ww*ww); wsum += wt;
                }
            }
        }
        return (wsum > 0.0) ? (sum / wsum) * (u_0 * u_0) : 0.0;  // [m2/s2] (velocities are /u_0)
    };
    // returns drift %/window (or -1.0 while still warming up, < 2 windows of samples);
    // abs_out receives the drift in the metric's own units (K for T, m2/s2 for KE).
    auto conv_drift = [&](std::vector<double>& hist, double metric, double& abs_out) -> double {
        hist.push_back(metric);
        const int n = (int)hist.size();
        if(n < 2*conv_window){ abs_out = -1.0; return -1.0; }
        double cur = 0.0, prev = 0.0;
        for(int i = n-conv_window;   i < n;             i++) cur  += hist[i];
        for(int i = n-2*conv_window; i < n-conv_window; i++) prev += hist[i];
        cur /= conv_window; prev /= conv_window;
        abs_out = std::fabs(cur - prev);
        return 100.0 * abs_out / std::max(1e-12, std::fabs(cur));
    };
    { std::ofstream cf(output_path + "/convergence.csv");         // (re)create with header
      if(cf) cf << "iter,mean_T_K,drift_T_pct,drift_T_K,mean_KE_m2s2,drift_KE_pct,converged\n"; }
    // --------------------------------------------------------------------------------

    for(iter_n = iter_start; iter_n <= nm; iter_n++){
        print_loop_3D_headings();

        // Inviscid spin-up: zero diffusion for the first `inviscid_spinup_iters`
        // cumulative iterations, then ramp linearly to full over `inviscid_ramp_iters`.
        // total_iter_count is per-RHS-call and persists across time slices.
        total_iter_count++;
        if(inviscid_spinup_iters <= 0){
            inviscid_phase = false;
            diffusion_ramp = 1.0;
        } else if(total_iter_count <= inviscid_spinup_iters){
            inviscid_phase = true;
            diffusion_ramp = 0.0;
        } else {
            inviscid_phase = false;
            int over = total_iter_count - inviscid_spinup_iters;
            diffusion_ramp = (inviscid_ramp_iters > 0)
                ? std::min(1.0, static_cast<double>(over) / inviscid_ramp_iters)
                : 1.0;
        }
        dt = inviscid_phase ? dt_inviscid : dt_visc;

        // Buoyancy ramp: linearly from 0 at iter 0 to 1 at buoyancy_ramp_iters.
        // Spreads the build-up of the thermal driving force across many iters so the
        // multi-sweep Jacobi can converge p_dyn against the growing forcing instead
        // of falling behind and ringing.
        buoyancy_ramp = (buoyancy_ramp_iters > 0)
            ? std::min(1.0, static_cast<double>(total_iter_count) / buoyancy_ramp_iters)
            : 1.0;

        cout << "      AGCM: iter = " << total_iter_count
             << "  inviscid_phase = " << (inviscid_phase ? "true" : "false")
             << "  diffusion_ramp = " << std::fixed << std::setprecision(3) << diffusion_ramp
             << "  buoyancy_ramp = " << std::fixed << std::setprecision(3) << buoyancy_ramp
             << "  dt = " << std::scientific << std::setprecision(4) << dt
             << std::defaultfloat << endl;



        // Pressure-solve cadence. EXPERIMENT (2026-06-10, user hypothesis): widen the stride
        // 2→4 so p_dyn is held fixed for 3 of every 4 global iterations while the velocity
        // evolves under that frozen pressure-gradient force. Rationale: the residual coastal
        // p_dyn oscillation (visible at i=0 over h=1 cells: large +/− jumps between adjacent
        // coastal columns) looks like a 2Δt ping-pong between the projection and the velocity
        // update; subcycling the projection breaks that resonance and lets the velocity
        // components "come to rest" between solves, while leaving the emergent i=0 vortices
        // (which we want to keep) intact. Trade-off: divergence accumulates over the 3 frozen
        // steps and is cleaned only on the 4th — watch whether the Pamir runaway worsens.
        // NOTE (2026-06-09): re-enabling 20-sweep Jacobi here was TESTED and made things WORSE
        // — converging the Poisson more faithfully drove pgr at the Tian Shan column from −2.6
        // to −5.2 and pegged p_dyn at its −10 ceiling sooner, proving the PRESSURE faithfully
        // balances a pathological DIVERGENCE SOURCE (pgr→u→divergence→pgr), not the root cause.
        // Single sweep retained. See [[project-upper-velocity-secular-growth]].
        constexpr int pressure_stride = 4;   // was 2
        if(iter_n % pressure_stride == 0){
            PressureSolverAtm(*this).run();
        }



        // Two cadences during inviscid spin-up:
        //   moist_stride    — PressureSolverAtm, SaturationAdjustment, ice scheme,
        //                     MoistConvection, ThermoAtm, turbulence sources. Throttled
        //                     because (a) the moist/thermo sources feed back unbounded
        //                     in the undamped (diffusion_ramp=0) flow, and (b) the
        //                     pressure solver is a single Jacobi sweep whose target
        //                     balances the full body-force field — calling it every step
        //                     races p_dyn toward that (very large) balance value.
        //   momentum_stride — all BCs (radial/theta/phi extrapolation, scalar surface,
        //                     solid ground). Run every step during spin-up so the
        //                     topographic walls are re-enforced before each RK4 cycle —
        //                     this is what lets the flow "see" the mountains and form
        //                     lee/windward circulations.
        // In the viscous phase moist_stride collapses back to 2 (original behaviour).
        // momentum_stride is kept at 1 in BOTH phases: with stride 2 the solid-ground
        // wall BC was re-enforced only every other RK4 cycle, leaving a 2Δt computational
        // mode pinned at steep high-latitude coasts (Antarctic, ×2.5 via 1/sinθ at 68°S)
        // that grew to a near-surface NaN seed by ~iter 142. Re-enforcing the walls every
        // step (the existing, tested BCs) suppresses it; cost is BCs running 2× as often.
        int moist_stride    = inviscid_phase ? 10  : 2;
        int momentum_stride = 1;

        if(iter_n % moist_stride == 0){
/*
            // Multi-sweep Jacobi: the original single-sweep cadence (one sweep every
            // moist_stride iters) couldn't keep up with the divergence accumulated by
            // a fully-spun-up global circulation, letting p_dyn drift to ±750 hPa by
            // iter 120 and the simulation crash into the velocity caps by iter 300.
            // 20 sweeps per call (~10× the previous rate) actually converges the
            // Poisson before the next time step, so dpdr/dpdthe/dpdphi in the RHS
            // reflect the real pressure each iter.  Verbose only on the first sweep.
            {
//                constexpr int n_pressure_sweeps = 20;
                constexpr int n_pressure_sweeps = 2;
                PressureSolverAtm solver(*this);
                for(int s = 0; s < n_pressure_sweeps; s++){
                    solver.run(s == 0);
                }
            }
*/
            // Moist-physics gate: the implicit microphysics (SaturationAdjustment,
            // ice scheme, MoistConvection) is too stiff during the initial velocity
            // transient — a single tropical maritime column repeatedly drove q_c, q_i
            // and S_s into runaway. Delay activation until the circulation has formed
            // on a dry field; the moist physics then sees a stable advection pattern
            // it was designed for. 0 disables the gate (always on).
            const bool moist_phys_active = (moist_phys_start_iter <= 0)
                                        || (total_iter_count > moist_phys_start_iter);

            if(moist_phys_active){
                // Water budget either side of the moist physics, so the creation this
                // reports can be ATTRIBUTED: everything before this point is RK4 transport,
                // everything after it is condensation, the ice scheme and the limiters.
                ThermoAtm(*this).waterBudget(iter_n % diagnosticStride() == 0, "post-RK4");

                SaturationAdjustment(*this).run();                      // based on the initial distribution, recomputation of the cloud water and cloud ice formation in case of saturated water vapour detected
                // Temperature 2Δt de-checkerboard. c/cloud/ice receive the same stride-2 moist
                // forcing AND damp_wiggles and stay stable; t got the forcing but NO damping —
                // the only prognostic without it — letting an undamped 2Δt computational mode
                // at the steep Pamir surface (t[0] sampled by densities() on its diverging odd
                // branch) run away to 53 K and crash the barometric formula. Mirror the scalar
                // treatment on t to close the asymmetry. See [[project_upper_velocity_secular_growth]].
                AtomUtils::damp_wiggles(t,     &i_topography, true, true, true);
                AtomUtils::damp_wiggles(ice,   &i_topography, true, true, true);
                AtomUtils::damp_wiggles(c,     &i_topography, true, true, true);
                AtomUtils::damp_wiggles(cloud, &i_topography, true, true, true);

                switch(CategoryIceScheme){                              // rain, snow graupel and precipitation production and reduction
                    case -1: cout << endl << endl << endl               // no CategoryIceScheme used
                        << "  no CategoryIceScheme used" << endl;
                        break;
                    case 0: ZeroCatIceScheme(*this).run();              // development of rain and snow fall, water vapour and cloud water
                        break;
                    case 1: OneCatIceScheme(*this).run();               // development of rain and snow fall, water vapour, cloud water and ice
                        break;
                    case 2: TwoCatIceScheme(*this).run();               // development of rain and snow fall, water vapour, cloud water and ice
                        break;
                    case 3: ThreeCatIceScheme(*this).run();             // development of rain and snow fall, water vapour, cloud water, ice and graupel
                        break;
                }

                AtomUtils::damp_wiggles(P_rain, &i_topography, true, true, true);
                AtomUtils::damp_wiggles(P_snow, &i_topography, true, true, true);

                MoistConvection(*this).run(iter_n);                     // rainfall from convecting clouds

                AtomUtils::damp_wiggles(P_conv, &i_topography, true, true, true);
                AtomUtils::damp_wiggles(E_u,    &i_topography, true, true, true);
                AtomUtils::damp_wiggles(M_u,    &i_topography, true, true, true);
                AtomUtils::damp_wiggles(q_v_u,  &i_topography, true, true, true);
                AtomUtils::damp_wiggles(q_c_u,  &i_topography, true, true, true);
                AtomUtils::damp_wiggles(E_d,    &i_topography, true, true, true);
                AtomUtils::damp_wiggles(M_d,    &i_topography, true, true, true);
                AtomUtils::damp_wiggles(q_v_d,  &i_topography, true, true, true);
                // The convective tendencies fed straight into rhs_t/rhs_v/rhs_w
                // (coeff_MC_*·MC in RHS_Atm) enter spatially rough over steep orography;
                // smooth them like the other MC fields so the momentum/heat source does
                // not inject a jagged divergence the pressure solver amplifies (Andes seed).
                AtomUtils::damp_wiggles(MC_t,   &i_topography, true, true, true);
                AtomUtils::damp_wiggles(MC_v,   &i_topography, true, true, true);
                AtomUtils::damp_wiggles(MC_w,   &i_topography, true, true, true);

                UtilsAtm(*this).precipitationSum();

                // Physical caps on the microphysics source terms.
                // The ice-scheme S-terms feed rhs_t (latent heat: S_c,S_r,S_i,S_s,S_g)
                // and rhs_c (moisture: S_v) UNCAPPED. They depend on c/cloud/ice, so
                // c,cloud,ice ↔ S_* ↔ rhs_t/rhs_c is a closed loop. Over the steep Andes
                // (≈50°S, 1/sin²θ-amplified orographic ascent) the condensation↔latent-
                // heat↔ascent feedback drives them to ~1e260 by iter ~309 — the leg of
                // the moist runaway the MoistConvection caps (rhsForcing) don't cover.
                // Capping S_* bounds rhs_{c,cloud,ice,t}, hence c/cloud/ice, hence S_*
                // next call. Healthy |S| ≲ 1e-4 (kg/kg)/s (e.g. S_c_c, S_s diagnostics),
                // so 1e-3 is ~10× headroom and clips only the runaway. NaN/Inf→0 via the
                // bit-level is_finite_safe (std::min/max are unreliable under -ffast-math).
                {
                    constexpr double S_max = 1.0e-3;                     // [(kg/kg)/s]
                    auto cap_S = [&](Array& S){
                        #pragma omp parallel for collapse(2) schedule(static)
                        for(int i = 0; i < im; i++){
                            for(int j = 0; j < jm; j++){
                                for(int k = 0; k < km; k++){
                                    double v = S.x[i][j][k];
                                    if(!AtomUtils::is_finite_safe(v)) S.x[i][j][k] = 0.0;
                                    else if(v >  S_max)               S.x[i][j][k] =  S_max;
                                    else if(v < -S_max)               S.x[i][j][k] = -S_max;
                                }
                            }
                        }
                    };
                    cap_S(S_v); cap_S(S_c); cap_S(S_i); cap_S(S_r);
                    cap_S(S_s); cap_S(S_g); cap_S(S_c_c);
                }

                // Microphysics state clamp.  ThermoAtm::densities computes
                // r_humid = scale·p_i / ((1 + 0.608·c − cloud − ice)·t_u);
                // if cloud + ice exceeds ~1 the denominator factor flips sign and
                // r_humid goes negative, inverting -∇p/ρ in the momentum equations
                // and producing the iter-356 NaN cascade at 62°N (Cook Inlet).
                // RK4 central-difference advection of a sharp surface cloud peak
                // also drives undershoots (negative cloud) at i=5,6 in the same
                // column.  Defensive clamp: c, cloud, ice ≥ 0 and cloud+ice ≤
                // CLOUD_ICE_MAX, applied EVERY moist iter before densities reads
                // them.  See [[project-cloud-runaway-cook-inlet]].
                {
                    // Same parameter as SaturationAdjustment's per-species cap, so the
                    // two cannot drift apart. This one bounds the SUM and is the binding
                    // one. See param.py (cloud_cap): on ATHAD it is not a backstop, it is
                    // setting the condensate.
                    const double CLOUD_ICE_MAX = cloud_cap;                 // [kg/kg]
                    #pragma omp parallel for collapse(2) schedule(static)
                    for (int i = 0; i < im; i++) {
                        for (int j = 0; j < jm; j++) {
                            for (int k = 0; k < km; k++) {
                                double cv  = c.x[i][j][k];
                                double cld = cloud.x[i][j][k];
                                double ic  = ice.x[i][j][k];
                                if (!AtomUtils::is_finite_safe(cv))  cv  = 0.0;
                                if (!AtomUtils::is_finite_safe(cld)) cld = 0.0;
                                if (!AtomUtils::is_finite_safe(ic))  ic  = 0.0;
                                if (cv  < 0.0) cv  = 0.0;
                                if (cld < 0.0) cld = 0.0;
                                if (ic  < 0.0) ic  = 0.0;
                                // ATHAD: water vapour's physical CEILING, enforced here
                                // because this is the last thing to touch c before
                                // densities() and the diagnostics read it. The three mass
                                // fractions sum to one and the background is the remainder
                                // 1 - c - co2, so c > 1 - co2 is a negative background mass.
                                // UtilsAtm::valueLimitationAtm carries the same bound but is
                                // commented out at both of its call sites, and the sub-cloud
                                // band was sitting at c = 1.0000 against co2 = 0.2053 — a
                                // composition summing to 1.21, invisible because
                                // AtmMixture::split() renormalises defensively and only the
                                // gas constant, pinned at 415.1, showed it. On Earth water is
                                // ~1 % of the mass, so 1 and 1 - co2 are the same number.
                                const double cv_max = std::max(0.0, 1.0 - co2.x[i][j][k]);
                                if (cv > cv_max) cv = cv_max;
                                const double sum = cld + ic;
                                if (sum > CLOUD_ICE_MAX) {
                                    const double scale = CLOUD_ICE_MAX / sum;
                                    cld *= scale;
                                    ic  *= scale;
                                }
                                c.x[i][j][k]     = cv;
                                cloud.x[i][j][k] = cld;
                                ice.x[i][j][k]   = ic;
                            }
                        }
                    }
                }
            }  // moist_phys_active

            // ATHAD: water vapour's physical CEILING, enforced EVERY iteration.
            //
            // The three mass fractions sum to one and the background is carried as the
            // remainder 1 - c - co2, so c > 1 - co2 is a negative background mass.
            // AtmMixture::split() renormalises defensively, so the violation never showed
            // as a crash — only as a gas constant pinned at 415.1 instead of moving with
            // the composition, while the sub-cloud band sat at c = 1.0000 against
            // co2 = 0.2053, a composition summing to 1.21.
            //
            // It lives here, outside the moist-physics block, because that block runs only
            // on moist iterations while the RK4 transport of c runs on every one, and
            // UtilsAtm::valueLimitationAtm — which carries the same bound, at the wrong
            // ceiling of 1 — is commented out at both of its call sites. On Earth water is
            // ~1 % of the mass and CO2 is ppm, so 1 and 1 - co2 are the same number.
            {
                long   n_ceiling = 0;
                double q_clipped = 0.0;      // mass-weighted water deleted this iteration
                double w_air     = 0.0;
                #pragma omp parallel for collapse(2) schedule(static) \
                        reduction(+:n_ceiling,q_clipped,w_air)
                for (int i = 0; i < im; i++)
                    for (int j = 0; j < jm; j++)
                        for (int k = 0; k < km; k++) {
                            if (c.x[i][j][k] < 0.0) c.x[i][j][k] = 0.0;
                            const double dp_Pa = (i < im - 1)
                                ? (p_stat.x[i][j][k] - p_stat.x[i+1][j][k]) * 100.0
                                :  p_stat.x[i][j][k] * 100.0;
                            const double u_air = (dp_Pa > 0.0) ? dp_Pa / g : 0.0;
                            w_air += u_air;
                            const double c_ceiling = std::max(0.0, 1.0 - co2.x[i][j][k]);
                            if (c.x[i][j][k] > c_ceiling) {
                                q_clipped += (c.x[i][j][k] - c_ceiling) * u_air;
                                c.x[i][j][k] = c_ceiling;
                                n_ceiling++;
                            }
                        }
                if (w_air > 0.0) m_q_h2o_clipped += q_clipped / w_air;
                // A ceiling that keeps biting is a signal, not a solution: water is pumped
                // down out of the one condensing level by sedimentation and evaporates into
                // the superheated band below with no return path, so c there grows until
                // something stops it. Do not trust a run where this count stays large.
                if (n_ceiling > 0 && iter_n % diagnosticStride() == 0)
                    cout << "      AGCM: water-vapour ceiling c = 1 - co2 hit in "
                         << n_ceiling << " cells" << endl;
            }

            // ATHAD: CO2 is PROGNOSTIC and is no longer re-imposed here.
            //
            // What was here: ThermoAtm::co2Atmosphere(), which fills the whole field with
            // the uniform co2_0*co2_scale. Calling it inside the time loop overwrote the
            // transported field every iteration — rhs_co2 is built by RHS_Atm_Turb and
            // carried through all four RK4 stages, and every bit of it was discarded — so
            // the CO2 in the output was exactly constant (min = max = 0.205300 in every
            // cell) and the "prognostic CO2" of CLAUDE.md was diagnostic. It now runs at
            // initialisation only: the field starts well mixed, which is what a
            // non-condensable gas below the homopause should be, and is then transported.
            //
            // co2Column() takes over the call site. It fills co2_total (the column CO2 mass
            // path, previously declared and never written) and reports the drift of the
            // global mass-weighted mean q_CO2 — a conserved quantity here, since the model
            // has no CO2 source or sink, so any drift is transport error.
            ThermoAtm(*this).co2Column(iter_n % diagnosticStride() == 0);
            ThermoAtm(*this).waterBudget(iter_n % diagnosticStride() == 0, "post-moist");
            ThermoAtm(*this).densities();
            ThermoAtm(*this).forces();
            ThermoAtm(*this).standAtm_DewPoint_HumidRel();              // International Standard Atmosphere temperature profile, dew point temperature, relative humidity profile
            ThermoAtm(*this).waterVapourEvaporation();                  // correction of surface water vapour by evaporation
            ThermoAtm(*this).latentSensibleHeat();                      // latent and sensible heat
            ThermoAtm(*this).vegetationLand();                          // vegetation on land

            if(turb_model != "laminar" && !inviscid_phase) {
                TurbulenceAtm(*this).run();                             // update turbulence sources (skipped in laminar mode / inviscid spin-up)

                AtomUtils::damp_wiggles(tke,        &i_topography, true, true, true);
                AtomUtils::damp_wiggles(dis,        &i_topography, true, true, true);
                AtomUtils::damp_wiggles(nue,        &i_topography, true, true, true);
                AtomUtils::damp_wiggles(prod,       &i_topography, true, true, true);
                AtomUtils::damp_wiggles(tke_source, &i_topography, true, true, true);
                AtomUtils::damp_wiggles(dis_source, &i_topography, true, true, true);
            }

//            UtilsAtm(*this).valueLimitationAtm();                       // value limitation prevents local formation of NANs
        }  // iter_n % moist_stride == 0

        if(iter_n % momentum_stride == 0){
            BC_Atm(*this).bcRadius();                                   // extrapolation in i-direction alomg grid boundaries
            if(turb_model != "laminar") TurbulenceAtm(*this).apply_wall_bc();  // reassert ω_wall at i=0 after bcRadius cubic extrapolation
            BC_Atm(*this).bcTheta();                                    // extrapolation in j-direction alomg grid boundaries
            BC_Atm(*this).bcPhi();                                      // extrapolation in k-direction alomg grid boundaries

            BC_Atm(*this).bcScalarSurfSur();                            // scalar variable at surfaces extrapolated by von Neumann
            BC_Atm(*this).bcSolidGround();                              // values inside mountains

            // ATHAD: text diagnostics on their own cadence — every 10 iterations for a
            // short run, every 100 for a long one. The equatorial profile shows the
            // column structure; the level summary answers whether ANY column has gone
            // non-physical, which a single profile cannot.
            if(iter_n % diagnosticStride() == 0){
                ThermoAtm probe(*this);
                std::string tag = "iter " + std::to_string(iter_n);
                probe.printColumnProfile((jm - 1) / 2, ("equator, " + tag).c_str());
                probe.printLevelSummary(tag.c_str());
                probe.printPlanetaryBalance(tag.c_str());
            }

            if(iter_n % checkpoint == 0){
                // Psi is FILLED by write_meridional_streamfunction and READ by
                // print_min_max_atm, so the fill goes first. Ported with ATHAD's item 2:
                // with the old order the printed Psi extrema were one checkpoint stale and
                // the Psi FIELD was identically zero at the first checkpoint — which is
                // exactly what the field port measured here before this swap (max Psi =
                // min Psi = 0.000000 at iteration 20). The CSV was always correct; only the
                // printed extrema and the 3D field were behind.
                write_meridional_streamfunction(iter_n);   // Hadley/Ferrel cell strength (zonal-mean v + Ψ) per vtk checkpoint
                print_min_max_atm();
                UtilsAtm(*this).writeFile(bathymetry_name, output_path, false);
                cout << endl << "      AGCM: write_file in run_3D_loop atm ......................." << endl;
            }
        }  // iter_n % momentum_stride == 0

        // ---- zonal-mean v momentum budget: snapshot vbar before RK4, then difference
        // it across each step (RK4 physics, polar filter, orographic Shapiro, radial
        // Shapiro) to attribute the Hadley/Ferrel spin-down. Checkpoint iters only.
        const bool do_vbudget = (iter_n % checkpoint == 0);
        std::vector<std::vector<double> > vb_dyn, vb_polar, vb_orog, vb_radial, vb_stage, vb_prev;
        auto vb_diff = [&](std::vector<std::vector<double> >& dst){
            zonal_mean_v(vb_stage);
            for(int i = 0; i < im; i++)
                for(int j = 0; j < jm; j++){ dst[i][j] = vb_stage[i][j] - vb_prev[i][j]; vb_prev[i][j] = vb_stage[i][j]; }
        };
        // Parallel ZONAL-wind (w) budget — same per-stage differencing, for the trades/Walker.
        std::vector<std::vector<double> > wb_dyn, wb_polar, wb_orog, wb_radial, wb_stage, wb_prev;
        auto wb_diff = [&](std::vector<std::vector<double> >& dst){
            zonal_mean_w(wb_stage);
            for(int i = 0; i < im; i++)
                for(int j = 0; j < jm; j++){ dst[i][j] = wb_stage[i][j] - wb_prev[i][j]; wb_prev[i][j] = wb_stage[i][j]; }
        };
        if(do_vbudget){
            vb_stage.assign(im, std::vector<double>(jm, 0.0));
            vb_prev .assign(im, std::vector<double>(jm, 0.0));
            vb_dyn  .assign(im, std::vector<double>(jm, 0.0));
            vb_polar.assign(im, std::vector<double>(jm, 0.0));
            vb_orog .assign(im, std::vector<double>(jm, 0.0));
            vb_radial.assign(im, std::vector<double>(jm, 0.0));
            zonal_mean_v(vb_prev);   // vbar before RK4
            wb_stage.assign(im, std::vector<double>(jm, 0.0));
            wb_prev .assign(im, std::vector<double>(jm, 0.0));
            wb_dyn  .assign(im, std::vector<double>(jm, 0.0));
            wb_polar.assign(im, std::vector<double>(jm, 0.0));
            wb_orog .assign(im, std::vector<double>(jm, 0.0));
            wb_radial.assign(im, std::vector<double>(jm, 0.0));
            zonal_mean_w(wb_prev);   // wbar before RK4
        }

        // ubud_* is read by the ParaView/Results block at the TOP of a checkpoint iteration,
        // which runs before this RK4 — so capturing only on checkpoint iterations would leave
        // the written field one whole checkpoint interval stale (ATHAD item 2: identically
        // zero at the first checkpoint). Capturing on the PRE-checkpoint iteration as well
        // makes what the vtk carries one iteration old instead. vbudget/wbudget deliberately
        // keep the checkpoint-only condition: their CSV is written after this RK4, so it is
        // already current.
        ubudget_capture = do_vbudget || ((iter_n + 1) % checkpoint == 0);
        vbudget_capture = do_vbudget;     // have rhs_v store its per-term split this RK4 (turbulent path)
        wbudget_capture = do_vbudget;     // have rhs_w store its per-term split this RK4 (turbulent path)
        // Single dynamical core (2026-07-08): the turbulent RHS reduces to LAMINAR when
        // use_turbulence_model is false (molecular diffusion, zeroed turbulence) and to EULER
        // when diffusion_ramp=0 (inviscid spin-up). The former separate laminar
        // solveRungeKutta_Atmosphere / RHS_Atm.cpp path was dropped — inviscid is now an
        // independent switch (diffusion_ramp), decoupled from the turbulence selection.
        solveRungeKutta_Atmosphere_Turb();
        ubudget_capture = false;
        vbudget_capture = false;
        wbudget_capture = false;
        if(do_vbudget){ vb_diff(vb_dyn); wb_diff(wb_dyn); }   // RK4 net (PGF+Coriolis+advection+diffusion+drag+MC)

        // Polar zonal (φ) filter — root-cause stabiliser for the high-latitude blow-up.
        // Near the poles the zonal cell width rm·sinθ·dφ → 0, so the explicit zonal CFL
        // diverges and short φ-waves in u,v,w explode within a single RK4 step (the clamp
        // in bcSolidGround only bounds RK4 inputs, never the in-step growth). This damps
        // those short zonal wavelengths with a latitude-dependent pass count; runs every
        // iteration on the freshly-solved field, and the filtered values flow into
        // un,vn,wn via storeIntermediateData3D below.
        // start=30°, gain=3 strengthens the high-mid-latitude band (≈45–67°) where the
        // dry near-surface w-mode grows (seed 52°N, explosion 60°N/150°W Gulf of Alaska);
        // gives ~3 passes at 52°, ~6 at 60°, while ≤35° (subtropical jet) stays at 0.
        AtomUtils::polar_zonal_filter(u, the.z, &i_topography, 30.0, 12, 3.0);
        AtomUtils::polar_zonal_filter(v, the.z, &i_topography, 30.0, 12, 3.0);
        AtomUtils::polar_zonal_filter(w, the.z, &i_topography, 30.0, 12, 3.0);
        if(do_vbudget){ vb_diff(vb_polar); wb_diff(wb_polar); }   // polar zonal (φ) filter

        // Orographic Shapiro filter — same purpose as the polar filter, but for the
        // steep-orography CFL hot-spots (Andes/Atlas/NZ Alps) the latitude-keyed polar
        // filter doesn't reach. Damps the near-surface 2Δ velocity oscillation at cliffs
        // that the bounded pressure projection alone could not contain. passes=3 (was 1)
        // restores near-surface damping after n_layers_above was trimmed 12→3 to keep
        // Tibet aloft; heavier configs (n_layers_above=5, passes=6) tested 2026-05-31
        // SPREAD the Gulf-of-Alaska coastal blow-up across more cells (each pass averages
        // the ±100 cap into neighbours) — wrong tool, the runaway is forcing-driven, not
        // a 2Δ checkerboard. Use this gentle setting and attack the MC forcing instead.
        AtomUtils::orographic_shapiro_filter(u, i_topography, /*steep_threshold=*/2, /*n_layers_above=*/3, /*passes=*/3);
        AtomUtils::orographic_shapiro_filter(v, i_topography, /*steep_threshold=*/2, /*n_layers_above=*/3, /*passes=*/3);
        AtomUtils::orographic_shapiro_filter(w, i_topography, /*steep_threshold=*/2, /*n_layers_above=*/3, /*passes=*/3);
        if(do_vbudget){ vb_diff(vb_orog); wb_diff(wb_orog); }   // orographic Shapiro filter

        // NOTE (2026-06-08): an orographic Shapiro filter on the SCALAR fields was
        // tested here and REVERTED. Energy-budget tracing of the cliff cell (5,28,209,
        // Cook Inlet, i_topography=5) showed its iter-357 cooling is a real zonal-advection
        // 2Δ checkerboard (−w·∂t/∂φ), and filtering t DID hold that cell (t≈0.999 vs
        // crashing to 0.56). BUT the cliff is NOT the first-to-fail cell: the first NaN is
        // a u-velocity blow-up at the 62°N / k=0 upper-troposphere meridian seam (i≈36,
        // ~11 km), reached at iter 357 while the cliff is still only t≈0.75. The t-filter
        // perturbed the global field enough to tip that seam ~12 iters EARLIER (357→343/345
        // for passes 3/1), a net regression. The cliff is a secondary mode; the seam is the
        // real target. See [[project-cook-inlet-moist-runaway]].

        // Radial (i) de-checkerboarding. The 2Δ grid-scale checkerboard that seeds the
        // near-surface CFL blow-up has its purest component on the radial axis, and no other
        // filter (polar=k, orographic=j+k steep-only) reaches it. The vertical velocity u gets
        // a plain 1-2-1 pass — it carries little resolved shear and this annihilates the i=1
        // high-latitude 2Δ spike. For v and w a 4th-order Shapiro is used instead: it still
        // removes the 2Δ mode completely (CFL guard) but preserves the resolved vertical shear.
        // A 1-2-1 on v,w run EVERY iteration measurably erodes the jet's thermal-wind shear and
        // the Hadley/Ferrel cells' two-branch v structure, spinning the circulation down despite
        // a steady baroclinic (equator-pole) temperature gradient; the 4th-order filter does not.
        // ⚠️ A/B KNOB 2026-07-21 (ATM_RADIAL_SHAPIRO_STRENGTH, default 1.0 = bit-identical).
        // The column-integrated w-momentum budget (scratchpad/column_radial_budget.py) showed
        // these radial Shapiro passes on v,w are the DOMINANT net sink of extratropical-jet
        // momentum (~85% of Mdot_net at 30-60°, ~6x the resolved dynamics; the ocean has no
        // radial velocity Shapiro at all). The 4th-order _ho variant was already chosen for v,w
        // to spare the resolved shear, but it still dominates. This scales the filter STRENGTH so
        // we can A/B whether easing it slows the jet spin-down WITHOUT re-triggering the radial 2Δ
        // checkerboard / near-surface CFL blow-up it guards against. 0 = filter off (CFL risk).
        static const double radial_shapiro_strength = [](){
            const char* e = getenv("ATM_RADIAL_SHAPIRO_STRENGTH"); return e ? atof(e) : 1.0; }();
        AtomUtils::radial_shapiro_filter   (u, i_topography, /*passes=*/2, radial_shapiro_strength);
        AtomUtils::radial_shapiro_filter_ho(v, i_topography, /*passes=*/2, radial_shapiro_strength);
        AtomUtils::radial_shapiro_filter_ho(w, i_topography, /*passes=*/2, radial_shapiro_strength);
        if(do_vbudget){ vb_diff(vb_radial); wb_diff(wb_radial);   // radial (vertical) Shapiro filter  [prime spin-down suspect]
            write_v_momentum_budget(iter_n, vb_dyn, vb_polar, vb_orog, vb_radial);
            write_w_momentum_budget(iter_n, wb_dyn, wb_polar, wb_orog, wb_radial);
        }

        // Radial de-checkerboarding for the SCALARS at orographic columns. The earlier
        // velocity radial filters do nothing for t/c/cloud/ice, which carry their own
        // single-cell 2Δz vertical checkerboard at the orographic condensation level — the
        // dry seed of the NZ-Alps blow-up (cold cell ~surface+7 over the steep Alps,
        // [[project_nz_alps_t_checkerboard]]). Extremum-gated + steep-column-confined so it
        // touches only the actual checkerboard cells and leaves the (healthy) bulk scalar
        // field bit-unchanged — avoiding the global perturbation that reverted the 2026-06-08
        // horizontal scalar filter. n_layers_above=10 reaches the surf+7 mode at NZ (surf=6).
        AtomUtils::orographic_radial_shapiro_filter(t,     i_topography, /*steep=*/2, /*n_layers_above=*/10, /*passes=*/2);
        AtomUtils::orographic_radial_shapiro_filter(c,     i_topography, /*steep=*/2, /*n_layers_above=*/10, /*passes=*/2);
        AtomUtils::orographic_radial_shapiro_filter(cloud, i_topography, /*steep=*/2, /*n_layers_above=*/10, /*passes=*/2);
        AtomUtils::orographic_radial_shapiro_filter(ice,   i_topography, /*steep=*/2, /*n_layers_above=*/10, /*passes=*/2);
        // ==============================================================================================

        // Soft steep-massif field smoothing (2026-06-10, user request). The init-time terrain
        // massif smoothing fixed the SURVIVAL (steep ramparts were the runaway driver), but the
        // residual steepest columns keep a locally unphysical field aloft (user: saturated p_dyn
        // block + bad u/v/w above the Himalaya, while smooth-contour South Pole behaves well).
        // Gently smooth u,v,w,p_dyn ONLY over the HIGH+STEEP contour columns (same narrow mask
        // as the terrain smoothing), reaching deep up the troposphere — so the near-orography
        // field is cleaned WITHOUT touching the broad real flow, plains/coasts, or the i=0
        // ocean vortices. Single gentle pass. See [[project-upper-velocity-secular-growth]].
        AtomUtils::steep_massif_field_smoothing(u,     i_topography, 6, 3, /*deep=*/30, /*passes=*/1);
        AtomUtils::steep_massif_field_smoothing(v,     i_topography, 6, 3, /*deep=*/30, /*passes=*/1);
        AtomUtils::steep_massif_field_smoothing(w,     i_topography, 6, 3, /*deep=*/30, /*passes=*/1);
        AtomUtils::steep_massif_field_smoothing(p_dyn, i_topography, 6, 3, /*deep=*/30, /*passes=*/1);

        // NOTE (2026-06-09): a meridional (j/θ) Shapiro filter was added here and REMOVED.
        // It bought only +8 iters against the iter-483 NaN (a symptom-only palliative — the
        // real cause is the pressure/divergence loop at the Tian Shan, see
        // [[project-upper-velocity-secular-growth]]) AND it smoothed away the emergent i=0
        // vortices that are a healthy-run validation signal ([[project-emergence-validation]]).
        // Net negative; do not re-add. Radial (i) + polar zonal (k) filtering is retained.

        // (Upper-troposphere Rayleigh sponge re-tested 2026-06-09 with i_start=30,
        // alpha_max=0.2 against the iter-483 upper-level velocity blow-up — REVERTED, net
        // negative. It did NOT touch the secular growth (diag62N max|u| peaks at i=25,
        // BELOW i_start=30, and was identical with/without the sponge: 0.230 vs 0.231 at
        // iter 320), and relaxing the zonal anomaly aloft perturbed the sensitive k=0 seam
        // mode → crash moved 483→479 and the seed RELOCATED from the 40°N/62°E Pamir column
        // to the 55°N/0°E meridian seam at i=36 (10.9 km). Confirms the rupture point is
        // whack-a-mole: the energy source is the secular velocity growth in the 55–62°N
        // band centred at i≈25; the seam is just the weakest point where it erupts. The
        // existing bcPhi 2-pass seam Shapiro (BC_Atm.h) already covers the seam and is not
        // enough. Attack the growth source, not the rupture point.)
        // AtomUtils::upper_rayleigh_sponge(u, 30, 0.2);
        // AtomUtils::upper_rayleigh_sponge(v, 30, 0.2);
        // AtomUtils::upper_rayleigh_sponge(w, 30, 0.2);

        // Near-surface coastal velocity HARD CAP — last-resort ceiling for the dry-seeded
        // velocity runaway at steep coasts (Gulf-of-Alaska 60°N/151°W). The slope cap
        // skips ocean↔land cliffs and bcVelSurfSur is init-only; without an in-loop RHS
        // friction term, coastal pressure-gradient force drives a linear velocity growth
        // that saturates at the ±100 m/s BC cap. At the bottom `layers` air cells of a
        // coastal/sloped column, if |u_nondim| > 3.0 (=30 m/s), clamp to ±3.0. Rate=1.0
        // selects the hard-clamp branch (Rayleigh form rounds 6-7 was too weak — per-step
        // RHS forcing >> %-damping). layers=8: the runaway is at i≈5 (298 m) above the
        // OCEAN side of the Cook Inlet coast, not at the land cliff itself — layers=3
        // (round 8) only reached i=0..2 over ocean columns and missed the actual blow-up.
        // 8 covers i_topo..i_topo+7 = 0..7 (≈440 m) above ocean / 298 m..1.5 km over land.
        AtomUtils::coastal_velocity_sponge(u, i_topography, 3.0, 1.0, 12);
        AtomUtils::coastal_velocity_sponge(v, i_topography, 3.0, 1.0, 12);
        AtomUtils::coastal_velocity_sponge(w, i_topography, 3.0, 1.0, 12);

        UtilsAtm(*this).findResiduumAtm();

        UtilsAtm(*this).storeIntermediateData3D(1.0);

        // Interactive radiation, per mode: (A) re-solve MLR on the current field to refresh the
        // t_eq relaxation target every teq_refresh_stride iters; (B) apply MLR as a direct
        // radiative heating on T every iteration. Either way the radiation tracks the evolving
        // climate (water-vapour feedback, CO2) instead of being frozen at the initial state.
        if (radiation_mode == 1 && iter_n % teq_refresh_stride == 0)
            refresh_radiative_teq();                            // A: refresh MLR -> t_eq target
        else if (radiation_mode == 2) {
            apply_radiative_heating();                          // B: MLR direct heating on T
            // MLR has just rebuilt albedo.y from the current condensate, so this is the
            // point at which the model's own albedo is known. Relax t_skin toward the
            // energy-balance value it implies; densities() picks the new value up on the
            // next iteration when it re-imposes the isothermal top.
            updateSkinTemperature(iter_n % diagnosticStride() == 0);
        }
        else if (radiation_mode == 4 && iter_n % teq_refresh_stride == 0)
            refresh_radiation_diag();                           // 4: refresh radiation.x diag (t unchanged)
        else if (radiation_mode == 5) {
            if (iter_n % teq_refresh_stride == 0) {
                apply_co2_perturbation(false);                  // 5: WV feedback — rebuild deep t_eq perturbation on the warmed field
                cloud_radiation_diag();                         // 5: refresh radiation.x with clouds + print cloud LW effect
            }
            apply_teq_relaxation();                             // 5: strong relax t -> Scotese+CO2 perturbation
        }

        ThermoAtm(*this).printDataAtm();

        // Convergence monitor: sample the slow integral metrics and report secular drift.
        if(conv_stride > 0 && total_iter_count > 0 && total_iter_count % conv_stride == 0){
            double dT_K = 0.0, dKE_abs = 0.0;
            const double mT  = domain_mean_T();
            const double mKE = domain_mean_KE();
            const double dT  = conv_drift(conv_histT,  mT,  dT_K);
            const double dKE = conv_drift(conv_histKE, mKE, dKE_abs);
            const bool warming = (dT < 0.0) || (dKE < 0.0);
            const bool converged = !warming && dT_K < conv_tol_T_K && dKE < conv_tol_pct;
            cout << "      [conv] iter " << total_iter_count
                 << "  T=" << mT << " K  KE=" << mKE << " m2/s2";
            if(warming) cout << "  (warming up)";
            else cout << "  drift: T=" << dT_K << " K (" << dT << "%)  KE=" << dKE
                      << "%/window" << (converged ? "  -> CONVERGED (both < tol)" : "");
            cout << std::endl;
            std::ofstream cf(output_path + "/convergence.csv", std::ios::app);
            if(cf) cf << total_iter_count << "," << mT << "," << (warming?0.0:dT) << ","
                      << (warming?0.0:dT_K) << "," << mKE << "," << (warming?0.0:dKE) << ","
                      << (converged?1:0) << "\n";
        }

        // Checkpoint the full 3D state once total_iter_count reaches the requested iter.
        if(checkpoint_save_iter >= 0 && total_iter_count == checkpoint_save_iter)
            save_state(checkpoint_save_iter, Ma);

        // Periodic restart checkpoints every 100 iters, but ONLY once the run is clean
        // (no non-finite cell in any serialized prognostic field) so we never write a
        // NaN state that would later poison restart_from_iter. Files land at
        // output_path/atm_restart_<Ma>Ma_<total_iter_count>.bin (Ma+iter stamped, so a
        // paleo chain of several time-slices never clobbers another slice's restart).
        {
            const int restart_save_stride = restart_stride;
            if(restart_save_stride > 0 && total_iter_count > 0
               && total_iter_count % restart_save_stride == 0){
                bool clean = true;
                for(Array* a : restart_arrays()){
                    for(int i = 0; i < im && clean; i++)
                        for(int j = 0; j < jm && clean; j++)
                            for(int k = 0; k < km && clean; k++)
                                if(!AtomUtils::is_finite_safe(a->x[i][j][k])) clean = false;
                }
                if(clean)
                    save_state(total_iter_count, Ma);
                else
                    cout << "      AGCM: restart checkpoint SKIPPED at iter "
                         << total_iter_count << " — non-finite cell present (not clean)" << endl;
            }
        }

    }  // end iter_n


    cout << endl << "      AGCM: run_3D_loop atm ended ..........................." << endl;
}
/*
*
*/
void cAtmosphereModel::Run(){

    std::time_t Run_start;

    struct tm * timeinfo_begin;

    std::time(&Run_start);

    timeinfo_begin = std::localtime(&Run_start);

    std::cout << std::endl << std::endl;
    std::cout << " ... AGCM: time and date at run time begin:   " 
        << std::asctime(timeinfo_begin);
    std::cout << std::endl << std::endl;




    mkdir(output_path.c_str(), 0777);

    if(!has_printed_welcome_msg)  print_welcome_msg();

    for(int i = time_start; i <= time_end; i+=time_step){
        RunTimeSlice(i);
    }

    print_final_remarks();




    std::time_t Run_end;

    struct tm * timeinfo_end;

    std::time(&Run_end);

    timeinfo_end = std::localtime(&Run_end);

    std::cout << " ... AGCM: time and date at run time end:   " 
        << std::asctime(timeinfo_end) << std::endl;

    int Run_total = Run_end - Run_start;
    int Run_total_minutes = Run_total/60;
    int Run_total_hours = Run_total_minutes/60;

    std::cout <<  " ... AGCM: computer time needed for Ma = " << time_start 
        << " to Ma = " << time_end << " in " << time_step << " Ma steps:     "
        << Run_total << " seconds" << std::endl
        << " ... compares to:" << std::endl << std::endl
        << setw(30) << setfill(' ') << Run_total_hours 
        << " hours" << std::endl
        << setw(30) << setfill(' ') << Run_total_minutes 
        << " minutes" << std::endl
        << setw(30) << setfill(' ') << Run_total%60 
        << " seconds" << std::endl << std::endl
        << std::endl;
}
/*
*
*/

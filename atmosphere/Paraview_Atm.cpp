/*
 * Atmosphere General Circulation Modell (AGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to write sequel, transfer and paraview files
*/

#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cstdlib>

#include "cAtmosphereModel.h"
#include "Utils.h"
#include "SaturationH2O.h"
#include "MixtureAtm.h"

using namespace std;
using namespace AtomUtils;

namespace ParaViewAtm{
    // Bit-level non-finite check.
    // The Makefile uses -ffast-math, which implies -ffinite-math-only; under that flag the
    // compiler assumes no NaN/Inf can occur and optimizes std::isfinite(v) into constant true.
    // That silently disables this filter and lets the literal text "nan" reach the VTK file,
    // tripping vtkDataReader with "Unsupported point attribute type: nan". Reading the IEEE-754
    // exponent bits via memcpy is invisible to the math optimizer: NaN and ±Inf both have all
    // 11 exponent bits set, so masking with 0x7FF0…ULL detects them regardless of -ffast-math.
    inline double safe_val(double v) {
        std::uint64_t bits;
        std::memcpy(&bits, &v, sizeof(bits));
        return ((bits & 0x7FF0000000000000ULL) == 0x7FF0000000000000ULL) ? 0.0 : v;
    }

    void dump_array(const string &name, Array &a, double multiplier, ofstream &f){
        f <<  "    <DataArray type=\"Float32\" Name=\"" << name << "\" format=\"ascii\">\n";
        for(int k = 0; k < a.km; k++){
            for(int j = 0; j < a.jm; j++){
                for(int i = 0; i < a.im; i++){
                    f << safe_val(a.x[i][j][k] * multiplier) << "\n";
                }
                f << "\n";
            }
            f << "\n";
        }
        f << "\n";
        f << "    </DataArray>\n";
        f.clear();
    }
/*
 * 
*/
    void dump_radial(const string &desc, Array &a, double multiplier, int i, ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << "\n";
        f << "LOOKUP_TABLE default" << "\n";
        for(int j = 0; j < a.jm; j++){
            for(int k = 0; k < a.km; k++){
                f << safe_val(a.x[i][j][k] * multiplier) << "\n";
            }
        }
    }
/*
 * 
*/
    void dump_radial_2d(const string &desc, Array_2D &a, double multiplier, ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << "\n";
        f << "LOOKUP_TABLE default" << "\n";
        for(int j = 0; j < a.jm; j++){
            for(int k = 0; k < a.km; k++){
                f << safe_val(a.y[j][k] * multiplier) << "\n";
            }
        }
    }
/*
 * 
*/
    void dump_zonal(const string &desc, Array &a, double multiplier, int k, ofstream &f){
        f <<  "SCALARS " << desc << " float " << 1 << "\n";
        f <<  "LOOKUP_TABLE default" << "\n";
        for(int i = 0; i < a.im; i++){
            for(int j = 0; j < a.jm; j++){
                f << safe_val(a.x[i][j][k] * multiplier) << "\n";
            }
        }
    }
/*
 * 
*/
    void dump_longal(const string &desc, Array &a, double multiplier, int j, ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << "\n";
        f << "LOOKUP_TABLE default" << "\n";
        for(int i = 0; i < a.im; i++){
            for(int k = 0; k < a.km; k++){
                f << safe_val(a.x[i][j][k] * multiplier) << "\n";
            }
        }
    }
}
/*
 * FluidPhase (2026-11-28; user: "there are oceans in the family but I can't see them in Paraview, add a type of
 * species ... that indicates the fluid phase and its location"). OUTPUT ONLY: computed from the state at each
 * write, read by no routine. The ocean / melt of this family is a boundary condition on the ground row
 * i = i_topography (= 0, invariant 1: no topography), never a field, so nothing in the output showed it.
 * One categorical value per point, from the model's own water physics (SaturationH2O.h):
 *   ground row (the surface, i == i_topography):
 *     10  liquid-water ocean   T_TRIPLE <= T_s < T_crit and p_s > p_sat(T_s)
 *     11  water ice            T_s < T_TRIPLE
 *     12  silicate melt        no liquid water possible (T_s >= T_crit or p_s <= p_sat) and T_s >= T_solidus
 *     13  solid rock           neither
 *   air (i > i_topography):
 *     0  vapour only           1  + liquid cloud      2  + ice / graupel cloud      3  + both (mixed)
 *     4  supercritical         T >= T_crit (isSupercritical: no liquid can form; the cloud test is skipped)
 *     5  condensate ON THE CAP cloud >= 0.99 cloud_cap or ice >= 0.99 cloud_cap (the per-field ceiling of
 *                              SaturationAdjustment.h; takes precedence over 1-3). Not an ocean -- the densest
 *                              cell is ~0.1 % liquid by volume -- but where condensate the model would otherwise
 *                              hold is clipped (register row COND-CAP-COLUMN).
 *   inside the ground (i < i_topography; never reached without topography): 13.
 * Condensate counts above 1e-9 kg/kg. T_solidus = ATM_PHASE_T_SOLIDUS [K], default 1393 K (dry peridotite at
 * low pressure, ~1120 C); an assumption, not a model constant -- nothing else reads it.
 * In ParaView: Threshold FluidPhase 10..13 shows the surface and what it is; 4 the supercritical column.
 * ATM_FLUID_PHASE=0 does not write the array (the files are then as before the field existed).
 */
static bool fluid_phase_on(){
    static const bool v = [](){ const char* e = getenv("ATM_FLUID_PHASE"); return !(e && atoi(e) == 0); }();
    return v;
}
static Array& fluid_phase_field(cAtmosphereModel& m){
    static Array fp;
    if(fp.im != m.im || fp.jm != m.jm || fp.km != m.km) fp.initArray(m.im, m.jm, m.km, 0.0);
    static const double T_sol = [](){ const char* e = getenv("ATM_PHASE_T_SOLIDUS");
                                      return (e && atof(e) > 0.0) ? atof(e) : 1393.0; }();
    constexpr double q_min = 1.0e-9;
    for(int j = 0; j < m.jm; j++)
        for(int k = 0; k < m.km; k++){
            const int ig = m.i_topography[j][k];
            for(int i = 0; i < m.im; i++){
                const double T = m.t.x[i][j][k] * m.t_0;               // K
                const double p = m.p_stat.x[i][j][k];                  // hPa
                double code;
                if(i < ig) code = 13.0;
                else if(i == ig){
                    if(T < SaturationH2O::T_TRIPLE) code = 11.0;
                    else {
                        const double e = SaturationH2O::saturationPressure(T);   // NO_SATURATION above T_crit
                        if(e != SaturationH2O::NO_SATURATION && e < p) code = 10.0;
                        else code = (T >= T_sol) ? 12.0 : 13.0;
                    }
                }
                else if(SaturationH2O::isSupercritical(T)) code = 4.0;
                else {
                    const bool liq = m.cloud.x[i][j][k] > q_min;
                    const bool sol = m.ice.x[i][j][k] + m.gr.x[i][j][k] > q_min;
                    const double cap = 0.99 * m.cloud_cap;
                    if(m.cloud_cap > 0.0 && (m.cloud.x[i][j][k] >= cap || m.ice.x[i][j][k] >= cap)) code = 5.0;
                    else code = (liq && sol) ? 3.0 : liq ? 1.0 : sol ? 2.0 : 0.0;
                }
                fp.x[i][j][k] = code;
            }
        }
    return fp;
}
/*
 * The eight species in the 3D files (2026-11-28; user: "except CO2 and H2O the other species are ... constant and
 * not represented in the RK scheme ... are they so invisible"). OUTPUT ONLY. The radial / zonal / longitudinal
 * slices have written all eight since README item 61 (`THE EIGHT SPECIES, ALL FROM split()` below); the panorama and
 * the sphere had only WaterVapour and the raw CO2 tracer, so N2, CH4, NH3, H2, CO, SO2 were invisible in 3D. Same
 * names, same source, same unit as the slices: AtmMixture::split() with the cell's condensate, H2O = q_v, CO2 = q_c
 * (carrier-renormalised), the rest q_b * f_bg[n] -- MASS FRACTIONS, dimensionless. The six background gases are not
 * prognostic: they are a fixed share of the background and move only with q_b, i.e. with the water and the CO2.
 * ATM_SPECIES_3D=0 does not write them.
 */
static bool species_3d_on(){
    static const bool v = [](){ const char* e = getenv("ATM_SPECIES_3D"); return !(e && atoi(e) == 0); }();
    return v;
}
static double species_q(cAtmosphereModel& m, int n, int i, int j, int k){
    double q_v, q_c, q_b;
    const double q_l = m.cloud.x[i][j][k] + m.ice.x[i][j][k] + m.gr.x[i][j][k];
    AtmMixture::split(m.c.x[i][j][k], m.co2.x[i][j][k], q_v, q_c, q_b, q_l);
    return (n == 0) ? q_v : (n == 1) ? q_c : q_b * m.m_comp.f_bg[n - 2];
}
/*
 * 
*/
void cAtmosphereModel::paraview_panorama_vts(string &Name_Bathymetry_File, int n){

    using namespace ParaViewAtm;

    string Atmosphere_panorama_vts_File_Name = output_path + "/"
        + Name_Bathymetry_File + "_Atm_panorama_" + std::to_string(n) + ".vts";
    ofstream Atmosphere_panorama_vts_File;

    Atmosphere_panorama_vts_File.precision(8);
    Atmosphere_panorama_vts_File.setf(ios::fixed);

    Atmosphere_panorama_vts_File.open(Atmosphere_panorama_vts_File_Name);

    if(!Atmosphere_panorama_vts_File.is_open()){
        cerr << "ERROR: could not open panorama_vts file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }

    ostringstream buf;
    buf.precision(8);
    buf.setf(ios::fixed);

    buf << "<?xml version=\"1.0\"?>\n\n"
        << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n\n"
        << " <StructuredGrid WholeExtent=\"" << 1 << " " << im << " " << 1 << " " << jm << " " << 1 << " " << km << "\">\n\n"
        << "  <Piece Extent=\"" << 1 << " " << im << " " << 1 << " " << jm << " " << 1 << " " << km << "\">\n\n"
        << "   <PointData Vectors=\"Velocity MagneticField\" Scalars=\"Topography u-component v-component w-component Temperature CondensationTemp EvaporationTemp Epsilon_3D PressureDynamic PressureStatic WaterVapour CloudWater CloudIce CO2-Concentration Q_Latent Rain RainSuper Ice PrecipitationRain PrecipitationSnow PrecipitationConv Updraft Downdraft\">\n\n"
        << "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n\n";

    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                buf << safe_val(u.x[i][j][k]) << " "
                    << safe_val(v.x[i][j][k]) << " " << safe_val(w.x[i][j][k]) << '\n';
            }
            buf << "\n\n";
        }
        buf << "\n\n";
    }

    buf << "\n\n    </DataArray>\n\n";
    Atmosphere_panorama_vts_File << buf.str();
    buf.str("");
    buf.clear();



    dump_array("Topography", h, 1.0, Atmosphere_panorama_vts_File);
    if(fluid_phase_on()) dump_array("FluidPhase", fluid_phase_field(*this), 1.0, Atmosphere_panorama_vts_File);

    dump_array("u-component", u, u_0, Atmosphere_panorama_vts_File);
    dump_array("v-component", v, u_0, Atmosphere_panorama_vts_File);
    dump_array("w-component", w, u_0, Atmosphere_panorama_vts_File);


    buf << "    <DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">\n\n";
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                buf << safe_val(t.x[i][j][k] * t_0 - t_0) << '\n';
            }
            buf << "\n\n";
        }
        buf << "\n\n";
    }
    buf << "\n\n    </DataArray>\n\n";
    Atmosphere_panorama_vts_File << buf.str();
    buf.str("");
    buf.clear();


    dump_array("Radiation", radiation, 1.0, Atmosphere_panorama_vts_File);

    dump_array("WaterVapour", c, 1e3, Atmosphere_panorama_vts_File);
    dump_array("CloudWater", cloud, 1e3, Atmosphere_panorama_vts_File);
    dump_array("CloudIce", ice, 1e3, Atmosphere_panorama_vts_File);
    dump_array("CloudGraupel", gr, 1e3, Atmosphere_panorama_vts_File);
    // Moist source/sink terms. S_c and S_r are the rhs_t latent-heat drivers
    // (RHS_Atm.cpp:397, coeff_energy*lv*(S_c+S_r)) that couple precip -> temperature
    // -> buoyancy -> velocity; dump them so the NE-Pacific (j=37,k=229) precip<->u
    // coupling is inspectable in ParaView. Were previously commented out AND all three
    // mistakenly dumped S_r — fixed to S_v/S_c/S_r.
    dump_array("S_c_c", S_c_c, 1e3, Atmosphere_panorama_vts_File);
    dump_array("S_v", S_v, 1e3, Atmosphere_panorama_vts_File);
    dump_array("S_c", S_c, 1e3, Atmosphere_panorama_vts_File);
    dump_array("S_r", S_r, 1e3, Atmosphere_panorama_vts_File);

    dump_array("Precipitation", Precipitation, 8.64e4, Atmosphere_panorama_vts_File);
    dump_array("PrecipitationRain", P_rain, 8.64e4, Atmosphere_panorama_vts_File);
    dump_array("PrecipitationSnow", P_snow, 8.64e4, Atmosphere_panorama_vts_File);
    dump_array("PrecipitationGraupel", P_graupel, 8.64e4, Atmosphere_panorama_vts_File);

//    dump_array("PressureStatic", p_stat, 1.0, Atmosphere_panorama_vts_File);
//    dump_array("PressureDynamic", p_dyn, p_0, Atmosphere_panorama_vts_File);
//    dump_array("r_humid", r_humid, 1.0, Atmosphere_panorama_vts_File);

    dump_array("CO2-Concentration", co2, co2_0, Atmosphere_panorama_vts_File);
    if(species_3d_on())                                                 // ATM_SPECIES_3D: the eight, as in the slices
        for(int n = 0; n < 8; n++){
            Atmosphere_panorama_vts_File << "    <DataArray type=\"Float32\" Name=\"" << AtmMixture::SPECIES_NAMES()[n]
                                         << "\" format=\"ascii\">\n";
            for(int k = 0; k < km; k++)
                for(int j = 0; j < jm; j++)
                    for(int i = 0; i < im; i++)
                        Atmosphere_panorama_vts_File << safe_val(species_q(*this, n, i, j, k)) << "\n";
            Atmosphere_panorama_vts_File << "    </DataArray>\n";
        }
//    dump_array("Q_Latent", Q_Latent, 1.0, Atmosphere_panorama_vts_File);
//    dump_array("Q_Radiation", radiation, 1.0, Atmosphere_panorama_vts_File);
//    dump_array("Q_Sensible", Q_Sensible, 1.0, Atmosphere_panorama_vts_File);

/*
    dump_array("u_u", u_u, 1.0, Atmosphere_panorama_vts_File);
    dump_array("u_d", u_d, 1.0, Atmosphere_panorama_vts_File);
    dump_array("M_u", M_u, 1.0, Atmosphere_panorama_vts_File);
    dump_array("M_d", M_d, 1.0, Atmosphere_panorama_vts_File);
*/
/*
    dump_array("q_v_u", q_v_u, 1e3, Atmosphere_panorama_vts_File);
    dump_array("q_v_d", q_v_d, 1e3, Atmosphere_panorama_vts_File);
    dump_array("q_c_u", q_c_u, 1e3, Atmosphere_panorama_vts_File);
*/
/*
    dump_array("g_p", g_p, 1e3, Atmosphere_panorama_vts_File);
    dump_array("c_u", c_u, 1e3, Atmosphere_panorama_vts_File);
    dump_array("e_d", e_d, 1e3, Atmosphere_panorama_vts_File);
    dump_array("e_p", e_p, 1e3, Atmosphere_panorama_vts_File);
*/
    dump_array("PrecipitationConv", P_conv, 8.64e4, Atmosphere_panorama_vts_File);

    if(turb_model != "laminar"){
        dump_array("TKE",           tke,        1.0, Atmosphere_panorama_vts_File);
        dump_array("Dissipation",   dis,        1.0, Atmosphere_panorama_vts_File);
        dump_array("EddyViscosity", nue,        1.0, Atmosphere_panorama_vts_File);
        dump_array("TurbProd",      prod,       1.0, Atmosphere_panorama_vts_File);
        dump_array("TKE_Source",    tke_source, 1.0, Atmosphere_panorama_vts_File);
        dump_array("Dis_Source",    dis_source, 1.0, Atmosphere_panorama_vts_File);
    }

    buf << "   </PointData>\n\n"
        << "   <Points>\n\n"
        << "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n\n";

    double dx = 0.1;
    double dy = 0.1;
    double dz = 0.1;

    for(int k = 0; k < km; k++){
        double z = k * dz;
        for(int j = 0; j < jm; j++){
            double y = j * dy;
            for(int i = 0; i < im; i++){
                buf << i * dx << " " << y << " " << z << '\n';
            }
            buf << "\n\n";
        }
        buf << "\n\n";
    }

    buf << "    </DataArray>\n\n"
        << "   </Points>\n\n"
        << "  </Piece>\n\n"
        << " </StructuredGrid>\n\n"
        << "</VTKFile>\n\n";
    Atmosphere_panorama_vts_File << buf.str();
    Atmosphere_panorama_vts_File.close();

    cout << "   File:  " << Atmosphere_panorama_vts_File_Name
        << "  has been written\n";
}
/*
*
*/
void cAtmosphereModel::paraview_sphere_vts(string &Name_Bathymetry_File, int n){
    // Hadean_Atm_sphere_<n>.vts (2026-11-27, giants' GIANTS-SPHERE-VTS fix ported; written with the panorama,
    // ATM_SPHERE=0 switches it off). Two changes against the writer that had been commented out:
    //  - the pole rows j = 0 and j = jm-1 are NOT written: there sin(the) = 0, all km longitudes are one point,
    //    and ParaView's Contour cannot invert the collapsed cells ("Cannot compute gradient of grid"). The
    //    sphere spans 1..179 deg (extent j = 2..jm-1); the pole values stay in the panorama .vts.
    //  - the Cartesian velocity is built in locals. The old writer stored it in aux_u / aux_v / aux_w, which the
    //    model uses as work arrays, so a write could have changed the run.
    using namespace ParaViewAtm;
    const string fname = output_path + "/" + Name_Bathymetry_File + "_Atm_sphere_" + std::to_string(n) + ".vts";
    ofstream f;
    f.precision(8);
    f.setf(ios::fixed);
    f.open(fname);
    if(!f.is_open()){
        cerr << "ERROR: could not open paraview_sphere_vts file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    const int j0 = 1, j1 = jm - 1;                                  // rows written: j0 .. j1-1
    ostringstream buf;
    buf.precision(8);
    buf.setf(ios::fixed);
    buf << "<?xml version=\"1.0\"?>\n\n"
        << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n\n"
        << " <StructuredGrid WholeExtent=\"" << 1 << " " << im << " " << j0 + 1 << " " << j1 << " " << 1 << " " << km << "\">\n\n"
        << "  <Piece Extent=\"" << 1 << " " << im << " " << j0 + 1 << " " << j1 << " " << 1 << " " << km << "\">\n\n"
        << "   <PointData Vectors=\"Velocity\" Scalars=\"Temperature\">\n\n"
        << "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n";
    // Cartesian components (m/s) kept for the three scalar arrays below
    const size_t npt = (size_t)im * (j1 - j0) * km;
    std::vector<float> cx(npt), cy(npt), cz(npt);
    size_t q = 0;
    for(int k = 0; k < km; k++){
        const double sinphi = sin(phi.z[k]), cosphi = cos(phi.z[k]);
        for(int j = j0; j < j1; j++){
            const double sinthe = sin(the.z[j]), costhe = cos(the.z[j]);
            for(int i = 0; i < im; i++, q++){
                const double ui = u.x[i][j][k], vi = v.x[i][j][k], wi = w.x[i][j][k];
                const double au = u_0 * (sinthe * cosphi * ui + costhe * cosphi * vi - sinphi * wi);
                const double av = u_0 * (sinthe * sinphi * ui + costhe * sinphi * vi + cosphi * wi);
                const double aw = u_0 * (costhe * ui - sinthe * vi);
                cx[q] = (float)safe_val(au); cy[q] = (float)safe_val(av); cz[q] = (float)safe_val(aw);
                buf << cx[q] << " " << cy[q] << " " << cz[q] << '\n';
            }
        }
    }
    buf << "    </DataArray>\n";
    f << buf.str();
    buf.str("");
    buf.clear();

    auto dump_ = [&](const char* name, Array& a, double mult, double offs){
        f << "    <DataArray type=\"Float32\" Name=\"" << name << "\" format=\"ascii\">\n";
        for(int k = 0; k < km; k++)
            for(int j = j0; j < j1; j++)
                for(int i = 0; i < im; i++)
                    f << safe_val(a.x[i][j][k] * mult + offs) << '\n';
        f << "    </DataArray>\n";
    };
    auto dump_v = [&](const char* name, const std::vector<float>& c){
        f << "    <DataArray type=\"Float32\" Name=\"" << name << "\" format=\"ascii\">\n";
        for(size_t p = 0; p < c.size(); p++) f << c[p] << '\n';
        f << "    </DataArray>\n";
    };
    dump_("Temperature", t, t_0, -t_0);                             // deg C, as the panorama
    dump_("Bathymetry", h, 1.0, 0.0);
    if(fluid_phase_on()) dump_("FluidPhase", fluid_phase_field(*this), 1.0, 0.0);
    dump_v("u-component", cx);                                      // Cartesian x, y, z in m/s
    dump_v("v-component", cy);
    dump_v("w-component", cz);
    dump_("PressureDynamic", p_dyn, p_0, 0.0);
    dump_("Pressure_bar", p_stat, 1.0e-3, 0.0);
    dump_("WaterVapour", c, 1e3, 0.0);                              // g/kg
    dump_("CloudWater", cloud, 1e3, 0.0);
    dump_("CloudIce", ice, 1e3, 0.0);
    dump_("CloudGraupel", gr, 1e3, 0.0);
    dump_("Density", r_humid, 1.0, 0.0);                           // kg/m3
    dump_("Precipitation", Precipitation, 8.64e4, 0.0);            // mm/d
    dump_("PrecipitationRain", P_rain, 8.64e4, 0.0);
    dump_("PrecipitationSnow", P_snow, 8.64e4, 0.0);
    if(species_3d_on())                                                 // ATM_SPECIES_3D: the eight, as in the slices
        for(int n = 0; n < 8; n++){
            f << "    <DataArray type=\"Float32\" Name=\"" << AtmMixture::SPECIES_NAMES()[n] << "\" format=\"ascii\">\n";
            for(int k = 0; k < km; k++)
                for(int j = j0; j < j1; j++)
                    for(int i = 0; i < im; i++)
                        f << safe_val(species_q(*this, n, i, j, k)) << '\n';
            f << "    </DataArray>\n";
        }

    buf << "   </PointData>\n\n"
        << "   <Points>\n\n"
        << "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(int k = 0; k < km; k++){
        const double sinphi_k = sin(phi.z[k]), cosphi_k = cos(phi.z[k]);
        for(int j = j0; j < j1; j++){
            const double sinthe_j = sin(the.z[j]), costhe_j = cos(the.z[j]);
            for(int i = 0; i < im; i++){
                const double r = rad.z[i] + 11.0;                  // display radius, as before
                buf << r * sinthe_j * cosphi_k << " " << r * sinthe_j * sinphi_k << " " << r * costhe_j << '\n';
            }
        }
    }
    buf << "    </DataArray>\n\n"
        << "   </Points>\n\n"
        << "  </Piece>\n\n"
        << " </StructuredGrid>\n\n"
        << "</VTKFile>\n\n";
    f << buf.str();
    f.close();
    cout << "   File:  " << fname << "  has been written\n";
}
/*
 * Hadean_Atm_cloudview_<n>.vtk (2026-11-27, ported from the giants' <Planet>_cloudview; written with the
 * panorama, ATM_CLOUDVIEW=0 switches it off): what a camera above the planet would see of the clouds. Every
 * radial layer i is a sheet holding the column mass ABOVE it, summed from the top down to and including i:
 *
 *     X_above(i, j, k) = sum_{i' = i .. im-2} q_X(i', j, k) * rho(i', j, k) * (z(i'+1) - z(i'))
 *
 * with q_X the mass fraction (kg/kg), rho = r_humid (kg/m3) and z = get_layer_height (m), so layer 0 is the
 * whole column (the photo) and stepping i upward strips the deeper cloud away. Cells inside the topography
 * (i < i_topography) count 0; the lid im-1 has no thickness and holds 0.
 *   WaterVapour_above_kgm2    precipitable water above the layer, kg/m2
 *   CloudWater_above_gm2 / CloudIce_above_gm2 / CloudGraupel_above_gm2   g/m2
 *   Condensate_above_gm2      the sum of the three
 * Grid: concentric spherical shells, r = 1 + ATM_CLOUDVIEW_SHELL (im-1-i)/(im-1) (default 0.05, a display
 * exaggeration), layer 0 OUTERMOST, so the globe seen from outside is the photo; Extract Subset / Clip on the
 * z index exposes the deeper layers. The pole rows are left out as on the sphere .vts.
 * PhotoRGB (colour by PhotoRGB, untick "Map Scalars"): a recipe, not radiative transfer -- a dark sea
 * (0.06, 0.12, 0.22) mixed toward white by the condensate column, b = its value scaled between the 2nd and 98th
 * percentile of LAYER 0 over |lat| <= 80 deg, so the scale is the same on every layer of one file (not between
 * files); the ice + graupel share tints it slightly blue-white.
 */
void cAtmosphereModel::paraview_cloudview(string &Name_Bathymetry_File, int n){
    const string fname = output_path + "/" + Name_Bathymetry_File + "_Atm_cloudview_" + std::to_string(n) + ".vtk";
    ofstream f;
    f.precision(4);
    f.setf(ios::fixed);
    f.open(fname);
    if(!f.is_open()){
        cerr << "ERROR: could not open cloudview file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    static const double cv_shell = [](){ const char* e = getenv("ATM_CLOUDVIEW_SHELL"); return e ? atof(e) : 0.05; }();
    const int j0 = 1, j1 = jm - 1;
    auto keep_ = [&](size_t q){ const int j_ = (int)((q / km) % jm); return j_ >= j0 && j_ < j1; };
    const size_t npt = (size_t)im * (j1 - j0) * km;
    f << "# vtk DataFile Version 3.0\n";
    f << "Cloudview_Data_Atmosphere_Circulation\n";
    f << "ASCII\n";
    f << "DATASET STRUCTURED_GRID\n";
    f << "DIMENSIONS " << km << " " << (j1 - j0) << " " << im << "\n";
    f << "POINTS " << npt << " float\n";
    {
        const std::streamsize prec_ = f.precision();
        f.precision(6);                                             // 4 decimals would merge shells 0.05/(im-1) apart
        for(int i = 0; i < im; i++)
            for(int j = j0; j < j1; j++)
                for(int k = 0; k < km; k++){
                    const double r = 1.0 + cv_shell * (double)(im - 1 - i) / (double)(im - 1);
                    f << r * sin(the.z[j]) * cos(phi.z[k]) << " " << r * sin(the.z[j]) * sin(phi.z[k])
                      << " " << r * cos(the.z[j]) << "\n";
                }
        f.precision(prec_);
    }
    f << "POINT_DATA " << npt << "\n";
    auto column = [&](Array& a, double unit){
        std::vector<double> col((size_t)im * jm * km, 0.0);         // col[i][j][k], layer im-1 stays 0
        for(int j = 0; j < jm; j++) for(int k = 0; k < km; k++){
            double s = 0.0;
            const int i_ground = i_topography[j][k];
            for(int i = im - 2; i >= 0; i--){
                if(i >= i_ground){
                    const double dz = (double)get_layer_height(i + 1) - (double)get_layer_height(i);
                    // safe_val: -ffast-math turns std::isfinite into constant true (see the top of this file)
                    const double v_ = ParaViewAtm::safe_val(a.x[i][j][k] * r_humid.x[i][j][k] * dz * unit);
                    if(v_ > 0.0) s += v_;
                }
                col[((size_t)i * jm + j) * km + k] = s;
            }
        }
        return col;
    };
    auto write_ = [&](const char* desc, const std::vector<double>& col){
        f << "SCALARS " << desc << " float 1\n";
        f << "LOOKUP_TABLE default\n";
        for(size_t q = 0; q < col.size(); q++) if(keep_(q)) f << col[q] << "\n";
    };
    const std::vector<double> c_vap = column(c, 1.0);               // kg/m2
    const std::vector<double> c_liq = column(cloud, 1.0e3);         // g/m2
    const std::vector<double> c_ice = column(ice, 1.0e3);
    const std::vector<double> c_gr  = column(gr, 1.0e3);
    std::vector<double> c_all(c_liq.size());
    for(size_t q = 0; q < c_all.size(); q++) c_all[q] = c_liq[q] + c_ice[q] + c_gr[q];
    write_("WaterVapour_above_kgm2", c_vap);
    write_("CloudWater_above_gm2", c_liq);
    write_("CloudIce_above_gm2", c_ice);
    write_("CloudGraupel_above_gm2", c_gr);
    write_("Condensate_above_gm2", c_all);
    {
        std::vector<double> v_;
        for(int j = 10; j < jm - 10; j++) for(int k = 0; k < km - 1; k++) v_.push_back(c_all[(size_t)j * km + k]);
        std::sort(v_.begin(), v_.end());
        const double lo = v_[v_.size() * 2 / 100], hi0 = v_[v_.size() * 98 / 100];
        const double hi = hi0 > lo ? hi0 : lo + 1.0;
        f << "COLOR_SCALARS PhotoRGB 3\n";
        for(size_t q = 0; q < c_all.size(); q++){
            if(!keep_(q)) continue;
            const double b = std::min(1.0, std::max(0.0, (c_all[q] - lo) / (hi - lo)));
            const double fi = c_all[q] > 0.0 ? (c_ice[q] + c_gr[q]) / c_all[q] : 0.0;
            const double sea[3] = { 0.06, 0.12, 0.22 };
            const double cl[3]  = { 1.0 - 0.08 * fi, 1.0 - 0.04 * fi, 1.0 };
            f << sea[0] + b * (cl[0] - sea[0]) << " " << sea[1] + b * (cl[1] - sea[1]) << " "
              << sea[2] + b * (cl[2] - sea[2]) << "\n";
        }
    }
    f.close();
    cout << "   File:  " << fname << "  has been written\n";
}
/*
 * 
*/
void cAtmosphereModel::paraview_vtk_radial(string &Name_Bathymetry_File, 
    int i_radial, int n){
    using namespace ParaViewAtm;
    string Atmosphere_radial_File_Name = output_path + "/"
        + Name_Bathymetry_File + "_Atm_radial_" + std::to_string(i_radial)
        + "_" + std::to_string(n) + ".vtk";
    ofstream Atmosphere_vtk_radial_File;
    Atmosphere_vtk_radial_File.precision(8);
    Atmosphere_vtk_radial_File.setf(ios::fixed);
    Atmosphere_vtk_radial_File.open(Atmosphere_radial_File_Name);
    if(!Atmosphere_vtk_radial_File.is_open()){
        cerr << "ERROR: could not open paraview_vtk file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }

    ostringstream buf;
    buf.precision(8);
    buf.setf(ios::fixed);

    buf << "# vtk DataFile Version 3.0\n"
        << "Radial_Data_Atmosphere_Circulation\n"
        << "ASCII\n"
        << "DATASET STRUCTURED_GRID\n"
        << "DIMENSIONS " << km << " " << jm << " " << 1 << '\n'
        << "POINTS " << jm * km << " float\n";

    double dx = 0.1;
    double dy = 0.1;

    for(int j = 0; j < jm; j++){
        double x = j * dx;
        for(int k = 0; k < km; k++){
            buf << x << " " << k * dy << " " << 0.0 << '\n';
        }
    }

    buf << "POINT_DATA " << jm * km << '\n';
    Atmosphere_vtk_radial_File << buf.str();
    Atmosphere_vtk_radial_File.flush();
    if(!Atmosphere_vtk_radial_File){
        cerr << "ERROR: write failed for " << Atmosphere_radial_File_Name
             << " after header (" << jm * km << " points, ~"
             << (buf.str().size() / 1024) << " KB). Stream state: "
             << "fail=" << Atmosphere_vtk_radial_File.fail()
             << " bad=" << Atmosphere_vtk_radial_File.bad() << "\n";
        abort();
    }
    buf.str("");
    buf.clear();
    dump_radial("u-Component", u, u_0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("v-Component", v, u_0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("w-Component", w, u_0, i_radial, Atmosphere_vtk_radial_File);

    buf << "SCALARS Temperature float " << 1 << '\n'
        << "LOOKUP_TABLE default\n";
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            buf << ParaViewAtm::safe_val(t.x[i_radial][j][k] * t_0 - t_0) << '\n';
        }
    }
    Atmosphere_vtk_radial_File << buf.str();
    buf.str("");
    buf.clear();
    dump_radial_2d("v-velocity_NASA", velocity_v_NASA, 1.0, Atmosphere_vtk_radial_File);
    dump_radial_2d("w-velocity_NASA", velocity_w_NASA, 1.0, Atmosphere_vtk_radial_File);

    dump_radial_2d("Temperature_NASA", temperature_NASA, 1.0, Atmosphere_vtk_radial_File);
//    dump_radial_2d("Temperature_Reconst", temp_reconst, 1.0, Atmosphere_vtk_radial_File);
    dump_radial_2d("Temperature_Landscape", temp_landscape, 1.0, Atmosphere_vtk_radial_File);

    dump_radial("Radiation", radiation, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("Epsilon", epsilon, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("TauAbove", tau_above, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("TauLayer", tau_layer, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("BruntVaisala_N2", brunt_N2, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("PsiMerid", Psi, 1.0, i_radial, Atmosphere_vtk_radial_File);
    // Ubud_* dropped from the VTK output 2026-08-20: six radial-momentum-budget fields
    // per slice was a third of the file for a diagnostic that is read as min/max in
    // Results_Atm, not as a field. The arrays are still computed and still printed
    // there; only the dumps are off. Uncomment to put them back.
//  dump_radial("Ubud_pgf", ubud_pgf, 1.0, i_radial, Atmosphere_vtk_radial_File);
//  dump_radial("Ubud_cor", ubud_cor, 1.0, i_radial, Atmosphere_vtk_radial_File);
//  dump_radial("Ubud_advv", ubud_advv, 1.0, i_radial, Atmosphere_vtk_radial_File);
//  dump_radial("Ubud_advh", ubud_advh, 1.0, i_radial, Atmosphere_vtk_radial_File);
//  dump_radial("Ubud_diff", ubud_diff, 1.0, i_radial, Atmosphere_vtk_radial_File);
//  dump_radial("Ubud_buoy", ubud_buoy, 1.0, i_radial, Atmosphere_vtk_radial_File);

    dump_radial_2d("Tropopause", Tropopause, 1.0, Atmosphere_vtk_radial_File);

    dump_radial_2d("Pressure_Landscape_bar", p_stat_landscape, 1.0e-3, Atmosphere_vtk_radial_File);

    dump_radial_2d("Evaporation", Evaporation, 1.0, Atmosphere_vtk_radial_File);

    dump_radial_2d("Albedo", albedo, 1.0, Atmosphere_vtk_radial_File);

    dump_radial("Topography", h, 1.0, i_radial, Atmosphere_vtk_radial_File);
    if(fluid_phase_on()) dump_radial("FluidPhase", fluid_phase_field(*this), 1.0, i_radial, Atmosphere_vtk_radial_File);

    dump_radial("WaterVapour", c, 1000.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("CloudWater", cloud, 1e3, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("CloudIce", ice, 1e3, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("CloudGraupel", gr, 1e3, i_radial, Atmosphere_vtk_radial_File);

    dump_radial("S_c_c", S_c_c, 1e3, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("S_v", S_v, 1e3, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("S_c", S_c, 1e3, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("S_i", S_i, 1e3, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("S_r", S_r, 1e3, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("S_s", S_s, 1e3, i_radial, Atmosphere_vtk_radial_File);

    dump_radial("Precipitation", Precipitation, 8.64e4, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("PrecipitationRain", P_rain, 8.64e4, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("PrecipitationSnow", P_snow, 8.64e4, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("PrecipitationGraupel", P_graupel, 8.64e4, i_radial, Atmosphere_vtk_radial_File);

    dump_radial_2d("PrecipitableWater", precipitable_water, 1.0, Atmosphere_vtk_radial_File);
    dump_radial_2d("Precipitation_NASA", precipitation_NASA, 1.0, Atmosphere_vtk_radial_File);

    dump_radial("Pressure_bar", p_stat, 1.0e-3, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("PressureDynamic", p_dyn, p_0, i_radial, Atmosphere_vtk_radial_File);

//    dump_radial("r_dry", r_dry, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("r_humid", r_humid, 1.0, i_radial, Atmosphere_vtk_radial_File);
//    dump_radial("HumidityRel", HumidityRel, 1.0, i_radial, Atmosphere_vtk_radial_File);
//    dump_radial("TempDewPoint", TempDewPoint, 1.0, i_radial, Atmosphere_vtk_radial_File);
//    dump_radial("TempStand", TempStand, 1.0, i_radial, Atmosphere_vtk_radial_File);

/*
    dump_radial_2d("Q_radiation_2D", Q_radiation, 1e-3, Atmosphere_vtk_radial_File);
    dump_radial_2d("Q_bottom_2D", Q_bottom, 1e-3, Atmosphere_vtk_radial_File);
    dump_radial_2d("Q_latent_2D", Q_latent, 1e-3, Atmosphere_vtk_radial_File);
    dump_radial_2d("Q_sensible_2D", Q_sensible, 1e-3, Atmosphere_vtk_radial_File);
*/

    dump_radial("BuoyancyForce", BuoyancyForce, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("CoriolisForce", CoriolisForce, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("CentrifugalForce", CentrifugalForce, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("PresGradForce", PresGradForce, 1.0, i_radial, Atmosphere_vtk_radial_File);
//    dump_radial("LorentzForce", LorentzForce, 1e9, i_radial, Atmosphere_vtk_radial_File);
//    dump_radial("Q_Radiation", radiation, 1.0, i_radial, Atmosphere_vtk_radial_File);

    dump_radial("Q_Latent", Q_Latent, 1e-3, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("Q_Sensible", Q_Sensible, 1e-3, i_radial, Atmosphere_vtk_radial_File);

    // ==================== THE EIGHT SPECIES, ALL FROM split() ====================
    // One source for all eight: AtmMixture::split, the routine the thermodynamics and the
    // radiation use, so a plotted mass fraction cannot drift from the integrated one.
    //   H2O  = q_v, the vapour;
    //   CO2  = q_c, the CARRIER-RENORMALISED CO2 — NOT the raw co2 array, which is the
    //          transported tracer and carries no dilution. That one is written separately
    //          below as CO2_tracer.
    //   rest = q_b*f_bg[n], the per-species background (ATHAD item 60).
    // All are MASS FRACTIONS, dimensionless, so the eight are directly comparable. The older
    // WaterVapour field above is the same water in g/kg, read from the raw c array.
    for (int n = 0; n < 8; n++) {
        Atmosphere_vtk_radial_File << "SCALARS " << AtmMixture::SPECIES_NAMES()[n] << " float " << 1 << "\n";
        Atmosphere_vtk_radial_File << "LOOKUP_TABLE default" << "\n";
        const double f_n = (n >= 2) ? m_comp.f_bg[n-2] : 0.0;
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                double q_v, q_c, q_b;
                const double q_l = cloud.x[i_radial][j][k] + ice.x[i_radial][j][k] + gr.x[i_radial][j][k];
                AtmMixture::split(c.x[i_radial][j][k], co2.x[i_radial][j][k], q_v, q_c, q_b, q_l);
                const double q_n = (n == 0) ? q_v : (n == 1) ? q_c : q_b * f_n;
                Atmosphere_vtk_radial_File << safe_val(q_n) << "\n";
            }
        }
    }
    // The transported tracer itself, undiluted: the array rhs_co2 advances.
    dump_radial("CO2_tracer", co2, 1.0, i_radial, Atmosphere_vtk_radial_File);

    if(turb_model != "laminar"){
        dump_radial("TKE",        tke,        1.0, i_radial, Atmosphere_vtk_radial_File);
        dump_radial("Dissipation", dis,        1.0, i_radial, Atmosphere_vtk_radial_File);
        dump_radial("EddyViscosity", nue,      1.0, i_radial, Atmosphere_vtk_radial_File);
        dump_radial("TurbProd",   prod,        1.0, i_radial, Atmosphere_vtk_radial_File);
        dump_radial("TKE_Source", tke_source,  1.0, i_radial, Atmosphere_vtk_radial_File);
        dump_radial("Dis_Source", dis_source,  1.0, i_radial, Atmosphere_vtk_radial_File);
    }

//    dump_radial("TempStandard", TempStand, 1.0, i_radial, Atmosphere_vtk_radial_File);
//    dump_radial("TempDewPoint", TempDewPoint, 1.0, i_radial, Atmosphere_vtk_radial_File);

//    dump_radial_2d("Evaporation_Dalton", Evaporation_Dalton, 1.0, Atmosphere_vtk_radial_File);
//    dump_radial_2d("Evaporation_Meyer", Evaporation_Meyer, 1.0, Atmosphere_vtk_radial_File);
//    dump_radial_2d("Evaporation_Rohwer", Evaporation_Rohwer, 1.0, Atmosphere_vtk_radial_File);
    dump_radial_2d("Evaporation", Evaporation, 1.0, Atmosphere_vtk_radial_File);

    dump_radial_2d("Vegetation", Vegetation, 1.0, Atmosphere_vtk_radial_File);
    dump_radial_2d("Landscape", Landscape, 1.0, Atmosphere_vtk_radial_File);
    dump_radial_2d("Tropopause", Tropopause, 1.0, Atmosphere_vtk_radial_File);

    buf << "VECTORS v-w-ATOM float\n";
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            buf << safe_val(v.x[i_radial][j][k])
                << " " << safe_val(w.x[i_radial][j][k]) << " " << 0.0 << '\n';
        }
    }

    buf << "VECTORS v-w-NASA float\n";
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            buf << velocity_v_NASA.y[j][k]
                << " " << velocity_w_NASA.y[j][k] << " " << 0.0 << '\n';
        }
    }
    Atmosphere_vtk_radial_File << buf.str();
    buf.str("");
    buf.clear();


    dump_radial("PrecipitationConv", P_conv, 8.64e4, i_radial, Atmosphere_vtk_radial_File);
//    dump_radial("MC_t", MC_t, 1.0, i_radial, Atmosphere_vtk_radial_File);
//    dump_radial("MC_q", MC_q, 1e3, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("MC_v", MC_v, 1e3, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("MC_w", MC_w, 1e3, i_radial, Atmosphere_vtk_radial_File);
//    dump_radial("M_u", M_u, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("M_d", M_d, 1.0, i_radial, Atmosphere_vtk_radial_File);
//    dump_radial_2d("CloudBase", i_Base, 1.0, Atmosphere_vtk_radial_File);
//    dump_radial_2d("CloudLFS", i_LFS, 1.0, Atmosphere_vtk_radial_File);
//    dump_radial("c_u", c_u, 1e3, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("e_d", e_d, 1e3, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("e_p", e_p, 1e3, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("q_v_u", q_v_u, 1e3, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("q_c_u", q_c_u, 1e3, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("q_v_d", q_v_d, 1e3, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("u_u", u_u, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("v_u", v_u, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("w_u", w_u, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("u_d", u_d, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("v_d", v_d, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("w_d", w_d, 1.0, i_radial, Atmosphere_vtk_radial_File);


    buf << "VECTORS v-w-Updraft float\n";
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            if(u.x[i_radial][j][k] < 0.0) u_u.x[i_radial][j][k] = 0.0;
                buf << safe_val(v_u.x[i_radial][j][k])
                    << " " << safe_val(w_u.x[i_radial][j][k]) << " " << 0.0 << '\n';
        }
    }
    buf << "VECTORS v-w-Downdraft float\n";
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
        if(u.x[i_radial][j][k] < 0.0) u_u.x[i_radial][j][k] = 0.0;
        buf << safe_val(v_d.x[i_radial][j][k])
            << " " << safe_val(w_d.x[i_radial][j][k]) << " " << 0.0 << '\n';
        }
    }
    Atmosphere_vtk_radial_File << buf.str();


    Atmosphere_vtk_radial_File.close();
    cout << "   File:  " << Atmosphere_radial_File_Name
        << "  has been written\n";
}
/*
 * 
*/
void cAtmosphereModel::paraview_vtk_zonal(string &Name_Bathymetry_File, 
    int k_zonal, int n){
    using namespace ParaViewAtm;
    string Atmosphere_zonal_File_Name = output_path + "/" + Name_Bathymetry_File
        + "_Atm_zonal_" + std::to_string(k_zonal) + "_" + std::to_string(n) + ".vtk";
    ofstream Atmosphere_vtk_zonal_File;
    Atmosphere_vtk_zonal_File.precision(8);
    Atmosphere_vtk_zonal_File.setf(ios::fixed);
    Atmosphere_vtk_zonal_File.open(Atmosphere_zonal_File_Name);
    if(!Atmosphere_vtk_zonal_File.is_open()){
        cerr << "ERROR: could not open vtk_zonal file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }

    // Buffer all output to reduce I/O syscalls
    ostringstream buf;
    buf.precision(8);
    buf.setf(ios::fixed);

    buf << "# vtk DataFile Version 3.0\n"
        << "Zonal_Data_Atmosphere_Circulation\n"
        << "ASCII\n"
        << "DATASET STRUCTURED_GRID\n"
        << "DIMENSIONS " << jm << " " << im << " " << 1 << '\n'
        << "POINTS " << im * jm << " float\n";

    double dx = 0.1;
    double dy = 0.05;

    // ==================================================================
    // THE VERTICAL AXIS IS TRUE HEIGHT, NOT LEVEL INDEX.
    //
    // x = i*dx puts the points at LEVEL INDICES, which on this exponentially stretched
    // grid is not a height axis at all: layer 0 is 1.2 km and layer 39 is 22.8 km
    // (zeta = 3, im = 41), so index space stretches the bottom of the atmosphere and
    // squashes the top by 18.6x. Everything drawn on it inherits that — contours, glyph
    // angles and streamline curvature alike — and the distortion VARIES with height, so
    // no single aspect-ratio setting in ParaView can undo it.
    //
    // Mapping the true height onto the same 0..(im-1)*dx span keeps the figure the size
    // it always was while making the axis linear in metres, so the vertical exaggeration
    // becomes ONE constant (~30x here) instead of a function of altitude.
    // ==================================================================
    const double h_top  = get_layer_height(im - 1);
    const double x_span = (im - 1) * dx;
    const double x_of_h = (h_top > 0.0) ? x_span / h_top : 0.0;   // plot units per metre

    for(int i = 0; i < im; i++){
        double x = get_layer_height(i) * x_of_h;
        for(int j = 0; j < jm; j++){
            buf << x << " " << j * dy << " " << 0.0 << '\n';
        }
    }

    buf << "POINT_DATA " << im * jm << '\n';
    Atmosphere_vtk_zonal_File << buf.str();
    buf.str("");
    buf.clear();

    dump_zonal("u-Component", u, u_0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("v-Component", v, u_0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("w-Component", w, u_0, k_zonal, Atmosphere_vtk_zonal_File);
    if(fluid_phase_on()) dump_zonal("FluidPhase", fluid_phase_field(*this), 1.0, k_zonal, Atmosphere_vtk_zonal_File);

    double inv_u_0 = 1.0 / u_0;

    // ==================================================================
    // TWO DEFECTS IN THE GLYPH/STREAMLINE VECTOR, AND ONLY ONE IS A SCALING TYPO.
    //
    // (1) The arrays hold NON-DIMENSIONAL velocity (u/u_0), so m/s is `* u_0` — which is
    //     what dump_zonal writes for the u/v/w SCALARS three lines above. The VECTORS
    //     field divided instead, giving u/u_0^2. Measured in a shipped file: the scalar
    //     v-Component and the vector's y differ by exactly 64 = u_0^2 at u_0 = 8 m/s.
    //     Direction is untouched by a uniform factor, so this never moved an arrow — it
    //     just made every magnitude 64x too small and inconsistent with the scalars
    //     plotted beside it.
    //
    // (2) THE ONE THAT MAKES GLYPHS DISAGREE WITH STREAMLINES. The geometry written
    //     above is INDEX SPACE: x = i*dx per LEVEL, y = j*dy per LATITUDE INDEX. A
    //     velocity in m/s does not live in that space. One level is 1.2 km at the surface
    //     and 22.8 km at the top (zeta = 3, im = 41) while one latitude index is a fixed
    //     ~111 km, so the two axes are compressed by wildly different and, in the radial
    //     case, height-DEPENDENT factors.
    //
    //     The physical flow here is ~5800:1 horizontal to vertical (measured: v = 3.25 m/s
    //     against u = 0.00056 m/s). Drawn as raw m/s in index space, every glyph lies flat
    //     along the latitude axis and the overturning is invisible — while Psi and the
    //     stream tracer, which integrate, resolve cells that close through exactly those
    //     tiny vertical velocities over a 300 km depth. That is the glyphs and the
    //     streamlines telling different stories about the same field.
    //
    // uv_plot converts the velocity INTO the plot's coordinates — plot-units per second —
    // so a closed cell is drawn closed and glyphs and streamlines agree by construction:
    //
    //     u_plot = u_phys * dx / dz_local        dz_local = the layer's true thickness
    //     v_plot = v_phys * dy / dl              dl = metres per latitude index
    //
    // u-v-Cell is KEPT (now correctly dimensional) because it is the honest m/s field and
    // some workflows want it; uv_plot is the one to glyph and stream-trace. Diagnostic
    // output only — no physics reads either.
    // ==================================================================
    const double dl = r_Earth * 1000.0 * M_PI / (double)(jm - 1);   // m per latitude index

    buf << "VECTORS u-v-Cell float\n";
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            buf << safe_val(u.x[i][j][k_zonal] * u_0) << " "
                << safe_val(v.x[i][j][k_zonal] * u_0) << " " << 0.0 << '\n';
        }
    }

    // Units are plot-units per DAY, not per second. The file is written ios::fixed at 8
    // decimals, and plot-units per second are ~1e-10 here, so a per-second field prints as
    // exactly 0.00000000 for most of the domain. 86400 puts it in a printable range and is
    // a real unit rather than a fudge: "how far a parcel moves across this plot in a day".
    // A uniform factor, so it cannot affect direction.
    const double per_day = 86400.0;

    // sx is now the SAME constant at every level, because the geometry above is linear in
    // height. On the old index-space geometry it had to be dx/dz_local, which varied 18.6x
    // across the column — correct for that geometry, but the geometry was the defect.
    const double sx = x_of_h;      // plot units per metre, vertical
    const double sy = dy / dl;     // plot units per metre, meridional

    buf << "VECTORS uv_plot float\n";
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            buf << safe_val(u.x[i][j][k_zonal] * u_0 * sx * per_day) << " "
                << safe_val(v.x[i][j][k_zonal] * u_0 * sy * per_day) << " " << 0.0 << '\n';
        }
    }

    buf << "SCALARS Temperature float " << 1 << '\n'
        << "LOOKUP_TABLE default\n";
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            buf << safe_val(t.x[i][j][k_zonal] * t_0 - t_0) << '\n';
        }
    }
    Atmosphere_vtk_zonal_File << buf.str();
    buf.str("");
    buf.clear();

    // height is invariant across j, hoist get_layer_height out of j-loop
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < im; i++){
        double height = get_layer_height(i);
        for(int j = 0; j < jm; j++){
            aux_t.x[i][j][k_zonal] = height;
        }
    }


    dump_zonal("Radiation", radiation, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Epsilon", epsilon, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("TauAbove", tau_above, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("TauLayer", tau_layer, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("BruntVaisala_N2", brunt_N2, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("PsiMerid", Psi, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    // Ubud_* dropped from the VTK output 2026-08-20: six radial-momentum-budget fields
    // per slice was a third of the file for a diagnostic that is read as min/max in
    // Results_Atm, not as a field. The arrays are still computed and still printed
    // there; only the dumps are off. Uncomment to put them back.
//  dump_zonal("Ubud_pgf", ubud_pgf, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
//  dump_zonal("Ubud_cor", ubud_cor, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
//  dump_zonal("Ubud_advv", ubud_advv, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
//  dump_zonal("Ubud_advh", ubud_advh, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
//  dump_zonal("Ubud_diff", ubud_diff, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
//  dump_zonal("Ubud_buoy", ubud_buoy, 1.0, k_zonal, Atmosphere_vtk_zonal_File);

    dump_zonal("Topography", h, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("height", aux_t, 1e-3, k_zonal, Atmosphere_vtk_zonal_File);

//    dump_zonal("TempStandard", TempStand, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
//    dump_zonal("TempDewPoint", TempDewPoint, 1.0, k_zonal, Atmosphere_vtk_zonal_File);

    dump_zonal("PressureDynamic", p_dyn, p_0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Pressure_bar", p_stat, 1.0e-3, k_zonal, Atmosphere_vtk_zonal_File);

//    dump_zonal("r_dry", r_dry, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("r_humid", r_humid, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("HumidityRel", HumidityRel, 1.0, k_zonal, Atmosphere_vtk_zonal_File);

    dump_zonal("WaterVapour", c, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("CloudWater", cloud, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("CloudIce", ice, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("CloudGraupel", gr, 1e3, k_zonal, Atmosphere_vtk_zonal_File);

    dump_zonal("Precipitation", Precipitation, 8.64e4, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("PrecipitationRain", P_rain, 8.64e4, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("PrecipitationSnow", P_snow, 8.64e4, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("PrecipitationGraupel", P_graupel, 8.64e4, k_zonal, Atmosphere_vtk_zonal_File);

    dump_zonal("S_v", S_v, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("S_c", S_c, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("S_i", S_i, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("S_g", S_g, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("S_r", S_r, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("S_s", S_s, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("S_c_c", S_c_c, 1e3, k_zonal, Atmosphere_vtk_zonal_File);


    dump_zonal("BuoyancyForce", BuoyancyForce, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("CoriolisForce", CoriolisForce, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("CentrifugalForce", CentrifugalForce, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("PresGradForce", PresGradForce, 1.0, k_zonal, Atmosphere_vtk_zonal_File);

    dump_zonal("Q_Latent", Q_Latent, 1e-3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Q_Sensible", Q_Sensible, 1e-3, k_zonal, Atmosphere_vtk_zonal_File);

    // ==================== THE EIGHT SPECIES, ALL FROM split() ====================
    // One source for all eight: AtmMixture::split, the routine the thermodynamics and the
    // radiation use, so a plotted mass fraction cannot drift from the integrated one.
    //   H2O  = q_v, the vapour;
    //   CO2  = q_c, the CARRIER-RENORMALISED CO2 — NOT the raw co2 array, which is the
    //          transported tracer and carries no dilution. That one is written separately
    //          below as CO2_tracer.
    //   rest = q_b*f_bg[n], the per-species background (ATHAD item 60).
    // All are MASS FRACTIONS, dimensionless, so the eight are directly comparable. The older
    // WaterVapour field above is the same water in g/kg, read from the raw c array.
    for (int n = 0; n < 8; n++) {
        Atmosphere_vtk_zonal_File << "SCALARS " << AtmMixture::SPECIES_NAMES()[n] << " float " << 1 << "\n";
        Atmosphere_vtk_zonal_File << "LOOKUP_TABLE default" << "\n";
        const double f_n = (n >= 2) ? m_comp.f_bg[n-2] : 0.0;
        for(int i = 0; i < im; i++){
            for(int j = 0; j < jm; j++){
                double q_v, q_c, q_b;
                const double q_l = cloud.x[i][j][k_zonal] + ice.x[i][j][k_zonal] + gr.x[i][j][k_zonal];
                AtmMixture::split(c.x[i][j][k_zonal], co2.x[i][j][k_zonal], q_v, q_c, q_b, q_l);
                const double q_n = (n == 0) ? q_v : (n == 1) ? q_c : q_b * f_n;
                Atmosphere_vtk_zonal_File << safe_val(q_n) << "\n";
            }
        }
    }
    // The transported tracer itself, undiluted: the array rhs_co2 advances.
    dump_zonal("CO2_tracer", co2, 1.0, k_zonal, Atmosphere_vtk_zonal_File);

    if(turb_model != "laminar"){
        dump_zonal("TKE",          tke,        1.0, k_zonal, Atmosphere_vtk_zonal_File);
        dump_zonal("Dissipation",  dis,        1.0, k_zonal, Atmosphere_vtk_zonal_File);
        dump_zonal("EddyViscosity", nue,       1.0, k_zonal, Atmosphere_vtk_zonal_File);
        dump_zonal("TurbProd",     prod,       1.0, k_zonal, Atmosphere_vtk_zonal_File);
        dump_zonal("TKE_Source",   tke_source, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
        dump_zonal("Dis_Source",   dis_source, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    }

    dump_zonal("Cloud_Base", CloudBase, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("LFS", LevelFreeSinking, 1.0, k_zonal, Atmosphere_vtk_zonal_File);

    dump_zonal("Deep_beg", Deep_beg, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Deep_end", Deep_end, 1.0, k_zonal, Atmosphere_vtk_zonal_File);

    dump_zonal("PrecipitationConv", P_conv, 8.64e4, k_zonal, Atmosphere_vtk_zonal_File);

    dump_zonal("E_u", E_u, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("D_u", D_u, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("E_d", E_d, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("D_d", D_d, 1e3, k_zonal, Atmosphere_vtk_zonal_File);

    dump_zonal("M_u", M_u, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("M_d", M_d, 1e3, k_zonal, Atmosphere_vtk_zonal_File);

    dump_zonal("MC_t", MC_t, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("MC_q", MC_q, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("MC_v", MC_v, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("MC_w", MC_w, 1e3, k_zonal, Atmosphere_vtk_zonal_File);

    dump_zonal("g_p", g_p, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("e_d", e_d, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("e_p", e_p, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("e_l", e_l, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("c_u", c_u, 1e3, k_zonal, Atmosphere_vtk_zonal_File);

    dump_zonal("q_c_u", q_c_u, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("q_v_u", q_v_u, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("q_v_d", q_v_d, 1e3, k_zonal, Atmosphere_vtk_zonal_File);

    dump_zonal("s", s, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("s_u", s_u, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("s_d", s_d, 1.0, k_zonal, Atmosphere_vtk_zonal_File);

    dump_zonal("u_u", u_u, u_0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("u_d", u_d, u_0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("v_u", v_u, u_0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("v_d", v_d, u_0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("w_u", w_u, u_0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("w_d", w_d, u_0, k_zonal, Atmosphere_vtk_zonal_File);


    buf << "VECTORS u-v-Updraft float\n";
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            if(u.x[i][j][k_zonal] < 0.0) u_u.x[i][j][k_zonal] = 0.0;
            buf << safe_val(u_u.x[i][j][k_zonal] * inv_u_0)
                << " " << safe_val(v_u.x[i][j][k_zonal] * inv_u_0) << " " << 0.0 << '\n';
        }
    }

    buf << "VECTORS u-v-Downdraft float\n";
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            if(u.x[i][j][k_zonal] < 0.0) u_d.x[i][j][k_zonal] = 0.0;
            buf << safe_val(u_d.x[i][j][k_zonal] * inv_u_0)
                << " " << safe_val(v_d.x[i][j][k_zonal] * inv_u_0) << " " << 0.0 << '\n';
        }
    }
    Atmosphere_vtk_zonal_File << buf.str();


    Atmosphere_vtk_zonal_File.close();
    cout << "   File:  " << Atmosphere_zonal_File_Name
        << "  has been written\n";
}
/*
 * 
*/
void cAtmosphereModel::paraview_vtk_longal(string &Name_Bathymetry_File, 
    int j_longal, int n){
    using namespace ParaViewAtm;
    string Atmosphere_longal_File_Name = output_path + "/" + Name_Bathymetry_File
        + "_Atm_longal_" + std::to_string(j_longal) + "_" + std::to_string(n) + ".vtk";
    ofstream Atmosphere_vtk_longal_File;
    Atmosphere_vtk_longal_File.precision(8);
    Atmosphere_vtk_longal_File.setf(ios::fixed);
    Atmosphere_vtk_longal_File.open(Atmosphere_longal_File_Name);
    if(!Atmosphere_vtk_longal_File.is_open()){
        cerr << "ERROR: could not open vtk_longal file " << __FILE__
            << " at line " << __LINE__ << "\n";
        abort();
    }

    ostringstream buf;
    buf.precision(8);
    buf.setf(ios::fixed);

    buf << "# vtk DataFile Version 3.0\n"
        << "Longitudinal_Data_Atmosphere_Circulation\n"
        << "ASCII\n"
        << "DATASET STRUCTURED_GRID\n"
        << "DIMENSIONS " << km << " " << im << " " << 1 << '\n'
        << "POINTS " << im * km << " float\n";

    double dx = 0.1;
    double dz = 0.025;

    // THE VERTICAL AXIS IS TRUE HEIGHT, NOT LEVEL INDEX -- the same repair the zonal writer
    // got in item 70, which fixed that writer only and left this one on index space. On an
    // exponentially stretched grid, x = i*dx stretches the bottom of the atmosphere and
    // squashes the top by the ratio of the thickest layer to the thinnest, and the
    // distortion VARIES with height, so no ParaView aspect setting can undo it. Mapping the
    // true height onto the same 0..(im-1)*dx span keeps the figure the size it always was
    // while making the axis linear in metres.
    const double h_top  = get_layer_height(im - 1);
    const double x_span = (im - 1) * dx;
    const double x_of_h = (h_top > 0.0) ? x_span / h_top : 0.0;   // plot units per metre

    for(int i = 0; i < im; i++){
        double x = get_layer_height(i) * x_of_h;
        for(int k = 0; k < km; k++){
            buf << x << " " << 0.0 << " " << k * dz << '\n';
        }
    }

    buf << "POINT_DATA " << im * km << '\n';
    Atmosphere_vtk_longal_File << buf.str();
    buf.str("");
    buf.clear();
    dump_longal("u-Component", u, u_0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("v-Component", v, u_0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("w-Component", w, u_0, j_longal, Atmosphere_vtk_longal_File);
    if(fluid_phase_on()) dump_longal("FluidPhase", fluid_phase_field(*this), 1.0, j_longal, Atmosphere_vtk_longal_File);

    buf << "SCALARS Temperature float " << 1 << '\n'
        << "LOOKUP_TABLE default\n";
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            buf << safe_val(t.x[i][j_longal][k] * t_0 - t_0) << '\n';
        }
    }
    Atmosphere_vtk_longal_File << buf.str();
    buf.str("");
    buf.clear();

    for(int i = 0; i < im; i++){
        double height = get_layer_height(i);
        for(int k = 0; k < km; k++){
            aux_t.x[i][j_longal][k] = height;
        }
    }

    dump_longal("Radiation", radiation, 1.0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("Epsilon", epsilon, 1.0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("TauAbove", tau_above, 1.0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("TauLayer", tau_layer, 1.0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("BruntVaisala_N2", brunt_N2, 1.0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("PsiMerid", Psi, 1.0, j_longal, Atmosphere_vtk_longal_File);
    // Ubud_* dropped from the VTK output 2026-08-20: six radial-momentum-budget fields
    // per slice was a third of the file for a diagnostic that is read as min/max in
    // Results_Atm, not as a field. The arrays are still computed and still printed
    // there; only the dumps are off. Uncomment to put them back.
//  dump_longal("Ubud_pgf", ubud_pgf, 1.0, j_longal, Atmosphere_vtk_longal_File);
//  dump_longal("Ubud_cor", ubud_cor, 1.0, j_longal, Atmosphere_vtk_longal_File);
//  dump_longal("Ubud_advv", ubud_advv, 1.0, j_longal, Atmosphere_vtk_longal_File);
//  dump_longal("Ubud_advh", ubud_advh, 1.0, j_longal, Atmosphere_vtk_longal_File);
//  dump_longal("Ubud_diff", ubud_diff, 1.0, j_longal, Atmosphere_vtk_longal_File);
//  dump_longal("Ubud_buoy", ubud_buoy, 1.0, j_longal, Atmosphere_vtk_longal_File);

    dump_longal("Topography", h, 1.0, j_longal, Atmosphere_vtk_longal_File);

    dump_longal("WaterVapour", c, 1.0e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("CloudWater", cloud, 1e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("CloudIce", ice, 1e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("CloudGraupel", gr, 1e3, j_longal, Atmosphere_vtk_longal_File);

    dump_longal("S_c_c", S_c_c, 1e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("S_v", S_v, 1e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("S_c", S_c, 1e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("S_i", S_i, 1e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("S_r", S_r, 1e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("S_s", S_s, 1e3, j_longal, Atmosphere_vtk_longal_File);

    dump_longal("Precipitation", Precipitation, 8.64e4, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("PrecipitationRain", P_rain, 8.64e4, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("PrecipitationSnow", P_snow, 8.64e4, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("PrecipitationGraupel", P_graupel, 8.64e4, j_longal, Atmosphere_vtk_longal_File);

    dump_longal("Pressure_bar", p_stat, 1.0e-3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("PressureDynamic", p_dyn, p_0, j_longal, Atmosphere_vtk_longal_File);

//    dump_longal("r_dry", r_dry, 1.0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("r_humid", r_humid, 1.0, j_longal, Atmosphere_vtk_longal_File);
//    dump_longal("HumidityRel", HumidityRel, 1.0, j_longal, Atmosphere_vtk_longal_File);
//    dump_longal("TempDewPoint", TempDewPoint, 1.0, j_longal, Atmosphere_vtk_longal_File);
//    dump_longal("TempStand", TempStand, 1.0, j_longal, Atmosphere_vtk_longal_File);

    dump_longal("BuoyancyForce", BuoyancyForce, 1.0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("CoriolisForce", CoriolisForce, 1.0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("CentrifugalForce", CentrifugalForce, 1.0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("PresGradForce", PresGradForce, 1.0, j_longal, Atmosphere_vtk_longal_File);

    dump_longal("Q_Latent", Q_Latent, 1e-3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("Q_Sensible", Q_Sensible, 1e-3, j_longal, Atmosphere_vtk_longal_File);
//    dump_longal("Q_Radiation", radiation, 1.0, j_longal, Atmosphere_vtk_longal_File);
//    dump_longal("TempStandard", TempStand, 1.0, j_longal, Atmosphere_vtk_longal_File);
//    dump_longal("TempDewPoint", TempDewPoint, 1.0, j_longal, Atmosphere_vtk_longal_File);

    // ==================== THE EIGHT SPECIES, ALL FROM split() ====================
    // One source for all eight: AtmMixture::split, the routine the thermodynamics and the
    // radiation use, so a plotted mass fraction cannot drift from the integrated one.
    //   H2O  = q_v, the vapour;
    //   CO2  = q_c, the CARRIER-RENORMALISED CO2 — NOT the raw co2 array, which is the
    //          transported tracer and carries no dilution. That one is written separately
    //          below as CO2_tracer.
    //   rest = q_b*f_bg[n], the per-species background (ATHAD item 60).
    // All are MASS FRACTIONS, dimensionless, so the eight are directly comparable. The older
    // WaterVapour field above is the same water in g/kg, read from the raw c array.
    for (int n = 0; n < 8; n++) {
        Atmosphere_vtk_longal_File << "SCALARS " << AtmMixture::SPECIES_NAMES()[n] << " float " << 1 << "\n";
        Atmosphere_vtk_longal_File << "LOOKUP_TABLE default" << "\n";
        const double f_n = (n >= 2) ? m_comp.f_bg[n-2] : 0.0;
        for(int i = 0; i < im; i++){
            for(int k = 0; k < km; k++){
                double q_v, q_c, q_b;
                const double q_l = cloud.x[i][j_longal][k] + ice.x[i][j_longal][k] + gr.x[i][j_longal][k];
                AtmMixture::split(c.x[i][j_longal][k], co2.x[i][j_longal][k], q_v, q_c, q_b, q_l);
                const double q_n = (n == 0) ? q_v : (n == 1) ? q_c : q_b * f_n;
                Atmosphere_vtk_longal_File << safe_val(q_n) << "\n";
            }
        }
    }
    // The transported tracer itself, undiluted: the array rhs_co2 advances.
    dump_longal("CO2_tracer", co2, 1.0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("height", aux_t, 1e-3, j_longal, Atmosphere_vtk_longal_File);

    if(turb_model != "laminar"){
        dump_longal("TKE",           tke,        1.0, j_longal, Atmosphere_vtk_longal_File);
        dump_longal("Dissipation",   dis,        1.0, j_longal, Atmosphere_vtk_longal_File);
        dump_longal("EddyViscosity", nue,        1.0, j_longal, Atmosphere_vtk_longal_File);
        dump_longal("TurbProd",      prod,       1.0, j_longal, Atmosphere_vtk_longal_File);
        dump_longal("TKE_Source",    tke_source, 1.0, j_longal, Atmosphere_vtk_longal_File);
        dump_longal("Dis_Source",    dis_source, 1.0, j_longal, Atmosphere_vtk_longal_File);
    }


    buf << "VECTORS u-w-Cell float\n";
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            // * u_0: the scalars written beside this vector by dump_longal are
            // dimensional, and this vector was not -- the same component/vector unit
            // mismatch item 70 found in the zonal writer, in the other direction.
            buf << safe_val(u.x[i][j_longal][k] * u_0)
                << " " << 0.0 << " " << safe_val(w.x[i][j_longal][k] * u_0) << '\n';
        }
    }
    Atmosphere_vtk_longal_File << buf.str();
    buf.str("");
    buf.clear();

    dump_longal("Cloud_Base", CloudBase, 1.0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("LFS", LevelFreeSinking, 1.0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("PrecipitationConv", P_conv, 8.64e4, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("M_u", M_u, 1.0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("M_d", M_d, 1.0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("MC_t", MC_t, 1.0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("MC_q", MC_q, 1e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("MC_v", MC_v, 1e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("MC_w", MC_w, 1e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("g_p", g_p, 1.0e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("e_d", e_d, 1.0e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("e_p", e_p, 1.0e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("e_l", e_l, 1.0e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("c_u", c_u, 1.0e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("q_c_u", q_c_u, 1.0e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("q_v_u", q_v_u, 1.0e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("q_v_d", q_v_d, 1.0e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("s", s, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("s_u", s_u, 1.0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("s_d", s_d, 1.0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("u_u", u_u, u_0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("u_d", u_d, u_0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("v_u", v_u, u_0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("v_d", v_d, u_0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("w_u", w_u, u_0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("w_d", w_d, u_0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("E_u", E_u, 1.0e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("D_u", D_u, 1.0e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("E_d", E_d, 1.0e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("D_d", D_d, 1.0e3, j_longal, Atmosphere_vtk_longal_File);


    buf << "VECTORS u-w-Updraft float\n";
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            if(u.x[i][j_longal][k] < 0.0) u_u.x[i][j_longal][k] = 0.0;
            buf << safe_val(u_u.x[i][j_longal][k])
                << " " << safe_val(w_u.x[i][j_longal][k]) << " " << 0.0 << '\n';
        }
    }
    buf << "VECTORS u-w-Downdraft float\n";
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            if(u.x[i][j_longal][k] < 0.0) u_u.x[i][j_longal][k] = 0.0;
            buf << safe_val(u_d.x[i][j_longal][k])
                << " " << safe_val(w_d.x[i][j_longal][k]) << " " << 0.0 << '\n';
        }
    }
    Atmosphere_vtk_longal_File << buf.str();


    Atmosphere_vtk_longal_File.close();
    cout << "   File:  " << Atmosphere_longal_File_Name
        << "  has been written\n";
}

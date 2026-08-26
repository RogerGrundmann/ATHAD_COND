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

#include "cAtmosphereModel.h"
#include "Utils.h"

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

    using namespace ParaViewAtm;

    string Atmosphere_panorama_vts_File_Name = output_path + "/" + Name_Bathymetry_File + "_Atm_sphere_" + std::to_string(n) + ".vts";

    ofstream Atmosphere_panorama_vts_File;

    Atmosphere_panorama_vts_File.precision(8);
    Atmosphere_panorama_vts_File.setf(ios::fixed);

    Atmosphere_panorama_vts_File.open(Atmosphere_panorama_vts_File_Name);

    if(!Atmosphere_panorama_vts_File.is_open()){
        cerr << "ERROR: could not open paraview_sphere_vts file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }

    ostringstream buf;
    buf.precision(8);
    buf.setf(ios::fixed);

    buf << "<?xml version=\"1.0\"?>\n\n"
        << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n\n"
        << " <StructuredGrid WholeExtent=\"" << 1 << " " << im << " " << 1 << " " << jm << " " << 1 << " " << km << "\">\n\n"
        << "  <Piece Extent=\"" << 1 << " " << im << " " << 1 << " " << jm << " " << 1 << " " << km << "\">\n\n"
        << "   <PointData Vectors=\"Velocity\" Scalars=\"Bathymetry Temperature PressureDyn Salinity\">\n\n"
        << "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n\n";

    for(int k = 0; k < km; k++){
        double sinphi = sin(phi.z[k]);
        double cosphi = cos(phi.z[k]);

        for(int j = 0; j < jm; j++){
            double sinthe = sin(the.z[j]);
            double costhe = cos(the.z[j]);
            double sinthe_cosphi = sinthe * cosphi;
            double costhe_cosphi = costhe * cosphi;
            double sinthe_sinphi = sinthe * sinphi;
            double costhe_sinphi = costhe * sinphi;

            for(int i = 0; i < im; i++){
                double ui = u.x[i][j][k];
                double vi = v.x[i][j][k];
                double wi = w.x[i][j][k];

                double au = sinthe_cosphi * ui + costhe_cosphi * vi - sinphi * wi;
                double av = sinthe_sinphi * ui + costhe_sinphi * vi + cosphi * wi;
                double aw = costhe * ui - sinthe * vi;

                aux_u.x[i][j][k] = au;
                aux_v.x[i][j][k] = av;
                aux_w.x[i][j][k] = aw;

                buf << au << " " << av << " " << aw << '\n';
            }
            buf << "\n\n";
        }
        buf << "\n\n";
    }

    buf << "\n\n    </DataArray>\n\n";
    Atmosphere_panorama_vts_File << buf.str();
    buf.str("");
    buf.clear();


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



    dump_array("Bathymetry", h, 1.0, Atmosphere_panorama_vts_File);

    dump_array("u-component", aux_u, u_0, Atmosphere_panorama_vts_File);
    dump_array("v-component", aux_v, u_0, Atmosphere_panorama_vts_File);
    dump_array("w-component", aux_w, u_0, Atmosphere_panorama_vts_File);

    dump_array("PressureDynamic", p_dyn, p_0, Atmosphere_panorama_vts_File);

    dump_array("WaterVapour", c, 1e3, Atmosphere_panorama_vts_File);
    dump_array("CloudWater", cloud, 1e3, Atmosphere_panorama_vts_File);
    dump_array("CloudIce", ice, 1e3, Atmosphere_panorama_vts_File);

    dump_array("Precipitation", Precipitation, 8.64e4, Atmosphere_panorama_vts_File);
    dump_array("PrecipitationRain", P_rain, 8.64e4, Atmosphere_panorama_vts_File);
    dump_array("PrecipitationSnow", P_snow, 8.64e4, Atmosphere_panorama_vts_File);
//    dump_array("PrecipitationGraupel", P_graupel, 8.64e4, Atmosphere_panorama_vts_File);


    buf << "   </PointData>\n\n"
        << "   <Points>\n\n"
        << "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n\n";

    for(int k = 0; k < km; k++){
        double sinphi_k = sin(phi.z[k]);
        double cosphi_k = cos(phi.z[k]);

        for(int j = 0; j < jm; j++){
            double sinthe_j = sin(the.z[j]);
            double costhe_j = cos(the.z[j]);

            for(int i = 0; i < im; i++){
                double r = rad.z[i] + 11.0;
                double r_sinthe = r * sinthe_j;

                buf << r_sinthe * cosphi_k << " "
                    << r_sinthe * sinphi_k << " "
                    << r * costhe_j << '\n';
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

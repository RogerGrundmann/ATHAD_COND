/*
 * Atmosphere General Circulation Modell(AGCM)applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
*/

#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include "Utils.h"
#include "cAtmosphereModel.h"

using namespace std;
using namespace AtomUtils;

void cAtmosphereModel::print_min_max_atm(){

    cout << endl << endl << endl << "      print_min_max_atm" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    cout << endl << endl << endl << " flow properties: " << endl << endl;


    cout << endl << " time step by the Courant-number       dt = 2.8284 * dr/u_0 = 2.8284 * 400/8 = 141.42 s:       dt = " << dt << endl;

    cAtmosphereModel::searchMinMax_3D(" max temperature ", " min temperature ", 
        " deg", t, 273.15, [](double i)->double{return i - 273.15;}, true);
    searchMinMax_3D(" max standard_temp ", " min standard_temp ", 
        " deg", TempStand, 1.0);


    cout << endl << " long wave radiation: " << endl;
    searchMinMax_3D(" max radiation ", " min radiation ", 
        " deg", radiation, 1.0);
    searchMinMax_3D(" max epsilon ", " min epsilon ", 
        " deg", epsilon, 1.0);
    searchMinMax_3D(" max tau_above ", " min tau_above ", 
        " 1 ", tau_above, 1.0);
    searchMinMax_3D(" max tau_layer ", " min tau_layer ", 
        " 1 ", tau_layer, 1.0);
    searchMinMax_3D(" max BruntVaisala_N2 ", " min BruntVaisala_N2 ", 
        " 1/s2 ", brunt_N2, 1.0);
    searchMinMax_3D(" max Psi ", " min Psi ", 
        " kg/s ", Psi, 1.0);
    searchMinMax_3D(" max ubud_pgf ", " min ubud_pgf ", 
        " nd ", ubud_pgf, 1.0);
    searchMinMax_3D(" max ubud_cor ", " min ubud_cor ", 
        " nd ", ubud_cor, 1.0);
    searchMinMax_3D(" max ubud_advv ", " min ubud_advv ", 
        " nd ", ubud_advv, 1.0);
    searchMinMax_3D(" max ubud_advh ", " min ubud_advh ", 
        " nd ", ubud_advh, 1.0);
    searchMinMax_3D(" max ubud_diff ", " min ubud_diff ", 
        " nd ", ubud_diff, 1.0);
    searchMinMax_3D(" max ubud_buoy ", " min ubud_buoy ", 
        " nd ", ubud_buoy, 1.0);
    searchMinMax_2D(" max epsilon_2D ", " min epsilon_2D ", 
         " ppm ", epsilon_2D, 1.0);
    searchMinMax_2D(" max albedo ", " min albedo ", 
         " ppm ", albedo, 1.0);



    cout << endl << " velocity components: " << endl;
    searchMinMax_3D(" max u-component ", " min u-component ", 
        "m/s", u, u_0);
    searchMinMax_3D(" max v-component ", " min v-component ", 
        "m/s", v, u_0);
    searchMinMax_3D(" max w-component ", " min w-component ", 
        "m/s", w, u_0);




    cout << endl << " pressures: " << endl;
    // ATHAD reports pressure in BAR, not hPa. A 250 bar atmosphere printed in hPa is
    // 250000 with the interesting range spanning eight orders of magnitude down to 4e-4;
    // bar is the unit the column profile and the level summary already use, so all four
    // diagnostics now agree. The ParaView fields are renamed Pressure_bar /
    // Pressure_Landscape_bar rather than silently rescaled, so an old saved colour scale
    // fails loudly instead of being wrong by a factor of 1000.
    searchMinMax_3D(" max pressure static ", " min pressure static ", 
        "bar", p_stat, 1.0e-3);
    searchMinMax_3D(" max pressure dynamic ", " min pressure dynamic ", 
        "bar", p_dyn, p_0 * 1.0e-3);




    cout << endl << " density: " << endl;
    searchMinMax_3D(" max density dry air ", " min density dry air ", 
        "kg/m3", r_dry, 1.0);
    searchMinMax_3D(" max density wet air ", " min density wet air ", 
        "kg/m3", r_humid, 1.0);




    cout << endl << " energies: " << endl;
//    searchMinMax_3D(" max radiation ",  " min radiation ", 
//        "kJ/kg", radiation, 1e-3);
    searchMinMax_3D(" max Q_Sensible ", " min Q_Sensible ", 
        "kJ/kg", Q_Sensible, 1e-3);
    searchMinMax_3D(" max Q_Latent ", " min Q_Latent ", 
        "kJ/kg", Q_Latent, 1e-3);




    cout << endl << " cloud thermodynamics: " << endl;
    searchMinMax_3D(" max water vapour ",  " min water vapour ", 
        "g/kg", c, 1e3);
    searchMinMax_3D(" max cloud water ", " min cloud water ", 
        "g/kg", cloud, 1e3);
    searchMinMax_3D(" max cloud ice ", " min cloud ice ", 
        "g/kg", ice, 1e3);
    searchMinMax_3D(" max cloud graupel ", " min cloud graupel ", 
        "g/kg", gr, 1e3);



    cout << endl << " --- parameterization: " << endl;
    searchMinMax_3D(" max S_c_c ", " min S_c_c ", "g/kgs", S_c_c, 1e3);
    searchMinMax_3D(" max S_v ", " min S_v ", "g/kgs", S_v, 1e3);
    searchMinMax_3D(" max S_c ", " min S_c ", "g/kgs", S_c, 1e3);
    searchMinMax_3D(" max S_i ", " min S_i ", "g/kgs", S_i, 1e3);
    searchMinMax_3D(" max S_g ", " min S_g ", "g/kgs", S_g, 1e3);
    searchMinMax_3D(" max S_r ", " min S_r ", "g/kgs", S_r, 1e3);
    searchMinMax_3D(" max S_s ", " min S_s ", "g/kgs", S_s, 1e3);



    cout << endl << " precipitation: " << endl;
    searchMinMax_3D(" max precipitation ", " min precipitation ", "mm/d", Precipitation, 8.64e4);
    searchMinMax_3D(" max precipitation rain ", " min precipitation rain ", "mm/d", P_rain, 8.64e4);
    searchMinMax_3D(" max precipitation snow ", " min precipitation snow ", "mm/d", P_snow, 8.64e4);
    searchMinMax_3D(" max precipitation graup ", " min precipitation graup ", "mm/d", P_graupel, 8.64e4);
    searchMinMax_3D(" max precipitation conv ", " min precipitation conv ", "mm/d", P_conv, 8.64e4);


    cout << endl << " moist convection: " << endl;
    cout << endl << " --- up/downdraft: " << endl;
    searchMinMax_3D(" max E_u ", " min E_u ", "g/m³s", E_u, 1e3);
    searchMinMax_3D(" max D_u ", " min D_u ", "g/m³s", D_u, 1e3);
    searchMinMax_3D(" max M_u ", " min M_u ", "g/m²s", M_u, 1e3);
    cout << endl;


    searchMinMax_3D(" max E_d ", " min E_d ", "g/m³s", E_d, 1e3);
    searchMinMax_3D(" max D_d ", " min D_d ", "g/m³s", D_d, 1e3);
    searchMinMax_3D(" max M_d ", " min M_d ", "g/m²s", M_d, 1e3);


    cout << endl << " --- source terms: " << endl;
    searchMinMax_3D(" max MC_t ", " min MC_t ", "K/s ", MC_t, 1.0);     //  to be complete: (kg/kg)K/s
    searchMinMax_3D(" max MC_q ", " min MC_q ", "g/kgs", MC_q, 1e3);
    searchMinMax_3D(" max MC_v ", " min MC_v ", " m/s²", MC_v, 1.0);
    searchMinMax_3D(" max MC_w ", " min MC_w ", " m/s²", MC_w, 1.0);


    cout << endl << " --- velocity components in the up- and downdraft: " 
        << endl;
    searchMinMax_3D(" max u_u ", " min u_u ", "m/s", u_u, 1.0);
    searchMinMax_3D(" max v_u ", " min v_u ", "m/s", v_u, 1.0);
    searchMinMax_3D(" max w_u ", " min w_u ", "m/s", w_u, 1.0);
    cout << endl;


    searchMinMax_3D(" max u_d ", " min u_d ", "m/s", u_d, 1.0);
    searchMinMax_3D(" max v_d ", " min v_d ", "m/s", v_d, 1.0);
    searchMinMax_3D(" max w_d ", " min w_u ", "m/s", w_d, 1.0);


    cout << endl 
        << " --- condensation and evaporation terms in the up- and downdraft: " 
        << endl;
    searchMinMax_3D(" max c_u ", " min c_u ", "g/kgs", c_u, 1e3);
    searchMinMax_3D(" max e_d ", " min e_d ", "g/kgs", e_d, 1e3);
    searchMinMax_3D(" max e_l ", " min e_l ", "g/kgs", e_l, 1e3);
    searchMinMax_3D(" max e_p ", " min e_p ", "g/kgs", e_p, 1e3);
    searchMinMax_3D(" max g_p ", " min g_p ", "g/kgs", g_p, 1e3);


    cout << endl 
        << " --- water vapour, cloud water and entropies in the up-and downdraft: " 
        << endl;
    searchMinMax_3D(" max q_v_u ", " min q_v_u ", "g/kg", q_v_u, 1e3);
    searchMinMax_3D(" max q_c_u ", " min q_c_u ", "g/kg", q_c_u, 1e3);
    searchMinMax_3D(" max q_v_d ", " min q_v_d ", "g/kg", q_v_d, 1e3);
    searchMinMax_3D(" max s ", " min s ", "./.", s, 1.0);
    searchMinMax_3D(" max s_u ", " min s_u ", "./.", s_u, 1.0);
    searchMinMax_3D(" max s_d ", " min s_d ", "./.", s_d, 1.0);


    cout << endl << endl << endl << " greenhouse gas: " << endl;
    searchMinMax_3D(" max co2 ", " min co2 ", "kg/kg", co2, 1.0);

    // ---- the background, split -----------------------------------------------------------
    // The six background gases are well mixed and source-free, so their mass ratios never
    // change: every per-species field is this constant times the local q_bg, and the split
    // collapses to one effective opacity. At this fork's composition the background is pure
    // N2, so the table below is one line of substance — which is exactly what makes the
    // lumped kappa_bg correct HERE and 1867x too small in ATHAD. Printed so that is visible
    // rather than lucky.
    {
        const double kap[6] = { kappa_N2, kappa_CH4, kappa_NH3, kappa_H2, kappa_CO, kappa_SO2 };
        const char* const* nm = AtmMixture::BG_NAMES();
        const double keff = kappaBackground();
        cout << endl << " background gases, split (q_bg = 1 - q_H2O - q_CO2 locally): " << endl;
        cout << "        species   frac. of bg    kappa [m2/kg]   share of kappa_bg_eff" << endl;
        for (int n = 0; n < 6; n++) {
            const double share = (keff > 0.0) ? 100.0 * m_comp.f_bg[n] * kap[n] / keff : 0.0;
            printf("        %-8s %10.4f %16.3e %18.1f %%\n",
                   nm[n], m_comp.f_bg[n], kap[n], share);
        }
        printf("        %-8s %10.4f %16.3e   (lumped kappa_bg %.3e, ratio %.0fx)\n",
               "TOTAL", 1.0, keff, kappa_bg, (kappa_bg > 0.0) ? keff / kappa_bg : 0.0);
    }
    searchMinMax_3D(" max epsilon ",  " min epsilon ", "%", epsilon, 1.0);


    cout << endl << endl << endl << " forces per unit volume: " << endl;
    searchMinMax_3D(" max presgrad force ", " min pressure force ", 
        "kN", PresGradForce, 1.0);
    searchMinMax_3D(" max buoyancy force ", " min buoyancy force ", 
        "kN", BuoyancyForce, 1.0);
    searchMinMax_3D(" max Coriolis force ", " min Coriolis force ", 
        "N", CoriolisForce, 1.0);
    searchMinMax_3D(" max centrifugal force ", " min centrifugal force ", 
        "N", CentrifugalForce, 1.0);





    cout << endl << endl << endl << " printout of maximum and minimum values of properties at their locations: latitude, longitude" 
        << endl <<
        " results based on two dimensional considerations of the problem" 
        << endl;
    cout << endl << " co2 distribution: " << endl;
    searchMinMax_2D(" max co2 column ", " min co2 column ", 
         " kg/m2 ", co2_total, 1.0);



    cout << endl << " precipitation: " << endl << endl;
    searchMinMax_3D(" max precipitation total ", " min precipitation total ", 
        "mm/d", Precipitation, 8.64e4);
    searchMinMax_2D(" max precipitable NASA ", " min precipitable NASA ", 
        "mm/d", precipitation_NASA, 1.0);
    searchMinMax_2D(" max precipitable water ", " min precipitable water ", 
        "mm", precipitable_water, 1.0);





    cout << endl << " energies at see level without convection influence: " 
        << endl;
//    searchMinMax_2D(" max 2D Q bottom ", " min 2D Q bottom heat ", 
//        "kJ/kg", Q_bottom_2D, 1e-3);
    searchMinMax_2D(" max 2D Q latent ", " min 2D Q latent heat ", 
        "kJ/kg", Q_latent_2D, 1e-3);
    searchMinMax_2D(" max 2D Q sensible ", " min 2D Q sensible heat ", 
        "kJ/kg", Q_sensible_2D, 1e-3);





    cout << endl << " secondary data: " << endl;
    searchMinMax_2D(" max Vegetation ", " min Vegetation ", 
        " ./.", Vegetation, 1.0);
    searchMinMax_2D(" max Evaporation Dalton ", " min Evaporation Dalton ",
        "mm/d", Evaporation_Dalton, 1.0);
    searchMinMax_2D(" max Evaporation Meyer ", " min Evaporation Meyer ",
        "mm/d", Evaporation_Meyer, 1.0);
    searchMinMax_2D(" max Evaporation Rohwer ", " min Evaporation Rohwer ",
        "mm/d", Evaporation_Rohwer, 1.0);
    searchMinMax_2D(" max Evaporation ", " min Evaporation ",
        "mm/d", Evaporation, 1.0);
/*
    searchMinMax_2D(" max 2D albedo ", " min 2D albedo ", 
        "%", albedo, 1.0);
    searchMinMax_2D(" max 2D epsilon_2D ", " min 2D epsilon_2D ", 
        "%", epsilon_2D, 1.0);

    cout << endl << " topography and the vertical cloud extension: " 
        << endl;
    searchMinMax_2D(" max 2D topography ", " min 2D topography ", 
        "m", Topography, 1.0);

    searchMinMax_2D(" max cloud base ",  " min cloud base ", 
        "./.", i_Base, 1);
    searchMinMax_2D(" max cloud top ",  " min cloud top ", 
        "./.", i_LFS, 1);
*/

    if(turb_model != "laminar"){
        cout << endl << " turbulence: " << endl;
        searchMinMax_3D(" max TKE ",          " min TKE ",          "m²/s²",     tke,        1.0);
        searchMinMax_3D(" max dissipation ",  " min dissipation ",  "m²/s³",     dis,        1.0);
        searchMinMax_3D(" max eddy viscosity "," min eddy viscosity ","m²/s",    nue,        1.0);
        searchMinMax_3D(" max turb production "," min turb production ","m²/s³", prod,       1.0);
        searchMinMax_3D(" max TKE source ",   " min TKE source ",   "m²/s³",     tke_source, 1.0);
        searchMinMax_3D(" max dis source ",   " min dis source ",   "m²/s⁴",     dis_source, 1.0);
    }

    cout << endl << endl;

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for print_min_max_atm\n", elapsed.count() * 1e-9);

    cout << "      print_min_max_atm ended" << endl;
}

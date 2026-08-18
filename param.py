# Given a parameter definition, generates necessary C++, Python and XML bindings
# coding=utf-8



def main():

    # read the input definition
    # name, description, datatype, default in the tuples


    PARAMS = {                                                          # dictionary{} (PARAMS) with keys ('common', etc.):[ and their tuples ('output_path', etc.)]
        'common': [
            # Relative to the directory the run is launched from. The ATOM line's convention
            # is to run from python/ and write to python/output_<name>/ — ATOM itself uses
            # output_ATOM/ — so a default of output_Hadean_condensation/ lands in
            # python/output_Hadean_condensation. (The giant-planet siblings hyphenate,
            # output-Uranus/, but ATHAD_COND is from the ATOM line and follows it.)
            ('output_path', 'directory where model outputs should be placed(must end in /)', 'string', 'output_Hadean_condensation/'),



            ('config_xml_path', 'directory where the configuration files are located', 'string', '../python'),

            ('verbose', 'some description of the module', 'bool', False),

#            ('paraview_panorama_vts_flag','flag to control if create paraview panorama', 'bool', False),
            ('paraview_panorama_vts_flag','flag to control if create paraview panorama', 'bool', True),

            ('Coriolis', 'Coriolis force', 'double', 1),
#            ('Coriolis', 'Coriolis force', 'double', 0),

            ('centrifugal', 'centrifugal force', 'double', 1),
#            ('centrifugal', 'centrifugal force', 'double', 0),

            ('buoyancy', 'buoyancy force', 'double', 1),
#            ('buoyancy', 'buoyancy force', 'double', 0),

            ('r_Earth', 'radius of the Earth in km', 'double', 6370.001),
            # ATHAD: the Hadean day was far shorter — tidal braking has been slowing the
            # Earth ever since. A 5.5 h day at ~4.4 Ga gives omega = 2*pi/19800 s =
            # 3.17e-4 rad/s, 4.35x the modern 7.292e-5. ASSUMPTION: estimates for the
            # early rotation range over roughly 4-6 h, and the Rossby number scales with
            # it, so this directly sets how many circulation cells the model can support.
            ('omega', 'ATHAD: rotation rate in rad/s (5.5 h Hadean day)', 'double', 3.17e-4),

            ('g', 'gravitational acceleration of the earth in m/s²', 'double', 9.8066),



            # ATHAD is ONE EPOCH. The paleo-reconstruction inputs that lived here — the
            # topography grids, the NASA surface fields, the Scotese temperature curves, the
            # pygplates reconstruction script, the Ma switches and the hydrosphere SST
            # coupling — are all gone: nothing reconstructs 4.4 Ga, and there is no
            # hydrosphere to couple to. time_start/end/step remain only because the
            # time-slice loop is still structural; a single slice is all that runs.

            ('time_start', 'start time', 'int', 0),
#            ('time_end', 'end time', 'int', 10),
            ('time_end', 'end time', 'int', 0),
            ('time_step', 'step size between timeslices', 'int', 10),



            ('use_earthbyte_reconstruction', 'control whether use earthbyte method to recontruct grids', 'bool', False),
#            ('use_earthbyte_reconstruction', 'control whether use earthbyte method to recontruct grids', 'bool', True),

            ('use_NASA_velocity', 'if use NASA velocity to initialise velocity', 'bool', False),
#            ('use_NASA_velocity', 'if use NASA velocity to initialise velocity', 'bool', True),

            ('use_NASA_temperature', 'if use NASA temperature to initialise surface temperature', 'bool', True),



#            ('CategoryIceScheme', 'number chooses Three(3)-Category Ice Scheme with rain, snow and graupel', 'int', 3),
            ('CategoryIceScheme', 'number chooses Two(2)-Category Ice Scheme with rain, snow', 'int', 2),
#            ('CategoryIceScheme', 'number chooses One(1)-Category Ice Scheme with rain, snow', 'int', 1),
#            ('CategoryIceScheme', 'number chooses Zero(0)-Category Ice Scheme with rain (Warm Rain Scheme)', 'int', 0),
#            ('CategoryIceScheme', 'number chooses no scheme(-1) no precipitation', 'int', -1),

            # ATHAD_COND: 60 bar = 60000 hPa, the midpoint of the 27-100 bar range. NOTE the
            # surface pressure is not actually taken from p_0 — ThermoAtm::densities() builds
            # it as 1e-2*(r_air*R_mix*T), i.e. from the reference density — so r_air below
            # must be consistent with this value. (The old comment here said R_Air, the
            # background gas constant; the code has used the full-mixture R_mix since ATHAD's
            # thermodynamics went in. Corrected rather than carried.)
            ('p_0', 'pressure at sea level in hPa', 'double', 60000.0),
            ('t_0', 'temperature in K compare to 0°C', 'double', 273.15),

            # ATHAD_COND: the sea-surface temperature is PRESCRIBED, as in ATHAD — no
            # paleo-temperature curve reaches this epoch. 513.15 K (240 C) is the midpoint of
            # the stated 230-250 C; the pole is the low end. The 10 K contrast is smaller than
            # ATHAD's 50 K because a liquid ocean is a far better heat reservoir than a magma
            # one, but it is an assumption, not a result.
            #
            # HARD CONSTRAINT, checked in test/cond_column_selftest.cpp: the surface must stay
            # below the boiling point at p_0 (548.7 K at 60 bar) or there is no liquid sea to
            # hold the surface layer at saturation, and the whole Phase 3 boundary condition
            # is void. At the 27 bar end of the input range water boils at 501.2 K, so a
            # 503 K surface there has no ocean at all.
            ('t_surf_equator', 'ATHAD_COND: prescribed sea-surface temperature at the equator in K', 'double', 513.15),
            ('t_surf_pole', 'ATHAD_COND: prescribed sea-surface temperature at the poles in K', 'double', 503.15),

            # ATHAD: the COSMO barometric profile T(h) = T0*sqrt(1 - 2*beta*g*h/(R*T0^2)) has a
            # near-surface lapse rate beta*g/(R*T0), so beta is derived rather than fixed:
            #
            #     beta = cosmo_lapse_fraction * R_mix * T_surf / cp
            #
            # cosmo_lapse_fraction is the lapse as a fraction of the DRY ADIABAT g/cp.
            #
            # Earth's value is 0.51: its observed ~5 K/km sits about half way to the 9.8 K/km
            # dry adiabat. That deficit is not a coincidence or a tuning constant — it is
            # LATENT HEAT. A rising saturated parcel condenses water and releases heat, which
            # partly offsets adiabatic cooling and flattens the lapse toward the moist adiabat.
            #
            # ATHAD's deep column condenses NOTHING: it is supercritical from the surface to
            # ~177 km, so no latent heat is released and there is nothing to flatten the lapse.
            # The physically consistent value is therefore 1.0, the dry adiabat of the mixture,
            # g/cp = 4.81 K/km. Carrying Earth's 0.51 across gave 2.46 K/km, a column that only
            # reached 1279 K at 300 km, never crossed the critical temperature, and so could
            # never condense — a self-fulfilling assumption, since the flattened lapse was
            # itself justified by condensation that the flattening then prevented.
            #
            # On the dry adiabat the column crosses 647 K at 177 km and reaches saturation
            # (p_H2O > p_sat) near 300 K at ~249 km and 0.04 bar, which is where a runaway
            # greenhouse should put its cloud deck.
            ('cosmo_lapse_fraction', 'ATHAD: near-surface lapse rate as a fraction of the dry adiabat; 1.0 = dry adiabat', 'double', 1.0),

            # ATHAD: temperature of the optically thin upper atmosphere, which is radiative
            # rather than convective. Set from the ENERGY BUDGET, not from a measured OLR:
            #
            #     sigma*T_skin^4 = (1 - albedo)*SW + geothermal
            #                    = 0.92*93.5 + 150 = 236.0 W/m2  ->  254.0 K
            #
            # An earlier value of 231.3 K came from T_skin = (OLR/2sigma)^(1/4) using a
            # previously MEASURED OLR, which is circular — the model then reproduced an OLR of
            # sigma*T_skin^4 and confirmed nothing but its own arithmetic.
            #
            # This is now the STARTING value of an iterated fixed point, not a fixed input.
            # With t_skin_relax > 0 the model re-derives it in the loop against its OWN mean
            # albedo — the clear-sky value here is only where the iteration begins.
            #
            # READ THIS BEFORE QUOTING THE OLR. The column's top is isothermal AT t_skin by
            # construction (ThermoAtm::densities re-imposes the profile every iteration) and
            # it is optically thick, so the emission level sits there and the reported OLR is
            # identically sigma*t_skin^4. Measured: t_skin = 254 K -> OLR 236.01 W/m2;
            # t_skin = 240 K -> OLR 188.13 W/m2; sigma*T^4 = 236.01 and 188.13. The OLR is an
            # INPUT wearing an output's clothes. Closing the fixed point below makes it equal
            # the absorbed flux — which is then true by construction, not by test. Making the
            # OLR a genuine prediction requires letting the top find its own temperature
            # radiatively instead of having it prescribed. See the README.
            #
            # It replaces the inherited t_00 = 236.15 K ("-37 C"), which was Earth's
            # tropopause temperature and landed near the right range here by coincidence.
            ('t_skin', 'ATHAD: starting radiative-equilibrium temperature of the optically thin top, in K', 'double', 254.0),

            # ATHAD: under-relaxation of the t_skin fixed-point iteration, per radiation call.
            #
            #     sigma*t_skin^4 = (1 - albedo_mean)*SW_mean + geothermal_flux
            #
            # with albedo_mean the model's OWN cos(latitude)-weighted albedo, built by
            # MultiLayerRadiation from the condensate the model actually made. Set to 0 to
            # hold t_skin at the configured value (the old behaviour, which left the budget
            # open by -32 W/m2). 0.25 converges in a handful of iterations and is well inside
            # the stability limit, the map being a gentle T^(1/4).
            ('t_skin_relax', 'ATHAD: relaxation of the t_skin fixed point per radiation call; 0 = hold t_skin fixed', 'double', 0.25),
            # ATHAD_COND reference density of the MIXTURE at the sea surface, not of dry air:
            # rho = p/(R_mix*T) = 6e6 Pa / (286.6 * 513.15 K) = 40.80 kg/m³ (Earth: 1.2041;
            # ATHAD: 42.97). Almost exactly ATHAD's, at a quarter of the pressure — the gas
            # is twice as heavy and three times as cold, and the two nearly cancel.
            # This is what sets the surface pressure, via p = 1e-2*(r_air*R_mix*T).
            ('r_air', 'ATHAD_COND: reference density of the atmospheric mixture at the sea surface in kg/m³', 'double', 40.80),
            ('r_0_water', 'reference density of fresh water in kg/m3', 'double', 997.0),
            ('t_equat_modern', 'mean temperature of the modern earth in °C', 'double', 15.4),

 
            # ATHAD insolation: the faint young Sun, at the TOP OF THE ATMOSPHERE.
            #
            # At 4.4 Ga the solar constant was about 0.71 of today's 1361 W/m2, i.e.
            # S = 966.3 W/m2, so the global mean incident flux is S/4 = 241.6 W/m2.
            #
            # The previous values, 116 and 71, were ATOM_Precipitation's 163.3 and 100.0
            # scaled by 0.71 — and those are Earth's absorbed-at-the-SURFACE shortwave,
            # already reduced by Earth's albedo and its atmosphere's absorption.
            # MultiLayerRadiation uses these as the flux INCIDENT at the top and then
            # applies ATHAD's own albedo, so Earth's albedo was being counted twice: the
            # cos(latitude)-weighted mean came to 107.5 W/m2, lighting the planet at 45 %
            # of its own insolation.
            #
            # The values below are a fit of the model's parabola short_wave_radiation[j] to
            # the annual-mean insolation of a ZERO-OBLIQUITY planet, Q(phi) = S/pi*cos(phi),
            # constrained so the cos(latitude)-weighted mean is exactly S/4. RMS error
            # 6 W/m2. Hadean obliquity is unknown; zero is the assumption, and it is also
            # the shape the parabola can actually represent — at Earth's 23.44 deg the true
            # curve flattens toward the pole (122 W/m2 there) in a way a parabola cannot
            # follow, giving twice the RMS error. ASSUMPTION — see README.
            ('rad_equator_short', 'ATHAD: TOA short wave insolation at the equator in W/m2', 'double', 298.0),
            ('rad_pole_short', 'ATHAD: TOA short wave insolation at the poles in W/m2', 'double', 0.0),

            # These two are the LONGWAVE boundary values of the inherited scheme. They are
            # present-day Earth fluxes and have no Hadean meaning; mode 2 computes the
            # longwave from the optical depth instead, so they survive only where the old
            # code still reads them.
            ('rad_equator', 'inherited Earth longwave boundary value in W/m2 (unused in radiation_mode 2)', 'double', 398.2),
            ('rad_pole', 'inherited Earth longwave boundary value in W/m2 (unused in radiation_mode 2)', 'double', 360.0),

            ('sigma', 'Stefan-Boltzmann constant W/(m²*K4)', 'double', 5.670280e-8),

            # ==================================================================
            # ATHAD longwave opacity — grey mass absorption coefficients [m2/kg].
            #
            # The layer optical depth is tau = SUM_s kappa_s * q_s * (dp/g) * (p/p_ref),
            # i.e. absorber mass times a pressure-broadening factor. These replace the
            # Bignami/Atwater-Ball emissivity regressions, which are fits to present-day
            # terrestrial columns and saturate to 1 the moment the water path exceeds an
            # Earth-like value.
            #
            # THESE THREE NUMBERS ARE THE BIGGEST SINGLE LEVER ON THE ANSWER and carry
            # roughly a factor-of-two uncertainty. kappa_H2O is the value conventionally
            # used in grey runaway-greenhouse models; CO2 absorbs less per unit mass; the
            # N2/CO/CH4/H2 background is nearly transparent in the thermal infrared. The
            # test of whether they are right is the outgoing longwave flux: a runaway
            # steam atmosphere should sit near the Nakajima / Komabayashi-Ingersoll limit
            # of ~280-310 W/m2, NOT at sigma*T_surf^4.
            ('kappa_H2O', 'ATHAD: grey longwave mass absorption of water vapour in m2/kg', 'double', 0.01),
            ('kappa_CO2', 'ATHAD: grey longwave mass absorption of CO2 in m2/kg', 'double', 0.001),
            ('kappa_bg', 'ATHAD: grey longwave mass absorption of the background gases in m2/kg', 'double', 1.0e-6),

            # ATHAD: geothermal / magma-ocean heat flux through the base of the atmosphere
            # in W/m2. A quenching magma ocean radiates far more than the modern Earth's
            # 0.09 W/m2, and at 1500 K this term is plausibly comparable to the absorbed
            # solar — it cannot be omitted. ASSUMPTION.
            ('geothermal_flux', 'ATHAD: heat flux from the molten surface into the atmosphere in W/m2', 'double', 150.0),

            # ATHAD_COND: near-surface Rayleigh (boundary-layer) drag, made CONFIGURABLE by
            # ATHAD's README item 44 and ported here 2026-08-18. Both were constexpr in
            # RHS_Atm_Turb.cpp with comments justifying them by EARTH'S GEOGRAPHY, which is the
            # pattern this fork was expected to bring and did:
            #
            #   rayleigh_kf   "baseline 1/day gave ~34 m/s eastward w off W-coast S-America;
            #                  10x cut the surface to ~28 m/s ... so 10x is the settled strength"
            #   drag_n_layers "that ~27 m/s max at 27S/71W is an Andes orographic feature"
            #
            # There is no South America and no Andes here either — invariant 1 makes is_land()
            # false everywhere. Worse, drag_n_layers is a COUNT OF CELLS, not a length, so its
            # physical depth is whatever the grid makes it. At this fork's grid:
            #
            #   ATOM_Precipitation (L_atm 400 m,    zeta 3.715, im 41)  5 cells =  236 m
            #   ATHAD_COND         (L_atm 6287.5 m, zeta 3.0,   im 61)  5 cells = 1786 m  -> 7.6x
            #   ATHAD              (L_atm 15719 m,  zeta 3.0,   im 41)  5 cells = 7152 m  -> 30.3x
            #
            # AND IT MOVES WITH im, which ATHAD's writeup did not have to state because it holds
            # im fixed at 41. This fork is at 61, so the SAME constant means a different depth
            # here than next door: 2861 m at im 41, 1786 m at 61, 1297 m at 81 — a 2.2x swing
            # from the level count alone. A grid-refinement study silently rescales the drag.
            # Same defect shape as init_tropopause_layers' round(h / L_atm): a grid index used
            # as a physical length.
            #
            # DEFAULTS REPRODUCE THE OLD constexpr EXACTLY, so this change alone is
            # bit-identical. They exist so the drag can be SCANNED. Reformulating the depth as a
            # LENGTH in metres is the follow-up; it is kept as a cell count here so the default
            # stays exactly the inherited baseline.
            ('rayleigh_kf', 'ATHAD_COND: near-surface Rayleigh drag rate in 1/s (was a constexpr 1/86400; ATHAD item 44)', 'double', 1.0/86400.0),
            ('drag_n_layers', 'ATHAD_COND: boundary-layer drag depth in AIR CELLS, not metres (was a constexpr 5.0; ATHAD item 44 - 5 cells is 1786 m here against 236 m on Earth, and moves with im)', 'double', 5.0),

            ('eps_residuum', 'relative error, end of iterations reached, 1% error  allowed', 'double', 1.0e-4),

#            ('turb_model', 'turbulence model: none, k_epsilon, k_omega, k_omega_SST', 'string', 'none'),
            ('turb_model', 'turbulence model: none, k_epsilon, k_omega, k_omega_SST', 'string', 'k_omega_SST'),
#            ('turb_model', 'turbulence model: none, k_epsilon, k_omega, k_omega_SST', 'string', 'k_omega'),
#            ('turb_model', 'turbulence model: none, k_epsilon, k_omega, k_omega_SST', 'string', 'k_epsilon'),

#            ('inviscid_spinup_iters', 'cumulative iterations to run inviscid (Euler + free-slip mountains) before viscous physics activates; 0 disables', 'int', 100),
#            ('inviscid_spinup_iters', 'cumulative iterations to run inviscid (Euler + free-slip mountains) before viscous physics activates; 0 disables', 'int', 40),
#            ('inviscid_spinup_iters', 'cumulative iterations to run inviscid (Euler + free-slip mountains) before viscous physics activates; 0 disables', 'int', 80),
#            ('inviscid_spinup_iters', 'cumulative iterations to run inviscid (Euler + free-slip mountains) before viscous physics activates; 0 disables', 'int', 300),
            ('inviscid_spinup_iters', 'cumulative iterations to run inviscid (Euler + free-slip mountains) before viscous physics activates; 0 disables', 'int', 0),
            ('inviscid_ramp_iters', 'iterations over which diffusion coefficient ramps from 0 to 1 after the inviscid phase', 'int', 20),

            # ATHAD_COND: 0 — the moist physics runs FROM THE FIRST ITERATION.
            #
            # ATHAD's 300 is defensible there: its column is supercritical from the ground to
            # ~177 km, so the saturation adjustment has almost nothing to do, and delaying it
            # buys a settled circulation before the stiff microphysics starts. Here
            # condensation is the model. A run shorter than 300 iterations with this set to
            # 300 never calls SaturationAdjustment, the ice scheme or MoistConvection at all,
            # and reports a dry-column result that looks like a moist one — the saturated
            # profile is built by densities() either way, so nothing in the printout says the
            # microphysics never ran.
            ('moist_phys_start_iter', 'ATHAD_COND: cumulative iterations before moist physics (SaturationAdjustment, ice scheme, MoistConvection) activates; 0 = always on, which is the ATHAD_COND default because condensation is the subject', 'int', 0),

            ('checkpoint_save_iter', 'dump the full 3D prognostic state to output_path/atm_restart_<iter>.bin when total_iter_count reaches this, for a fast debug restart; -1 disables', 'int', 300),
#            ('checkpoint_save_iter', 'dump the full 3D prognostic state to output_path/atm_restart_<iter>.bin when total_iter_count reaches this, for a fast debug restart; -1 disables', 'int', 200),
            ('restart_from_iter', 'load output_path/atm_restart_<iter>.bin and resume from it, skipping the dry spin-up (debug shortcut); -1 disables', 'int', -1),

            # ATHAD: iterations between full 3D restart dumps. Was a constexpr 100 buried in
            # the iteration loop, next to two configurable siblings (checkpoint for the vtk
            # slices, panorama_print for the vts panoramas) — so the one output you cannot
            # switch off was the one that writes 925 MB per dump. 0 or negative disables.
            ('restart_stride', 'ATHAD: iterations between full 3D restart dumps (.bin, ~925 MB each); 0 disables', 'int', 100),
#            ('restart_from_iter', 'load output_path/atm_restart_<iter>.bin and resume from it, skipping the dry spin-up (debug shortcut); -1 disables', 'int', 300),

#            ('dt_visc', 'non-dimensional time step used in the viscous (production) phase', 'double', 0.001),
#            ('dt_visc', 'non-dimensional time step used in the viscous (production) phase', 'double', 0.0005),
            # The stretched radial grid is fine at the surface (dr ~ 0.014 vs the old
            # constant 0.025); once the radial finite differences correctly use that
            # spacing, the explicit diffusion CFL limit at the surface tightens to
            # dt ~ dr^2/(2D) ~ 1e-4. dt_visc=5e-4 was ~5x over it and blew up the
            # near-surface cells at iter 174. 1e-4 is CFL-safe (validated: passes 174).
            # ATHAD: im 41 -> 61 shrinks the radial step dr from 1/40 to 1/60, and the
            # explicit diffusion CFL limit goes as dr^2 — a factor (40/60)^2 = 0.44. The
            # inherited 1e-4 was validated at im=41 and would be ~2.2x over the limit here,
            # so it is scaled to 4e-5. Verify against a long run before trusting it.
            ('dt_visc', 'non-dimensional time step used in the viscous (production) phase', 'double', 0.00004),
            ('dt_inviscid', 'non-dimensional time step used during the inviscid spin-up phase (smaller to absorb the missing diffusive damping)', 'double', 0.000008),
#            ('dt_inviscid', 'non-dimensional time step used during the inviscid spin-up phase (smaller to absorb the missing diffusive damping)', 'double', 0.0005),
#            ('dt_inviscid', 'non-dimensional time step used during the inviscid spin-up phase (smaller to absorb the missing diffusive damping)', 'double', 0.0003),
        ],


        'atmosphere': [

#            ('nm', 'the maximum number of iterations', 'int', 4),
#            ('nm', 'the maximum number of iterations', 'int', 100),
            ('nm', 'the maximum number of iterations', 'int', 400),
            ('checkpoint', "control when to write output files", 'int', 20),

            # Cadence of the TEXT diagnostics (column profile, level summary, OLR), separate
            # from `checkpoint`, which writes ParaView files. 0 = automatic: every 10
            # iterations for a short run (nm <= 100), every 100 for a longer one.
            ('diagnostic_stride', 'ATHAD: iterations between text diagnostics; 0 = auto (10 short / 100 long)', 'int', 0),
            ('panorama_print', "control when to write panorama files", 'int', 100),


            ('coeff_Dalton', "diffusion coefficient in evaporation by Dalton", 'double', 0.7),

#            ('convection_perturbation', 'convective trigger perturbation: 0=fixed, 1=Bechtold (2008) surface-flux-based for shallow/fixed for deep', 'int', 0),
            ('convection_perturbation', 'convective trigger perturbation: 0=fixed, 1=Bechtold (2008) surface-flux-based for shallow/fixed for deep', 'int', 1),

#            ('convection_mode', 'convection type: 0=deep only (precipitating), 1=deep+shallow (non-precipitating if p_diff<p_stat_diff), 2=deep+shallow+midlevel (also non-precipitating for cloud base above p_stat_midlevel/700 hPa)', 'int', 0),
            ('convection_mode', 'convection type: 0=deep only (precipitating), 1=deep+shallow (non-precipitating if p_diff<p_stat_diff), 2=deep+shallow+midlevel (also non-precipitating for cloud base above p_stat_midlevel/700 hPa)', 'int', 1),

#            ('iter_prec', 'precipitation sub-iteration count: min 3 for evaporation (e_d, e_p) to act on non-zero P_conv; check convergence at 4-5', 'int', 4),
            ('iter_prec', 'precipitation sub-iteration count: min 3 for evaporation (e_d, e_p) to act on non-zero P_conv; check convergence at 4-5', 'int', 3),

            ('evap_model', "evaporation formula driving surface humidity update: Dalton, Meyer, or Rohwer", 'string', 'Meyer'),



            # ==================================================================
            # ATHAD vertical grid.
            #
            # L_atm is the AMPLITUDE of the exponential stretch, NOT the shell thickness
            # and NOT a layer spacing. The shell is (exp(zeta) - 1) * L_atm:
            #     Earth:      (exp(3.715) - 1) *   400.0 =  16.0 km
            #     ATHAD:      (exp(3.000) - 1) * 15719.0 = 300.0 km
            #     ATHAD_COND: (exp(3.000) - 1) *  6287.5 = 120.0 km
            #
            # ATHAD_COND: 120 km, and the reasoning below is ATHAD's, re-run.
            #
            # The scale height here is 15.0 km at the sea against ATHAD's 59.3, so the same
            # number of pressure decades fits in a quarter of the height. The saturated
            # troposphere reaches 66 km (measured, not assumed — it is where the moist
            # adiabat meets t_skin), the grey radiating level tau ~ 1 sits near 0.06 bar at
            # kappa_CO2 = 0.001, which the profile reaches at ~62 km, and 120 km puts the lid
            # at ~1e-5 bar — 6e-7 of the surface pressure, against the 1.5e-6 that made
            # ATHAD's lid transparent.
            #
            # The grid also lands better. With a 300 km shell and im = 61 the troposphere got
            # 33 levels and 27 were spent on isothermal vacuum; at 120 km it gets ~49, and the
            # first layer is 0.32 km instead of 0.81.
            #
            # The shell is set by where the column reaches the radiating level, and the
            # profile has moved under it twice: 300 km on the inherited COSMO shape, 230 km
            # once the proper dry adiabat turned out to compress the atmosphere far more.
            # Finishing the saturation conversion moved it back — with the upper column
            # staying steam it keeps steam's high cp and cools along a shallower adiabat, so
            # 230 km now tops out at 0.29 bar, three times ABOVE the 0.1 bar radiating level,
            # with a top layer of optical depth ~6. An opaque lid emits sigma*T_lid^4 straight
            # out of the domain, which is why the model's OLR is not an output (see t_skin).
            #
            # Deepening it used to be impossible: at 260 km and 300 km the old tridiagonal
            # radiation solve produced NaN across the whole field in its first call, its rows
            # degenerating as eps -> 0. That is fixed (MultiLayerRadiation is now two-stream
            # flux sweeps), and 300 km is what the fix buys:
            #
            #     shell   top p      lid eps   OLR      sigma*T_lid^4
            #     230 km  0.29  bar  1.0000    795 W/m2  787 W/m2   <- OLR == lid, an input
            #     260 km  0.017 bar  0.0610    519       271
            #     300 km  3.8e-4 bar 0.0000    581       271        <- OLR is a real integral
            #
            # At 300 km the top layer is transparent, the isothermal skin is resolved from
            # 256 km up, and the outgoing flux is no longer the boundary temperature read
            # back out. The OLR still varies with the shell depth (519 vs 581), because im is
            # fixed at 61 and a deeper shell is a coarser grid — it is not yet converged.
            #
            # zeta 3.0 with im = 61 keeps the top cell at ~1.7 local scale heights.
            ('L_atm', 'ATHAD_COND: amplitude of the radial stretch in m; shell = (exp(zeta)-1)*L_atm = 120 km', 'double', 6287.5),
            ('zeta', 'ATHAD: radial coordinate-stretching factor (was a hard-coded 3.715)', 'double', 3.0),

            # ATHAD: the radiative-convective boundary of a runaway steam atmosphere sits
            # far higher than Earth's. At 250 km the column is still at 0.68 bar and at
            # 290 km at 0.086 bar, so nearly the whole shell convects. ASSUMPTION — these
            # should be derived from the lapse rate once the radiation is right (Phase 5),
            # not prescribed.
            # ATHAD: the convective column now ends where the adiabat meets the skin
            # temperature, ~207 km. Above that the atmosphere is isothermal and radiative.
            ('tropopause_pole', 'ATHAD: top of the convective column at the poles in m', 'double', 195000.0),
            ('tropopause_equator', 'ATHAD: top of the convective column at the equator in m', 'double', 207000.0),


            # ATHAD_COND circulation-cell layout. Ported from ATHAD (its items 31-38) with
            # the MECHANISM copied and every NUMBER re-derived, because the two models are
            # not in the same regime:
            #
            #                  scale height   dtheta/theta      Ro_T     Held-Hou edge
            #     Earth            8.4 km     45/288  = 0.156   5.96e-2     18.1 deg
            #     ATHAD           59.3 km     50/1500 = 0.033   4.75e-3      5.1 deg
            #     ATHAD_COND      15.0 km     10/513  = 0.0195  7.03e-4      2.0 deg
            #
            # Ro_T = g*H*(dtheta/theta)/(omega^2 a^2). Condensation cost this atmosphere
            # both terms: the sea-surface contrast is 10 K rather than 50, and the column
            # collapsed from 59.3 km of scale height to 15.0. Rotation is unchanged, so
            # Ro_T falls a further 6.8x below ATHAD's and 85x below Earth's, and the direct
            # cell is narrower again -- 2.0 deg against ATHAD's 5.1.
            #
            #   cell_lat_scale = 0.11, from 1.96/18.06 (Held-Hou edge, this model against
            #   Earth's). ATHAD's 0.33 is ATHAD's number and would be 3x too wide here.
            #   CAVEAT: it puts the Hadley EDGE at 30*0.11 = 3.3 deg, which is three points
            #   on a 1 deg grid. The direct cell is marginally resolved and that is a grid
            #   statement, not a physical one -- if it matters, the fix is resolution, not
            #   a wider prescribed cell.
            #
            #   n_cells_hemisphere = 5, from THIS model's own emergent wind. Rhines
            #   L_beta = pi*sqrt(2U/beta) with beta = 2*omega*cos(phi)/a and the measured
            #   U = 20.2 m/s (max zonal at iteration 100) gives a 21.4 deg band at 45 deg,
            #   so (90 - 3.3)/21.4 = 4.1 extratropical bands plus the direct cell -> 5.1.
            #   ATHAD reaches the same count from U = 21.4 and a 9.9 deg Hadley edge, which
            #   is a coincidence of two different roads, not a shared derivation.
            #
            #   cell_amp_mode = 0, as in ATHAD. Unmeasured here.
            #
            # THE SAME CAVEAT APPLIES AS IN ATHAD, and if anything harder: n sets the
            # structure the model is HANDED, not one it can generate. Nothing here
            # maintains an indirect cell. Set cell_lat_scale 1.0 and n_cells_hemisphere 3
            # to recover Earth's layout, which is what everything in this repo before this
            # commit was measured on.
            ('cell_lat_scale', 'ATHAD_COND: Hadley-edge latitude as a fraction of Earth\'s; 0.11 = Held-Hou edge for this rotation rate and contrast, 1.0 = Earth', 'double', 0.11),
            ('cell_amp_mode', 'ATHAD_COND: scale the prescribed velocity amplitudes with cell_lat_scale; 0 = off, 1 = v,w by s (continuity), 2 = v by s and w by s^2 (angular momentum)', 'int', 0),
            ('n_cells_hemisphere', 'ATHAD_COND: prescribed circulation cells per hemisphere; 5 = this regime (Rhines ~21 deg bands), 3 = Earth', 'int', 5),


            # ATHAD albedo. These replace the inherited albedo_pole/albedo_equator, which
            # were INERT — MultiLayerRadiation built its own albedo from bare literals and
            # never read them, so the pole/equator pair was configuration theatre.
            #
            # albedo_surface: a quenching silicate melt is dark. Measured basaltic-melt
            # albedos are 0.05-0.10. There is no land, no snow and no sea ice, so this is a
            # single global clear-sky value, not Earth's latitude parabola.
            #
            # albedo_cloud: the SW albedo of an optically thick cloud top, composited over
            # the surface value by refl = tau/(tau+2) on the condensate path. IT IS
            # CURRENTLY THE WHOLE ANSWER. The reflectivity saturates the moment any
            # condensate is present, and ATHAD's cloud deck is enormously thick, so the
            # model's mean planetary albedo IS this number to four decimal places. It is
            # therefore the second-biggest lever after the opacities, and it is an
            # ASSUMPTION: 0.50 is a thick terrestrial water cloud. A deep, cold, slowly
            # sedimenting Hadean deck could plausibly be brighter.
            ('albedo_surface', 'ATHAD: clear-sky albedo of the molten silicate surface', 'double', 0.08),
            ('albedo_cloud', 'ATHAD: shortwave albedo of an optically thick cloud top', 'double', 0.50),

            ('epsilon_equator', 'emissivity and absorptivity caused by other gases than water vapour/(by Häckel)', 'double', 0.48),
            ('epsilon_pole', 'emissivity and absorptivity caused by other gases than water vapour at the poles', 'double', 0.45),
            ('epsilon_tropopause', 'emissivity and absorptivity caused by other gases than water vapour in the tropopause', 'double', 0.001),

            ('re', 'Reynolds number for laminar flows: ratio viscous to inertia forces, Re = u * L/nue', 'double', 1000.0),
            ('sc_WaterVapour', 'Schmidt number of water vapour, Sc = nue/D', 'double', 0.61),
            ('sc_CO2', 'Schmidt number of CO2', 'double', 0.96),
            ('pr', 'Prandtl number of air for laminar flows', 'double', 0.7179),
            ('pr_turb', 'turbulent Prandtl number for temperature transport in turbulent flows', 'double', 0.9),
            # ATHAD: the boundary layer scales with the scale height, which is 59 km here
            # against Earth's 8.4 km — a ~7x ratio applied to Earth's 1500 m. At the chosen
            # grid this puts ~12 cells inside the ABL. ASSUMPTION.
            ('abl_height', 'ATHAD: physical depth of the atmospheric boundary layer in m', 'double', 10000.0),
            # ==================================================================
            # ATHAD_COND atmospheric composition — MOLE fractions AT THE SEA SURFACE.
            #
            # These are NOT the numbers the literature quotes for this epoch. The quoted
            # composition — H2O 0.4-2 %, CO2 89-95 %, N2 5-20 % — is the DRY atmosphere,
            # what survives above the cold trap. The air in contact with a 240 C ocean
            # carries water at its own vapour pressure instead: p_sat(513.15 K) = 33.47 bar
            # of a 60 bar column, so x_H2O = 0.5578 at the surface, fifty-five times the
            # dry figure, and CO2 and N2 are diluted in their dry ratio to fit.
            #
            # The two are consistent once separated by height, and the column runs between
            # them; test/cond_column_selftest.cpp asserts both mixtures and the arithmetic
            # that connects them. Change these and that test fails, which is the point.
            #
            # The five trace gases of ATHAD's reducing mixture are gone: the epoch this
            # models is oxidised and degassed, and CH4/NH3/H2/CO/SO2 are not part of the
            # stated composition. They stay as parameters, at zero, because MixtureAtm's
            # background pseudo-species still needs somewhere to put anything added later.
            # ==================================================================
            ('x_H2O', 'ATHAD_COND: mole fraction of H2O at the sea surface (saturated: p_sat(T_surf)/p_0)', 'double', 0.5578),
            ('x_CO2', 'ATHAD_COND: mole fraction of CO2 at the sea surface', 'double', 0.4109),
            ('x_N2',  'ATHAD_COND: mole fraction of N2 at the sea surface',  'double', 0.0313),
            ('x_CH4', 'ATHAD_COND: mole fraction of CH4', 'double', 0.0),
            ('x_NH3', 'ATHAD_COND: mole fraction of NH3', 'double', 0.0),
            ('x_H2',  'ATHAD_COND: mole fraction of H2',  'double', 0.0),
            ('x_CO',  'ATHAD_COND: mole fraction of CO',  'double', 0.0),
            ('x_SO2', 'ATHAD_COND: mole fraction of SO2', 'double', 0.0),

            # ep = R_background / R_H2O = 296.8/461.5. NOTE: the dilute approximation this
            # constant serves, q_sat = ep*E/(p-(1-ep)*E), is INVALID here for the same reason
            # it was invalid in ATHAD and NOT for the opposite one: water is 0.4 % by mole
            # above the cold trap, where the dilute form would be fine, but 56 % at the sea
            # surface, where it is not. A form that is valid over part of the column is worse
            # than one that is valid over none, because it looks right in the printouts.
            # The exact mass-fraction form is used everywhere; ep remains only where a
            # genuine gas-constant ratio is wanted.
            ('ep', 'ATHAD_COND: ratio of the background-mixture to water-vapour gas constants', 'double', 0.6431),
            ('hp', 'water vapour pressure at T = 0°C: E = 6.1 hPa', 'double', 6.1078),

            # ATHAD_COND: "Air" means the NON-CONDENSABLE BACKGROUND (everything but H2O and
            # CO2). With the five reducing trace gases gone that background is N2 alone:
            # M_bg = 28.014 g/mol -> R_bg = 296.8 J/(kg K). It is not air, and the fact that
            # it now coincides with the gas constant of N2 is arithmetic, not a return to
            # Earth. The quantity that matters for saturation is M_nonwater (CO2 + background,
            # 42.88 g/mol here), which is a different number and a different function.
            ('R_Air', 'ATHAD_COND: specific gas constant of the non-condensable background in J/(kg*K)', 'double', 296.8),
            ('R_WaterVapour', 'specific gas constant of water vapour in J/(kg*K)', 'double', 461.5),
            ('r_water_vapour', 'density of saturated water vapour in kg/m³ at 10°C', 'double', 0.0094),
            ('R_co2', 'specific gas constant of CO2 in J/(kg*K)', 'double', 188.9),
            ('lv', 'specific latent evaporation heat(condensation heat) in J/kg', 'double', 2.52e6),
            ('ls', 'specific latent vaporisation heat(sublimation heat) in J/kg', 'double', 2.83e6),

            # ATHAD_COND: cp of the MIXTURE at the sea surface, 1349 J/(kg K) from
            # MixtureAtm::cp_of at 513 K on the saturated composition — one third of ATHAD's
            # 2040 and, by coincidence, close to Earth's 1005. It varies by 25 % across the
            # column, from 1349 at the wet surface to 1028 in the dry CO2 air aloft, so this
            # constant is only the fallback; cp_of() gives the local value.
            # cv_l = cp_l - R_mix = 1349 - 286.6.
            ('cp_l', 'ATHAD_COND: specific heat capacity of the mixture at constant pressure in J/(kg K)', 'double', 1349.0),
            ('cv_l', 'ATHAD_COND: specific heat capacity of the mixture at constant volume in J/(kg K)', 'double', 1062.4),
            ('lamda', 'heat transfer coefficient of air in W/(m K)', 'double', 0.0262),
            ('r_co2', 'density of CO2 in kg/m³ at 25°C', 'double', 0.0019767),
            ('gam', 'constant slope of temperature    gam = 6.5 K/1000 m', 'double', 0.0065),

            ('u_0', 'annual mean of surface wind velocity in m/s, 8 m/s compare to 28.8 km/h', 'double', 8.0),
            # ATHAD: an upper PHYSICAL bound on the prognostic temperature, replacing the
            # hard-coded 333.15 K (60 °C) literal in SaturationAdjustment. That literal was
            # justified as "well above any physical surface temperature" — true on Earth,
            # but a factor of 4.5 BELOW ATHAD's 1500 K surface, and it was written back into
            # the prognostic field, collapsing 250 bar to 30 bar in one iteration.
            # 2000 K leaves headroom over the surface and stays inside the cp Shomate fits.
            # ATHAD: ceiling on the condensate mass fraction, applied in two places that
            # were both bare literals — SaturationAdjustment::clampAndFade caps cloud and
            # ice SEPARATELY at this value, and the moist-physics block in RunTimeSlice caps
            # their SUM at it before densities() reads them. Kept as one number so they
            # cannot drift apart; the sum-cap is the binding one.
            #
            # The inherited comment called 0.05 "~50x the largest physical cloud/ice mixing
            # ratio, so it never clips a real cloud — it only stops a runaway". That is an
            # EARTH statement: terrestrial cloud water is ~1 g/kg. ATHAD's atmosphere is
            # 67 % water by mass and its polar column pegs this cap exactly (49.999996 g/kg
            # measured at 80N, 185.6 km), so here it is not a runaway backstop at all — it
            # is setting the condensate, and with it the optical depth and the albedo.
            # Raise it and re-measure before trusting a cloud field that sits on it.
            ('cloud_cap', 'ATHAD: ceiling on the condensate mass fraction (cloud, ice, and their sum) in kg/kg', 'double', 0.05),

            ('t_max_phys', 'ATHAD: upper physical bound on the prognostic temperature in K', 'double', 2000.0),

            ('t_00', 'temperature in K compare to -37°C', 'double', 236.15),
            ('t_000', 'temperature in K compare to -20°C', 'double', 253.15),
            ('s_0', 'entropy at 0°C, cp_l * t_0 in J/kg', 'double', 274515.75),
            # ATHAD_COND: the water-vapour scale is the SEA-SURFACE mass fraction
            # q_H2O = 0.3464 (saturated at 513 K / 60 bar), not Earth's 0.035 trace and not
            # ATHAD's 0.6724. c_0 is a normalisation in the RHS energy/moisture coefficients
            # (RHS_Atm_Turb.cpp: coeff_energy, coeff_MC_q, coeff_L) and, unlike in ATHAD, it
            # is NOT the value the field is initialised to: the water field is a profile
            # here, running from c_0 at the sea to c_h2o_dry_top above the cold trap.
            ('c_0', 'ATHAD_COND: reference water vapour mass fraction at the sea surface in kg/kg', 'double', 0.3464),

            # ATHAD_COND: the water that survives above the cold trap, as a MASS fraction.
            # This is the model's link to the quoted composition: the literature's "H2O
            # 0.4-2 % by mole" is the DRY atmosphere, x_H2O = 0.010 of a CO2/N2 mixture of
            # 42.63 g/mol, i.e. q = 0.00423 kg/kg. initWaterWapour() uses it as the floor
            # under the saturation profile, and the height at which the two cross IS the
            # cold trap — a computed quantity, not a prescribed level.
            ('c_h2o_dry_top', 'ATHAD_COND: water vapour mass fraction above the cold trap in kg/kg', 'double', 0.00423),

            # ATHAD_COND: the CO2 field is a MASS FRACTION, not ppm. At 62% by mass at the
            # sea surface, ppm is meaningless. This epoch has no biosphere and no vegetation;
            # it does have an ocean, and therefore in reality a carbonate sink — but that sink
            # is a hydrosphere process and there is no hydrosphere here, so CO2 remains well
            # mixed and conserved, and that is an ASSUMPTION this model makes rather than a
            # result it derives. q_CO2 = 0.4109*44.010/29.009 = 0.6233 at the surface.
            ('co2_0', 'ATHAD_COND: reference CO2 mass fraction at the sea surface in kg/kg', 'double', 0.6233),
            ('co2_scale', 'multiplier applied to the whole CO2 field for sensitivity experiments (1.0 = field as built; 2.0 = doubled CO2)', 'double', 1.0),

            # ATHAD: no land, so no land/ocean humidity split — a single surface relative
            # humidity applies everywhere. Kept as one knob rather than two identical ones.
            ('c_ocean', 'ATHAD: surface water vapour as a fraction of the saturation value, in %', 'double', 70),

        ],
    }




    XML_READ_FUNCS = {                                                  # dictionary (XML_READ_FUNCS) with keys ('"string", etc.) and one list element ("FillStringWithElement")
        "string": "FillStringWithElement",
        "double": "FillDoubleWithElement",
        "int": "FillIntWithElement",
        "bool": "FillBoolWithElement"
    }



    # functions begin

    def write_cpp_defaults(filename, classname, sections):

        with open(filename, 'w') as f:

            f.write("// header files\n")
            f.write("// THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n")
            f.write("// ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n")
            f.write("\n")
            f.write("void %s::SetDefaultConfig() {\n" % classname)

            for section in sections:
                f.write('\n  // %s section\n' % section)

                for slug, desc, ctype, default in PARAMS[section]:
                    rhs = default

                    if ctype == 'string':
                        rhs = '"%s"' % default
                    elif ctype == 'bool':

                        if default:
                            rhs = 'true'
                        else:
                            rhs = 'false'

                    f.write('  %s = %s;\n' %(slug, rhs))

            f.write("}")



    def write_cpp_load_config(filename, classname, sections):

        with open(filename, 'w') as f:

            f.write("// config files\n")
            f.write("// THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n")
            f.write("// ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n")
            f.write("\n")

            for section in sections:
                f.write('\n  // %s section\n' % section)
                element_var_name = 'elem_%s' % section
                f.write('\n  if(%s) {\n' %(element_var_name))

                for slug, desc, ctype, default in PARAMS [section]:
                    func_name = XML_READ_FUNCS [ctype]
                    f.write('    Config::%s(%s, "%s", %s);\n' %(func_name, element_var_name, slug, slug))

                f.write("  }\n")



    def write_cpp_headers(filename, sections, is_extern = False):

        with open(filename, 'w') as f:

            f.write("// header files\n")
            f.write("// THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n")
            f.write("// ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n")
            f.write("\n")

            if is_extern:
                f.write("#include<string>\n\n")
                f.write("using namespace std;\n")

                #if 'atmosphere' in filename:
                #    f.write("namespace AtmParameters{\n")
                #else:
                #    f.write("namespace HydParameters{\n")

            for section in sections:
                f.write('\n// %s section\n' % section)

                for slug, desc, ctype, default in PARAMS [section]:
                    if is_extern:
                        f.write('   extern %s %s;\n' %(ctype, slug))
                    else:
                        f.write('%s %s;\n' %(ctype, slug))
           
            if is_extern:
                f.write("}\n")



    def write_pxi(input_filename, output_filename, substitutions):

        data = open(input_filename, 'r').read()
        indent = '    '

        for key, classname, sections in substitutions:
            rep = ''

            for section in sections:
                rep += '%s# %s section\n' %(indent, section)

                for slug, desc, ctype, default in PARAMS[section]:
                    rep += '%sproperty %s:\n' %(indent, slug)
                    rep += '%s    def __get__(%s self):\n' %(indent, classname)
                    rep += '%s        self._check_alive()\n' % indent
                    rep += '%s        return self._thisptr.%s\n' %(indent, slug)
                    rep += '%s\n' % indent
                    rep += '%s    def __set__(%s self, value):\n' %(indent, classname)
                    rep += '%s        self._check_alive()\n' % indent
                    rep += '%s        self._thisptr.%s = <%s> value\n' %(indent, slug, ctype)
                    rep += '%s\n' % indent

            data = data.replace('{{ %s }}' % key, rep)

        with open(output_filename, 'w') as f:

            f.write("""# pxi files\n""")
            f.write("# THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n")
            f.write("# ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n")
            f.write(data)



    def write_pxd(filename, model, sections):

        with open(filename, 'w') as f:
            # Sadly, Cython docs are incorrect on usage of 'include', so we must include a whole lot of boilerplate

            f.write("""# pxd files\n""")
            f.write("""# THIS FILE IS AUTOMATICALLY GENERATED BY param.py
# ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME
from libcpp.vector cimport vector
cdef extern from "c%sModel.h":
    cppclass c%sModel:
        c%sModel() except +  # NB! std::bad_alloc will be converted to MemoryError
        void LoadConfig(const char *filename)
        void Run()
        void RunTimeSlice(int time_slice)
        vector[float] get_layer_heights()
""" %(model, model, model))

            for section in sections:
                f.write('        # %s section\n' % section)

                for slug, desc, ctype, default in PARAMS [section]:
                    f.write('        %s %s\n' %(ctype, slug))



    def write_config_xml(filename, sections):

        with open(filename, 'w') as f:
            f.write("""<!-- THIS FILE IS GENERATED AUTOMATICALLY BY param.py. DO NOT EDIT. -->""")
            f.write('<atom>')

            for section in sections:
                f.write('    <%s>\n' % section)

                for slug, desc, ctype, default in PARAMS [section]:
                    if ctype == 'bool':
                        default = str(default).lower()                  # Python uses True/False, C++, uses true/false
                    f.write('        <%s>%s</%s>  <!-- %s(%s) -->\n' %(slug, default, slug, desc, ctype))

                f.write('    </%s>\n' % section)
 
            f.write('</atom>')

    # functions end



    atmosphere_sections = ['common', 'atmosphere']

    for filename, classname, sections in [
        ('atmosphere/cAtmosphereDefaults.cpp.inc', 'cAtmosphereModel', atmosphere_sections)
    ]:
        write_cpp_defaults(filename, classname, sections)


    for filename, classname, sections in [
        ('atmosphere/AtmosphereLoadConfig.cpp.inc', 'cAtmosphereModel', atmosphere_sections)
    ]:
        write_cpp_load_config(filename, classname, sections)


    for filename, sections in [
        ('atmosphere/AtmosphereParams.h.inc', atmosphere_sections)
    ]:
        write_cpp_headers(filename, sections)


    write_pxi('python/pycond.pyx.template', 'python/pycond.pyx', [
        ('atmosphere_params', 'Atmosphere', atmosphere_sections)])


    for filename, model, sections in [
        ('python/atmosphere_pxd.pxi', 'Atmosphere', atmosphere_sections)
    ]:
        write_pxd(filename, model, sections)


    for  filename, sections in [
        ('python/config_cond.xml', atmosphere_sections)
    ]:
        write_config_xml(filename, sections)


    for  filename, sections in [
        ('cli/config_cond.xml', atmosphere_sections)
    ]:
        write_config_xml(filename, sections)



if __name__ == '__main__':
    main()

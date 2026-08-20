// Standalone self-test for the ATHAD_COND regime — the post-condensation Hadean column.
//
// ATHAD's numbers were checked once and then lived in CLAUDE.md; here they are asserted,
// because ATHAD_COND inverts the thermodynamic regime and almost every one of them moves.
// This test owns the composition arithmetic that the config values are derived FROM, so a
// change to either has to be made in both places and the disagreement shows up here.
//
// Build and run:
//   g++ -std=c++17 -Iatmosphere -Ilib -o test/cond_column_selftest test/cond_column_selftest.cpp
//   ./test/cond_column_selftest
//
// or from the project root:  make test-cond-column

#include "MixtureAtm.h"
#include "SaturationH2O.h"

#include <cstdio>
#include <cmath>

static int failures = 0;

static void check(const char* what, double got, double want, double rel_tol)
{
    const double rel = std::fabs(got - want) / std::fabs(want);
    const bool ok = (rel <= rel_tol);
    if (!ok) failures++;
    std::printf("  %-54s got %12.6g   want %12.6g   %s\n",
                what, got, want, ok ? "ok" : "FAIL");
}

static void check_true(const char* what, bool cond)
{
    if (!cond) failures++;
    std::printf("  %-54s %s\n", what, cond ? "ok" : "FAIL");
}

// The nominal run point. These four numbers are the model's inputs; everything else in
// this file is derived from them, and the config must agree.
static const double P_SURF_BAR = 60.0;
static const double T_SURF     = 513.15;      // 240 C, the midpoint of 230-250 C
static const double T_POLE     = 503.15;      // 230 C

// The DRY composition — what the atmosphere is above the cold trap, and the figure the
// literature quotes: H2O 0.4-2 %, CO2 89-95 %, N2 5-20 %. Midpoints, renormalised.
//
// These stay as the literature REFERENCE and as the source of the CO2:N2 ratio the surface
// mixture inherits. They are no longer the whole dry atmosphere: since 2026-08-20 the five
// trace gases are present at ATHAD's 1.4 mole-% each, and above the cold trap — where the
// water is gone — they are 15.7 % of what is left, which pushes CO2 to 77.4 % and OUT of the
// quoted 89-95 % band. That is a deliberate departure from the quoted composition; see the
// composition block in param.py.
static const double X_H2O_DRY  = 0.010;
static const double X_CO2_DRY  = 0.920;
static const double X_N2_DRY   = 0.070;

// The five trace gases, at ATHAD's mole fraction. Non-condensable and well mixed, so the
// same fraction applies at every level; what changes with height is everyone else's share.
static const double X_TRACE    = 0.014;       // each of CH4, NH3, H2, CO, SO2
static const double X_TRACE_5  = 5.0 * X_TRACE;

// The SURFACE composition, copied from param.py. These are the config's x_H2O/x_CO2/x_N2
// and r_air; the checks below re-derive them from T_SURF and P_SURF_BAR and fail if the
// config has drifted from the physics it is supposed to encode.
static const double X_H2O_CFG  = 0.5578;
static const double X_CO2_CFG  = 0.3459;
static const double X_N2_CFG   = 0.0263;
static const double R_AIR_CFG  = 39.08;       // [kg/m3]

int main()
{
    std::printf("\nATHAD_COND column self-test  (%.0f bar, %.2f K)\n\n", P_SURF_BAR, T_SURF);

    // ------------------------------------------------------------------
    // 1. The regime itself. ATHAD's invariant 2 is that water is supercritical; this
    //    model exists because it is not. If these fail, nothing below means anything.
    // ------------------------------------------------------------------
    std::printf(" the regime is subcritical and the sea is liquid\n");
    check_true("T_surf below the critical temperature",
               T_SURF < AtmMixture::T_CRIT_H2O);
    check_true("p_surf below the critical pressure",
               P_SURF_BAR * 1000.0 < AtmMixture::P_CRIT_H2O);

    // Boiling point at the surface pressure: the sea must not be boiling, or there is no
    // liquid ocean to hold the surface layer at saturation.
    const double T_boil = SaturationH2O::dewPoint(P_SURF_BAR * 1000.0);
    check("boiling temperature at 60 bar",            T_boil, 548.7, 2e-3);
    check_true("the sea is below boiling at the equator",  T_SURF < T_boil);
    check_true("the sea is below boiling at the pole",     T_POLE < T_boil);
    // For the record, and because the 27 bar corner of the input range fails this:
    std::printf("      (at 27 bar water boils at %.1f K, so a 503 K surface there is dry)\n",
                SaturationH2O::dewPoint(27000.0));

    // ------------------------------------------------------------------
    // 2. The surface is SATURATED. This is the modelling decision that separates
    //    ATHAD_COND from a dry-surface variant: there is an ocean, so the air in contact
    //    with it carries water at its own vapour pressure — not at the dry-atmosphere
    //    0.4-2 %, which is what survives above the cold trap.
    // ------------------------------------------------------------------
    std::printf("\n the surface mole fraction is set by the sea, not by the dry composition\n");
    const double p_sat_surf = SaturationH2O::saturationPressure(T_SURF);        // [hPa]
    check("p_sat(T_surf) [bar]",                      p_sat_surf / 1000.0, 33.470, 5e-3);

    const double x_H2O_surf = p_sat_surf / (P_SURF_BAR * 1000.0);
    check("saturated x_H2O at the surface",           x_H2O_surf, 0.5578, 5e-3);
    check_true("which is 50x the dry-atmosphere value",
               x_H2O_surf > 50.0 * X_H2O_DRY);

    // The rest of the surface column is the five trace gases at their fixed 1.4 % each,
    // and then CO2 and N2 in the dry ratio within whatever is left. The order matters and
    // is the physics: x_H2O is not free — it is p_sat(T_surf)/p_0, set by the sea — so
    // anything added to the mixture displaces CO2 and N2, never water.
    const double x_rest      = 1.0 - x_H2O_surf - X_TRACE_5;
    const double dry_non_h2o = X_CO2_DRY + X_N2_DRY;
    const double x_CO2_surf  = x_rest * X_CO2_DRY / dry_non_h2o;
    const double x_N2_surf   = x_rest * X_N2_DRY  / dry_non_h2o;
    check("saturated x_CO2 at the surface",           x_CO2_surf, 0.3459, 1e-2);
    check("saturated x_N2  at the surface",           x_N2_surf,  0.0263, 2e-2);

    // ...and the config carries exactly those, to the four decimals it stores them in.
    check("config x_H2O matches the saturated value", X_H2O_CFG, x_H2O_surf, 1e-3);
    check("config x_CO2 matches the saturated value", X_CO2_CFG, x_CO2_surf, 2e-3);
    check("config x_N2  matches the saturated value", X_N2_CFG,  x_N2_surf,  2e-3);
    check_true("config mole fractions sum to 1",
               std::fabs(X_H2O_CFG + X_CO2_CFG + X_N2_CFG + X_TRACE_5 - 1.0) < 1e-9);

    // ------------------------------------------------------------------
    // 3. The two mixtures, from the model's own resolve(). The config's x_* and r_air are
    //    the SURFACE values; the dry mixture is what the column tends to aloft.
    // ------------------------------------------------------------------
    // Resolved from the CONFIG values, not the exact ones, so what this test reports is
    // what the model will actually run with.
    std::printf("\n mixture properties at the sea surface\n");
    const AtmMixture::Composition S =
        AtmMixture::resolve(X_H2O_CFG, X_CO2_CFG, X_N2_CFG,
                            X_TRACE, X_TRACE, X_TRACE, X_TRACE, X_TRACE);
    check_true("mole fractions sum to 1",             S.valid);
    check("M_mean [g/mol]",     S.M_mean * 1000.0,    27.79,  1e-2);
    check("R_mix [J/(kg K)]",   S.R_mix,              299.2,  1e-2);
    check("q_H2O [kg/kg]",      S.q_H2O,              0.3616, 2e-2);

    // The background is no longer N2: it is N2 plus the five traces, and by MASS the traces
    // are 71 % of it. R_bg lands within 0.3 % of ATHAD's 317.26 for that reason.
    check("M_bg [g/mol]",       S.M_bg * 1000.0,      26.138, 1e-3);
    check("R_bg [J/(kg K)]",    S.R_bg,               318.10, 1e-3);
    check_true("N2 is a minority of the background by mass", S.f_bg[0] < 0.35);
    check_true("the six background fractions sum to 1",
               std::fabs(S.f_bg[0]+S.f_bg[1]+S.f_bg[2]+S.f_bg[3]+S.f_bg[4]+S.f_bg[5]-1.0) < 1e-9);

    const double cp_surf = AtmMixture::cp_of(S.q_H2O, S.q_CO2, T_SURF, S.M_bg);
    check("cp(T_surf) [J/(kg K)]", cp_surf,           1372.6, 3e-2);

    // r_air is the column's anchor: ThermoAtm::densities() sets p_stat[0] = 1e-2*r_air*
    // R_mix*T_surf, so r_air and R_mix together ARE the 60 bar. Getting this wrong does not
    // fail loudly — it quietly runs a different atmosphere.
    const double rho_surf = P_SURF_BAR * 1.0e5 / (S.R_mix * T_SURF);
    check("rho_surf [kg/m3]",                         rho_surf, 39.08, 2e-2);
    check("config r_air reproduces it",               R_AIR_CFG, rho_surf, 1e-3);
    check("...and so the anchor gives back p_0 [hPa]",
          1e-2 * R_AIR_CFG * S.R_mix * T_SURF, P_SURF_BAR * 1000.0, 1e-3);
    const double H_surf = S.R_mix * T_SURF / 9.81;
    check("scale height [m]",   H_surf,               15650.0, 3e-2);

    std::printf("\n mixture properties above the cold trap (the dry composition)\n");
    // Built from the SURFACE mixture with the water taken out, not from the literature
    // midpoints: the trace gases do not condense, so above the cold trap they keep their
    // partial pressures while water loses all but x = 0.010 of its own, and everyone else's
    // fraction rises to fill the gap. CO2 ends at 77.4 %, below the quoted 89-95 %.
    const double dry_scale = 0.990 / (1.0 - X_H2O_CFG);
    const AtmMixture::Composition D =
        AtmMixture::resolve(X_H2O_DRY, X_CO2_CFG * dry_scale, X_N2_CFG * dry_scale,
                            X_TRACE * dry_scale, X_TRACE * dry_scale, X_TRACE * dry_scale,
                            X_TRACE * dry_scale, X_TRACE * dry_scale);
    check_true("mole fractions sum to 1",             D.valid);
    check("M_mean [g/mol]",     D.M_mean * 1000.0,    39.897, 1e-3);
    check("R_mix [J/(kg K)]",   D.R_mix,              208.40, 1e-3);
    check("q_H2O [kg/kg]",      D.q_H2O,              0.00452, 2e-2);
    check("dry x_CO2 is out of the quoted 0.89-0.95 band",
          X_CO2_CFG * dry_scale, 0.7744, 1e-2);

    // ------------------------------------------------------------------
    // 4. Why the dry adiabat cannot be used here.
    //
    //    ATHAD sets cosmo_lapse_fraction = 1.0 and integrates dT/dz = -g/cp, and that is
    //    correct there BECAUSE nothing condenses. The saturated lapse carries the extra
    //    heat capacity L*dq_sat/dT, and at 35 % water by mass that term is an order of
    //    magnitude larger than cp itself. This check is the justification for the moist
    //    adiabat in ThermoAtm::densities(); if it ever stops holding, that change is no
    //    longer needed.
    // ------------------------------------------------------------------
    std::printf("\n the saturated lapse rate is nothing like the dry one\n");
    // The "other" molar mass in the saturation formula is everything that is NOT water,
    // i.e. CO2 AND the background — M_nonwater, not M_bg. Passing M_bg here (N2 only,
    // 28.014) makes the sea surface come out at q_sat = 0.448 instead of 0.346, because
    // it weighs the non-water gas as if the 41 % CO2 were not there. That is a 29 % error
    // in the single quantity this whole model turns on, and it is the same confusion
    // CLAUDE.md records M_nonwater itself having had.
    const double M_nw_surf = AtmMixture::M_nonwater(S.q_H2O, S.q_CO2, S.M_bg);
    check("M_nonwater at the surface [g/mol]",        M_nw_surf * 1000.0, 40.118, 1e-2);

    const double q_sat_surf = SaturationH2O::saturationMassFractionAt(
                                  T_SURF, P_SURF_BAR * 1000.0, M_nw_surf);
    check("q_sat at the sea surface [kg/kg]",         q_sat_surf, S.q_H2O, 5e-2);

    const double L_surf = SaturationH2O::latentHeat(T_SURF);
    check("latent heat at 513 K [J/kg]",              L_surf, 1.72e6, 5e-2);

    // NOTE: dqSatdT is the DILUTE Clausius-Clapeyron form, q_sat*L/(R_v T^2). At q_sat =
    // 0.35 it is not the derivative of the exact mass fraction — see the exact form below,
    // which differs by M_other/(x*M_H2O + (1-x)*M_other) ~ 1.45 here. Both are printed so
    // the gap is visible; Phase 2 makes the exact one the one the model uses.
    const double dqdT_dilute = q_sat_surf * L_surf / (AtmMixture::R_H2O * T_SURF * T_SURF);
    const double dqdT_exact  = SaturationH2O::dqSatdT(T_SURF, P_SURF_BAR * 1000.0, M_nw_surf);
    std::printf("      dq_sat/dT  dilute %.5g /K   exact %.5g /K   ratio %.3f\n",
                dqdT_dilute, dqdT_exact, dqdT_exact / dqdT_dilute);
    check("the exact/dilute ratio is M_other/den",
          dqdT_exact / dqdT_dilute,
          M_nw_surf / (x_H2O_surf * AtmMixture::M_H2O + (1.0 - x_H2O_surf) * M_nw_surf), 1e-2);

    // The lapse rate itself. NOTE the trap this check was written wrong the first time:
    // the familiar g/(cp + L dq/dT) keeps only the denominator of the moist-static-energy
    // balance and gives 0.7 K/km here, which is seven times too small. The pressure
    // dependence of q_sat is not negligible when water is 56 % by mole.
    const double gamma_dry     = 9.81 / cp_surf;
    const double gamma_trunc   = 9.81 / (cp_surf + L_surf * dqdT_exact);
    const double gamma_moist   = SaturationH2O::moistLapse(T_SURF, P_SURF_BAR * 1000.0,
                                     M_nw_surf, cp_surf, S.R_mix, 9.81);
    std::printf("      dry %.3f K/km   saturated %.3f K/km   (truncated form would give %.3f)\n",
                gamma_dry * 1000.0, gamma_moist * 1000.0, gamma_trunc * 1000.0);
    check("dry lapse [K/km]",                         gamma_dry * 1000.0,   7.15, 1e-2);
    check("saturated lapse [K/km]",                   gamma_moist * 1000.0, 4.84, 2e-2);
    check_true("the saturated lapse is well below the dry one",
               gamma_moist < 0.75 * gamma_dry);
    check_true("...and well above what the truncated form claims",
               gamma_moist > 4.0 * gamma_trunc);

    // Earth, as a control on the formula itself: the standard saturated adiabatic lapse
    // rate at 300 K / 1000 hPa is ~3.6-3.9 K/km. A formula that only ever runs at 60 bar
    // has nothing to be checked against; this is the one point where it does.
    const double gamma_earth = SaturationH2O::moistLapse(300.0, 1000.0, 0.02896,
                                                          1005.0, 287.0, 9.81);
    check("Earth control: saturated lapse at 300 K / 1 bar [K/km]",
          gamma_earth * 1000.0, 3.79, 3e-2);

    // And the dry limit. Above the critical point there is no saturation curve at all, so
    // moistLapse must return g/cp EXACTLY — this is the path every unsaturated layer of the
    // column takes, and a formula that only approximately reduces to the dry adiabat would
    // put a slow bias through all of it.
    check("dry limit: moistLapse above T_crit returns g/cp",
          SaturationH2O::moistLapse(700.0, 60000.0, M_nw_surf, cp_surf, S.R_mix, 9.81),
          9.81 / cp_surf, 1e-15);
    // At 200 K there IS a saturation curve (over ice) but the vapour pressure is ~1e-4 hPa
    // of a 60 bar column, so the correction is real and utterly negligible. Asserted as
    // "negligible" rather than "zero" because it is not zero.
    check("cold limit: the moist correction at 200 K is below 1e-5 relative",
          SaturationH2O::moistLapse(200.0, 60000.0, M_nw_surf, cp_surf, S.R_mix, 9.81),
          9.81 / cp_surf, 1e-5);

    // ------------------------------------------------------------------
    // 5. Known limits of the imported fits at this regime, asserted so they cannot be
    //    forgotten rather than silently relied on.
    // ------------------------------------------------------------------
    std::printf("\n limits of the imported fits\n");
    check_true("Shomate cp fits cover the sea surface",
               T_SURF > AtmMixture::T_FIT_MIN && T_SURF < AtmMixture::T_FIT_MAX);
    check_true("but NOT the cold trap: cp below 298 K is clamped, not extrapolated",
               AtmMixture::T_FIT_MIN > 250.0);

    std::printf("\n%s  (%d failure%s)\n\n",
                failures ? "self-test FAILED" : "self-test passed",
                failures, failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}

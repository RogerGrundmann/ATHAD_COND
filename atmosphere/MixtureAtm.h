#ifndef MIXTUREATM_H
#define MIXTUREATM_H

// ============================================================================
// MixtureAtm — thermodynamic properties of the Hadean atmospheric mixture
// ============================================================================
//
// The inherited Earth model could treat R and cp as constants because water vapour was
// a trace (~1 % by mass) riding on a fixed N2/O2 carrier. ATHAD cannot: H2O is 67 % of
// the mass and CO2 another 21 %, so the gas constant and heat capacity of a parcel are
// set by its own composition and vary as the water field evolves.
//
// This header provides the local mixture properties as functions of the two prognostic
// mass fractions (c = H2O, co2 = CO2), with everything else lumped into a fixed
// non-condensable background.
//
//   R_of  (c, co2)     specific gas constant       [J/(kg K)]
//   cp_of (c, co2, T)  specific heat, T-dependent  [J/(kg K)]
//   M_of  (c, co2)     mean molar mass             [kg/mol]
//
// MASS weighting, not mole weighting: for specific (per-kilogram) quantities the
// mixture value is the mass-fraction-weighted sum. Mole-weighting them is a real and
// easy mistake — ATNEPT c116d71 ("Weight the mixture properties by mass fraction") is
// the family's precedent for getting it wrong first.
//
// cp is strongly temperature-dependent over the 300–1500 K range ATHAD spans: H2O rises
// from ~1.86 to ~2.6 kJ/(kg K), CO2 from ~0.85 to ~1.33. A constant cp misplaces the
// lapse rate everywhere, so cp_i(T) uses Shomate-form fits (NIST/JANAF) per species.

#include <cmath>
#include <cstdlib>
#include <algorithm>

namespace AtmMixture {

    // ---- water critical point -----------------------------------------------
    // Above these, liquid and vapour are one phase: there is no saturation vapour
    // pressure, no condensation and no latent heat. ATHAD's surface sits at 1500 K and
    // p_H2O = 200 bar, i.e. deeply supercritical, so every condensation path must be a
    // genuine NO-OP there — not a clamp. See CLAUDE.md, invariant 2.
    constexpr double T_CRIT_H2O = 647.096;       // [K]
    constexpr double P_CRIT_H2O = 220640.0;      // [hPa] = 22.064 MPa

    // ---- universal constant and molar masses [kg/mol] ----------------------
    constexpr double R_STAR  = 8.314462618;      // [J/(mol K)]

    constexpr double M_H2O   = 0.018015;
    constexpr double M_CO2   = 0.044010;
    constexpr double M_N2    = 0.028014;
    constexpr double M_CH4   = 0.016043;
    constexpr double M_NH3   = 0.017031;
    constexpr double M_H2    = 0.002016;
    constexpr double M_CO    = 0.028010;
    constexpr double M_SO2   = 0.064066;

    // Specific gas constants [J/(kg K)] — R*/M.
    constexpr double R_H2O   = R_STAR / M_H2O;   // 461.5
    constexpr double R_CO2   = R_STAR / M_CO2;   // 188.9

    // ------------------------------------------------------------------------
    // Composition, resolved once from the configured mole fractions.
    //
    // Everything that is neither H2O nor CO2 is collapsed into one "background"
    // pseudo-species, because ATHAD transports only H2O and CO2; the rest stay well
    // mixed and can be represented by their aggregate molar mass and heat capacity.
    // ------------------------------------------------------------------------
    // Names of the six background gases, in the order f_bg[] uses.
    inline const char* const* BG_NAMES() {
        static const char* n[6] = {"N2", "CH4", "NH3", "H2", "CO", "SO2"};
        return n;
    }

    // All eight species in one order: the two prognostic ones (water vapour, CO2) followed by
    // the six background species in BG_NAMES order. The ParaView writers use this so the
    // plotted composition is ONE list from ONE source (split), rather than the three naming
    // conventions it grew. Ported from ATHAD (README item 61).
    inline const char* const* SPECIES_NAMES() {
        static const char* n[8] = {"H2O", "CO2", "N2", "CH4", "NH3", "H2", "CO", "SO2"};
        return n;
    }

    struct Composition {
        double M_mean   = 0.0;   // mean molar mass of the full mixture   [kg/mol]
        double R_mix    = 0.0;   // gas constant of the full mixture      [J/(kg K)]

        double q_H2O    = 0.0;   // mass fractions of the full mixture    [kg/kg]
        double q_CO2    = 0.0;
        double q_bg     = 0.0;

        double M_bg     = 0.0;   // background pseudo-species             [kg/mol]
        double R_bg     = 0.0;   //                                       [J/(kg K)]

        // The six background gases, as mass fractions OF THE BACKGROUND (they sum to 1).
        // Fixed: all six are well mixed and source-free, so their ratios never change and any
        // per-species quantity is this constant times the local q_bg. At this fork's shipped
        // composition only N2 is non-zero. Order is N2, CH4, NH3, H2, CO, SO2 — see BG_NAMES.
        double f_bg[6]  = {0,0,0,0,0,0};

        bool   valid    = false; // mole fractions summed to 1
        double x_sum    = 0.0;
    };

    // ------------------------------------------------------------------------
    // THE NON-WATER CARRIER AND ITS DILUTION (ported from ATHAD, its README items 57 and 59).
    //
    // WORSE HERE THAN THERE, because of what this fork is. ATHAD_COND is the post-condensation
    // atmosphere: CO2 is 54.78 % of the mass and water only 36.16 % (62.33 / 34.64 before the
    // trace gases were restored on 2026-08-20), and the water field spans 20.1 to 339.0 g/kg —
    // a 16.8x range against ATHAD's 1.4x. With the old split, which read co2 straight from the
    // field and gave the BACKGROUND every change in the water:
    //
    //     at the sea surface (q_v = 0.339)  q_CO2 = 0.5478 (pinned)   q_bg = 0.1132
    //     aloft              (q_v = 0.020)  q_CO2 = 0.5478 (pinned)   q_bg = 0.4322
    //
    // The background — 3 % of this atmosphere at the reference — was being handed 0.31 of the
    // mass aloft, a 9.5x swing, all of it stolen from CO2. R_mix aloft comes out 240.1 where
    // the carrier-consistent answer is 200.2: a 20 % error in the gas constant, hence in the
    // scale height. And kappa_CO2 = 1e-3 against kappa_bg = 1e-6, so the same misallocation
    // understates the CO2 optical depth aloft by ~50 %, at the levels the OLR comes from.
    //
    // The reference carrier fraction 1 - q_H2O_ref, set once by resolve(). 0 = unset, which
    // selects the legacy behaviour so nothing can depend on call order.
    // Set ONCE by cAtmosphereModel::initComposition(), from c_0 — the water mass fraction the
    // stored co2 value is quoted at, which is what makes q_CO2_of an identity at the reference.
    // 0 = unset, which selects the legacy no-dilution path rather than a wrong number.
    //
    // IT IS NOT SET INSIDE resolve(), AND THE COND SELF-TEST IS WHY. It was, briefly: resolve()
    // is the one function that knows the configured composition, so it looked like the natural
    // home. But resolve() is a pure function that any caller may invoke with any composition,
    // and cond_column_selftest.cpp calls it with a near-dry one — which left carrierRef() at
    // 0.9957 instead of 1 - c_0 and moved M_nonwater by a quarter. A global written
    // by whoever called last is not a reference. The test caught it on the first run.
    inline double& carrierRef() { static double v = 0.0; return v; }

    // ATM_CO2_DILUTE=0 restores the old split. Default ON.
    inline bool co2_dilute()
    {
        static const bool v = [](){
            const char* e = std::getenv("ATM_CO2_DILUTE"); return !(e && std::atoi(e) == 0); }();
        return v;
    }

    // The LOCAL CO2 mass fraction. The transported value is read as the CO2 mass fraction AT
    // THE REFERENCE WATER CONTENT and scaled to the local carrier:
    //
    //     q_CO2 = co2_stored * (1 - q_v) / (1 - q_H2O_ref)
    //
    // At the reference this is the identity, so co2_0, the initial condition and every printed
    // value keep their meaning; away from it q_CO2 and q_bg scale together and the carrier's
    // composition is invariant, which is the physical statement. CO2 stays prognostic.
    inline double q_CO2_of(double c, double co2)
    {
        const double q_c0 = std::min(std::max(co2, 0.0), 1.0);
        const double ref  = carrierRef();
        if (!co2_dilute() || ref <= 0.0) return q_c0;
        const double q_v  = std::min(std::max(c, 0.0), 1.0);
        return q_c0 * (1.0 - q_v) / ref;
    }

    // Build the composition from the eight configured mole fractions.
    inline Composition resolve(double x_H2O, double x_CO2, double x_N2,
                               double x_CH4, double x_NH3, double x_H2,
                               double x_CO,  double x_SO2)
    {
        Composition C;

        C.x_sum = x_H2O + x_CO2 + x_N2 + x_CH4 + x_NH3 + x_H2 + x_CO + x_SO2;
        C.valid = std::fabs(C.x_sum - 1.0) < 1.0e-9;

        // Mean molar mass: sum of x_i * M_i.
        const double m_H2O = x_H2O * M_H2O;
        const double m_CO2 = x_CO2 * M_CO2;
        const double m_bg  = x_N2  * M_N2  + x_CH4 * M_CH4 + x_NH3 * M_NH3
                           + x_H2  * M_H2  + x_CO  * M_CO  + x_SO2 * M_SO2;

        C.M_mean = m_H2O + m_CO2 + m_bg;
        C.R_mix  = R_STAR / C.M_mean;

        C.q_H2O = m_H2O / C.M_mean;
        C.q_CO2 = m_CO2 / C.M_mean;
        C.q_bg  = m_bg  / C.M_mean;

        const double x_bg = x_N2 + x_CH4 + x_NH3 + x_H2 + x_CO + x_SO2;
        C.M_bg = (x_bg > 0.0) ? (m_bg / x_bg) : M_N2;
        C.R_bg = R_STAR / C.M_bg;

        // Per-species split of the background, by MASS within the background.
        if (m_bg > 0.0) {
            C.f_bg[0] = x_N2  * M_N2  / m_bg;
            C.f_bg[1] = x_CH4 * M_CH4 / m_bg;
            C.f_bg[2] = x_NH3 * M_NH3 / m_bg;
            C.f_bg[3] = x_H2  * M_H2  / m_bg;
            C.f_bg[4] = x_CO  * M_CO  / m_bg;
            C.f_bg[5] = x_SO2 * M_SO2 / m_bg;
        }

        // NOTE: carrierRef() is deliberately NOT set here. See its declaration above — it is
        // set once by cAtmosphereModel::initComposition() from c_0, because resolve() is a
        // pure function that anything may call with any composition, and a global written by
        // "whoever called last" is not a reference.

        return C;
    }

    // ------------------------------------------------------------------------
    // Specific heat capacities [J/(kg K)] as a function of temperature.
    //
    // Shomate form  cp_molar = A + B*t + C*t^2 + D*t^3 + E/t^2   with t = T/1000 [K],
    // cp_molar in J/(mol K); divided by the molar mass to get the specific value.
    // Coefficients are the NIST/JANAF gas-phase fits valid over roughly 300–2000 K,
    // which is the range ATHAD spans. T is clamped to that range rather than
    // extrapolated — beyond it the polynomials diverge fast.
    // ------------------------------------------------------------------------
    constexpr double T_FIT_MIN = 298.0;
    constexpr double T_FIT_MAX = 2000.0;

    inline double shomate(double T, double A, double B, double Cc, double D, double E)
    {
        const double t  = std::min(std::max(T, T_FIT_MIN), T_FIT_MAX) * 1.0e-3;
        const double t2 = t * t;
        return A + B * t + Cc * t2 + D * t2 * t + E / t2;      // [J/(mol K)]
    }

    // H2O (gas), NIST 500–1700 K fit.
    inline double cp_H2O(double T) {
        return shomate(T, 30.09200, 6.832514, 6.793435, -2.534480, 0.082139) / M_H2O;
    }
    // CO2, NIST 298–1200 K fit.
    inline double cp_CO2(double T) {
        return shomate(T, 24.99735, 55.18696, -33.69137, 7.948387, -0.136638) / M_CO2;
    }
    // Background: dominated by N2 and CO (both 28 g/mol, near-identical cp), with CH4,
    // NH3, H2 and SO2 minor by mole. The N2 fit (100–500 K extended) is representative
    // to a few per cent, which is well inside the uncertainty of the composition itself.
    inline double cp_bg(double T, double M_background) {
        return shomate(T, 28.98641, 1.853978, -9.647459, 16.63537, 0.000117) / M_background;
    }

    // ------------------------------------------------------------------------
    // Local mixture properties from the two prognostic mass fractions.
    //
    // c and co2 are mass fractions in [0,1]; the background takes up the remainder.
    // Both are clamped and the background floored at zero, so a transport overshoot
    // degrades the properties smoothly instead of producing a negative gas constant.
    // ------------------------------------------------------------------------
    // q_cond is the SUSPENDED CONDENSATE mass fraction (cloud + ice + graupel), defaulting to 0
    // so every existing call compiles unchanged. The carrier is what is left of a kilogram of
    // parcel once the water is taken out, and "the water" is vapour AND condensate — a droplet
    // is still in the parcel and still weighs something.
    //
    // THE RETURNED FRACTIONS SUM TO 1 - q_cond, NOT TO 1, deliberately: they are per unit TOTAL
    // parcel mass, so R_of returns (1 - q_cond)*R_gas and p/(R_of*T) is the TOTAL density
    // directly. That identity is what lets densities()'s separate `water_factor` divisor be
    // deleted — this model already had the correction, in one place, applied only to r_humid,
    // spelled 1 - cloud - ice (no graupel) and floored at 0.5.
    inline void split(double c, double co2, double& q_v, double& q_c, double& q_b,
                      double q_cond = 0.0)
    {
        q_v = std::min(std::max(c,   0.0), 1.0);

        if (!co2_dilute() || carrierRef() <= 0.0) {      // legacy path, unchanged
            q_c = std::min(std::max(co2, 0.0), 1.0);
            const double sum = q_v + q_c;
            if (sum > 1.0) {                  // renormalise rather than go negative
                q_v /= sum;
                q_c /= sum;
            }
            q_b = std::max(0.0, 1.0 - q_v - q_c);
            return;
        }

        const double q_l  = std::min(std::max(q_cond, 0.0), 1.0);
        const double carr = std::max(0.0, 1.0 - q_v - q_l);          // the gas carrier
        q_c = std::min(std::max(co2, 0.0), 1.0) * carr / carrierRef();
        if (q_c > carr) q_c = carr;                                  // cannot exceed the carrier
        q_b = carr - q_c;
    }

    // Specific gas constant of the local mixture [J/(kg K)].
    inline double R_of(double c, double co2, double R_background, double q_cond = 0.0)
    {
        double q_v, q_c, q_b;
        split(c, co2, q_v, q_c, q_b, q_cond);
        return q_v * R_H2O + q_c * R_CO2 + q_b * R_background;
    }

    // Specific heat at constant pressure of the local mixture [J/(kg K)].
    inline double cp_of(double c, double co2, double T, double M_background, double q_cond = 0.0)
    {
        double q_v, q_c, q_b;
        split(c, co2, q_v, q_c, q_b, q_cond);
        return q_v * cp_H2O(T) + q_c * cp_CO2(T) + q_b * cp_bg(T, M_background);
    }

    // Mean molar mass of the local mixture [kg/mol]. 1/M = sum(q_i / M_i).
    inline double M_of(double c, double co2, double M_background)
    {
        double q_v, q_c, q_b;
        split(c, co2, q_v, q_c, q_b);
        const double inv = q_v / M_H2O + q_c / M_CO2 + q_b / M_background;
        return (inv > 0.0) ? 1.0 / inv : M_background;
    }

    // Mole fraction of CO2 given the local mass fractions — what an optical-depth or
    // partial-pressure calculation needs. x_i = q_i * M_mix / M_i.
    inline double x_CO2_of(double c, double co2, double M_background)
    {
        double q_v, q_c, q_b;
        split(c, co2, q_v, q_c, q_b);
        return q_c * M_of(c, co2, M_background) / M_CO2;
    }

    // Mean molar mass of everything EXCEPT water [kg/mol] — what the exact saturation
    // mass-fraction conversion needs as its "other" carrier. Built from the local CO2 and
    // background mass fractions, renormalised to exclude H2O.
    //
    // It needs the water mass fraction to do that. The first version took only (co2,
    // M_background) and set q_b = 1 - q_c, so the renormalisation by (q_c + q_b) was
    // identically 1 and the water's share of the mass was silently handed to the
    // background gas. At the reference composition that returned 28.58 g/mol instead of
    // 35.11 — the carrier reported 19 % too light, hence a q_sat some 23 % too large
    // wherever the conversion is dilute. Harmless on Earth, where water is ~1 % of the
    // mass and the misallocation is invisible; here water is 67 %.
    inline double M_nonwater(double c, double co2, double M_background)
    {
        double q_v, q_c, q_b;
        split(c, co2, q_v, q_c, q_b);
        const double sum = q_c + q_b;
        if (!(sum > 0.0)) return M_background;
        const double inv = (q_c / sum) / M_CO2 + (q_b / sum) / M_background;
        return (inv > 0.0) ? 1.0 / inv : M_background;
    }

    // Mole fraction of H2O, likewise.
    inline double x_H2O_of(double c, double co2, double M_background)
    {
        double q_v, q_c, q_b;
        split(c, co2, q_v, q_c, q_b);
        return q_v * M_of(c, co2, M_background) / M_H2O;
    }

}   // namespace AtmMixture

#endif

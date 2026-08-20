#ifndef SATURATIONH2O_H
#define SATURATIONH2O_H

// ============================================================================
// SaturationH2O — water saturation physics valid to the critical point
// ============================================================================
//
// The inherited Earth model used the Magnus formula,
//     E_sat = hp * exp(a*(T - T_0)/(T - b))
// with (a, b) = (17.2694, 35.86) over water and (21.8746, 7.66) over ice. That is an
// empirical fit calibrated on roughly 200–320 K, and the code that used it had to cap T
// at 60 °C to stay inside its range. ATHAD's surface is 1500 K and its column spans
// 1500 K down to ~200 K, so Magnus is unusable across most of the domain — evaluated at
// 1500 K it returns ~1.2e7 hPa, fifty times the total pressure of the atmosphere, which
// flips the sign of every q_sat denominator that divides by (p - E).
//
// This header replaces it with the IAPWS formulations, which are the reference
// correlations for water and are valid over exactly the range ATHAD needs:
//
//   saturationPressure(T)   Wagner & Pruss (IAPWS-95) saturation line, 273.16 K to the
//                           critical point 647.096 K / 22.064 MPa. Returns NO_SATURATION
//                           above T_crit, where liquid and vapour are one phase and the
//                           quantity does not exist.
//   sublimationPressure(T)  IAPWS (2011) ice-Ih sublimation curve, below 273.16 K.
//   latentHeat(T)           Watson correlation, anchored at 2.501e6 J/kg at 273.15 K and
//                           correctly vanishing at the critical point.
//   saturationMassFraction  the EXACT mass-fraction conversion — see below.
//
// Why the exact conversion matters. The Earth code computed
//     q_sat = ep * E / (p - (1 - ep) * E)
// which is the dilute approximation: it assumes the condensable is a trace riding on a
// carrier of fixed molar mass. Here water IS the bulk gas (67 % by mass, 80 % by mole),
// so the mean molar mass of the mixture changes as water condenses out and the
// approximation has no small parameter left. The exact conversion from a mole fraction
// x = E/p to a mass fraction is used instead:
//
//     q = x*M_H2O / ( x*M_H2O + (1 - x)*M_other )
//
// which reduces to the dilute form when x << 1 and stays correct when it does not.

#include "MixtureAtm.h"

#include <cmath>
#include <algorithm>

namespace SaturationH2O {

    // Sentinel returned when no saturation state exists (T at or above the critical
    // temperature). Callers must test for it rather than using the value.
    constexpr double NO_SATURATION = -1.0;

    constexpr double T_TRIPLE   = 273.16;        // [K]
    constexpr double P_TRIPLE   = 6.11657;       // [hPa]  611.657 Pa
    constexpr double P_CRIT_hPa = 220640.0;      // [hPa]  22.064 MPa

    inline bool isSupercritical(double T) { return T >= AtmMixture::T_CRIT_H2O; }

    // ------------------------------------------------------------------------
    // Saturation vapour pressure over liquid water [hPa].
    //
    // Wagner & Pruss (IAPWS-95) saturation-line equation:
    //     ln(p/p_c) = (T_c/T) * sum_i a_i * theta^n_i ,   theta = 1 - T/T_c
    //
    // Valid 273.16 K -> 647.096 K. Reference points (checked in test/saturation_selftest):
    //     273.16 K -> 6.11657 hPa      373.15 K -> 1013.5 hPa (approx 1 atm)
    //     647.096 K -> 220640 hPa      (the critical point)
    // ------------------------------------------------------------------------
    inline double saturationPressure(double T)
    {
        if (isSupercritical(T)) return NO_SATURATION;
        if (T <= 0.0)           return 0.0;

        const double Tc    = AtmMixture::T_CRIT_H2O;
        const double theta = 1.0 - T / Tc;

        static constexpr double a1 = -7.85951783;
        static constexpr double a2 =  1.84408259;
        static constexpr double a3 = -11.7866497;
        static constexpr double a4 =  22.6807411;
        static constexpr double a5 = -15.9618719;
        static constexpr double a6 =  1.80122502;

        const double t15 = theta * std::sqrt(theta);          // theta^1.5
        const double t3  = theta * theta * theta;
        const double t35 = t3 * std::sqrt(theta);             // theta^3.5
        const double t4  = t3 * theta;
        const double t75 = t4 * t3 * std::sqrt(theta);        // theta^7.5

        const double sum = a1 * theta + a2 * t15 + a3 * t3
                         + a4 * t35   + a5 * t4  + a6 * t75;

        return P_CRIT_hPa * std::exp((Tc / T) * sum);
    }

    // ------------------------------------------------------------------------
    // Sublimation pressure over ice Ih [hPa], IAPWS (2011).
    //     ln(p/p_t) = [ sum_i a_i * theta^b_i ] / theta ,   theta = T/T_t
    // Valid roughly 50 K -> 273.16 K.
    //
    // NOTE the summand is a_i * theta^b_i, NOT a_i * (1 - theta^b_i). Writing the
    // (1 - ...) form gives a result that is the exact negative of the right answer in the
    // exponent, i.e. p_subl(250 K) = 49 hPa instead of 0.76 hPa — a factor of 65 the wrong
    // way, and increasing as T falls. The self-test caught it.
    // ------------------------------------------------------------------------
    inline double sublimationPressure(double T)
    {
        if (T <= 0.0) return 0.0;

        const double theta = T / T_TRIPLE;

        static constexpr double a1 = -0.212144006e2, b1 = 0.333333333e-2;
        static constexpr double a2 =  0.273203819e2, b2 = 0.120666667e1;
        static constexpr double a3 = -0.610598130e1, b3 = 0.170333333e1;

        const double sum = a1 * std::pow(theta, b1)
                         + a2 * std::pow(theta, b2)
                         + a3 * std::pow(theta, b3);

        return P_TRIPLE * std::exp(sum / theta);
    }

    // Saturation pressure over whichever phase is stable at T [hPa], or NO_SATURATION.
    inline double saturationPressureAuto(double T)
    {
        if (isSupercritical(T)) return NO_SATURATION;
        return (T >= T_TRIPLE) ? saturationPressure(T) : sublimationPressure(T);
    }

    // ------------------------------------------------------------------------
    // Latent heat of vaporisation [J/kg], Watson correlation:
    //     L(T) = L_ref * ((T_c - T)/(T_c - T_ref))^0.38
    //
    // The inherited model used a CONSTANT lv = 2.52e6 J/kg. Latent heat is not constant
    // and, more importantly, it goes to ZERO at the critical point — there is no phase
    // change left to release it. A constant lv near T_crit injects heat that the physics
    // does not have, which is exactly the runaway the old Magnus cap was there to stop.
    // ------------------------------------------------------------------------
    inline double latentHeat(double T)
    {
        const double Tc = AtmMixture::T_CRIT_H2O;
        if (T >= Tc) return 0.0;

        constexpr double L_ref = 2.501e6;    // [J/kg] at 273.15 K
        constexpr double T_ref = 273.15;     // [K]

        return L_ref * std::pow((Tc - T) / (Tc - T_ref), 0.38);
    }

    // Latent heat of sublimation [J/kg]: vaporisation plus fusion (3.337e5 J/kg).
    inline double latentHeatSublimation(double T)
    {
        constexpr double L_fusion = 3.337e5;                 // [J/kg]
        return latentHeat(T) + L_fusion;
    }

    // ------------------------------------------------------------------------
    // EXACT saturation mass fraction.
    //
    //   E        saturation partial pressure of H2O   [hPa]
    //   p        total pressure                       [hPa]
    //   M_other  mean molar mass of everything else   [kg/mol]
    //
    // Returns kg H2O per kg of mixture. When E >= p the whole column would be vapour, so
    // the mass fraction saturates at 1 rather than dividing by a negative denominator —
    // which is precisely what the dilute form did, and where its runaway condensation
    // came from.
    // ------------------------------------------------------------------------
    inline double saturationMassFraction(double E, double p, double M_other)
    {
        if (E == NO_SATURATION) return 1.0;      // supercritical: no condensation limit
        if (!(p > 0.0))         return 0.0;
        if (E >= p)             return 1.0;

        const double x = E / p;                              // mole fraction at saturation
        const double num = x * AtmMixture::M_H2O;
        const double den = num + (1.0 - x) * M_other;
        return (den > 0.0) ? num / den : 1.0;
    }

    // ------------------------------------------------------------------------
    // EXACT INVERSE of saturationMassFraction: vapour partial pressure [same units as p]
    // from a water MASS fraction q and the molar mass of everything that is not water.
    //
    // This is the partner the file was missing, and its absence is why
    // ThermoAtm::waterVapourEvaporation() computed q -> e with the dilute q*p/ep while
    // computing e -> q exactly two lines earlier. Inverting
    //     q = x*M_H2O / ( x*M_H2O + (1-x)*M_other )
    // for the mole fraction x gives
    //     x = q*M_other / ( M_H2O*(1-q) + q*M_other )
    // and e = x*p. It reduces to the dilute e = q*p/ep when q << 1 and stays correct when
    // it does not — at q = 0.3616 (ATHAD_COND's sea surface) the dilute form is low by
    // 1150 hPa out of 33470, which reads as a saturation deficit at a surface that is
    // saturated by construction. See README, "What the surface state still breaks".
    inline double vapourPressureFromMassFraction(double q, double p, double M_other)
    {
        if (!(p > 0.0))   return 0.0;
        if (!(q > 0.0))   return 0.0;
        if (q >= 1.0)     return p;                          // pure vapour
        const double den = AtmMixture::M_H2O * (1.0 - q) + q * M_other;
        if (!(den > 0.0)) return p;
        const double x = q * M_other / den;                  // mole fraction
        return x * p;
    }

    // ------------------------------------------------------------------------
    // Dew-point temperature [K] from a vapour partial pressure e [hPa].
    //
    // Magnus can be inverted in closed form; IAPWS cannot, so this bisects the
    // (strictly monotone) saturation curve instead. 60 iterations over 150-647 K give
    // better than 1e-15 K of bracket width, and the cost is irrelevant — this is a
    // diagnostic, evaluated once per cell per output step.
    //
    // Returns T_CRIT if e exceeds the critical pressure: there is no dew point above it.
    inline double dewPoint(double e)
    {
        if (!(e > 0.0))       return 0.0;
        if (e >= P_CRIT_hPa)  return AtmMixture::T_CRIT_H2O;

        double lo = 150.0, hi = AtmMixture::T_CRIT_H2O;
        if (saturationPressureAuto(lo) >= e) return lo;      // drier than the table reaches

        for (int it = 0; it < 60; it++) {
            const double mid = 0.5 * (lo + hi);
            if (saturationPressureAuto(mid) < e) lo = mid; else hi = mid;
        }
        return 0.5 * (lo + hi);
    }

    // Latent heat of whichever phase change is stable at T [J/kg] — the partner of
    // saturationPressureAuto, so a derivative taken along the Auto curve uses the L that
    // belongs to it. Getting these two out of step is a 20 % error in the ice range.
    inline double latentHeatAuto(double T)
    {
        return (T >= T_TRIPLE) ? latentHeat(T) : latentHeatSublimation(T);
    }

    // ------------------------------------------------------------------------
    // The saturation-derivative factor A [dimensionless], shared by dq_sat/dT and the
    // moist lapse rate so the two cannot disagree.
    //
    //     q_sat = x*Mw / (x*Mw + (1-x)*Mo),        x = E(T)/p
    //     dq/dx = Mw*Mo / den^2,                   den = x*Mw + (1-x)*Mo
    //     A     = x * Mw*Mo / den^2
    //
    // so that  (dq/dT)_p = A * L/(R_v T^2)  and  (dq/dp)_T = -A/p.
    //
    // WHY THIS IS NOT THE FAMILIAR q_sat*L/(R_v T^2): that form is the x -> 0 limit. A
    // reduces to q_sat exactly when x is small (den -> Mo, A -> x*Mw/Mo = q_sat), which is
    // why Earth never noticed. At the ATHAD_COND sea surface x = 0.558 and A/q_sat =
    // Mo/den = 1.48, so the dilute derivative understates the real one by a third. Water is
    // 35 % of the mass there; it is not a trace and cannot be differentiated like one.
    // ------------------------------------------------------------------------
    // Phase given explicitly, for callers that need the liquid and the ice branch at the
    // same temperature (the saturation adjustment weighs both).
    inline double satDerivFactorAt(double E, double p, double M_other)
    {
        if (E == NO_SATURATION || !(p > 0.0) || E >= p) return 0.0;

        const double x   = E / p;
        const double den = x * AtmMixture::M_H2O + (1.0 - x) * M_other;
        if (!(den > 0.0)) return 0.0;

        return x * AtmMixture::M_H2O * M_other / (den * den);
    }

    // Phase chosen by temperature, the common case.
    inline double satDerivFactor(double T, double p, double M_other)
    {
        return satDerivFactorAt(saturationPressureAuto(T), p, M_other);
    }

    // ------------------------------------------------------------------------
    // dq_sat/dT [1/K] at constant pressure, exact for a mixture of any water content.
    // Follows the IAPWS curve because L(T) does; the Magnus analytic derivative it
    // replaces carried the fit's own coefficients and is wrong wherever the fit is.
    inline double dqSatdTFrom(double E, double L, double T, double p, double M_other)
    {
        if (!(T > 0.0)) return 0.0;
        return satDerivFactorAt(E, p, M_other) * L / (AtmMixture::R_H2O * T * T);
    }

    inline double dqSatdT(double T, double p, double M_other)
    {
        return dqSatdTFrom(saturationPressureAuto(T), latentHeatAuto(T), T, p, M_other);
    }

    // ------------------------------------------------------------------------
    // SATURATED (moist) adiabatic lapse rate [K/m], positive downward-decreasing.
    //
    // From conservation of moist static energy along a saturated ascent,
    //     cp dT + L dq + g dz = 0,   q = q_sat(T, p(z)),   dp = -rho g dz,
    // and with (dq/dT)_p = A L/(R_v T^2) and (dq/dp)_T = -A/p, rho = p/(R_mix T):
    //
    //     dT/dz = -g * [1 + L*A/(R_mix*T)] / [cp + L^2*A/(R_v*T^2)]
    //
    // BOTH terms matter and the numerator is the one that gets dropped. Keeping only the
    // denominator — the familiar g/(cp + L dq/dT) — gives 0.7 K/km at the ATHAD_COND sea
    // surface; the full expression gives 5.04 K/km against a dry 7.27. A factor of seven
    // apart, and the truncated form is the one that looks like the textbook.
    //
    // Checks out against Earth: at 300 K / 1000 hPa with M_other = 0.02896 it returns
    // 3.79 K/km, the standard saturated adiabatic value. Returns the DRY lapse g/cp wherever
    // no condensation is possible, so a caller can use it unconditionally.
    // ------------------------------------------------------------------------
    inline double moistLapse(double T, double p, double M_other,
                             double cp, double R_mix, double g)
    {
        if (!(cp > 0.0) || !(T > 0.0) || !(R_mix > 0.0)) return 0.0;

        const double A = satDerivFactor(T, p, M_other);
        if (!(A > 0.0)) return g / cp;                       // nothing to condense

        const double L   = latentHeatAuto(T);
        const double num = 1.0 + L * A / (R_mix * T);
        const double den = cp  + L * L * A / (AtmMixture::R_H2O * T * T);

        return (den > 0.0) ? g * num / den : g / cp;
    }

    // Convenience: saturation mass fraction directly from T, p and the background molar
    // mass, using whichever phase is stable. Returns 1.0 where supercritical (i.e. no
    // constraint), so callers that compare q_v > q_sat naturally find nothing to condense.
    inline double saturationMassFractionAt(double T, double p, double M_other)
    {
        return saturationMassFraction(saturationPressureAuto(T), p, M_other);
    }

}   // namespace SaturationH2O

#endif

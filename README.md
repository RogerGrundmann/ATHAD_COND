# ATHAD_COND

**The Hadean atmosphere after condensation.**

An atmospheric general-circulation model of the epoch that follows
[ATHAD](https://github.com/RogerGrundmann/ATHAD): the magma ocean has quenched, the steam
atmosphere has rained out into a liquid ocean, and what remains above it is a
CO₂-dominated atmosphere of **27–100 bar over a 230–250 °C sea**. Same solver — a
finite-difference Navier–Stokes solver on a spherical shell with RK4 time integration and
vertical coordinate stretching — applied to the opposite thermodynamic regime.

ATHAD_COND is forked from ATHAD at `29ca2f9`, **carrying its full history**, so fixes
cherry-pick in both directions. ATHAD in turn is the atmosphere half of
[ATOM_Precipitation](https://github.com/RogerGrundmann/ATOM_Precipitation) (forked at
`1e3f319`).

**Where it differs from ATHAD, in one sentence:** ATHAD's water is supercritical from the
ground to ~177 km and *nothing condenses*; here water is subcritical everywhere and
condensation is live from the sea surface up, so the code paths ATHAD spent seventeen
defect-fixes making inert are the ones this model depends on.

## Repository layout

```
atmosphere/   C++ source — the atmosphere model
cli/          Command-line interface source and config file
lib/          Array classes, config reader, shared utilities
python/       Cython bindings, run scripts, output directory
test/         Standalone self-tests
tinyxml2/     Vendored XML parser
```

## Build, run, test

```bash
make cond                                       # -> cli/cond
make test                                       # self-tests, run these first
cd python && OMP_NUM_THREADS=8 ../cli/cond config_cond.xml
```

Output lands in `python/output_Hadean_condensation/`.

`make` regenerates the parameter bindings from `param.py` first. The generated files are
tracked on purpose, so a `param.py` change shows its full effect in the diff.

## Conditions

The stated inputs are **T = 230–250 °C, p = 27–100 bar, H₂O 0.4–2 %, CO₂ 89–95 %,
N₂ 5–20 % by mole**. The nominal run point is the midpoint, **60 bar / 513.15 K**.

**The stated composition and the stated temperature describe different heights, and this
is the single most important thing to understand about the configuration.** The quoted
0.4–2 % H₂O is the *dry* atmosphere, what survives above the cold trap. Air in contact
with a 240 °C sea holds water at its own vapour pressure: `p_sat(513.15 K) = 33.47 bar` of
a 60 bar column, so `x_H₂O = 0.5578` at the surface — fifty-five times the dry figure. Both
are true, at different heights, and the column runs between them.

| | at the sea | above the cold trap |
|---|---|---|
| x = (H₂O, CO₂, N₂) | .5578 / .4109 / .0313 | .010 / .920 / .070 |
| M_mean [g/mol] | 29.01 | 42.63 |
| R_mix [J/(kg·K)] | 286.6 | 195.0 |
| cp [J/(kg·K)] | 1349 | 1028 |
| q_H₂O [kg/kg] | 0.3464 | 0.00423 |
| ρ [kg/m³] | 40.80 | — |
| scale height | 15.0 km | 10.2 km |

Water boils at 548.7 K at 60 bar, so the sea is 35 K below boiling. (At the 27 bar end of
the stated range water boils at 501.2 K, so a 503 K surface there has **no liquid ocean at
all** — that corner is not a colder version of this one, it is a different model.)

`test/cond_column_selftest.cpp` asserts every number in that table, plus the arithmetic
connecting the config to it. The config's `x_*` and `r_air` are derived quantities; change
one without the other and the test fails, which is the point.

## Four invariants — do not silently break these

1. **There is no topography, and the planet is hemispherically symmetric.** `h ≡ 0`,
   `i_topography ≡ 0` everywhere, so `AtomUtils::is_land()` is false at every point — now a
   global *ocean* rather than a global magma surface, but the same statement in code.
   `LandOceanFraction()` throws if a land point ever appears. Any north–south asymmetry in
   the output is a defect: the surface temperature is a symmetric parabola, the insolation
   is mirrored, and there is no obliquity and no seasonal cycle.

2. **Water is subcritical everywhere, and condensation is live.** *This is ATHAD's
   invariant 2, inverted.* Critical point 647.096 K / 220.64 bar; this column runs
   511 K / 60 bar downward, so every condensation path is a path that runs. What carries
   over unchanged is *which formulas are allowed*: **never reintroduce Magnus** (calibrated
   to ~320 K, returns 1.2e7 hPa at 513 K) and **never the dilute `q_sat = ep·E/(p−E)`**.
   The dilute form is invalid here for the opposite reason it was invalid in ATHAD — water
   is 0.4 % by mole above the cold trap, where it would be fine, and 56 % at the sea, where
   it is not. A form valid over part of the column is worse than one valid over none,
   because it looks right in the printouts.

3. **Radiation runs in mode 2** (direct σT⁴). Modes 0/1/3/4/5 lean on the Scotese snapshot
   or the 280 ppm CO₂ reference; neither exists in a 92 % CO₂ atmosphere at 4.4 Ga.

4. **The column is on its own adiabat, integrated not fitted — and the adiabat is MOIST
   where the air is saturated.** `SaturationH2O::moistLapse`, hydrostatic on the layer-mean
   T, isothermal above `t_skin`. Do **not** restore the COSMO `T = T₀√(1−coeff·h)` form, and
   do **not** revert to the dry adiabat: `cosmo_lapse_fraction = 1.0` is justified in ATHAD
   *because nothing condenses there*, and that justification is void here.

## What the model currently says

Everything below is measured. Items refer to the commits.

**1. The moist adiabat, and the number that nearly went in wrong.** The saturated lapse
rate is

    dT/dz = −g·[1 + L·A/(R_mix·T)] / [cp + L²·A/(R_v·T²)],   A = x·Mw·Mo/den²

and the familiar `g/(cp + L·dq_sat/dT)` is only its denominator. That truncation gives
**0.7 K/km** at this sea surface; the full expression gives **5.04 K/km** against a dry
7.27. A factor of seven, and the truncated form is the one that looks like the textbook.
The neglected term is the pressure dependence of `q_sat`, negligible for a trace gas and
not for 56 % water. Checked against Earth at 300 K / 1 bar, where it returns 3.81 K/km.

**2. `dq_sat/dT` was the dilute limit.** `q_sat·L/(R_v·T²)` is the `x → 0` form; the exact
mass-fraction derivative is larger by `M_other/(x·Mw + (1−x)·Mo)` = **1.478** at the sea
surface. This fed the saturation adjustment's Newton damping `ω = 1/(1+Gain)`, so the
damping was a third too weak in exactly the high-`q_sat` regime it exists for.

**3. The initial state is a fixed point of three fields, and iterating it does not
converge.** `c` depends on T and p, `co2` depends on `c`, and T and p depend on both.
Iterating them closes a loop over the runaway-greenhouse feedback: 3 passes gave 446.7 K at
8.9 km, 8 passes gave 463.8 K, still climbing. That is the correct physics arriving through
the wrong numerics. `densities(true)` integrates the water profile in the same upward sweep
as T and p — one pass, exact.

**4. The measured column.** A saturated troposphere **62 km deep**: 511 K at the sea, 487 K
at 4.5 km, 429 K at 18 km, 367 K at 34 km, 280 K at 60 km, with the first 10 km at
**5.02 K/km** against the 5.04 predicted. `q_H₂O` tracks `q_sat` to four digits throughout.
Precipitable water **133.5 m**.

**5. The cold trap is eight times too wet.** It sits at the inherited `t_skin` and 0.02 bar
and passes **0.0349 kg/kg** of water, against the **0.00423** the stated dry composition
implies. The stated "H₂O 0.4–2 % by mole" and the stated 230–250 °C sea are mutually
consistent only if the cold trap is near **217 K**, which is what `p_sat/p = 0.0099` at
0.02 bar requires. So `c_h2o_dry_top` is currently a floor the column never reaches.

**6. The OLR is not yet an output.** At 20 iterations:

```
mean planetary albedo    = 0.4997        absorbed SW + geothermal = 120.84 + 150.00 = 270.84 W/m2
outgoing long wave (OLR) = 268.10 W/m2   imbalance                = 2.74 W/m2
lid temperature / eps    = 262.22 K / 0.0000
sigma*T_lid^4            = 268.10 W/m2
```

**OLR and σT_lid⁴ are equal to six figures.** This is exactly the condition ATHAD's item 10
diagnosed and the diagnostic was built to flag. The cause here is different from ATHAD's,
though: the lid is genuinely transparent (eps = 0.0000) and the shell is deep enough. The
grey optical depth reaches 1 at ~0.07 bar, which is *where the isothermal skin begins*, so
the emission to space comes from a layer whose temperature is prescribed — and `t_skin` is
a fixed point relaxing to close the very balance being reported. The balance closes to
2.74 W/m² because it was constructed to.

Deepening the shell will not fix this. It is the isothermal-skin prescription that has to
go, which is the same open task as ATHAD's invariant 3.

## Assumptions vs. results — read this before quoting any number

| Parameter | Value | Status |
|---|---|---|
| `geothermal_flux` | 150 W/m² | **Inherited from ATHAD's magma ocean and almost certainly wrong here.** It is 55 % of the entire energy input, against modern Earth's 0.087 W/m². A quenched surface under a liquid ocean does not deliver this. **Biggest unexamined input.** |
| `kappa_H2O` / `kappa_CO2` / `kappa_bg` | 0.01 / 0.001 / 1e-6 m²/kg | Factor-of-2 uncertain. `kappa_CO2` now dominates: the CO₂ column is 5.8e5 kg/m² against water's 2.6e3 |
| `t_surf_equator` / `t_surf_pole` | 513.15 / 503.15 K | **Prescribed, not solved.** The 10 K contrast is an assumption |
| `t_skin` | fixed-point iterate | Currently *sets* the OLR — see item 6 |
| `albedo_cloud` | 0.50 | IS the model's planetary albedo (0.4997 measured); saturates the moment condensate exists, and condensate is now everywhere |
| insolation | 0.71 S₀ | Faint young Sun |
| `omega` | 3.17e-4 (5.5 h day) | Inherited from ATHAD; a post-condensation Earth is later and slower |
| `delta_i_c` | 500 s over ocean | Bechtold's convective adjustment timescale, Earth-calibrated. A time, not a pressure — no unit error to correct, and no evidence here to replace it |
| Shomate `cp` below 298 K | clamped | The cold trap and skin sit below the fit range |

A grey scheme cannot represent the window regions that set the real runaway limit.

## The pattern this fork was expected to bring, and did

ATHAD's README records seventeen defects of one shape: *Earth's numbers as bare literals
inside physics kernels*. At 250 bar those are wrong by orders of magnitude and obvious. **At
60 bar they are wrong by a factor, and plausible** — which is harder to see. Two examples,
both found here:

- `MoistConvection`'s triggers (1000/970/900/800 hPa). At ATHAD's 250 bar the whole column
  sits above 1000 hPa, every trigger fails, and deep convection is inactive — wrong, but
  loudly wrong. At 60 bar `p_stat <= 1000 hPa` does **not** fail: it selects everything
  above ~53 km, so the scheme fires at the top of the saturated column instead of at the
  sea, producing a plausible-looking deep convection keyed to the stratosphere. Now
  fractions of the local surface pressure.
- `initCloudIce`'s critical-humidity parabola, roots at `p = 0` and `p_crit = 1000 hPa`. At
  60 bar `x = p/p_crit` exceeds 1 through the whole troposphere and `x·(1−x)` goes
  **negative** — the critical humidity stops reducing the cloud threshold and starts
  amplifying it, precisely where the cloud is. Now anchored to `p_0`.

And one of the same shape in a field rather than a constant: ATHAD holds the CO₂ **mass**
fraction uniform, which is right there because its water is uniform too. Here water runs
0.346 → 0.004, so a uniform mass fraction makes the background absorb all the mass the
water vacates, and R settles at 230 J/(kg·K) aloft instead of 195 — an 18 % error through
the whole upper atmosphere, from the one field that is supposed to have no structure. What
is well mixed is the CO₂:background ratio *within the dry air*.

## Remaining work

- **The OLR is not independent of `t_skin`** (item 6). This is the first thing to fix and
  it is the same task as ATHAD's invariant 3: the temperature profile is prescribed, so
  radiation cannot set it.
- **`geothermal_flux = 150 W/m²` is inherited and probably an order of magnitude or more
  too large.** It is over half the energy budget. Nothing quantitative here survives it.
- **The cold-trap water content is 8× the stated composition** (item 5). Testable
  prediction: a 217 K cold trap reconciles them.
- **Rain has nowhere to go.** The sea is a boundary condition, not a reservoir, so once
  the moist physics runs (from iteration 0 — `moist_phys_start_iter` is 0 here, not ATHAD’s 300) the water budget stops being a closed-system test.
  It needs an evaporation − precipitation − storage closure.
- **The initial water profile is built on the provisional temperature** that
  `initTemperatureData` leaves, because that routine still runs before the composition
  exists. `densities(true)` repairs it, but it should be built on the real mixture.
- **Not yet run beyond 20 iterations.** No stability run, no grid-convergence check at
  100/120/140 km, no thread-determinism check. ATHAD's experience says 400 iterations is a
  stability check and not a convergence check — its meridional wind was still in free
  acceleration at 400, with the pressure gradient at 1.8 % of the Coriolis term.
- **There is no balanced initial state here, and the one next door must not be ported as
  first written.** ATHAD found that its prescribed circulation is buried within ~5 iterations
  by an unopposed Coriolis torque and built `initBalancedState` against it (its items 26-28);
  this model has neither the fix nor a measurement of whether it has the problem. **The
  caution is specific**: item 27's version balances the θ-momentum equation alone, and with
  the switch defaults this model shares — `coriolis_nontraditional()` and
  `metric_curvature()` both false, `buoyancy_ramp` = 0 at iteration 0 — the radial equation
  is `rhs_u = −dp_dyn/dr·exp_rm` and nothing else, so its radial gradient is an unopposed
  vertical force. In ATHAD it drove `max |u|` from 0.114 to 11.17 m/s over 200 iterations
  with the tropics sinking, while the meridional streamfunction, built from `v` alone,
  reported the cell as healthy. Port item 28's two-component version with its residual
  diagnostic, or nothing. See CLAUDE.md, *Relationship to the family*.
- **A term written in `RHS_Atm_Turb.cpp` is not necessarily a term the model applies.** The
  above is the first instance found of this shape — the Earth-constant pattern one level up,
  a code path gated off by a default switch — and it is worth expecting again, because a
  balance, a diagnostic or a budget derived from the source as written will silently disagree
  with the model as configured.
- **Inherited unchanged from ATHAD**: the column air mass is not conserved (`p_stat.x[0]`
  re-anchored every iteration); the `c ≤ 1 − co2` ceiling deletes water; `ATM_ANELASTIC`
  ships default-off and has not been measured in this regime, where the density span is
  ~3 orders of magnitude instead of ~5 — which makes this the better testbed for whether it
  matters at all.

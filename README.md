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

## The background opacity, split per species — a no-op that stopped being one

Ported from ATHAD (its README item 60), where `kappa_bg` = 1e-6 m2/kg — the value of a
radiatively inert diatomic — was being applied to a background containing NH3 and CH4 at
1.4 mole-% each, and splitting it raised the effective background opacity **1867x**.

**AT THIS FORK'S COMPOSITION IT CHANGED NOTHING — UNTIL 2026-08-20, WHEN THE COMPOSITION
CHANGED.** The background here was pure N2, `x_CH4` through `x_SO2` all zero, so
`sum(f_i*kappa_i)` = `kappa_N2` = 1e-6 and the lumped value was exact: OLR at iteration 20 was
**275.29 W/m2 in both arms**, identical, as a no-op must be. Then the five trace gases were
given ATHAD's mole fractions (see *Remaining work*), and the same code produced:

| species | fraction of bg | kappa [m2/kg] | share of kappa_bg_eff |
|---|---|---|---|
| N2 | 0.2927 | 1.000e-06 | 0.0 % |
| CH4 | 0.0892 | 3.000e-03 | 13.8 % |
| **NH3** | 0.0947 | 1.000e-02 | **48.7 %** |
| H2 | 0.0112 | 1.000e-05 | 0.0 % |
| CO | 0.1558 | 1.000e-04 | 0.8 % |
| **SO2** | 0.3563 | 2.000e-03 | **36.7 %** |
| **TOTAL** | 1.0000 | **1.944e-03** | **1944x the lumped value** |

**That is the whole case for the split, made twice over.** It was right by luck rather than by
construction while the background was nitrogen, and nothing in the code said which; the same
constant was exact here and 1867x too small in the sibling. Under the lumped `kappa_bg` the
composition change would have run five new absorbers at nitrogen's opacity and printed nothing
to say so. `ATM_BG_LUMPED=1` restores the single value and every background opacity that
predates the composition change.

`q_N2 … q_SO2` are written to all three VTK slices here as well, computed on the fly from
`AtmMixture::split` (five of them identically zero). `ATM_BG_LUMPED=1` forces the single
`kappa_bg`, which at this composition is the same number either way.

## The sweep count: p_dyn's amplitude and Psi are DECOUPLED, and the 28 % was the initial projection

**The first version of this section attributed a 28 % drop in `Psi_max` to the time-loop solve.
That was wrong, and the confound noted at the bottom of it turned out to be the entire effect.**
`ATM_PRESS_SWEEPS` was doing two jobs — `project_initial_velocity` makes 200 passes through
`run()`, so raising the knob multiplied the startup projection by the same factor.
`ATM_PROJ_SWEEPS` now separates them. Three arms, iteration 20, one binary:

| | `p_dyn` min | `p_dyn` max | radial j=45 | latitudinal | `Psi_max` |
|---|---|---|---|---|---|
| loop 1, proj 200x1 — shipped | -0.00566 | 0.00854 | 0.00191 | 0.00814 | 40959.17 |
| loop 10, proj 200x1 — **loop alone** | -0.04454 | 0.05899 | 0.01587 | 0.08211 | **40959.80** |
| loop 10, proj 200x10 — both | -0.04454 | 0.05898 | 0.01588 | 0.08209 | **29436.77** |

**THE LOOP SWEEPS GROW `p_dyn` TENFOLD AND LEAVE THE CIRCULATION ALONE.** Isolated, ten sweeps
per iteration take the dynamic pressure from a range of 0.014 to 0.104 — the shipped field is
about a tenth of a converged one, which is real — and `Psi_max` moves by **+0.0015 %**,
40959.17 to 40959.80.

**THE 28 % IS THE INITIAL PROJECTION, ALL OF IT.** The only difference between the second and
third arms is 2000 startup relaxations instead of 200, and it costs `Psi_max` 11 523 units.
Their `p_dyn` fields at iteration 20 are identical to four digits, so the projection is not
leaving a different pressure behind — it leaves a different VELOCITY. `project_initial_velocity`
applies `v <- v - grad p` and then clears `p_dyn`; with ten times the relaxation it removes more
of the initial divergence, the circulation starts weaker, and `Psi` carries that difference
forward while the pressure re-equilibrates to the same amplitude either way.

**So the two counts answer two different questions and neither substitutes for the other.** The
time-loop solve sets the amplitude of `p_dyn` and, at least over 20 iterations, nothing else.
The initial projection sets the circulation. **Every `Psi` number in this fork therefore depends
on how well the initial velocity field was projected, and 200 x 1 was never chosen — it is the
default of a routine written for a different purpose.**

That also puts the earlier reading of `p_dyn`'s radial structure in its place: the
radial/latitudinal ratio is 0.235 at one sweep and 0.193 at ten, so more solver work builds the
whole field in proportion. Under-convergence is an amplitude problem, not a shape one, and it is
not what makes `p_dyn` two-dimensional here — the absence of a balanced initial state is.

**Limits.** 20 iterations, one arm each, and ten sweeps is not converged either: a 50-sweep arm
was abandoned when 200 x 50 startup passes proved impractical, which is what forced the knobs
apart in the first place. What is established is the DECOUPLING and the attribution, not where
either count saturates.

## The CO2 distribution work ported from ATHAD (its README items 56-57, 59)

**This fork had already found half of it, and the half it found was the right half.**
`ThermoAtm::co2Atmosphere()` did not fill a uniform mass fraction. It stored

    co2(i) = (1 - c(i)) * f_CO2,     f_CO2 = q_CO2/(q_CO2 + q_bg) at the sea = 0.9536

with a comment giving exactly the right reason: the water runs from 0.346 at the sea to 0.004
above the cold trap, so a uniform MASS fraction would make the background absorb 0.34 of the
mass and put R aloft at 230 J/(kg K) instead of 195 — an 18 % error through the whole upper
atmosphere. **ATHAD had the same defect and did not notice it for months**, because its water
field is nearly uniform (537-739 g/kg, a 1.4x range) and there the two statements coincide.

**WHAT IT COULD NOT DO.** It baked the distribution in ONCE, against the initial water field.
The moment the water evolves the stored field is stale and the background silently resumes
absorbing every change — the same defect, deferred by one initialisation. `AtmMixture::q_CO2_of`
now applies `q_CO2 = co2_stored*(1 - q_v)/(1 - c_0)` CONTINUOUSLY, so the invariant is
maintained rather than imprinted, and the stored field becomes uniform at `co2_0`.

**THE TWO AGREE EXACTLY AT t = 0**, which is the check that the port preserves what this fork
got right: `co2_0/(1 - c_0)` = 0.6233/0.6536 = **0.9536** = `f_CO2`. Measured at iteration 10,
stored field against effective field:

| at | old stored `co2` | new effective `q_CO2` |
|---|---|---|
| sea, c = 0.3399 | 0.6296 | 0.6233 x 0.6601/0.6536 = **0.6294** |
| 70 km, c = 0.0214 | 0.9333 | 0.6233 x 0.9786/0.6536 = **0.9331** |

**AND THE GAS-MASS NORMALISATION IS WHAT MOVES THE ANSWER HERE.** `split()` now takes the
suspended condensate (cloud + ice + graupel) so the fractions sum to `1 - q_cond`, and
`densities()`'s `water_factor` divisor — the same correction, applied to `r_humid` alone,
without graupel, floored at 0.5 — is deleted in the same edit so it is not counted twice.
In a condensing atmosphere with a sea this is not the second-order effect it is in ATHAD:

40 iterations, 24 threads, same config, one binary generation apart:

| iter | OLR before | OLR after | imbalance before | imbalance after |
|---|---|---|---|---|
| 10 | 281.77 | **272.95** (-3.1 %) | -10.99 | **-2.17** |
| 20 | 283.81 | **275.29** (-3.0 %) | -13.03 | **-4.51** |
| 30 | 283.82 | 275.39 (-3.0 %) | -13.04 | -4.61 |
| 40 | 283.72 | **275.36** (-2.9 %) | -12.94 | **-4.58** |

| | before | after |
|---|---|---|
| `Psi_max` @ 20 / @ 40 | 42218.93 / 42186.38 | 40959.17 / 40927.25 (-3.0 %) |
| co2 column average | 438312.257 | 438259.583 (-0.012 %) |
| max cloud water | 47.084 / 47.002 g/kg | 46.912 / 46.860 |
| max water vapour | 339.932 g/kg | 339.943 |
| mean albedo | 0.5000 | 0.5000 |

**The OLR shift is stable at -3.0 % across all four diagnostics**, so it is a property of the
change and not of the first iteration. ATHAD's equivalent was -0.11 %. **The 28x difference is
the condensate**: ATHAD's is 12-47 g/kg confined to a thin band, this fork condenses throughout.

**The imbalance is NOT a claim that the budget now closes.** It improves from -12.94 to
-4.58 W/m2 at iteration 40, but both arms are 40-iteration transients, and this file's own
history (ATHAD items 25, 43, 46) is a series of imbalances that looked like convergence and were
not. What is measured is a difference between two arms, not an equilibrium.

**A DESIGN FLAW THE COND SELF-TEST CAUGHT ON THE FIRST RUN.** `carrierRef()` was set inside
`AtmMixture::resolve()` — the one function that knows the configured composition, so it looked
like the natural home. But `resolve()` is a pure function any caller may invoke with any
composition, and `cond_column_selftest.cpp` calls it with a near-dry one. That left
`carrierRef()` at 0.9957 instead of 0.6536 and moved `M_nonwater` from 42.88 to 36.27 g/mol.
**A global written by whoever called last is not a reference.** It is now set once, explicitly,
by `initComposition()` from `c_0` — the water content `co2_0` is quoted at — and the same repair
went back into ATHAD, where the two happened to coincide and the bug was invisible.

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

**Since 2026-08-20 the five trace gases carry ATHAD's mole fractions** — CH₄, NH₃, H₂, CO
and SO₂ at 0.014 each — where they used to be zero. `x_H₂O` is fixed by the sea and cannot
give way, so they displace CO₂ and N₂ within the 44.22 % that is not water; above the cold
trap they are **15.7 % of the dry atmosphere and CO₂ falls to 77.4 %, below the quoted
89–95 %**. That is a deliberate departure from the stated composition. Note also that
matching ATHAD's *fractions* does not match its *amounts*: 1.4 % of ATHAD's 250 bar is
3.5 bar per gas, 1.4 % of this 60 bar column is 0.84 bar, so the same numbers imply four
fifths of each trace gas left with the water.

| | at the sea | above the cold trap |
|---|---|---|
| x = (H₂O, CO₂, N₂) | .5578 / .3459 / .0263 | .010 / .7744 / .0589 |
| x(CH₄, NH₃, H₂, CO, SO₂) | .014 each | .0313 each |
| M_mean [g/mol] | 27.79 | 39.90 |
| R_mix [J/(kg·K)] | 299.2 | 208.4 |
| cp [J/(kg·K)] | 1372.6 | 1042.5 † |

† **The 1042.5 is the dry composition at the SEA-SURFACE temperature**, not at the cold
trap. `cp_of(dry, 513.15 K)` = 1042.5; at cold-trap conditions it is **886.4** (the
Shomate fits clamp below 298 K). This table warns two paragraphs above that the stated
composition and the stated temperature describe different heights — and then does it
itself, for this one entry. The 1372.6 is sound: it is `cp_of` at the sea, and it is what
`cp_l` is set to.

| q_H₂O [kg/kg] | 0.3616 | 0.00452 |
| ρ [kg/m³] | 39.08 | — |
| scale height | 15.7 km | 10.9 km |

The background is no longer N₂ alone: M_bg = 26.138 g/mol, **R_bg = 318.10 J/(kg·K)**
(was 296.8), within 0.3 % of ATHAD's 317.26 because in both forks the traces are ~71 % of
the background by mass. `M_nonwater` at the sea is **40.12 g/mol** (was 42.88). And the
background's grey opacity Σf_i·κ_i goes **1e-6 → 1.944e-3 m²/kg**, so item 60's per-species
split — an exact no-op at the old composition — now carries NH₃ (48.7 %), SO₂ (36.7 %) and
CH₄ (13.8 %).

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

Three more of the shape, ported from ATHAD 2026-08-18 and **measured here rather than
assumed** — the fork inherited all three silently:

- **`exp_rm` is not the Jacobian of the radial stretch.** It is documented as one in two
  files that agree with each other and with the variable's name (`TurbulenceAtm.h`;
  `PressureSolverAtm.h` writes `dp/dr_physical = exp_rm * dp/d(rad.z)`), and it is the
  Jacobian of a *quadratic* stretch while `init_layer_heights` builds an *exponential* one.
  `checkRadialMetric()` now prints the spread every startup: **12.29×** here, the core's
  radial length unit running 9 837 m at the surface to 120 850 m at the top. So the dynamics
  and the radiation do not agree where the levels are. CLAUDE.md carried "~12×" as an
  inference from `zeta = 3.0`; the inference was right and is now a measurement. The test is
  unit-free, so it is not a units convention. **`im` is not the lever** — the spread is a
  function of `zeta` and `im` only, and 61 → 41 moves it to 11.77×. This is the pattern's
  worst form yet: not an Earth literal but a *wrong formula with a confident comment beside
  it*, which is why the check is print-only and permanent rather than a one-off measurement.
- **The surface drag is Earth's twice over.** `rayleigh_kf` was tuned against a jet off the
  west coast of South America and `drag_n_layers` defended by an Andes orographic feature;
  invariant 1 makes `is_land()` false everywhere here. Both are config parameters now,
  bit-identically (verified: a 5-iteration run against the pre-change binary differs in
  timestamps, wall clock and the new metric-check lines, and in no number).
- **And `drag_n_layers` is a cell COUNT, so the drag depth is set by the grid.** 5 cells is
  236 m on `ATOM_Precipitation`'s grid and **1 786 m here**. ATHAD reported this against its
  own 7 152 m; what its writeup did not have to say is that the depth moves with **`im`** as
  well as `zeta`, because ATHAD holds `im` at 41 and this fork is at 61. The same constant
  means 2 861 m at `im` 41, 1 786 m at 61, 1 297 m at 81 — **a 2.2× swing in a physical
  momentum sink from the level count alone.** A grid-refinement study here silently rescales
  the drag, so `im` and `drag_n_layers` cannot be varied independently until the depth is
  reformulated as a length in metres.

And one of the same shape in a field rather than a constant: ATHAD holds the CO₂ **mass**
fraction uniform, which is right there because its water is uniform too. Here water runs
0.346 → 0.004, so a uniform mass fraction makes the background absorb all the mass the
water vacates, and R settles at 230 J/(kg·K) aloft instead of 195 — an 18 % error through
the whole upper atmosphere, from the one field that is supposed to have no structure. What
is well mixed is the CO₂:background ratio *within the dry air*.

## Reproducibility

**Run-to-run bit-identical at a fixed thread count; thread-count dependent at ~1e-8.**
Measured at 24 threads, three one-iteration runs identical in every number, and at 20
iterations 24 threads against 8 agree in every wind extremum and differ in the last digit
of `residuum_atm` (0.77593489 against 0.77593490).

**It was not always.** Until the race fix ported from ATHAD (its item 18), two runs of the
*same* binary at the *same* thread count gave different answers — max u-component 0.080795
against 0.080751 at 24 threads, the extremum wandering in longitude and hemisphere. Two
defects, both inherited:

- **The Poisson loop wrote `p_dyn` in place** under `collapse(2) schedule(dynamic, 4)` over
  (i,j), while the stencil reads `p_dyn[i±1][j±1]` — across the very two indices it writes.
  `schedule(dynamic)` made it worse than a thread-count dependence: which thread got which
  chunk varied with timing, so the same binary at the same thread count varied run to run.
  Now two red-black passes over a checkerboard colouring of (i+j+k), `schedule(static)`.
- **`UtilsAtm::findResiduumAtm` wrote the shared `m.residuum_old`** from inside every
  thread's loop. It is read three times — the "declining" vs "too high" message and both
  reported errors — so a raced value drove the line a human reads to decide whether the run
  is converging. Now captured once before the parallel region and set once after the
  reduction, with a lexicographic tie-break so the reported error *location* stops depending
  on thread arrival order too.

The varying cell sat at 37 905 m, i ≈ 39 of 61 — interior, which is where an in-place
stencil race puts it and not where a boundary defect would.

**Red-black does not reproduce the old single-thread answer**: it is a different sweep order
from lexicographic Gauss–Seidel and converges to the same solution by a different path.
Every number in this file measured before the fix moves in its last digits.

**Not fixed — floating-point reduction order.** OpenMP combines partial sums in a
thread-count-dependent order and `+` is not associative. Curing it needs ordered reductions,
as in ATHAD.

**What reproducibility bought immediately.** The cell-structure port could not be verified
when it landed — there was nothing to compare against, because two runs of the same binary
already disagreed at ~1e-3. With the races cured the check is possible and it passes: a
pre-port binary *carrying the same race fixes* against the ported one, both at Earth
settings (`cell_lat_scale` 1.0, `n_cells_hemisphere` 3), 20 iterations, 24 threads, agree in
**every number of the entire log**. The only lines that differ are the three diagnostics the
port adds. That is the verification the port's own commit message said it could not make,
and it is a fair illustration of what a race costs: not just wrong answers, but the loss of
the instrument you would use to detect them.

## What the surface state still breaks — audit, 2026-08-18 (items 1-3 RETRACTED same day)

The surface moved from ATHAD's **250 bar / 1500 K / supercritical** to **60 bar / 513.15 K /
subcritical**. This section first claimed three live defects in
`ThermoAtm::waterVapourEvaporation()`. **Two of the three were not live and the third could not
happen. The retraction is kept above the surviving findings because the way it was got wrong is
the more useful result.**

### RETRACTED: the evaporation routine was already repaired, in the branch that runs

The claim was that `waterVapourEvaporation()` used a scalar `ep = 0.6431` that excludes CO₂
(giving q_sat = 0.4479 against the correct 0.3464, +29.3 %), that it inverted q → e with the
dilute form while computing e → q exactly, and that the active Meyer formula turned the
resulting spurious 1151 hPa deficit into ~317 mm/day of evaporation from a saturated surface.

**The arithmetic is all correct and describes code that never executes.** `ThermoAtm.h:212-257`
is an ATHAD_COND-specific branch that already does the right thing — per-cell
`AtmMixture::M_nonwater`, the exact q → e inverse written out inline, and a relaxation of
`c.x[0]` toward `q_sea` instead of an empirical flux — and it ends in `continue` for every
subcritical cell. The supercritical guard above it catches the rest. **No cell reaches the
`ep` code.** Its own comment says so in capitals: *"ATHAD_COND KEEPS THIS GUARD … and adds the
branch below, which is the one that runs."* The self-test's warning at
`cond_column_selftest.cpp:158-165` had been acted on; what it protects is that branch.

**And the third claim was impossible, not merely inert.** `evap_model` was said to be one flag
from a Rohwer sign flip. The live branch ends `m.Evaporation.y[j][k] = m.Evaporation_Dalton...`
— it **hardcodes Dalton and never reads `evap_model` at all**. Measured: the config says
`Meyer`, and `max Evaporation` = `max Evaporation Dalton` = 62.833923 mm/d exactly, while
`Evaporation Meyer` = 58.636186. **`evap_model` is a config parameter that does nothing**,
which is a real finding, and not the one that was claimed.

**How it was got wrong: liveness was inferred from a grep instead of read from the control
flow** — the exact failure CLAUDE.md names ("a term written in `RHS_Atm_Turb.cpp` is not
necessarily a term the model applies"), committed while quoting that lesson. The grep found
`m.ep` at nine sites in one routine; reading 80 lines further back would have found the
`continue`. **The generalisation that does survive is narrower and still worth having: an
assertion in a test is not a check on the code it describes** — `cond_column_selftest.cpp`
states this error to the digit but asserts against a value it computes itself, so it would not
have caught the dead path had it been live.

Rohwer's Earth regression is still worth recording as arithmetic, since it is printed every
run: `(1.465 − 0.000732·p_mmHg)` is +0.909 at Earth sea level, **−31.478** here and −135.796 in
ATHAD, and `Evaporation_Rohwer_average` duly prints **−2497.795 mm/d**. It is a diagnostic that
cannot be selected, not a latent switch.

### What was actually fixed

**The precipitable-water diagnostic was on the dilute inverse, and it is live.**
`ThermoAtm.h:699` computed `e = q·p/ep` inside the column integral. Replaced by the exact
inverse on per-cell `M_nonwater`. Measured: **135 791 → 153 529 mm, +13.1 %**, and it is the
*only* number in the whole model output that moves — everything else is bit-identical over a
full run.

`SaturationH2O::vapourPressureFromMassFraction(q, p, M_other)` is new: the exact inverse of
`saturationMassFraction`, which the file was missing. The live evaporation branch had written
that expression out inline; it now calls the helper, bit-identically, so the two cannot drift.
The dead Earth path was converted to the same exact pair — **a trap removal, not a fix, and it
changes no output.**

### Still open, and these two are verified live

**A. `p_stat.x[0]` is re-anchored to `r_air·R_mix·T_surf`** (`ThermoAtm.h:1456,1491`;
`InitValues_Atm.cpp:542`). **Measured, not read**: the printed `max pressure static` is
59.919 bar **at 0°N**, so the surface pressure carries the equator-to-pole temperature contrast
— 58.75 bar at the pole against 59.92 at the equator, a 1.95 % latitudinal variation imposed by
construction. Under ATHAD's item 19 it would be `p_0` everywhere.

The design point is consistent to 0.007 % (1e-2·40.8·286.6·513.15 = 60004 hPa against
`p_0` = 60000), so porting item 19 would not move the initial state, only the drift. **But it
must not be ported blindly.** ATHAD's argument is "if no mass enters or leaves, `p_s` is a
constant", and here mass *does* leave — rain reaches a sea that is a boundary condition, not a
reservoir. A constant `p_0` would assert a conservation this configuration does not have.
Decide what `p_s` should do under net precipitation first.

**B. Stale ATHAD arithmetic in comments on live code.** `ThermoAtm.h:1452-1455` justifies
`R_mix` by "using `R_Air` here instead yields 204 bar rather than the intended 250";
`MultiLayerRadiation.h:137,152` explains pressure broadening as "at 250 bar … a factor of 250
over the 1 bar reference". Both describe ATHAD; here it is 60 bar. The code is right.
`p_ref = 1.0e5` Pa is a genuine broadening reference and is fine.

**C. Dormant.** `OneCatIceScheme.h:239,242` tests `p_stat <= 500.0` hPa absolute;
`CategoryIceScheme` is 2, so it never runs — verified the same way this time.

**Already repaired, verified in passing**: the deep-convection triggers (fractions of
`p_stat.x[0]`), `initCloudIce`'s `H_crit` parabola (anchored to `p_0`), and
`latentSensibleHeat`'s sea-surface humidity (per-cell `M_nonwater`).

## Remaining work

- **WHY Psi DOES NOT CLOSE AT THE GROUND, AND THE TROPOPAUSE CONSTANTS THAT CAUSE IT**
  (2026-08-24). `Psi(i=0)` must be zero: `u` at the surface is 0 and `Psi` at the lid is 0, so
  the column-integrated meridional MASS flux is forced to vanish. It does not, and the reason
  is not the solver, not the diagnostic, and not the projection.

  **Psi(ground) can only vanish if the two branches of the cell carry equal and opposite
  MASS.** Compared at the same latitude (+2 deg) with ATHAD:

  | | v reverses at | **mass above the reversal** | \|net\|/gross |
  |---|---|---|---|
  | ATHAD | 49.6 km | **54.9 %** | **0.097** |
  | this fork (was) | 59.6 km | **0.04 %** | 0.9992 |
  | ATHAD_PERID (was) | 59.6 km | **0.00 %** | 1.0000 |

  In ATHAD the branches nearly cancel. Here the return branch existed in VELOCITY and not in
  MASS FLUX, so `Psi(ground)` was the whole circulation by construction. In scale heights the
  reversal sat at **0.84 H in ATHAD, 3.85 H here, 7.7 H in ATHAD_PERID**.

  **TWO DEFECTS, AND THEY ARE NOT THE SAME SIZE.**

  **(A) `tropopause_equator/pole` were ATHAD's 207/195 km, ABOVE THIS MODEL'S 120 km LID.**
  Inherited on the fork commit and never re-sized. `height_to_level(207000)` returns 70.5 and
  clamps to `im-1 = 60`, which the startup line has been printing all along:

      tropopause_pole=195000  pole_index=60.000000  (height 120000.3 m of 120000.3 m)   <- here
      tropopause_pole=195000  pole_index=35.000000  (height 201275.1 m of 300005.6 m)   <- ATHAD

  Two consequences: the initial wind ramp was stretched over the whole shell, and
  `VelocityInitializer::init_v_or_w_above_tropopause` opens with `if (tl >= m.im - 1) return;`
  so **the taper above the tropopause never executed in either fork.** Same defect class as
  `init_tropopause_layers`'s `round(h / L_atm)`: a height constant that is right in the parent
  and becomes a wrong level index in a fork with a different shell.

  **(B) `init_v_or_w` is linear in GEOMETRIC HEIGHT**, so the reversal sits at a fixed fraction
  of the tropopause height while the mass is distributed exponentially. ATHAD's cancellation is
  a coincidence of its 59 km scale height; this fork's is 15.5 km.

  **MEASURED, 40 iterations, 24 threads, config-only arms:**

  | tropopause | tl | v reverses | mass above | RMS Psi(0) | interior Psi | ratio | Psi_max at | div rms |
  |---|---|---|---|---|---|---|---|---|
  | 195/207 km (ATHAD) | **60 = lid** | 59.6 km | 0.04 % | 8.79e12 | 3.58e12 | **2.456** | 0 m | 7.570e-02 |
  | 70/74.3 km (thermal) | 50 | 40.2 km | 0.80 % | 6.55e12 | 3.03e12 | 2.164 | 0 m | 6.744e-02 |
  | **22/23.4 km (mass)** | **30** | **7.0 km** | **79.8 %** | 1.03e13 | **1.49e13** | **0.693** | **7.0 km** | **2.825e-02** |

  **Fixing (A) alone is not enough here.** This model's own thermal tropopause is 70 km -- its
  profile goes isothermal at 221.1 K from level 50 -- and putting it there still leaves only
  **0.80 %** of the mass above the reversal. It is (B) that governs: placing the tropopause by
  MASS (p = 0.1 p_surf, 21.9 km) puts 79.8 % above it and

  - the closure ratio falls **2.456 -> 0.693**;
  - **the global maximum of `Psi` LEAVES THE GROUND for the first time in this fork** (2.74e13
    at 7.0 km) -- the signature ATHAD item 68 used to argue `Psi_max` had stopped reporting the
    defect. Here it was still reporting it until now;
  - the real interior circulation grows **4.1x**;
  - **`div(rho u)/rho` rms falls 7.570e-02 -> 2.825e-02**, BELOW the pre-port branch's
    4.864e-02 -- so most of the divergence rise attributed to `ATM_PROJ_SWEEPS` was a symptom
    of an initial condition carrying a huge mass-flux imbalance, not a cost of the extra sweeps;
  - the OLR does not move (196.44 -> 196.42), the family's null yet again.

  **22.0/23.4 km IS THE DEFAULT SINCE 2026-08-24**; `tropopause_pole=195000
  tropopause_equator=207000` restore the old branch exactly. **HONEST CAVEAT, recorded in
  `param.py` too**: this column is still condensing at 60 km, so 22 km is a good initial-wind
  scale and a BAD description of the convective top. The parameter does double duty and the two
  jobs disagree here. **The real repair is to express the cell profile in a mass coordinate**
  and let the parameter mean what its name says; until then it is set for the job it does.
  ATHAD_PERID needs no such compromise -- its thermal tropopause IS very nearly a mass one.

- **THE ONE-KNOB DECOMPOSITION OF THE PORT** (2026-08-24, 40 iterations, 24 threads, one binary,
  each arm turning on exactly ONE knob from the old branch):

  | arm | knob | OLR | delta | div rms | Psi_max |
  |---|---|---|---|---|---|
  | old | -- | 172.69 | -- | 4.864e-02 | 40 516 |
  | P | `ATM_PROJ_SWEEPS=10` | 172.69 | **0.00** | **7.624e-02** | **29 065** |
  | R | `ATM_RAD_DIRECT` | 175.42 | **+2.73 (+1.6 %)** | 4.863e-02 | 40 501 |
  | S | `ATM_SAT_SUPERHEAT` | 172.69 | **0.00 (exact null)** | 4.864e-02 | 40 516 |
  | B | `ATM_PRECIP_BANDS` | **197.56** | **+24.87 (+14.4 %)** | 4.825e-02 | 40 514 |
  | new | all four | 196.44 | +23.75 | 7.570e-02 | 29 054 |

  **Every effect is single-valued.** The divergence rise is entirely `ATM_PROJ_SWEEPS`; the
  other three are null to four digits. `ATM_SAT_SUPERHEAT` is an exact null in every printed
  field, confirming the census end to end. **And the OLR is `ATM_PRECIP_BANDS`, not the
  radiation solver** -- an earlier reading of the four-knob pair as a solver-convergence effect
  scaling with optical thickness is **WITHDRAWN**: `ATM_RAD_DIRECT` is 1.6 %, one ninth of it.

  **The bands result is two instruments coming into agreement, not "the OLR rose."** The log
  carries two independent OLR calculations with the same formula and different domains:
  `printColumnProfile` at the **equator only** (j = 90, k = 0) and `printPlanetaryBalance` as
  the **global cos-weighted mean**.

  | | equator column | global mean |
  |---|---|---|
  | old | **196.6** | 172.69 |
  | bands on | **196.6** | 197.56 |

  The equatorial value is identical in every arm to 0.1 W/m2. What changed is **everything
  except the equator**: the extratropics were radiating 12 % below the tropics and now radiate
  the same. Supporting evidence: the equatorial column is unchanged to four digits, and
  `max precipitation total` relocates from 69 deg N at 66.6 km to 61 deg N at the ground.
  **The microphysical path is NOT established** -- the leading explanation is that the old
  bands destroyed snow outside 253.15-273.15 K, so condensate accumulated in the colder
  extratropical upper column instead of falling out, but the per-latitude condensate
  comparison has not been done and two earlier guesses at this mechanism were wrong.

- **THE ATHAD PORT OF 2026-08-24 (its README items 68-80): four defaults changed, seven knobs
  added, and the measurement is OPEN.** The two trees' physics files are kept name-identical so
  fixes cherry-pick in both directions; this is the second batch to travel (the first was item
  67's grey skin factor, 2026-08-21). **What transfers is each repair's argument. What does not
  transfer is its size** — ATHAD is 250 bar over a 1500 K melt.

  **The four that change results.** Restore all of them with
  `ATM_PROJ_SWEEPS=1 ATM_RAD_DIRECT=0 ATM_SAT_SUPERHEAT=0 ATM_PRECIP_BANDS=0`.

  1. **`ATM_PROJ_SWEEPS` 1 -> 10** (item 68). The initial pressure projection ran one relaxation
     sweep per pass. `Psi(ground)` must be identically zero — `u` at the surface is 0 and `Psi` at
     the lid is 0, so the column-integrated meridional mass flux is forced to vanish — and in
     ATHAD it was **2.09x the interior circulation**. Ten sweeps remove 52.5 % of that and then
     plateau, so ~44 % is structural and unexplained (ATHAD item 72: the projection has CONVERGED
     to a fixed point that is not divergence-free, and 64x the sweeps changes it by nothing).
     The cost is a one-time startup expense. **In ATHAD this also changed what `Psi_max` MEANS**:
     at one sweep the global maximum of `Psi` WAS the spurious surface flux, so `Psi_max` was
     reporting the defect rather than the circulation, and the two maxima only coincide at 10.
  2. **`ATM_RAD_DIRECT`, new and ON** (items 30, 71). The Lambda iteration is Jacobi on an
     `im`-link chain: information moves one layer per sweep, so it needs O(N^2) sweeps, and the
     inherited `n_lambda = 4` is an Earth constant — a loop bound nobody reads as a physical
     assumption. The system does not need iterating: `a_i + b_i = 1` makes the net flux constant
     with height, which closes it in two O(N) passes, exact. ATHAD measured the Lambda iteration
     converging MONOTONICALLY onto the closed form (243.43 -> 227.55 -> 212.54 -> **211.57** W/m2
     at 4/64/512/direct), which is the test that matters: the closed form is the answer the
     sweeps are trying to reach, not an alternative to them. It is also free — 281 s against
     277 s for four sweeps and 467 s for 512. `ATM_N_LAMBDA` arrives with it.
  3. **`ATM_SAT_SUPERHEAT`, new and ON** (item 75). `SaturationAdjustment::clampAndFade` tested
     whether a cell could hold a condensed phase, condensed, added the latent heat **that makes
     its own answer false**, and never re-tested; `IceSchemeCommon::evaporateWhereImpossible`
     then deleted the result, correctly, every iteration. Work done and undone, invisible because
     the diagnostics print after the whole moist block. The guard re-tests after the write-back
     and rejects the step whole, so what a rejected cell keeps is its supersaturation — the
     honest state of a parcel that cannot condense — rather than a manufactured phase. **This is
     the largest single effect ATHAD has recorded: OLR -49.4 %, photosphere +44 km, max cloud
     water 0.000000 -> 40.09 g/kg.** It is also what refuted that tree's "microphysics is
     unmeasurable" wall: five consecutive null repairs had been measuring an annihilation.
  4. **`ATM_PRECIP_BANDS`, new and ON** (item 76). Each precipitation category was written
     `(band) ? (inherited + produced) : 0`, which conflates "can this phase be PRODUCED here"
     with "can a flux PASS THROUGH here". The second has no temperature bound — falling ice does
     not cease to exist because the air it passes through is cold — and the `: 0.0` DESTROYS a
     flux arriving from the level above. The bands bottom out at Earth's 236.15 K and 253.15 K.
     Fixed: the inherited flux always passes, only production is gated, and snow loses its
     -20 C floor. Graupel keeps its -37 C floor, which is physical (riming needs supercooled
     liquid); rain keeps 273.15 K for the same kind of reason.

  **Default-off, off-branch bit-identical:** `ATM_CELL_ALTERNATE` (item 69 — `centreAmp` fills
  every middle cell with a Ferrel copy and the polar template already carries the Ferrel sense,
  so four of five prescribed cells turn the same way and their mass fluxes ADD; `edgeRadialCoeff`
  already assumes the parity the cores do not impose), `ATM_METRIC_EXACT` (item 80 —
  `exp_rm = 1/(rm+1)` is the Jacobian of a QUADRATIC stretch applied to an exponential grid,
  documented as the transformation's Jacobian in two files that agree with each other and with
  the variable's name; `metricExpRm()` replaces all eleven sites and `metricCurv()` adds the
  `-(J'/J)f'` curvature term the Laplacian omits under BOTH metrics. **Default off because in
  ATHAD the correct Jacobian made the projection WORSE** — `div(rho u)/rho` rms 2.7e-02 ->
  7.7e-02 — with the OLR unmoved), `ATM_PRECIP_CAP` (item 77).

  **Print-only:** the staged `ATM_ICE_CENSUS` (entry / post-adjustment / post-`damp_wiggles` /
  post-ice-scheme / leaving the block / at the diagnostic), its `canCondense`-violation counter
  and sample dump, the `P_rain` cap probe and the `S_r` term decomposition (items 74, 77, 78).
  Plus the zonal ParaView writer's **true-height vertical axis** (item 70): it wrote the vertical
  coordinate as LEVEL INDEX, which on an exponentially stretched grid distorts by ~18x and varies
  with altitude, so no ParaView aspect setting could undo it — contours, glyph angles and
  streamline curvature all inherited it. The glyph vector was also scaled `1/u_0` where the
  scalars beside it used `u_0` (64x, direction unaffected) and was raw m/s on an index-space
  geometry; `uv_plot` now carries the field in plot units per day. **No physics reads these
  fields**, so no computed result changes — but any conclusion drawn by eye from an older zonal
  plot should be re-examined. `psicheck.py`, `survival.py` and `mksweeps.py` come with it,
  adapted to this tree's `im = 61` and config name.

  **Deliberately NOT ported.** `ATM_SKIN_TAU` (item 73) was measured in ATHAD and **withdrawn**:
  `T_rad(tau)` exceeds the prescribed profile at all 41 levels there, so there is no
  radiative-convective crossing to switch at, and the whole-column form runs away to 15 599 K in
  four iterations. It also presumes ATHAD's isothermal lid over a prescribed DRY adiabat and its
  `ATM_PROGNOSTIC_T` knob; this tree has neither — invariant 4 here is a MOIST, state-dependent
  adiabat. `moist_phys_start_iter = 0` (item 74) has been this tree's default since the fork.

  **Both trees build clean and `make test` passes with 0 failures.**

  ### THE MEASUREMENT, RUN 2026-08-24 — 40 iterations, 24 threads, one binary, env-only arms

  | | old branch | new defaults | |
  |---|---|---|---|
  | OLR | 172.69 W/m2 | **196.44** | **+13.8 %**, and still rising (182.5 -> 192.1 -> 196.3 -> 196.4) against an old arm flat from iteration 10 |
  | imbalance (in - out) | +98.08 | +74.34 | |
  | photosphere | 65.7 km / 261.03 K | 65.8 km / **264.72 K** | +3.7 K, height unmoved |
  | `t_skin` | 262.88 K | 262.88 K | unmoved, as in every ATHAD arm |
  | mean albedo | 0.5000 | 0.5000 | saturated on presence, as before |
  | max cloud water | 49.03 g/kg | **50.000000** | **AT `cloud_cap`** — see below |
  | max cloud ice | 2.475 g/kg | 3.245 | +31 % |
  | `P_snow` mean | **0.000e+00 mm/a** | **9.444e+04** | the band repair, and it is not subtle |
  | total precip mean | 9.461e+04 mm/a | 2.025e+05 | x2.14 |
  | `Psi_max` | 40 516 (1e9 kg/s) | 29 054 | **-28 %, and still AT z = 0 m** |
  | `div(rho u)/rho` rms | 4.864e-02 | **7.570e-02** | **+56 %, from the first diagnostic on** |
  | `residuum_atm` | 0.094 | 0.047 | -50 % |

  **THREE RESULTS THAT ARE NOT ATHAD'S, AND ONE THAT IS.**

  **(1) `ATM_SAT_SUPERHEAT` IS INERT HERE. The largest effect ATHAD has ever recorded is a
  NO-OP in this tree, and that is measured, not assumed.** 10 iterations with
  `ATM_ICE_CENSUS=1`: the guard is REACHED 1 140 760 times in the first call and 3.4 million
  by the third, and it REJECTS **zero** every time. The staged census agrees from the other
  side — `condensate where canCondense is false` is **0 cells, 0 kg/kg, at all four stages**
  (entry, post-`adjustSaturation`, post-`applyTopography`, post-`clampAndFade`). Nothing in
  this column ever condenses into a state the ice scheme would have to undo. **So ATHAD's
  item 75 is to this fork what `alpha_entry` was** (see the saturation-adjustment entry above):
  the same code, the same call, and a defect that cannot fire because the cell it needs is
  40-225 K away. It is kept ON so the two trees stay identical and so the guard is present if
  the column ever reaches those states; it costs one `saturationMassFraction` call per Newton
  pass and changes no number here. **Corollary: none of the OLR change above is item 75's.**

  **(2) THE PRECIPITATION FLUX IN THIS MODEL IS THE CAP, BY SEVEN ORDERS OF MAGNITUDE.** This
  is what the item-77 probe was ported to ask, and the answer here is the opposite of ATHAD's
  (there the cap binds on the first ice-scheme call and never again). Measured, every call:

      CAP PROBE: P_rain cells at the cap 6 127 404,
                 largest value the recurrence wanted 2.795e+04 kg/(m2 s)   (cap 3.000e-03)

  **A factor of 9.3 million.** `P_rain` and `P_snow` are both pinned at 9.461e+04 mm/a — mean
  equal to max, i.e. a uniform field sitting on `P_max_flux` = 3.0e-3 kg/(m2 s) ~ 260 mm/d,
  Earth's "well above any physical precip". So the x2.14 in total precipitation is not a
  measurement of the band repair's size: the repair unlocked a SECOND capped category, and
  both categories are at the ceiling. **No precipitation number from this tree means anything
  until `P_max_flux` is sized for a 60 bar atmosphere over a 513 K sea.** `ATM_PRECIP_CAP`
  scales it and arrived with the same port. This is item 52's rule — *look at what the caps
  are holding back* — and this cap is holding back everything.

  **(3) ATHAD'S `Psi_max` RESULT DOES NOT REPRODUCE. `Psi_max` is still reporting the
  defect.** In ATHAD, ten projection sweeps moved the global maximum of `Psi` off the ground
  and into the interior, which is what licensed reading `Psi_max` as the circulation again.
  Here the maximum is **at z = 0 m in BOTH arms** — 40 516 at the ground, 29 054 at the
  ground. The surface flux falls 28 % and still dominates the field. **Do not import ATHAD's
  "the two maxima coincide now" sentence into this tree.**

  **(4) AND `div(rho u)/rho` GOT WORSE, IN BOTH SIBLINGS, BY ABOUT THE SAME FRACTION.**
  4.864e-02 -> 7.570e-02 here (+56 %) and 7.068e-02 -> 1.168e-01 in ATHAD_PERID (+62 %),
  present at the FIRST diagnostic and persisting to the last. **The suspect is
  `ATM_PROJ_SWEEPS`, and the evidence is PERID's**: its OLR moved -0.5 % across the same four
  knobs while its divergence rose 62 %, so the rise cannot be a radiative side-effect there.
  **IT IS NOT ATTRIBUTED — no one-knob arm has been run** — and ATHAD never measured `div`
  against this knob either (its item 72 varied `ATM_PRESS_SWEEPS`, the TIME LOOP, and found
  the residual bit-identical under 64x). If it holds up, converging the initial projection
  makes the running divergence worse while making `Psi(ground)` better, which is the same
  shape as ATHAD item 80's metric result and belongs to the same open question.

  **What the OLR change is NOT.** Not item 75 (inert, measured above). The remaining
  candidates are `ATM_RAD_DIRECT`, `ATM_PRECIP_BANDS` (through the long-wave cloud opacity —
  cloud ice is +31 %, and ATHAD item 75 established that path responds to AMOUNT), and
  `ATM_PROJ_SWEEPS`. The split across the two forks is suggestive — **+13.8 % at 60 bar
  against -0.5 % at 12.5 bar**, which is the direction a solver-convergence effect would take,
  since four Lambda sweeps fall further short the thicker the column — but that is an
  inference from two points. **Note also that the SIGN differs from ATHAD**, where converging
  the solver LOWERED the OLR by 15 %. **The one-knob decomposition is the next run**: four
  arms, ~6.6 min each.

- **The saturation-adjustment knobs are ported, and this fork is the control that settles what
  `alpha_entry` is worth** (ATHAD README item 64). `ATM_SAT_TRACE=1` (print-only) and
  `ATM_SAT_NO_ALPHA=1` (default off); the traced levels default to 20 (10.8 km) and 48 (63 km,
  the 148.9 % tail) here against ATHAD's 37/38, settable with `ATM_SAT_TRACE_I1/I2`.

  | | level 20 — 27.8 bar, 462 K | level 48 — 0.075 bar, 279 K |
  |---|---|---|
  | entry `q_v/q_sat` | 1.00 | 1.00 → 1.07 |
  | `alpha_entry` | **1.000000** | **0.9998** |
  | latent heating over the call | +0.05 K | +0.78 K |
  | condensed | 0.0459 kg/kg | 0.000398 |
  | exit `q_v` vs `q_sat` | equal to six digits | equal to six digits |

  **`alpha_entry` never engages here** — the −37 °C ice threshold it doubles is 40 to 225 K below
  this fork's condensing column, so the master gain is 1 and the defect is inert. And the
  adjustment converges cleanly to `q_v = q_sat`: none of ATHAD's pathology appears, because
  condensing 0.046 kg/kg at 27.8 bar releases **0.05 K** rather than the 76 K that throws
  ATHAD's level 38 past boiling and inverts its target. **Same code, same call, opposite
  behaviour — set by the pressure and the temperature of the cell, not by the scheme.** That is
  the control ATHAD item 64 asked for, and it confirms its conclusion: ATHAD's supersaturation
  is a temperature problem, not a saturation-adjustment problem.

- **`M_max` SIZED, AND THE GEOPOTENTIAL IS ON BY DEFAULT** (2026-08-20). Two changes that
  belong together, because the first was suppressing the second by 4.3x.

  **`mc_M_max` is a config parameter now, at 100 kg/(m²s).** It was a bare 3.0 — Earth's, whose
  own comment reads "~10x any realistic value" with "healthy ~0.3", true at 1.2 kg/m³ — written
  in **two** places: the namespace constant `clamp_M` reads, and a second shadowing `constexpr`
  inside `rhsForcing`. Both 3.0, agreeing by luck. Measured, `M_u` sat at **exactly 3.0000 from
  2.2 km to 21.9 km in up to 177 of 181 columns**: the vertical structure of the convective mass
  flux was not computed, it was the cap.

  Sized two ways that agree. **Density**: `M = ρ·σ·w` and only ρ changes between the planets;
  Earth's 3.0 at 1.2 kg/m³ implies a ceiling `σw` = 2.5 m/s, which at this model's cloud-base
  density (36.2 at 2.2 km) is 90 and at the sea surface (39.3–42.1) is 98–105. **What the scheme
  wants**: with the cap lifted to 1e9, `|M_u|` is bounded and steady — 13.7 kg/(m²s) at iteration
  0, 13.0 at 40, no runaway — i.e. `σw` ≈ 0.4 m/s, the same `σw` Earth's *healthy* convection has.
  Earth's cap sits 10x above its healthy value; 10x of 13.7 is 137. **100.0 is 33x Earth's, the
  density ratio, with 7.3x headroom.** Verified: at 100 every number matches the uncapped run to
  the digit, so the cap is a backstop again rather than the answer.

  **`ATM_MC_GEOPOTENTIAL` is on by default here** (`=0` restores the old behaviour), and the pair
  re-measured at the new cap — one binary, env-only, 40 iterations:

  | | off | on |
  |---|---|---|
  | **`max c_u`** | 0.000000 g/kg/s | **3.116e-03 @ 24 854 m** |
  | `max s_u`, and where | 2.545962 at **661 m** | 2.911126 at **24 854 m** |
  | `max MC_t` | 4.383e-03 K/s | 5.080e-03 (+15.9 %) |
  | `max M_u` | 13035.640441 g/m²s | 13041.235334 (+0.04 %) |
  | OLR / photosphere / `Psi_max` / precipitable water | 273.70 / 65.9 km / 40515.99 / 160626.497 | **all identical** |

  **`c_u` = 3.12e-03 is 4.3x the 7.19e-04 first reported, and 4.3x is exactly what the cap was
  suppressing the mass flux by** (3.0 against the natural 13.0) — `c_u = dcond·M_u/(ρ·step)`, so
  the condensation rate scales with the flux. **The earlier figure is superseded.** The rest holds:
  condensation at 24.9 km, `max s_u` moving from cloud base to column top, and nothing integrated
  moving — the albedo wall. `max M_u` differing by 0.04 % between arms confirms this is a genuine
  pair and not a re-measurement of a cap. Runs: `run_mmax_free`, `run_mmax100`, `run_pair_{on,off}`.

  **ATHAD keeps the geopotential OFF**: its updraft is one grid level deep and has no ascent for
  it to act on. The code is identical in both trees; only the default differs, deliberately.

- **`ATM_MC_GEOPOTENTIAL` PAYS HERE, AND IT IS THE FIRST UPDRAFT CONDENSATION IN EITHER FORK**
  (2026-08-20, ATHAD README item 63). Item 53 predicted that adding `g·z` to the static energy
  would let the rising parcel cool and condense; ATHAD tried it and got nothing, because its
  updraft is one grid level deep. This fork's is **20 km deep** — cloud base 1786–2200 m at the
  sea, LFS 20 517–28 130 m, 5218 of 11041 cells in the zonal slice convecting — because the
  convective triggers here have been fractions of surface pressure since the fork was cut.

  **MEASURED, 40 iterations, ONE binary, env-only A/B (`run_alb006` vs `run_gz_cond`):**

  | | off | on |
  |---|---|---|
  | **`max c_u`** | **0.000000 g/kg/s** | **7.19e-04 @ 24 854 m, 37°S** |
  | `max s_u`, and where | 2.5460 at **661 m** | 2.9111 at **24 854 m** |
  | `max MC_t` | 1.030e-03 K/s | 1.175e-03 (+14 %) |
  | `max q_c_u` | 7.0537 g/kg at 1017 m | 7.1392 at 1785 m |
  | OLR / photosphere / `Psi_max` / precipitable water | 273.70 / 65.9 km / 40515.96 / 160626.497 | **all identical** |

  The condensation appears **25 km up**, which is where a parcel cooling at `g/cp_l` = 4.81 K/km
  for 23 km first saturates, and `max s_u` moves from the cloud base to the top of the column —
  the same statement read off the other field, since with the geopotential in, static energy is
  largest where `g·z` is. **Nothing integrated moves**: the albedo wall again, as with every
  microphysics correction since item 51.

  **Still open, and now the loudest thing in this scheme: `max M_u` is 3.0000 kg/(m²s), PINNED at
  `M_max`.** Half the convective mass flux in this model is a cap, not a result — and `M_max` = 3.0
  is another Earth constant ("~10× any realistic value" on a 1 bar planet). Size it before reading
  anything quantitative off the convection here. That is the third time a cap has turned out to be
  standing where a measurement should be (`MCt_max` in item 53, `cc_factor`'s reference in the same
  item, this one now).

- **The startup energy-balance check was reading 62 % of the insolation, and its albedo was a
  literal** (ATHAD README item 62, repaired in both trees 2026-08-20). `cAtmosphereModel.cpp`
  took the planetary mean of the insolation parabola as `0.5·(equator + pole)` = 149.0 W/m²,
  where the cos(latitude)-weighted mean is `equator·8/π² + pole·(1 − 8/π²)` = **241.55 W/m²**,
  and it used a bare `0.08` commented *"molten surface"* instead of `albedo_surface`. The
  printed estimate goes 287.08 → **372.23 W/m²** and 266.75 → **284.64 K**, so the gap it
  reports against the configured `t_skin` = 254 K goes 12.75 → 30.64 K. **Diagnostic only** —
  `planetaryShortWave()` always weighted correctly, so no result moves.
- **`albedo_surface` is a water value now — 0.06 — and it is not a lever at all** (2026-08-20).
  It was 0.08, justified in `param.py` by *"a quenching silicate melt is dark, measured
  basaltic-melt albedos are 0.05–0.10"*: ATHAD's magma ocean, inherited unchanged into a model
  whose surface is a 240 °C sea. 0.06 is water — a calm sea is 0.03–0.06 broadband at small
  zenith angle, rising toward 0.06–0.08 once high-latitude incidence and whitecaps enter, and
  the low-latitude weighting is the right one because a zero-obliquity planet puts 81 % of its
  insolation (8/π², the same factor as above) in the tropics and subtropics. The
  `MultiLayerRadiation` comment block was ATHAD's magma text verbatim, down to *"at 1500 K
  nothing within 1200 K of the ice thresholds"*; it now describes this sea, and the variable is
  `alb_surface_clear` in both trees rather than `alb_surface_molten`.

  **MEASURED — 40 iterations, one binary, config-only A/B (`run_alb006` vs `run_alb008`), and
  the result is a null so complete it is worth stating precisely:**

  | | 0.08 | 0.06 |
  |---|---|---|
  | albedo field, max and min | 0.499996 / 0.499996 | **0.499996 / 0.499996** |
  | mean planetary albedo | 0.5000 | 0.5000 |
  | OLR at iteration 20 / 40 | 273.47 / 273.70 W/m² | **273.47 / 273.70** |
  | photosphere, skin %, `Psi_max`, cloud water | 65.9 km, 58.6 %, 40515.96, 49.080208 g/kg | **all identical** |
  | max `tau_above` | 2696177.796942 | 2696177.800091 (1.2e-09) |

  **A 25 % change in the sea's albedo is invisible, and the mechanism is exact.** The cloud bump
  composites as `alpha_eff = alpha_surf + (alpha_cloud − alpha_surf)·refl` with
  `refl = tau/(tau+2)` on the condensate path, and this fork's deck saturates it: back out `refl`
  from the printed field and it is **0.9999909**. The surface therefore contributes
  `(alpha_cloud − alpha_surf)(1 − refl)` ≈ 4e-06 of the answer, and the two arms differ by
  ~2e-07 in albedo — which is why everything downstream agrees to nine significant figures. The
  albedo field is also *uniform*: max equals min to six decimals, so there is no latitude
  gradient left for a surface value to show through.

  **So this is `albedo_cloud` = 0.50 IS the planetary albedo, restated from the other side.**
  Fixing `albedo_surface` was worth doing because a constant carrying the parent's conditions is
  a defect whatever it is worth numerically — but nobody should expect a number from it while
  the deck is total. **The lever is `albedo_cloud`, and the question underneath it is whether a
  240 °C sea really maintains an unbroken optically thick deck at every latitude.**

- **THE FIVE TRACE GASES ARE IN, AT ATHAD'S MOLE FRACTIONS** (2026-08-20). CH₄, NH₃, H₂, CO
  and SO₂ carry 0.014 each here now, where they were zero and documented as "gone: the epoch
  is oxidised and degassed". **Every number in this file measured before this date was measured
  on the old composition.** What changed and why it is not a free parameter swap:

  | | without traces | with traces |
  |---|---|---|
  | x(H₂O) / x(CO₂) / x(N₂) at the sea | 0.5578 / 0.4109 / 0.0313 | **0.5578** / 0.3459 / 0.0263 |
  | M_mean [g/mol] | 29.009 | 27.789 |
  | R_mix [J/(kg·K)] | 286.61 | 299.20 |
  | M_bg / R_bg | 28.014 / 296.80 | 26.138 / **318.10** |
  | q_H₂O / q_CO₂ at the sea | 0.3464 / 0.6233 | 0.3616 / 0.5478 |
  | ρ_surf [kg/m³] | 40.80 | 39.08 |
  | cp at the sea [J/(kg·K)] | 1349.2 | 1372.6 |
  | M_nonwater [g/mol] | 42.888 | 40.118 |
  | scale height | 15.0 km | 15.7 km |
  | dry / saturated lapse [K/km] | 7.27 / 5.04 | 7.15 / 4.84 |
  | Σf_i·κ_i (background opacity) | 1.000e-06 | **1.944e-03 (1944×)** |

  **`x_H₂O` does not move, and that is the whole shape of the change.** It is
  `p_sat(513.15 K)/p_0` — set by the sea, not by the mixture — so the 7 % the traces take comes
  out of CO₂ and N₂, which keep their dry 0.920:0.070 ratio inside what is left. The surface is
  still saturated afterwards (`q_sat` = 0.361635 against `q_H₂O` = 0.361611), which is this
  fork's defining condition and the self-test's first real check.

  **TWO COSTS, BOTH STATED IN `param.py` RATHER THAN BURIED:**
  1. **Above the cold trap the traces are 15.7 % of the dry atmosphere and CO₂ falls to 77.4 %,
     below the quoted 89–95 % band.** The stated composition for this epoch does not have room
     for ATHAD's trace inventory; adopting it is a departure from that composition, not a
     refinement of it.
  2. **The fractions match ATHAD, the amounts do not.** 1.4 mole-% of ATHAD's 250 bar column is
     3.5 bar of each gas; 1.4 % of this 60 bar column is 0.84 bar. Matching the numbers assumes
     four fifths of each trace gas left with the water. Matching the *amounts* instead would put
     them at ~5.8 mole-% each, 29 % of the column — a different atmosphere again, and the open
     question if what was wanted was one inventory across the two epochs.

  **THE PER-SPECIES OPACITY SPLIT WAS WRITTEN FOR EXACTLY THIS AND IT PAID.** Item 60's port
  was recorded here as "a verified no-op at this composition" — true while the background was
  pure N₂, and it stopped being true the moment the composition changed. `Σf_i·κ_i` is now
  1.944e-03 m²/kg, **1944× the lumped `kappa_bg` = 1e-6**, with NH₃ 48.7 %, SO₂ 36.7 % and CH₄
  13.8 % while N₂ — still 29 % of the background by mass — contributes 0.02 %. Under the lumped
  constant the model would have run five new absorbers at nitrogen's opacity and printed nothing
  to say so. *A no-op guard is worth keeping precisely because it is a no-op until it isn't.*

  **MEASURED, 40 iterations, 24 threads, ONE binary, config-only A/B** (`run_traces40` against
  `run_notraces40`, so nothing but the composition differs — see ATHAD item 61 for why the
  binary and the thread count have to be held fixed):

  | | without traces | with traces | |
  |---|---|---|---|
  | OLR | 275.36 W/m² | 273.70 | −0.60 % |
  | photosphere (τ=1) | 62.5 km, 268.23 K | **65.9 km**, 267.29 K | **+3.4 km** |
  | emission from the isothermal skin | 41.4 % of columns | **58.6 %** | +17.2 points |
  | max τ_above (surface) | 2.520e6 | 2.696e6 | +7.0 % |
  | mean planetary albedo | 0.5000 | 0.5000 | 0 |
  | `Psi_max` | 40927.25 | 40515.96 | −1.0 % |
  | max water vapour | 338.66 g/kg | 354.74 | +4.7 % |
  | max cloud water | 46.86 g/kg (6373 m) | 49.08 (7023 m) | +4.7 %, deck +650 m |
  | precipitable water | 153 564 mm | 160 626 mm | +4.6 % |
  | evaporation (Dalton) | 4426 mm/a | **6329 mm/a** | **+43 %** |
  | max surface T | 239.05 °C | 239.19 °C | +0.13 K |

  **Read the water column, not the OLR.** The OLR moves 0.60 %, which is the same size as the
  ~0.5 % ATHAD item 61 measured for a 40-iteration OLR across a bare thread-count change, and
  this fork emits **58.6 % of its columns from the isothermal skin**, so its OLR is `t_skin`
  restated even more directly than ATHAD's. What is well clear of any noise floor is the
  structure: the photosphere climbs 3.4 km, the cloud deck 650 m and the water column 4.6 %,
  all of them following the +4.4 % scale height that a lighter mixture (M_mean 29.01 → 27.79)
  buys. The **+43 % evaporation** is the largest single response and comes from the sea-surface
  saturation mass fraction rising 0.3464 → 0.3616 against an unchanged surface temperature —
  the same water, in a mixture that holds more of it per kilogram.

  **The six κ are now assumptions of the same standing as `kappa_CO2` and `kappa_H2O`**, and a
  grey scheme cannot represent the window regions that would decide them. `ATM_BG_LUMPED=1`
  restores the single `kappa_bg` and, with it, every pre-2026-08-20 background opacity.

- **The eight species have one name each now, and the `N2` collision is gone** (ported from
  ATHAD README item 61, 2026-08-20). The VTK files carried `N2` for the **Brunt-Vaisala
  frequency squared** and `q_N2` for **nitrogen**, in the same file; the array is now
  `brunt_N2` and the field `BruntVaisala_N2`. All eight species are written as bare names —
  `H2O CO2 N2 CH4 NH3 H2 CO SO2` — as mass fractions from `AtmMixture::split()`, the routine
  the thermodynamics and the radiation use, so what is plotted is what is integrated. The raw
  transported CO2 array is still written, as `CO2_tracer`; the two answer different questions
  (`CO2` is what the physics reads, `CO2_tracer` is what the transport equation advances).
  **Unverified in this fork** — the change is textual and was checked in ATHAD, but no COND run
  has been made against it.
- **`ATM_MC_GEOPOTENTIAL` is ported and UNMEASURED here** (ATHAD item 61). It adds the missing
  `g·z` to `s`, `s_u`, `s_d`, so a rising parcel cools on the dry adiabat instead of carrying
  its cloud-base temperature upward. Default off. **In ATHAD it did not restore `c_u`, and the
  reason was geometry**: the updraft there is one grid level deep — cloud base and LFS on
  adjacent levels, so the recurrence loop `for(i = i_base+1; i <= i_LFS-1)` never executes.
  **Check the base/LFS separation here BEFORE expecting the `c_u` this fork lost in item 53 to
  come back**, because if it is one level here too, the geopotential cannot be the fix.

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
- **The drag has not been scanned here.** It is now scannable (`rayleigh_kf`,
  `drag_n_layers`), and ATHAD's four-arm 200-iteration template — baseline, drag ×0.1, ×10,
  and depth 1 cell — ports directly. Do it *after* deciding the `im` question below, since
  the arms are only comparable at fixed `im`.
- **`im` = 61 here against ATHAD's 41, and the two are not freely interchangeable.** `dr` is
  derived (`cAtmosphereModel.cpp:57`, `1/(im-1)`), so the classic `dr = 0.025` coupling is
  already repaired and `im` can be changed by editing one line. But it is not a neutral
  resolution knob: it rescales the drag depth 1.6× (above), it does **not** fix the metric
  (12.29 → 11.77×), and ATHAD's own 41 was chosen for wall clock (`0ece6c7`: 1.44× faster,
  0.68× the memory) while explicitly costing accuracy at the top — after which its item 39
  found the photosphere reduced to a single grid cell of `dτ = 55` there. Two `im = 41`
  remnants are still in the tree and are cosmetic only: `Paraview_Atm.cpp:855`'s
  `dz = 0.025` (a VTK plot spacing, already wrong at 61) and `Array_1D.cpp:80`'s
  `if(mm == 41)` (suppresses one print).
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

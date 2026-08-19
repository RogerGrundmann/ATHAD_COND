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

## ATM_PRESS_SWEEPS: the elliptic solve is under-converged by ~10x, and Psi pays 28 % for it

The knob ported from ATHAD (its README item 54) answers here a question it could not answer
there. In ATHAD 99.97 % of `p_dyn` is a prescribed balanced initial state, so the solver's
contribution is swamped whatever it does. **This fork has no `initBalancedState`: `p_dyn` starts
identically zero and every bit of it is built by the flow through the Poisson solve**, which
makes the sweep count the only thing between the divergence and the pressure that removes it.

Iteration 20, same config, one binary, the knob the only difference:

| | `p_dyn` min | `p_dyn` max | radial range j=45 | latitudinal range | ratio | `Psi_max` |
|---|---|---|---|---|---|---|
| 1 sweep (shipped) | -0.00566 | 0.00854 | 0.00191 | 0.00814 | 0.235 | 40959.17 |
| 10 sweeps | -0.04454 | 0.05898 | **0.01588** | **0.08209** | 0.193 | **29436.77** |

**`p_dyn` grows by an order of magnitude and `Psi_max` falls 28 %.** The shipped one-sweep field
is not a converged pressure — it is about a tenth of one — and the meridional circulation is
being held up in part by a projection that has not finished removing the divergence. **That is a
larger effect on Psi than anything else measured in either tree this week**: the metric terms
(ATHAD item 58) net to zero at 400 iterations, the CO2 dilution moves it 3 %.

**The radial structure was never what the solver was failing to build.** The radial/latitudinal
ratio barely moves, 0.235 -> 0.193: ten sweeps build the WHOLE field roughly in proportion. So
under-convergence is an amplitude problem, not a shape problem, and item 54's radial question and
this one are genuinely separate.

**TWO THINGS THIS DOES NOT ESTABLISH, and the second is a confound in the table above.**

- **10 sweeps is not converged either.** A 50-sweep arm was launched and abandoned: see below.
  Nothing here says where the amplitude saturates, only that it is still climbing at 10.
- **THE KNOB ALSO MULTIPLIES THE INITIAL PROJECTION.** `project_initial_velocity` calls `run()`
  200 times, so `ATM_PRESS_SWEEPS=10` makes the startup 2000 relaxation passes instead of 200 —
  the 10-sweep arm therefore differs from the 1-sweep arm in its INITIAL STATE as well as in its
  per-iteration solve, and the 28 % cannot be attributed cleanly between them. Separating them
  needs a second knob. This is also why the 50-sweep arm was abandoned: 200 x 50 = **10 000**
  startup sweeps, which had not finished the projection when the run was killed. The cost of the
  knob is not linear in the count as its comment says — it is linear in the count TWICE.

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

| | at the sea | above the cold trap |
|---|---|---|
| x = (H₂O, CO₂, N₂) | .5578 / .4109 / .0313 | .010 / .920 / .070 |
| M_mean [g/mol] | 29.01 | 42.63 |
| R_mix [J/(kg·K)] | 286.6 | 195.0 |
| cp [J/(kg·K)] | 1349 | 1028 † |

† **The 1028 is the dry composition at the SEA-SURFACE temperature**, not at the cold
trap. `cp_of(dry, 513.15 K)` = 1028.1; at cold-trap conditions it is **856.8** (the
Shomate fits clamp below 298 K). This table warns two paragraphs above that the stated
composition and the stated temperature describe different heights — and then does it
itself, for this one entry. The 1349 is sound: it is `cp_of` at the sea, and it is what
`cp_l` is set to.

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

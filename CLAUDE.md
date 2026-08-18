# ATHAD_COND — the Hadean atmosphere after condensation

An atmospheric general-circulation model of the Earth in the epoch that follows
[ATHAD](../ATHAD): the magma ocean has quenched, the steam atmosphere has rained out into a
liquid ocean, and what remains above it is a **CO₂-dominated atmosphere of 27–100 bar over a
230–250 °C sea**. Nominal run point: **60 bar / 513 K**.

**Forked from `ATHAD` @ `29ca2f9` (2026-08-13), carrying its full history**, so fixes
cherry-pick in both directions (`git remote athad` points at the sibling working copy).
ATHAD in turn came from `ATOM_Precipitation` @ `1e3f319`. There is still no paleogeography
and no time-slice series — ATHAD_COND is one epoch. There is now an **ocean**, but as a
boundary condition only: no prognostic hydrosphere.

**The defining inversion.** ATHAD's central fact is that water is supercritical from the
ground to ~177 km, so nothing condenses and half the code exists to make condensation a
genuine no-op. Here water is subcritical everywhere and **condensation is live from the sea
surface up**. Every path ATHAD made inert is one this model depends on. See "Four
invariants" below — invariant 2 is inverted, the other three are not.

## Build, run, test

```bash
make cond                                       # -> cli/cond
make test                                       # IAPWS + column self-tests, run these first
cd python && OMP_NUM_THREADS=8 ../cli/cond config_cond.xml
```

**`atom_log.txt` is reserved by the model — never redirect a run into it.**
`lib/Utils.cpp:45` opens `atom_log.txt` in the run directory with `std::ofstream::out`,
which truncates it. A shell redirect to the same name loses the whole console output when
the model's handle closes; a 100-iteration run's printouts were destroyed this way. Use any
other filename.

Output lands in `python/output_Hadean_condensation/` — the ATOM line's convention
(`output_<name>/`, underscore, relative to the run directory; the giant-planet
siblings hyphenate instead).

`make` regenerates the parameter bindings from `param.py` first. The generated files
(`atmosphere/*.inc`, `python/atmosphere_pxd.pxi`, `python/pycond.pyx`,
`cli/config_cond.xml`, `python/config_cond.xml`) are **tracked on purpose**, so a
`param.py` change shows its full effect in the diff. Regenerate and commit them together.

Removing a parameter from `param.py` deletes a C++ member, so it must be removed
together with its uses or the build breaks.

## Atmospheric composition

Mole fractions are the input; the model works in mass fractions. **The configured mole
fractions are the SEA-SURFACE values, not the quoted dry composition** — see below.

| | at the sea (config) | above the cold trap |
|---|---|---|
| x(H₂O) / x(CO₂) / x(N₂) | 0.5578 / 0.4109 / 0.0313 | 0.010 / 0.920 / 0.070 |
| M_mean [g/mol] | 29.01 | 42.63 |
| R_mix [J/(kg·K)] | 286.6 | 195.0 |
| cp [J/(kg·K)] | 1349 | 1028 |
| q_H₂O / q_CO₂ [kg/kg] | 0.3464 / 0.6233 | 0.00423 / 0.9498 |

- **p_surf = 60 bar**, **T_surf = 513.15 K** (prescribed), **ρ_surf = 40.80 kg/m³**
- **scale height 15.0 km** at the sea (ATHAD: 59.3 km) — the shell is 120 km, not 300
- **Background** (everything except H₂O and CO₂) is N₂ alone: M_bg = 28.014 g/mol,
  **R_bg = 296.8 J/(kg·K)**. This is what `R_Air` means. The quantity that matters for
  saturation is **`M_nonwater` = 42.88 g/mol** (CO₂ + background) — a different number and
  a different function. Passing `M_bg` where `M_nonwater` belongs puts the sea surface at
  `q_sat` = 0.448 against the true 0.346.
- ATHAD's five reducing trace gases (CH₄/NH₃/H₂/CO/SO₂) are zero here.

**Why the config does not carry the quoted composition.** The literature figure — H₂O
0.4–2 %, CO₂ 89–95 %, N₂ 5–20 % — is the *dry* atmosphere. Air over a 240 °C sea holds
water at its own vapour pressure: `p_sat(513.15 K) = 33.47 bar` of a 60 bar column, so
`x_H₂O = 0.5578` at the surface. Both figures are true at different heights.
`test/cond_column_selftest.cpp` asserts the arithmetic connecting them; `r_air` and the
`x_*` are derived and guarded there.

Only **H₂O (`c`) and CO₂ (`co2`) are prognostic**, both as mass fractions. CO₂ is well
mixed **in the dry air** — `q_CO2 = (1 − c)·f_CO2` — not as a uniform mass fraction, which
is what ATHAD does and is only equivalent when the water field is uniform too.

## Where the physics lives

| File | What it owns |
|---|---|
| `MixtureAtm.h` | Composition → mass fractions, `R_of`, `cp_of` (Shomate), `M_of`, `M_nonwater`, water's critical point |
| `SaturationH2O.h` | IAPWS saturation + sublimation, Watson `latentHeat`, exact `saturationMassFraction`, `dewPoint`, `satDerivFactor`, exact `dqSatdT`, **`moistLapse`** |
| `MultiLayerRadiation.h` | Grey optical depth from column mass with pressure broadening; surface energy balance |
| `ThermoAtm.h` | The column integration (`densities`), the sea saturation boundary (`waterVapourEvaporation`), diagnostics |
| `test/saturation_selftest.cpp` | IAPWS reference points |
| `test/cond_column_selftest.cpp` | The ATHAD_COND regime: both mixtures, the config's own arithmetic, the moist lapse (with an Earth control) |

## Four invariants — do not silently break these

1. **There is no topography, and the planet is hemispherically symmetric.** `h ≡ 0`,
   `i_topography ≡ 0`, `is_land()` false everywhere — now a global *ocean*, the same
   statement in code. `LandOceanFraction()` throws if a land point appears. Any
   north–south asymmetry is a defect: symmetric surface parabola, mirrored insolation, no
   obliquity, no seasons.

2. **Water is subcritical everywhere and condensation is LIVE.** *This is ATHAD's
   invariant 2, inverted, and it is why this fork exists.* Every path ATHAD made a genuine
   no-op is a path that now runs. What carries over is which formulas are allowed: **never
   Magnus** (calibrated to ~320 K, returns 1.2e7 hPa at 513 K) and **never the dilute
   `q_sat = ep·E/(p−E)`** — invalid here for the opposite reason it was invalid in ATHAD:
   water is 0.4 % by mole aloft, where the dilute form would be fine, and 56 % at the sea,
   where it is not. A form valid over part of the column is worse than one valid over
   none, because it looks right in the printouts.

3. **Radiation runs in mode 2** (direct σT⁴). Modes 0/1/3/4/5 lean on the Scotese snapshot
   or the 280 ppm CO₂ reference; neither exists in a 92 % CO₂ atmosphere at 4.4 Ga.
   Radiation must *set* the profile rather than be handed one — **and it does not yet**:
   `densities()` re-imposes the adiabat every iteration, and the OLR currently equals
   `σT_lid⁴` exactly (README item 6).

4. **The column is on its own adiabat, integrated not fitted, and the adiabat is MOIST
   where the air is saturated.** `SaturationH2O::moistLapse`, hydrostatic on the layer-mean
   T, isothermal above `t_skin`. Do **not** restore the COSMO `T = T₀√(1−coeff·h)` form.
   Do **not** revert to `dT/dz = −g/cp`: `cosmo_lapse_fraction = 1.0` is justified in ATHAD
   *because nothing condenses there*, and that justification is void here.
   And do not use `g/(cp + L·dq_sat/dT)` — that is only the denominator of the moist-static-
   energy balance, and it gives 0.7 K/km here against the correct 5.04.

## Assumptions vs. results — read this before quoting any number

| Parameter | Value | Status |
|---|---|---|
| `geothermal_flux` | 150 W/m² | **Inherited from ATHAD's magma ocean and almost certainly wrong.** 55 % of the whole energy input, against modern Earth's 0.087 W/m². Biggest unexamined input |
| `kappa_CO2` / `kappa_H2O` / `kappa_bg` | 0.001 / 0.01 / 1e-6 m²/kg | Factor-of-2 uncertain. `kappa_CO2` now dominates (CO₂ column 5.8e5 kg/m² vs water 2.6e3) |
| `t_surf_equator` / `t_surf_pole` | 513.15 / 503.15 K | **Prescribed, not solved** |
| `t_skin` | fixed-point iterate | Currently *sets* the OLR rather than following it |
| `albedo_cloud` | 0.50 | IS the planetary albedo (0.4997 measured); saturates wherever condensate exists, which is now everywhere |
| `omega` | 3.17e-4 (5.5 h day) | Inherited; this epoch is later and slower |
| `delta_i_c` | 500 s | Bechtold convective timescale, Earth-calibrated. A time, not a pressure — no unit error, and no evidence here to replace it |

## What the model currently says

Measured; see the README for the full items.

- **Saturated troposphere 62 km deep.** 511 K at the sea, 429 K at 18 km, 280 K at 60 km;
  first 10 km at **5.02 K/km** against the 5.04 the moist lapse predicts. Precipitable
  water 133.5 m. `q_H₂O` tracks `q_sat` to four digits.
- **The cold trap is 8× too wet**: 0.0349 kg/kg against the 0.00423 the stated dry
  composition implies. Reconciling them needs a cold trap near **217 K**. Testable.
- **The OLR is not yet an output**: OLR = σT_lid⁴ = 268.10 W/m² to six figures, with
  eps = 0.0000 and a transparent lid. τ ≈ 1 falls at ~0.07 bar, which is where the
  isothermal skin begins, so the emission comes from a prescribed temperature and `t_skin`
  relaxes to close the balance it is reporting. Deepening the shell will not fix it.
- Not yet run beyond 20 iterations: no stability run, no grid convergence, no thread
  determinism check.

## Relationship to the family

Siblings live beside this directory: `ATHAD` (the parent, ~4.4 Ga magma ocean),
`ATOM_Precipitation` (modern Earth), `ATJUP`, `ATSAT`, `ATURAN`, `ATNEPT` (giants),
`ASTIM` (impacts).

C++ class, file and function names are kept **identical to ATHAD**, which keeps them
identical to `ATOM_Precipitation`; only the outer shell is renamed (`libcond.a`,
`cli/cond`, `config_cond.xml`, `pycond`). Preserve that. The git history is ATHAD's, and
`git remote athad` points at the sibling working copy, so `git cherry-pick` works in both
directions — use it rather than re-implementing a fix that already exists next door.

**The pattern to expect, and why it is worse here than in ATHAD.** ATHAD found seventeen
defects of one shape: *Earth's numbers as bare literals inside physics kernels.* At 250 bar
those are wrong by orders of magnitude and obvious. **At 60 bar they are wrong by a factor,
and plausible.** `p_stat <= 1000 hPa` fails at every level of ATHAD's column and so does
nothing; here it selects everything above 53 km and fires the deep-convection scheme in the
stratosphere. When something behaves oddly, look for a constant that was true at 1 bar and
288 K — and do not assume that a branch which now *runs* is therefore correct.

**⚠️ DO NOT CHERRY-PICK ATHAD'S BALANCED INITIAL STATE AS IT WAS FIRST WRITTEN.** ATHAD
items 26-28 build a `initBalancedState` that this model does not have. Item 27's version
balances the **θ-momentum equation only**, and porting it here would reproduce a failure
ATHAD measured in full: it drove `max |u|` from 0.114 to **11.17 m/s** over 200 iterations,
reversed the tropics to sinking over the hottest surface on the planet, and the meridional
streamfunction — **which is built from `v` alone and therefore cannot see a radial failure**
— reported the circulation as healthy throughout. The streamfunction here was just repaired
(`a85122c`, the density inside the integral); that fixes what it measures, not what it is
blind to.

The reason it fails is a **switch**, and this model has the identical one. At `u = v = 0`:

| | ATHAD | here |
|---|---|---|
| `rhs_u` | `RHS_Atm_Turb.cpp:1001` | line 1001, identical |
| `AtomUtils::coriolis_nontraditional()` (gates `coriolis_rad`) | default **false** | `lib/Utils.h:52`, default **false** |
| `AtomUtils::metric_curvature()` (gates the `−(v²+w²)/r` term) | default **false** | `lib/Utils.h:90`, default **false** |
| `buoyancy_ramp` at iteration 0 | **0** (`buoyancy_ramp_iters` = 300) | **0**, same 300 |
| `ATM_METRIC_RADIUS` | default on | default on |

So the radial momentum equation at iteration 0 is `rhs_u = −dp_dyn/dr·exp_rm` **and nothing
else**, and any radial gradient a balance writes is a pure unopposed vertical force. The
same switch governs the `w²cotθ` term item 27's θ-balance is *built on*, so with the
defaults that balance is also computed against a force the model does not apply — though
that part is worth only 0.3 % of `F_θ` (both balances are ~99.7 % Coriolis).

**And turning the switch on does not fix it**, which was worth testing before assuming:
with `ATOM_METRIC_CURVATURE=1` the radial equation gains a force of 0.075 rms against the
~293 a θ-only balance creates — **~4000× too small**. The θ-force's shape varies with height
because `w` does, and no radial force of comparable size exists to pay for that variation.
Completing the metric does not rescue a one-component balance; only balancing both does.

**The lesson generalises past the balance, and is the reason this note is here rather than
in a commit message: a term written in `RHS_Atm_Turb.cpp` is not necessarily a term the
model applies.** Read the switch, not just the source. This is the family's constant trap
one level up — not an Earth number in a kernel, an Earth-tested code path that never runs.

If the balance is ported later, port **ATHAD item 28's two-component version** (it minimises
the unbalanced acceleration in both equations at once and reads `F_r`/`F_θ` switch by switch,
so it follows whatever this model is configured to do) **together with its residual
diagnostic**, which prints what each component is left holding. There is no version of "port
item 27 now and fix it after". Whether a balanced initial state helps here **at all** is a
separate question needing its own A/B in this regime — 120 km shell, saturated sea — and
ATHAD's own 200-iteration confirmation was still running when this was written (2026-08-15).

Fixes made here that belong upstream in ATHAD: the exact `dqSatdT` (ATHAD's dilute form is
correct only in the `x → 0` limit, and ATHAD's water is not dilute either — it is the bulk
gas, so the error there is larger, not smaller); `moistLapse`; and the `M_bg`/`M_nonwater`
distinction at the saturation call sites.

Traps already solved elsewhere in the family — check before re-deriving:
Coriolis/centrifugal signs (ATURAN `8b284cb`, `4201957`; ATNEPT `024c37f`, `e412b1b`);
mass- not mole-weighted mixture properties (ATNEPT `c116d71`); in-place Gauss–Seidel as a
threading defect (ATURAN `ffd0e0e`); report failures and limits in the README (ATURAN
`74b4ded`, ATNEPT `34286b8`).

## Open risks

- **The dynamics and the radiation do not agree where the levels are — measured here now,
  not inherited as an estimate** (ported from ATHAD README item 39, 2026-08-18).
  `exp_rm = 1/(rm+1)` is documented in two places as the Jacobian of the radial coordinate
  transformation (`TurbulenceAtm.h`; `PressureSolverAtm.h` writes it as
  `dp/dr_physical = exp_rm * dp/d(rad.z)`). It is the Jacobian of a **quadratic** stretch,
  `z ~ (rm+1)²/2`, while `init_layer_heights` builds an **exponential** one. `checkRadialMetric()`
  now prints the spread at every startup, and it is **12.29×** here: the core's radial length
  unit is 9 837 m at the surface and 120 850 m at the top, while `get_layer_height()` — which
  the radiation, `ConvectiveAdjustment` and every diagnostic use — has the true heights.
  CLAUDE.md previously carried "~12×" as an inference from `zeta = 3.0`; that inference was
  right, and it is now a measurement. The test is unit-free
  (`ratio(i) = [inv_2dr·exp_rm]/[1/(z[i+1]−z[i−1])]` must be constant in `i`), so it cannot be
  argued away as a units convention. **`im` is not the lever**: the spread depends only on
  `zeta` and `im`, and moving `im` 61 → 41 changes it 12.29 → 11.77. `zeta = ln(1.5) = 0.405`
  would make `exp_rm` correct by construction. `ATM_METRIC_CHECK=1` adds the per-level table.
  The diagnostic is print-only and self-silencing below a spread of 1.05.
- **The surface drag is Earth's, twice over, and its depth moves with the grid**
  (ported from ATHAD README item 44, 2026-08-18). `rayleigh_kf` and `drag_n_layers` were
  `constexpr` in `RHS_Atm_Turb.cpp` with comments justifying them by Earth's geography — a jet
  off the west coast of South America for the strength, an Andes orographic feature for the
  depth. Invariant 1 makes `is_land()` false everywhere here too. **`drag_n_layers = 5` is a
  count of AIR CELLS, not a length**, so its physical depth is whatever the grid makes it:
  236 m on `ATOM_Precipitation`'s grid, **1 786 m here**, 7 152 m in ATHAD. Both are config
  parameters now, bit-identically. **The part ATHAD's writeup did not have to state**: the
  depth moves with `im` as well as `zeta`, and ATHAD holds `im` fixed at 41 while this fork is
  at 61 — the same constant means 2 861 m at `im` 41, 1 786 m at 61 and 1 297 m at 81, a 2.2×
  swing from the level count alone. **A grid-refinement study here silently rescales the drag**,
  so `im` and `drag_n_layers` cannot be varied independently until the depth is reformulated as
  a length in metres. Not yet scanned here; ATHAD's four-arm scan is the template.
- **`Psi`, `tau_above`, `tau_layer` and `N2` are fields now, and `N2` VALIDATES the moist
  adiabat** (ported from ATHAD's item 42, 2026-08-18). All four were absent: `Psi` existed only
  as a CSV of zonal means plus the scalar `Psi_max`, and the other three were computed **nowhere
  in this model**. Measured at 20 iterations:

  | field | value |
  |---|---|
  | `tau_above` | 0 at the lid → **2.55e6** at the surface |
  | `tau_layer` | up to **1.11e5** per layer, at 7.7 km |
  | `N2` | **4.30e-5** 1/s² at the surface, 4.25e-4 at 114 km |
  | `Psi` | 4.2212365e13 kg/s, symmetric about the equator |

  **`N2` is the result worth having.** A column on the MOIST adiabat is DRY-stable, so the dry
  N² must be positive and equal to `(g/T)(Γ_dry − Γ_moist)`. With the self-test's own lapse
  rates that predicts **4.270e-5**, against a measured **4.30e-5** — agreement to **0.7 %**.
  So invariant 4 in its moist form is now *measured* rather than asserted, and the saturated
  sweep in `densities()` is integrating the lapse rate it claims to. Note this differs from
  ATHAD by design: there N² ≈ 0 through the column because the adiabat is dry.
  **`Psi` also caught a defect** — it read identically zero until the `write_meridional_
  streamfunction` / `print_min_max_atm` order was swapped (ATHAD's item 2, which had not been
  ported with the rest). The CSV was always right; the printed extrema and the field were not.
  **Not ported deliberately**: the `.vts` panorama carries none of these in either fork.
  **The photosphere is measured now too** (block ported 2026-08-18, 20 iterations, 24
  threads): **59.7 km, T there = 278.6 K, and 0.0 % of columns emitting from within 1 K of
  `t_skin`.** Against ATHAD's 239.4 km / 325.2 K with a skin fraction rising 1.1 → 39.2 %.
  **That difference matters more than the numbers.** ATHAD's `t_skin` pinning works by the
  radiating level migrating into the isothermal top, and here it has not: the photosphere sits
  at half the shell depth, 16 K above `t_skin`, in **zero** columns of the skin. So this fork's
  "the OLR is not independent of `t_skin`" risk cannot be the *same* mechanism as ATHAD's, and
  the OLR is correspondingly not sitting on it — 283.81 W/m² against σT_lid⁴ = 270.62 and
  absorbed = 270.78, a **−13.03 W/m² imbalance that is not identically the distance to a
  constant**. Whether it closes onto one is now checkable rather than assumable, which is what
  the instrument was for. Note the lid emissivity is still 0.0000, so the *lid* is transparent
  while the photosphere is deep — the two are decoupled here, unlike ATHAD at 230 km.
- **The radial momentum balance is the pressure gradient and nothing else, and the buoyancy
  body force is absent** (`ubud_*` ported from ATHAD's item 42, 2026-08-18, measured here at
  20 iterations). Ψ is built from `v` alone, so it cannot see a radial failure; `vbud_*`/`wbud_*`
  covered two components of three and this closes the set.

  | term | value | vs `ubud_pgf` |
  |---|---|---|
  | `ubud_pgf` | 0.100816 | — |
  | `ubud_cor` | **0.000000** | exactly zero |
  | `ubud_advv` | 0.002092 | 2.1 % |
  | `ubud_advh` | 0.002844 | 2.8 % |
  | `ubud_diff` | 0.000445 | 0.44 % |
  | `ubud_buoy` | **0.000000** | < 5e-7, so >2e5 down |

  Same structure as ATHAD's item 42 (pgf 1.15, cor exactly 0, advh 0.0085, buoy 1e-6) at ~10×
  smaller magnitude. Two of these are self-checks that now pass rather than being asserted:
  **`ubud_cor` is identically zero**, so `coriolis_nontraditional()` is genuinely off and not
  merely documented off — the CLAUDE.md warning that "a term written in `RHS_Atm_Turb.cpp` is
  not necessarily a term the model applies" is now instrumented rather than trusted. And
  **`ubud_buoy` is effectively zero**, which measures ATHAD item 34's *other* extra-`*dt` term —
  the half **not** repaired by the `surf_drag` fix in `a609d2b`. The buoyancy body force in a
  model whose density spans ~3 orders of magnitude is currently not acting.
- **The OLR is not independent of `t_skin`.** First thing to fix, and the same task as
  invariant 3: the profile is prescribed, so radiation cannot set it.
- **`geothermal_flux` is over half the energy budget** and is ATHAD's magma-ocean number.
  Nothing quantitative survives it being wrong.
- **The surface temperature is prescribed, not solved.** Every result is conditional on it.
- **Rain has nowhere to go.** The sea is a boundary condition, not a reservoir, so the
  water budget is not a closed-system test: `moist_phys_start_iter` is 0 here, not ATHAD’s 300, because condensation is the subject.
- **Boussinesq.** `ATM_ANELASTIC` is ported from ATHAD and ships default-off; it has not
  been measured in this regime. The density span here is ~3 orders of magnitude instead of
  ~5, which makes this the better testbed for whether it matters at all.
- **The column air mass is not conserved** (inherited): `p_stat.x[0]` is re-anchored every
  iteration to `r_air·R_mix·T_surf`.
- `time_start/end/step` remain because the time-slice loop is still structural, though only
  one slice ever runs.

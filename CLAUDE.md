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

**The five trace gases carry ATHAD's mole fractions since 2026-08-20** (they were zero
before). x(H₂O) is not free — it is `p_sat(T_surf)/p_0` — so they displace CO₂ and N₂,
and above the cold trap they are 15.7 % of the dry atmosphere, which puts CO₂ at 77.4 %,
**below the quoted 89–95 % band**. That is a deliberate departure, not a refinement.

| | at the sea (config) | above the cold trap |
|---|---|---|
| x(H₂O) / x(CO₂) / x(N₂) | 0.5578 / 0.3459 / 0.0263 | 0.010 / 0.7744 / 0.0589 |
| x(CH₄) = x(NH₃) = x(H₂) = x(CO) = x(SO₂) | 0.014 each (7 % total) | 0.0313 each (15.7 %) |
| M_mean [g/mol] | 27.79 | 39.90 |
| R_mix [J/(kg·K)] | 299.2 | 208.4 |
| cp [J/(kg·K)] | 1372.6 | 1042.5 † |

† **The 1042.5 is the dry composition at the SEA-SURFACE temperature**, not at the cold
trap — `cp_of(dry, 513.15 K)`; at cold-trap conditions it is **886.4** (the Shomate fits
clamp below 298 K). This table warns two paragraphs above that the stated composition and
the stated temperature describe different heights — and then does it itself, for this one
entry. The 1372.6 is sound: it is `cp_of` at the sea, and it is what `cp_l` is set to.

| q_H₂O / q_CO₂ [kg/kg] | 0.3616 / 0.5478 | 0.00452 / 0.8542 |

- **p_surf = 60 bar**, **T_surf = 513.15 K** (prescribed), **ρ_surf = 39.08 kg/m³**
- **scale height 15.7 km** at the sea (ATHAD: 59.3 km) — the shell is 120 km, not 300
- **Background** (everything except H₂O and CO₂) is N₂ **plus the five traces**:
  M_bg = 26.138 g/mol, **R_bg = 318.10 J/(kg·K)**, within 0.3 % of ATHAD's 317.26 because
  in both forks the traces are ~71 % of the background by mass. This is what `R_Air`
  means. The quantity that matters for saturation is **`M_nonwater` = 40.12 g/mol**
  (CO₂ + background) — a different number and a different function. Passing `M_bg` where
  `M_nonwater` belongs puts the sea surface at the wrong `q_sat`.
- **The background's grey opacity went 1e-6 → 1.944e-3 m²/kg with them** (Σ f_i·κ_i), so
  item 60's per-species split, which was an exact no-op at the old composition, is now
  carrying NH₃ (48.7 %), SO₂ (36.7 %) and CH₄ (13.8 %). Had the composition changed with
  the lumped κ still in place, five new absorbers would have run at nitrogen's opacity.
- **The trace INVENTORY is not conserved from ATHAD**, only the fractions: 1.4 % of 250 bar
  is 3.5 bar per gas, 1.4 % of 60 bar is 0.84 bar. Matching the numbers assumes four fifths
  of each trace gas left with the water.

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
| `kappa_CO2` / `kappa_H2O` | 0.001 / 0.01 m²/kg | Factor-of-2 uncertain. `kappa_CO2` dominates (CO₂ column 5.8e5 kg/m² vs water 2.6e3) |
| `kappa_N2/CH4/NH3/H2/CO/SO2` | 1e-6 / 3e-3 / 1e-2 / 1e-5 / 1e-4 / 2e-3 m²/kg | ATHAD's values, and they now matter here: Σf_i·κ_i = **1.944e-3**, 1944× the lumped `kappa_bg` = 1e-6 that was exact while the background was pure N₂. Six assumptions of the same standing as the two above |
| trace-gas mole fractions | 0.014 each (CH₄/NH₃/H₂/CO/SO₂) | ATHAD's numbers, adopted 2026-08-20. **The fractions match, the inventory does not** (1.4 % of 250 bar ≠ 1.4 % of 60 bar), and they push dry CO₂ to 77.4 %, below the quoted 89–95 % |
| `t_surf_equator` / `t_surf_pole` | 513.15 / 503.15 K | **Prescribed, not solved** |
| `t_skin` | fixed-point iterate | Currently *sets* the OLR rather than following it. **Its fixed point is `sigma*T^4 = F/2` since 2026-08-21** (ATHAD item 67); it used to be `= F`, which made the budget close by construction. `ATM_SKIN_GREY=0` restores that. Effect size here UNMEASURED |
| `albedo_cloud` | 0.50 | IS the planetary albedo (0.4997 measured); saturates wherever condensate exists, which is now everywhere |
| `albedo_surface` | **0.06** | A **water** value since 2026-08-20; it was 0.08 with ATHAD's basaltic-melt justification attached. Clear-sky only, and the cloud bump overwrites it almost everywhere, so its reach is small — measured, not assumed |
| `omega` | 3.17e-4 (5.5 h day) | Inherited; this epoch is later and slower |
| `mc_M_max` | **100 kg/(m²s)** | Sized 2026-08-20 from the density ratio (Earth's 3.0 implies σw = 2.5 m/s; ×36–42 kg/m³ here) and confirmed against the uncapped scheme's own 13.0–13.7. Was Earth's 3.0 and **pinned from 2.2 to 21.9 km in 177 of 181 columns** |
| `delta_i_c` | 500 s | Bechtold convective timescale, Earth-calibrated. A time, not a pressure — no unit error, and no evidence here to replace it |

## What the model currently says

> **THE GREY SKIN FACTOR IS IN, AND IT IS ON BY DEFAULT (ported from ATHAD item 67, 2026-08-21).**
> `t_skin` was solving `sigma*T^4 = F` where `F` = absorbed SW + geothermal. That is the planet's
> EFFECTIVE EMISSION temperature; the top layer of a grey atmosphere is at `T_eff/2^(1/4)`, because
> it sees no downward flux and re-emits half of what passes up through it. **The precondition was
> checked in THIS tree, not assumed**: `SW_abs` appears in exactly one place, `MultiLayerRadiation`'s
> surface energy balance, and nowhere in the column, so the classical skin result `sigma*T^4 = F/2`
> is exact here rather than a fit. `MultiLayerRadiation` already contained the right relation twice.
>
> **Consequence: `t_skin` is now ~2^(1/4) colder, and the OLR follows it down.** Because the reported
> OLR descends onto `sigma*t_skin^4`, setting `sigma*t_skin^4 = F` made "OLR -> absorbed" an
> identity. In ATHAD, correcting it turned a -1.41 W/m2 imbalance at 200 iterations into **+123.76**.
>
> **EVERY OLR AND IMBALANCE FIGURE RECORDED BELOW PREDATES THIS AND WAS MEASURED WITH THE FACTOR
> OFF.** None is retracted -- they are correct measurements of a branch still reachable with
> `ATM_SKIN_GREY=0` -- but none describes the shipped model.
>
> **THE SIZE OF THE EFFECT IN THIS TREE IS NOT MEASURED.** ATHAD is 250 bar and this is a different
> column; the defect and its repair transfer, the magnitude does not. Run the A/B before quoting a
> number here.


> **THE 2026-08-24 PORT FROM ATHAD (its README items 68-80). FOUR DEFAULTS CHANGED, AND THE
> A/B IS RUN.** Every number recorded in this file and the README predates the port; they are
> not retracted -- each of the four is restorable with one environment variable -- but none of
> them describes the shipped model. **Measured, 40 iterations, 24 threads, env-only arms: OLR
> 172.69 -> 196.44 W/m2 (+13.8 %), `Psi_max` -28 %, `div(rho u)/rho` rms +56 %, `P_snow` 0 ->
> 9.44e4 mm/a.** Three things the README's port entry says at length and this table cannot:
> **`ATM_SAT_SUPERHEAT` is INERT here** (reached 3.4 million times, rejects zero, and zero
> forbidden cells at all four census stages -- ATHAD's largest effect is a no-op in a
> subcritical column); **the precipitation flux is the `P_max_flux` cap by a factor of 9.3
> MILLION**, so no precip number here means anything yet; and **`Psi_max` is still at z = 0 m
> in both arms**, so ATHAD's "the maximum has moved off the ground" must NOT be imported.
>
> | knob | was | is | ATHAD's finding |
> |---|---|---|---|
> | `ATM_PROJ_SWEEPS` | 1 | **10** | the initial pressure projection was under-converged: `Psi(ground)`, which must be 0, was **2.09x the interior circulation**, and 10 sweeps take 52.5 % of that off. A knee, not a cure (item 68) |
> | `ATM_RAD_DIRECT` | absent | **on** | the closed-form two-stream equilibrium solve, exact in two O(N) passes, replaces the 4-sweep Lambda iteration. `n_lambda = 4` is an Earth constant: **15 % high** on ATHAD's column, and the exact answer costs nothing (items 30, 71). `ATM_N_LAMBDA` is new alongside it |
> | `ATM_SAT_SUPERHEAT` | absent | **on** | `clampAndFade` tested whether a cell could hold a condensed phase, condensed, added the latent heat that makes the answer FALSE, and never re-tested. The ice scheme deleted the result every iteration. **The largest single effect ATHAD has recorded: OLR -49.4 %** (item 75) |
> | `ATM_PRECIP_BANDS` | absent | **on** | the precipitation phase bands were Earth's, and they did two jobs at once: `(band) ? (inherited + produced) : 0` DESTROYS a flux arriving from above rather than only withholding production. Snow's -20 C floor goes; the inherited flux always passes (item 76) |
>
> **Why they ship on rather than off.** Each is a repair whose ARGUMENT is universal -- a
> projection that has not converged, a relaxation that has not converged, a routine disagreeing
> with `IceSchemeCommon::canCondense` about where a phase can exist, a flux destroyed by the
> temperature of the air it falls through. **The SIZES were ATHAD's and did not transfer, which
> is now measured rather than warned about**: one of the four turned out to be inert here and
> another turned out to be capped. Restore the old branch with
> `ATM_PROJ_SWEEPS=1 ATM_RAD_DIRECT=0 ATM_SAT_SUPERHEAT=0 ATM_PRECIP_BANDS=0`.
> **THE ONE-KNOB DECOMPOSITION IS RUN (2026-08-24) and every effect is single-valued**:
> the OLR is `ATM_PRECIP_BANDS` (+14.4 %), `ATM_RAD_DIRECT` is +1.6 %, `ATM_PROJ_SWEEPS` and
> `ATM_SAT_SUPERHEAT` are 0.00; the divergence rise is entirely `ATM_PROJ_SWEEPS`. And the
> bands effect is TWO INSTRUMENTS AGREEING, not a rise -- the equator-only OLR is 196.6 in
> every arm while the global mean goes 172.69 -> 197.56, i.e. the extratropics were radiating
> 12 % below the tropics. See the README.
>
> **Default-off knobs ported with it, off-branch bit-identical:** `ATM_CELL_ALTERNATE` (item 69 --
> four of five prescribed cells turn the same way, and `edgeRadialCoeff` already assumes the
> parity it does not impose), `ATM_METRIC_EXACT` (item 80 -- `exp_rm = 1/(rm+1)` is a QUADRATIC
> stretch's Jacobian applied to an exponential grid; the repair also adds the `-(J'/J)f'`
> curvature term the Laplacian omits in both metrics. **Default off because in ATHAD the correct
> Jacobian made the projection WORSE**), `ATM_PRECIP_CAP` (item 77). Print-only: the `ATM_ICE_CENSUS`
> stages, the `canCondense`-violation counter, the `P_rain` cap probe and the `S_r` decomposition
> (items 74, 77, 78), and the zonal ParaView writer's **true-height vertical axis** (item 70 --
> it was level index, which on this stretched grid distorts by ~18x and varies with altitude, so
> no aspect setting could undo it; no physics reads those fields).
>
> **Deliberately NOT ported.** `ATM_SKIN_TAU` (item 73): measured in ATHAD and **withdrawn** --
> `T_rad(tau)` exceeds the prescribed profile at all 41 levels there, so there is no
> radiative-convective crossing to switch at. It also presumes ATHAD's isothermal lid over a
> prescribed dry adiabat and its `ATM_PROGNOSTIC_T` knob, and this tree has neither: its adiabat
> is moist and state-dependent. `moist_phys_start_iter = 0` (item 74) is already the default here
> and has been since the fork.

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

**THE BUILD CAN SHIP A STALE OBJECT, AND THE `-MMD` REPAIR DOES NOT COVER IT** (ported
2026-08-26 from ATOM_Precipitation, where it actually bit). `-MMD` tracks headers only for
objects ALREADY compiled with it; one built before the flag has no `.d`, so make sees only its
`.cpp` and never rebuilds it. Upstream's `cli/atm.o` sat four weeks stale through a large port,
and because it holds the model as a **stack local**, adding members to the class made `main`
reserve the OLD `sizeof` while the library constructed the NEW one -- the constructor ran off
the end of main's frame onto the stack canary. Every run aborted with
`*** stack smashing detected ***` AFTER completing successfully, with no compiler warning, and
**AddressSanitizer saw nothing** (an ODR size mismatch is not an out-of-bounds access to
anything it tracks). Found by disassembling `main`. The cure is one line -- **every object
depends on the Makefile itself** -- and it is PROPHYLACTIC here: every C++ object in this tree
already had its `.d`, checked before applying.

**The constant-density streamfunction claim in ATHAD's CLAUDE.md was already out of date about
this tree** -- the density-weighted version is present here, and ATOM_Precipitation got it on
2026-08-26 (2.28x overweight at its cell core, `Psi_max` 851.68 -> 373.84 (1e9 kg/s), location
unmoved). Nothing to do here; recorded so it is not re-derived.

**The longitudinal ParaView slice was on level index here too**, fixed the same day: item 70
put the ZONAL writer on true height and fixed only that one, so every longitude-height figure
kept the distortion. Plot-only; any conclusion drawn by eye from an older longitudinal plot
should be re-examined.

## Open risks

- **THE CONVECTION SCHEME IS A MOIST CONVECTIVE ADJUSTMENT SINCE 2026-10-10; TIEDTKE IS NOT CALLED**
  (README, first section; the user's instruction). `ATM_CONV_SCHEME=tiedtke` restores it, and
  **every number recorded before that date, and every `ATM_MC_*` knob and instrument below, needs
  it** — with the default they act on a routine that does not run. **In this tree `adjust` adds
  no routine**: `densities()` already rebuilds `t` on the moist adiabat every moist step
  (invariant 4) and `SaturationAdjustment` + the ice scheme condense and rain; `adjust` stops
  calling `MoistConvection` and zeroes `MC_t`, `MC_q`, `MC_v`, `MC_w`, `P_conv` and the updraft /
  downdraft fields (`convectionStep()` in `cAtmosphereModel.cpp`). Measured, 20 iterations,
  8 threads: 120 of 180 printed extrema identical, the scheme's own fields go to zero, **everything
  else moves in the 6th-8th digit** (max T 238.881236 -> 238.881216 C, precipitation average
  199 163.86 -> 199 163.62 mm/a). **The long run, 400 iterations, three arms
  (`run_long.sh`; two Tiedtke runs byte-identical, so the noise floor is zero)**: 104 of 180
  extrema identical, the rest below 0.03 %; fields differ by at most 0.0023 K and 0.32 g/kg of
  vapour (1e-3 of its range) — **but that difference grows linearly, 0.08 g/kg per 100
  iterations, and has not levelled off**, so this says nothing about a spun-up state. ATHAD and ATHAD_PERID still call
  Tiedtke.
- **THREE PRINT-ONLY INSTRUMENTS FROM UPSTREAM, 2026-10-10: `ATM_SATADJ_DIAG`, `ATM_MC_DIAG`,
  `ATM_MC_CAP_DIAG`** (README, first section). Default off; byte-identical off and on (20/20
  files, 4 iterations, 1 thread). **Use them before and after any moist-physics change here.**
  Their first readings, 20 iterations, 8 threads, none followed up: **`cloud_cap` deletes
  ~19-22 mm of column water per moist call** (−202 mm in eleven calls, 0.13 % of the
  precipitable water, 99 % of what `SaturationAdjustment` removes); the cold fade −0.29 mm per
  call; the phase split creates +0.005 mm per call (`ATM_SAT_PHASE_AVAIL`'s defect, sized);
  **convective rain generation `g_p` is exactly zero while 100 % of columns are active and
  deep**; mean cloud-base `M_u` is 11.4 kg/(m2 s); of the four tendency caps only `MC_q` binds
  (1.11 % of cells, never > 2x).
- **THE CONVECTION RAINS NOWHERE BECAUSE THE UPDRAFT DETRAINS THE ENVIRONMENT'S CLOUD, AND THE TWO
  KNOBS THAT REPAIR IT ARE BOTH DEFAULT 0** (2026-10-10, README). The condensate recurrence
  subtracts `step*D_u*cloud` — the ENVIRONMENT's cloud, 36.5 g/kg — from an updraft carrying
  7.1 g/kg, at `step*D_u/M_u` = 0.89: 4.6x what it holds, in 64 261 of 64 261 columns, so
  `q_c_u` is zero from the second level up and `g_p = K_p*q_c_u` never forms.
  **`ATM_MC_QC_DETRAIN=1`** (upstream `f90c462`) detrains `q_c_u` instead and rain forms,
  5.2e6 mm/a — **AND MUST NOT BE SET ALONE**: the evaporation demand is then 13 432x the
  generation, applied in full by `rhsForcing()`, and `MC_t` sits on its cap in 51 % of all cells,
  all cooling. **`ATM_MC_EVAP_LIMIT=1`** (upstream `e8d06ce`, the default there) charges
  `e_d + e_p` against the rain that exists: alone a byte-identical null, with `QC_DETRAIN` it
  takes the unmet demand and the `MC_t` truncation to zero. **With both, the downdraft still
  re-evaporates 100 % of the generation and nothing reaches the ground.**
  **`ATM_MC_QVD`** (upstream `6b49292`; the downdraft humidity was built from the cloud-WATER
  array; 1 = environment vapour, 2 = saturated downdraft; default 0) **was the predicted cure and
  is not one**: it cuts the unbacked demand from ~11 000x the generation to 235x (mode 1) or
  467x (mode 2), and `e_d` still takes 100 %. Its one visible effect under the limit is `MC_q`
  off its cap (1.108 % -> 0). **What remains is the evaporation law `e_d` itself** (`a_ev`,
  `b_ev`, charged on a 60 bar layer mass against a `q_sat` of hundreds of g/kg) — a reading, not
  a diagnosis. 20 iterations only; the flip of any of these defaults is the user's call.
- **THE CONVECTIVE RAIN AT THE GROUND DEPENDS ON THE PARITY OF `iter_prec`, AND THE DEFAULT 3 IS THE
  PHASE WHERE ALL OF IT EVAPORATES** (2026-10-10, README). `downdraftEntrainment` forms `e_d` and
  `e_p` from `P_conv` as the PREVIOUS sub-iteration pass left it, so the passes alternate: no
  evaporation, total evaporation, none, total. Measured with `ATM_MC_QC_DETRAIN=1
  ATM_MC_EVAP_LIMIT=1 ATM_MC_ED_AREA=1 ATM_MC_QVD=1`: `iter_prec` 3 and 5 give **0 mm/a** at the
  ground with `e_d` taking 100 %; 4 and 6 give **~3e7 mm/a** with `e_d` taking 3-12 %. Neither is
  the answer (3e7 mm/a is 82 m of rain a day). **Every convective-rain number recorded on
  2026-10-10 is the odd phase, including "the downdraft re-evaporates everything".**
  `ATM_MC_ED_AREA` (upstream `c8f1489`, default 0: `e_d` weighted by the rain-shaft area like
  `e_p`) is ported and cannot show while that is so.
  **`ATM_MC_INLINE` (default 0, not an upstream knob) removes the cycle**: 1 evaluates `e_d` / `e_p`
  at the flux arriving at each level inside `downdraftRecurrence`'s walk, 2 also forms `g_p` from
  the same pass's `q_c_u` inside `updraftRecurrence`. With 2 the result is identical to every
  printed digit for `iter_prec` 3, 4, 5 and 6. **The consistent answer: 2.33e7 mm/a of convective
  rain is formed and none reaches the ground — the downdraft evaporates 91 %, the sub-cloud layer
  9 %, all of it backed.** The probe's "0.3x" is refuted. 2 iterations only.
  **Do not tune a convective constant with `ATM_MC_INLINE` unset, and read no convective-rain
  number from the default `iter_prec` = 3 as a result.**
- **Field names changed 2026-08-20** (ported from ATHAD README item 61). The Brunt-Vaisala
  frequency squared is `brunt_N2` in the code and `BruntVaisala_N2` in the VTK files — it used to
  be `N2`, in the same file as nitrogen's `q_N2`. The eight species are written under bare names
  (`H2O CO2 N2 CH4 NH3 H2 CO SO2`), all mass fractions from `AtmMixture::split()`, so what is
  plotted is what is integrated; the raw transported CO2 array is `CO2_tracer`. **Unverified in
  this fork — no COND run has been made against the renamed writers.**
- **`ATM_MC_GEOPOTENTIAL` PAYS HERE — measured 2026-08-20, and it is the first updraft
  condensation in either fork.** `max c_u` 0 → **3.116e-03 g/kg/s at 24 854 m** (the 7.19e-04 first measured was at the old `M_max` = 3.0, which suppressed the mass flux 4.3× — superseded), with `max s_u`
  moving from the cloud base (661 m) to the top of the column (24 854 m). The reason it works here
  and not in ATHAD is geometry: this fork's convective column is **20 km deep** (base 1786–2200 m,
  LFS 20.5–28.1 km, 47 % of columns convecting) because its triggers have been fractions of surface
  pressure since the fork was cut, while ATHAD's updraft is one grid level deep and has no ascent
  for the geopotential to act on. Nothing integrated moves — OLR, photosphere, Ψ and the water
  column are identical — which is the albedo wall, not a null result about the physics.
  **ON BY DEFAULT here since 2026-08-20** (`ATM_MC_GEOPOTENTIAL=0` restores the old behaviour);
  ATHAD keeps it off, because its updraft is one grid level deep. `M_max` is sized now (100, see
  the table above), so the mass flux is no longer a cap.
- **The 40-iteration OLR is reproducible to only ~0.5 % across a thread-count change in ATHAD**
  (item 61), and this fork shares the convective trigger set that amplifies it. Ψ, albedo and the
  photosphere are unmoved at 1e-6. **State the thread count with any OLR number**, and treat COND's
  own 0.004-0.03 % OLR deltas as inside the envelope until the same null pair is run here.

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


## `ATM_PROGNOSTIC_T` and `ATM_CONV_ADJ_MOIST` (ported from ATHAD 2026-11-26, default 0)

User: "port ATM_PROGNOSTIC_T to ATHAD_PERID and COND", with the moist option. **Default: the temperature is
prescribed.** `ThermoAtm::densities()` rebuilds `t` every iteration as `max(t_skin, adiabat from the surface)`
(moist where saturated), so everything above the point where the adiabat reaches `t_skin` is pinned to
`t_skin` = (F / 2 sigma)^(1/4) = 221 K (`updateSkinTemperature`, `ATM_SKIN_GREY` 1) and has no dynamics.

- `ATM_PROGNOSTIC_T=1` -- in the TIME LOOP `densities(false, true)` keeps the temperature the dynamics,
  radiation and microphysics computed (floored at 180 K, ATHAD's floor) and only re-integrates `p_stat`
  hydrostatically on it. The initial calls (`densities(true)`, the one after `initCO2`) still build the
  prescribed profile, so both modes START from the same state. `ThermoAtm::prognosticT()` reads the knob.
- `atmosphere/ConvectiveAdjustment.h` -- ATHAD's Manabe-Strickler adjustment (whole unstable segments, local
  cp, per-layer critical drop on the stretched grid, dp-weighted enthalpy conserved), copied from ATHAD's
  working tree. It runs ONLY under `ATM_PROGNOSTIC_T=1`, immediately before `densities()`: with a free profile
  it is the only thing holding the lapse rate stable. `ATM_CONV_ADJ_MOIST=1` takes the saturated lapse
  (`SaturationH2O::moistLapse`) where the lower cell is saturated (c >= 0.99 q_sat), the dry one elsewhere.
  Further knobs: `ATM_CONV_ADJ_LAPSE` (scale of the critical lapse, 1), `ATM_CONV_ADJ_PASSES` (sweep cap, 64).

Checks 2026-11-26 (`python/run_progt_check.sh`): knob unset, 2 it, 1 thread = the binary before the port,
every output file byte-identical; `ATM_PROGNOSTIC_T=1 ATM_CONV_ADJ_MOIST=1`, 4 it, 2 threads: exit 0, the
adjustment acts in every column (max 2 sweeps, max dT ~1.3 K, enthalpy drift ~1e-14). 

FIRST LONG RUN 2026-11-26 (run_progt_long, killed at it 50) and two FIXES, both under the knob only (default
byte-identical, `python/run_fix_check.sh` 13/13): the lid ran away cold (204 K at it 50, below the grey skin limit
221 K). A per-stage trace showed why: `apply_radiative_heating()` (radiation_mode 2) changes `t` AFTER
`storeIntermediateData3D()` copied it into `tn`, and the RK step builds the next `t` from `tn`, so every interior
level lost its radiative change every iteration while the lid -- not in the RK loop (i = 1 .. im-2) -- kept it
(-4.7 K per call) and `bcRadius` pinned it to `densities()`' snapshot of itself. (1) Under the knob the heated `t` is
copied into `tn`, so radiation reaches the interior. (a) Under the knob `bcRadius` pins the lid to `t_skin`.
Traced 6 it after the fixes: lid 221.079 K everywhere; the levels below now keep the radiation and COOL ~4 K per
iteration below t_skin (i-1 203 -> 195 K, i-5 208 -> 203 K by it 6), which grey radiative equilibrium does not
allow -- OPEN: suspects the explicit step against the radiative time scale aloft and the 7-40x excess clear-sky
optical depth on the pressure grid (TauLayer vs kappa*q*dp/g*p/p_ref). ATHAD's items 45-46, 65-67 and 79 are
the record of the free top there.


## `ATM_IM` -- the number of radial levels is a build option (2026-11-26), DEFAULT 41 (was 61)

**Default 41 since 2026-11-26** (user: "make im=41 the default for COND and PERID"), as in ATHAD; `make clean && make
cond IM=61` restores the grid of every result before that date. Byte checks (`python/run_imdefault.sh`, 2 it, 1 thread):
the default build = the IM=41 test binary, and IM=61 = the 61-level binary before, 13/13 files each. With the legacy
height grid (no `ATM_GRID_PRESSURE`) im 41 gives ~490 m bottom layers, so `drag_n_layers` 5 reaches ~2.9 km instead
of ~1.8 km; the pressure-grid settings below keep it near 1.8 km. `survival.py` / `psicheck.py` default to `ATM_IM` 41.

User: "is it finally possible to have im=41 in COND and PERID like in ATHAD", then "add the im option". `static const
int im = ATM_IM` in `cAtmosphereModel.h` (default 61, `static_assert` 5..81 = `lib/Array.cpp` MAXI); `make clean &&
make cond IM=41` builds it (objects of different im must not be linked together, hence the clean). Unset = every
output file byte-identical to the build before (`python/run_im41.sh`, 2 it, 1 thread, 13/13).

Use it with the pressure grid and a `ATM_GRID_BETA` that keeps the bottom layer ~300 m, so `drag_n_layers` (a count of
CELLS) keeps its depth: here `ATM_GRID_PRESSURE=1 ATM_GRID_PTOP=1.68e-5 ATM_GRID_BETA=4.09`, `t_skin` 221.05, lid 3.0e-4 bar / 87.86 km
(first five layers ~350-380 m). Test 2026-11-26, 40 it from scratch, im 41 vs the same lid at im 61 (`run_gridp1`):
OLR 159.1 vs 169.3 W/m2 (-6 %), photosphere 66.0 vs 67.4 km (265 vs 256 K), div rms 0.0391 vs 0.0437, RH max 6015 vs 2573 % -- the coarse top (1.09 in ln p) shows here; 22 min at 4 threads. Output and restarts of different im do not mix; `survival.py` / `psicheck.py` read
`ATM_IM`.

---

## Sphere `.vts` and cloudview, written with the panorama (2026-11-27), `ATM_SPHERE` / `ATM_CLOUDVIEW` default 1

Every panorama write (`paraview_panorama_vts_flag`, `panorama_print`) now also writes, in `Paraview_Atm.cpp`:

- `Hadean_Atm_sphere_<n>.vts`: the old `paraview_sphere_vts`, whose call had been commented out, with the giants'
  GIANTS-SPHERE-VTS fix ported: the pole rows j = 0 / jm-1 are left out (extent j = 2..jm-1, 1..179 deg;
  `sin(the) = 0` collapses every longitude and ParaView's Contour stopped with "Cannot compute gradient of grid"),
  and the Cartesian velocity (m/s) is built in locals instead of `aux_u / aux_v / aux_w`, which the model uses as
  work arrays. Fields: Velocity, Temperature (deg C), Bathymetry, u/v/w-component, PressureDynamic, Pressure_bar,
  WaterVapour / CloudWater / CloudIce / CloudGraupel (g/kg), Density (`r_humid`, kg/m3), Precipitation / Rain /
  Snow (mm/d). Display radius `rad.z + 11` as before. ~0.64 GB per write at im 41.
- `Hadean_Atm_cloudview_<n>.vtk` (ported from the giants' `<Planet>_cloudview`): each radial layer i holds the
  column mass ABOVE it, `sum_{i'=i..im-2} q * r_humid * (z(i'+1) - z(i'))` with z = `get_layer_height`, cells
  inside the topography counted 0: `WaterVapour_above_kgm2`, `CloudWater / CloudIce / CloudGraupel /
  Condensate_above_gm2`, and `PhotoRGB` (colour by it, untick "Map Scalars"): a recipe, not radiative transfer --
  a dark sea mixed toward white by the condensate column scaled between its 2nd and 98th percentile of layer 0
  over |lat| <= 80 deg (per file, so two files are not on one scale). Concentric shells, layer 0 outermost,
  r = 1 + `ATM_CLOUDVIEW_SHELL` (im-1-i)/(im-1), default 0.05; pole rows left out. ~0.24 GB per write.

`ATM_SPHERE=0` / `ATM_CLOUDVIEW=0` switch either off. Both only read the model state: byte check 2 it, 1 thread,
repo-default config with every writer at 2, HEAD against the change: 15 / 15 existing output files identical in
ATHAD, ATHAD_COND and ATHAD_PERID, both new files open and contour in VTK.

## `ATM_SFC_DIAG=<n>` -- surface energy balance terms (2026-11-27, from ATHAD_PERID 5a209c0), print-only, default 0 = off

`MultiLayerRadiation.h`. Every n-th call of the radiation, one `[sfc]` line per row at k = 0 for 90 N, 60 N,
30 N and the equator: T_s0 and the T_s the balance sets (the change per call), T_air1 as it entered (the
`densities()` adiabat + the RK step: `ad`) and as the radiative equilibrium makes it (`re`), SW_abs, the
geothermal flux, the down-welling long-wave from both air profiles, sigma T_s0^4, the bulk turbulent term
15 (T_s0 - T_air1) for both, and net = SW + geo + Ldown - sigma T_s0^4 - H for both (`re` = what the balance
uses). Reads only. Byte check 2 it, 1 thread, every writer at 2: 17 / 17 output files identical unset and =1.
First reading on this model (iteration 1, default config): net -19 W/m2 at the pole, +64 / +109 / +125 at
60 N / 30 N / equator, -0.43 .. +2.74 K per call before the omega_rad = 0.05 blend.
On ATHAD_PERID it traced the run_1000 cooling of the ground to the moist block's radial `damp_wiggles(t)` at the
ground row (T0 += (T1 - T0)/4 per moist step); the fix is `ATM_T_DAMP_GROUND` (below, default 1 since `5d11a36`).

## `ATM_T_DAMP_GROUND` -- the vertical Shapiro pass on t leaves the ground row alone (`5d11a36`, 2026-11-28), DEFAULT 1, =0 restores

**This changes the default: every number in this file from before `5d11a36` belongs to the old branch.**
Same change as ATHAD `c27ae39`. `cAtmosphereModel.cpp`, helper `damp_wiggles_t`, at both `damp_wiggles(t)` sites
(set-up and moist block). The 1-2-1 de-checkerboard ran along i from the fluid floor, and with the clamped floor
its update at the ground row is T0 += (T1 - T0)/4: every call pulled the ground a quarter of the way to the first
air layer. The ground row is not an air cell -- it carries `MultiLayerRadiation`'s surface balance and
`densities()` hangs the column and the surface pressure from it. With 1 the horizontal passes are unchanged
(ground row included) and the vertical pass starts one row higher, so level `i_topography+1` is the clamped
floor; same k, j, i order as the single call. `=0` makes the old single call, so it restores the old branch by
construction (no byte check run). Measured on ATHAD_PERID (-1.30 K pole / -1.36 K equator per moist step, the column falling 4.5 K per 100
iterations); **not measured on this model**.

## `ATM_SFC_COUPLE` -- the ground with a heat capacity against the actual level-1 air (`fc5e5ee`, 2026-11-28), default 0 = off

Same change as ATHAD `f9e9e5a`. `MultiLayerRadiation.h` (`sfc_couple_on`) and `apply_radiative_heating` in
`cAtmosphereModel.cpp`. WHY: the ground had no working energy balance in the time loop. The Newton step solves
the balance against the RADIATIVE-EQUILIBRIUM air (which `densities()` then overwrites), has no heat capacity,
and the omega_rad = 0.05 blend keeps 5 % of it per call, so the ground relaxed on an iteration count, not on a
time. With 1, in the time loop's call only (other callers of `run()` keep the Newton step):

    ground   C dT_s/dt = SW + geo + Ldown - sigma T_s^4 - H,   H = c_H (T_s - T_1)
             T_1 = the ACTUAL level-1 temperature as it entered; linearised backward Euler over
             tau = dt L_atm/u_0 * ACCEL (stable for any ACCEL; ACCEL -> infinity is the Newton step
             with the actual air)
    level 1  + H dt L_atm/u_0 / (rho_1 cp_1 dz_1)   (real time, not accelerated)

Both are written AFTER the blend, so the blend no longer touches the ground. `ATM_SFC_HEATCAP` [J/m2/K]
default 2.1e8 (a 50 m water layer), `ATM_SFC_ACCEL` default 1. At ACCEL 1 the ground is practically constant
over any affordable run (one iteration is ~0.03 s on `L_atm` here; +180 W/m2 moves it ~1e-8 K per iteration),
so the ACCEL value is part of the decision. Prints its settings once. Not measured yet; next:
`ATM_SFC_COUPLE=1` with `ATM_SFC_DIAG`.

# Moved from README.md on 2026-11-26

The sections below were the technical part of README.md and were moved here verbatim (user: "CLAUDE.md everywhere"); the README keeps each heading with a one-sentence summary.

## The convection scheme is a moist convective adjustment since 2026-10-10 — the Tiedtke scheme is no longer called

**`ATM_CONV_SCHEME=<adjust|tiedtke>`, default `adjust`, at the user's instruction** ("replace
the Tiedtke scheme with a simple moist convective adjustment"). `ATM_CONV_SCHEME=tiedtke`
restores the mass-flux scheme. **Every number in this file recorded before this section was made
with Tiedtke running**; the sections below it, written the same day, are the measurements that
led here and all of them need `ATM_CONV_SCHEME=tiedtke` to be reproduced.

**What `adjust` is in this tree — and it is less than the name suggests.** A moist convective
adjustment holds a column on the moist adiabat where it is saturated and on the dry one where it
is not, and lets the resulting condensate rain out. Both halves already run here on every moist
step, unconditionally:

- `ThermoAtm::densities()` rebuilds `t` from the surface up with `SaturationH2O::moistLapse`
  (invariant 4) — the saturated adiabat wherever condensation is possible, `g/cp` elsewhere.
  That is the HARD limit of an adjustment: no superadiabatic layer is left for a separate
  routine to find, and one called before `densities()` would be overwritten by it.
- `SaturationAdjustment` condenses what exceeds saturation, with its latent heat, and the ice
  scheme turns the condensate into precipitation.

So **`adjust` adds no routine**. It stops calling `MoistConvection`, and zeroes the four
tendencies that scheme feeds the right-hand sides (`MC_t`, `MC_q`, `MC_v`, `MC_w`), its
precipitation `P_conv` and its updraft / downdraft fields, so that nothing stale is read or
written out. `convectionStep()` in `cAtmosphereModel.cpp` is the one place.

**Why the mass-flux scheme goes.** It is `ATOM_Precipitation`'s Tiedtke scheme with Earth
constants on a 60 bar column, and the sections below measured it the same day: its updraft
condensate is annihilated by detraining the environment's cloud, so it forms no rain; with
that repaired its sub-iteration cycles with period 2 and the answer depends on the parity of
`iter_prec`; made consistent it forms 2.3e7 mm/a of rain and re-evaporates all of it.

**What the replacement changes, measured** — 20 iterations from scratch, 8 threads
(`python/run_ca_tied`, `run_ca_adj`, `run_conv.sh`), last print:

- **120 of 180 printed extrema are identical.** Of the 60 that differ, 37 are the scheme's own
  fields going to zero (`M_u` max 13 210, `MC_q` max 0.02 g/kg/s — its cap —, `MC_t`
  +0.0052 / −0.0044 K/s, `MC_v` ±3.9e-5, `MC_w` −2.5e-4 m/s2, and so on).
- **The other 23 move in the 6th to 8th digit**: max temperature 238.881236 → 238.881216 °C,
  max water vapour 352.600544 → 352.600387 g/kg, `Psi` 2.90858128e13 → 2.90857955e13, max
  radiation 4028.3136 → 4028.3130, precipitation average 199 163.86 → 199 163.62 mm/a; the
  largest relative changes are the evaporation maxima, +0.004 %. (`max presgrad force` reads
  −1e-6 → −8e-6, a field whose printed maximum is at the last digit of the format.)
- `cloud_cap`'s deletion in the last call is −19.156 → −19.169 mm; the saturation budget closes
  in both.
- No NaN, both exit 0; 7 min 21 s against 7 min 42 s wall.

**So at 20 iterations the mass-flux scheme was contributing nothing the large-scale state could
see**, and the replacement is safe at that length. **That is all this measurement says.** 20
iterations are minutes of model time; the run that would show whether the missing convective
transport matters is a long one; it follows.

**The long run, 400 iterations (2026-10-10)** — from scratch, 8 threads, three arms side by
side: `adjust`, `tiedtke`, and `tiedtke` a second time as the noise floor (`python/run_long.sh`,
`run_long_adj`, `run_long_tied`, `run_long_tied2`; 2 h 32 min wall for the three). `cmp_extrema.py` reads the last printed extrema, `cmp_fields.py` the section VTKs of every
checkpoint; the largest difference is taken over the three sections and given as a fraction of
the field's range.

- **The noise floor is zero**: the two Tiedtke runs agree in 180 of 180 extrema and their VTKs
  are byte-identical. Every difference below is the scheme's.
- **104 of 180 extrema identical** between `tiedtke` and `adjust`; of the 76 that differ, 39 are
  the scheme's own fields. The rest stay below 0.03 %: max temperature 238.788879 → 238.788659 °C,
  max water vapour 351.871236 → 351.869503 g/kg, max cloud ice −0.004 %, evaporation maxima
  +0.023 %. (`max presgrad force` −6e-6 → −3.8e-5 is the last digit of its format again.)
- **Fields, at iteration 400**: temperature differs by at most 0.0023 K (rms 0.0007 K), water
  vapour by 0.32 g/kg at 6.7 bar (rms 0.029; 1e-3 of the range), cloud water by 0.037 g/kg,
  `v` by 1.7e-4 and `w` by 1.1e-3 m/s, the dynamic pressure by 3e-4 of its range.
- **The differences GROW LINEARLY and have not levelled off**: water vapour 1.2e-4 of its range
  at iteration 50, 4.8e-4 at 200, 9.6e-4 at 400; cloud water, temperature, `u` and `w` likewise.
  Only `v` and the dynamic pressure flatten.
- The mass-weighted mean water drifts +0.6443 % with `tiedtke` and +0.6480 % with `adjust`.
- One transient: at iteration 350, four cells of the zonal section at 0.015 bar differ by 33 K
  (−52.1 against −18.6 °C); at 400 they agree again to 0.002 K. Not explained.
- All arms exit 0, no NaN.

**So over 400 iterations the convective transport does not matter here** — three to five digits
below the fields themselves. **It is not a statement about a spun-up state**: the difference in
water vapour grows by about 0.08 g/kg per 100 iterations without slowing, so a run a hundred
times longer could reach tens of g/kg if nothing bounds it. That run has not been made.

**What is lost, knowingly:** convective transport of moisture (`MC_q`, which sat on its cap in
1.1 % of cells) and of momentum (`MC_v`, `MC_w`); and `MC_t`, which `densities()` discarded in
any case. A convective adjustment that also MIXES moisture and momentum over the adjusted depth
is the natural next step if a long run asks for it. It is not written.

**Checks.** `ATM_CONV_SCHEME=tiedtke` against the binary from before all of 2026-10-10's changes:
20 of 20 output files byte-identical (2 iterations, 1 thread), and 20 of 20 against the first
run of the day at 20 iterations and 8 threads (`run_dg_20`). Control: `adjust` against
`tiedtke`, 2 of 20.

**Not touched:** ATHAD and ATHAD_PERID still call Tiedtke. ATHAD has a prognostic-temperature
branch (`ATM_PROGNOSTIC_T`) with a DRY Manabe-Strickler `ConvectiveAdjustment.h`; there a moist
adjustment would be a real routine (the same sweep with `moistLapse` where the air is
saturated). This tree has no such branch.

## Three instruments from upstream, and what they read on their first call (2026-10-10)

`ATM_SATADJ_DIAG`, `ATM_MC_DIAG` and `ATM_MC_CAP_DIAG`, all **print-only and default off**, ported
from `ATOM_Precipitation` (`03c8451`, `11ec280`, and `15fd99a` / `e8d06ce` / `f6f5642`), where
between 2026-08-31 and 2026-09-24 they attributed that tree's moist-physics defects one by one.
They are here first because the repairs that followed them upstream are not yet ported, and no
repair should be ported into a tree that cannot measure it.

| knob | where | what it prints |
|---|---|---|
| `ATM_SATADJ_DIAG=1` | `SaturationAdjustment.h` | the routine's water budget per call, clip by clip, in mm of column water: `entry_clip`, `phase_split`, `cold_delete` in `adjustSaturation`; `neg_clip`, `ceiling`, `impossible`, `supersat`, `cap`, `fade` in `clampAndFade`; the total, the unattributed remainder, and the running sum |
| `ATM_MC_DIAG=1` | `MoistConvection.h` | the convective rain budget as an identity, generation − applied `e_d` − applied `e_p` = ground, beside the model's own `P_conv`; the evaporative demand that met no rain; the deep / shallow / midlevel census; `q_c_u`, `M_u`, `sigma_p`, base and LFS heights; sub-cloud RH |
| `ATM_MC_CAP_DIAG=1` | `MoistConvection.h` | how many cells each of the four tendency caps truncates and by how much (> 2x, > 10x, > 100x), `MC_t` split into its mass-flux and latent halves, and where the over-cap cells sit |

**What had to change against upstream.** The shallow / midlevel tests are fractions of the
surface pressure (`f_stat_*`); saturation is `SaturationH2O`'s with the exact mass fraction, and
a cell with `q_sat >= 1` has no relative humidity and is left out of the RH means; the rain walk
ends at level 1 because this tree sets `P_conv[0] = P_conv[1]`; the saturation budget has three
buckets upstream does not need (`ceiling`, `impossible`, and graupel in the total, since
`evaporateWhereImpossible` empties it); the cap census is a SECOND pass after `rhsForcing()`
instead of counters inside it, so that `rhsForcing()` is untouched but for reading its three
caps from one definition (`McCaps`). Upstream's later probes under `ATM_MC_DIAG` split by land
and ocean and by Earth regions and are not carried.

**Print-only, verified.** 4 iterations from scratch, 1 thread, `config_port_new.xml` with `nm` = 4
and `checkpoint` = 2 (`python/run_dg_pre`, `run_dg_new`, `run_dg_on`): the binary before the
change against the new one with the knobs unset, **20 of 20 output files byte-identical**; the
new one unset against all three knobs set, **20 of 20**. No control arm was run (nothing shows
that this comparison would have caught a difference), and the printed budget figures are OpenMP
reductions of doubles, so their last digits may move with the thread count.

**First readings** — 20 iterations from scratch, 8 threads, all three on (`python/run_dg_20`,
eleven moist calls, exit 0). They are readings, not diagnoses; none has been followed up.

1. **`cloud_cap` deletes about 19-22 mm of column water in every call but the first**, −22.23
   falling to −19.16 mm, **−202 mm over eleven calls** against a precipitable water of
   160 626 mm (0.13 %). The cap (`cloud_cap` = 0.05 kg/kg) is 99 % of everything this routine
   removes. Its comment says it "SETS the condensate rather than bounding a pathology"; the
   budget puts a number on what that costs.
2. **The cold fade removes 0.29 mm per call**, steady from the second call — the sigmoid upstream
   measured as "a global exponential decay on cloud water" (`dfcfd7c`), present here unchanged.
3. **The phase split creates water**, +0.0012 rising to +0.0056 mm per call (+0.054 mm in eleven
   calls): the defect `ATM_SAT_PHASE_AVAIL` repairs, default off, now with a size. It is 3 500x
   smaller than the cap's deletion.
4. `entry_clip`, `cold_delete`, `neg_clip` and `ceiling` are exactly zero; `impossible` and
   `supersat` are zero to 1e-17. Both budgets close: unattributed ≤ 8e-15 mm.
5. **Convective rain generation is exactly zero** — `g_p` 0.0000 mm/a in all eleven calls, and
   the model's own `P_conv` at the surface 0.0000 — **although every column is active and
   classified deep** (100 %, none shallow or midlevel), every column reaches a level with
   `K_p > 0`, the updraft carries condensate (mean `q_c_u` 0.29 g/kg over 26 cells per column, max
   7.8 g/kg) and base and LFS are never on adjacent levels (mean base 1.8 km, LFS 24.7 km).
   So the scheme convects everywhere and rains nowhere. **The cause was found the same day —
   next subsection.** (The two candidates first written here, the `t_u >= t_0` gate and the
   `is_land || t_u <= t_00` reset, are both EXCLUDED: every cell is above `t_0`, 388-500 K.)
6. **Mean cloud-base mass flux is 11.4 kg/(m2 s)** (12.9 at the first call), where the
   inherited comments call ~0.3 healthy on Earth; `mc_M_max` is 100 here.
7. **Of the four caps only `MC_q` binds: 43 439 cells, 1.11 % of the fluid, and never beyond
   2x** (59 235 at the first call). `MC_t`, `MC_v`, `MC_w` are truncated in 0 cells, in both
   halves of `MC_t`. So ATHAD's item 53 concern — `MCt_max` as an Earth constant bounding this
   column — does not show at 20 iterations with `g_p` at zero; it has to be read again once
   the convection rains.

### Why the convection rains nowhere: the updraft detrains the ENVIRONMENT's cloud (2026-10-10)

`[MC-GP]`, a census added to `ATM_MC_DIAG` for this question, counts the three things
`g_p = K_p * q_c_u` needs in one cell — `t >= t_0`, a level more than `delta_i_c` above cloud
base, and updraft condensate — and their intersections, stage by stage through the convective
passes (the stage-by-stage form was a temporary probe; one line at the end of the call is kept).

- **Before the first `updraftRecurrence` the condensate is there**: `q_c_u > 0` in 1 735 765 of
  1 800 026 cells (the seed `initUpdraft` writes), and all three conditions hold in 1 671 504.
- **After it, `q_c_u > 0` in 64 261 cells — exactly one per active column — and all three
  conditions hold in 0.** The surviving cell is the first level above cloud base, the one level
  where `K_p` = 0 (within `delta_i_c` = 500 m of the base). Nothing changes in passes 2 and 3.
- **The line that does it** is the condensate recurrence,
  `M_u*q_c_u + step*(-D_u*cloud - rho*g_p)`: it detrains `cloud`, **the environment's cloud
  water**, out of the updraft's condensate, where a detraining updraft loses its own `q_c_u`.
  Measured at the step base+1 → base+2: mean `step*D_u/M_u` = **0.89** (one level detrains 89 %
  of the mass flux — `eps_u` = 0.2/`R_cloud` = 2e-3 /m on ~450 m layers), mean environment cloud
  **36.5 g/kg** against **7.1 g/kg** carried, so the detrained term is **4.6x** what the updraft
  holds, **in 64 261 of 64 261 columns**. `q_c_u` goes negative and is clipped to zero from
  there up.
- **The parcel does condense** — `c_u > 0` in every cell above the first — and that condensate
  is among what is wiped.

**Upstream has the same line and repaired it as `ATM_MC_QC_DETRAIN`** (`f90c462`, 2026-10-01,
default off there). On Earth the environment's cloud is ~0.009 g/kg, so the same error
UNDER-detrains and the condensate grows; here the environment holds five times what the updraft
does and the error annihilates it. One defect, opposite sign of consequence — the knob is
ported under its upstream name, **default 0**.

**THE A/B, AND WHY THE KNOB MUST NOT BE SET ALONE.** `ATM_MC_QC_DETRAIN=1`, 20 iterations from
scratch, 8 threads (`python/run_el_q`): `q_c_u > 0` in 1 554 466 cells, `g_p > 0` in 445 529,
**generation 5.2 million mm/a** (7.7 million at the first call). And with rain forming, the
downdraft and sub-cloud evaporation `e_d + e_p` — which `rhsForcing()` applies in full to `MC_t`
and `MC_q` whether or not the rain exists — **demand 13 432x the generation** (UNMET
7.0e10 mm/a): **`MC_t` truncated in 51.1 % of all cells, its latent half > 100x the cap in
32.4 %, 100 % of it cooling** (largest term `e_d` in 80.8 %, `e_p` in 19.2 %), and `MC_q`
truncated in 52.2 %, > 100x in 9.8 %. Half the atmosphere on a cap.

### `ATM_MC_EVAP_LIMIT` (2026-10-10, default 0)

Ported from `ATOM_Precipitation` `e8d06ce`, where it has been the default since 2026-09-22 and
the knob is retired. In `downdraftRecurrence`, before `P_conv` is formed: the level's `e_d + e_p`
is charged against the rain arriving plus the rain made there and both are scaled together; at
`i = 0`, where the flux passes through unchanged, both are zeroed. Same site and same
arithmetic as upstream.

Four arms, 20 iterations from scratch, 8 threads, all three instruments on, final call:

| | default (`run_dg_20`) | `EVAP_LIMIT` alone (`run_el_e`) | `QC_DETRAIN` alone (`run_el_q`) | both (`run_el_qe`) |
|---|---|---|---|---|
| generation `g_p`, mm/a | 0 | 0 | 5 216 971 | 5 203 328 |
| UNMET evaporative demand, mm/a | 0 | 0 | 7.0e10 (13 432x) | **0** |
| `MC_t` truncated, % of cells | 0 | 0 | 51.12 (> 100x: 32.33) | **0** |
| `MC_q` truncated, % of cells | 1.108 | 1.108 | 52.23 (> 100x: 9.83) | **1.108** |
| share of generation taken by `e_d` | – | – | 100 % | 100 % |
| convective rain at the ground, mm/a | 0 | 0 | 0 | **0** |
| mean `q_c_u`, g/kg | 0.289 | 0.289 | 1.355 | 1.355 |
| precipitation average (ice scheme), mm/a | 199 163.86 | 199 163.86 | 198 629.24 | 199 163.91 |
| `cloud_cap` deletion, final call, mm | −19.156 | −19.156 | −19.319 | −19.158 |
| max T, °C | 238.881236 | 238.881236 | 238.883411 | 238.881248 |

- **Alone it is a null, as it must be**: with `g_p` = 0 there is no demand to limit. Every printed
  figure equals the default's, and at 2 iterations, 1 thread, the output is byte-identical
  (20 of 20 files).
- **With `ATM_MC_QC_DETRAIN` it removes the whole eruption**: unmet demand to zero, `MC_t` off its
  cap everywhere, `MC_q` back to the default's 1.108 %. At the first call 39 849 cells (1.02 %)
  still exceed the `MC_t` cap, all cooling, all from `e_l` (the detrained condensate) and none
  beyond 2x; by the last call there are none.
- **AND THE CONVECTION STILL DOES NOT RAIN AT THE GROUND.** The downdraft evaporation takes
  **100 % of the generation**, now within what exists. So the two knobs together turn "no rain
  is formed" into "all rain formed is re-evaporated by the downdraft", and the large-scale state
  barely notices (precipitation average 199 163.86 → 199 163.91 mm/a, max T in the 8th digit).
  **(Tested the same day, next subsection: it is NOT enough.)**
  Upstream met the same wall — "`e_d` takes 100 % of it" — and its repair is `ATM_MC_QVD`
  (`6b49292`): the downdraft humidity is built as `0.5*(cloud + 0.98*q_sat)`, **half of the
  cloud-WATER array**, where Tiedtke's mix is environment VAPOUR. That line is here unchanged
  (`MoistConvection.h`, twice in `downdraftEntrainment`) and with 36 g/kg of cloud it is not a
  small error. NOT PORTED, not measured here: it is the next candidate, not a diagnosis.

**Byte checks** (2 iterations from scratch, 1 thread, `python/run_vqcd.sh`, `run_el.sh`):
the binary before all of this against the new one with every knob unset, **20 of 20 files**;
the three instruments on, 20 of 20; `ATM_MC_EVAP_LIMIT=1` alone, 20 of 20; controls
`ATM_MC_QC_DETRAIN=1` against unset 7 of 20, and both knobs against `QC_DETRAIN` alone 2 of 20.

**Not measured:** anything longer than 20 iterations; either knob in ATHAD or ATHAD_PERID;
what happens to the water the parcel condenses (`c_u` debits the environment's vapour through
`MC_q`) while its condensate is wiped — on the default branch that is a candidate water sink,
and the column budget instrument upstream has for it (`ColumnWaterBudget.h`) is not ported.

### `ATM_MC_QVD` (2026-10-10, default 0): it cuts the demand 24-47x and the downdraft still takes everything

Ported from `ATOM_Precipitation` `6b49292` in upstream's form (`qvDown()`, two call sites in
`downdraftEntrainment`): the downdraft humidity `e_d` and `e_p` evaluate against is
`0.5*(cloud + 0.98*q_sat)` — half the cloud-WATER array — by default; `=1` uses the
environment's vapour, `0.5*(c + 0.98*q_sat)`; `=2` a saturated downdraft, `0.98*q_sat`.

**The prediction written above this section was that this is "the next candidate" for why no
convective rain reaches the ground. It is not enough, in either mode.** 20 iterations from
scratch, 8 threads, `ATM_MC_QC_DETRAIN=1 ATM_MC_EVAP_LIMIT=1`, final call
(`python/run_el_qe`, `run_qv1`, `run_qv2`):

| | `QVD=0` | `QVD=1` | `QVD=2` |
|---|---|---|---|
| generation `g_p`, mm/a | 5 203 328 | 5 218 599 | 5 218 565 |
| share taken by `e_d` | 100 % | **100 %** | **100 %** |
| convective rain at the ground, mm/a | 0 | **0** | **0** |
| `MC_t` truncated, cells | 0 | 0 | 0 |
| `MC_q` truncated, % of cells | 1.108 | **0** | **0** |
| precipitation average (ice scheme), mm/a | 199 163.91 | 199 163.81 | 199 163.81 |
| max T, °C | 238.881248 | 238.881255 | 238.881255 |

**How far off it is, measured with the limit OFF** so that the raw demand shows (2 iterations,
8 threads, `ATM_MC_QC_DETRAIN=1`; `run_el_q`, `run_qvraw1`, `run_qvraw2`; first and second
call): the evaporative demand that no rain backs is

    QVD=0   10 996x and 13 005x the generation     MC_t truncated in 54.8 % of cells
    QVD=1      235x and    272x                                      34.7 %
    QVD=2      467x and    553x                                      50.7 %

So the cloud-water term was worth a factor of 24-47 in demand — it IS a large error here, as
expected with 36 g/kg of cloud — **and a demand 235x the supply remains with it gone.** Mode 1
demands half of what mode 2 does. The one thing it changes in the limited run is `MC_q`: the
43 439 cells (1.108 %) that sat on that cap in every run today are off it in both modes.

**What is left is the evaporation law itself, and that is a reading, not a diagnosis.**
`e_d = a_ev*(1 + b_ev*P^(1/6))*(q_sat - q_v_d)*P^(4/9)` with `a_ev` = 1e-3 and `b_ev` = 5.9, and
the demand it is charged as is `step*rho*e_d`. Two things in that product are not Earth's: the
layer mass `step*rho` (60 bar, ~450 m layers) and `q_sat`, which is hundreds of g/kg in a
388-500 K column, so even the 2 % deficit of a saturated downdraft is several g/kg. Whether
the coefficients, the deficit or the downdraft's area weighting (upstream's `ATM_MC_ED_AREA`,
`ATM_MC_ED_ABOVE_BASE`, both unported) is the one to move has NOT been looked at.

**Byte check** (`python/run_qvd.sh`): the binary before any of today's changes against the new
one with every knob unset, 20 of 20 files (2 iterations, 1 thread); controls `QVD=1` against
`QVD=0` 2 of 20, `QVD=2` against `QVD=1` 2 of 20 (20 iterations, 8 threads).

### The downdraft evaporation law, factor by factor: the law is not the outlier, the layer is (2026-10-10)

A probe build (`[MC-ED]`, not kept in the source; logs `python/run_ed0`, `run_ed1`, `run_ed2`)
re-evaluated the law on every level that carries convective rain, with
`ATM_MC_QC_DETRAIN=1 ATM_MC_EVAP_LIMIT=1`, 2 iterations from scratch, 8 threads, second call.
No code changed; this is a reading of what is there.

**The law** is `e_d = a_ev*(1 + b_ev*P^(1/6))*(q_sat - q_v_d)*P^(4/9)`, `a_ev` = 1e-3, `b_ev` = 5.9:
the Kessler rain-evaporation rate written on the precipitation flux, as in COSMO. It is a rate
per kg of AIR, and `downdraftRecurrence` charges it as `step*rho*e_d` — to the whole grid layer.

| mean over the ~590 000 levels with rain | `QVD=0` | `QVD=1` | `QVD=2` |
|---|---|---|---|
| layer mass `step*rho`, kg/m2 | 16 062 | 16 114 | 16 066 |
| temperature / `q_sat` | 443.7 K / 235 g/kg | 442.8 K / 233 g/kg | 443.0 K / 234 g/kg |
| deficit `q_sat - q_v_d`, g/kg | **103.5** | **2.20** | **4.68** |
| rain available, kg/(m2 s) | 0.0161 | 0.0163 | 0.0163 |
| law rate per unit deficit, 1/s | 5.30e-4 | 5.34e-4 | 5.36e-4 |
| levels whose demand ≥ ALL the rain arriving | 100 % | 79.9 % | 100 % |
| law demand at the rain actually there / generation | **53.2x** | **1.1x** | **2.4x** |
| the same law AREA-WEIGHTED, `sigma*E(P/sigma)` / generation | 15.9x | **0.3x** | **0.7x** |
| Tiedtke's bound, `sum |M_d|*0.98*dq_sat` / generation | 2.95x | 2.95x | 2.95x |

(At the first call the `QVD=1` figures are 1.2x, 0.4x and 2.49x.)

Three things, in the order they matter:

1. **The deficit under the default is 44 % of `q_sat`** — 103.5 g/kg — because half of the
   downdraft "humidity" is the cloud-water array. `ATM_MC_QVD` takes it to 2.2 or 4.7 g/kg, a
   factor 22-47, which is the factor measured in the unlimited runs above. That part is done.
2. **The layer is 16 000 kg/m2.** The law's rate is ordinary — 5.3e-4 /s per unit deficit at
   58 mm/h of rain — but it is applied to every kilogram of a layer that holds 20-50x what an
   Earth layer does (18.6 kg/m3 over ~1070 m; the ~450 m quoted earlier in this section is the
   spacing near cloud base only), with the deficit held fixed while it evaporates. So even at
   a 2 % deficit one level asks for all the rain arriving in 80 % of the levels. **The rain falls
   in a shaft, not in the whole cell**: `e_p` already carries the precipitating area fraction
   `sigma_p` (0.08 here), `e_d` does not. Weighted the same way — which is exactly upstream's
   `ATM_MC_ED_AREA` (`c8f1489`), unported — the demand is **0.3x the generation with `QVD=1` and
   0.7x with `QVD=2`**. That is the first combination today under which convective rain would
   reach the ground; it is a number computed from the fields, NOT a run. **(Run the same day, next
   subsection: the model does not do what this computation assumes.)**
3. **Independently of this law, a downdraft this strong would take the rain anyway.** Tiedtke's
   own requirement — evaporate what keeps the descending air saturated — is `|M_d|*dq_sat` per
   level: `|M_d|` = 3.4 kg/(m2 s) (`gam_d` = 0.3 of a cloud-base `M_u` of 11) against a `q_sat` that
   rises 4.8 g/kg per level on the way down, **2.95x the generation**. So the scheme's rain
   production is small against the water its own downdraft can absorb, whatever law is used.
   That points at the mass-flux closure (`c_mb`, `gam_d`, both Earth / COSMO constants on a
   60 bar column) and at how little of the updraft's water becomes rain (`g_p = K_p*q_c_u` with
   `q_c_u` capped at 10 g/kg), and none of that has been examined.

**So: the coefficients `a_ev`, `b_ev` are not the thing to move.** The cheap next step is
`ATM_MC_ED_AREA` with `ATM_MC_QVD=1`; the larger question behind it is point 3.

### `ATM_MC_ED_AREA` (2026-10-10, default 0), and what it uncovered: THE ANSWER DEPENDS ON THE PARITY OF `iter_prec`

Ported from `ATOM_Precipitation` `c8f1489` in upstream's form: `e_d` is evaluated in the rain
shaft and weighted by its area, `sigma_p*E(P/sigma_p)`, as `e_p` already is.

**The subsection above computed that with this and `ATM_MC_QVD=1` the demand would be 0.3x the
generation and rain would reach the ground. Run, it does not** — 20 iterations from scratch,
8 threads, `ATM_MC_QC_DETRAIN=1 ATM_MC_EVAP_LIMIT=1 ATM_MC_ED_AREA=1` (`python/run_ea0`,
`run_ea1`, `run_ea2`):

| | `QVD=0` | `QVD=1` | `QVD=2` |
|---|---|---|---|
| generation `g_p`, mm/a | 5 203 328 | 5 218 620 | 5 218 569 |
| share taken by `e_d` | 100 % | 100 % | 100 % |
| convective rain at the ground, mm/a | 0 | 0 | 0 |

With `QVD=0` the output is **byte-identical** to the run without the knob (20 of 20 files against
`run_el_qe`, at 8 threads): where the demand exceeds the supply before and after, the limiter
returns the same numbers and the knob cannot show.

**Why the computed 0.3x was not what the model does.** The probe evaluated the law at the rain
ACTUALLY ARRIVING at each level, inside one downward walk. The model does not:
`downdraftEntrainment` forms `e_d` and `e_p` from `P_conv` **as the previous `iter_prec` pass left
it**, and `downdraftRecurrence` then walks the column with that fixed demand. So pass 2 sees
pass 1's flux (zero: `g_p` is forced to 0 on the first pass) and evaporates nothing; pass 3 sees
pass 2's UN-EVAPORATED flux, demands accordingly, and evaporates everything; a fourth pass would
see pass 3's zero flux and evaporate nothing again.

**Tested, and it is exactly that** — same knobs with `QVD=1`, 2 iterations from scratch, 8 threads,
only `<iter_prec>` changed (`python/run_ip4`, `run_ip5`, `run_ip6`; first and second call):

| `iter_prec` | generation, mm/a | share taken by `e_d` | convective rain at the ground, mm/a |
|---|---|---|---|
| 3 (the default) | 7.68e6, 5.45e6 | 100 %, 100 % | 0, 0 |
| 4 | 3.17e7, 3.07e7 | 5.2 %, 2.9 % | **3.00e7, 2.98e7** |
| 5 | 9.92e6, 7.30e6 | 100 %, 100 % | 0, 0 |
| 6 | 3.01e7, 2.93e7 | 11.6 %, 6.5 % | **2.66e7, 2.74e7** |

**The convective rain at the ground is 0 for an odd pass count and ~3e7 mm/a for an even one.**
The sub-iteration is a period-2 oscillation, not a convergence, and the shipped `iter_prec` = 3
samples the phase in which everything evaporates. The generation swings with it (5-10e6 against
3e7), because the updraft's condensate is debited by the previous pass's `g_p` in the same lagged
way. `param.py`'s note on `iter_prec` — "min 3 for evaporation to act on non-zero `P_conv`; check
convergence at 4-5" — describes the lag and did not see the cycle.

**Neither phase is the answer.** 3e7 mm/a is 82 metres of rain a day, 150x the ice scheme's whole
199 164 mm/a; 0 is the other extreme. Every convective-rain figure in this section, and the
"e_d takes 100 %" wall of the three subsections above, is the ODD phase of this cycle.
The self-consistent value needs `e_d` and `e_p` evaluated at the flux arriving at each level,
inside `downdraftRecurrence`'s own downward walk — a restructuring, not a port (upstream has the
same lag, and there one level does not usually ask for all the rain, so it shows less;
not checked in that tree). NOT DONE. **(Done the same day as `ATM_MC_INLINE`, next subsection.)**

**Byte check** (`python/run_eda.sh`): the binary before all of today's changes against the new
one with every knob unset, 20 of 20 files (2 iterations, 1 thread); control `ED_AREA=1` with
`QVD=1` against `QVD=1` alone, 2 of 20 (20 iterations, 8 threads).

### `ATM_MC_INLINE` (2026-10-10, default 0): the sub-iteration made pass-count independent

Not an upstream knob. `=1` evaluates `e_d` and `e_p` inside `downdraftRecurrence`'s downward
walk, at the flux arriving at the level plus the rain made there, with the law, the deficit
(`qvDown`) and the area weighting exactly as `downdraftEntrainment` writes them. `=2` also forms
`g_p` inside `updraftRecurrence`'s upward walk, from the `q_c_u` the SAME pass has just computed
on the level below, before it is debited from the level above. The gates are unchanged
(nothing on the first pass, nothing where `is_land || t <= t_00`).

`ATM_MC_QC_DETRAIN=1 ATM_MC_EVAP_LIMIT=1 ATM_MC_ED_AREA=1 ATM_MC_QVD=1`, 2 iterations from
scratch, 3 threads, only `<iter_prec>` changed (`python/run_in<mode>_<iter_prec>`); generation
in mm/a at the first and second call, and where the rain goes at the second:

| `iter_prec` | lagged (previous subsection) | `ATM_MC_INLINE=1` | `ATM_MC_INLINE=2` |
|---|---|---|---|
| 3 | 7.68e6, 5.45e6 → ground 0 | 7.68e6, 5.45e6 → ground 0 | **2.470e7, 2.329e7 → ground 0** |
| 4 | 3.17e7, 3.07e7 → ground 2.98e7 | 3.17e7, 3.07e7 → ground 0 | **2.470e7, 2.329e7 → ground 0** |
| 5 | 9.92e6, 7.30e6 → ground 0 | 9.92e6, 7.30e6 → ground 0 | **2.470e7, 2.329e7 → ground 0** |
| 6 | 3.01e7, 2.93e7 → ground 2.74e7 | 3.01e7, 2.93e7 → ground 0 | **2.470e7, 2.329e7 → ground 0** |

- **Mode 2 is identical for 3, 4, 5 and 6 passes to every printed digit** — generation
  24 702 532.5162 and 23 289 005.7756 mm/a, `e_d` 23 126 676.4427 and 21 243 803.7100, `e_p`
  1 575 856.0735 and 2 045 202.0655, the cap census cell for cell. The cycle is gone.
- **Mode 1 removes the evaporation's flip and leaves the generation's**: the ground is 0 at every
  pass count, but the generation still alternates 5-10e6 / 3e7, because `g_p` is still formed
  from the previous pass's condensate. That is why mode 2 exists.
- **THE SELF-CONSISTENT ANSWER IS THAT NO CONVECTIVE RAIN REACHES THE GROUND**: of 2.33e7 mm/a
  generated, the downdraft evaporates **91.2 %** and the sub-cloud layer **8.8 %**, all of it
  backed by rain that exists (unmet demand 0). So the lagged scheme's "0" was the right number
  for the wrong reason, and its "3e7" was not a number at all.
- **The 0.3x computed two subsections above is refuted a second time, now by the consistent
  scheme.** It was evaluated on the odd phase's state, where the generation is 4-5x smaller than
  the consistent one, and it did not include `e_p`. What does hold from that subsection is its
  point 3: a downdraft carrying 0.3 of this updraft's mass flux can absorb about three times
  the rain formed.
- `MC_t` is truncated in 78 627 cells (2.01 %), none beyond 2x; `MC_q` in none.

**Byte check** (`python/run_inl.sh`): the binary before all of today's changes against this
build with every knob unset, 20 of 20 files (2 iterations, 1 thread).

**Not done:** nothing longer than 2 iterations was run with this knob; the
large-scale effect is not measured; `e_l`, the detrainment `D_u` at the LFS and
`CloudBaseUpdate` still depend on the pass number by construction (first pass against later
ones), which the identity across 3-6 passes shows does not matter from the third pass on.

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

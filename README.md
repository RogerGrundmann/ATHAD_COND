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

> **Status: under construction.** The sections below are inherited from ATHAD and are being
> replaced phase by phase. Anything not yet marked for ATHAD_COND still describes the
> 250 bar / 1500 K parent and should not be quoted for this model.

## Repository layout

```
atmosphere/   C++ source — the atmosphere model
cli/          Command-line interface source and config file
lib/          Shared utilities (arrays, geometry, …)
python/       Python bindings (Cython)
test/         Self-tests
tinyxml2/     Embedded TinyXML-2 XML parser
param.py      Parameter definitions (generates .inc, .pxi and .xml at compile time)
Makefile      Top-level build file
CLAUDE.md     Project constants, invariants and open risks
```

## Prerequisites

| Dependency | Notes |
|---|---|
| g++ ≥ 7 or clang++ ≥ 6 | C++17 required |
| OpenMP | Typically ships with the compiler |
| Cython, Python ≥ 3.6 | Only for the Python interface |

## Build

```bash
make had          # the CLI executable, cli/had
make all          # CLI + Python extension
```

## Run

```bash
cd python && OMP_NUM_THREADS=8 ../cli/had config_athad.xml
```

Output lands in `python/output_Hadean/` as ParaView `.vtk` slices and `.vts`
panoramas — the ATOM line's convention (`output_<name>/`, relative to the run
directory). Run from `cli/` instead and it writes `cli/output_Hadean/`; the
`<output_path>` entry in the config controls it either way.

## Configuration

`cli/config_athad.xml` is generated from `param.py` and carries every parameter with
its default and a description. You need only include the entries you want to change;
everything else falls back to the compiled default.

## Hadean conditions

Surface pressure 250 bar, surface temperature ≈ 1500 K, mean molar mass 21.43 g/mol,
R_mix = 387.9 J/(kg·K), surface density 43 kg/m³. Composition by mole fraction:
80 % H₂O, 10 % CO₂, 3 % N₂, and 1.4 % each of CH₄, NH₃, H₂, CO and SO₂. H₂O and CO₂
are prognostic; the rest are a fixed well-mixed background.

Full constants, the three model invariants, and the assumptions that are still open
are in [CLAUDE.md](CLAUDE.md).

## Status and limitations

Under construction. This section records what has actually been measured — including
the measurements that did not work out — rather than what is intended.

1. **Bootstrap (done).** The tree builds and links as `libathad.a` + `cli/had`.
   Hydrosphere, paleogeography and the NASA/Scotese data pipeline are removed.

2. **Flat surface (done).** `init_topography()` prescribes `h ≡ 0`,
   `i_topography ≡ 0`; the run reports 65341 surface points, all water. The surface
   temperature is prescribed from `t_surf_equator` / `t_surf_pole` (1500 K / 1450 K,
   area mean 1483 K) instead of read from a paleo curve. A 2-iteration run completes
   and writes 20 ParaView files, **bit-identically at 1, 4 and 8 OpenMP threads**.

   Three defects had to be fixed to get there, all inherited and all latent on Earth:

   - `MoistConvection::findCloudBaseLFS` indexed `t.x[-1][j][k]`. Its scan for
     `p_stat <= 1000 hPa` leaves the index at its `-1` sentinel when no level
     qualifies; the two other consumers of that index guard the sentinel, this one
     did not. Latent on Earth (p_surf ≈ 1013 hPa always fires by level 1), live on a
     column that sits entirely above 1000 hPa.
   - `AtomUtils::GetMean_2D/3D` built the global `m_node_weights` lazily with
     `clear()` + `push_back()`. `printDataAtm()` calls them from ~9 concurrent OpenMP
     sections, so the first call raced several threads through the same reallocation
     and corrupted the heap. It was masked in the parent only because an earlier
     single-threaded `GetMean_2D(temperature_NASA)` happened to run first.
   - `get_temperatures_from_curve()` dereferenced `begin()` and decremented `end()`
     *before* its own `size() < 2` guard — undefined behaviour on an empty map, which
     is what every call became once the Scotese curves were gone.

   The Makefile also gained `-MMD -MP` header dependencies; without them, edits to
   the headers where nearly all the physics lives did not trigger a rebuild.

   **What is not yet right.** The physics is still Earth physics on a 16 km shell.
   The column mean settles near 332 K and the prescribed 1500 K survives only at the
   domain top, because radiation still relaxes toward an Earth-like target, the
   saturation formula is still Magnus (capped ~101 °C), and the shell is ~20× too
   shallow for a 250 bar atmosphere. Deep convection is inactive: its trigger
   thresholds (1000/970/900/800 hPa) are absolute Earth surface pressures and never
   fire here. Phases 2–5 address these.

3. **Composition and mixture thermodynamics (done).** `MixtureAtm.h` derives the mass
   fractions, mean molar mass and gas constants from the configured mole fractions and
   supplies the *local* mixture properties `R_of(c, co2)`, `cp_of(c, co2, T)` (Shomate
   fits, 298–2000 K) and `M_of(c, co2)`. The constant `R_Air` / `cp_l` were replaced at
   the sites that set the pressure, density and Poisson exponent. CO₂ became a well-mixed
   mass fraction (ppm is meaningless at 20 % by mass) and water is initialised well-mixed
   rather than as a fraction of a saturation value that does not exist.

   Measured at the equator, initial state — every target hit exactly:

   | quantity | target | measured |
   |---|---|---|
   | surface pressure | 250 bar | 249.947 bar |
   | surface temperature | 1500 K | 1499.5 K |
   | surface density | 42.97 kg/m³ | 42.9695 kg/m³ |
   | mixture gas constant | 387.9 J/(kg·K) | 387.915 |
   | mean molar mass | 21.434 g/mol | 21.4337 |
   | H₂O mass fraction | 0.6724 | 0.6724 |

4. **Domain and non-dimensionalisation (done).** The shell went from 16 km to 300 km
   (`L_atm` is the stretch *amplitude*, so the shell is `(exp(zeta)−1)·L_atm`), `zeta`
   3.715 → 3.0 and `im` 41 → 61. The top cell now spans 1.65 local scale heights instead
   of 2.92. The column reaches **0.0237 bar** at the equator and 0.0051 bar at the pole,
   both below the 0.1 bar radiating level. `dt_visc` was scaled by the `dr²` CFL ratio,
   `abl_height` by the scale-height ratio, and the COSMO lapse parameter `beta` is now
   *derived* (`cosmo_lapse_fraction · R_mix · T_surf / cp`) rather than Earth's 42 K —
   which would have given a 0.71 K/km, essentially isothermal, 300 km column.

   Diagnostics added: a per-column profile and a per-level global summary
   (`ThermoAtm::printColumnProfile` / `printLevelSummary`), printed at init and every
   checkpoint. The level summary is what located two of the defects below.

   Five more inherited defects surfaced, all latent on Earth:

   - **`dr` was hard-coded to 0.025**, silently tied to `im = 41` (0.025 × 40 = 1). With
     61 levels the radial span became 1.5, one rad.z unit was read as 932 km instead of
     300, and the domain top landed at 1399 km. Now derived as `1/(im−1)`.
   - **`SaturationAdjustment` capped the prognostic temperature at 333.15 K and wrote the
     cap back** — "well above any physical surface temperature" is an Earth statement, and
     it collapsed 250 bar to 30 bar in one iteration. Replaced by a physical bound plus
     the correct statement: above 647.096 K water is supercritical and there is nothing to
     condense, so the adjustment is a genuine no-op there.
   - **The COSMO profile drives T → 0 at finite height** (305 km at the 1500 K equator,
     285 km at the 1450 K pole, since the coefficient goes as 1/T²). Temperature was
     floored at `t_00` but pressure kept falling at a rate its own temperature no longer
     justified. Now continues isothermally above the floor.
   - **`initCloudIce` manufactured cloud from a negative `q_sat`.** Magnus `E_sat` at
     1500 K (~1.2 × 10⁷ hPa) exceeds the 250 bar column, so `ep·E_sat/(p−E_sat)` goes
     negative and `c − H_crit·q_sat` *adds* water. It produced a condensate mass fraction
     of 0.47, inflating density by 1.87×. Same root cause drove the surface evaporation
     scheme to write a water mass fraction of 21.
   - **`bcRadius` extrapolated `p_stat`, `r_humid` and `r_dry` cubically at the lid.** The
     file already documents this stencil overshooting through zero for velocity and
     amplifying concavity for turbulence — but the hydrostatic quantities were left on it.
     On a 300 km shell it drove `p_stat` to −36 hPa across ~11 500 top-level cells. They
     now use log-linear extrapolation, which is exact for an isothermal layer and positive
     by construction.

   Also fixed: `co2Atmosphere()` ran *after* `densities()`, so the density was built with
   `R_of(c, 0) = 414.2` instead of 387.9 — a 7 % error through the whole column. Harmless
   on Earth, where CO₂ was ppm and never entered the density.

   **Stability.** A 20-iteration run at 8 threads completes cleanly: no NaN, no
   non-positive pressure at any of the 5 checkpoints, 34 ParaView files written. The
   equatorial surface settles rather than drifting away —

   | iter | T [K] | p [bar] | ρ [kg/m³] |
   |---|---|---|---|
   | 0  | 1499.5 | 249.947 | 42.9695 |
   | 5  | 1497.1 | 249.615 | 42.9222 |
   | 10 | 1496.0 | 249.360 | 42.8784 |
   | 15 | 1495.6 | 249.307 | 42.8692 |
   | 20 | 1495.5 | 249.280 | 42.8646 |

   with the per-checkpoint change decaying 2.4 → 1.1 → 0.4 → 0.1 K. The domain top rises
   from 0.024 to 0.049 bar over the run and stays below the radiating level.

   **What is not yet right.** Radiation is still the inherited Earth scheme (Bignami /
   Atwater–Ball emissivity, `radiation_mode = 5` relaxing toward a target that no longer
   exists), so the small residual drift above has no physical meaning — it is an Earth
   parameterisation being evaluated far outside its calibration, not a Hadean climate
   settling. Deep convection remains inactive: its trigger thresholds are absolute Earth
   surface pressures. Between 373 K and 647 K the saturation curve is still Magnus, far
   outside its validity. These are Phases 4 and 5.

5. **Saturation physics (done).** `SaturationH2O.h` replaces the Magnus formula
   everywhere — all ~25 call sites across 8 files — with the IAPWS correlations:

   - Wagner & Pruss (IAPWS-95) saturation line, 273.16 K to the critical point;
   - IAPWS (2011) ice-Ih sublimation curve below the triple point;
   - Watson latent heat, correctly vanishing at the critical point;
   - the **exact** mass-fraction conversion `q = x·M_H₂O / (x·M_H₂O + (1−x)·M_other)`,
     replacing the dilute `ep·E/(p−E)` that assumed water is a trace;
   - dew point by bisection (IAPWS has no closed-form inverse) and dq_sat/dT from
     Clausius–Clapeyron rather than the Magnus fit's own derivative.

   `make test` runs `test/saturation_selftest.cpp` against published reference points.
   All pass: the liquid curve to ~1e-5 relative (6.11657 hPa at 273.16 K, 1013.25 hPa at
   373.124 K, 220640 hPa at 647.095 K), sublimation to ~6e-3, plus monotonicity,
   [0,1]-boundedness, and agreement with the dilute form in the dilute limit.

   Writing the test first paid for itself: it caught a sign error in the sublimation
   formula (the summand is `a_i·θ^b_i`, not `a_i·(1−θ^b_i)` — the wrong form gives
   49 hPa instead of 0.76 hPa at 250 K) and two wrong expectations of my own.

   **Result: there is no condensation anywhere in the model.** The equatorial column runs
   1499.5 K at the surface to 676.7 K at 243 km and never crosses 647.096 K, so every
   level is supercritical — no cloud, no rain, no latent heat. The column printout now
   reports the saturation state per level and says so explicitly.

   That is physically coherent for a runaway greenhouse, but it is *not yet a result*: it
   follows from the prescribed lapse rate, and `cosmo_lapse_fraction = 0.5108` is still
   Earth calibration. At the dry adiabat (4.81 K/km) the column would cross 647 K near
   177 km and condensation would begin there. Which is right is for the radiation to
   decide, not for a prescribed lapse — see Phase 5.

6. **Radiation (done).** The Bignami/Atwater–Ball emissivity is replaced by a grey
   optical depth built from the absorber mass each layer actually holds, with pressure
   broadening:

   ```
   tau_i = SUM_s kappa_s * q_s * (dp_i/g) * (p_i/p_ref)
   ```

   The pressure-broadening factor is what no Earth-calibrated emissivity fit contains and
   what matters most here — at 250 bar it is a factor of 250 over the 1 bar reference.
   The old fit could not work in principle: `eps = 0.684 + 0.0056*e_surf` saturates to
   0.999 once `e_surf` exceeds ~56 hPa, and here it is ~2×10⁵ hPa, so every layer became
   a blackbody and the scheme carried no information about composition at all.

   Also changed: `radiation_mode` 5 → **2** (direct σT⁴ heating; modes 0/1/3/4/5 all lean
   on the Scotese snapshot or the 280 ppm CO₂ reference, neither of which exists at
   4.4 Ga); insolation rescaled to the faint young Sun (0.71 S₀); a **geothermal flux**
   added to the surface energy balance, since a molten surface supplies heat from below;
   albedo made a single cloud-deck value, Earth's ice/ocean contrast having no subject
   here; and the hard-coded `287.0` J/(kg·K) in the cloud-density calculation replaced by
   the local mixture gas constant.

   **The headline measurement.** Outgoing longwave flux at the top of the atmosphere:

   | | value |
   |---|---|
   | OLR | **324.5 W/m²** |
   | σ·T_surf⁴ | 286 629 W/m² |
   | suppression | **×883** |
   | transmitted fraction | 0.000 |

   and it stays at 324.5 W/m² at *every* checkpoint while σ·T_surf⁴ drifts from 286 629 to
   286 462. That near-total independence from surface temperature is the defining
   signature of a runaway greenhouse — the emission level sits inside the optically thick
   water column, not at the ground. The value is ~5 % above the canonical Nakajima /
   Komabayashi–Ingersoll limit of 280–310 W/m², which is close for a grey scheme whose
   three κ coefficients were chosen a priori rather than fitted.

   **What this does and does not establish.** It establishes that the opacity is the right
   order of magnitude and that the scheme reproduces the correct qualitative behaviour. It
   does not establish that 324.5 W/m² is the right number: `kappa_H2O`, `kappa_CO2` and
   `kappa_bg` carry roughly a factor-of-two uncertainty, and they are the biggest single
   lever on the answer. A grey scheme also cannot represent the window regions that set
   the real limit. Treat the value as a consistency check that passed, not as a
   prediction.

   Still open: the column remains supercritical throughout (Phase 4), so there is no
   cloud deck to justify the 0.4 albedo that is being assumed — the two are inconsistent,
   and resolving it needs the upper atmosphere to cool enough to condense.

7. **Lapse rate and albedo made self-consistent (done).** This closed the inconsistency
   Phase 5 left standing — a 0.4 cloud albedo asserted over a column that condensed
   nothing.

   The root cause was `cosmo_lapse_fraction = 0.5108`, carried over from Earth. That
   number is not a tuning constant: Earth's ~5 K/km sits half way to its 9.8 K/km dry
   adiabat *because of latent heat release*. ATHAD's deep column condenses nothing, so
   there is nothing to flatten the lapse, and the consistent value is the dry adiabat.
   Carrying Earth's value across was self-fulfilling — the flattened lapse was justified
   by condensation that the flattening then prevented.

   Setting the fraction to 1.0 alone broke the model (surface NaN, 30 bar instead of 250),
   which exposed a deeper problem: the inherited COSMO profile is `T = T₀√(1−coeff·h)`, a
   **sqrt** in height. Matching its near-surface slope to the adiabat does not make it an
   adiabat — it plunges to zero at 156 km, inside the domain. So the fitted profile was
   replaced by the physics it was standing in for, integrated layer by layer:

   ```
   dry adiabat     dT/dz = -g/cp          (cp local: follows composition and T)
   hydrostatic     dp/dz = -p*g/(R*T)     (R local, on the layer-mean T)
   isothermal top  T = t_skin             where the adiabat falls below it
   ```

   No tuned constant, exact for a constant-cp adiabat, and it cannot produce the zero
   temperature the sqrt form did.

   Consequences measured:

   - The atmosphere is far more compressed on the true adiabat, so the domain came back
     from 300 km to **230 km**; 300 km put the top 90 km into near-vacuum (6e-8 bar).
   - **A cloud deck now forms**, from ~207 km (0.09 bar) upward, where p_H₂O finally
     exceeds p_sat. The top two levels report `CONDENSING`. The albedo is no longer
     asserted: `MultiLayerRadiation` uses the dark molten-surface value (0.08) and the
     existing condensate-driven cloud bump raises it where the model actually makes cloud.
     (The `albedo_pole`/`albedo_equator` parameters turned out to be **inert** — the
     radiation module builds its own albedo and never read them, so the earlier 0.4 was
     never in effect.)
   - `t_skin` is now derived from the **energy budget**, σT⁴ = (1−α)·SW + geothermal =
     236.0 W/m² → 254.0 K, replacing a value taken from a previously measured OLR, which
     was circular.

   **Energy balance closes:** OLR = 236.0 W/m², exactly the absorbed SW plus geothermal,
   with the surface suppressed by ×1214.

   **And that is the interesting result.** 236 W/m² is *below* the 280–310 W/m² runaway
   limit. At 0.71 S₀ with the assumed 150 W/m² geothermal flux, the planet does not absorb
   enough to sustain a runaway greenhouse — so the 1500 K surface is being held by fiat
   (it is prescribed), not by the budget. For a 1500 K surface to be self-consistently in
   runaway, the geothermal flux would have to be **≥ ~195 W/m²** (absorbed ≥ 280). That is
   not implausible for a genuine magma ocean, but it is a prediction the model now makes
   rather than an assumption it was given, and `geothermal_flux` should be revisited
   against magma-ocean cooling estimates.

   Still not closed: `t_skin` is a one-shot estimate using the clear-sky albedo, while the
   cloud deck raises it (albedo 0.4 would give 245.5 K). Making it a true fixed point means
   iterating `t_skin` against the model's own albedo. `initComposition()` prints both values
   and warns when they diverge.

8. **Rotation, paleo-time removal, Python bindings (done).**

   **Rotation.** ω = 3.17e-4 rad/s — a 5.5 h Hadean day, 4.35× modern, reported at
   startup. Checked against the family's four Coriolis/centrifugal sign fixes rather than
   re-derived. Two findings:

   - The **dynamical** Coriolis signs in `RHS_Atm_Turb.cpp` already agree with ATURAN
     `8b284cb` / ATNEPT `024c37f` once the opposite storage convention is accounted for
     (they store −a because their RHS subtracts; this file stores +a because its RHS adds).
     ATOM_Precipitation had fixed this independently. No change needed.
   - The **diagnostic** in `ThermoAtm::forces()` did not match. It carried a
     +2Ω·sinθ·u term the momentum equations drop under the traditional approximation, so
     the ParaView "Coriolis force" field showed a force that was never applied. It now
     follows the same `ATOM_CORIOLIS_NONTRAD` switch as the dynamics.
   - The centrifugal diagnostic used `(1 + |sinθ|)`, which is **maximal at the pole** where
     the true value is zero. Replaced by the correct decomposition about the rotation axis,
     `a_r = ω²r sin²θ`, `a_θ = ω²r sinθcosθ` (ATURAN `4201957` / ATJUP `8649675`).
     Diagnostic only — ATHAD's RHS carries no centrifugal term, the force being curl-free
     and absorbed by the pressure projection, which is what ATURAN found.

   **Paleo-time removal.** Nineteen parameters deleted (topography grids, NASA surface
   fields, Scotese curves, the pygplates reconstruction script, `Ma_switch`, `Ma_max`,
   `t_paleo_max`, `co2_paleo`, …), plus `read_Hydrosphere_SST()` and its two parameters —
   the atmosphere half of an atm↔ocean Picard loop with no ocean to couple to.
   `time_start/end/step` remain: the time-slice loop is still structural, and a single
   slice is all that runs.

   **Python bindings.** `pyatom` → `pyathad`, hydrosphere halves stripped. `model.py` was
   driving *both* spheres and — at the bottom — actually running the ocean rather than the
   atmosphere; it now drives the one model that exists.

   **Diagnostic cadence** (separate from `checkpoint`, which writes ParaView files):
   every **10** iterations for a short run (`nm ≤ 100`), every **100** for a longer one,
   overridable with `diagnostic_stride`. The chosen value is reported at startup.
   Verified: nm=20 prints at iters 10 and 20; nm=400 selects the 100 stride.

   A 20-iteration run completes clean — no NaN, no non-positive pressure, OLR steady at
   236.0 W/m².

9. **The saturation conversion finished, and 60 km of manufactured cloud removed (done).**

   Phase 5 replaced Magnus with IAPWS everywhere, which fixed the saturation *pressure*.
   It did not finish the *conversion* from that pressure to a mass fraction. About twenty
   sites still carried the dilute form `q_sat = ep·E/(p − E)` that invariant 2 in
   CLAUDE.md forbids, and the equatorial column printout was showing the consequence
   without anyone reading it that way: water vapour oscillating between 0 and 0.88 above
   140 km, and a density that *rose* with height at 142 km.

   Three distinct defects, all latent on Earth:

   - **`SaturationAdjustment::clampAndFade` had the superheated branch inverted.** It read
     `q_sat = (p > E_sat) ? ep·E_sat/(p − E_sat) : ep·1e-5`. When p_sat(T) exceeds the
     local pressure the vapour is superheated and **cannot condense at all**, so the
     saturation limit is 1 — all of the water stays vapour. The fallback asserted the
     opposite, a limit of ~7e-6, and the block below it therefore dumped the **entire
     0.67 water mass fraction into cloud** and released its latent heat, in exactly the
     layers where nothing can condense. On this column that is every level between ~373 K
     and the critical point: a **60 km slab of manufactured cloud from ~140 to ~200 km**,
     which set the planetary albedo and, smeared downward by `damp_wiggles`, inflated the
     density in the supercritical layers beneath it through the `(1 − cloud − ice)`
     loading term. Latent on Earth because no terrestrial cell is ever above 373 K, so the
     branch never ran — the same shape of defect as the 333.15 K temperature cap.
   - **The fix from Phase 5 was applied to the entry `q_sat` but not to the copy inside
     the same function's Newton loop**, which went on pulling q_v toward `ep·1e-5`.
   - **`AtmMixture::M_nonwater` could not do what its name and comment claim.** It took
     `(co2, M_background)` and set `q_b = 1 − q_c`, so its renormalisation "to exclude
     H₂O" was identically 1 and the water's share of the mass was silently handed to the
     background gas. It returned 28.58 g/mol where the reference composition gives 35.11 —
     the carrier 19 % too light, hence a q_sat some 23 % too large in every exact
     conversion in the model. Invisible on Earth, where water is ~1 % of the mass. It now
     takes the local water mass fraction as well.

   All remaining dilute sites went with them: `initCloudIce` (guarded above the critical
   point but not in the subcritical superheated band, where its q_sat goes negative and
   the `c − H_crit·q_sat` subtraction manufactures cloud), all four ice schemes,
   `RHS_Atm_Turb`'s `q_Rain`/`q_Ice` thresholds — where a negative threshold switches the
   latent-heat term on unconditionally — and the dead `init_vapour_cloud`, which is one
   uncomment away from being live.

   Diagnostics added, because the defect was visible for weeks in output nobody could
   read: the column profile now prints `q_cld` and `q_ice` and flags a level `CLOUD?!`
   when it carries condensate that cannot exist there, and a new
   `ThermoAtm::printPlanetaryBalance` prints the cos(latitude)-weighted mean albedo,
   absorbed shortwave, OLR and the implied skin temperature every diagnostic step.

   **Measured, 20 iterations, 8 threads, clean (no NaN, pressure positive everywhere):**

   | | before | after |
   |---|---|---|
   | condensate at 142–200 km | up to 0.42 | ≤ 0.05, and falling |
   | q_H₂O at 186 km | 0.0000 | 0.9177 |
   | T at 218 km | 254.0 K | 427.3 K |
   | domain top | 0.028 bar | 0.295 bar |
   | mean planetary albedo | — (not measured) | 0.4999 |

   The upper column is **much warmer**, and that is the point: with the bogus condensation
   gone it stays steam, keeps steam's high cp, and cools along a shallower adiabat. The
   234 km shell no longer reaches the radiating level.

   **Two things this breaks, both now open:**

   - **The shell is too shallow again.** The top is at 0.295 bar, above the 0.1 bar
     radiating level. The 300 → 230 km reduction in Phase 7 was made against the profile
     the manufactured latent heat produced, so it has to be re-derived.
   - **The planetary albedo is pinned at 0.4999**, which is the hard-coded thick-cloud
     asymptote `alpha_cloud = 0.50` in `MultiLayerRadiation`. Once any condensate is
     present the reflectivity `tau/(tau+2)` saturates, so the model is not *computing* an
     albedo, it is *reporting a constant* — and that constant is a bare literal inside a
     physics kernel, not a parameter. It has to be exposed and justified before the OLR or
     the runaway claim that depends on it means anything.

   With that albedo the budget no longer closes: absorbed SW + geothermal = 203.8 W/m²
   against an OLR of 236.0, an imbalance of **−32.3 W/m²**, because `t_skin` is still the
   configured 254 K while the energy balance now implies 244.8 K. That is the fixed point
   below, and it is no longer optional.

10. **The insolation, the albedo parameters, the t_skin fixed point — and what they
    exposed: the OLR is not an output (done, and it is bad news).**

    Three fixes, and then the thing they uncovered.

    **The insolation was Earth's surface flux used as top-of-atmosphere insolation.**
    `rad_equator_short` / `rad_pole_short` were 116 / 71 W/m², which is
    ATOM_Precipitation's 163.3 / 100.0 scaled by 0.71 — and those are Earth's
    absorbed-at-the-*surface* shortwave, already reduced by Earth's albedo and its
    atmosphere's absorption. `MultiLayerRadiation` uses them as the flux incident at the
    top and then applies ATHAD's own albedo, counting Earth's albedo twice. The
    cos(latitude)-weighted mean was **107.5 W/m² where 0.71 S₀ delivers 241.6**: the planet
    was lit at 45 % of its own insolation. Replaced by a fit of the model's parabola to the
    annual-mean insolation of a zero-obliquity planet, `Q(φ) = S/π·cos φ`, constrained so
    the weighted mean is exactly S/4 — 298.0 at the equator, 0.0 at the pole, RMS 6 W/m².
    Measured mean after the change: **241.56 W/m²**. (Hadean obliquity is unknown. Zero is
    the assumption and also the only shape the parabola can follow; at Earth's 23.44° the
    true curve flattens toward the pole in a way it cannot.)

    **The albedo constants became parameters.** `alpha_cloud = 0.50` and the molten-surface
    0.08 were bare literals in `MultiLayerRadiation` while the config carried an
    `albedo_pole`/`albedo_equator` pair that nothing read. They are now `albedo_cloud` and
    `albedo_surface`, and the inert pair is gone. This matters more than tidiness:
    the cloud deck's condensate path is ~10⁵ g/m² against a `cwp_tau` of 100, so
    `refl = τ/(τ+2)` saturates and **every cloudy column returns `albedo_cloud` to four
    decimals**. The measured mean planetary albedo is 0.4999. The model is not computing an
    albedo; it is reporting that constant.

    **`t_skin` is now an iterated fixed point** against the model's own albedo,
    `σ·t_skin⁴ = (1−ᾱ)·S̄ + F_geo`, relaxed by `t_skin_relax` (0.25) after each radiation
    call. It converges from the configured 254 K to **262.85 K** against a target of
    262.88 K within 20 iterations. Setting `t_skin_relax = 0` restores the old fixed value.

    **What all of that exposed.** The model's outgoing longwave flux is not a computed
    quantity. It is σT⁴ of the topmost layer, which is prescribed.

    | run | lid T | lid ε | reported OLR | σ·T_lid⁴ |
    |---|---|---|---|---|
    | `t_skin` = 254 K | 254.0 K | 1.0000 | 236.01 W/m² | 236.01 |
    | `t_skin` = 240 K | 240.0 K | 1.0000 | 188.13 W/m² | 188.13 |
    | after these fixes | 343.98 K | 1.0000 | 802.32 W/m² | 793.84 |

    The first two rows are a deliberate experiment: changing one configured number moved
    the "measured" OLR to σ times that number to five significant figures. The lid is
    optically thick (τ ≈ 6 at 0.29 bar), so the emission level sits on the domain boundary
    and the model reports the boundary condition back.

    **The 236.0 W/m² of Phase 7, and the "energy balance closes" that went with it, was
    this.** It was worse than a coincidence: `bcRadius` pins the lid temperature to a
    snapshot `t_top_init` taken once from the initial condition, and that snapshot came
    from `initTemperatureData`, which builds its adiabat *before* the water and CO₂ fields
    exist and therefore with a background-only cp and a much steeper lapse. `densities()`
    rebuilt the whole column on the correct mixture cp every iteration — and could never
    move the lid, because bcRadius pinned it straight back to the stale value. The lid
    stayed at 254.0 K, the OLR stayed at σ·254⁴, and the agreement with the energy budget
    was the budget being read back out of the number it had been used to set. `densities()`
    now refreshes the snapshot, which is what makes the 343.98 K above visible.

    **So the honest current state is a −531 W/m² imbalance and no meaningful OLR**, and
    that is a better description of the model than the closed budget it replaces.

    **The shell cannot currently be deepened to fix it.** 230 km tops out at 0.29 bar,
    three times above the 0.1 bar radiating level. Going deeper — 260 km (top 0.017 bar)
    and 300 km (top 3.2e-4 bar) — produces a finite, sane initial state and then **NaN
    across the entire field in the first `MultiLayerRadiation` call**: the field is finite
    in the diagnostic immediately before it and NaN in the one immediately after, in both
    runs. The scheme's tridiagonal (Thomas) assembly is built from products of the layer
    emissivity and its rows degenerate as ε → 0 — which is exactly the condition at the top
    of any domain that reaches the radiating level. The radiation has to be reformulated
    before the shell can grow, so 230 km stays.

    Also fixed in passing: `python/PythonStream.cpp` still included `pyatom.h`, so the
    default `make` target had not built the Python bindings since the fork. It does now.

11. **The radiation scheme rewritten as two-stream flux sweeps; the shell deepened to
    300 km; the OLR is a computed quantity for the first time (done).**

    The inherited scheme assembled a tridiagonal system whose every entry was a product of
    the layer emissivity — sub-diagonal `ε_{i-1}σT_{i-1}⁴`, diagonal `−2ε_iσT_i⁴`,
    super-diagonal `ε_{i+1}σT_{i+1}⁴` — with a right-hand side built from differences of
    cumulative transmitted sums, solved by Thomas for a correction to `σT⁴`. Optical depth
    goes as p² through the pressure broadening, so near the top adjacent rows differ by
    orders of magnitude while the right-hand side is a difference of two nearly equal
    sums. It did not survive a domain that reaches the radiating level (item 10).

    Replaced by the plane-parallel grey transfer it was standing in for, integrated
    directly:

    ```
    up[i] = up[i-1]·(1 − ε_i) + ε_i·σT_i⁴        surface → top,  up[0] = σT_surf⁴
    dn[i] = dn[i+1]·(1 − ε_i) + ε_i·σT_i⁴        top → surface,  nothing enters from space
    ```

    Radiative equilibrium of a layer is then the statement that it absorbs what it emits,
    `ε_i(up[i-1] + dn[i+1]) = 2ε_iσT_i⁴`, **and the emissivity cancels**:

    ```
    σT_i⁴ = (up[i-1] + dn[i+1]) / 2
    ```

    Nothing divides by ε. The ε → 0 limit is exactly right rather than merely survivable —
    a transparent layer passes both streams and takes the mean of what goes by — and at
    the top, where `dn → 0`, it reduces to `σT⁴ = up/2`, the classical skin temperature
    that this model has until now been *prescribing* as `t_skin`. The temperature and the
    fluxes are iterated against each other (Lambda iteration, 4 passes); convergence is
    slow in optically thick layers, the known weakness of the method, but there the update
    is nearly a no-op anyway, and the fluxes are exact for the current profile at every
    pass.

    **Measured, 20 iterations, 8 threads, no NaN, pressure positive everywhere:**

    | shell | top p | lid ε | OLR | σ·T_lid⁴ | |
    |---|---|---|---|---|---|
    | 230 km | 0.29 bar | 1.0000 | 795 W/m² | 787 W/m² | OLR *is* the lid — an input |
    | 260 km | 0.017 bar | 0.0610 | 519 W/m² | 271 W/m² | decoupled |
    | 300 km | 3.8e-4 bar | 0.0000 | 581 W/m² | 271 W/m² | **a real column integral** |

    The 230 km row is the regression check: where the old scheme worked, the new one
    reproduces it (795 against the old 802 W/m²). The other two rows are runs that
    previously NaN'd in the first radiation call.

    **`L_atm` is now 300 km.** The column tops out at 3.8e-4 bar with a transparent lid,
    the isothermal skin is resolved from 256 km up, condensation begins at 242.8 km, and
    the outgoing flux is no longer the boundary temperature read back out.

    **And the first thing the model says with it is that its own opacity is too low.**
    OLR = 581 W/m² against 271 W/m² absorbed plus geothermal: the atmosphere radiates away
    more than twice what it takes in, so it cannot hold the prescribed 1500 K surface. That
    is now a statement about `kappa_H2O` = 0.01 m²/kg rather than an artefact of the
    boundary. Note also that the OLR is **not yet grid-converged** — 519 W/m² at 260 km
    against 581 at 300 km, because `im` is fixed at 61 and a deeper shell is a coarser one.

    Two consequential changes came with it:

    - **`radiation.x` is now the upward long-wave flux at the top of each layer**, not
      `σT⁴` of that layer. It is a diagnostic field (nothing feeds it back into the
      dynamics), it is continuous across the surface by construction — so the 1-2-1
      smoothing pass that used to hide the surface kink is gone — and `radiation.x[im-1]`
      is the OLR, which is what the mode-5 cloud diagnostic already assumed it was.
    - **`bcRadius` no longer pins the radiation lid** to `σ·(t·t_0)⁴`. That pin was
      justified by "radiation.x = σ·(t·t_0)⁴", which has stopped being true; keeping it
      would have thrown away the one number the column integration exists to produce.

    **A separate defect, found by asking whether the vertical stretch reached everything:**
    `init_tropopause_layers` converted a height to a level index with `round(h / L_atm)`.
    That is a level index only for uniformly spaced layers, and the grid is exponentially
    stretched — `height(i) = (exp(zeta·i/(im−1)) − 1)·L_atm`, so `L_atm` is the *amplitude*
    of the stretch, not a spacing. The exact inverse is
    `i = (im−1)·ln(1 + h/L_atm)/zeta`. At 300 km the pole's 195 km convective top sits at
    level **52**; the old formula returned **12**, which is 13 km. `VelocityInitializer`
    builds its entire jet profile between the surface and that level and applies a linear
    taper to zero from it up to the domain top, so the initial wind structure was
    compressed into the bottom 4 % of the atmosphere and the remaining 96 % got the taper.
    Wrong on Earth too (28 of 41 levels = 4.9 km, not the intended 11 km), but wrong there
    in a way that still landed inside the troposphere, so it never showed.

12. **CO₂ made genuinely prognostic, and a 100-iteration run (done).**

    The CO₂ in the ParaView output was *exactly* constant — `min co2 = max co2 = 0.205300`
    in every cell at every checkpoint. `ThermoAtm::co2Atmosphere()` fills the whole field
    with the uniform `co2_0·co2_scale`, and it was being called **inside the time loop**,
    immediately before `densities()`. Meanwhile the model does carry a full CO₂ transport
    equation — `RHS_Atm_Turb` builds `rhs_co2`, `RungeKutta_Atm_Turb` carries it through
    all four stages — and every iteration the result was discarded. CLAUDE.md's "H₂O and
    CO₂ are prognostic" was false for CO₂.

    `co2Atmosphere()` is now the initial condition only. In its place in the loop,
    `ThermoAtm::co2Column()`:

    - fills `co2_total`, which was declared ("areas of higher co2 concentration") and
      **never written**, so its min/max and the `co2_average` derived from it both read
      0.000. It now carries the column CO₂ mass path in kg/m², the CO₂ analogue of
      precipitable water: 519 108 kg/m² global mean, 505 205 (pole) to 522 343 (equator);
    - reports the drift of the global mass-weighted mean q_CO₂. There is no CO₂ source or
      sink anywhere in the model, so that mean is conserved and any drift is transport
      error. Measured over 20 iterations: **−0.0000 %**.

    Unit labels corrected with it: the 3-D field is kg/kg not ppm, the column is kg/m², and
    the startup banner's `co2_0=0.205 ppm` is kg/kg.

    **The field still comes out uniform — and that is now the answer rather than the
    assumption.** A passive tracer with no sources and no gradients has `∇q = 0`, so
    advection and diffusion both vanish and uniform is the exact solution. The difference
    is that the model now computes it, monitors it, and will transport any structure that
    does appear instead of erasing it every step.

    **100-iteration run**, 300 km shell, 24 threads, no NaN, ParaView every 10 iterations.
    The model relaxes monotonically toward radiative balance:

    | iter | 10 | 20 | 40 | 60 | 80 | 100 |
    |---|---|---|---|---|---|---|
    | OLR [W/m²] | 617.8 | 581.0 | 522.0 | 477.0 | 441.7 | 413.4 |
    | imbalance [W/m²] | −347.0 | −310.2 | −251.2 | −206.2 | −170.9 | −142.6 |

    The imbalance decays by about 9 % per 10 iterations and the rate is slowing, so
    reaching balance (OLR → 271 W/m²) needs of order 400 iterations. Note the OLR passes
    *through* the 280–310 W/m² Nakajima / Komabayashi–Ingersoll band on the way down —
    which is the first time this model has approached that limit from a computed flux
    rather than a prescribed one. `t_skin` has converged to 262.88 K and the lid emissivity
    is 0.0000, so the outgoing flux is a column integral throughout.

13. **The Hadley cells were not symmetric, and nothing here can make them so (done).**

    Spotted in the iteration-0 ParaView output. `VelocityInitializer::compute()` set the
    meridional wind at 15°N to a surface coefficient of **4.0** and at 15°S to **3.0**:

    ```cpp
    init_v_or_w(m.v,  75, -3.0,  4.0);   // lat:  15N
    init_v_or_w(m.v, 105, -3.0,  3.0);   // lat:  15S
    ```

    Every other mirror pair in the whole u/v/w initialisation is identical — 0/180, 15/165,
    30/150, 45/135, 60/120 — and the `form_diagonals` spans are mirrored too, so this was
    the only asymmetry in the velocity initial condition. The commented-out lines that sat
    directly above each showed the values had once been the other way round (3.0 north,
    4.0 south), so it had never been symmetric; someone had swapped which hemisphere won.

    Measured in `meridional_streamfunction_10.csv`, at the level of maximum |Ψ|:

    | lat | Ψ(+lat) | Ψ(−lat) | sum (0 if symmetric) |
    |---|---|---|---|
    | 30° | 74 187 | −72 049 | 2 138 |
    | 15° | −278 138 | 366 213 | **88 075** |
    | 10° | −198 131 | 260 698 | 62 567 |

    (10⁹ kg/s.) The southern cell is stronger than the northern by 366/278 = **1.32**
    against the coefficients' 4.0/3.0 = 1.33, and the global antisymmetry error is 11.4 %.
    At 30°, away from the injected asymmetry, the mismatch is 2.9 %. It does not wash out:
    still 8.3 % at iteration 100.

    On Earth this asymmetry is physical — the ITCZ sits north of the equator because of the
    land–sea distribution. **ATHAD cannot have it.** There is no land, no topography, the
    prescribed surface temperature is a symmetric parabola, the insolation profile is
    explicitly mirrored (`short_wave_radiation[j] = short_wave_radiation[j_max-j]`), and
    there is no obliquity and no seasonal cycle. Nothing in the model can sustain a
    hemispheric asymmetry, so all of it was inherited from these two numbers.

    Both are now **3.5**, their mean, which removes the asymmetry and leaves the total
    initial Hadley mass flux unchanged. Measured after the fix, same diagnostic:

    | lat | Ψ(+lat) | Ψ(−lat) | sum |
    |---|---|---|---|
    | 30° | 73 094 | −73 138 | −44 |
    | 15° | −322 169 | 322 181 | 12 |
    | 10° | −229 399 | 229 418 | 19 |

    **Global antisymmetry error 11.38 % → 0.0152 %**, the residue being floating-point and
    OpenMP reduction noise rather than structure.

14. **No condensed phase where none can exist (done) — and what enforcing it exposed.**

    The cloud deck in the ParaView output sat at 166–230 km, and every level of it was
    already flagged `CLOUD?!` by the model's own column diagnostic: condensate in cells
    whose `q_sat` is 1.0000, which is what "no saturation limit" means. The genuine
    condensation level is 242.8 km.

    Cause: **none of the four ice schemes, nor MoistConvection, contained a single
    critical-point check** — `grep T_CRIT` across all five files returned nothing.
    `TwoCatIceScheme` ran its full warm and cold microphysics from the ground up, and the
    precipitation maxima landed at i=50 (175.8 km) and i=54 (218.2 km), inside the band.
    Condensate formed legitimately at 243 km, sedimented into air at 350–700 K where it
    should flash to vapour instantly, and the scheme kept it and shuffled it between rain,
    snow, cloud and ice. Latent on Earth, where no cell is near 647 K or above its own
    saturation pressure, so the question never arises.

    `IceSchemeCommon::canCondense()` now states the two conditions — supercritical
    (T ≥ 647.096 K, where liquid and vapour are one phase) and superheated (p_sat(T) > p,
    where the vapour cannot reach saturation whatever its abundance) — and
    `evaporateWhereImpossible()` enforces them: condensate goes back to the vapour with its
    latent heat (no latent heat above the critical point, where there was never a separate
    phase), and the sources and precipitation fluxes are cleared. It is called from all
    four ice schemes and from `SaturationAdjustment::clampAndFade`, whose "supercritical
    cells carry no condensate and need no fade" `continue` was exactly the wrong response —
    they carry none only if something removes it, and nothing did.

    **Result: zero `CLOUD?!` levels anywhere in the run**, from every level between 166 and
    230 km before.

    **What it exposed.** With the condensate returned to the vapour, `q_H2O` in the
    sub-cloud band went to **1.27** — a mass fraction above one. It turned out the field had
    been unphysical all along, just less visibly: before the fix it sat at 0.9971 against
    `co2` = 0.2053, a composition summing to 1.20. Nothing bounded it. The three mass
    fractions must sum to one and the background is carried as the remainder `1 − c − co2`,
    so `c > 1 − co2` is a **negative background mass** — and `AtmMixture::split()`
    renormalises defensively, so the only symptom was a gas constant pinned at 415.1 instead
    of moving with the composition. `UtilsAtm::valueLimitationAtm` carries such a bound at
    the wrong ceiling of 1, and is commented out at both of its call sites.

    The ceiling `c ≤ 1 − co2` is now enforced every iteration, outside the moist-physics
    block (which runs only on moist iterations, while the RK4 transport of `c` runs on all
    of them). The column is physical again: `c` = 0.7947 through the band, R = 405.6.

    **But the ceiling bites in 282 302 cells, and that is an open problem, not a fix.**
    Water is pumped downward out of the single condensing level by sedimentation and
    evaporates into the superheated band with no return path, so `c` there grows until the
    clamp stops it — which deletes water. The count is printed at every diagnostic step.
    A run in which it stays large is not to be trusted, and this one does.

15. **The water is not lost, it is created — and the moist physics was not running (done:
    diagnosis and instrumentation).**

    ⚠️ **The headline of this item is wrong, and item 17 corrects it.** The water was not
    being created: the mass-weighted mean divides by an air mass that `densities()` rebuilds
    every iteration, and it was the air mass that moved. Measured against frozen weights the
    water field drifts +0.0011 % where this item reports +0.17 %. The instrumentation and
    the two experiments below stand; the attribution to the transport does not. Read this
    item for how the measurement was built and item 17 for what it turned out to measure.

    Item 14 left the `c ≤ 1 − co2` ceiling deleting water in 282 302 cells per iteration,
    with the working hypothesis that sedimentation pumps water out of the condensing level
    and nothing returns it. That hypothesis is wrong.

    `ThermoAtm::waterBudget()` now reports the global mass-weighted mean of
    vapour + cloud + ice + graupel, against its initial value. ATHAD has no water source
    and no water sink — the surface is supercritical, so `waterVapourEvaporation()` returns
    at once — so that mean is conserved exactly, and any drift is scheme error. Nothing had
    ever measured it. Measured over 20 iterations:

    ```
    iter  0:  q = 0.672751   drift +0.0000 %   deleted by the ceiling 0.000000
    iter 10:  q = 0.673903   drift +0.1712 %   deleted by the ceiling 0.007173
    iter 20:  q = 0.673893   drift +0.1698 %   deleted by the ceiling 0.014293
    ```

    Water is **created** at ~0.0007 kg/kg per iteration — about 0.11 % of the total per
    iteration — and the ceiling removes almost exactly as much as appears. The clamp is not
    starving the band; it is bailing out a leak.

    Two experiments locate it:

    - **`CategoryIceScheme = -1`** (no ice scheme at all): identical, +0.1718 % against
      +0.1712 %. Not the microphysics.
    - **`moist_phys_start_iter = 300`.** The moist physics — SaturationAdjustment, the ice
      schemes, MoistConvection — **does not run at all until iteration 300.** Every
      20-iteration diagnostic in this README was therefore a *dry* run, and the water is
      created with no moist physics executing. What remains is the RK4 tracer transport.

    **And the CO₂ conservation of item 12 does not contradict this — it is worthless as a
    test.** CO₂ drifts +0.0000 % because it is *uniform*, and a uniform tracer is conserved
    by any consistent advection scheme, however non-conservative, since `u·∇q = 0`. Water
    carries structure (put there by the one-off moist physics in the initialisation), and
    only a field with structure can expose the error. The earlier claim that the CO₂ result
    validated the transport was wrong.

    **Where the drift is.** Per-level attribution in `waterBudget()` puts it at levels
    42–47, **113 to 149 km** — deep in the supercritical column, far below both the cloud
    deck at 243 km and the band the `c` ceiling clamps at 185–230 km. It is static between
    iterations 10 and 20 (+2.83e-3 against +2.82e-3), so the redistribution happens early
    and then settles into a steady state in which creation balances the ceiling's deletion.
    Water moved *downward* out of the low-density band into air one to two orders of
    magnitude denser, and the mass-weighted total rose accordingly.

    **The cause is that the flow does not satisfy the continuity the tracer equation
    assumes.** `RHS_Atm_Turb` integrates `∂q/∂t = −u·∇q + …`, which conserves `∫ρq dV` only
    when `∇·(ρu) = 0`. The pressure projection enforces `∇·u = 0`, and ATHAD's density spans
    five orders of magnitude across the column.

    **The flux-form correction proposed here first is the wrong fix, and was not applied.**
    Expanding `∂(ρq)/∂t + ∇·(ρuq) = 0` gives `q[∂ρ/∂t + ∇·(ρu)] + ρ[∂q/∂t + u·∇q] = 0`, and
    the first bracket *is* continuity — so for a **mass fraction** the advective form is
    already exactly right, conditional on continuity. Adding `−(q/ρ)·∇·(ρu)` against a fixed
    `ρ` and a velocity field with `∇·(ρu) ≠ 0` would restore the global integral by making a
    *uniform* tracer develop structure: CO₂ would stop being well mixed, which is precisely
    the failure the uniformity test in the scope below exists to catch. It would trade a
    measured 0.17 % global error for an unphysical local one.

    So the defect is upstream, in the flow, and the fix belongs in `PressureSolverAtm`:
    anelastic continuity, `∇·(ρu) = 0`. This is the **Boussinesq risk in CLAUDE.md arriving
    in concrete form**. It changes the dynamical core and every result in this README, which
    is why it is scoped below rather than applied.

    Two further consequences worth knowing: the cloud deck discussed in items 9 and 14 is
    an *initialisation-time* state that the first 300 iterations only advect, and the
    ParaView output of a 400-iteration run has active moist physics only over its last 100
    iterations.

16. **The cloud came back: an undamped Newton step, and precipitation mass returned
    (done). And step 1 of the anelastic scope, measured and falsified.**

    The output showed water vapour everywhere and **no cloud, ice or rain anywhere in the
    domain** — while the column at 243 km sat ninefold supersaturated, `q_H2O` = 0.72
    against `q_sat` = 0.076, flagged `CONDENSING`. Two separate causes.

    **The precipitation mass was being destroyed.** `evaporateWhereImpossible()` (item 14)
    zeroed `P_rain` and `P_snow` without returning their mass to the vapour. The deck
    condensed, the ice scheme autoconverted it to rain, the rain fell one level into the
    superheated band, and it vanished. A downward flux `P` falling at speed `v` corresponds
    to a mass fraction `P/(v·ρ)` in the air it passes through; that is now converted before
    the flux is cleared, using the fall speeds the schemes' own residence times are built
    from (1.6 m/s rain, 0.96 m/s snow).

    **The first Newton step of the saturation adjustment was undamped**, and that is what
    actually suppressed the cloud. `adjustSaturation` initialised `q_v_hyp = q_sat`, a jump
    straight to the saturation value, and applied its `omega = 1/(1+Gain)` damping only
    from the second pass. On Earth this is harmless — condensing the ~0.01 kg/kg a
    terrestrial parcel holds releases about 12 K. Here the undamped step condenses
    **0.65 kg/kg in one go and releases 800 K** of latent heat. `T` is then clamped at the
    critical point, where `p_sat` = 220 bar against a local 0.085 bar, so `q_sat` becomes 1,
    the target inverts, and the next pass evaporates everything back. The iteration
    flip-flops between fully condensed and fully evaporated and finishes at zero.

    This one is not an inherited constant, it is an inherited **assumption**: that latent
    heating is a perturbation. At 67 % water by mass, condensation is a bulk phase change of
    the atmosphere. Starting the iteration from `q_v_b` makes the first pass compute a
    properly damped target. The equilibrium it should find is modest — condensing ~0.04
    kg/kg warms the parcel ~47 K, after which the vapour is superheated and nothing more
    can condense.

    **Measured after the fix** (20 iterations, moist physics on from iteration 1):

    ```
    max cloud water = 15.687 g/kg  @ 73°S, 242 774 m
    max cloud ice   = 14.400 g/kg  @ 68°S, 256 027 m
    ```

    A cloud deck at the condensation level, in both hemispheres, at the high latitudes
    where the column is coldest — and none anywhere it cannot exist.

    **Step 1 of the anelastic scope, measured.** `PressureSolverAtm` now reports the
    divergence of the *actual* velocity after the projection — nothing had ever measured
    it, so "the projection enforces ∇·u = 0" was an assumption about the code rather than
    an observation of it. A/B on the existing `ATM_POISSON_METRIC_FIX` knob:

    | | `div(u)` rms | water drift, 20 iters | ceiling deletions |
    |---|---|---|---|
    | metric fix off | 2.625e-02 | +0.1698 % | 0.014293 |
    | metric fix on | **2.153e-02** | +0.1698 % | 0.014293 |

    The consistent metric reduces the residual divergence by 18 % — real, and worth
    keeping — but the water drift is **bit-identical**. So the metric inconsistency is not
    the cause of the mass error, and step 1 of the scope is falsified as an explanation
    while remaining valid as a repair. Note also the absolute number: an rms `∇·u` of
    2.6e-02 is not a small residual. The projection is leaving a great deal of divergence
    behind, which is consistent with the anelastic diagnosis and makes steps 2–7 the
    remaining candidate.

17. **The anelastic projection, built and measured — and the water was never being created
    (done).**

    Steps 2–7 of the scope below are implemented behind `ATM_ANELASTIC` (default 0,
    bit-identical when unset): a one-dimensional base state `ρ̄(z)` rebuilt by
    `ThermoAtm::densities()` alongside the profile it averages; the divergence source
    `D = ∇·u* + u*_r·dlnρ̄/dr`; the same `dlnρ̄/dr` as a first-derivative term in the Poisson
    stencil, so source and operator stay adjoint — the lesson of step 1; and `ρ̄u_r = 0` at
    the surface and the lid in place of the `c43/c13` extrapolation, which permitted a
    through-wall mass flux. Step 7 needed no change: `t_ref_level[i]` is already the same
    horizontal mean, of the same prescribed profile, that `ρ̄` is.

    **It works, and it does not fix the water.** A/B at 20 iterations, metric fix on in
    both, `ATM_ANELASTIC` the only difference:

    | | off | on |
    |---|---|---|
    | `∇·u` rms | 2.153e-02 | 2.670e-02 |
    | `∇·(ρ̄u)/ρ̄` rms | 3.359e-02 | **2.607e-02** |
    | max radial wind after the initial projection | 0.1937 m/s | **0.0972 m/s** |
    | max meridional wind, iter 20 | 3.562967 m/s | 3.652312 m/s |
    | water drift, 20 iters | +0.0630 % | **+0.0630 %** |
    | drift by level, top 6 | i=44 +2.82e-03, i=45 +2.80e-03, … | **identical to 3 s.f.** |
    | CO₂ drift; min vs max | 0.0000 %; equal | 0.0000 %; equal |

    The anelastic residual falls 22 %, the Boussinesq one rises — it is no longer the
    enforced quantity — and half the spurious radial wind the old projection left in the
    initial field is gone. **Step 6 turns out to be unnecessary**: the `p_dyn_cap` source
    clamp, which the scope suspected of clipping the projection before it could act, binds
    in **0 of 3 791 399 fluid cells**. That is now printed every solve rather than assumed.

    And the water drift does not move — the third repair in a row to leave it identical.
    A tracer error indifferent to a velocity field that changed by a factor of two in `u`
    is not an advection error, which rules out the flow exactly as the algebra of item 15
    ruled out the tracer equation.

    **The drift is in the diagnostic's denominator.** `waterBudget()` reports water mass
    over air mass, and *both* come from `p_stat`, which `densities()` re-integrates
    hydrostatically every iteration on the local `R` and `cp` — which follow the
    composition and the surface temperature. Nothing separated the two. Measured against
    weights frozen at the reference time, with the column air mass those weights carry
    reported beside it:

    ```
    iter 10:  q_mean drift +0.0621 %   q against FROZEN weights +0.0011 %   column air mass -0.1258 %
    iter 20:  q_mean drift +0.0630 %   q against FROZEN weights +0.0025 %   column air mass -0.2128 %
    ```

    The water field moved by **+0.0025 % over 20 iterations**, twenty-five times less than
    the number this README has been quoting. What moved is the air: the column is losing
    about 0.01 % of its mass per iteration, steadily and without sign of stopping, so
    water-over-air rises. The per-level attribution at 113–149 km is the same artefact —
    those are the levels where `dp` changed most, not where water arrived.

    So **item 15's headline is wrong and is corrected here**: water was not being created
    at 0.11 % per iteration. The transport error is ~0.0001 % per iteration, and what the
    budget was measuring is the hydrostatic column being re-weighed. Two consequences:

    - The **column air mass is not conserved**, because nothing in this model makes it a
      conserved quantity. `p_stat.x[0]` is re-anchored every iteration to
      `r_air·R_mix·T_surf`, so the mass of a "250 bar atmosphere" follows the surface
      temperature — which the model evolves (1496.6 K against the prescribed 1500) — and
      the whole column is then re-integrated on the new anchor. It is a prescription
      defect, not a transport defect, and it is the open one.
    - The `c ≤ 1 − co2` ceiling is still deleting water (0.000164 kg/kg over 20 iterations),
      and that deletion is real. It is now the larger of the two mass errors.

    **A fourth repair, correct and not the cause.** `ATM_TRACER_DIFF_FLUX` (default 0) adds
    the missing term of the diffusive flux: for a mass fraction the conservative form is
    `∂(ρq)/∂t = ∇·(ρK∇q)`, i.e. `∂q/∂t = K∇²q + K∇lnρ̄·∇q`, and the second term was absent —
    the tracer diffusion was conserving `∫q dV` rather than `∫ρq dV`, in the same way the
    advective form does but *without* continuity to repair it. Unlike the advective
    flux-form correction of item 15 this one is proportional to `∇q`, so a uniform tracer
    stays uniform and CO₂ stays well mixed. Worth having and kept; measured effect at 20
    iterations is the fifth decimal of `residuum_atm` and nothing else. Note also that
    `diff_prec_re_inv` is **not** multiplied by `diffusion_ramp`, so moisture diffusion runs
    at full strength through a spin-up in which heat and momentum diffusion are ramped from
    zero. That asymmetry is undocumented and probably unintended.

    **`ATM_ANELASTIC` is left off by default.** It is the right continuity for this column
    and it measurably improves the projection, but 20 iterations is not evidence of
    stability, and the scope's own warning — every result in this README moves with it —
    stands. Flip it after a 400-iteration run, the way the metric fix was flipped.

## Scope: the anelastic continuity fix

**Status: steps 1–7 are implemented and measured (items 16 and 17). The premise below —
that the tracer mass error comes from the velocity field — is false; the error was in the
budget's denominator, not in the flow. The scope is kept because the anelastic projection
is right on its own terms and the reasoning is what the measurements were made against.**

The tracer mass error of item 15 comes from a velocity field that does not satisfy
continuity for the prescribed density. `PressureSolverAtm` projects onto `∇·u = 0`; a
column whose density spans five orders of magnitude needs `∇·(ρ̄u) = 0`. This is the scope
of that change, written before touching it because it is the dynamical core and every
result in this README moves with it.

**What the solver does now.** `PressureSolverAtm::run()` builds a provisional velocity
(`aux_u/v/w`), takes its divergence with single-power metric factors, solves a Poisson
equation for `p_dyn` with a variable-coefficient Laplacian, and corrects the velocity by
the pressure gradient. The file documents its own defect: the Laplacian coefficients use
`exp_2_rm`, `inv_rm2`, `inv_rm2sinthe2` while the divergence source and the gradient
correction are single-power, so `div` and `grad` are **not discretely adjoint** and the
projection does not exactly remove the divergence it measured. Three stabilisers
(`p_dyn_cap = 2.0`, the `p_dyn_ceiling`, the topography Dirichlet pins) were calibrated
against that inconsistency.

**The change, in order of dependency.**

1. **Resolve the metric inconsistency first.** Until `div` and `grad` are adjoint, no
   projection — Boussinesq or anelastic — removes the divergence. This is the file's own
   open TODO and it is independently testable: measure `∇·u` before and after the
   projection and require it to drop. *Do this step alone first and re-measure the water
   drift; it may account for much of it.*
2. **Introduce a base-state density `ρ̄(z)`** — one-dimensional, time-independent, the
   horizontal average of the prescribed profile. It must be a function of height only:
   a fully three-dimensional ρ makes the elliptic operator time-varying and costs
   solvability. Rebuild it when `densities()` rebuilds the profile.
3. **Anelastic divergence source.** `D = (1/ρ̄)∇·(ρ̄u*) = ∇·u* + u*_r · dln ρ̄/dr`. Because
   ρ̄ depends on r only, this is one extra term on the radial component.
4. **Anelastic Poisson operator.** `∇·(ρ̄∇φ) = ρ̄∇²φ + (dρ̄/dr)(∂φ/∂r)`, so the existing
   7-point stencil gains a first-derivative term in `r` proportional to `dln ρ̄/dr`. No new
   solver is needed — the change is to the stencil coefficients.
5. **Boundary conditions.** Impose `ρ̄u_r = 0` at the surface and the lid. `aux_u` is
   currently `c43/c13`-extrapolated at both ends, which does not impose zero normal mass
   flux.
6. **Re-derive the stabilisers.** `p_dyn_cap` and the ceiling were tuned against the
   inconsistent operator and will otherwise clip the anelastic projection before it acts.
7. **Check the buoyancy reference.** `BuoyancyForce = coeff_buoy·(t − t_ref_level[i])` is
   already base-state-referenced; confirm the reference is the same ρ̄.

**Verification — the point of the diagnostics already in place.**

| test | required | result (item 17, `ATM_ANELASTIC=1`, 20 iters) |
|---|---|---|
| `waterBudget()` drift, 20 iters | ≈ 0 | +0.0630 %, **unchanged** — and the wrong test: see item 17 |
| water drift against frozen weights | ≈ 0 | +0.0025 %, the number that was wanted all along |
| water-vapour ceiling hits | ≈ 0 | still deleting 0.000164 kg/kg / 20 iters |
| CO₂ max − min (uniformity) | still 0 | 0 exactly |
| `∇·(ρ̄u)/ρ̄` after projection | measured, and small | 3.359e-02 → **2.607e-02** |
| `∇·u` after projection | measured | 2.153e-02 → 2.670e-02 (no longer the enforced one) |
| source clamped at `p_dyn_cap` | not clipping the projection | **0 of 3 791 399 cells** |
| Ψ_max over 400 iters | bounded | not yet run at 400 |
| bit-identical at 1/4/8 threads | still passes | not re-checked |

The CO₂ uniformity test is the one that catches an over-correction: any scheme that makes a
well-mixed tracer develop structure is wrong, whatever it does for the mass budget. It is
what ruled out the advective flux-form correction, and it passes here because the anelastic
change is in the flow rather than in the tracer equations.

**Effort and risk.** The code is modest — a few dozen lines across `PressureSolverAtm` and
the divergence source. The verification is the work, and every number in this README
changes. Step 1 is separable, independently valuable, and should be measured before
steps 2–7 are started.

That estimate held: about forty lines, and the verification took four 20-iteration runs.
What it did not anticipate is that all of it would be measured against a diagnostic whose
denominator was moving. Three repairs — step 1, steps 2–7, and the diffusive flux term —
each left the water drift identical, and *that* is what finally identified the diagnostic
rather than any one of them. A repair that changes nothing measurable is evidence about
the measurement.

## Remaining work







- **The temperature profile is still prescribed, not solved.** `ThermoAtm::densities()`
  re-imposes the adiabat + isothermal top on `t` every iteration, so what the dynamics and
  the radiation compute is overwritten before it can matter. Invariant 3 in CLAUDE.md says
  radiation must *set* the profile rather than nudge it toward a prescribed one; the
  prescription still wins. The OLR is now a genuine integral **over a prescribed profile** —
  a real improvement, but not yet a prediction.
- **The column air mass is not conserved** (item 17). `p_stat.x[0]` is re-anchored every
  iteration to `r_air·R_mix·T_surf`, so a "250 bar atmosphere" loses mass whenever the
  surface temperature moves — 0.126 % over 20 iterations, which is what the water budget
  was reading as water creation. The tracer transport itself drifts ~0.0001 % per
  iteration. Fixing this means anchoring the column to a mass rather than to a prescribed
  surface density, and it is entangled with the prescribed profile below.
- **The `c ≤ 1 − co2` ceiling still deletes water**, 0.000164 kg/kg over 20 iterations.
  With the transport exonerated this is now the larger of the two mass errors, and unlike
  the other one it is a straightforward deletion with no compensating source.
- **`moist_phys_start_iter = 300`** means a 400-iteration run is dry for three quarters of
  its length. Deliberate (it lets the circulation form before the stiff microphysics
  starts), but it must be stated whenever a run is quoted.
- **The run is not converged**: 100 iterations leaves a −143 W/m² imbalance, still decaying
  ~9 % per 10 iterations. Of order 400 iterations are needed. Run it.
- **The OLR is not grid-converged either**: 519 W/m² at a 260 km shell against 581 at
  300 km, with `im` fixed at 61. Refine vertically and check.
- **The opacity is too low to hold the surface.** OLR 581 W/m² against 271 absorbed. Since
  `kappa_H2O`/`kappa_CO2`/`kappa_bg` carry a factor-of-two uncertainty and are the biggest
  lever, this is the first quantity worth testing against the new scheme.


- **`albedo_cloud = 0.50` IS the model's planetary albedo** and is an assumption. It is
  now the second-biggest lever after the opacities. A deep, cold, slowly sedimenting
  Hadean deck could plausibly be brighter.
- **`geothermal_flux` is the open number.** The ≥ ~195 W/m² figure was derived with the
  clear-sky albedo and the too-low insolation, so it has to be redone once the radiation
  can produce an OLR worth comparing to. Check against magma-ocean cooling estimates.
- **`initTemperatureData` builds its column before the composition exists** — no water, no
  CO₂, so a background-only cp and too steep a lapse. `densities()` overwrites it, so the
  only thing that ever escaped was the lid snapshot (item 10), which is now refreshed. It
  should still be built on the real mixture.
- **`initCloudIce`'s H_crit parabola is keyed to absolute pressure** (`p_crit = 1000` hPa,
  `p_mid = 550`), an Earth surface pressure. It should be a fraction of the local surface
  pressure, like the deep-convection triggers.
- **The three κ opacities** carry factor-of-two uncertainty and are the biggest lever on
  the OLR. A grey scheme cannot represent the window regions that set the real limit.
- **The surface temperature is prescribed, not solved.** Everything above is conditional
  on that.
- **Boussinesq** remains untested against a column whose density spans two orders of
  magnitude (see CLAUDE.md).

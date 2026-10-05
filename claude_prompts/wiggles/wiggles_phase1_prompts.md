# WATERHYPERNET -- Wiggles Phase 1: Twin experiment

## Goal

Know, before touching real data, which interpolation method removes which
wiggle and what each does to real spectral features.  A high-resolution
irradiance (TSIS-1 HSRS × HAPI O₂/H₂O × O₃ on a 0.01 nm grid) multiplied by
smooth OSOAA fields is observed through a HYPSTAR instrument model with the
VEIT FWHM(λ) curves and the real E/L grids, then ρw is recomputed with
`linear`, cubic/sinc (null), `ruddick2023` and `srf` and judged by ρw''.

**Gate G1** (plan §4): the `srf` method reduces the line-region ρw'' error to
the away-from-lines level (target > 80 % reduction in cases (ii)-(iii)),
degrades gracefully under (iv)-(vi), and leaves the fluorescence and Raman
controls unchanged to within the noise floor of case (vii).  If `ruddick2023`
already achieves this in case (iii) at the VEIT SRF difference, H2 is a
refinement rather than a requirement and Phase 2 defaults to eq. (14).

## Conventions

- `ocean14`.  Run via `conda run -n ocean14 python ...`; `conda activate` fails
  non-interactively.
- **Model.** Use Opus 5.5 (`claude-opus-5-5`) for this work, including any
  subagents (pass `model: opus`).
- **JXP runs git.**  Claude does not run any state-changing git command -- in
  this repo or in the OSOAA fork.
- Run scripts **from the repository root** so `hypernet` imports resolve.
- Run `pytest -q` after each step where relevant.  58 tests on 2026-10-03; the
  archive-dependent ones skip themselves when `$OS_COLOR` is not mounted.  OSOAA-
  dependent tests must skip when `OSOAA_ROOT` is unset or the exe is missing.
- **Tier 2** -- the VEIT grids and noise levels come from `$OS_COLOR`.  Do
  **not** unset `$OS_COLOR`.
- Any calculation goes into a script on disk, not into the chat.
- Data intermediates go **outside** the repo: reference spectra and line
  databases under `$OS_COLOR/hypernet/wiggles/ref/`, OSOAA runs and twin-
  experiment outputs under `$OS_COLOR/hypernet/wiggles/phase1/`.  Only
  figures and small tables are committed.
- Every `.md` in `docs/` is published.  Gate reports stay out of `docs/`.
- OSOAA runs take 30 s to 15 min per wavelength: run grids in the background
  and start with one case before launching the full set.
- Ask questions in the Q&A section below; log completed work under `## Logs`.

## Context

### Where things are

- **Plan:** `docs/wiggles_planning.md` §3 (the three methods and the
  generalised eq. 14), §4 "Phase 1", §6 (reference spectra and RT tools),
  §9 (Emod line shapes, angular effects, Ring).
- **Setup doc:** `claude_prompts/wiggles/wiggles_prompts.md` -- Q&A rounds 2-3
  and Logs on OSOAA (unbuilt, gas-free, monochromatic; wrapper copies in the
  fork's notebooks and tests; no `LUM_Advanced_Down.txt` parser; Level
  numbering conflict) and HAPI.
- **OSOAA fork:** `/Users/xavier/Oceanography/python/RadiativeTransferCode-OSOAA`
  -- `CLAUDE.md` (build: `cp gen/Makefile_OSOAA.gfortran Makefile; make all`
  → `exe/OSOAA_MAIN.exe`), `src/OSOAA_MAIN.F`, `doc/OSOAA-V2.0_UserManual-V2.0.pdf`,
  `docs/user_guide/output_files.rst` (levels listed as TOA, sea bottom, 0+,
  0−, user-defined; `input_parameters.rst` gives the example `-1`),
  `tests/test_water_Ed.py` (the `OSOAASimulation` class: `get_default_params`
  with `'OSOAA.View.Level': 4  # 0-`, `build_command`, `run`,
  `_parse_flux_file`), `notebooks/*.ipynb`.  Homebrew gfortran 16.2.0 at
  `/opt/homebrew/bin/gfortran`.  **Built 2026-10-03** (task 2):
  `exe/OSOAA_MAIN.exe`, ~10 s per wavelength.  `OSOAA.View.Level`: 1 = TOA,
  2 = sea bed, 3 = 0+, 4 = 0−, 5 = user `-OSOAA.View.Z` (the fork's docs
  saying ±1 are wrong).
- **HAPI:** `/Users/xavier/Oceanography/python/HAPI`, installed editable in
  `ocean14` (hitran-api 1.3.0.0).  HITRAN downloads need network access; cache
  the tables under `$OS_COLOR/hypernet/wiggles/ref/hitran/`.
- **Reference spectra to fetch:** TSIS-1 HSRS v2 (LASP; highest native
  resolution), Serdyuchenko et al. 2014 O₃ cross-sections.  Record URLs and
  checksums in the Logs.  *HSRS is already fetched (Phase 0 task 3b,
  2026-09-29):* `hypernet/refspec.py` (`fetch_hsrs`, `load_hsrs(wmin, wmax,
  step, frame='vac'|'air')`, URL and SHA-256 in the module) and the file
  `$OS_COLOR/hypernet/ref/hybrid_reference_spectrum_p005nm_resolution_c2022-11-30_with_unc.nc`
  (0.005 nm resolution, 0.001 nm sampling, vacuum wavelengths).  The Emod
  builder should import it, not download again.
- **From Phase 0a:** `hypernet/data/veit_srf_model.json` (E, Ld, Lu FWHM(λ)
  models from the HSRS template fit of task 3b: the SRF widths themselves,
  E ≈ 2.2-2.4 nm, Ld ≈ 2.9 → 2.6 nm over 400-700 nm; valid 390-680 and
  850-870 nm, unconstrained 680-850 nm; read with `srf.load_srf_models`;
  the empirical line-width models are in `veit_srf_model_empirical.json`),
  the per-scan noise levels (Phase 0 task 6), the module names chosen there
  (`hypernet/srf.py`, `hypernet/whn_l1a.py`, `hypernet/whn_srf.py`,
  `hypernet/refspec.py`).  Ld is the reference L channel (Phase 0 Q11).  The real E (1536 px) and L (1538 px) grids come from the
  VEIT sample at `$OS_COLOR/WATERHYPERNET/Wavelengths/`; HYPSTAR geometry
  there (vza is measured from nadir): Lu at `viewing_zenith_angle` ≈ 40°
  (water, 40° off nadir), Ld at ≈ 140° (sky, 40° off zenith), SZA 36.7°.
- **Paper:** `context/papers/ruddick2023.pdf` (eqs. 11-15, §3 simulations).

### Dependencies

- No RBINS data needed.  Tasks 2-5 need only the fork, HAPI and the network.
- Tasks 7-9 need the Phase 0a SRF model: it exists, as the template SRFs in
  `hypernet/data/veit_srf_model.json` (122304, narrow-E).  Per Q&A Q6, task
  7 also builds `hypernet/data/release2_srf_models.json` from Phase 0 task
  8c so that every case runs for two instruments: 122304 (narrow-E,
  FWHM_Ld − FWHM_E ≈ 0.5 nm) and 122305 (E ≈ L).  See the Phase 0 Q&A
  "Build #8c" and "Build #8f".
- **The one-function interpolator** (`interpolate_ed_to_l`) is needed here
  and hardened in Phase 2.  Recommendation: write a minimal, tested version
  in task 8 of this doc and let Phase 2 add validation, caching and
  uncertainty.  JXP confirms in Setup Q&A, together with its module name.

## Prompts

### Setup

1. Read the plan sections and setup-doc material listed in Context, this doc,
   the fork's `CLAUDE.md` and `tests/test_water_Ed.py`, and the two tutorial
   notebooks.  Put your questions in Q&A: the interpolator's home and name,
   the module names for the Emod builder and the twin-experiment scene code
   (defaults `hypernet/emod.py`, `hypernet/twin.py`), whether OSOAA is built
   in the fork in place, and which aerosol and water cases to run first.  Do
   not build or run anything yet.  Log your work.

### Build

2. **Task 0: build OSOAA.**  In the fork: `cp gen/Makefile_OSOAA.gfortran
   Makefile; make all`, adding `-std=legacy` if gfortran 16.2 rejects the
   fixed-form sources.  Set `OSOAA_ROOT` (document how in the Logs).  **First
   settle `OSOAA.View.Level`:** read its handling in `src/OSOAA_MAIN.F` and
   the user manual, and record the mapping (the fork's docs say ±1, the
   wrappers use 1 = TOA and 4 = 0−; we need 0+).  Smoke-test with the fork's
   `tests/test_water_Ed.py`.  Do not commit anything in the fork.  Log your
   work, including build fixes and timings per wavelength.

3. **`hypernet/rt/osoaa.py`.**  A merged wrapper, imported from nowhere in the
   fork: `osoaa_root()`, `default_params(wavelength_nm, sza, ...)`,
   `build_command(params)`, `run(params, work_dir)` (subprocess, one
   wavelength per call), parsers `parse_flux(path)` (direct/diffuse/total Ed
   per level), `parse_lum_vsvza(path)` (upward radiance vs VZA at the chosen
   level) and the **new** `parse_lum_advanced_down(path)` (downward radiance
   vs VZA and azimuth, for Ld).  Save one real run's output files as fixtures
   under `hypernet/tests/data/osoaa/` and test the three parsers on them;
   add a `needs_osoaa` skip marker for one live 550 nm run.  Log your work.

4. **Emod builder.**  Module per Q&A (default `hypernet/emod.py`):
   HSRS via `hypernet.refspec.load_hsrs` (already downloaded), `gas_transmittance(lam,
   airmass, pwv_mm, pressure_hpa, temperature_k)` using HAPI `absorptionCoefficient_Voigt`
   for O₂ and H₂O on the 0.01 nm grid 380-1000 nm, `ozone_transmittance(lam,
   airmass, ozone_du)` from Serdyuchenko, and `build_emod(sza, pwv_mm=15,
   ozone_du=300, pressure_hpa=1013)` returning λ, F₀, T_direct (air mass for
   the SZA) and T_diffuse (an effective diffuse air mass, documented).  Cache
   results as npz under `ref/emod/`.  Tests: transmittances in [0, 1]; O2-A
   depth at 760.6 nm within a stated range; grid spacing 0.01 nm; cache hit
   reproduces the computed array.  Log your work.

5. **OSOAA smooth fields.**  `hypernet/wiggles/phase1_osoaa_fields.py`: on a 5 nm grid
   380-1000 nm, for SZA 30/50/70°, two aerosol loads and a small
   chlorophyll/sediment grid, run OSOAA and collect direct and diffuse Ed at
   0+, Ld at the HYPSTAR sky geometry (VZA 40°, relative azimuth 90°: the VEIT
   sample and all 224 requested sequences) and Lw at 0+.  Run one SZA/aerosol/water case first, time it,
   then the rest in the background.  Save
   `$OS_COLOR/hypernet/wiggles/phase1/osoaa_fields.npz` and a figure
   `hypernet/wiggles/figs/phase1/osoaa_fields.png`.  As a check, compare
   the modelled clear-sky Ld(750)/Ed(750) with the observed clear floor of
   the sky index (0.011-0.015 sr⁻¹ at every SZA; Phase 0 Context → Sky
   index).  Log your work.

6. **Scene composer and ρw library.**  Module per Q&A (default
   `hypernet/twin.py`): `compose_scene(emod, fields, case)` returning
   high-resolution Ed (direct + diffuse, separate air masses), Ld and Lw;
   `rhow_library(fields)` from the OSOAA water grid, with a chlorophyll
   fluorescence control (Gaussian at 683 nm, ~25 nm FWHM, amplitude scaled to
   the clear-water and turbid cases) and a water-Raman control (Ed
   redistributed by the ~3357 cm⁻¹ shift with a simple Raman efficiency);
   `rhow_true(scene, srf_L, grid_L)` -- the truth is Lw and Ed both observed
   with the L SRF on the L grid, then π Lw / Ed.  Tests: flat Emod and flat
   fields give flat ρw; the fluorescence control adds a peak at 683 nm.  Log
   your work.

7. **Instrument model.**  In the same module: `observe(lam_hr, spec_hr, grid,
   srf_model, seed=None, noise=None)` -- Gaussian SRF with FWHM(λ) from
   `hypernet/data/veit_srf_model.json`, the real 1536/1538 grids saved once to
   `hypernet/data/hypstar_grids.npz`, and a `Case` table for (i) grid offset
   only, equal SRFs; (ii) SRF-width mismatch, no offset; (iii) both; (iv)
   wavelength error ±0.1 and ±0.3 nm; (v) wrong SRF width in the correction
   ±0.3 nm; (vi) wrong Emod (SZA off by 10°, water vapour ×2, no aerosol);
   (vii) photon noise at the VEIT per-scan level.  Tests: the SRF integrates
   to 1; a delta line observed with FWHM 3 nm returns FWHM 3 nm through
   `hypernet.srf.fit_line`.  Log your work.

8. **Minimal interpolator.**  In the module agreed in Q&A:
   `interpolate_ed_to_l(wav_irr, irradiance, wav_rad, *, emod=None,
   srf_irr=None, srf_rad=None, method="srf")` implementing plan §3 for
   `method in {"linear", "ruddick2023", "srf"}`, plus `"cubic"` (CubicSpline)
   and `"sinc"` as nulls.  Tests in `hypernet/tests/test_edinterp.py`: with
   constant Emod every method equals `np.interp` to machine precision; with
   `srf_rad = srf_irr` the `srf` output equals a direct implementation of
   eq. (14) to machine precision; a spectrum equal to Emod_E returns Emod_L.
   Log your work.

### Analysis

9. **Cases × methods.**  `hypernet/wiggles/phase1_twin.py`: for every case, SZA,
   aerosol and water type, and every method, compute ρw on the L grid, then
   ρw'' with h = 1 nm and h = 5 nm second differences as in the paper.
   Metric: rms of ρw'' − ρw''_true within ±5 nm of the ten lines and away
   from them.  Write `$OS_COLOR/hypernet/wiggles/phase1/twin_metrics.parquet`,
   a committed `hypernet/wiggles/phase1_twin_metrics.csv`, and figures under
   `hypernet/wiggles/figs/phase1/` (ρw'' at O2-A, Hα and Ca H/K per method; the
   reduction per case).  Log your work.

10. **Controls and sensitivity.**  Extend the script: the change in the
    fluorescence and Raman controls' second derivative per method against the
    case (vii) noise floor; degradation curves for (iv), (v) and (vi).
    Table `hypernet/wiggles/phase1_controls.csv` and a figure.  Log your work.

11. **Prediction for Phase 3.**  From case (iii), tabulate the expected
    reduction of the line-region ρw'' excess for `ruddick2023` and `srf` as a
    function of ΔFWHM = FWHM_L − FWHM_E (0 to 1 nm) and line depth -- the
    number Gate G3 will be judged against.  `hypernet/wiggles/phase1_prediction.csv`
    and a figure.  Log your work.

### Gate

12. **Gate G1.**  Evaluate the criteria with the numbers from tasks 9-11:
    the reduction in cases (ii)-(iii) against 80 %, graceful degradation in
    (iv)-(vi), the controls against the (vii) floor, and whether
    `ruddick2023` alone passes case (iii) at the VEIT SRF difference.  Write
    a short gate report (proposed location `claude_prompts/wiggles/gate_G1.md`;
    JXP to confirm -- not under `docs/`) with the verdict, the Phase 2 default
    method and open concerns.  Propose the corresponding edits to
    `docs/wiggles_planning.md` §4 in Q&A and apply them only if JXP approves.
    Log your work.

## Q&A

### Setup #1 -- 2026-10-03 (Opus 5.5)

**What has changed since this doc was written** (from Phase 0, tasks 3b,
8c, 8e):

- **Two instrument types.**  On the 224 Release 2 sequences:
  - 122302, 122304 and 120242 have an E SRF ~0.5 nm narrower than L
    ("narrow-E"; E ≈ 2.2-2.4 nm, Ld ≈ 2.7-2.9 nm).
  - 121222, 121231, 122303 and 122305 have E ≈ L (within 0.07 nm).
  - Both types are stable to 0.01-0.06 nm over SZA, season and
    recalibration.
- **ρw follows that split.**  The narrow-E instruments show the
  line-selective ρw excess H2 predicts, at about the predicted amplitude.
  The E ≈ L instruments still carry an unexplained 0.5-2 % line-residual
  floor in ρw.
- **The wavelength scale is air, within 0.07 nm, on every instrument**, so
  case (iv)'s ±0.1/±0.3 nm errors bracket reality generously.
- **Kevin's sky index is in** (Phase 0 Context).  Task 5 now uses RAA 90°
  and checks against its clear floor.
- **Stale here:** "21 tests today" (now 58), and "until [the Phase 0a SRF
  model] exists use the §2.3 table" (it exists: template SRFs in
  `hypernet/data/veit_srf_model.json`).

**What I learned from the fork** (`tests/test_water_Ed.py`, both tutorial
notebooks; nothing run):

- One `OSOAASimulation` class: a keyword dict → `OSOAA_MAIN.exe -KEY value`,
  run in a work dir.
  - `Flux.txt` (in `Advanced_outputs/`) holds direct/diffuse/total Ed and Eu
    at all 108 levels, normalised to π at TOA.  So one run gives Ed at 0+
    whatever `View.Level` is.
  - `LUM_vsVZA.txt` gives the radiance (Stokes I = πL/Esun) vs VZA at the one
    level set by `OSOAA.View.Level`.  Both notebooks use 4 = 0−.
  - The water-leaving notebook also writes `LUM_Advanced.txt`; nothing
    parses it, nor the downward (sky) radiance we need for Ld.
- **A bug to fix before task 5: `'AER.Waref': wavelength_um`.**  The aerosol
  reference wavelength is set to the simulated wavelength, so AOT equals
  `AOTref` at *every* λ and the aerosol's spectral shape (its Ångström
  behaviour) is lost.  For spectral fields, Waref must be fixed (e.g.
  0.55 µm) and OSOAA left to scale AOT(λ) with its Mie model.
- **Lw at 0+ is not what the notebooks compute.**  Level 4 is Lu(0−).  Either
  transmit it (Lw = Lu(0−)·(1 − ρ_F)/n²), or run Level 3 (0+) twice, with
  and without the water body, and difference them to remove the surface
  reflection.
- **Run count.**  Each OSOAA call is one wavelength and one radiance level.
  Ld (downward at 0+) and Lu (0−) may need separate runs, unless the
  Advanced output holds both directions at the chosen level (task 2 checks).
  At 30 s-15 min per wavelength, the task 5 grid (3 SZA × 2 aerosol ×
  water cases × 125 wavelengths × 1-2 levels) could take days.

**Questions for JXP** (each has a default; say "defaults" to accept them all)

1. **The interpolator's home and name.**  *Default:
   `interpolate_ed_to_l` in a flat module, `hypernet/edinterp.py`, next to
   `srf.py` and `refspec.py`.  Task 8 writes the minimal tested version;
   Phase 2 hardens it (validation, caching, punpy).*
>A. Use your default
2. **Module names.**  *Default: `hypernet/emod.py` (Emod builder),
   `hypernet/twin.py` (scene composer, ρw library, instrument model and the
   `Case` table), and `hypernet/rt/osoaa.py` (a new `hypernet/rt/`
   subpackage).  Phase 1 scripts go in `hypernet/wiggles/phase1_*.py`,
   figures in `hypernet/wiggles/figs/phase1/`.*
>A. Use your default
3. **Building OSOAA.**  *Default:*
   - *Build in place in the fork (`cp gen/Makefile_OSOAA.gfortran Makefile;
     make all`, adding `-std=legacy` if needed), commit nothing there, and
     leave the copied Makefile and the build products untracked.*
   - *`hypernet.rt.osoaa.osoaa_root()` reads `$OSOAA_ROOT`, falling back to
     `/Users/xavier/Oceanography/python/RadiativeTransferCode-OSOAA`.  I will
     not edit your shell profile; add `export OSOAA_ROOT=...` yourself if you
     want it set globally.*
>A. Use your default
4. **Which cases first, and how big a grid.**  *Default:*
   - *Start with one VEIT-like case: SZA 40° (the sample is 36.7°), RAA 90°,
     AOT(550) = 0.1 with `AER.Waref` fixed at 0.55 µm (the fix above), and
     Chl 1 mg m⁻³.*
   - *Time it per wavelength.*
   - *Then the full grid: SZA 30/50/70 × AOT(550) {0.05, 0.25} × water
     {Chl 0.1, 1, 10 mg m⁻³; plus one turbid case with sediment}.*
   - *If a wavelength takes > 1 min, use a 10 nm grid.  The fields are
     smooth, and a 10 nm spline loses nothing at 3 nm resolution.  Run in
     the background, as the doc says.*
>A. Use your default
5. **Lw at 0+.**  *Default: transmit Lu(0−) to 0+ with
   (1 − ρ_F(θ))/n², one Level-4 run per wavelength.  It is smooth in λ,
   which is all the twin experiment needs.  The difference method (two
   Level-3 runs) is kept as a check at a few wavelengths.*
>A. Use your default
6. **Two instruments, not one.**  The doc models only the VEIT (narrow-E)
   SRFs.  *Default:*
   - *Add a small step to task 7 that writes
     `hypernet/data/release2_srf_models.json`: per instrument × calibration
     period, the median template FWHM(λ) for E and Ld over the Phase 0c
     sequences.*
   - *Run every case for two instruments: 122304 (narrow-E, ΔFWHM ≈ 0.5 nm)
     and 122305 (E ≈ L).*
   - *Case (i) on 122305 then directly asks whether linear interpolation
     (H1) explains the 0.5-2 % floor seen on E ≈ L instruments.*
   - *Task 11's ΔFWHM sweep (0-1 nm) is unchanged.*
>A. Use your default
7. **O₃ cross-sections.**  *Default: Serdyuchenko et al. (2014), from the
   IUP Bremen spectroscopy database, at 293 K, cached under
   `$OS_COLOR/hypernet/wiggles/ref/`, with the URL and checksum in the Logs.
   Fallback: HITRAN's UV O₃ cross-sections via HAPI.*
>A. Use your default
8. **Housekeeping.**  *Default: I update this doc's stale lines (the test
   count, the "until it exists" note, a pointer to the two instrument
   types).  In the Phase 0 doc, the two tasks numbered 8e (mine, "ρw
   wiggles by instrument class", done; and the Sky-index task added in the
   Setup #6 session) need distinct numbers.  I propose renaming mine to
   **8f** and updating its Q&A/Log headers and the script docstring, once
   that session is no longer editing the file.*
>A. Use your default

### Gate #12 -- 2026-10-04 (Opus 5.5)

The verdict is in `claude_prompts/wiggles/gate_G1.md`: **pass for narrow-E
instruments; Phase 2 default `srf` above ΔFWHM 0.15 nm, `linear` below; H2 is
required** (`ruddick2023` removes 0 % of the mismatch).

**Questions for JXP** (each has a default; say "defaults" to accept them all)

1. **Report location.**  *Default: keep `claude_prompts/wiggles/gate_G1.md`
   (outside `docs/`, where the Phase 2 doc already expects it).*
>A. Yes, that is fine
2. **Plan §4 edits.**  *Default: apply these four to
   `docs/wiggles_planning.md` once you approve:*
   - (a) *Phase 1, the G1 paragraph: replace "leaves the fluorescence and
     Raman controls unchanged to within the noise floor of case (vii)" with
     "changes the fluorescence and Raman controls' ρw'' by no more than the
     linear method does, and by ≪ the controls themselves in case (iii)
     (the controls' ρw'' is itself below the single-sequence noise)".*
   - (b) *Phase 1, after the G1 paragraph, add the outcome: "**G1 outcome
     (2026-10-04):** passed for narrow-E instruments.  `srf` removes ~100 %
     of the mismatch error in (ii)-(iii); `ruddick2023` removes 0 %, so H2 is
     required.  Tolerances for 80 %: E wavelength scale ≤ 0.05 nm relative to
     Emod, FWHM_E ≲ 0.07 nm, FWHM_L ≲ 0.1 nm; Emod errors are harmless.  On
     E ≈ L instruments the correction gains nothing and is fragile.  See
     `claude_prompts/wiggles/gate_G1.md`."*
   - (c) *Phase 2: add two bullets.*
     - *"Default method `srf`, applied where the measured ΔFWHM ≥ 0.15 nm;
       `linear` otherwise."*
     - *"E wavelength scale: fit a shift and linear stretch of E against
       the HSRS-based Emod (per sequence or per calibration period) to
       ≤ 0.05 nm, and carry its uncertainty."*
   - (d) *Phase 3, the G3 paragraph: after "by the amount Phase 1 predicts
     for the measured SRF difference", add "(the noise row of
     `hypernet/wiggles/phase1_prediction.csv`, observable excess
     √(rms²_line − rms²_away) of ρw'' at h = 1 nm; e.g. 60-75 % at
     ΔFWHM 0.5 nm, ≈ 0 below 0.1 nm)".*
>A. That is all ok
3. **The threshold.**  *Default: 0.15 nm, from the twin's break-even under
   realistic knowledge (gate report, "Why 0.15 nm").  Phase 2 revisits it
   once the L1A E SRFs and the punpy uncertainties exist.*
>A. Ok
4. **A non-Gaussian SRF check** (gate concern 3).  *Default: add it as a
   small Phase 2 task (twin truth with the Phase 0 pvoigt/supergauss shape,
   Gaussian correction), not a new Phase 1 task.*
>A. Agreed

## Logs

### 2026-10-03 -- Setup #1 (Opus 5.5)

- Read:
  - plan §3, §4 (Phase 1), §6 and §9;
  - the setup doc's OSOAA/HAPI rounds and the new Sky-index round
    (Setup #6, by another session);
  - this doc, and the fork's `CLAUDE.md` and `tests/test_water_Ed.py` (the
    `OSOAASimulation` class: default params, `build_command`, `run`, the
    `Flux.txt` parser);
  - both tutorial notebooks, through their code cells and outputs.
- Nothing was built or run.
- Learned:
  - `Flux.txt` covers every level, so Ed at 0+ comes free with any run, but
    `LUM_vsVZA.txt` is a single level per run.
  - Both notebooks output Lu(0−) (Level 4), not Lw at 0+.
  - `AER.Waref` is tied to the simulated wavelength, which flattens the
    aerosol spectrum: a bug for spectral fields.
  - Nothing parses the downward-radiance output we need for Ld.
- Wrote the Setup #1 Q&A: what changed since Phase 0 (two instrument types,
  the air scale, the sky index), the fork findings, and 8 questions with
  defaults.  They cover module names, the OSOAA build, the first cases and
  grid, Lw at 0+, modelling two instruments, the O₃ source and
  housekeeping, including the duplicate "8e" in the Phase 0 doc.
- No package code changed; `pytest -q` not re-run (58 passed at the end of
  Phase 0 task 8e).

### 2026-10-03 -- Build #2: OSOAA built (Opus 5.5)

- JXP accepted all 8 Setup defaults.  The Q8 housekeeping is done: the
  Phase 0 ρw-wiggles task is now 8f (the Sky-index task keeps 8e), and the
  stale lines here are fixed (the test count; the Phase 0a SRF model exists;
  two instruments).
- **Build.**  In the fork, `mkdir obj; cp gen/Makefile_OSOAA.gfortran
  Makefile; OSOAA_ROOT=<fork> make all`.
  - It compiled the 12 sources and linked `exe/OSOAA_MAIN.exe` (787 kB) in
    17 s with gfortran 16.2.0.
  - **No `-std=legacy` and no source fixes were needed**, and the build log
    has no warnings.
  - `git status` in the fork shows only the untracked `Makefile`; `exe/`
    and `obj/` are not reported.  Nothing was committed there.
- **`OSOAA_ROOT`.**  For a session: `export
  OSOAA_ROOT=/Users/xavier/Oceanography/python/RadiativeTransferCode-OSOAA`.
  To make it permanent, add that line to `~/.zshrc` yourself; I did not edit
  shell profiles (Q3).  `hypernet.rt.osoaa.osoaa_root()` (task 3) will fall
  back to that path.
- **`OSOAA.View.Level`, settled from the source.**
  - `src/OSOAA_MAIN.F` (the comment block at lines 1026-1034, the parser at
    3073, the 1-5 range check at 4174) and `src/OSOAA_SOS.F` (2077-2082,
    which set `Z_OUT`) agree: **1 = TOA, 2 = sea bed, 3 = sea surface 0+,
    4 = sea surface 0−, 5 = user altitude/depth (`-OSOAA.View.Z`)**.
  - The fork's Sphinx docs (`input_parameters.rst`,
    `parameter_reference.rst`, `faq.rst`) giving 1/−1 are wrong, and their
    `-1` examples would fail the range check.
  - So **we need 3 (0+) for Ld and 4 (0−) for Lu**.  The source files are
    ISO-8859-1 (French comments), so read them with `LC_ALL=C grep -a`.
- **Smoke test.**  `pytest tests/test_water_Ed.py` in the fork, with
  `OSOAA_ROOT` set, `-p no:cacheprovider` and `PYTHONDONTWRITEBYTECODE=1`
  so that nothing is written into the fork: **33 passed in 346 s**.
- **Timing.**  Every single-wavelength run took 10.2-11.1 s (the default
  setup: 100 m sea, Junge phytoplankton, mono-modal aerosol, SZA 30°); the
  four-SZA test took 22 s for 2 runs.
  - So the task 5 grid is cheap at 5 nm: 125 wavelengths ≈ 21 min per case
    and level.
  - Q4's default (24 cases: 3 SZA × 2 AOT × 4 water) is ~8.5 h for one
    level, or ~17 h if Ld (Level 3) and Lu (Level 4) need separate runs.
  - Runs are single-threaded, so 4-8 in parallel would cut that to a few
    hours.
- No hypernet package code changed for task 2; `pytest -q` (hypernet): 60
  passed.

### 2026-10-03 -- Build #3: `hypernet/rt/osoaa.py` (Opus 5.5)

- New subpackage `hypernet/rt/` (Q2), with `osoaa.py`:
  - `osoaa_root(path=None)`: the argument, then `$OSOAA_ROOT`, then the
    fork's default path (Q3); raises if the exe is missing.  Also
    `exe_path()`.
  - `default_params(wavelength_nm, sza, work_dir, mie_dir=None, phi=90,
    level=3, aot550=0.1, aer_waref_um=0.55, chl=1, csed=0, ys440=0,
    det440=0, wind=5, sea_depth=100, pressure=1013)`.
    - Its aerosol, hydrosol and surface settings are the fork's defaults,
      but **`AER.Waref` is fixed at 0.55 µm** (the Setup bug: the fork tied
      it to the run wavelength, which made AOT spectrally flat).
    - `mie_dir` lets one Mie/surface-matrix database be reused across runs.
    - It always requests `LUM_vsVZA.txt` and both Advanced files.
  - `build_command(params)`, and `run(params, work_dir=None)`, which runs
    one wavelength in a subprocess with `OSOAA_ROOT` in its environment,
    raises `RuntimeError` on failure, and returns the parsed outputs plus
    `seconds`.
  - Parsers:
    - `parse_flux` (`Flux.txt`, every level);
    - `parse_lum_advanced` (Up or Down: level, z, signed VZA, scattering
      angle, I/Q/U, polarisation; the header's 0+/0− level numbers and the
      azimuths for ±VZA);
    - `parse_lum_advanced_down` (checks the direction);
    - `parse_lum_vsvza` (the standard per-level file, with I and REFL);
    - `radiance_at(adv, level, vza)` (VZA interpolation; level `'0+'`/`'0-'`
      allowed).
- **A finding that simplifies task 5.**  The keywords
  `-OSOAA.ResFile.Adv.Up` and `-OSOAA.ResFile.Adv.Down` (`OSOAA_MAIN.F`
  1067-1085, parser 3102-3110) write the upward and downward radiance **at
  every level and every VZA** for the chosen azimuth.  So **one run per
  wavelength gives Ld at 0+ (Down, level 26), Lu at 0− (Up, level 27) and
  Lu at 0+**, and no separate Level-3/Level-4 runs are needed.  That halves
  the task 5 budget: ~21 min per case on the 5 nm grid; 24 cases ≈ 8.5 h
  serial.
- **Units.**  Radiances are I = πL/E_sun and the TOA flux is π, so L/E =
  I/F, and REFL in the vsVZA file = πI/F(level).  At 550 nm (SZA 40°, AOT
  0.1) the fixture gives Ld/Ed(0+) in 0.01-0.1 sr⁻¹ at the HYPSTAR sky view;
  it is a test bound.  ~10 s per run.
- Fixtures in `hypernet/tests/data/osoaa/` come from one real run (550 nm,
  SZA 40°, RAA 90°, AOT(550) 0.1, Chl 1).  They are `Flux.txt`,
  `LUM_vsVZA.txt`, `ListParam.txt` and the two Advanced files, **trimmed to
  levels 0, 25, 26, 27 and 28** (1.26 MB → 60 kB each, headers kept), plus a
  README.
- New `hypernet/tests/test_osoaa.py`, 6 tests:
  - the flux parser (TOA = π cos SZA, direct + diffuse = total, Ed falling
    with depth);
  - the advanced parsers (direction, 0± levels, azimuths, zero downward
    radiance at TOA, a direction check);
  - sky and water radiances at the HYPSTAR geometry, with Lu at the
    refracted angle;
  - vsVZA REFL = πI/F;
  - the aerosol-reference fix;
  - a live 550 nm run (`needs_osoaa`), which reproduces the fixture fluxes
    to 1e-4.
- `pytest -q` (hypernet): see the next log entry.
- **Task 4 head start** (from Phase 0 task 8d): `hypernet/emod.py` already
  has the HITRAN line lists (O₂ and H₂O, 379-1000 nm, cached under
  `$OS_COLOR/hypernet/wiggles/ref/hitran`), `cross_section`, `columns`,
  `optical_depth` and `gas_transmittance(lam, airmass, pwv_mm,
  pressure_hpa, temperature_k)`, with tests in `hypernet/tests/test_emod.py`.
  Task 4 still needs ozone (Serdyuchenko), `build_emod`, the diffuse air
  mass, the 0.01 nm full-range cache and its tests.
- Also from 8d, for the twin experiment: the irradiance SRF is mildly
  flat-topped (super-Gaussian p ≈ 2.2), the radiance SRF Gaussian.  The L
  wavelength scale has a dispersion term (−0.03 to −0.05 nm per 100 nm), so
  case (iv) should include a *linear* wavelength error (±0.1 nm at the ends),
  not only a rigid shift.
- `pytest -q` (hypernet): 71 passed.

### 2026-10-04 -- Build #4: Emod builder (Opus 5.5)

- Completed `hypernet/emod.py`, which 8d started:
  - `build_emod(sza, pwv_mm=15, ozone_du=300, pressure_hpa=1013.25,
    temperature_k=296, ozone_temperature_k=293, diffuse_airmass=1.66,
    lam_min=380, lam_max=1000, step=0.01, cache=True, cache_dir=None)`
    returns `lam` (a regular 0.01 nm **air** grid, 62,001 points), `F0`
    (TSIS-1 HSRS, W m⁻² nm⁻¹, bin-averaged), `T_direct`, `T_diffuse`,
    `T_gas_direct`, `T_o3_direct`, the air masses and the inputs.
  - It is cached as npz in `$OS_COLOR/hypernet/wiggles/ref/emod/`, keyed by
    every input.  A first build takes 3.2 s; a cache hit is instant and
    identical.
  - `gas_transmittance(lam, airmass, pwv_mm, pressure_hpa, temperature_k)`
    now forms exp(−mτ) on the fine HITRAN grid (0.02 cm⁻¹) and
    **bin-averages it** onto `lam` (flux-conserving `_bin_average`).  The
    8d interpolation is kept as `bin_average=False`.
  - `ozone_transmittance(lam, airmass, ozone_du, temperature_k=293)` and
    `ozone_cross_section` read Serdyuchenko et al. (2014).
    - URL: `https://www.iup.uni-bremen.de/gruppen/molspec/downloads/serdyuchenkogorshelev5digits.dat`.
    - SHA-256: `4dfbf021b746512c192df5f0d43c54cee6ea3b4365bb490bcf6ed347f0ce7092`.
    - 12.5 MB, 213-1100 nm at 0.01 nm, vacuum, 11 temperatures 193-293 K.
    - Stored in `$OS_COLOR/hypernet/wiggles/ref/ozone/`, with an npz copy.
  - Also `fetch_ozone`, `kasten_young_airmass` and `f0_on_grid` (HSRS
    bin-averaged in air).
- **Documented choices** (module docstring):
  - The direct air mass is Kasten & Young (1989).
  - The diffuse air mass is fixed at **1.66** (the two-stream diffusivity
    factor: skylight is scattered above most of the H₂O and arrives from
    all directions); case (vi) bounds it.
  - One gas layer at surface p and T.
  - O₃ at 293 K (Q7).
  - O₂-O₂ CIA (477/577/630 nm) is not in HITRAN line lists, so it is not
    modelled (broad, a few per cent).
  - No Rayleigh or aerosol: OSOAA supplies them, smooth, in task 5.
- **Full-range cross-sections** (379-1001 nm, 819,761 wavenumber points)
  are cached.  O₂ is instant; H₂O took 21 s for 163,805 lines.
- **VEIT geometry check** (SZA 36.7°, air mass 1.246):
  - T_direct(760.6 nm) = 0.128 and T_diffuse = 0.065, as 0.01 nm bin means.
  - Band means: O₂-A 759-770 0.60; O₂-B 686-695 0.85; H₂O 925-960 0.48.
  - O₃ T(600) = 0.950, which matches 300 DU × σ × m by hand.
  - F0(550) = 2.0 W m⁻² nm⁻¹.
- Tests in `hypernet/tests/test_emod.py`, now 6:
  - the Kasten-Young values;
  - `_bin_average` conserves flux (analytic sin bins, 1e-5);
  - O₃ transmittance at the Chappuis peak;
  - `build_emod`: 0.01 nm spacing, 380-1000 nm, T in [0, 1], O₂-A
    T_direct(760.6) in the **stated range 0.05-0.25** (air mass 1.25,
    0.01 nm bins), T_diffuse < T_direct there, F0(550) in 1.5-2.2, and a
    cache hit reproducing every array exactly (in a tmp dir).
  - The HITRAN-, HSRS- and O₃-dependent tests skip when the files are
    absent.
- `pytest -q`: 75 passed.

### 2026-10-04 -- Build #5: OSOAA smooth fields (Opus 5.5)

- New `hypernet/wiggles/phase1_osoaa_fields.py` (`--veit`, `--grid
  --workers N`, `--assemble`).
  - It runs one OSOAA call per wavelength on the 5 nm grid 380-1000 nm (125
    wavelengths).
  - Extracted, in units of E_sun:
    - `ed_dir` and `ed_dif` at 0+ (Flux.txt / π);
    - `ld`, the sky radiance at 0+ for the HYPSTAR sky view (40° from
      zenith, RAA 90°);
    - `lu0m`, Lu(0−) at the refracted angle (28.65°);
    - `lw` = lu0m·(1 − ρ_F)/n² (Q5);
    - `lu0p`, Lu(0+) at 40°, for checks.
  - Each case is saved on its own (resumable) under
    `$OS_COLOR/hypernet/wiggles/phase1/osoaa_cases/`, with its own
    Mie/surface database (`osoaa_db/<case>/`, 1.1 GB in all; deletable), so
    that parallel cases never write the same file.
  - Run directories are deleted after parsing.
  - The assembled file is `$OS_COLOR/hypernet/wiggles/phase1/osoaa_fields.npz`
    (25 cases × 125 wavelengths, plus the `raw_*` arrays, an
    `extrapolated` mask and the case metadata).  The figure is
    `hypernet/wiggles/figs/phase1/osoaa_fields.png`.
- **OSOAA wrapper fixes** (`hypernet/rt/osoaa.py`):
  - When `OSOAA.Wa` ≠ `AER.Waref`, OSOAA needs the aerosol index at the
    reference wavelength.  `AER.MMD.MRwaref`/`MIwaref` (1.45, −0.001; no
    dispersion) are now set.  The 550 nm live test never hit this.
  - `default_params(csed>0)` now adds the required `SED.JD.*` Junge
    parameters (slope 4, relative index 1.15 − 0.001i, rate 1).
- **Cases (Q4).**
  - The VEIT-like case: SZA 40°, AOT(550) 0.10, Chl 1.
  - The grid: SZA {30, 50, 70} × AOT(550) {0.05, 0.25} × water {Chl 0.1, 1,
    10 mg m⁻³; turbid = Chl 1 + 5 mg L⁻¹ sediment + YS(440) 0.1 m⁻¹}.
- **Timing.**  Only a case's first wavelength is slow (~10 s, for the
  SZA/wind surface matrices and the first Mie tables); later wavelengths
  take 0.4-0.7 s.  So:
  - the VEIT case took 71 s;
  - the 24-case grid took **6.0 min on 12 workers** (92-218 s per case);
  - every wavelength succeeded (no NaNs).

  The day-long budget estimated in task 2 was far too pessimistic: that
  smoke test built a fresh database for every run.
- **Below 400 nm.**  OSOAA's phytoplankton absorption table
  (`fic/OSOAA_SEA_PHYT_COEFFS.txt`, Bricaud) starts at 400 nm; the
  pure-water table starts at 200 nm.
  - At 380-395 nm OSOAA runs without pigment absorption: ρw is 0.13-0.24 in
    the Chl cases, and through coupling Ed and Ld step by 2-5 %.  The
    second difference of ln Ed is ±0.03-0.06 there, against −0.0005 above
    400 nm.
  - Every quantity below 405 nm is replaced by a **quadratic-in-ln
    extrapolation fitted over 405-445 nm**, which keeps the curvature
    continuous (junction second differences ±0.0003).  The raw values are
    kept, and the fork is not edited.
  - Ca K (393 nm) and Ca H (397 nm) are therefore on extrapolated smooth
    fields, which is fine for multipliers that only need to be smooth.
- **Clear-floor check** (the sky index at 750 nm, observed clear floor
  0.011-0.015 sr⁻¹):
  - The VEIT-like case gives **0.0147**.
  - AOT 0.05 gives 0.0095 (SZA 30 and 50) and 0.0143 (SZA 70).
  - AOT 0.25 gives 0.030-0.039: hazy but still "clear" (< 0.05).
  - The observed floor is bracketed by AOT 0.05-0.10.  The water case does
    not change the sky (to 4 decimals), as it should not.
- **Other checks.**
  - Diffuse fraction at 750 nm: 0.04-0.29.
  - The effective surface reflectance ρ_eff = (Lu(0+) − Lw)/Ld lies within
    0.017-0.036 over all cases and wavelengths (Mobley ~0.028).
  - The turbid case gives ρw(550) ≈ 0.056.
  - ρw in the NIR falls to 10⁻⁶-10⁻⁴ in the Chl cases: OSOAA's Lw is
    elastic only, with no Raman.
- `pytest -q`: 75 passed (no new tests; the npz is a data product).

### 2026-10-04 -- Build #6: scene composer and ρw library (Opus 5.5)

- New `hypernet/twin.py` (Q2):
  - `convolve_to_grid(lam_hr, spec_hr, grid, fwhm)`: a Gaussian SRF whose
    FWHM may be a constant, one value per pixel, or a callable (e.g.
    `SRFModel.fwhm`).  It preserves constants and line equivalent widths.
    Task 7's `observe` will reuse it.
  - `load_fields()` reads `osoaa_fields.npz`.
  - `rhow_library(fields)`: per water type (chl0.1 / chl1 / chl10 /
    turbid), the elastic ρw on 5 nm from the OSOAA case SZA 50°, AOT 0.05,
    plus the control parameters.
  - `compose_scene(emod, fields, case, water=None, controls=('fl','raman'))`
    gives the 0.01 nm scene:
    - Ed = F0(ed_dir·T_direct + ed_dif·T_diffuse);
    - Ld = F0·ld·T_diffuse;
    - Lw = ρ_el·Ed/π + Lw_fl + Lw_raman;
    - Lu = Lw + ρ_eff·Ld (ρ_eff from OSOAA).

    The smooth fields go 5 nm → 0.01 nm by a cubic spline in ln.
  - `rhow_true(scene, srf_L, grid_L)` = πLw/Ed, with both observed with the
    L SRF on the L grid.
- **Controls.**
  - **Fluorescence:** a Gaussian at 683 nm, 25 nm FWHM.  Its peak ρ is
    `fluorescence_amplitude(chl)` = 3e-4·Chl^0.7 (6e-5 / 3e-4 / 1.5e-3 for
    Chl 0.1 / 1 / 10; turbid = Chl 1).  It is applied to Ed smoothed with a
    5 nm Gaussian, so it carries no Fraunhofer structure.
  - **Water Raman** (`raman_radiance`): Ed at the excitation wavelength
    (3357 cm⁻¹ shift), broadened by a 200 cm⁻¹ Gaussian band.  The
    efficiency (`raman_efficiency`) is C·b_R(exc)/(2a_t(exc) + a_t(λ)), with
    b_R ∝ exc^−5.5, a_t = a_w + A·Chl^E + a_ys (OSOAA's pure-water and
    Bricaud tables), and C set so that pure water gives 5e-4 at 550 nm.
- **Two corrections made while checking the figure**
  (`hypernet/wiggles/phase1_scene_check.py` →
  `hypernet/wiggles/figs/phase1/scene_check.png`):
  1. The fluorescence level first used F0 without the gas transmittance, so
     the peak read 19 % high near O₂-B/H₂O.  It now uses the smoothed
     actual Ed.
  2. A constant Raman efficiency gave NIR Raman *larger* than the elastic ρw
     (Δρw ≈ 0.001 at 950 nm).  The absorption-limited model fixes the red,
     and the water-dependent a_t stops the blue being too large in greener
     water.  Final Raman shares of ρw at 550 / 650 nm: 9 % / 11 % (Chl 0.1),
     3 % / 4 % (Chl 1), 1 % / 1 % (Chl 10), ~0 (turbid).
- **A real feature, kept on purpose.**  Fluorescence and Raman are emitted
  in the water, so they do not cross the O₂/H₂O path, while Ed does.  So
  πLw/Ed genuinely peaks in the gas bands (e.g. the fluorescence bump at
  O₂-B, 687 nm).  This is the "true structure at low-transmittance
  wavelengths" of Kevin's doc, which a correction must not erase.
- Tests: `hypernet/tests/test_twin.py`, 6:
  - `convolve_to_grid` keeps constants and equivalent widths, and accepts a
    callable;
  - flat Emod and flat fields give flat ρw (1e-10);
  - the fluorescence control peaks at 683 nm with the stated amplitude and
    is nil outside;
  - the Raman excitation is 464.4 nm for 550, and a line reappears at its
    emission;
  - the Raman efficiency falls in the red and NIR and with Chl;
  - a real VEIT scene (Tier 2): every term finite and ≥ 0, the elastic truth
    = the smooth ρw to 2 %, and the Raman share in the green between 0.5 %
    and 50 %.
- `pytest -q`: 81 passed.

### 2026-10-04 -- Build #7: instrument model (Opus 5.5)

- `hypernet/twin.py`, added:
  - `gaussian_srf(dlam, fwhm)`, a density that integrates to 1.
  - `observe(lam_hr, spec_hr, grid, srf_model, seed=None, noise=None,
    true_grid=None, dfwhm=0)`.  It takes an SRFModel, a callable, an array
    or a constant FWHM.  `true_grid` is where the pixels really sample, for
    a wavelength error; `noise` is a relative 1-σ per pixel, seeded;
    `dfwhm` is a width offset.
  - `load_grids()` and `load_instrument_srfs()`.
  - The `Case` dataclass and `case_table(instrument)`, 16 cases per
    instrument:
    - (i) `i_offset`: equal SRFs (both the L SRF), real grids;
    - (ii) `ii_mismatch`: the E and L SRFs, E on the L grid;
    - (iii) `iii_both`;
    - (iv) rigid E-grid errors of ±0.1 and ±0.3 nm, **plus a linear
      stretch of ±0.1 nm at 390/870 nm** (Phase 0 8d);
    - (v) the correction told FWHM_E or FWHM_L ±0.3 nm;
    - (vi) a wrong Emod: SZA + 10°, water vapour × 2, or low-aerosol fields;
    - (vii) VEIT per-scan noise for a 6-scan mean.
- New `hypernet/wiggles/phase1_instrument.py` writes two committed data
  files.  `setup.py` `package_data` now also ships `data/*.npz`.
  - **`hypernet/data/hypstar_grids.npz`** (75 kB):
    - `grid_E_122304` (1536) and `grid_L_122304` (1538) from the VEIT L1A;
    - `grid_L_122305` from a Release 2 L2B.  **122305 has 1536 L pixels,
      not 1538: pixel counts are instrument-specific.**  It has no E grid
      without its L1A, so the cases use 122304's;
    - `noise_{E,Ld,Lu}_scan`, the flattened per-scan relative scatter (Phase
      0 Q9).  Medians over 400-900 nm: E 0.41 %, Ld 0.43 %, Lu 0.75 %.
  - **`hypernet/data/release2_srf_models.json`**: E and Ld SRFModels per
    instrument × calibration period (keys `HYPSTAR_<id>@<IRR cal date>`,
    plus `HYPSTAR_<id>` for the latest period; Q6).
    - Built from **clear-sky** Release 2 sequences only (sky index < 0.05;
      the fixed-veil Ld widths carry a Ring dependence on sky, per 8d).
    - Per 10 nm window: the median FWHM and dlam over sequences, with error
      1.25 MAD/√n, then `SRFModel.from_lines` with the covariance × χ²_ν.
    - 5-50 sequences per group.  FWHM_E/FWHM_Ld at 450 nm: 122304 (2024-11)
      2.35/2.81; 122305 2.68/2.77; 122302 2.19/2.72; 120242 2.17/2.66;
      121222 (2025-03) 2.70/2.73; 121231 2.83/2.88; 122303 2.67/2.76.
    - χ²_ν is 10-200, because the median errors are small and the windows
      differ by more than a quadratic allows.  Meta notes the proxy E grids
      and the ~0.1 nm E bias of the L2 route.
- **The twin runs every case for 122304 (narrow-E) and 122305 (E ≈ L)**
  (Q6).
- Tests added to `test_twin.py` (5):
  - `gaussian_srf` integrates to 1;
  - **a ~delta line observed at FWHM 3 nm returns 3.00 nm** (±0.01) through
    `srf.fit_line`;
  - a +0.3 nm true-grid error shifts the fitted centroid by −0.3 nm; seeded
    noise is reproducible and of the right size; `dfwhm` adds in
    quadrature;
  - the grids, noise and instrument SRFs load, and keep the narrow-E / E ≈ L
    split (122304 Ld − E > 0.3 nm, 122305 within 0.15 nm);
  - the case table contains (i)-(vii), with equal SRFs in (i), no offset in
    (ii), a stretch of ±0.1 nm at the ends, and so on.
  - Note: noise-free, exactly Gaussian synthetic lines leave `curve_fit`'s
    covariance undefined (zero residuals), so `fit_line` returns NaN.  The
    tests add 1e-4 seeded noise.
- `pytest -q`: 86 passed.

### 2026-10-04 -- Build #8: minimal interpolator (Opus 5.5)

- New `hypernet/edinterp.py` (Q1): `interpolate_ed_to_l(wav_irr,
  irradiance, wav_rad, *, emod=None, srf_irr=None, srf_rad=None,
  method='srf')`, with `METHODS = ('linear', 'ruddick2023', 'srf', 'cubic',
  'sinc')`.
  - `linear` is `np.interp`.
  - `ruddick2023` is Emod_E(λ_L)·interp(E/Emod_E(λ_E)) (eq. 14, the E SRF on
    both sides).
  - `srf` is Emod_L(λ_L)·interp(E/Emod_E(λ_E)) (plan §3, the L SRF in the
    numerator).
  - `cubic` is a scipy CubicSpline.  `sinc` is a Lanczos-windowed sinc
    (a = 8) in fractional-pixel space of the E grid (the grids are not
    uniform).  Both are nulls.
  - `emod` is `(lam, flux)` or a dict with `lam` and `E`/`Ed`/`Emod` (e.g. a
    twin scene).  It is convolved to each grid by
    `twin.convolve_to_grid` (`model_on_grid`).
  - The SRF arguments accept a constant, a per-pixel array on their own
    grid (resampled when `ruddick2023` evaluates the E SRF at the L
    pixels), a callable or an `SRFModel`.
  - Outside `wav_irr` the linear-weight methods hold the end values, as
    `np.interp` does.
- Tests: `hypernet/tests/test_edinterp.py`, 6, on synthetic 0.01 nm spectra
  with ~900 narrow lines and HYPSTAR-like grids with a drifting phase:
  - **with a constant Emod, linear = ruddick2023 = srf = `np.interp` to
    1e-12**;
  - **with srf_rad = srf_irr, `srf` equals a pixel-by-pixel direct
    implementation of eq. 14 to 1e-12** (and `ruddick2023` to 1e-14);
  - **a spectrum equal to Emod_E returns Emod_L** (and α·Emod_E returns
    α·Emod_L) to 1e-12;
  - with E at FWHM 2.3 and L at 2.8, the line residual of truth_L/Ed_L is
    < 1e-6 for `srf`, and below 1 % of linear's;
  - the cubic and sinc nulls reproduce smooth functions to 2e-3;
  - per-pixel SRF arrays give the same results as constants;
  - error paths (no emod, unknown method).
- `pytest -q`: 92 passed.

### 2026-10-04 -- Analysis #9: cases × methods (Opus 5.5)

- New `hypernet/wiggles/phase1_twin.py` (`--workers`, `--scenes`).
  - The grid: 2 instruments (122304 narrow-E, 122305 E ≈ L) × 16 cases
    (`twin.case_table`) × 25 OSOAA scenes × 5 methods (linear, cubic, sinc,
    ruddick2023, srf) × h = 1 and 5 nm: **8,000 metric rows in about 1
    minute on 12 workers** (~5 s per scene).
  - Per scene: the 0.01 nm scene with the fluorescence and Raman controls;
    the truth (`rhow_true`, L SRF on the L grid); Lu and Ld through the L SRF.
  - Per case: E through the case's E SRF at its true grid, reported on the
    nominal grid; for case vii, VEIT per-scan noise/√6 on E, Lu and Ld.  The
    correction gets the case's SRFs (±0.3 nm in v) and an Emod: the scene's
    own high-resolution Ed (a perfect model), or a case-vi variant.
  - ρw = π(Lu − ρ_eff·Ld)/Ed_L, then ρw'' by h-second differences; the metric
    is rms(ρw'' − ρw''_true) within ±5 nm of the plan's ten lines and away
    from them, over 400-900 nm, plus per line.
  - Outputs: `$OS_COLOR/hypernet/wiggles/phase1/twin_metrics.parquet`;
    `hypernet/wiggles/phase1_twin_metrics.csv` (medians over scenes, with
    `reduction` = 1 − rms/rms_linear and `line_over_away`); figures
    `figs/phase1/twin_rho2_lines.png` and `twin_reduction.png`.
- **Results, 122304 (narrow-E), h = 1 nm**, medians over the 25 scenes:
  - The true line structure (controls and real features) is rms ρw''
    2.6e-5.  The linear method's line error is **1.45e-4, about 6× the true
    structure**.
  - **(ii)/(iii) SRF mismatch: srf removes ~100 %** (1.45e-4 → 7e-8, the
    numerical floor).  `ruddick2023`, `cubic` and `sinc` remove **0 %** (or
    6 % worse).  Eq. 14 cannot fix a width mismatch, as plan §3 predicted.
  - **(i) offset only:** the linear error is only **5.6e-6, ~4 % of the
    mismatch error** (matching Phase 0's ~3 %).  `ruddick2023` and `srf`
    remove 99 % of it, and even `cubic` 98 %.
  - **(iv) E wavelength error is the main weakness:** a rigid ±0.1 nm shift
    leaves srf a 62 % reduction (error 5.7e-5); ±0.3 nm only 12 %; a ±0.1 nm
    stretch, 69 %.  So the E wavelength scale (relative to the model)
    must be known to ≪ 0.1 nm, i.e. Phase 2 must fit a shift (and
    stretch) along with the SRF.
  - **(v) SRF told ±0.3 nm wrong:** 6-55 % reduction (E −0.3 worst).
  - **(vi) wrong Emod: robust**, ≥ 99 % (SZA + 10°, PWV × 2, low aerosol).
    The model's line *shapes* matter, not its broad level, as plan §3
    predicted.
  - **(vii) noise:** 45 % (h = 1) and 59 % (h = 5).  The srf line error
    (1.1e-4) sits at the away-from-lines noise level (9.7e-5;
    line/away 1.14).
- **Results, 122305 (E ≈ L):**
  - The linear line error in (iii) is 1.23e-5, **already below the true
    structure** (2.6e-5).
  - srf still removes 99.5 % in (ii)/(iii), but it is fragile: a 0.3 nm SRF
    error makes it **5-7× worse than linear** (v), and a ±0.1 nm wavelength
    error leaves it no better than linear (iv).
  - `ruddick2023` is 2.3× *worse* than linear in (iii): linear
    interpolation's own smoothing partly mimics a wider L SRF, and removing
    it exposes the remaining mismatch.
- **What this means for G1** (task 12 decides):
  - The 80 % target in (ii)-(iii) is met by srf (~100 %) and not by
    ruddick2023 (0 %), so H2 is a requirement on narrow-E instruments.
  - Graceful degradation holds for (vi), but (iv) and (v) show that the
    SRF width (to ≪ 0.3 nm) and the E wavelength scale (to ≪ 0.1 nm) must
    be known.  Phase 0 found per-instrument SRF stability of 0.01-0.06 nm,
    so the table route is viable; wavelength errors call for fitting a
    shift and stretch.
  - On E ≈ L instruments the correction gains little and risks much, so
    Phase 2 might apply srf only where |FWHM_L − FWHM_E| is above a
    threshold.
- The controls (fluorescence, Raman, the real band peaks) are in the
  truth; task 10 quantifies how much each method changes them.
- `pytest -q`: 92 passed.

### 2026-10-04 -- Analysis #10: controls and sensitivity (Opus 5.5)

- Extended `hypernet/wiggles/phase1_twin.py` with `--task 10`
  (`run_scene_t10`, `sweep_cases`, `summarise_controls`,
  `summarise_degradation`, `figures_t10`).  All 25 scenes × 2 instruments in
  **22 s** on 12 workers (17,000 control rows, 14,700 sweep rows).
  - **Controls.**  Each scene is composed three times: all controls, no
    fluorescence, no Raman.  Ed (so every method's Ed_L) is identical across
    the three, so a method's view of a control is ρ_m(all) − ρ_m(without),
    and the truth's likewise.  Fluorescence is judged over 665-705 nm
    (`FL_WINDOW`, it peaks in O₂-B); Raman within ±5 nm of the ten lines.
    Every case (i)-(vii) and all five methods.  In case vii the three
    variants share one noise draw.
  - **Noise floor:** rms[(ρ_m(vii) − ρ_m(iii))''] in the same region, per
    method; it is method-independent to a few %.
  - **Degradation sweeps** (`SWEEP`, case iii with one perturbation; linear,
    ruddick2023, srf): rigid E shift ±0.3 nm (11 points), stretch ±0.2 nm
    (7), correction FWHM_E and FWHM_L error ±0.3 nm (9 each), Emod SZA
    offset ±10° (7), Emod PWV factor 0.5-3 (6).  `main_t10` pre-builds the
    extra Emods (~3 s each, cached under `ref/emod/`) so the workers only
    read them.
  - Outputs: `$OS_COLOR/hypernet/wiggles/phase1/twin_controls.parquet`,
    `twin_degradation.parquet`; committed `hypernet/wiggles/phase1_controls.csv`
    (medians over scenes: `ctl_amp`, `ctl_err`, `ctl_level_err`, `tot_err`,
    `noise_floor`, `ctl_err_over_noise`, `ctl_err_over_amp`,
    `ctl_amp_over_noise`) and `phase1_degradation.csv` (`rms_line`,
    `reduction` vs linear at the same perturbation, `reduction_vs_iii_linear`);
    figures `figs/phase1/twin_controls.png`, `twin_degradation.png`.
- **Controls (the G1 criterion).**
  - **No method changes a control by more than 0.20 of the case (vii)
    noise floor**, in any case (i)-(vii), on either instrument, at h = 1 or
    5 nm (worst: fluorescence, 122304, h = 5, ±0.3 nm shift).  In case (iii)
    srf changes neither control at all (ratio 1e-13); linear, cubic, sinc and
    ruddick2023 change them by 0.02-0.06 of the floor.
  - **The criterion is weak as worded.**  The controls' own ρw'' is itself
    below the single-sequence noise: fluorescence 0.17 and Raman 0.05 of the
    floor at h = 1 nm, 1.3 and 0.4 at h = 5 nm.  A sharper measure is
    the error relative to the control (`ctl_err_over_amp`).  In case (iii),
    122304, h = 1 nm: linear distorts the fluorescence ρw'' by 23 % and the
    Raman ρw'' by 44 %, ruddick2023 by 31 % and 46 %, srf by 0.  The
    distortion is Ed_L's line error multiplying the control (the controls
    *peak* in O₂-B/O₂-A, where that error is largest; `twin_controls.png`
    top row).  ρ-level errors are < 1 % of the control for every method.
    srf's worst case over (iv)-(vi) is 48 % (fluorescence, ±0.3 nm shift),
    no worse than linear's 51 %.  So the corrections never *remove* a real
    feature; srf restores it, and the other methods distort it.
  - Suggest G1 reads the control criterion as "control error ≤ linear's, and
    ≪ the control itself in case (iii)", or judges it on a multi-sequence
    floor (√N lower) for Phase 3 composites.
- **Degradation curves (122304, h = 1 nm; reduction vs linear):**
  - E wavelength, rigid shift: ±0.02 nm → 92 %, ±0.05 → 80 %, ±0.1 → 62 %,
    ±0.2 → 32 %, ±0.3 → 12 %.  Stretch (at the ends): ±0.05 → 84 %, ±0.1 →
    69 %, ±0.2 → 44 %.  **The 80 % target needs the E-relative-to-model
    wavelength scale to ≤ 0.05 nm.**
  - Correction FWHM error: FWHM_E ±0.05 → 86-87 %, ±0.1 → 72-74 %, ±0.2 →
    40-50 %, ±0.3 → 6-29 %; FWHM_L is less sensitive (±0.1 → 81-83 %, ±0.3 →
    38-55 %).  The response is linear in |ΔFWHM| near 0, so **80 % needs
    FWHM_E to ≲ 0.07 nm and FWHM_L to ≲ 0.1 nm**: within Phase 0's 0.01-0.06 nm
    SRF stability, but only with per-instrument SRFs.
  - Emod: SZA ±10° and PWV × 0.5-3 all keep ≥ 98 % (errors ≤ 3e-6, against
    1.45e-4).  The Emod's atmosphere is not a concern.
  - ruddick2023 sits on linear (0 % or slightly worse) along every curve.
- **122305 (E ≈ L):** srf beats linear only for |shift| ≲ 0.05 nm
  and |ΔFWHM| ≲ 0.05 nm.  Beyond that it is no better than linear (shift)
  or several times worse (FWHM: −6.6× at −0.3 nm).  ruddick2023 is 2.3×
  worse than linear at the nominal point, as in task 9.  This supports
  applying srf only above a ΔFWHM threshold (task 11 will give the number).
- `pytest -q`: 92 passed (no new tests; the task 10 code is a script).

### 2026-10-04 -- Analysis #11: prediction for Phase 3 (Opus 5.5)

- New `hypernet/wiggles/phase1_prediction.py` (reuses `phase1_twin`'s masks,
  second differences and constants).  Case (iii) on the real 122304 E/L
  grids with the 122304 Ld SRF model as L and a synthetic E SRF,
  FWHM_E(λ) = FWHM_L(λ) − ΔFWHM, ΔFWHM ∈ {0, 0.05, 0.1, 0.2, …, 1.0} nm; 25
  scenes; linear, ruddick2023, srf; h = 1 and 5 nm.  54,450 rows in 11 s.
  - Three scenarios: `ideal` (noise-free, correction told the truth),
    `noise` (case vii noise), `realistic` (noise + a 0.05 nm E shift + the
    correction told FWHM_E + 0.05 nm, i.e. the edges of Phase 0's
    wavelength-scale and SRF-stability findings).
  - Two metrics:
    - the twin metric (`reduction`: of rms(ρw'' − ρw''_true) near the
      lines, vs linear);
    - **an observable one for G3**, because Phase 3 has no truth:
      excess = √(rms²_line(ρw'') − rms²_away(ρw'')), with `reduction_obs`
      vs linear.  The truth has an excess of its own (real line filling,
      gas-band peaks), so `reduction_obs_max` = 1 − excess_true /
      excess_linear is the most a perfect correction can show.  The
      away-from-lines rms carries the noise, so the quadrature removes it
      to first order.
  - Line depth = 1 − min(Ed_L, ±1 nm)/max(Ed_L, ±5 nm) of the true Ed
    through the L SRF; per line and in depth bins (rows `line =
    'depth[a,b)'`).
  - Outputs: `$OS_COLOR/hypernet/wiggles/phase1/twin_prediction.parquet`;
    committed `hypernet/wiggles/phase1_prediction.csv` (3,168 rows: medians
    over scenes per scenario × ΔFWHM × h × method × line or depth bin);
    `figs/phase1/twin_prediction.png`.
- **Results, ten lines together, h = 1 nm** (srf; ruddick2023 is −4 to −8 %
  for every ΔFWHM > 0 and every scenario, i.e. never helps):

  | ΔFWHM (nm) | 0 | 0.1 | 0.2 | 0.3 | 0.5 | 0.7 | 1.0 |
  |---|---|---|---|---|---|---|---|
  | twin, ideal | 99 | 100 | 100 | 100 | 100 | 100 | 100 |
  | twin, noise | 0 | 2 | 14 | 27 | 51 | 66 | 80 |
  | **observable, noise** | −1 | 11 | 31 | 52 | **74** | 84 | 90 |
  | observable, realistic | −4 | 3 | 21 | 39 | **62** | 73 | 80 |
  | observable max (perfect), noise | 69 | 71 | 80 | 87 | 92 | 94 | 96 |

  - Linear's error grows ~linearly with ΔFWHM above 0.1 nm (rms 2.2e-5 at
    0.1, 1.6e-4 at 0.5, 4.2e-4 at 1.0); the true excess is 2.0e-5.  At
    ΔFWHM ≲ 0.1 nm the wiggle excess is below the true line structure and
    there is little to remove.
  - **For a narrow-E instrument (ΔFWHM ≈ 0.5 nm) the G3 prediction for a
    single sequence is a 60-75 % reduction of the observable excess**
    (realistic-ideal knowledge), with ~90 % the ceiling.  For E ≈ L
    instruments (ΔFWHM < 0.1 nm) it is ≈ 0-10 %: G3 should not expect a
    reduction there.
  - **h = 5 nm is a poor G3 metric:** the true excess dominates (Ca H/K
    wings, broad gas bands), so even a perfect correction shows only 28 %
    at ΔFWHM 0.5 (srf reaches that ceiling).
- **Per line (noise scenario, ΔFWHM 0.5):** Ca K/H 76-78 %, G band 65 %,
  H β 43 %, Mg b/Na D 33 %, O₂-B 36 %, O₂-A 24 %, H α 13 %, H₂O 1 %.  With
  noise the reduction follows the line's SNR (depth × signal / noise), not
  depth alone.  The deep O₂-A (0.69) and H₂O (0.43) bands sit where Lu
  noise is largest (dark water in the NIR), so the depth bins are
  confounded with wavelength; the figure's panel (c) is therefore per line.
  Noise-free, every line and depth bin is ≥ 93 % for ΔFWHM ≥ 0.2.
- **Implications for G3 (task 12):**
  - Judge on the observable excess at h = 1 nm, against the noise row for
    the instrument's measured ΔFWHM, per sequence, with
    `reduction_obs_max` as the ceiling.
  - Expect the clearest signal at Ca H/K and the G band; O₂-A, H α and H₂O
    are noise-limited for one sequence.
  - Averaging N sequences (noise /√N) moves the result toward the ideal
    row.
- `pytest -q`: 92 passed.

### 2026-10-04 -- Gate #12: G1 (Opus 5.5)

- Wrote `claude_prompts/wiggles/gate_G1.md` (the proposed location; Q&A
  "Gate #12" Q1 asks JXP to confirm).  Verdict: **pass for narrow-E
  instruments** with two knowledge requirements (E wavelength ≤ 0.05 nm,
  SRF widths ≲ 0.07-0.1 nm); **H2 is required** (`ruddick2023` 0 % in (ii),
  −6 % in (iii)); Phase 2 default `srf` above ΔFWHM 0.15 nm, `linear`
  below.
  - Criterion 1: pass (~100 %).
  - Criterion 2: pass on narrow-E, conditional on the tolerances; fail on
    E ≈ L.
  - Criterion 3: pass, but weak as worded (rewording proposed).
  - Criterion 4: `ruddick2023` alone does not pass.
- Seven open concerns: the E wavelength scale (binding), E SRFs from L1A,
  SRF shape, physics outside the twin (Ring, the Lu NIR excess), E ≈ L
  instruments and their unexplained floor, the control criterion, and
  single-sequence noise.
- New finding pulled from the task 9 table: on 122305 the twin's rms
  relative ρw error near the lines is 0.08 % (offset only) and 0.46 %
  (offset + residual mismatch).  Phase 0 measured a 0.5-2 % floor, so H1
  does not explain it.
- Proposed §4 edits to `docs/wiggles_planning.md` in the Q&A (Gate #12 Q2):
  the G1 control wording, the G1 outcome, two Phase 2 bullets (the default
  method and threshold, and the E shift/stretch fit), and the G3 yardstick.
  **Not applied** pending JXP's approval.
- No code changes; `pytest -q` not rerun (92 passed at task 11).

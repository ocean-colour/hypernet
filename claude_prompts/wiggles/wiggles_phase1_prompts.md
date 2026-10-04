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

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
- Run `pytest -q` after each step where relevant.  21 tests today; the archive-
  dependent ones skip themselves when `$OS_COLOR` is not mounted.  OSOAA-
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
  `/opt/homebrew/bin/gfortran`.  No exe yet.
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
- Tasks 7-9 need the Phase 0a SRF model; until it exists use the §2.3 table.
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

5. **OSOAA smooth fields.**  `wavecal/phase1_osoaa_fields.py`: on a 5 nm grid
   380-1000 nm, for SZA 30/50/70°, two aerosol loads and a small
   chlorophyll/sediment grid, run OSOAA and collect direct and diffuse Ed at
   0+, Ld at the HYPSTAR sky geometry (VZA 40°, relative azimuth from the
   VEIT sample) and Lw at 0+.  Run one SZA/aerosol/water case first, time it,
   then the rest in the background.  Save
   `$OS_COLOR/hypernet/wiggles/phase1/osoaa_fields.npz` and a figure
   `wavecal/figs/phase1/osoaa_fields.png`.  Log your work.

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

9. **Cases × methods.**  `wavecal/phase1_twin.py`: for every case, SZA,
   aerosol and water type, and every method, compute ρw on the L grid, then
   ρw'' with h = 1 nm and h = 5 nm second differences as in the paper.
   Metric: rms of ρw'' − ρw''_true within ±5 nm of the ten lines and away
   from them.  Write `$OS_COLOR/hypernet/wiggles/phase1/twin_metrics.parquet`,
   a committed `wavecal/phase1_twin_metrics.csv`, and figures under
   `wavecal/figs/phase1/` (ρw'' at O2-A, Hα and Ca H/K per method; the
   reduction per case).  Log your work.

10. **Controls and sensitivity.**  Extend the script: the change in the
    fluorescence and Raman controls' second derivative per method against the
    case (vii) noise floor; degradation curves for (iv), (v) and (vi).
    Table `wavecal/phase1_controls.csv` and a figure.  Log your work.

11. **Prediction for Phase 3.**  From case (iii), tabulate the expected
    reduction of the line-region ρw'' excess for `ruddick2023` and `srf` as a
    function of ΔFWHM = FWHM_L − FWHM_E (0 to 1 nm) and line depth -- the
    number Gate G3 will be judged against.  `wavecal/phase1_prediction.csv`
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

## Logs

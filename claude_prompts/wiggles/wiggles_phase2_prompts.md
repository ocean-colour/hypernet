# WATERHYPERNET -- Wiggles Phase 2: Code

## Goal

Turn the Phase 0 SRF model and the Phase 1 interpolator into the deliverable
code in `hypernet/`: one function `interpolate_ed_to_l` with
`method = {"linear", "ruddick2023", "srf"}`, Fraunhofer-line SRF
self-calibration, a cached Emod builder, `punpy` uncertainty propagation into
`u_rel_random` / `u_rel_systematic` in `obsarray` conventions, a per-instrument
SRF table, and a drop-in measurement-function snippet in the style of
`hypernets_processor`'s `interpolate_wav_linear.py`.  We will not run
`hypernets_processor` end to end.

**Gate:** the plan sets no numbered gate for Phase 2.  The exit criterion is
plan §7 items 1, 2 and 4 delivered with tests, in particular: the two special
cases reproduce `np.interp` and eq. (14) to machine precision on synthetic
data, and the G0 and G1 verdicts are applied (table vs self-calibration,
wavelength recalibration or not, default method).

**G1 has closed (2026-10-04, `claude_prompts/wiggles/gate_G1.md`):**

- Default `srf` where the measured ΔFWHM = FWHM_L − FWHM_E ≥ 0.15 nm;
  `linear` otherwise.
- `srf` needs the E wavelength scale to ≤ 0.05 nm relative to Emod (a shift +
  stretch fit; task 4b) and FWHM_E/FWHM_L to ≲ 0.07/0.1 nm.
- A non-Gaussian SRF twin check is added (task 9b).

**G0 has not closed** (Phase 0 task 14 waits on Kevin's L1A delivery); its
preliminary findings are in Context.

## Conventions

- `ocean14`.  Run via `conda run -n ocean14 python ...`; `conda activate` fails
  non-interactively.
- **Model.** Use Opus 5.5 (`claude-opus-5-5`) for this work, including any
  subagents (pass `model: opus`).
- **JXP runs git.**  Claude does not run any state-changing git command.
- Run scripts **from the repository root** so `hypernet` imports resolve.
- Run `pytest -q` after each step where relevant.  92 tests on 2026-10-05; the
  archive-dependent ones skip themselves when `$OS_COLOR` is not mounted, the
  OSOAA ones when the exe is missing.
- **Tier 2** -- the VEIT check needs `$OS_COLOR`.  Do **not** unset it.
- Any calculation goes into a script on disk, not into the chat.
- Data intermediates go **outside** the repo under
  `$OS_COLOR/hypernet/wiggles/ref/` (Emod cache) and
  `$OS_COLOR/hypernet/wiggles/phase2/`.  Only figures and small tables are
  committed.
- Every `.md` in `docs/` is published.  Nothing unpublished goes there.
- Ask questions in the Q&A section below; log completed work under `## Logs`.

## Context

### Where things are

- **Plan:** `docs/wiggles_planning.md` §3 (the equations), §4 "Phase 2", §7
  (deliverables and uncertainty), §9 (Gaussian shape, Emod line shapes,
  processor drift).
- **Gate reports:** `claude_prompts/wiggles/gate_G1.md` (done; read it
  before task 2).  `gate_G0.md` does not exist yet (Phase 0 task 14).  Until
  it does, use the preliminary G0 findings from the Phase 0 Logs (Build #8c,
  #8d, #8f):
  - **G0(a):** E ≠ L on the narrow-E instruments (122302, 122304, 120242:
    Ld − E ≈ 0.5 nm); E ≈ L within 0.07 nm on 121222, 121231, 122303 and
    122305.  H2 stays.
  - **G0(b), preliminary pass:** SRFs are stable to 0.01-0.06 nm over SZA,
    season and recalibration (L channel from Release 2).
  - **G0(c):** air wavelength scale, mean offset +0.007 to +0.047 nm, so no
    relative (L − E) recalibration is needed.  But there is a dispersion trend
    in L (−0.03 to −0.05 nm per 100 nm) that exceeds 0.1 nm at the ends.
    Together with G1's 0.05 nm tolerance, this motivates task 4b.
- **Phase 0 code:** `hypernet/srf.py`:
  - the line fits: `LINES`, `fit_line`, `fit_lines`;
  - the FWHM model: `fit_fwhm_model`, `SRFModel` (JSON);
  - **the HSRS template fit**, the standard SRF measurement since Phase 0
    task 3b: `fit_srf_template` with `via_grid`, `shape` in `SHAPES`
    (gauss/pvoigt/supergauss) and `absorber`, and `fit_template_windows`.
    Gaussian line fits measure SRF ⊗ the line's intrinsic width, so they are
    diagnostics.

  Readers are in `hypernet/whn_l1a.py` (incl. `load_l2b`, `sky_index`),
  `hypernet/whn_srf.py` and `hypernet/refspec.py` (HSRS).
- **Phase 0 SRF products:**
  - `hypernet/data/veit_srf_model.json` (template; the empirical line-width
    version is `veit_srf_model_empirical.json`);
  - `hypernet/data/release2_srf_models.json`: E and Ld per instrument × cal
    period, clear-sky median template fits over the 224 Release 2 sequences
    (E from L2 reads ~0.1 nm wide; E grids are a 122304 proxy except
    122304@2024-11);
  - `$OS_COLOR/hypernet/wiggles/phase0/` (`srf_models.parquet` = VEIT only
    until the delivery; `release2_template.parquet`, `template_fits.parquet`).
- **Phase 1 code:**
  - the minimal `interpolate_ed_to_l` in the flat module
    **`hypernet/edinterp.py`** (Phase 1 Q&A Q1), `METHODS = (linear,
    ruddick2023, srf, cubic, sinc)`, where cubic/sinc are twin nulls;
  - `hypernet/emod.py` (`build_emod`, `ref/emod/` cache);
  - `hypernet/twin.py` (`observe`, `Case`, `case_table`,
    `load_instrument_srfs`, `load_grids`) for synthetic test data;
  - `hypernet/rt/osoaa.py`;
  - the twin scripts `hypernet/wiggles/phase1_twin.py` (`--task 9|10`) and
    `phase1_prediction.py`;
  - `hypernet/tests/test_edinterp.py`.
- **Phase 1 numbers** (in `hypernet/wiggles/`): `phase1_twin_metrics.csv`,
  `phase1_controls.csv`, `phase1_degradation.csv` (the tolerances), and
  `phase1_prediction.csv` (G3's yardstick).
- **hypernets_processor:** public repo `github.com/HYPERNETS/hypernets_processor`,
  file `hypernets_processor/interpolation/measurement_functions/interpolate_wav_linear.py`.
  Kevin's idea doc links it at commit
  `9a12819a3ffb224de9e5f761b1217d130fa69f83` (the version he points to; the
  public repo is not current with Release 2).  Its interface has not yet been
  re-read; Task 2 fetches it.
- **CoMet toolkit:** `punpy` 1.1.0, `comet_maths` 1.0.10, `obsarray` 1.0.3 in
  `ocean14` and in `requirements.txt`.  De Vis et al. (2024) describes the
  processor's use of them.
- **Data:** VEIT sample at `$OS_COLOR/WATERHYPERNET/Wavelengths/`; the full
  request (if delivered by now) alongside it, indexed by Phase 0b at
  `$OS_COLOR/hypernet/wiggles/phase0/request_index.parquet`.
- **Paper:** `context/papers/ruddick2023.pdf`.

### Dependencies

- Depends on G0 (SRF table vs self-calibration; wavelength recalibration) and
  G1 (default method).  If the full data have not arrived, G0 is provisional
  on the VEIT sample: build both the table loader and the self-calibration
  path, and leave the default to be set when G0 closes.
- G1 is closed; G0 is provisional (above).  Tasks 4-6 build both the table
  and the self-calibration paths so that G0 can set the default later.
- Module naming: the phase scripts live in the subpackage `hypernet/wiggles/`
  (run as `python -m hypernet.wiggles.<script>`), and the reusable code in
  flat modules.  `interpolate_ed_to_l` already lives in `hypernet/edinterp.py`
  (Phase 1 Q1), so the default is to harden it there and put the
  self-calibration, wavelength fit and SRF table next to it.  The plan's
  `hypernet/wavecal/` name is superseded unless JXP prefers a subpackage
  (Setup Q&A).

## Prompts

### Setup

1. Read the plan sections, gate reports and Phase 0/1 modules listed in
   Context, and this doc.  Put your questions in Q&A:
   - module layout (default: flat modules next to `hypernet/edinterp.py`);
   - the snippet's location (proposal `hypernet/snippets/`);
   - how to proceed while G0 is provisional;
   - how the method is chosen (where the ΔFWHM ≥ 0.15 nm rule lives, and
     ΔFWHM over which band);
   - the granularity of the E wavelength fit (per sequence or per
     calibration period);
   - whether to target Kevin's pinned commit or the current default branch
     of `hypernets_processor`;
   - the obsarray variable names to target.

   Do not write code yet.  Log your work.

### Build

2. **Read the processor's interface.**  Fetch
   `interpolate_wav_linear.py` from `github.com/HYPERNETS/hypernets_processor`
   at commit `9a12819a3ffb224de9e5f761b1217d130fa69f83` (and diff against the
   default branch if JXP asks) and the module that registers
   measurement functions; save a copy under `context/hypernets_processor/`
   with the commit hash in a `SOURCE.txt`.  Record in the Logs the class name,
   `function` signature, argument names, return convention, and how
   uncertainties and `obsarray` attributes flow through it.  Log your work.

3. **Harden `interpolate_ed_to_l`.**  Final signature per plan §4:
   `interpolate_ed_to_l(wav_irr, irradiance, wav_rad, *, emod, srf_irr,
   srf_rad=None, method="srf", u_rel_random=None, u_rel_systematic=None)`.
   Add input validation, NaN and out-of-range L pixels handled explicitly,
   2-D `irradiance` (wavelength, scan), a precomputed SRF-convolution matrix
   on the 0.01 nm grid cached per `SRFModel`, and full docstrings with the
   equations.
   - Keep `cubic`/`sinc` only as test nulls (or move them to the twin), as
     agreed in Q&A.
   - Add the G1 selection rule, e.g. `choose_method(srf_irr, srf_rad,
     threshold=0.15)` → `"srf"` if the median FWHM_L − FWHM_E over the
     agreed band is ≥ threshold, else `"linear"`.
   - Extend `test_edinterp.py`: the two machine-precision special cases
     remain, plus shape, NaN, edge, 2-D and threshold tests, and a timing test
     that a 1536 → 1538 call runs under a stated budget.

   Log your work.

4. **SRF self-calibration.**  `calibrate_srf(wav, spectrum, channel="E", ...)`
   → `SRFModel` (centroid offset and FWHM(λ) quadratic with covariance).
   - Build it on the HSRS template fit (`fit_template_windows`, then
     `SRFModel.from_lines`), the Phase 0 standard; `fit_lines` stays a
     diagnostic.
   - Add a quality flag for low blue signal (the low-sun risk in plan §9).
   - Test: an `observe()`d synthetic spectrum with a known FWHM(λ) is
     recovered within its reported uncertainty; a noisy one returns the
     flag.

   Log your work.

4b. **E wavelength-scale fit (G1).**  `fit_wavelength_scale(wav, spectrum,
   srf, emod=None)` → shift and linear stretch (and covariance) of a
   spectrum's wavelength scale relative to the HSRS-based model, fitted
   across the template windows.
   - Apply it as a wavelength correction to E before interpolation, at the
     granularity agreed in Q&A.
   - Requirement: ≤ 0.05 nm (G1); report the achieved uncertainty.
   - Tests: a synthetic E with a known shift/stretch (twin case iv) is
     recovered within its uncertainty and to < 0.02 nm at VEIT noise.
   - Run it on the VEIT sample and the Release 2 E spectra.  Compare with
     Phase 0 8d's dispersion trend, in a small table
     `hypernet/wiggles/phase2_wavescale.csv`.

   Log your work.

5. **Emod cache.**  `emod_for(sza, pwv_mm, ozone_du, pressure_hpa)` returning
   the nearest cached entry on a stated grid (SZA bins of 5°, a few water
   vapour columns), populated on demand under `ref/emod/`, plus cached
   E-SRF and L-SRF convolutions per `SRFModel` id.  Test: two calls in the
   same bin return the same array; a cache miss builds and stores.  Log your
   work.

6. **Per-instrument SRF table.**  A committed `hypernet/data/srf_table.csv`
   (instrument, channel, cal period, quadratic coefficients, centroid offset,
   uncertainties, n_sequences, ΔFWHM and the method the G1 rule selects)
   written by `hypernet/wiggles/phase2_srf_table.py`.
   - Source: Phase 0b's `srf_models.parquet` once the delivery has been
     fitted.  Until then use `hypernet/data/release2_srf_models.json`
     (Release 2 template fits, clear sky) plus the VEIT model, flagged as
     provisional.
   - Provide `load_srf_table()` and `srf_for(instrument, channel, date)`.  If G0(b) failed, the default path
   is `calibrate_srf` per sequence and the table is the fallback.  If G0(c)
   fired, apply the centroid offset as a wavelength recalibration before
   interpolation.  Tests on the CSV.  Log your work.

7. **Uncertainty with punpy.**  Propagate with `punpy.MCPropagation`: the
   SRF parameters (centroid and FWHM, with the fit covariance), the Emod
   inputs (F₀ scale, water vapour, ozone, SZA) and the wavelength grids,
   together with the input `u_rel_random` / `u_rel_systematic`, into
   `u_rel_random` and `u_rel_systematic` on Ed_L following the `obsarray`
   conventions found in task 2.  Tests: seeded reproducibility; for
   `method="linear"` with zero SRF/Emod uncertainty the output equals the
   analytic linear propagation of the input uncertainties; systematic and
   random components separate as expected for a common scale error.  Log
   your work.

8. **Drop-in snippet.**  In the location agreed in Q&A: a measurement-
   function class in the style of `interpolate_wav_linear.py` -- same
   `function` and argument-name interface -- taking the SRF table and Emod
   cache as inputs and calling `interpolate_ed_to_l`.  A test imports it and
   runs it on the Phase 1 synthetic scene.  A short module docstring says
   what RBINS must wire up; no end-to-end processor run.  Log your work.

9. **VEIT consistency check.**  `hypernet/wiggles/phase2_veit_check.py`: run the
   three methods on the VEIT L1A_IRR with the sample's SRF models; confirm
   `linear` reproduces the L1C `irradiance` to the 2 × 10⁻⁵ of plan §2.3, and
   show Ld/Ed_L for the three methods with their `u_rel_*` bands.  Figure
   `hypernet/wiggles/figs/phase2/veit_three_methods.png`; table
   `hypernet/wiggles/phase2_veit_check.csv`.  Log your work.

9b. **Non-Gaussian SRF check (G1 concern 3).**
    `hypernet/wiggles/phase2_shape_check.py`: rerun twin case (iii) (and the
    Phase 1 prediction at a few ΔFWHM) with the truth observed through the
    Phase 0 pvoigt/supergauss SRF shapes (`srf.SHAPES`, fitted parameters
    from the VEIT template fits) and the correction still Gaussian.
    - Report how much of the reduction survives, and whether the residual
      correlates with Ed''/Ed (G3's trigger).
    - Table `hypernet/wiggles/phase2_shape_check.csv` and a figure under
      `figs/phase2/`.

    Log your work.

### Gate

10. **Phase 2 exit check.**  Against plan §7 items 1, 2 and 4 and the G0/G1
    verdicts: list each deliverable, its module, its tests and the
    `pytest -q` count; confirm the machine-precision special cases; state the
    default `method` (and the ΔFWHM rule), the SRF source and wavelength-
    scale treatment the snippet ships with, the task 4b and 9b results, and
    what remains provisional while G0 is open.  Write a short report (proposed
    location `claude_prompts/wiggles/phase2_exit.md`; JXP to confirm -- not
    under `docs/`).  Propose the corresponding edits to
    `docs/wiggles_planning.md` §4 and §7 in Q&A and apply them only if JXP
    approves.  Log your work.

## Q&A

### Setup #1 -- 2026-10-07 (Opus 5.5)

**What I read:** plan §3, §4 (Phase 2, and the G1 outcome now in Phase 1),
§7 and §9; `gate_G1.md`; this doc; `hypernet/edinterp.py`, `srf.py`
(`SRFModel`), `emod.py` and `twin.py`; and the VEIT L1C file's irradiance and
uncertainty variables.  Nothing was fetched or run except two read-only
checks: the L1C variable list, and the time of one `srf` call.

**Findings that bear on the questions**

- **The L1C has three uncertainty components, not two.**
  `u_rel_random_irradiance` (err_corr random × random),
  `u_rel_systematic_indep_irradiance` and
  `u_rel_systematic_corr_rad_irr_irradiance` (both `err_corr_matrix` along
  wavelength, systematic along scan; the matrices are named in
  `err_corr_1_params`).  Units are **%**, dims (wavelength, scan); the
  radiances follow the same pattern.  The plan's `u_rel_random` /
  `u_rel_systematic` signature does not match this.  Phase 0 also found
  every L1C `u_rel_random` value to be 0, a placeholder.
- **`edinterp.py` imports `convolve_to_grid` from `twin.py`**, so the
  deliverable depends on the twin-experiment module.
- **Speed is not an issue.**  One 1536 → 1538 `srf` call with an SRFModel
  takes 0.02 s today.  A cached convolution matrix (task 3) is for punpy's
  Monte Carlo (hundreds of draws), not for a single call.
- **`SRFModel` already carries what punpy needs**: `coeffs` with `cov`
  (FWHM quadratic), and `offset` with `offset_err`.
- **G1 tested only E wavelength errors.**  `srf` also evaluates Emod at the
  L pixels, so an L wavelength error misplaces the model's lines relative to
  the measured L.  The twin has not measured this.

**Questions for JXP** (each has a default; say "defaults" to accept them all)

1. **Module layout.**  *Default: flat modules.*
   - *Harden `interpolate_ed_to_l` in `hypernet/edinterp.py`; add
     `choose_method` there.*
   - *Move `convolve_to_grid` (and the cached convolution matrix) into
     `hypernet/srf.py`, so `edinterp` no longer imports `twin`.  `twin`
     re-exports it, so the Phase 1 scripts are unchanged.*
   - *New `hypernet/srfcal.py`: `calibrate_srf` (task 4) and
     `fit_wavelength_scale` (task 4b).*
   - *New `hypernet/srftable.py`: `load_srf_table`, `srf_for` (task 6).*
   - *`emod_for` in `hypernet/emod.py` (task 5).*
   - *The punpy propagation in `edinterp.py`, behind the `u_rel_*`
     arguments (task 7).*

   *The plan's `hypernet/wavecal/` name is dropped.*

   >A. Use the default
2. **The nulls.**  *Default: public `METHODS = ("linear", "ruddick2023",
   "srf")`; `cubic` and `sinc` move to a private `NULL_METHODS` in the same
   module, still accepted by `interpolate_ed_to_l` so that
   `phase1_twin.py` reproduces its tables.*
   >A. Use the default
3. **The snippet.**  *Default: `hypernet/snippets/interpolate_wav_srf.py`,
   a class named after the processor's linear one (e.g.
   `InterpolateWavSRF`; the exact name follows task 2).*
   - *`hypernets_processor` is not installed, so the test uses a minimal
     stand-in for its base class.*
   - *`hypernet/snippets/` is not imported by `hypernet/__init__`.*
   >A. Use the default
4. **The processor version.**  *Default: target Kevin's pinned commit
   `9a12819a…`.  Diff it against the current default branch in task 2 for
   information only, and record any interface change in the Logs.*
   >A. Use the default
5. **Uncertainty components.**  *Default: match the L1C, not the plan's two
   names.*
   - *`interpolate_ed_to_l` takes and returns the three components,
     `u_rel_random`, `u_rel_systematic_indep` and
     `u_rel_systematic_corr_rad_irr`, in % on (wavelength[, scan]), with
     their err_corr matrices.*
   - *The input components are interpolated (linear weights) and
     propagated.*
   - *The new terms (SRF parameters, Emod inputs, the wavelength fit) are
     common to every scan of an instrument, so they are added to
     `u_rel_systematic_indep`, with a punpy-derived wavelength error-
     correlation matrix.*
   - *Plan §4/§7's signature gets updated at the Phase 2 exit (task 10).*
   >A. Use the default
6. **While G0 is open.**  *Default:*
   - *the snippet ships with the SRF table as its source (G0(b) passes
     preliminarily), and `calibrate_srf` per sequence is an option;*
   - *the table's Release 2 rows are flagged `provisional`;*
   - *no L − E centroid recalibration (G0(c) relative passes); E's absolute
     scale is handled by task 4b;*
   - *revisit all three when `gate_G0.md` exists.*
   >A. Use the default
7. **The method rule.**  *Default:*
   - *ΔFWHM = the median of FWHM_Ld(λ) − FWHM_E(λ) over 400-680 nm, where
     both template models are constrained (the VEIT model is unconstrained
     at 680-850 nm).*
   - *It is computed per instrument × calibration period and stored in the
     task 6 table, with the selected method.*
   - *`choose_method(srf_irr, srf_rad, threshold=0.15)` does the same for
     any pair of models; the snippet reads the table and allows an
     override.*
   >A. Use the default
8. **The wavelength-scale fit (task 4b).**  *Default:*
   - *Fit per sequence, on the scan-mean E, against the HSRS-based Emod
     through the E SRF, for the sequence's SZA.*
   - *Fall back to the calibration-period median when the fit is flagged
     (low blue signal) or its uncertainty exceeds 0.03 nm.*
   - *Fit L (Ld) the same way.*
   - *Add an **L** wavelength-error case to the task 9b twin check, since
     G1 never measured it.*
   >A. Use the default
9. **The Emod cache grid (task 5).**  *Default:*
   - *SZA bins of 5° (each entry builds in ~3.5 s);*
   - *water vapour snapped to {7.5, 15, 30} mm;*
   - *ozone 300 DU and pressure 1013 hPa fixed.*

   *G1 found the correction insensitive to SZA ± 10° and water vapour
   × 0.5-3, so this grid is generous.*
   >A. Use the default
10. **Edges.**  *Default: L pixels whose SRF reaches beyond the Emod range
    (380-1000 nm), or beyond the E grid, get `linear` and a flag; NaN
    inputs propagate as NaN at the affected pixels only.*
   >A. Use the default
11. **The exit report.**  *Default:
    `claude_prompts/wiggles/phase2_exit.md`, as proposed in task 10.*
   >A. Use the default

## Logs

### 2026-10-05 -- Doc refreshed after G1 (Opus 5.5)

- Brought the doc up to date before Setup #1, from the Phase 1 session:
  - **Goal:** the G1 outcome, and G0 still open.
  - **Context:**
    - the preliminary G0(a)/(b)/(c) findings;
    - the HSRS template fit as the standard SRF measurement;
    - `release2_srf_models.json`;
    - `hypernet/edinterp.py` as the interpolator's home;
    - the Phase 1 scripts and tables.
  - The test count, now 92.
- **Tasks:**
  - Setup #1: new questions on the method rule and the wavelength-fit
    granularity.
  - Task 3: the `choose_method` rule (ΔFWHM ≥ 0.15 nm).
  - Task 4: built on the template fit.
  - **New task 4b:** the E wavelength shift + stretch fit (G1, ≤ 0.05 nm).
  - Task 6: the Release 2 source and the ΔFWHM/method columns.
  - **New task 9b:** the non-Gaussian SRF twin check.
  - Task 10: states the new items.
- Nothing else changed; no code.

### 2026-10-07 -- Setup #1 (Opus 5.5)

- Read the plan (§3, §4 Phase 2 and the G1 outcome, §7, §9), `gate_G1.md`,
  this doc and the Phase 0/1 modules.  Wrote 11 questions in Q&A "Setup #1",
  each with a default.  No code written.
- **Learned:**
  - The L1C carries three relative-uncertainty components in %
    (`u_rel_random_`, `u_rel_systematic_indep_`,
    `u_rel_systematic_corr_rad_irr_` + variable), with err_corr matrices
    named in `err_corr_1_params`.  The plan's two-name signature needs
    updating (Q5).
  - `edinterp` depends on `twin` through `convolve_to_grid` (Q1).
  - One `srf` call takes 0.02 s.
  - `SRFModel` already holds `cov` and `offset_err` for punpy.
  - G1 never tested an **L** wavelength error, which `srf` is exposed to
    through Emod_L (Q8).
- Read-only checks only: the VEIT L1C variable list and one timed call.


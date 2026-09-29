# WATERHYPERNET -- Wiggles Phase 0: Diagnosis

## Goal

Measure, from the data themselves, the per-channel (E, L) and per-instrument
spectral response functions (Gaussian, FWHM quadratic in λ) and wavelength
offsets of the HYPSTAR radiometers, and test whether they are stable across
SZA, sky condition, season and recalibration.  Phase 0a does this now on the
one VEIT sequence in hand and produces the reusable fitting code; Phase 0b
repeats it over the ~150 requested sequences once the L1A data arrive.

**Gate G0** (plan §4): (a) if FWHM_E(λ) = FWHM_L(λ) within the fit uncertainty
on most instruments, H2 is dropped and eq. (14) suffices; (b) if the SRFs are
stable per instrument (spread across SZA, sky and season below ~0.2 nm, no
unexplained jump at recalibration), Phase 2 ships a per-instrument SRF table,
otherwise Phase 2 must self-calibrate per sequence; (c) if centroid offsets
exceed ~0.1 nm systematically on some instrument, a wavelength recalibration
goes into the function first.

## Conventions

- `ocean14`.  Run via `conda run -n ocean14 python ...`; `conda activate` fails
  non-interactively.
- **Model.** Use Opus 5.5 (`claude-opus-5-5`) for this work, including any
  subagents (pass `model: opus`).
- **JXP runs git.**  Claude does not run any state-changing git command.
- Run scripts **from the repository root** so `hypernet` imports resolve
  (`python -m hypernet.whn_explore 1`, `python -m hypernet.wiggles.phase0a_veit`).
- Run `pytest -q` after each step where relevant.  55 tests today; the archive-
  dependent ones skip themselves when `$OS_COLOR` is not mounted.
- **Tier 2** -- steps that read the VEIT sample or the archive need `$OS_COLOR`
  mounted.  Do **not** unset `$OS_COLOR`.
- Any calculation goes into a script on disk, not into the chat.
- Data intermediates (parquet/npz) go **outside** the repo under
  `$OS_COLOR/hypernet/wiggles/phase0/`; only figures and small tables are
  committed.
- **Paths** (Q&A Setup #1, Q5; moved into the package 2026-09-29): the
  Phase 0 scripts are the subpackage `hypernet/wiggles/`
  (`hypernet/wiggles/phase0a_*.py`, `phase0b_*.py`), run as modules from the
  repository root, e.g. `python -m hypernet.wiggles.phase0a_veit`.  Their
  committed tables go in `hypernet/wiggles/phase0_*.csv` and figures in
  `hypernet/wiggles/figs/phase0/`; the path constants (`WIGGLES_DIR`,
  `FIGDIR`, `DATA_DIR`, `OUT`) live in `hypernet/wiggles/__init__.py`.  The
  earlier exploratory scripts stay in top-level `wavecal/`.  The SRF JSON
  goes in `hypernet/data/`.
- Every `.md` in `docs/` is published.  Gate reports and other unpublished
  material stay out of `docs/`.
- Ask questions in the Q&A section below; log completed work under `## Logs`.

## Context

### Where things are

- **Plan:** `docs/wiggles_planning.md` §2.3 (what the VEIT sequence says, with
  the line table), §4 "Phase 0", §5 (data request), §9 (risks: SRF stability,
  Ring effect, Gaussian shape).
- **Setup doc:** `claude_prompts/wiggles/wiggles_prompts.md` -- Q&A decisions
  (Gaussian + quadratic FWHM(λ); widths not depths; stability is the gate) and
  Logs.
- **Existing scripts** (exploratory, to be promoted into the package):
  - `wavecal/line_fits_veit.py` -- `model()` (Gaussian absorption on a linear
    continuum), `fit()` (`curve_fit` within ±4 nm), the ten-line list `lines`,
    Ld/Lu split by `viewing_zenith_angle < 90`.
  - `wavecal/sanity_checks_veit.py` -- L1C = `np.interp` check, linear-vs-cubic
    error, `hf_power()`, rigid-shift scan, `d2()` and the ρw'' correlations.
  - `wavecal/instrument_timeline.py`, `wavecal/select_sequences.py` -- the
    instrument periods and the data request.
- **Package:** `hypernet/whn_explore.py` (`whn_root()`, `out_root()`,
  `load_spectrum`) is the model for path resolution and Tier-2 tests;
  `hypernet/tests/test_whn_explore.py` shows the `needs_whn` skip pattern.
- **Data:**
  - VEIT sample (in hand): `$OS_COLOR/WATERHYPERNET/Wavelengths/` --
    `HYPERNETS_W_VEIT_L1A_IRR_20260604T0845_*.nc`, `..._L1A_RAD_...`,
    `..._L1C_ALL_..._090_...`, `..._L2A_REF_..._090_...`.  Variables used so
    far: `wavelength`, `irradiance` (wavelength, scan), `radiance`,
    `viewing_zenith_angle`, `downwelling_radiance`, `upwelling_radiance`,
    `reflectance`, `reflectance_nosc`.
  - Full request: `docs/wiggles_data_request.csv` (224 rows; columns `row,
    priority, optional, site, instrument, cal_dates_rad_irr, sequence_time,
    azimuth, sza, sky, file`).  **Not yet delivered.**  Expected as
    L1A_IRR/L1A_RAD/L1C_ALL/L2A_REF per sequence plus the cal files
    `HYPERNETS_CAL_HYPSTAR_*_{RAD,IRR}_v2.3.nc`, mirrored from
    `AIOcean:data/Color/WATERHYPERNET/Wavelengths/` (rclone) to
    `$OS_COLOR/WATERHYPERNET/Wavelengths/`.
  - `$OS_COLOR/hypernet/wavecal/instrument_timeline.csv` -- instrument periods.
  - `$OS_COLOR/WATERHYPERNET/RELEASE_2` -- L2A/L2B only; no L1A.
- **Paper:** `context/papers/ruddick2023.pdf`.

### Dependencies

- 0a needs nothing beyond the VEIT sample.  Its FWHM(λ) curves feed the
  Phase 1 instrument model, so finish tasks 2-8 early.
- 0b (tasks 9-13) and the gate (task 14) wait on the data.
- Naming clash: the exploratory scripts live in top-level `wavecal/`, so a
  package module `hypernet/wavecal/` would shadow it in the mind if not on the
  path.  Candidate names: `hypernet/srf.py` (SRF model + line fits) and
  `hypernet/whn_l1a.py` (L1A/L1C readers).  JXP chooses in Setup Q&A.

## Prompts

### Setup

1. Read the plan sections and the setup doc listed in Context, this doc, and
   `wavecal/line_fits_veit.py` and `wavecal/sanity_checks_veit.py` in full.
   Then put your questions in the Q&A section below: at least the module names
   (`hypernet/srf.py` and `hypernet/whn_l1a.py`, or a subpackage), whether
   Ca H/K is fitted as a joint two-Gaussian blend, whether the 936 nm H₂O band
   stays in the line list, and where committed figures and small tables go
   (proposal: `wavecal/figs/phase0/` and `hypernet/data/`).  Do not write code
   yet.  Log your work.

1b. See my answers to your questions in the Q&A section below and respond accordingly.
   Then move on to task #2 under Build.  Log your work.

### Build (0a, VEIT sample)

2. **Line fitting module.**  Promote `line_fits_veit.py` into the package module
   agreed in Q&A (default `hypernet/srf.py`): a `LINES` table (name, λ_lab,
   half-window, blend group) for the ten lines of plan §2.3 with Ca H/K as a
   blended pair; `fit_line(wav, spec, lam, half=4.0)` returning centroid,
   sigma, depth, continuum and their `curve_fit` uncertainties; `fit_lines(wav,
   spec, lines=LINES)` returning a DataFrame with FWHM = 2.3548 σ.  Tests in
   `hypernet/tests/test_srf.py`: a synthetic Gaussian line on a sloped
   continuum recovers centroid and FWHM to 0.01 nm; the blended pair recovers
   both centroids; a line that fails returns NaN, not an exception.  Log your
   work.

3. **FWHM(λ) model.**  Add `fit_fwhm_model(lam, fwhm, err)` (weighted quadratic
   least squares, returning coefficients and covariance), `fwhm_at(lam,
   coeffs)`, and an `SRFModel` dataclass (channel, coefficients, covariance,
   centroid offset and its error, wavelength range) with `to_json`/`from_json`.
   Tests: an exact quadratic is recovered; the covariance shrinks with smaller
   `err`; JSON round-trip.  Log your work.

3b. **Template (HSRS) SRF fit** (Q&A Setup #1, Q2).  Fetch TSIS-1 HSRS
   (brought forward from Phase 1; record URL and checksum in the Logs) and add
   `fit_srf_template(wav, spec, hsrs, ...)` to the SRF module: per line
   window, fit σ, shift, scale and continuum against HSRS convolved with the
   Gaussian SRF (O2 bands only with HAPI O2).  Its `SRFModel` is the one
   written to `hypernet/data/veit_srf_model.json`; the empirical fit of tasks
   2-3 stays as a cross-check.  Tests: a synthetic spectrum made from HSRS
   with a known σ recovers it.  Log your work.

4. **L1A/L1C readers.**  New module (default `hypernet/whn_l1a.py`):
   `wavelengths_root()` resolving `$OS_COLOR/WATERHYPERNET/Wavelengths`;
   `load_l1a_irr(path)` (grid, per-scan and mean Ed); `load_l1a_rad(path)`
   splitting Lu (water view, `viewing_zenith_angle < 90`) and Ld (sky,
   ≥ 90; vza is measured from nadir -- see Q&A Setup #1), per scan and
   mean, each checked against the L1C variables; `load_l1c(path)` and `load_l2a(path)` returning `irradiance`,
   `downwelling_radiance`, `upwelling_radiance`, `reflectance`,
   `reflectance_nosc` and the geometry.  Tier-2 tests with a `needs_wavelengths`
   skip: 1536/1538 pixels, L1C `irradiance` equals `np.interp` of the L1A mean
   to 1e-4, 6 scans in L1A_IRR and L1C, 12 in L1A_RAD (6 Ld + 6 Lu).  Log
   your work.

5. **VEIT fits, with uncertainties.**  `hypernet/wiggles/phase0a_veit.py`: fit E from
   the mean L1A_IRR and L from the mean Ld and Lu separately, then fit
   FWHM(λ) per channel.  Write the line table to
   `$OS_COLOR/hypernet/wiggles/phase0/veit_lines.parquet`, a committed copy to
   `hypernet/wiggles/phase0_veit_lines.csv`, the three `SRFModel`s to
   `hypernet/data/veit_srf_model.json` (small, committed -- Phase 1 reads it),
   and a figure `hypernet/wiggles/figs/phase0/veit_fwhm_vs_lambda.png` (FWHM points
   with error bars, quadratic fits, E vs Ld vs Lu).  Reproduce the plan §2.3
   table to the quoted precision and note any differences in the Logs.  Log
   your work.

6. **Ld = Lu check and per-scan scatter.**  Extend the script (or add
   `hypernet/wiggles/phase0a_scans.py`): fit every scan separately (6 E scans, the Ld
   and Lu scans), report the scan-to-scan spread of centroid and FWHM per line
   against the `curve_fit` errors, and test Ld = Lu per line (difference over
   combined error).  Line depths in E, Ld and Lu are tabulated as a
   diagnostic only (Ring effect).  Figure
   `hypernet/wiggles/figs/phase0/veit_scan_scatter.png`.  Log your work.

7. **Preliminary H1/H2 budget.**  `hypernet/wiggles/phase0a_budget_veit.py`: per line,
   (H1) the residual line amplitude in Ld/Ed produced by linear interpolation
   alone (linear minus cubic on the same Ed, as in `sanity_checks_veit.py`) and
   (H2) the residual predicted from the FWHM difference,
   depth × (FWHM_L² − FWHM_E²)/FWHM², against the measured residual line
   amplitude in Ld/Ed and the ρw'' power within ±5 nm of the line.  Write
   `hypernet/wiggles/phase0_veit_budget.csv` and
   `hypernet/wiggles/figs/phase0/veit_budget.png`.  This script is the template for
   task 12.  Log your work.

8. **0a interim report.**  In Q&A, summarise: the FWHM(λ) coefficients and
   uncertainties handed to Phase 1, whether E ≠ L survives the uncertainties,
   the centroid offsets against the 0.1 nm threshold of G0(c), and anything
   that argues against the Gaussian.  Log your work.

8b. **0a interim report slides.**  Create a slideset in `docs/slides/wiggles_phase0_report.pptx` with a description of what you have accomplished so far.  Include the figures and tables from tasks 5-7.  Log your work.  When using text, use no font size smaller than 20pt

### Build (0b, full request -- after the data arrive)

9. **Ingest and index.**  Mirror the delivery with `rclone` to a subfolder of
   `$OS_COLOR/WATERHYPERNET/Wavelengths/` (name per the delivery).
   `hypernet/wiggles/phase0b_index.py`: match every L1A/L1C/L2A file to a row of
   `docs/wiggles_data_request.csv` by site, `sequence_time` and azimuth,
   carry `instrument`, `sza`, `sky` and the cal period, record which
   requested sequences are missing and which spares were substituted, and
   check the `instrument_calibration_file_rad` attribute bug.  Write
   `$OS_COLOR/hypernet/wiggles/phase0/request_index.parquet` and a committed
   summary `hypernet/wiggles/phase0_delivery_summary.csv`.  Log your work.

10. **Batch line fits.**  `hypernet/wiggles/phase0b_fit_all.py`: `fit_lines` on the
    mean E, Ld and Lu of every sequence (per scan optional, behind a flag) and
    `fit_fwhm_model` per sequence and channel.  Outputs
    `$OS_COLOR/hypernet/wiggles/phase0/line_fits.parquet` and
    `srf_models.parquet` keyed by site, instrument, cal period, sequence.
    Add a Tier-2 test that one delivered sequence fits without error.  Log
    your work.

11. **Stability.**  `hypernet/wiggles/phase0b_stability.py`: FWHM(λ) at 400, 500, 600,
    700 nm and the centroid offset (L − E) versus SZA, `sky`, month and time,
    per instrument; the spread per instrument against the 0.2 nm criterion;
    122304 and 121222 before/after recalibration; 122302 at BEFR vs THFR
    (site vs instrument); the two-instruments-one-site pairs at VEIT, GAIT and
    MAFR.  Figures under `hypernet/wiggles/figs/phase0/` and a committed table
    `hypernet/wiggles/phase0_stability.csv`.  Log your work.

12. **H1/H2 error budget per sequence.**  Generalise task 7 into
    `hypernet/wiggles/phase0b_budget.py` over the whole index: per sequence and line,
    the H1 and H2 contributions versus the measured ρw'' excess at the line;
    summarise by instrument and water type.  Table and figure as above.  Log
    your work.

13. **Calibration files and lab data.**  `hypernet/wiggles/phase0b_calfiles.py`: open
    each `HYPERNETS_CAL_HYPSTAR_*_{RAD,IRR}_v2.3.nc`, list the variables and
    attributes (wavelength coefficients, `bandwidth`, anything SRF-like), and
    compare with the fitted centroids and FWHM(λ); read any lab line-spread
    data Kevin sends.  Record in the Logs what the cal files do and do not
    contain.  Log your work.

### Gate

14. **Gate G0.**  Evaluate (a), (b) and (c) against the plan's criteria using
    the outputs of tasks 10-13, with numbers.  Write a short gate report
    (proposed location `claude_prompts/wiggles/gate_G0.md`; JXP to confirm --
    not under `docs/`) covering: the verdict on each of (a)-(c), the
    per-instrument SRF summary, what it means for Phase 2 (table vs
    self-calibration, wavelength recalibration or not), and open concerns.
    Propose the corresponding edits to `docs/wiggles_planning.md` §4 in Q&A,
    and apply them only if JXP approves.  Log your work.

## Q&A

### Setup #1 -- 2026-09-28 (Opus 5.5)

**What I found reading the inputs**

1. **The Ld/Lu split in the scripts is reversed.**  `line_fits_veit.py` sets
   Ld = `viewing_zenith_angle < 90`, but the convention is measured from
   nadir.  L1A_IRR, which looks up, has vza = 179.9.  In L1A_RAD the 6 scans
   at vza = 40.19 are the **water** view (Lu).  The 6 at 140.03 / 139.94 are
   the **sky** (Ld), taken as 3 before and 3 after the water scans.
   `wavecal/check_vza_convention_veit.py` shows it:
   - The vza < 90 mean equals L1C `upwelling_radiance` exactly (median
     relative difference 0).
   - The vza ≥ 90 mean matches L1C `downwelling_radiance` to 2 × 10⁻³.
   - The sky is 5-11x brighter than the water view.

   Consequences:
   - The "FWHM Ld" and "FWHM Lu" columns of the plan §2.3 table are swapped.
   - The "centroid L − E" column came from the water view, not the sky.
   - The conclusion stands (E sharper than L in the blue, the two L views
     agree), but the labels are wrong.
   - The same slip is in task 4 of this doc, in the Phase 1 doc (Context,
     "Ld at `viewing_zenith_angle` ≈ 40°, Lu at ≈ 140°") and possibly in the
     note to Kevin, which does not quote vza.
   - L1A_RAD has **12** scans, not 6, so the task 4 test needs changing.
   - `sanity_checks_veit.py` is unaffected: it takes Ld and Lu from the L1C
     variables.
2. **A Gaussian fitted to a solar feature measures the SRF convolved with the
   line, not the SRF itself.**  Several entries in the list are not single
   narrow lines:
   - Ca H/K have broad cores and damping wings reaching past the ±4 nm
     window, so the linear "continuum" sits inside the wing.
   - The G band is a CH molecular band.
   - Mg b is a triplet (516.7 / 517.3 / 518.4 nm).
   - Na D is a doublet 0.6 nm apart.
   - Hα and Hβ have wide wings.
   - O2-A and O2-B are rotational bands.  Their envelope shape depends on air
     mass, which differs between the direct beam and the sky path.

   For the E-vs-L and stability tests (G0 a/b) the intrinsic profile is common
   to both channels and largely cancels.  For the **absolute** FWHM(λ) handed
   to Phase 1 it does not: the empirical widths are biased high by a
   line-dependent amount, and a quadratic will absorb that bias as spurious
   shape.  The clean fix is a forward-model fit: convolve TSIS-1 HSRS (plus
   HAPI O2 for the bands) with a Gaussian SRF, then fit σ, shift, scale and
   continuum.  Blends and bands come out right automatically.  The cost is
   bringing the HSRS download forward from Phase 1.
3. **Air vs vacuum.**  The λ_lab values in the list are air wavelengths
   (Hα 656.28), whereas TSIS-1 HSRS is on a vacuum scale.  The difference is
   0.1-0.2 nm, the same size as the G0(c) threshold.  We do not know which
   scale the HYPSTAR wavelength calibration uses.
4. **G0(c) is ambiguous.**  "Centroid offsets exceed ~0.1 nm" could mean L − E
   (relative, which drives H1) or measured − laboratory (absolute, which
   drives a recalibration).  An absolute centroid means little on a blended
   feature.
5. **The red lever arm is thin.**  Above 660 nm the list has only Hα, the two
   O2 bands and the 936 nm H₂O band, and none of the bands is a clean SRF
   probe.  The Ca II infrared triplet (849.8 / 854.2 / 866.2 nm) is strong,
   solar, fairly free of telluric lines and inside the 350-1100 nm range.
   The water view is faint in the red (3.0 vs 34 at 650 nm), so red fits in
   Lu will be noisier.
6. **The uncertainties are not yet real.**  `curve_fit` with no `sigma`
   rescales the covariance by the residuals, so model mismatch leaks into the
   errors, and L1A carries no uncertainty variables.  Task 6 compares scan
   scatter against these errors, so the choice matters.
7. **Packaging.**  `setup.py` has no `package_data`, so `hypernet/data/*.json`
   would not ship with a non-editable install.

**Questions for JXP** (each has a default; say "defaults" to accept them all)

1. **Module names.**  *Default: `hypernet/srf.py` (SRF model + line fits,
   instrument-agnostic, hence no `whn_` prefix) and `hypernet/whn_l1a.py`
   (L1A/L1C/L2A readers).  No subpackage.*
>A. Defaults
2. **Empirical Gaussian or HSRS forward model?**  *Default: both.*
   - Tasks 2-3 as written (empirical Gaussian): fast, and adequate for the
     relative tests (E vs L, stability).
   - A task 3b: `fit_srf_template(wav, spec, hsrs, ...)`, fitting σ(λ) and the
     shift against HSRS convolved with the SRF.  Its `SRFModel` is the one
     written to `hypernet/data/veit_srf_model.json` for Phase 1.
   - The empirical fit stays as a cross-check.

   This pulls the HSRS download (a Phase 1 task) into Phase 0.
>A. Defaults
3. **Ca H/K.**  *Default: a joint fit over ~389-401 nm with two Gaussians,
   separate σ and depth, one linear continuum.  Both lines enter the
   FWHM(λ) fit but are flagged `blend`.  With the template fit (Q2) the issue
   goes away.*
>A. Defaults
4. **Line list.**  *Default:*
   - *936 nm H₂O stays as a diagnostic row (`use_for_srf=False`), out of the
     FWHM(λ) fit.*
   - *O2-A and O2-B are treated the same way in the empirical fit, and used in
     the template fit only with HAPI O2.*
   - *Add the Ca II infrared triplet (849.8, 854.2, 866.2 nm) to anchor the
     red.*
   - *The `LINES` table carries `lam_air`, `lam_vac`, `half`, `group` and
     `use_for_srf`.*
>A. Defaults
5. **Where outputs go.**  *Default: figures in `wavecal/figs/phase0/`, small
   CSVs in `wavecal/`, and the SRF JSON in `hypernet/data/` with
   `package_data={'hypernet': ['data/*.json']}` added to `setup.py`.*
>A. Use `wiggles` instead of `wavecal`
6. **Fit uncertainties.**  *Default:*
   - *Per-pixel σ = scan-to-scan std/√N for the mean spectra.*
   - *For single scans, σ = that std × √(N/(N−1)).*
   - *Fit with `absolute_sigma=True`.*
   - *Report χ²_ν per line: a large χ²_ν is the first sign that the Gaussian
     is the wrong shape (plan §9).*
>A. Defaults
7. **Fixing the Ld/Lu slip.**  *Default:*
   - *The reader in task 4 uses Lu = vza < 90 and Ld = vza ≥ 90, and checks
     both against L1C.*
   - *The task 4 test expects 6 IRR scans and 6 + 6 RAD scans.*
   - *Correct the Phase 1 Context line.*
   - *Leave `line_fits_veit.py` as a historical exploratory script, with a
     one-line comment pointing to the fix.*
   - *Correct the §2.3 column headers in `docs/wiggles_planning.md` after
     task 5 re-derives the table.*

   `docs/` is published, so I will not touch it without your OK.
>A. Defaults
8. **Air or vacuum, and G0(c).**  *Default:*
   - *Ask Kevin which scale the HYPSTAR wavelength calibration uses; assume
     air until he answers.*
   - *G0(c) tests both the relative (L − E) and the absolute (vs laboratory)
     offset against 0.1 nm.*
   - *The absolute offset is tested on unblended lines only (Hα, Hβ, Ca II
     IR), or on the template-fit shift.*
>A. Defaults

### Setup #1b -- 2026-09-28 (Opus 5.5): what I did with the answers

All defaults were accepted, except that Q5 uses `wiggles` instead of
`wavecal`.

- **Q5.**
  - Read as: new Phase 0 scripts, figures and small CSVs go in a new
    top-level `wiggles/` (`wiggles/phase0a_*.py`, `wiggles/figs/phase0/`,
    `wiggles/phase0_*.csv`).
  - The existing exploratory scripts stay in `wavecal/`, and
    `$OS_COLOR/hypernet/wavecal/` is untouched.
  - The task paths above and the Conventions now say so.
  - If you meant "rename `wavecal/` to `wiggles/`", that is a `git mv` for you
    (plus path edits in the setup doc, the plan and the Phase 1-3 docs); say
    so and I will do the edits.
- **Q2.**  Added task 3b (HSRS template fit) to Build.
- **Q7.**
  - Task 4 now has the correct Lu/Ld split and scan counts.
  - The Phase 1 Context line is corrected.
  - `wavecal/line_fits_veit.py` has a comment flagging the swapped labels.
  - The §2.3 table in `docs/wiggles_planning.md` waits for task 5.
- **Q8.**  Added the air/vacuum question as item 5 of the requests in
  `correspondence/wiggles_note_to_kevin.md`, which is still unsent.
- **Q5, packaging.**  Added `package_data={'hypernet': ['data/*.json']}` to
  `setup.py`.
- **Q1, Q3, Q4, Q6.**  Implemented in `hypernet/srf.py` (task 2, see Logs).
  - One extension to Q4: Ca II 849.8 and 854.2 nm are only 4.4 nm apart, so
    their ±4 nm windows overlap.  They are fitted jointly, as Ca H/K are
    (group `CaIR12`), but are not flagged `blend`.
  - `group` means "fitted jointly".  `blend` means "the feature is itself a
    blend or band", i.e. unsuitable for absolute centroids.

### Build #4 -- 2026-09-28 (Opus 5.5): what the VEIT files say, and one question

From `wiggles/phase0a_l1c_consistency.py`:

- **How L1C is built.**
  - Lu is the six L1A water scans, unchanged.
  - Ld is **the two sky-series means interpolated linearly in time** to the
    water-view time: it matches to 4.5 × 10⁻⁵, whereas the plain 6-scan mean
    is off by 3 × 10⁻³.
  - E also has two series (before and after); here they bracket the water
    view symmetrically (weight 0.500), so the time interpolation equals the
    plain mean.
  - L1C `irradiance` = `np.interp` of that E to ≤ 4 × 10⁻⁵ over 400-900 nm.
    It departs at low-signal pixels: 2 × 10⁻⁴ at the 935 nm H₂O band and
    2-3 % within 10 nm of the grid ends.  This suggests the processor handles
    some pixels differently when averaging (e.g. outlier rejection); it
    doesn't matter for the lines.
- **No usable uncertainties in L1C.**  `u_rel_random_*` and
  `std_downwelling_radiance` are all zero.
- **The scan-to-scan scatter is mostly broadband.**  Relative std/mean,
  median over 400-900 nm:

  | | raw | after removing a smooth (20 px) scale per scan |
  |---|---:|---:|
  | E | 6.3e-3 | 4.1e-3 |
  | Lu | 5.1e-2 | 7.5e-3 |
  | Ld (both series) | 5.9e-2 | 4.3e-3 |
  | Ld (within one series) | 3.7-4.0e-3 | 3.7-3.9e-3 |

  The sky changed by about 5 % between its two series, and the water view
  jitters by about 5 % scan to scan (glint and waves).  A line fitted on a
  free continuum sees only the flattened scatter.

**Question 9 -- revise Q6?**  The agreed Q6 default (raw std/√N) would
overstate the per-pixel errors by 7x for Lu and 14x for Ld.  That makes
χ²_ν ≪ 1 and the `curve_fit` errors too large, which would bias the task 6
comparison of scan scatter against fit errors.
*Default: use `srf.scan_errors(scans, flatten_px=20)` (new option; the mean
is unchanged) for the errors in tasks 5, 6 and 10, and report the raw scatter
alongside.*

>A. Ok, we should avoid overstating the errors for the reasons you mentioned.
Let's go with your suggestion and we may revisit it later.

### Build #8 -- 2026-09-28 (Opus 5.5): Phase 0a interim report

One sequence: VEIT, 2026-06-04 08:45, HYPSTAR 122304, post-recalibration
(IRR 2024-11-12, RAD 2024-11-08).  The numbers come from
`wiggles/phase0a_interim.py`, which reads the committed task 5-7 products.
The fits use flattened scan-scatter errors (Q9).  Task 6 showed that the scan
spread matches these errors (median ratio 0.9).

**1. FWHM(λ) handed to Phase 1** (`hypernet/data/veit_srf_model.json`)

The model is FWHM = c₀ + c₁x + c₂x², with x = (λ − 600 nm)/100 nm.  It is fitted
to 10 lines (Ca H/K, G band, Hβ, Mg b, Na D, Hα and the Ca II IR triplet)
over 393-866 nm.  The covariance is inflated by χ²_ν.

| channel | c₀ (nm) | c₁ | c₂ | χ²_ν | FWHM at 400 / 500 / 600 / 700 / 850 nm |
|---|---|---|---|---:|---|
| E | 2.80 ± 0.12 | +0.03 ± 0.04 | −0.02 ± 0.03 | 26 | 2.66 / 2.75 / 2.80 / 2.80 / 2.73 (± 0.06-0.14) |
| Ld | 3.17 ± 0.11 | −0.04 ± 0.04 | −0.04 ± 0.03 | 13 | 3.07 / 3.16 / 3.17 / 3.08 / 2.79 (± 0.05-0.19) |
| Lu | 3.28 ± 0.11 | −0.10 ± 0.06 | −0.08 ± 0.04 | 5 | 3.16 / 3.30 / 3.28 / 3.11 / 2.54 (± 0.06-0.32) |

Caveats for Phase 1:
- **These are Gaussian widths of solar features, not SRF widths.**  They are
  biased high by the intrinsic line widths, by different amounts per line;
  that is what χ²_ν = 5-26 measures.  The per-line scatter about the curves
  is 0.2-0.3 nm.
- **Only the E − L difference is reliable.**  The absolute level (the SRF
  itself) needs the HSRS template fit of task 3b (see Q10).
- **Curvature is poorly constrained.**  Every c₂ is within 2σ of zero, and
  above ~700 nm the curves rest on Hα and the Ca II IR triplet only.  Use
  them only inside 393-866 nm.

**2. Does E ≠ L survive the uncertainties?  Yes, in the blue and green; not
in the red.**

- **Per line.**  FWHM_Ld − FWHM_E is 0.25-0.59 nm at every SRF line from
  Ca K to Na D, each at 5-15σ.  At Hα it is 0.16 nm (2.7σ), and at the
  Ca II IR triplet 0.09-0.16 nm (< 1.5σ).  Lu − E is similar: 0.32-0.74 nm
  at 5-12σ from Ca K to Na D, 0.31 at Hα, and ≈ 0 ± 0.2 nm in the Ca II IR.
- **Model level.**  Ld − E is 0.41-0.42 ± 0.09 nm at 400-500 nm (4-5σ),
  0.37 ± 0.16 at 600 nm (2.3σ) and 0.14 ± 0.19 at 800 nm.  The difference
  fades to the red.  Model-level significance is lower than per-line, because
  the χ²_ν inflation charges both curves for the intrinsic-width scatter,
  which largely cancels in the difference.
- **Hedge.**  Of the lines where E ≠ L is strongest, five of seven are blends
  (Ca H/K, G band, Mg b, Na D).  Among the unblended lines, only Hβ is
  individually decisive (0.49 nm, 7.5σ; Hα 2.7σ).  The difference is common
  to blends and clean lines, so it is not a blend artefact.  But a clean-line
  confirmation at more S/N (more sequences, task 10) or through the template
  fit would make it airtight.
- **The residual structure matches H2** (task 7).  α₂ = 0.93 ± 0.06 in Ld/Ed
  and 1.08 ± 0.14 in ρw.  H1 has the wrong sign and is ~9× too small.  This
  check is partly circular; see the Build #7 log.
- **One L SRF.**  Ld ≈ Lu to −0.06 ± 0.03 nm at the clean lines (task 6).
  The exceptions are Ca H and the G band, where Lu is 0.16-0.33 nm wider,
  plausibly from Lw content.  Ld is the better L-channel probe (Q11).
- **Stable within the sequence.**  The two sky series, ~50 s apart, agree to
  |z| < 1.7.

**Verdict:** G0(a) is not met on this sequence, i.e. H2 stays.  It is one
sequence; task 11 decides.

**3. Centroid offsets vs the 0.1 nm threshold of G0(c)**

- **Relative (L − E), which drives H1: below threshold.**
  - Weighted mean over the clean lines (Hβ, Hα, Ca II IR): Ld − E = **+0.050
    ± 0.011 nm** (χ²_ν 0.9, max 0.075 nm) and Lu − E = +0.038 ± 0.018 nm.
  - Including the blends: +0.055 ± 0.015 nm.  Only Ca K exceeds 0.1 nm
    (+0.131 ± 0.013), a blend.
  - The offset is significant (4.5σ) but half the threshold.  It agrees with
    the joint-fit shift of task 7 (−0.04 ± 0.01 nm in that sign convention).
- **Absolute (measured − lab), which would drive a recalibration:
  undetermined.**
  - Weighted means over the clean lines: +0.15 ± 0.05 nm against air and
    −0.05 to 0.00 nm against vacuum, the same in all three channels.
  - The lines disagree by ±0.15 nm, i.e. 10× their errors: Hα sits on the
    air scale (−0.02), while Hβ (+0.20) and the Ca II IR (+0.10 to +0.37)
    sit closer to vacuum.
  - So the 0.1 nm test can't be applied until Kevin says which scale HYPSTAR
    uses (note to Kevin, item 5) and the template fit gives unbiased
    centroids.
  - If the scale is air, E, Ld and Lu are all offset by about +0.15 nm, which
    exceeds the threshold.  But they are offset *together*, so it would not
    feed H1.
- **Verdict:** G0(c) relative: no recalibration is needed for H1.  G0(c)
  absolute: open.

**4. Evidence against the Gaussian**

- Per-scan χ²_ν (task 6): clean lines 1.2-1.8, blends 1.6-3.2, bands 30+.
  The Gaussian is adequate for isolated lines at single-scan S/N.  It fails
  for features that aren't single lines, which says more about the line
  model than about the SRF shape.
- On the mean spectra (√6 more S/N), χ²_ν per line is 3-7 (task 5), so even
  clean lines show structure beyond a Gaussian-on-a-line.  This can't yet
  separate a non-Gaussian SRF from non-Gaussian intrinsic profiles (Hα and
  Hβ have broad wings, Ca II IR has damping wings).
- There are no obvious asymmetric residuals in the task 7 profiles.  The
  measured Ld/Ed bumps match the H2 shape, apart from Na D (weaker,
  α₂ ≈ 0.3-0.9) and Hα (stronger, α₂ ≈ 1.2-1.5).
- **Nothing here argues against a Gaussian SRF.**  The direct test is the
  template fit (task 3b): with the true line profiles supplied by HSRS, the
  residuals left over are the SRF shape.

**Other findings for Kevin / later phases**

- L1C `u_rel_random_*` are zeros, and L2A `std_reflectance` overstates the
  per-pixel noise (it includes broadband scan changes).
- Lu appears to carry an additive, line-free component in the red (EW_Lu/EW_E
  = 0.64-0.85 at O2-B and the Ca II IR, against ≈ 1 in Ld): stray light or a
  dark residual?
- In the red the ρw structure at the lines is not Ed-borne (task 7).

**Question 10 -- template fit before Phase 1?**  The widths handed over are
line widths, not SRF widths, so a Phase 1 twin experiment built on them would
use an L SRF too wide by an unknown ~0.1-0.5 nm per channel.  The E − L
difference is right to first order.
*Default: do task 3b (HSRS template fit) next, before Phase 1 reads the JSON,
and overwrite `veit_srf_model.json` with the template models.  Keep the
empirical file as `hypernet/data/veit_srf_model_empirical.json`.  If Phase 1
needs to start first, it uses the empirical file and flags every result as
provisional.*
>A. Follow your default

**Question 11 -- which channel is "L"?**  Ld and Lu share one SRF to ~0.1 nm,
but Lu is noisier, faint in the red, and shows possible Lw/Raman and
additive-light effects.
*Default: Ld is the reference L channel for the SRF (tasks 10-11 and
Phase 2), Lu is carried as a check, and the gate compares E with Ld.*
>A. Follow your default.

### Build #3b -- 2026-09-29 (Opus 5.5): update to the interim report

The HSRS template fit (task 3b) supersedes items 1, 3 and 4 of the Build #8
report.  Item 2 stands, and is now stronger.  Everything comes from
`wiggles/phase0a_template_veit.py`.

- **SRF handed to Phase 1** (`hypernet/data/veit_srf_model.json`, template
  fit, veil 0).

  | | 400 | 500 | 600 | 700 | 850 nm |
  |---|---:|---:|---:|---:|---:|
  | E | 2.35 | 2.23 | 2.21 | 2.28 | 2.55 |
  | Ld | 2.88 | 2.74 | 2.64 | 2.58 | 2.57 |
  | Lu | 2.92 | 2.73 | 2.63 | 2.62 | 2.76 |

  - Errors are ± 0.04-0.15 nm, with the covariance inflated by χ²_ν = 12-21.
  - The windows are 27 × 10 nm over 390-880 nm, with the telluric bands
    excluded, so the curves are **constrained over 390-680 and 850-870 nm
    only**.
  - The empirical line widths were 0.3-0.5 nm too wide: E was 2.7-2.8 nm.
- **E ≠ L, now without the blend caveat.**  Ld − E is 0.53 ± 0.04 nm at
  400 nm, 0.50 ± 0.06 at 500, 0.42 ± 0.08 at 600, 0.30 ± 0.09 at 700 and
  0.02 ± 0.15 at 850 nm.  16 of 27 windows are individually at z > 3.
- **Centroid offsets.**
  - Relative Ld − E: +0.045 ± 0.011 nm, which confirms the empirical +0.050.
  - Absolute, as robust medians over the windows:
    - E: −0.02 nm against air, −0.16 against vacuum;
    - Ld: +0.09 / −0.09;
    - Lu: +0.06 / −0.08.

    **E points to an air-scale calibration**, and on the air scale every
    channel is within 0.1 nm.  The L channels can't tell air from vacuum.
    Kevin's answer still settles it.
- **Gaussian shape.**
  - The additive veil (Ring filling-in or stray light) is zero within
    errors, at −0.03 ± 0.05 in every channel.
  - The per-window χ²_ν is 2.5-4.3.  That is not a perfect fit, but most of
    the misfit sits in windows near weak telluric lines (580-600 and
    660-680 nm).
  - A direct test of the SRF shape (e.g. a Gaussian + Lorentzian SRF) is a
    possible extension, not yet done.

## Logs

### 2026-09-28 -- Setup #1 (Opus 5.5)

- Read this doc, plan §2.3/§4/§5/§9 of `docs/wiggles_planning.md`, the setup
  doc `wiggles_prompts.md` (Q&A and Logs), and `wavecal/line_fits_veit.py` and
  `wavecal/sanity_checks_veit.py` in full.  Inspected the VEIT sample's
  dimensions and geometry:
  - L1A_IRR: 1536 px × 6 scans.
  - L1A_RAD: 1538 px × 12 scans.
  - L1C: 6 scans, vza = 40.19.
  - L2A: 1 series.
  - Both grids span 350-1100 nm.
- Found that the Ld/Lu split used so far is reversed: vza = 40 is the water
  view and vza = 140 is the sky.  Wrote
  `wavecal/check_vza_convention_veit.py` to show it against the L1C
  variables.  This is the only file added; no package code was written.
- Wrote 7 findings and 8 questions with defaults to Q&A.  They cover:
  - module names;
  - empirical vs HSRS-template SRF fitting;
  - the Ca H/K blend;
  - the line list (936 nm, O2 bands, Ca II IR triplet);
  - output locations and `package_data`;
  - fit weighting;
  - fixing the Ld/Lu slip;
  - air/vacuum and the meaning of G0(c).
- Learned:
  - `setup.py` has no `package_data`.
  - `sanity_checks_veit.py` reads Ld/Lu from L1C, so it is unaffected by the
    swap.
  - The Phase 1 prompt doc repeats the swapped vza labels in its Context
    section.
- No pytest run (no package code changed).

### 2026-09-28 -- Setup #1b and Build #2 (Opus 5.5)

- Applied JXP's Setup #1 answers (all defaults; Q5 uses `wiggles/`); details
  in Q&A Setup #1b.
  - Edited this doc: Conventions, the task 4 spec, new task 3b, and
    `wavecal/` → `wiggles/` for the new outputs in tasks 5-13.
  - Edited `wiggles_phase1_prompts.md` (the geometry line),
    `wavecal/line_fits_veit.py` (a comment), `setup.py` (`package_data`) and
    `correspondence/wiggles_note_to_kevin.md` (the air/vacuum question).
- **Task 2.**  New module `hypernet/srf.py`:
  - `LINES`: 13 rows, columns `name, lam_air, lam_vac, half, group, blend,
    use_for_srf`.  They are the ten lines of plan §2.3 plus the Ca II IR
    triplet.  Ca H/K and Ca II 849.8/854.2 are joint-fit groups.  H₂O and
    O2-A/B have `use_for_srf=False`.
  - `air_to_vac()`: the inverse Morton (2000) relation, as in VALD.
  - `scan_errors()`: the mean spectrum plus its error s/√N and the
    single-scan error s (ddof = 1).
  - `fit_blend()`: n Gaussians times a linear continuum, fitted by
    `curve_fit` with `absolute_sigma=True` when `err` is given.  It returns
    `chi2_nu` and gives NaN rows with `ok=False` on failure (too few pixels,
    no convergence, bad covariance, result outside the window, or depth /
    width out of range).
  - `fit_line()`: the single-line wrapper.
  - `fit_lines()`: a DataFrame with FWHM = 2.3548 σ and
    `dmu = mu − lam_air`.
  - Two differences from `line_fits_veit.py`:
    - the continuum is referenced to a fixed wavelength, not the fitted μ,
      which decouples the parameters;
    - `c0` is reported at the nominal line centre.
- New `hypernet/tests/test_srf.py`, 8 tests:
  - air → vacuum conversion against NIST Hα and Na D2;
  - the table schema;
  - a sloped-continuum single line recovering centroid and FWHM to 0.01 nm;
  - a weighted noisy fit whose errors and χ²_ν are consistent;
  - the Ca H/K-like pair recovering both centroids and σ to 0.01 nm;
  - NaN on failure (off-grid, all-NaN, emission);
  - `fit_lines` on a grid truncated at 800 nm, giving NaN red rows;
  - `scan_errors`.
- `pytest -q`: **29 passed** (21 + 8).
- Not yet run on the VEIT sample; that is task 5, after the readers in task
  4.

### 2026-09-28 -- Build #3 (Opus 5.5)

- Added the FWHM(λ) model to `hypernet/srf.py`:
  - `fit_fwhm_model(lam, fwhm, err, deg=2, lam_ref=600, lam_scale=100,
    scale_cov=False)`.
    - A weighted least-squares polynomial in x = (λ − 600 nm)/100 nm, so
      `coeffs[0]` is the FWHM at 600 nm and the fit is well conditioned.
    - Returns `(coeffs, cov)` with cov = (AᵀWA)⁻¹, taking `err` as absolute.
    - `scale_cov=True` inflates cov by χ²_ν when χ²_ν > 1.  This is useful
      because blend-biased widths will scatter more than their errors.
    - Drops non-finite points and points with error ≤ 0, and returns NaN when
      fewer than deg + 1 points remain.
  - `fwhm_at(lam, coeffs)` and `fwhm_err_at(lam, cov)` (the propagated
    1-sigma error), both taking the same `lam_ref` / `lam_scale` keywords.
  - The `SRFModel` dataclass:
    - Fields: `channel, coeffs, cov, offset, offset_err, lam_min, lam_max,
      lam_ref, lam_scale, chi2_nu, npts, instrument, frame, meta`.
    - Methods: `fwhm()`, `fwhm_err()`, `sigma()`, `to_dict/from_dict`, and
      `to_json(path=None)` / `from_json(str_or_path)`.  NaN is written as
      JSON `null`, so the files are standard JSON.
  - Additions beyond the spec, for task 5:
    - `SRFModel.from_lines(channel, fits)`: fits FWHM(λ) to the `ok &
      use_for_srf` rows of a `fit_lines` table.  `offset` is the weighted mean
      of `dmu` over the unblended lines (Hα, Hβ, the Ca II IR triplet; Q&A
      Q8), with its error inflated by √χ²_ν when the lines disagree.
    - `save_srf_models(path, models, meta)` and `load_srf_models(path)`: one
      JSON file holding several channels, as
      `hypernet/data/veit_srf_model.json` will.
- `offset` is absolute (measured − laboratory, in `frame` = air by default).
  The relative L − E offset is the difference of two models' offsets, or is
  taken per line in task 5.
- Tests added to `hypernet/tests/test_srf.py`:
  - an exact quadratic recovered to 1e-10;
  - the covariance scaling as err² (err halved gives cov/4);
  - bad points dropped, and too few points giving NaN;
  - JSON round trips (string, file, several models, NaN ↔ null);
  - `from_lines` on a synthetic 13-line spectrum recovering the quadratic to
    0.02 nm and a +0.04 nm offset.
- `pytest -q`: **34 passed**.

### 2026-09-28 -- Build #4 (Opus 5.5)

- New module `hypernet/whn_l1a.py`, following `whn_explore`'s style (dict
  returns, numpy docstrings):
  - `wavelengths_root(path=None)`: explicit path, then
    `$OS_COLOR/WATERHYPERNET/Wavelengths`; raises `FileNotFoundError`.
  - `parse_name`, `find_products(root, recursive=True)` (a DataFrame
    index; recursive for the 0b delivery subfolders) and
    `sequence_files(site, seq_time)` (product → path; raises on ambiguity).
  - `load_l1a_irr(path)`: `wave`, `scans`, `mean`, `bandwidth`, the per-scan
    geometry (`sza, saa, vza, vaa, paa, time, series_id, quality_flag`) and
    `meta` (attributes: `system_id`, cal files and dates, ...).
  - `load_l1a_rad(path)`:
    - `Lu` (vza < 90, water) and `Ld` (vza ≥ 90, sky) sub-dicts, each with
      `scans`, `mean` and geometry.
    - `Lu['t_mean']`, and `Ld['at_lu_time']`, the sky interpolated in time
      as L1C does.
  - `load_l1c(path)` and `load_l2a(path)`: `irradiance`,
    `downwelling_radiance`, `upwelling_radiance`, `water_leaving_radiance`,
    `reflectance` and `reflectance_nosc`, plus `rhof`, `epsilon`, the
    geometry and, for L2A, the `std_*` variables and scan counts.  Arrays
    keep the file layout.
  - `series_means`, `interp_in_time`, and `check_against_l1c(irr, rad, l1c)`
    (maximum relative E / Lu / Ld differences against L1C).
- `hypernet/srf.py`:
  - new `flatten_scans()`, which divides out a smooth per-scan scale;
  - new `scan_errors(..., flatten_px=None)` option (Q&A Q9).
- New `wiggles/phase0a_l1c_consistency.py`, the diagnostics behind the Build
  #4 Q&A: per-scan metadata, how L1C is built, the zero L1C uncertainties,
  raw vs flattened scatter, and where L1C E departs from `np.interp`.
- New `hypernet/tests/test_whn_l1a.py`, 8 tests:
  - Tier 1: filename parsing, time interpolation, explicit root.
  - Tier 2 (`needs_wavelengths`), the VEIT sample:
    - 1536/1538 px; 6 IRR, 12 RAD (6 Lu + 6 Ld) and 6 L1C scans;
    - the vza split (L1C vza = the water view; the sky is > 5x brighter);
    - L1C E = `np.interp` of the L1A mean to 1e-4 over 400-900 nm (1e-3
      over 380-1000);
    - `check_against_l1c`: E < 1e-4, Lu exact, Ld < 1e-3;
    - L2A sanity.
  - The spec's "1e-4" holds only over 400-900 nm; see the Q&A.
- Also in `test_srf.py`: a test that flattening recovers 0.4 % pixel noise
  under 5 % broadband scan changes.
- `pytest -q`: **43 passed**.

### 2026-09-28 -- Build #5 (Opus 5.5)

- JXP accepted Q9: errors from `scan_errors(..., flatten_px=20)`, with the
  raw scatter reported alongside.
- New `wiggles/phase0a_veit.py`:
  - Loads the VEIT sequence through `whn_l1a` and fits `srf.LINES` in:
    - E: the mean of the 6 L1A_IRR scans;
    - Ld: the mean of the 6 sky scans (both series, **not** the
      time-interpolated L1C Ld, for S/N);
    - Lu: the mean of the 6 water scans.
  - Each channel is fitted twice, with flattened errors and with raw errors
    (`err_kind`).
  - `SRFModel.from_lines` per channel on the flattened fits, with
    `scale_cov=True`.
  - Writes:
    - `$OS_COLOR/hypernet/wiggles/phase0/veit_lines.parquet`;
    - `wiggles/phase0_veit_lines.csv` (78 rows = 13 lines × 3 channels × 2
      error kinds);
    - `hypernet/data/veit_srf_model.json` (E, Ld, Lu, with meta: method
      "empirical ... to be superseded by task 3b", errors, cal dates);
    - `wiggles/figs/phase0/veit_fwhm_vs_lambda.png`.
  - The figure has two panels: FWHM per line with the quadratic ±1σ band,
    and FWHM_L − FWHM_E per line with the model difference.  Palette slots
    1-3 plus distinct markers; filled = clean, open = blend/band, faded =
    diagnostic.
  - Stdout: the re-derived §2.3 table next to the plan's (relabelled), and
    the absolute offsets in air and vacuum.
- **Errors.**  The median relative error of the mean over 400-900 nm, raw
  vs flattened:
  - E: 2.6e-3 vs 1.7e-3;
  - Ld: 2.4e-2 vs 1.8e-3;
  - Lu: 2.1e-2 vs 3.1e-3.

  With raw errors the per-line χ²_ν is 0.03 (Ld) and 0.09 (Lu), which
  confirms Q9.  With flattened errors it is 3-7: the Gaussian-on-linear-
  continuum model is not adequate to the S/N, which argues for the template
  fit.
- **Reproducing plan §2.3** (after swapping its mislabelled Ld/Lu columns):
  - G band, Hβ, Mg b, Na D, Hα (E, Ld) and O2-B reproduce to ≤ 0.05 nm.
  - Hα Lu: 3.08 vs 2.96 (reweighting; Lu is the noisiest channel).
  - **Ca H/K change a lot with the joint fit.**
    - FWHM_E 2.67 / 2.77 (plan 1.83 / 1.91); Ld 3.08 / 3.02; Lu 2.99 / 3.35.
    - The E-L excess shrinks from 0.6-0.9 nm to 0.25-0.6 nm.
    - The Lu − E centroids, which had opposite signs (+0.20 / −0.14), become
      +0.08 / +0.01.  So the sign flip was a blend artefact.
    - Neither version is a good SRF measure, because the H/K wings extend
      past the window.
  - O2-A: Ld 3.74 vs 3.42, Lu 3.28 vs 3.46.  The band shape depends on the
    path and the weights, which is why it is diagnostic only.
- **What survives:** FWHM_L − FWHM_E ≈ 0.4-0.6 nm over 400-600 nm (G band,
  Hβ, Mg b, Na D, each at > 5σ), ≈ 0.2-0.3 nm at Hα, and ≈ 0 ± 0.2 nm at the
  Ca II IR triplet.
  - G0(a) E ≠ L therefore holds on this sequence, with the difference fading
    to the red.
  - The earlier "0.5-0.9 nm" for the blue edge rested on the Ca H/K single
    fits.
  - Ld and Lu agree to within 0.1-0.2 nm, but Lu is systematically wider at
    Ca H, the G band and Hα.  Task 6 tests this.
- **The FWHM(λ) quadratics (flattened errors, cov × χ²_ν):**
  - E: flat, 2.66 / 2.75 / 2.80 / 2.80 / 2.73 nm at 400 / 500 / 600 / 700 /
    850 nm, χ²_ν = 26.
  - Ld: 3.07 / 3.16 / 3.17 / 3.08 / 2.79, χ²_ν = 13.
  - Lu: 3.16 / 3.30 / 3.28 / 3.11 / 2.54, χ²_ν = 5.
  - The large χ²_ν says the per-line widths scatter about any smooth curve
    by ~0.2-0.3 nm, far above their 0.02-0.05 nm errors.  That is the
    line-dependent intrinsic broadening (blends, G band 2.49 vs Mg b 3.06 in
    E).
  - So the empirical quadratic is a fair relative description, not an
    absolute SRF.  **Recommend doing task 3b (HSRS template) before Phase 1
    consumes `veit_srf_model.json`.**
- **Absolute wavelength scale (G0(c)): inconclusive.**
  - The clean lines disagree among themselves by ±0.15 nm, i.e. 10× their
    0.015 nm errors:
    - Hα: −0.02 nm vs air, −0.20 vs vacuum;
    - Hβ: +0.20 vs air, +0.06 vs vacuum;
    - Ca II 849.8 / 854.2 / 866.2: +0.37 / +0.10 / +0.25 vs air.
  - The weighted means are +0.15 nm vs air and −0.05 to 0.00 nm vs vacuum in
    all three channels.  That leans towards vacuum but rests on Hβ and the
    Ca II IR lines, which the empirical fit handles worst.
  - Kevin's answer on air/vacuum, and the template fit, are both needed.
  - The relative L − E centroids are small and consistent: +0.02 to +0.1 nm
    at the clean lines.
- **Plan §2.3 corrected** (Q7 default, approved): the headers now read
  `FWHM Lu | FWHM Ld | centroid Lu − E`, followed by a dated correction note
  (the swap, the joint Ca H/K widths, and the rows that differ).  The table
  values are unchanged.
- `pytest -q`: 43 passed (no package code changed in this step).

### 2026-09-28 -- Build #6 (Opus 5.5)

- No new JXP answers since Q9.
- New `wiggles/phase0a_scans.py`.  It reuses `phase0a_veit.load()` and
  `fit_all()` and adds:
  - every line fitted in every single scan (6 E, 6 Ld, 6 Lu), weighted by
    the flattened single-scan error;
  - the two Ld sky series fitted separately (3 scans each);
  - per line and channel, the scan std of centroid, FWHM and depth over the
    median `curve_fit` error, plus the error of the mean from the scan
    scatter;
  - the Ld − Lu test on the mean fits, with z taken over the larger of the
    fit and scan errors;
  - depths and weak-line equivalent widths (EW = depth × σ × √(2π)) as a
    diagnostic.
- Outputs:
  - `$OS_COLOR/hypernet/wiggles/phase0/veit_scan_fits.parquet` (all
    per-scan fits);
  - `wiggles/phase0_veit_scans.csv` (a 39-row summary);
  - `wiggles/figs/phase0/veit_scan_scatter.png`, in four panels: the FWHM
    and centroid ratios (log axis), Ld − Lu FWHM, and L/E depth and EW.
- **The errors are right (Q9 vindicated).**
  - Scan std / median fit error is 0.4-1.9, median ≈ 0.9, for FWHM and
    centroid in all three channels.
  - The mean-spectrum fit error and the scan-scatter error of the mean agree
    to within ~×1.5 (e.g. Hα E FWHM 0.041 vs 0.025 nm, Ld Hβ 0.051 vs
    0.066).
  - The outlier is Lu Ca II 849.8, with ratios of 2.3 for FWHM and 4.0 for
    the centroid: the faint red water view, in the joint 849.8/854.2 fit.
- **The Gaussian is adequate at single-scan S/N for isolated lines only.**
  Median per-scan χ²_ν:
  - 1.1-1.5 at Hβ, Mg b, Na D and Hα;
  - 0.9-1.8 at the Ca II IR lines;
  - 1.6-5 at Ca H/K;
  - 16 at the G band in E;
  - 4 at O2-B, and 60-170 at O2-A in E and Ld.

  The misfits are exactly the blends and bands.  This is the strongest
  argument yet for the template fit (task 3b).
- **The widths are stable over the sequence.**  Ld series A (before the
  water view) vs B (after), ~50 s apart with the sky ~5 % brighter or
  dimmer: every |z| < 1.7 in FWHM and < 1.8 in centroid (G band −0.10 ±
  0.07, Hβ +0.17 ± 0.10 nm).
- **Ld = Lu.**
  - FWHM over the 10 SRF lines: χ² = 26.5 for 10 dof.  Two lines have
    |z| > 2: Ca H (−0.33 ± 0.09 nm, 3.9σ) and the G band (−0.16 ± 0.07,
    2.2σ).
  - Weighted mean Ld − Lu = −0.06 ± 0.03 nm, i.e. Lu slightly wider.
  - Centroids: χ² = 13.1 for 10 dof, no |z| > 2; mean Ld − Lu = +0.020 ±
    0.010 nm.
  - O2-A: +0.46 ± 0.04 nm (11σ).  That is expected for a band whose shape
    depends on path: Lu adds the water-leaving and surface-reflected paths.
  - Reading: one L SRF to ~0.1 nm at the clean lines.  The blue blend
    excess in Lu is plausibly Lw content (water Raman filling of Ca H and
    the G band in clear water) rather than optics.  **The sky (Ld) is the
    cleaner L-SRF probe, and Lu the check.**
- **Depths and EW (diagnostic, Ring).**
  - Depths in Ld and Lu are 10-17 % shallower than in E from 390 to 660 nm,
    but EW_Ld / EW_E = 0.91-1.09 there.
  - So the shallower lines are mostly the wider L SRF (EW is conserved),
    not Ring filling-in.  Ld's lowest EW ratios, Ca H 0.91 and Hα 0.94,
    leave room for a few-per-cent filling.
  - In Lu, EW_Lu / EW_E falls to 0.64-0.85 at O2-B and the Ca II IR lines,
    where the water view is faintest, while Ld stays at 1.0-1.1.  That
    suggests an **additive, line-free component in Lu in the red** (stray
    light or a dark/offset residual), which dilutes relative depths.  Worth
    raising with Kevin, and a reason not to use Lu red lines for the SRF.
- `pytest -q`: 43 passed (no package code changed).

### 2026-09-28 -- Build #7 (Opus 5.5)

- No new JXP answers since Q9.
- New `wiggles/phase0a_budget_veit.py`.  `budget(irr, rad, l2a, fits_E,
  fits_L)` is the per-sequence function that task 12 will loop over; `main()`
  runs it on VEIT.
  - **Windows.**  One per group (Ca H/K and Ca II 849.8/854.2 jointly),
    ±5 nm, with a linear continuum through the pixels 5-7 nm out.
    Residual = y/continuum − 1.
  - **Measured:** Ld/Ed_L and Lu/Ed_L, with Ld interpolated in time to the
    water view and Ed_L = `np.interp` of E, i.e. the processor's own
    quantities; plus ρw from L2A.
  - **H1:** Ed_cubic/Ed_linear − 1.
  - **H2:** the E line re-observed at the Ld width with its equivalent width
    conserved, over the E line.  It uses the task 5 E depth/σ/μ and Ld σ;
    the measured L depths are not used, so the Ring effect is kept out.
  - **Metrics:**
    - the rms of each profile and of the noise within ±5 nm;
    - ρw in absolute units against ρw × H2 and ρw × H1, with noise
      `std_reflectance/√n`;
    - weighted-LS projections: onto H2 alone (α₂), H1 alone (α₁), and jointly
      onto (H1, H2, shift = d ln E_line/dλ, in nm);
    - rms ρw''.
  - Outputs: `wiggles/phase0_veit_budget.csv` (11 rows) and
    `wiggles/figs/phase0/veit_budget.png`.  The figure shows four example
    profiles, the rms per line for the ratios, the rms for ρw (absolute), and
    α₂ per line.
- **H2 accounts for the line structure; H1 doesn't.**
  - rms over ±5 nm, measured Ld/Ed vs H2 vs H1: Ca H/K 2.7e-2 / 2.4e-2 /
    2.8e-3; G band 1.7e-2 / 1.6e-2 / 1.3e-3; Hβ 7.7e-3 / 8.0e-3 / 5e-5; Mg b
    5.8e-3 / 5.7e-3 / 1.8e-4; Na D 5.8e-3 / 5.1e-3 / 5.7e-4.
  - At Hα H2 gives half the measured value (5.0e-3 vs 2.7e-3).  At the Ca II
    IR and O2-B the measured residual is close to the noise.
  - The median rms_H2/rms_H1 over the SRF lines is 9 (range 1.4-160).
  - α₂ (H2 alone), weighted over the 10 SRF lines: **0.93 ± 0.06** in
    Ld/Ed, **1.01 ± 0.07** in Lu/Ed and **1.08 ± 0.14** in ρw.  Per line,
    α₂ is 0.8-1.5 from Ca H/K to Hα, with Na D in Lu and ρw low (0.3 ± 0.5).
  - The H1 profile has the **opposite sign** to the measured feature.  Linear
    interpolation overestimates Ed at a line core, so the ratio dips there;
    the measured ratio peaks, as H2 predicts.  So α₁ alone is large and
    negative (−3.7 ± 0.3 in Ld/Ed), and the H1 and H2 bases are
    anti-correlated (r = −0.6 to −0.85).
  - The joint (H1, H2, shift) fit is therefore degenerate: α₂ drops to 0.5
    and α₁ goes to −3 to −6, which is unphysical.  It is reported but not
    used.
  - Its shift term is small, −0.04 ± 0.01 nm (Ld/Ed): a residual L − E
    centroid offset of the size task 5 found.
- **ρw.**
  - The measured ρw residual rms matches ρw × H2 from Ca H/K to Na D
    (e.g. G band 2.4e-4 vs 1.8e-4, Hβ 1.5e-4 vs 1.2e-4).  ρw × H1 is 10-100×
    smaller.
  - In the red (O2, Ca II IR, H₂O) the ρw residual stays at 0.3-3e-4 while
    both predictions fall below 1e-5.  That structure is not Ed-borne; there
    ρw → 0, so it comes from the Lu − ρ·Ld subtraction (Lu's red additive
    component, task 6; O2 path differences).
  - The L2A `std_reflectance/√n` noise is ≥ the measured residual at several
    lines.  Like the raw L1A scatter (Q9), it includes broadband scan-to-scan
    changes and overstates the per-pixel noise.
- **Caveat: the test is partly circular.**  H2 is built from Gaussian fits
  to the same E and Ld spectra whose ratio it predicts.  α₂ ≈ 1 does not
  independently prove the SRF difference.  What it does establish:
  - (i) the ratio residual is what the measured width difference produces
    **with EW conserved**, so Ring filling-in is not needed to explain it;
  - (ii) the ρw wiggle at the Fraunhofer lines has the H2 shape and full
    amplitude;
  - (iii) linear interpolation has the wrong sign and is 10× too small.

  The independent tests come from Phase 1 (twin experiment) and Phase 3
  (correcting the SRF and seeing whether ρw'' drops).
- **The task's closed form**, depth × (FWHM_L² − FWHM_E²)/FWHM_L², is
  1.5-2× the H2 central value.  To first order the central value is
  a(1 − σ_E/σ_L) ≈ a(σ_L² − σ_E²)/(2σ_L²).  Both are in the CSV
  (`h2_closed_form`, `h2_central`).  Task 12 should use the profile.
- `pytest -q`: 43 passed.

### 2026-09-28 -- Build #8 (Opus 5.5)

- No new JXP answers since Q9.
- New `wiggles/phase0a_interim.py`: the report's numbers (model FWHM with
  errors, E − L differences and significance, per-line and weighted L − E
  centroids, χ²_ν by line class), all read from the committed products.
- Fixed a precision bug: the committed CSVs were written with `%.5g`, which
  rounds centroids near 400-900 nm to 0.01 nm.  The writers now use `%.7g`
  (`%.6g` for the budget), and I re-ran `phase0a_veit.py`, `phase0a_scans.py`
  and `phase0a_budget_veit.py`.  The JSON and figures were unaffected (they
  are computed from the full-precision fits).
- Wrote the 0a interim report to Q&A (Build #8), covering four points:
  - the FWHM(λ) coefficients and their caveats;
  - E ≠ L survives at 5-15σ per line from Ca K to Na D, and fades to
    < 1.5σ in the Ca II IR, so G0(a) is not met and H2 stays;
  - the relative L − E centroid is +0.050 ± 0.011 nm, below 0.1 nm; the
    absolute offset is open (air/vacuum);
  - the Gaussian is adequate for clean lines, and it is the line model, not
    the SRF shape, that fails at blends and bands.
- Added Q10 (task 3b before Phase 1) and Q11 (Ld as the reference L
  channel), both with defaults.
- `pytest -q`: 43 passed.

### 2026-09-29 -- Build #3b and the 0b scripts (Opus 5.5)

- JXP accepted the Q10 and Q11 defaults: task 3b first, with the template
  models in `veit_srf_model.json` and the empirical ones in
  `veit_srf_model_empirical.json`; Ld is the reference L channel.  JXP has
  sent the data request to Kevin, and
  `correspondence/wiggles_note_to_kevin.md` is now marked "sent".
- **HSRS.**
  - URL: `https://lasp.colorado.edu/lisird/resources/lasp/hsrs/v2/hybrid_reference_spectrum_p005nm_resolution_c2022-11-30_with_unc.nc`
    (TSIS-1 HSRS v2, created 2022-11-30).
  - SHA-256: `dd9f62fb9b39433631013ebf052429f4daddb2bd7e0d970a6d292be6026f3e20`.
  - The file is 60.7 MB and is stored at `$OS_COLOR/hypernet/ref/`.
  - It covers 202-2730 nm at 0.005 nm resolution, sampled every 0.001 nm,
    in W m⁻² nm⁻¹, on **vacuum** wavelengths (variable `Vacuum
    Wavelength`).
  - The p01nm/p025nm variants at the same path 302 or 200 as well; only
    p005nm was fetched.
- **New `hypernet/refspec.py`:** `fetch_hsrs(force, verify)` with a
  checksum check, `hsrs_available()`, and `load_hsrs(wmin, wmax, step=5,
  frame='vac'|'air')`, which block-averages to 0.005 nm by default.
- **`hypernet/srf.py` additions:**
  - `vac_to_air()`, the fixed-point inverse of `air_to_vac`;
  - `convolve_gaussian()`, a direct weighted mean over ±6σ;
  - `fit_srf_template(wav, spec, ref_wave, ref_flux, lo, hi, err, veil)`,
    with model (c₀ + c₁x)·(R_σ(λ − dlam) + veil·⟨R_σ⟩).  Here `dlam` is
    measured − reference (the sign of `dmu`).  Bounds are σ 0.15-3.5 nm,
    |dlam| < 1 nm and veil −0.3 to 0.9; a parameter at a bound means
    failure, and the function never raises;
  - `template_windows()` and `TELLURIC`, the exclusion bands;
  - `fit_template_windows()`, which returns a table shaped for
    `SRFModel.from_lines`.
- **Telluric exclusions** are 570-580 (O2-O2 at 577), 626-634, 640-650,
  685-750, 757-773, 780-845 and > 870 nm.  They were widened after the first
  VEIT run, whose windows at 575, 645, 745, 795 and 875 nm gave offsets of
  −0.5 to −1 nm and one FWHM of 5.6 nm.
- **New `hypernet/whn_srf.py`:**
  - `fit_sequence(irr, rad, ref, template, empirical, veil, per_scan)`
    returns the empirical and template tables and six `SRFModel`s;
  - `fit_files()`;
  - `models_table()`, one row per model: coefficients, errors, the
    flattened covariance, and FWHM at 400/450/500/600/700/850 nm;
  - `load_reference()`.
- **New `wiggles/phase0a_template_veit.py` (task 3b):**
  - Fits the 10 nm template windows with the veil free and fixed.
  - Writes `veit_template.parquet`, `wiggles/phase0_veit_template.csv`,
    `hypernet/data/veit_srf_model.json` (template, veil 0 for all channels)
    and `wiggles/figs/phase0/veit_template_fwhm.png`.
  - The veil comes out at −0.03 ± 0.05 and freeing it makes about a third of
    the windows fail, so it is shipped fixed at 0.
  - Results are in the Build #3b Q&A.  `phase0a_veit.py` now writes
    `veit_srf_model_empirical.json`, and `phase0a_interim.py` reads it.
- **The 0b scripts, all run on the VEIT sequence:**
  - `wiggles/phase0b_index.py` (task 9) uses the new
    `whn_l1a.sequence_table()` and `match_request()`.
    - It matches by site code, `sequence_time` and azimuth, and reads the
      L1A attributes for the instrument, cal-date and cal-attribute-bug
      checks.
    - It counts primaries and spares per request row, and substitutions of
      spares for missing primaries.
    - It writes `request_index.parquet` and
      `wiggles/phase0_delivery_summary.csv`.
    - On VEIT: 224 missing, 1 extra (the sample is not in the request), and
      the attribute bug is detected.
  - `wiggles/phase0b_fit_all.py` (task 10):
    - loops over `whn_srf.fit_files`, with template veil 0 and flags
      `--per-scan`, `--no-template` and `--limit`;
    - keeps going past failures, recorded in `fit_status.csv`;
    - writes `line_fits`, `template_fits`, `srf_models` and `scan_fits`
      `.parquet`;
    - takes 1 s per sequence and reproduces the VEIT numbers.
  - `wiggles/phase0b_stability.py` (task 11):
    - per instrument × cal period × channel, the mean, std and range of FWHM
      at 400-700 nm against the 0.2 nm criterion, plus the Ld − E width and
      offset;
    - slopes against SZA, season (cos of month) and year, and means by sky
      class;
    - the plan's pairs, with 122304 vs 122305 at VEIT split by 122304's
      calibration period, since 122305 ran between them;
    - writes three CSVs and `stability_vs_{sza,time}.png`.
    - A synthetic 30-sequence smoke test (scratch, not committed) recovered
      an injected 0.3 nm recalibration jump as 0.28 nm (12σ) and a null
      site pair as z = −0.5.
    - On VEIT alone: one group, no pairs.
  - `wiggles/phase0b_budget.py` (task 12):
    - loops `phase0a_budget_veit.budget()` over sequences with L1A + L2A,
      using the E and Ld empirical fits from `line_fits.parquet`;
    - adds a ρw'' excess, over the 25th percentile of rms ρw'' in 10 nm
      tiles over 400-700 nm (the median tile is not line-free);
    - summarises by instrument and water type;
    - writes `budget_lines.parquet`, `wiggles/phase0_budget{,_summary}.csv`
      and `budget_all.png`.
    - On VEIT it reproduces task 7 (α₂ = 0.93 in Ld/Ed, 1.08 in ρw).  The
      ρw'' excess is 8.3× at Ca H/K, 5.1× at the G band, 2.6-3.2× at Hβ and
      Mg b, and 0.6-1.2× in the red.
  - `wiggles/phase0b_calfiles.py` (task 13):
    - finds `HYPERNETS_CAL_HYPSTAR_*_{RAD,IRR}_v*.nc` and any
      line-spread-like files;
    - lists every variable and attribute and flags the SRF-like ones;
    - compares the cal wavelengths with the L1A grids, and
      `bandwidth`/`fwhm` with the fitted template FWHM.
    - No cal files are present yet, so it says so and exits.  A smoke test
      on an L1A file as stand-in found `bandwidth` and `wavelength` and
      compared 3.0 nm against the fitted E 2.2-2.6 nm.
- **Tests.**
  - `test_srf.py` gained 8: vac↔air inversion, window tiling, template
    recovery of σ / dlam / veil (parametrised), graceful failure, template →
    `SRFModel`, and HSRS recovery (Tier 2).
  - `test_whn_l1a.py` gained a synthetic test of `sequence_table` and
    `match_request`.
  - New `test_whn_srf.py`, Tier 2: the HSRS checksum, and the Hα minimum at
    656.46 nm (vacuum) and 656.28 nm (air); plus `fit_sequence` on VEIT,
    which stands in for task 10's delivered-sequence test.
  - `pytest -q`: **53 passed**.
- Updated `wiggles_phase1_prompts.md` Context: HSRS already fetched via
  `hypernet/refspec.py`; the JSON now holds template SRFs with its valid
  range; the module names; Ld as the L reference.
- A slip to report: I ran `git mv -n` (a dry run) while checking a rename.
  It changed nothing (`git status` confirmed it), but it should not have
  been run.
- Still waiting on the delivery: running tasks 9-13 for real, and the gate
  (task 14).

### 2026-09-29 -- Moved `wiggles/` into the package (Opus 5.5)

- At JXP's request, the top-level `wiggles/` directory is now the subpackage
  `hypernet/wiggles/`, so it installs with the repository (`find_packages()`
  picks it up).  I moved it with plain `mv`; JXP runs git, and `git add -A`
  should record renames.
- New `hypernet/wiggles/__init__.py` holds the path constants: `WIGGLES_DIR`
  (committed CSVs), `FIGDIR` (`figs/phase0/`), `DATA_DIR` (`hypernet/data/`),
  `REPO` and `OUT` (`$OS_COLOR/hypernet/wiggles/phase0`, falling back to
  `.` when `$OS_COLOR` is unset).
- The scripts:
  - The `sys.path` hacks are gone.  Sibling imports are now
    `from hypernet.wiggles import phase0a_veit as p5`, and so on.
  - They run as modules: `python -m hypernet.wiggles.<script>`.
  - `phase0a_l1c_consistency` and `phase0a_interim` were flat scripts; they
    are now wrapped in `main()`, so nothing runs at import.
  - `phase0a_l1c_consistency` resolves its data directory through
    `whn_l1a.wavelengths_root()`.
- Paths updated:
  - the docstrings in `srf.py`, `whn_l1a.py`, `whn_srf.py` and
    `test_whn_l1a.py`;
  - one path in the `docs/wiggles_planning.md` §2.3 correction note;
  - everything above `## Q&A` in this doc, including a new Paths
    convention.
  - The Q&A and Logs keep their original `wiggles/...` paths as a
    historical record; read them as `hypernet/wiggles/...`.
  - `CLAUDE.md` layout now lists `srf.py`, `whn_l1a.py`, `whn_srf.py`,
    `refspec.py`, `data/` and `wiggles/`.
- Verification:
  - I re-ran all 11 scripts through `python -m`.  The three CSVs tracked at
    HEAD (`phase0_veit_{lines,scans,budget}.csv`) are byte-identical.
  - The JSON `script` provenance now reads `hypernet/wiggles/...`.
  - No top-level `wiggles/` is recreated.
- New tests in `test_import.py`: every `hypernet.wiggles` module imports and
  has `main()`, including in a subprocess with `$OS_COLOR` unset.
  `pytest -q`: **55 passed**.
- Note: the committed CSVs and PNGs sit inside the package directory, but
  `package_data` ships only `data/*.json`, so a non-editable install carries
  the code, not the figures and tables.

### 2026-09-29 -- Build #8b: Phase 0a report slides (Opus 5.5)

- Wrote `docs/slides/wiggles_phase0_report.pptx`, 14 slides on a 16:9 wide
  canvas (13.33 × 7.5 in).  **No text below 20 pt**, tables included; I
  checked the slide XML, whose minimum `sz` is 2000.  Order:
  1. title;
  2. H1 vs H2;
  3. what 0a built;
  4. four things the VEIT files taught us;
  5. the task 5 figure, then 6. its table;
  7. the task 6 figure, then 8. its table;
  9. the task 7 figure, then 10. the α₂ callouts, then 11. its table;
  12. the task 3b template SRF, figure and table;
  13. the gate G0 status;
  14. next steps.

  Every slide has speaker notes.
- The numbers are not hand-typed.  `docs/slides/wiggles_phase0_report_data.py`
  reads the committed products (`hypernet/wiggles/phase0_veit_{lines,scans,
  budget}.csv`, `hypernet/data/veit_srf_model.json`) and the figure aspect
  ratios, and writes `wiggles_phase0_report_data.json`.  Then
  `docs/slides/build_wiggles_phase0_report.js` (pptxgenjs) builds the deck.
  Rebuild with:

      python docs/slides/wiggles_phase0_report_data.py
      NODE_PATH=$(conda run -n slides npm root -g) conda run -n slides node docs/slides/build_wiggles_phase0_report.js

- Tooling: this Mac had no node, LibreOffice, PowerPoint or Keynote.
  - I created a separate conda env `slides` (nodejs, pptxgenjs installed
    with `npm -g`, python 3.12, lxml, defusedxml, pillow), leaving `ocean14`
    untouched.
  - With JXP's OK I installed LibreOffice with `brew install --cask
    libreoffice`, for the render QA.
- QA:
  - The pptx skill's `validate.py` passes.
  - I rendered every slide through LibreOffice and inspected it.  Fixes from
    the first render: the task 3b title wrapped onto the figure (shortened);
    the table footnotes floated far below the tables (fixed row height,
    notes placed under the table); the gate cards were too tall; and the
    "in time" callout became "interpolated".
  - The embedded figure PNGs keep their own axis labels, which are smaller
    than 20 pt on the slide.  They are images, not slide text.
- Design: navy/teal with a coral accent for the key numbers; a circle badge
  motif (H1/H2, steps 1-4, gate a-c); Cambria headings and Calibri body.
- `docs/slides/` holds no `.md`, so nothing there is published by Sphinx.

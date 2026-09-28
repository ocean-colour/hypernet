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
  (`python -m hypernet.whn_explore 1`, `python wiggles/phase0a_veit.py`).
- Run `pytest -q` after each step where relevant.  34 tests today; the archive-
  dependent ones skip themselves when `$OS_COLOR` is not mounted.
- **Tier 2** -- steps that read the VEIT sample or the archive need `$OS_COLOR`
  mounted.  Do **not** unset `$OS_COLOR`.
- Any calculation goes into a script on disk, not into the chat.
- Data intermediates (parquet/npz) go **outside** the repo under
  `$OS_COLOR/hypernet/wiggles/phase0/`; only figures and small tables are
  committed.
- **Paths** (Q&A Setup #1, Q5): new Phase 0 scripts, figures and small
  tables go in the top-level `wiggles/` directory (`wiggles/phase0a_*.py`,
  `wiggles/figs/phase0/`, `wiggles/phase0_*.csv`).  The earlier exploratory
  scripts stay in `wavecal/`.  The SRF JSON goes in `hypernet/data/`.
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

5. **VEIT fits, with uncertainties.**  `wiggles/phase0a_veit.py`: fit E from
   the mean L1A_IRR and L from the mean Ld and Lu separately, then fit
   FWHM(λ) per channel.  Write the line table to
   `$OS_COLOR/hypernet/wiggles/phase0/veit_lines.parquet`, a committed copy to
   `wiggles/phase0_veit_lines.csv`, the three `SRFModel`s to
   `hypernet/data/veit_srf_model.json` (small, committed -- Phase 1 reads it),
   and a figure `wiggles/figs/phase0/veit_fwhm_vs_lambda.png` (FWHM points
   with error bars, quadratic fits, E vs Ld vs Lu).  Reproduce the plan §2.3
   table to the quoted precision and note any differences in the Logs.  Log
   your work.

6. **Ld = Lu check and per-scan scatter.**  Extend the script (or add
   `wiggles/phase0a_scans.py`): fit every scan separately (6 E scans, the Ld
   and Lu scans), report the scan-to-scan spread of centroid and FWHM per line
   against the `curve_fit` errors, and test Ld = Lu per line (difference over
   combined error).  Line depths in E, Ld and Lu are tabulated as a
   diagnostic only (Ring effect).  Figure
   `wiggles/figs/phase0/veit_scan_scatter.png`.  Log your work.

7. **Preliminary H1/H2 budget.**  `wiggles/phase0a_budget_veit.py`: per line,
   (H1) the residual line amplitude in Ld/Ed produced by linear interpolation
   alone (linear minus cubic on the same Ed, as in `sanity_checks_veit.py`) and
   (H2) the residual predicted from the FWHM difference,
   depth × (FWHM_L² − FWHM_E²)/FWHM², against the measured residual line
   amplitude in Ld/Ed and the ρw'' power within ±5 nm of the line.  Write
   `wiggles/phase0_veit_budget.csv` and
   `wiggles/figs/phase0/veit_budget.png`.  This script is the template for
   task 12.  Log your work.

8. **0a interim report.**  In Q&A, summarise: the FWHM(λ) coefficients and
   uncertainties handed to Phase 1, whether E ≠ L survives the uncertainties,
   the centroid offsets against the 0.1 nm threshold of G0(c), and anything
   that argues against the Gaussian.  Log your work.

### Build (0b, full request -- after the data arrive)

9. **Ingest and index.**  Mirror the delivery with `rclone` to a subfolder of
   `$OS_COLOR/WATERHYPERNET/Wavelengths/` (name per the delivery).
   `wiggles/phase0b_index.py`: match every L1A/L1C/L2A file to a row of
   `docs/wiggles_data_request.csv` by site, `sequence_time` and azimuth,
   carry `instrument`, `sza`, `sky` and the cal period, record which
   requested sequences are missing and which spares were substituted, and
   check the `instrument_calibration_file_rad` attribute bug.  Write
   `$OS_COLOR/hypernet/wiggles/phase0/request_index.parquet` and a committed
   summary `wiggles/phase0_delivery_summary.csv`.  Log your work.

10. **Batch line fits.**  `wiggles/phase0b_fit_all.py`: `fit_lines` on the
    mean E, Ld and Lu of every sequence (per scan optional, behind a flag) and
    `fit_fwhm_model` per sequence and channel.  Outputs
    `$OS_COLOR/hypernet/wiggles/phase0/line_fits.parquet` and
    `srf_models.parquet` keyed by site, instrument, cal period, sequence.
    Add a Tier-2 test that one delivered sequence fits without error.  Log
    your work.

11. **Stability.**  `wiggles/phase0b_stability.py`: FWHM(λ) at 400, 500, 600,
    700 nm and the centroid offset (L − E) versus SZA, `sky`, month and time,
    per instrument; the spread per instrument against the 0.2 nm criterion;
    122304 and 121222 before/after recalibration; 122302 at BEFR vs THFR
    (site vs instrument); the two-instruments-one-site pairs at VEIT, GAIT and
    MAFR.  Figures under `wiggles/figs/phase0/` and a committed table
    `wiggles/phase0_stability.csv`.  Log your work.

12. **H1/H2 error budget per sequence.**  Generalise task 7 into
    `wiggles/phase0b_budget.py` over the whole index: per sequence and line,
    the H1 and H2 contributions versus the measured ρw'' excess at the line;
    summarise by instrument and water type.  Table and figure as above.  Log
    your work.

13. **Calibration files and lab data.**  `wiggles/phase0b_calfiles.py`: open
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

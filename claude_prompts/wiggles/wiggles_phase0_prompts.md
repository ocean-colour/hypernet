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
- Run `pytest -q` after each step where relevant.  43 tests today; the archive-
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

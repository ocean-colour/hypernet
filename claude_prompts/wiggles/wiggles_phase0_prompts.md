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
  (`python -m hypernet.whn_explore 1`, `python wavecal/phase0a_veit.py`).
- Run `pytest -q` after each step where relevant.  21 tests today; the archive-
  dependent ones skip themselves when `$OS_COLOR` is not mounted.
- **Tier 2** -- steps that read the VEIT sample or the archive need `$OS_COLOR`
  mounted.  Do **not** unset `$OS_COLOR`.
- Any calculation goes into a script on disk, not into the chat.
- Data intermediates (parquet/npz) go **outside** the repo under
  `$OS_COLOR/hypernet/wiggles/phase0/`; only figures and small tables are
  committed.
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

4. **L1A/L1C readers.**  New module (default `hypernet/whn_l1a.py`):
   `wavelengths_root()` resolving `$OS_COLOR/WATERHYPERNET/Wavelengths`;
   `load_l1a_irr(path)` (grid, per-scan and mean Ed); `load_l1a_rad(path)`
   splitting Ld (`viewing_zenith_angle < 90`) and Lu (≥ 90), per scan and
   mean; `load_l1c(path)` and `load_l2a(path)` returning `irradiance`,
   `downwelling_radiance`, `upwelling_radiance`, `reflectance`,
   `reflectance_nosc` and the geometry.  Tier-2 tests with a `needs_wavelengths`
   skip: 1536/1538 pixels, L1C `irradiance` equals `np.interp` of the L1A mean
   to 1e-4, six scans in each file.  Log your work.

5. **VEIT fits, with uncertainties.**  `wavecal/phase0a_veit.py`: fit E from
   the mean L1A_IRR and L from the mean Ld and Lu separately, then fit
   FWHM(λ) per channel.  Write the line table to
   `$OS_COLOR/hypernet/wiggles/phase0/veit_lines.parquet`, a committed copy to
   `wavecal/phase0_veit_lines.csv`, the three `SRFModel`s to
   `hypernet/data/veit_srf_model.json` (small, committed -- Phase 1 reads it),
   and a figure `wavecal/figs/phase0/veit_fwhm_vs_lambda.png` (FWHM points
   with error bars, quadratic fits, E vs Ld vs Lu).  Reproduce the plan §2.3
   table to the quoted precision and note any differences in the Logs.  Log
   your work.

6. **Ld = Lu check and per-scan scatter.**  Extend the script (or add
   `wavecal/phase0a_scans.py`): fit every scan separately (6 E scans, the Ld
   and Lu scans), report the scan-to-scan spread of centroid and FWHM per line
   against the `curve_fit` errors, and test Ld = Lu per line (difference over
   combined error).  Line depths in E, Ld and Lu are tabulated as a
   diagnostic only (Ring effect).  Figure
   `wavecal/figs/phase0/veit_scan_scatter.png`.  Log your work.

7. **Preliminary H1/H2 budget.**  `wavecal/phase0a_budget_veit.py`: per line,
   (H1) the residual line amplitude in Ld/Ed produced by linear interpolation
   alone (linear minus cubic on the same Ed, as in `sanity_checks_veit.py`) and
   (H2) the residual predicted from the FWHM difference,
   depth × (FWHM_L² − FWHM_E²)/FWHM², against the measured residual line
   amplitude in Ld/Ed and the ρw'' power within ±5 nm of the line.  Write
   `wavecal/phase0_veit_budget.csv` and
   `wavecal/figs/phase0/veit_budget.png`.  This script is the template for
   task 12.  Log your work.

8. **0a interim report.**  In Q&A, summarise: the FWHM(λ) coefficients and
   uncertainties handed to Phase 1, whether E ≠ L survives the uncertainties,
   the centroid offsets against the 0.1 nm threshold of G0(c), and anything
   that argues against the Gaussian.  Log your work.

### Build (0b, full request -- after the data arrive)

9. **Ingest and index.**  Mirror the delivery with `rclone` to a subfolder of
   `$OS_COLOR/WATERHYPERNET/Wavelengths/` (name per the delivery).
   `wavecal/phase0b_index.py`: match every L1A/L1C/L2A file to a row of
   `docs/wiggles_data_request.csv` by site, `sequence_time` and azimuth,
   carry `instrument`, `sza`, `sky` and the cal period, record which
   requested sequences are missing and which spares were substituted, and
   check the `instrument_calibration_file_rad` attribute bug.  Write
   `$OS_COLOR/hypernet/wiggles/phase0/request_index.parquet` and a committed
   summary `wavecal/phase0_delivery_summary.csv`.  Log your work.

10. **Batch line fits.**  `wavecal/phase0b_fit_all.py`: `fit_lines` on the
    mean E, Ld and Lu of every sequence (per scan optional, behind a flag) and
    `fit_fwhm_model` per sequence and channel.  Outputs
    `$OS_COLOR/hypernet/wiggles/phase0/line_fits.parquet` and
    `srf_models.parquet` keyed by site, instrument, cal period, sequence.
    Add a Tier-2 test that one delivered sequence fits without error.  Log
    your work.

11. **Stability.**  `wavecal/phase0b_stability.py`: FWHM(λ) at 400, 500, 600,
    700 nm and the centroid offset (L − E) versus SZA, `sky`, month and time,
    per instrument; the spread per instrument against the 0.2 nm criterion;
    122304 and 121222 before/after recalibration; 122302 at BEFR vs THFR
    (site vs instrument); the two-instruments-one-site pairs at VEIT, GAIT and
    MAFR.  Figures under `wavecal/figs/phase0/` and a committed table
    `wavecal/phase0_stability.csv`.  Log your work.

12. **H1/H2 error budget per sequence.**  Generalise task 7 into
    `wavecal/phase0b_budget.py` over the whole index: per sequence and line,
    the H1 and H2 contributions versus the measured ρw'' excess at the line;
    summarise by instrument and water type.  Table and figure as above.  Log
    your work.

13. **Calibration files and lab data.**  `wavecal/phase0b_calfiles.py`: open
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

## Logs

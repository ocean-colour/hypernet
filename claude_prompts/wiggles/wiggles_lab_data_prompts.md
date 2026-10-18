# WATERHYPERNET -- Wiggles: the lab wavelength-calibration data

## Goal

Fold the HYPSTAR laboratory wavelength-calibration data, sent by Joel on
2026-10-17, into the wiggles work while the L1A delivery is still
outstanding:

1. a reader and committed tables (the lab half of Phase 0 task 13);
2. a revised per-instrument SRF table, with each instrument's own E grid,
   122303 reclassified and an E-SRF shape column;
3. a test of whether the lab's common polynomial error explains the 0.13 nm
   common window pattern;
4. a validation of recovering native-grid E from the L2 irradiance (no
   Phase 3 run);
5. a provisional gate G0 report;
6. an update of the interim report.

Scope and decisions: `wiggles_prompts.md`, Q&A "New data #1"; JXP's answers
are recorded there.  **Phase 3 waits for the L1A delivery** (JXP,
2026-10-18).  Do not run the in-situ test on reconstructed E.

## Conventions

- `ocean14`.  Run via `conda run -n ocean14 python ...`; `conda activate` fails
  non-interactively.
- **Model.** Use Opus 5.5 (`claude-opus-5-5`) for this work, including any
  subagents (pass `model: opus`).
- **JXP runs git.**  Claude does not run any state-changing git command.
- Run scripts **from the repository root** so `hypernet` imports resolve.
- Run `pytest -q` after each step where relevant.  **137 passed on
  2026-10-18.**  The archive-dependent tests skip themselves when `$OS_COLOR`
  is not mounted, the OSOAA ones when the exe is missing.
- **Tier 2** -- most tasks need `$OS_COLOR`.  Do **not** unset it.
- Any calculation goes into a script on disk, not into the chat.
- Data intermediates go **outside** the repo under
  `$OS_COLOR/hypernet/wiggles/lab_data/`.  Only figures and small tables are
  committed.
- Phase scripts go in `hypernet/wiggles/` with a `lab_` prefix (run as
  `python -m hypernet.wiggles.lab_<name>`).  Reusable code goes in flat
  modules.  Committed figures go in `hypernet/wiggles/figs/lab/`.
- **Keep what you replace.**  When a committed product changes
  (`srf_table.csv`, `hypstar_grids.npz`, `release2_srf_models.json`,
  `wavescale_window_pattern.csv`), keep the previous version alongside as
  `<name>_pre_lab.<ext>` and report the differences row by row.
- Every `.md` in `docs/` is published.  Nothing unpublished goes there.
- Ask questions in the Q&A section below; log completed work under `## Logs`.

## Context

### Where things are

- **Decisions and findings:** `claude_prompts/wiggles/wiggles_prompts.md`,
  Q&A "New data #1" and its Log entry (2026-10-17).  Read these first.
- **Lab data:** `$OS_COLOR/hypernet/wiggles/wiggles_lab_data.zip`, unzipped
  to `$OS_COLOR/hypernet/wiggles/lab_data/hypstar_<sn>/calibration/wavelength/<YYYY_MM>/`.
  - 7 instruments and 13 epochs, VNIR only:
    - 120242: 2020_09, 2022_05, 2022_06, 2024_05;
    - 121222: 2021_07, 2021_07_tst;
    - 121231: 2021_07;
    - 122302: 2023_01, 2024_09;
    - 122303: 2023_01;
    - 122304: 2023_01, 2024_09;
    - 122305: 2023_01.
  - Four files per epoch:
    - `fwhm_vnir.dat`: lab FWHM_L and FWHM_E at the lamp lines.
    - `*_wl_coefs_*.dat`: 4th-order pixel→λ polynomials for VNIR_L and
      VNIR_E, with px counted from 0 over 2048 px; the header gives the
      residual σ.
    - `lamp_vnir.dat`: L/E column pairs, one per lamp (Hg-Ar, Xe, and Kr
      where present), each peak-normalised.  Its `wl` column is a fixed
      nominal per-instrument grid, probably the factory one.
    - `deltawl_vnir.dat`: reference λ − nominal-grid λ at the fitted line
      centre.
- **Exploratory scripts** (do not import them; promote what you need):
  - `wavecal/lab_data_examine.py` (Fable, 2026-10-17): inventory, file
    semantics, grid crops, FWHM comparison, line shapes, the common
    polynomial residual, stability and wings.  Its outputs are in
    `$OS_COLOR/hypernet/wiggles/lab_data/derived/` (`summary.txt` plus CSVs
    and PNGs).
    - ⚠ Its `egrid_122304_2024-11-08_2024-11-12.csv` uses the 2024_09
      polynomial.  The processor uses 2023_01 for that period, as it does
      for every 122304 period.
  - `wavecal/invert_l2_irradiance.py`: native-E recovery from the L1C
    irradiance on VEIT.
- **Code this touches:**
  - `hypernet/srf.py`: `fit_srf_template` with `via_grid`, `shape` in
    `SHAPES` (gauss / pvoigt / supergauss), `fit_template_windows`,
    `SRFModel`.
  - `hypernet/srftable.py`: `load_srf_table`, `row_for`, `model_from_row`,
    `srf_for`, `scale_for`, `method_for`, `TABLE_IS_DEFAULT`,
    `APPLY_CENTROID_OFFSET`, `SOURCE_RANK`.
  - `hypernet/srfcal.py`: `calibrate_srf`, `fit_wavelength_scale`,
    `load_window_pattern`, `correct_wavelengths`, `scale_or_fallback`.
  - `hypernet/edinterp.py`: `interpolate_ed_to_l`, `choose_method`,
    `delta_fwhm`.
  - `hypernet/twin.py`: `load_grids` reads `hypernet/data/hypstar_grids.npz`.
  - `hypernet/whn_l1a.py`: `load_l1a_irr`, `load_l1c`, `load_l2b`,
    `release2_path`.
- **Scripts that use the 122304 proxy E grid** (task 3 replaces it):
  - `hypernet/wiggles/phase1_instrument.py` writes `hypstar_grids.npz` and
    `release2_srf_models.json`;
  - `hypernet/wiggles/phase0c_release2_srf.py`;
  - `hypernet/wiggles/phase2_srf_table.py`;
  - `hypernet/wiggles/phase2_wavescale.py`.
- **Checks to re-run:** `hypernet/wiggles/phase2_veit_check.py` and
  `phase2_shape_check.py` (task 9b).
- **Committed products:**
  - `hypernet/data/srf_table.csv`, `hypstar_grids.npz`,
    `release2_srf_models.json`, `wavescale_window_pattern.csv`;
  - `hypernet/wiggles/phase0c_*.csv`, `phase2_*.csv`.
- **Gate and exit reports:** `claude_prompts/wiggles/gate_G1.md` and
  `phase2_exit.md`.  `gate_G0.md` does not exist yet (task 6 writes a
  provisional one).
- **Interim report:** `reports/wiggles_interim_report.md` (outside `docs/`,
  unpublished), for the HYPERNETS team, in JXP's voice.
- **Plan:** `docs/wiggles_planning.md` §4 (gates G0-G3) and §9 (risks).
- **Note to Joel and Kevin** about the lab data (plain text, for email):
  `correspondence/wiggles_note_lab_data.txt`.  It is drafted; JXP sends it.

### Key findings from the lab data (2026-10-17, to be re-derived in task 2)

- **ΔFWHM = L − E**, median over the 400-680 nm lines (lab | field table):

  | instrument | lab | field |
  |---|---|---|
  | 120242 | +0.44 | 0.39 |
  | 122302 | +0.47 | 0.46 |
  | 122304 | +0.39 | 0.39-0.49 |
  | 121222 | −0.04 | −0.05 |
  | 121231 | +0.02 | 0.08 |
  | 122305 | −0.06 | 0.01 |
  | 122303 | −0.02 | 0.14 |

  The split matches, and 122303 is `linear`.  The field FWHM(λ) models agree
  with the lamp widths to within about 0.1 nm.
- **E-SRF shape.**  The narrow-E units (120242, 122302, 122303, 122304,
  122305) have a super-Gaussian E with p ≈ 2.5-2.9.  L is Gaussian
  (p ≈ 2.1) on every unit, and E is too on 121222 and 121231.  There are no
  Lorentzian wings, and a slight blue-side asymmetry.  The VEIT field
  template fit gave p_E = 2.23.
- **Grids.**
  - The processor's grid is the lab polynomial on the pixels with
    350 ≤ λ < 1100 nm.  This reproduces the VEIT L1A E and L grids to 2e-5
    nm, and the Release 2 L grid of every instrument.  120242 has one
    red-end pixel in doubt.
  - **One epoch per instrument:**
    - 120242 ← 2024_05;
    - 121222 ← 2021_07;
    - 121231 ← 2021_07;
    - 122302, 122303, 122304, 122305 ← 2023_01.
  - The E epoch is verified only for 122304.  For the others it is inferred
    from the L grid.
- **The 122304 proxy E grid** is off by 0.15-1.3 nm (up to more than 2
  pixels) on the other instruments.
- **A common polynomial residual** appears on all seven instruments, in
  both channels: +0.09 nm at Hg 404.66, −0.09 nm at Hg 546.07, and ~±0.1 nm
  in the NIR.  It matches the window pattern near 545 and 585 nm but not at
  405 and 435 nm.  Between 440 and 690 nm only 546.07 constrains it (plus Kr
  587.09 on 121222, 121231 and 120242 2020_09).
- **The scale is air.**
- **Lab stability:** ≤ 0.12 nm in FWHM and ≤ 0.09 nm in the polynomial
  across epochs.
- **Stray light:** the floor is ≤ 3e-4.
- **L2 inversion.**  E_L = M·E with M fixed by the two grids.  On VEIT,
  native E is recovered to 1.5e-5 (second differences to 0.2 %) with the
  right grid.  With the 2024_09 grid the second-difference error is 18 %.

### Dependencies

- None of tasks 2-7 needs the L1A delivery.
- G0 stays provisional until Phase 0b fits the L1A.  When the delivery
  arrives, Phase 0b's `srf_models.parquet` outranks everything here
  (`SOURCE_RANK`).
- Field widths remain the source of truth for the SRF table.  The lab data
  are validation and a prior on shape.  They become the source only where
  the field has nothing, and such rows must be flagged.
- Task 3 depends on task 2 (the grids).  Task 4 depends on task 2 (the
  residual table).  Task 5 depends on task 2.  Task 6 depends on tasks 2-4.
  Task 7 depends on all of them.

## Prompts

### Setup

1. Read the Context above, Q&A "New data #1" in `wiggles_prompts.md`,
   `phase2_exit.md`, and `wavecal/lab_data_examine.py` and
   `wavecal/invert_l2_irradiance.py` in full.  Then put your questions in
   Q&A below.  Cover at least:
   - the module name (default `hypernet/labcal.py`);
   - the schema of the committed lab tables (task 2);
   - whether the SRF table should take the E shape p from the lab or from a
     field supergauss fit at the lab p's starting value;
   - whether 120242's ambiguous red-end pixel matters.

   Do not write code yet.  Log your work.

### Build

2. **Lab reader and tables (New data #1 item 1; the lab half of Phase 0
   task 13).**  `hypernet/labcal.py`:
   - **Readers** for the four file types.
   - **Committed coefficients** in `hypernet/data/lab_wl_coefs.csv`: one row
     per instrument × epoch × channel, with c0-c4, σ, the coefficient-file
     date and the lines used.
   - **`PROCESSOR_EPOCH`**: the epoch the processor uses per instrument,
     and whether that was verified (L and E) or inferred (L only).
   - **`grid_for(instrument, channel, epoch=None)`**: computes the grid from
     the coefficients with the 350 ≤ λ < 1100 crop.  JXP's choice: compute
     the grids from committed coefficients; do not commit grid CSVs.
   - **A committed line table** `hypernet/data/lab_srf.csv`: per
     instrument × epoch × channel × line, with the lab FWHM, our Gaussian
     and super-Gaussian refits (p and FWHM), S/N, and blend / low-S/N flags.
   - Promote only the fitting you need from `wavecal/lab_data_examine.py`.
   - **Tests:**
     - `grid_for` reproduces the VEIT L1A E and L grids to ≤ 1e-4 nm
       (Tier 2);
     - it reproduces every instrument's Release 2 L grid (Tier 2, one L2B
       file each);
     - the coefficient CSV round-trips (Tier 1);
     - the crop rule and pixel counts hold (Tier 1).
   - **Comparison:** `hypernet/wiggles/lab_compare.py` writes the committed
     `hypernet/wiggles/lab_fwhm_compare.csv` (lab vs `srf_table.csv` at the
     line wavelengths and ΔFWHM per instrument) and
     `figs/lab/fwhm_lab_vs_table.png`.
   - Record in the Phase 0 doc's Logs that task 13's lab half is done.  The
     cal `.nc` files have still not arrived.

   Log your work.

3. **Revise the SRF table (item 2).**
   - Regenerate `hypernet/data/hypstar_grids.npz` with every instrument's
     own E and L grids from `labcal.grid_for`, keeping the old file as
     `_pre_lab`.  Update `phase1_instrument.py` and `twin.load_grids` so the
     key names cover all instruments.  Note that the old file's
     `grid_L_122305` has 1536 px; check it.
   - Re-run `phase0c_release2_srf.py` with the own E grid in `via_grid` for
     every instrument; the `e_grid` column becomes `own`.  Report how the E
     FWHM(λ) and ΔFWHM change per row.
   - Re-estimate the "E from L2 reads ~0.06-0.09 nm wide" bias with the
     right grids.  The lamp does not support it.
   - Re-run `phase2_srf_table.py`:
     - 122303 must come out `linear` from the data.  If it does not, say so,
       and set it from the lab with a note rather than forcing it silently.
     - Add the columns `shape` and `shape_p` (and `shape_source`):
       super-Gaussian for E on the narrow-E instruments, with p per the
       Setup answer; Gaussian elsewhere.
     - Make `srftable.model_from_row` and `srf_for` honour the shape, so
       that `interpolate_ed_to_l` convolves with it.
     - Add a lab-validation column (lab − table FWHM at 405/436/546 nm).
   - Re-run `phase2_veit_check.py` and `phase2_shape_check.py` with the lab
     E shape.  Report how much of `srf`'s reduction remains on VEIT and in
     the twin, against the Gaussian rows and task 9b's 76-84 %.
   - Update `phase2_wavescale.py`'s E grids and re-run it, because the
     proxy also fed the E scale fits.
   - Run `pytest -q`; tests that pinned old numbers get updated, with the
     reason in the Log.

   Log your work.

4. **The common window pattern (item 3).**
   - Compare `hypernet/data/wavescale_window_pattern.csv` (the empirical
     0.13 nm rms per-window pattern of Phase 2 task 4b) with the lab's
     common polynomial residual, evaluated at the window centres.  Use the
     per-line residuals; between lines, use only interpolation you can
     defend, and refit the polynomial with and without the Kr lines.
   - Quantify the rms of the pattern before and after subtracting the lab
     residual, per instrument and channel, and where they disagree (405,
     435 nm).
   - If the lab residual explains a significant part, add a per-instrument
     lab correction option to `srfcal`, e.g.
     `load_window_pattern(source='empirical'|'lab'|'both')` (keep the
     empirical one as default until JXP decides).  Re-run
     `phase2_wavescale.py` with it.
   - Table `hypernet/wiggles/lab_window_pattern.csv` and a figure
     `figs/lab/window_pattern_vs_lab.png`.
   - State plainly whether Phase 2 open concern #1 is resolved, partly
     resolved, or not.

   Log your work.

5. **Validate the L2 inversion (item 4; no Phase 3).**
   - Promote the recovery to a function, e.g.
     `recover_native_e(irradiance_L, wav_L, wav_E)` in `hypernet/edinterp.py`
     (or per the Setup answer).  Handle the L pixels outside the E grid
     (the last pixel at 1099.9 nm extrapolates) and NaNs explicitly.
   - **Tests:**
     - recovery on VEIT reproduces L1A E to ≤ 1e-4 over 400-900 nm with
       `grid_for` (Tier 2);
     - the round trip M·recover(M·E) = M·E holds (Tier 1);
     - a 0.05 nm grid error is detected, e.g. as rising second-difference
       power (Tier 1, synthetic).
   - Then, on 3-5 Release 2 L2B sequences per instrument, check that the
     recovered E is physically sensible: no pixel-alternating ringing,
     smooth continuum, positive.  Fit the HSRS template on the recovered
     native-grid E, compare the widths with the `via_grid` route of task 3
     and with the lab FWHM, and compare the second-difference power with
     VEIT's.
   - This is the test of whether the processor used the inferred E epoch on
     instruments other than 122304.  A wrong grid shows up as excess
     second-difference power.
   - **Do not** recompute ρw or run any Phase 3 metric: JXP chose to wait for
     the L1A.
   - Table `hypernet/wiggles/lab_inversion_check.csv`.

   Log your work.

### Gate

6. **Provisional gate G0 (item 5).**  Write `claude_prompts/wiggles/gate_G0.md`
   (not under `docs/`), headed **provisional: Release 2 L2 + lab data, no
   L1A**.
   - Evaluate the plan's criteria (§4), with numbers from tasks 2-5, the
     Phase 0 Logs (8c-8f) and `phase2_exit.md`:
     - **G0(a):** FWHM_E vs FWHM_L, field and lab, per instrument.
     - **G0(b):** stability across SZA, sky index, season and
       recalibration (field), plus lab epochs, against 0.2 nm.
     - **G0(c):** centroid offsets; the air scale; the common polynomial
       residual; the one-epoch-per-instrument finding.
   - Give the per-instrument SRF summary (class, ΔFWHM, shape, method).
   - Say what it means for Phase 2: table vs self-calibration, any
     wavelength recalibration, and the `TABLE_IS_DEFAULT` /
     `APPLY_CENTROID_OFFSET` settings.
   - List exactly what the L1A delivery must still confirm.
   - Propose the corresponding edits to `docs/wiggles_planning.md` §4 in
     Q&A; apply them only if JXP approves.

   Log your work.

### Report

7. **Update the interim report (item 6).**  Edit
   `reports/wiggles_interim_report.md` in place, in the same voice and for
   the same readers (the HYPERNETS team).  Keep the previous version as
   `reports/wiggles_interim_report_2026-10-15.md`.
   - **Summary:** add the lab confirmation of the two instrument classes,
     the flat-topped E SRF and what it costs `srf`, the air scale, and the
     window-pattern verdict.
   - **Section 3:** a new subsection on the lab data, with a figure (lab vs
     field FWHM) and a table (ΔFWHM lab vs field, E shape p, processor
     epoch).
     - Update Table 1: 122303 moves to E ≈ L / `linear`; add a shape column.
     - Update any numbers that tasks 3-5 changed, with a note when a
       number moved.
   - **Section 4:** what comes next.  Phase 0b, G0 and Phase 3 still wait on
     the L1A, now with a provisional G0.  The SRF table is now validated
     against the lab.
   - **Section 5:** the asks.
     - #2 (lab data) is received; thank Joel.
     - #7 (air scale) is answered.
     - Add the new questions from `correspondence/wiggles_note_lab_data.txt`
       (the processor's wavelength epoch after recalibration, a 2025 lab
       epoch for 121222, the nominal `wl` grid, how the lamp illuminates the
       E diffuser).
     - Keep the L1A and cal files as ask #1.
   - Check every figure link resolves and every number against the Logs.

   Log your work.

## Q&A

### Setup #1 -- 2026-10-18 (Opus 5.5): questions before building

Read: this doc, `wiggles_prompts.md` Q&A "New data #1" and its two Logs,
`phase2_exit.md`, `wavecal/lab_data_examine.py` and
`wavecal/invert_l2_irradiance.py` in full; skimmed `srf.py`, `srftable.py`,
`srf_table.csv`, `twin.load_grids`, `edinterp.py`, `phase0c_release2_srf.py`,
`phase1_instrument.py`, the lab tree and `derived/summary.txt`.  Context
items that the code or data contradict are listed in the Log entry below;
the ones that change the work are in the questions.  Each question has a
default; it stands unless you say otherwise.

1. **Module name and contents.**  `hypernet/labcal.py` with:
   - `lab_root()` (`$OS_COLOR/hypernet/wiggles/lab_data`, the unzipped tree;
     the zip is never read);
   - readers for the four file types;
   - CSV I/O for the two committed tables;
   - `PROCESSOR_EPOCH` and `grid_for`;
   - a pixel-space lamp-line fitter (Gaussian + super-Gaussian), promoted
     from `fit_line_px` in `lab_data_examine.py`.

   The pseudo-Voigt is dropped: η = 0 on every line and every epoch.
   `grid_for` takes `channel` in {`E`, `L`, `Ld`, `Lu`}, since Ld and Lu
   share the L grid and the SRF table says `Ld`.
   *Default: yes, as listed.*
>A Use your default

2. **Schema of `hypernet/data/lab_wl_coefs.csv`.**  The epoch directory name
   has to be the key.  The coefficient-file date is not unique: 121222
   `2021_07` and `2021_07_tst` both hold `hypstar_121222_wl_coefs_210714.dat`,
   but the coefficients differ by up to 0.018 nm, and only `2021_07` matches
   Release 2.  Proposed columns, 26 rows (13 epochs × L/E):
   - `instrument` (`HYPSTAR_<sn>`), `sn`, `epoch`, `channel` (L/E);
   - `c0`…`c4`, written with `%.17g`.  c4 ≈ 3e-12 multiplies px⁴ ≈ 1.8e13,
     so a short float format breaks the 1e-4 nm test;
   - `sigma` (header value, 2 decimals);
   - `coef_file`, `coef_date` (the yymmdd in the file name), `file_created`
     (header line, blank for 120242 2020_09, 121222 2021_07 and 121231
     2021_07, which have none);
   - `lamps` (`HgAr+Xe`, `HgAr+Kr` for 120242 2020_09, `HgAr+Kr+Xe` for
     121222/121231 2021_07);
   - `n_lines`, `lines` (see Q3);
   - `processor` (bool), `processor_status` (`verified` / `inferred_L` /
     `no`; see Q8).

   The `_tst` epoch is kept with `processor = False` and is excluded from
   every comparison.  The zero SWIR columns are dropped.
   *Default: as above.*
>A. Use your default

3. **"The lines used" in `lab_wl_coefs.csv`.**  The coefficient files do not
   list the lines Joel fitted.  The evidence suggests he used all the lines
   in `deltawl_vnir.dat`:
   - our Gaussian-centre residuals through each file polynomial have rms
     0.06-0.09 nm, which matches the stated σ;
   - on the three Kr epochs, the file polynomial agrees with our all-line
     refit to ≤ 0.014 nm at the window wavelengths, so Kr was included;
   - 120242 2020_09 (σ = 0.03) has Kr and no Xe.

   Options:
   - (a) record the `deltawl_vnir.dat` line list (13/14/19 lines, including
     Hg 296.7 nm, which lies outside the 350-1100 crop), with a
     `lines_source = assumed_all_deltawl` column;
   - (b) record only the lines our refit finds at S/N ≥ 10;
   - (c) also ask Joel, as a sixth question in
     `correspondence/wiggles_note_lab_data.txt`.

   *Default: (a), plus (c) if you are happy to add one line to the unsent
   note.*
>A. (a) and (c).  Please add the line

4. **Schema of `hypernet/data/lab_srf.csv`.**  One row per instrument ×
   epoch × channel × line, for all 13 epochs (≈ 400 rows).  Columns:
   - `instrument`, `sn`, `epoch`, `channel`, `lam_ref`, `species`, `lamp`;
   - `lab_fwhm`, `lab_deltawl` (from the files);
   - `mu_px`, `mu_px_err`, `disp_nm_px`, `lam_fit`, `resid` (= poly(mu) −
     ref, the dlam sign; task 4 needs this);
   - `fwhm_gauss`, `fwhm_gauss_err`, `fwhm_halfmax`;
   - `fwhm_sg`, `sg_p`, `sg_p_err`;
   - `asym`, `asym_err` (the asymmetric Gaussian is the best model on 121222,
     121231, 122303 and 122305; see Q7);
   - `chi2_gauss`, `chi2_sg`, `chi2_asym`, `snr`;
   - flags: `blend` (Hg 404.66, Xe 823.16, which the script already treats
     as raised by a neighbour), `low_snr` (S/N < 10, mostly Ar 696/763/826,
     whose lab FWHM reads up to 1 nm high), `outside_grid` (Hg 296.7).

   **The errors are not usable as they stand.**  χ²_ν is 30-560 on the Hg
   lines but 0.1-0.2 on Ar.  The noise comes from line-free windows, so
   photon noise on the bright lines is missing.  *Default: store the formal
   errors multiplied by √χ²_ν (column suffix `_err` = scaled) and say so in
   the module docstring.*
>A. Use your default

5. **Wavelength scale for the `lab_srf.csv` refits.**  The exploratory code
   fits in pixel space and converts with the local dispersion of that
   epoch's own polynomial.  It uses the lamp files' nominal `wl` column only
   to choose the noise windows.  Options:
   - (a) keep that;
   - (b) use the nominal `wl` column.  It is rounded to 0.1 nm and sits
     0.08-0.39 nm rms off the epoch polynomials, which would bias `resid`;
   - (c) use the `PROCESSOR_EPOCH` polynomial for every epoch of an
     instrument.

   FWHM barely depends on the choice, because the dispersions agree to
   < 0.1 %.  The centres and residuals of task 4 do depend on it.
   *Default: (a).*
>A. Use your default

6. **Where the E shape p comes from in the SRF table.**  The `srf_table`
   FWHMs are Gaussian template fits.  For a flat-topped SRF, the Gaussian-fit
   FWHM is not the true FWHM: on the super-Gaussian units our Gaussian refit
   of the lamp lines is 0.03-0.10 nm below the lab value.  Pairing the lab p
   with a Gaussian-fit FWHM would therefore describe an SRF that nobody
   measured.  The code has a gap: `fit_srf_template(shape='supergauss')`
   always fits p, starting at 2.0 with bounds [1, 8], and has no way to fix
   p or set its starting value.  Options:
   - (a) take p from the lab and keep the field Gaussian FWHM(λ);
   - (b) fit field super-Gaussians with p free, started at the lab p;
   - (c) fit field super-Gaussians with p **fixed** at the lab value (new
     `shape_par` / `fix_shape_par` arguments), and refit FWHM(λ) with it.

   Option (b) is ill-conditioned: on L2 E, p trades off against width and
   against the `via_grid` resampling.  The VEIT field fit gave 2.23 (Phase 0
   weighted mean) and 2.54 (median), against a lab 2.79 for 122304.
   *Default: (c), with `shape_source = lab`.  The free-p field fit is kept
   as a diagnostic column `shape_p_field` only.*
>A. Use your default

7. **Which rows get a super-Gaussian E.**  The Context lists five
   "narrow-E" units with p ≈ 2.5-2.9: 120242, 122302, 122303, 122304 and
   122305.  But 122303 and 122305 are E ≈ L (lab ΔFWHM −0.02 and −0.06), so
   they are not narrow-E.  On their E lines the super-Gaussian hardly
   improves χ² (525 vs 562, 238 vs 233).  The asymmetric Gaussian wins
   instead (asym −0.38 and −0.21, wider on the blue side), as it does on
   121231 (−0.19).  The super-Gaussian clearly wins (χ²_ν down 2-9×) only on
   120242, 122302 and 122304.  Options:
   - (a) `supergauss` for E on 120242, 122302 and 122304 only;
   - (b) all five;
   - (c) a rule: super-Gaussian where median χ²_sg < 0.5 χ²_g, which gives
     (a).

   L is `gauss` on every unit (p 2.03-2.22, no χ² gain).  The asymmetry
   stays in `lab_srf.csv` only: `asym` is not in `SHAPES`, and it does not
   matter for the `linear` units.  *Default: (c) → (a).*
>A. Use your default

8. **122304 (and 122302): 2024_09 vs 2023_01.**  The grids come from
   `PROCESSOR_EPOCH` = 2023_01 for every period.  For widths and shape,
   the epoch the lab actually measured before the 2024-11 calibration is
   2024_09.  `lab_data_examine.py` uses 2024_09 for both grids and widths,
   because its `BEST_EPOCH` matches by date.  The FWHM difference between
   the two epochs is ≤ 0.004 nm.
   *Default:*
   - `grid_for` defaults to `PROCESSOR_EPOCH`, and `epoch=` is explicit;
   - `lab_fwhm_compare.csv` and the shape p use the latest non-`_tst` epoch
     before the period's calibration date, in a `lab_epoch` column;
   - tests:
     - the VEIT L1A E and L grids are reproduced by
       `grid_for('122304', 'E'/'L')` to ≤ 1e-4 nm;
     - a guard test: `epoch='2024_09'` misses them by > 0.02 nm rms
       (observed 0.023 E, 0.029 L);
   - nothing reads `derived/egrid_*.csv`.

   `processor_status`: `verified` for 122304 L and E; `inferred_L` for the
   other instruments' E.  Release 2 matches their L grid only.
>A. Use your default

9. **120242's red-end pixel: does it matter?**  The Release 2 L grid has
   1540 px (px 411-1950, ending at 1099.105 nm).  The 350 ≤ λ < 1100 rule
   on the 2024_05 polynomial gives 1541 px, because px 1951 = 1099.575 nm.
   One explanation fits the summary numbers: the pixel window was set with
   120242's **2020_09** polynomial, whose px 1951 is 1100.30 nm, and the
   2024_05 wavelengths were then evaluated on that window.  This is not yet
   checked for E, where the rule on 2024_05 gives 1541 px ending at
   1099.83.
   - **Science:** it does not matter.  The pixel lies above 1099 nm, outside
     every band we use (Emod 380-1000, FWHM 395-865, Phase 3 400-900).
   - **Code:** it matters in three places:
     - array lengths in `hypstar_grids.npz`;
     - the exact-length Release 2 test;
     - `recover_native_e`: the last E pixel would have no L pixel beyond
       it, so M has a zero column there.

   *Default:*
   - in task 2, test the 2020_09-window hypothesis on L.  If it holds, apply
     it to E as a documented per-instrument `CROP_EPOCH` exception, and mark
     120242 E `inferred_L`;
   - if it does not hold, add a one-pixel L exception and leave E on the
     rule;
   - either way, `recover_native_e` returns NaN for E pixels that no L pixel
     constrains.
>A. Use your default

10. **Where `recover_native_e` lives.**  `invert_l2_irradiance.interp_matrix`
    re-implements the M that `edinterp.linear_weights` already defines.
    *Default:*
    - put it in `hypernet/edinterp.py`, with `linear_weights` as the single
      definition of M, beside a small `second_difference_power(E, wav, band)`
      helper for task 5's ringing test;
    - edge rules:
      - L pixels beyond the E range (122304's last L pixel, 1099.93 nm,
        against an E maximum of 1099.70) are clamped rows.  They stay in
        the system, because they are exact equations of `np.interp`;
      - unconstrained E pixels → NaN;
      - NaN rows in E_L are dropped before `lsqr`;
    - it accepts 1-D or (n_L, n_seq) input.

    The alternative is `labcal.py`, which would make `edinterp` import nothing
    new but separate M from its inverse.
>A. Use your default

11. **`_pre_lab` files and the `hypstar_grids.npz` keys.**
    *Default:*
    - keep `<name>_pre_lab.<ext>` beside each replaced product in
      `hypernet/data/`.  No loader ever defaults to a `_pre_lab` file; they
      exist only for the row-by-row diffs;
    - in the new npz:
      - keep the existing keys with their meaning.  `grid_E_122304` and
        `grid_L_122304` are already the 2023_01 crop, and `grid_L_122305`
        (1536 px) is already the correct 122305 L crop (px 412-1947, which
        matches Release 2 exactly), so the old "check it" is answered;
      - add `grid_E_<sn>` / `grid_L_<sn>` for all seven instruments;
      - leave the noise keys unchanged;
    - a Tier-1 test that every npz grid equals `labcal.grid_for` to 1e-9 nm,
      so the two sources cannot drift.

    The alternative is to have `twin.load_grids` compute the grid keys from
    `grid_for` and keep only the noise in the npz.  That leaves one source of
    truth but changes what `phase1_instrument.py` writes.
>A. Use your default

12. **Test tiers, and committing Joel's numbers.**  The repo uses per-module
    `skipif` decorators (`needs_veit`, `needs_whn`, `needs_hsrs`), not
    registered markers.
    *Default:*
    - add `needs_lab`, which checks that the lab tree exists;
    - the readers are Tier 2;
    - Tier 1 covers the CSV round trip, the crop rule, the pixel counts,
      the npz = `grid_for` check and the synthetic `recover_native_e`
      tests;
    - no raw lab files go into `hypernet/tests/data/`.

    **A question first:** `lab_wl_coefs.csv` and `lab_srf.csv` would put
    Joel's unpublished calibration coefficients and FWHMs into a repository
    whose docs are public.  Is that cleared, or should we ask Joel in the
    note?  *Default: ask him in the note, and commit locally meanwhile; you
    decide what to push.*
>A. I will do the push.  Add to the note that I have done so.

13. **The Kr lines in task 4.**  Only three epochs have Kr: 121222 and
    121231 2021_07, which are both `PROCESSOR_EPOCH`, and 120242 2020_09,
    which is not.  Between 440 and 690 nm, the four 2023_01 instruments and
    120242 2024_05 have only Hg 546.07.  The pattern's largest windows (415
    −0.20, 505 −0.23, 605 +0.25, 665 −0.19 nm) are therefore unconstrained
    by any lab line on those units.  One lead: 120242 2020_09, fitted to
    Hg/Ar + Kr with no Xe, has σ = 0.03 nm against 0.08-0.09 at its later
    epochs.  That suggests the ±0.09 nm "common residual" comes from the
    4th-order fit of this particular line set, with Xe 823-1084 pulling it.
    It may not be instrument physics.
    *Default:*
    - for the verdict, compare only at the lines (405, 436, 546, and 587
      where available);
    - as the one defended interpolation, transfer the all-line − standard-set
      difference curve from the two Kr processor epochs to the other
      instruments.  The residual is common to ±0.01 nm across instruments,
      so the transfer is justified;
    - test the transfer on 120242 (2020_09 vs 2024_05);
    - refit each Kr epoch with and without Kr, and without Xe.
>A. Use your default

14. **Shape in `SRFModel`.**  To make `srf_for` honour the shape, the shape
    has to travel with the model.  `SRFModel` is Gaussian-only (no shape
    field), and `srf.convolution_matrix` and `srf_key` are Gaussian-only, so
    `interpolate_ed_to_l` cannot convolve with a super-Gaussian today.
    *Default:*
    - add `shape='gauss'` and `shape_par=nan` fields to `SRFModel`.  Old
      JSON loads as Gaussian, so the change is backward-compatible;
    - make `convolution_matrix` use `srf_kernel` for the non-Gaussian
      shapes;
    - add the shape to `srf_key`, so cached Gaussian matrices are never
      reused for a super-Gaussian;
    - add a Tier-1 test that `shape='supergauss', p=2` equals the Gaussian
      matrix.
>A. Use your default

## Logs

### 2026-10-18 -- Doc created

- Written from `wiggles_prompts.md` Q&A "New data #1" and JXP's answers:
  - proceed with items 1-5 in order, and add item 6 (update the interim
    report);
  - Phase 3 waits for the L1A;
  - grids are computed from committed coefficients.
- Items 1-6 map to tasks 2-7, with Setup as task 1.

### 2026-10-18 -- Setup #1 (Opus 5.5)

- **Read in full:**
  - this doc;
  - `wiggles_prompts.md` Q&A "New data #1" and its Logs of 2026-10-17 and
    2026-10-18;
  - `phase2_exit.md`;
  - `wavecal/lab_data_examine.py` (965 lines);
  - `wavecal/invert_l2_irradiance.py`.
- **Skimmed:**
  - `srf.py`: `srf_kernel`, `SHAPES`, `fit_srf_template`, `SRFModel`,
    `convolution_matrix`;
  - `srftable.py` and `srf_table.csv` (the columns, and the `e_grid` values);
  - `twin.load_grids`, `edinterp.py`;
  - `phase0c_release2_srf.py`, `phase1_instrument.py`;
  - the test skip conventions;
  - the lab tree (heads of the four files for 122304 2023_01, plus all
    coefficient headers);
  - `lab_data/derived/summary.txt`.
- Wrote 14 questions in Q&A "Setup #1".  **No code written**, no package
  file changed, and `pytest` was not run.
- **File semantics learned:**
  - `fwhm_vnir.dat` and `deltawl_vnir.dat` have four columns
    (`wl, L, wl, E`).  The two `wl` columns are identical in every file.
    Lines are not sorted: the Kr and Xe lines are appended after Hg 1013.98.
  - There are 13 standard lines (Hg 296.7 to Xe 1083.8), 14 on 120242 2020_09
    (Hg/Ar + Kr, **no Xe**) and 19 on 121222 and 121231 2021_07 (Hg/Ar + Kr
    + Xe).
  - Coefficient files are 5 rows × `VNIR_L VNIR_E SWIR_L SWIR_E` (SWIR all
    0) and give σ to 2 decimals.
    - Three files have no "File was created" line.
    - The 2024_09 files (named 240926) were created 2024-11-07, the day
      before the 2024-11-08 radiometric calibration, yet the processor kept
      2023_01.
    - The `_tst` coefficient file has the same name as `2021_07`'s but
      different coefficients (created 2022-04-12).  The epoch directory is
      the only unique key.
  - `lamp_vnir.dat` has 2048 rows (one per pixel) and 2 or 3 L/E pairs,
    each column peak-normalised.
- **Where the reusable pieces are in `lab_data_examine.py`:**
  - `read_coefs`, `read_two_col`, `read_lamp`;
  - `poly`, `dpoly`, `px_of`;
  - `crop_rule`, which is the processor rule;
  - `fit_line_px` with `gauss` / `sgauss` / `skew_gauss`, `noise_of` and
    `halfmax_fwhm`;
  - `lamp_pairs` / `line_pair`, which identify the lamp of each pair by its
    brightest line.

  The script runs top to bottom at import (no `main`) and must not be
  imported.  `invert_l2_irradiance.interp_matrix` duplicates
  `edinterp.linear_weights`.
- **Inconsistencies between the Context and the code or data** (none of
  them changes the headline findings):
  1. **The proxy.**  The script's "proxy" is the 122304 **2024_09** E
     polynomial (`proxy_c`).  The proxy that phase0c and phase1 actually
     used is `grid_E_122304`, the VEIT L1A E grid, which equals the
     **2023_01** crop exactly.  The two differ by 0.023 nm rms (0.053 max),
     so the "0.15-1.3 nm off" numbers need re-deriving against the right
     proxy in task 3.
  2. **122304's 2023-01 E row.**  `srf_table.csv` labels the 122304
     2023-01 E row `proxy (122304)`, but it already used its own grid (the
     same 2023_01 polynomial).  `phase0c`'s `own` flag tests for cal
     2024-11-12 only.
  3. **The "narrow-E" group.**  The Context's "narrow-E units (120242,
     122302, 122303, 122304, 122305)" mixes two groups.  122303 and 122305
     are E ≈ L.  Their E is best fitted by an asymmetric Gaussian, not a
     super-Gaussian (Q7).  The super-Gaussian wins only on 120242, 122302
     and 122304, with χ²_ν down 2-9× (not "2-5×": 122302 is 9×).
  4. **The asymmetry.**  "Slight blue-side asymmetry" understates it on
     122303 E (−0.38), 122305 E (−0.21) and 121231 (−0.19 in both
     channels).
  5. **Lab polynomial stability** is ≤ 0.10 nm, not ≤ 0.09: 120242
     2022_05 − 2020_09 L at 700 nm is −0.101.
  6. **The Kr lamp.**  "Hg-Ar, Xe, and Kr where present" misses that
     120242 2020_09 has Kr **instead of** Xe.  Its σ = 0.03 nm, against
     0.07-0.11 for every other epoch, is a lead for task 4 (Q13).
  7. **The 120242 red-end pixel.**  The R2 L grid (px 411-1950) is exactly
     the 350 ≤ λ < 1100 window of the **2020_09** polynomial, carrying
     2024_05 wavelengths.  This is a candidate explanation, still to be
     checked (Q9).
  8. **The pattern rms.**  The committed `wavescale_window_pattern.csv` has
     rms 0.119 nm over all windows; "0.13 nm" is the figure quoted from
     Phase 2.
  9. **`grid_L_122305` (1536 px)** is correct: the lab 2023_01 crop gives
     px 412-1947, matching Release 2.  The task 3 "check it" is answered.
  10. **The lab ΔFWHM.**  ΔFWHM over 400-680 nm rests on three lines (405,
      436, 546, plus 587 on the Kr epochs), and Hg 404.66 is a blend.  The
      Ar 696/763/826 lab FWHMs are inflated by up to 1 nm at S/N 2-30.
  11. **The line-fit χ²_ν** is 30-560 on Hg, so the formal p errors are
      not meaningful (Q4).
- **Code gaps the build will hit:**
  - `fit_srf_template` cannot fix p or set its starting value;
  - `SRFModel`, `convolution_matrix` and `srf_key` are Gaussian-only
    (Q6, Q14).

### 2026-10-18 -- Build #2 (Opus 5.5)

Task 2 (lab reader and tables), with JXP's Setup answers: the defaults
throughout, plus Q3 (c) and Q12 (the note now says the tables are in the
public repository).  Tasks 3 and 5 were not started.

- **New `hypernet/labcal.py`:**
  - `lab_root`, `lab_available`, `epochs`, `epoch_dir`;
  - readers `read_coefs`, `read_line_file` / `read_fwhm` / `read_deltawl`,
    `read_lamp`, `read_epoch`;
  - `poly`, `dpoly`, `px_of`, `crop_window`;
  - `PROCESSOR_EPOCH` (epoch plus L/E status), `CROP_EPOCH`, `coefs`,
    `grid_for(instrument, channel, epoch=None, return_px=False)`;
  - CSV I/O: `load_wl_coefs` / `write_wl_coefs` and `load_lab_srf` /
    `write_lab_srf`;
  - the pixel-space line fitter: `gauss`, `supergauss`, `asym_gauss`,
    `fit_line_px`, `fit_epoch_lines`, `lamp_names`, `noise_of`,
    `halfmax_fwhm`;
  - `coef_rows`, and `physical_epoch` (Q8).

  The module docstring documents the file semantics and the
  √χ²_ν-scaled errors.  No `wavecal/` script is imported.
- **New `hypernet/data/lab_wl_coefs.csv`** (26 rows) and
  **`hypernet/data/lab_srf.csv`** (364 rows = 182 lines × 2 channels; every
  fit converged), both written by the new
  `python -m hypernet.wiggles.lab_tables`.  Two runs gave byte-identical
  files.  The script's grid check goes to
  `$OS_COLOR/hypernet/wiggles/lab_data/lab_tables/`
  (`grid_check.csv`, `crop_hypothesis_120242.csv`).
- **New `hypernet/wiggles/lab_compare.py`.**  It writes the committed
  `lab_fwhm_compare.csv` (per table row × line, with `lab_epoch` and
  `processor_epoch`), `lab_fwhm_compare_summary.csv` (per table row) and
  `figs/lab/fwhm_lab_vs_table.png`.  The per-row summary is a second file
  because its schema differs from the per-line table.  Intermediates
  `lab_dfwhm_epochs.csv` and `lab_shape_clean.csv` go to the same
  `$OS_COLOR` directory.
- **Q9 verdict: the 2020_09-window hypothesis holds.**
  - The 350 ≤ λ < 1100 window of 120242's 2020_09 L polynomial is
    px 411-1950.  With the 2024_05 wavelengths on it, it reproduces the
    Release 2 L grid exactly (1540 px, max |Δ| 3.9e-5 nm).
  - The 2022_05, 2022_06 and 2024_05 windows give 1542, 1541 and 1541 px.
  - `CROP_EPOCH = {'120242': '2020_09'}` is applied to the processor epoch
    of both channels.  120242 E thus has 1539 px (px 412-1950, against
    1541 by the plain rule); it is unverified and marked `inferred_L`.
- **Grid checks** (`grid_for` against the processor grids):
  - VEIT L1A: E max 2.5e-5 nm, L 9e-6 nm.
  - Release 2 L, every instrument × period in the request: ≤ 3.9e-5 nm.
  - The 2024_09 guard: 0.023 nm rms on E.  On L, the 2024_09 crop has
    1537 px (not 1538), and the polynomial is 0.029 nm rms off on the
    processor pixels.
- **Lab vs table ΔFWHM** (lab = median L − E over the 400-680 nm lines
  with S/N ≥ 10 at `lab_epoch`; "lab − table" = median over the row's lines,
  E / Ld):

  | instrument | period | lab epoch | lab ΔFWHM | lab, no 404.66 | table ΔFWHM | lab − table FWHM, E / Ld | method (both) |
  |---|---|---|---|---|---|---|---|
  | 120242 | 2024-11 | 2024_05 | +0.436 | +0.430 | 0.390 | +0.079 / +0.035 | srf |
  | 121222 | 2021-10 | 2021_07 | −0.039 | −0.045 | −0.054 | +0.023 / +0.031 | linear |
  | 121222 | 2025-03 | 2021_07 | −0.039 | −0.045 | −0.025 | −0.016 / +0.018 | linear |
  | 121231 | 2021-10 | 2021_07 | +0.018 | +0.016 | 0.082 | −0.025 / −0.059 | linear |
  | 122302 | 2023-01 | 2023_01 | +0.469 | +0.428 | 0.464 | +0.081 / +0.007 | srf |
  | 122303 | 2023-01 | 2023_01 | −0.024 | −0.043 | 0.138 | +0.081 / +0.102 | linear |
  | 122304 | 2023-01 | 2023_01 | +0.395 | +0.376 | 0.491 | +0.086 / +0.003 | srf |
  | 122304 | 2024-11 (R2) | 2024_09 | +0.389 | +0.370 | 0.390 | +0.027 / −0.008 | srf |
  | 122304 | 2024-11 (VEIT L1A) | 2024_09 | +0.389 | +0.370 | 0.477 | +0.155 / +0.037 | srf |
  | 122305 | 2023-01 | 2023_01 | −0.062 | −0.047 | 0.011 | +0.127 / −0.015 | linear |

  - The G1 method from the lab agrees with the table on every row.
  - The lab E is wider than the table E by +0.02 to +0.13 nm on most rows.
    The VEIT L1A E row is the widest gap (+0.155), and 122303 Ld drops to
    −0.068 without Hg 404.66.
- **E shape on the processor epochs**, clean lines with S/N ≥ 10, median
  (χ²_sg / χ²_g in brackets):

  | instrument | p | χ²_sg / χ²_g |
  |---|---|---|
  | 120242 | 2.53 | 0.49 |
  | 122302 | 2.93 | 0.09 |
  | 122304 | 2.76 | 0.50 |
  | 122303 | 3.00 | 0.93 |
  | 122305 | 2.41 | 1.02 |
  | 121222 | 2.13 | 1.07 |
  | 121231 | 2.09 | 1.10 |

  ⚠ for task 3: 120242 (0.491) and 122304 (0.496) sit just under Q7's
  0.5 threshold.  The rule gives the intended answer (a), but only
  narrowly.  Dropping the weak Ar lines also moves p by up to 0.1 from
  `lab_data_examine.py`'s values (e.g. 122304: 2.76 vs 2.79).
- **Tests:** new `hypernet/tests/test_labcal.py`, with 7 tests.
  - Tier 1 (4 tests): the coefficient table's contents and
    `PROCESSOR_EPOCH`; a bit-exact CSV round trip; the crop rule and
    pixel counts per instrument × channel; the `lab_srf.csv` contents and
    the two-class split at Hg 546.
  - Tier 2 (3 tests): the readers against the committed coefficients
    (`needs_lab`); VEIT L1A E/L ≤ 1e-4 nm plus the 2024_09 guard
    (`needs_veit`); one Release 2 L2B L grid per instrument (`needs_whn`).
  - `pytest -q`: **144 passed** (137 before).
- **Correspondence:** `correspondence/wiggles_note_lab_data.txt` gained
  question 6 (which lines Joel fitted) and a paragraph saying the two
  derived tables are in the public repository, so he can object.  It is
  still unsent.
- **Phase 0 doc:** the Logs record that task 13's lab half is done and that
  the cal `.nc` files have not arrived.
- **Learned:**
  - pandas' default CSV float parser is not exact.  `load_wl_coefs` needs
    `float_precision='round_trip'` for the coefficients to survive a write
    and read bit for bit.  The writer uses Python `repr`.
  - Release 2 fixed 120242's pixel window at its first lab calibration
    (2020_09) while updating the wavelengths, so "one polynomial per
    instrument" is not the whole processor rule.
  - Not changed in this task: `srf_table.csv`, `hypstar_grids.npz`,
    `release2_srf_models.json`, `wavescale_window_pattern.csv`.

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

### Gate #6 -- 2026-10-18 (Opus 5.5): proposed §4 edits

From the provisional G0 (`claude_prompts/wiggles/gate_G0.md`).  None is
applied: `docs/wiggles_planning.md` is unchanged until you approve.  Each
edit quotes the current text of §4 and gives the replacement; answer per
edit (yes / no / changed text).

1. **G0(a): make it per instrument.**  As written, 4 of 7 instruments are
   E ≈ L ("most"), so the criterion would drop H2, which is required on the
   other three.

   *Current:*
   > (a) If FWHM_E(λ) = FWHM_L(λ) within the fit uncertainty on most
   > instruments, H2 is dropped and eq. (14) suffices.

   *Proposed:*
   > (a) Per instrument: where ΔFWHM = FWHM_L − FWHM_E (median over 400–680
   > nm, Gaussian-equivalent widths) is below 0.15 nm (gate G1), H2 is not
   > needed and the processor's linear interpolation is kept; where it is
   > above, the `srf` method applies.  Lab line-spread data, where they
   > exist, cross-check the class; a disagreement is resolved from the L1A
   > E, not silently.
>A. Ok use your Proposed

2. **G0(b): add the lab epochs, and bound the recalibration jump.**  The
   cal files may never "explain" a jump.  Every observed jump is
   ≤ 0.03 nm.

   *Current:*
   > (b) If the SRFs are stable per instrument (spread across SZA, sky and
   > season below ~0.2 nm, and no jump at recalibration that the cal files do
   > not explain), Phase 2 ships a per-instrument SRF table.

   *Proposed:*
   > (b) If the SRFs are stable per instrument (spread across SZA, sky and
   > season below ~0.2 nm; any jump at recalibration below ~0.1 nm or
   > explained by the cal files; lab calibration epochs, where available,
   > agreeing to ~0.2 nm), Phase 2 ships a per-instrument SRF table.

   The rest of (b) is unchanged.
>A. Ok use your Proposed

3. **G0(c): the dispersion, not only the centroid.**  The offsets are
   ≤ 0.10 nm, but shift + stretch reaches 0.1-0.24 nm at the ends of
   390-870 nm.  The processor grid is now known from the lab polynomials.

   *Current:*
   > (c) If centroid offsets exceed ~0.1 nm systematically on some
   > instrument, a wavelength recalibration goes into the function before
   > anything else.

   *Proposed:*
   > (c) If centroid offsets, or the wavelength-scale shift and stretch,
   > exceed ~0.1 nm systematically anywhere in 390–870 nm on some
   > instrument, a wavelength correction (per-sequence shift + stretch
   > against the model irradiance, with a calibration-period fallback) goes
   > into the function before anything else.  The processor's grids are the
   > lab wavelength polynomials (`hypernet.labcal.grid_for`); the correction
   > is relative to them.
>A. Ok use your Proposed

4. **Phase 0: the SRF shape.**  The lab shows flat-topped E SRFs on the
   narrow-E units.

   *Current* (first bullet of Phase 0):
   > Parameterise each channel's SRF as a Gaussian whose FWHM is a quadratic
   > in λ, fitted to ~10 clean solar/telluric lines per spectrum

   *Proposed:*
   > Parameterise each channel's SRF as a Gaussian whose FWHM is a quadratic
   > in λ (with a super-Gaussian shape where laboratory line-spread data
   > show one, its exponent from the lab and checked against the field),
   > fitted to ~10 clean solar/telluric lines per spectrum

   The rest of the bullet is unchanged.
>A. Ok use your Proposed

5. **A status line under *Gate G0*.**

   *Proposed* (new last sentence of the Gate G0 paragraph):
   > A provisional evaluation from Release 2 L2 and the laboratory data
   > (`claude_prompts/wiggles/gate_G0.md`, 2026-10-18) passes (a) per
   > instrument and (b), and triggers (c), handled by the wavelength-scale
   > correction; the final G0 waits for the L1A.
>A. Ok use your Proposed

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

### 2026-10-18 -- Build #3 (Opus 5.5)

Task 3 (revise the SRF table), with the Setup defaults (Q6, Q7, Q8, Q11,
Q14).  Tasks 4 and 5 were not started; no Phase 3 metric was run and ρw was
not recomputed.

**Code**

- `hypernet/srf.py`:
  - `SRFModel` gains `shape='gauss'` and `shape_par=nan`.  JSON written
    before these fields loads as Gaussian, so the change is
    backward-compatible.
  - `ShapedFWHM`: an uncached callable that carries a shape.
  - `srf_shape(spec)` returns a spec's shape.
  - `srf_key` includes the shape.
  - `convolution_matrix` and `convolve_to_grid` use `srf_kernel` for
    non-Gaussian shapes (`shape=` overrides the spec), and the cache key
    includes the shape.
  - `fit_srf_template` and `fit_template_windows` take `shape_par` (a
    starting value) and `fix_shape_par` (hold p).
- `edinterp._fwhm_fn` and the snippet's MC `_srf` keep the shape through
  `ShapedFWHM`, so the MC draws convolve with it.
- `srftable.model_from_row` reads `shape` / `shape_p`, so `srf_for` and
  `interpolate_ed_to_l` convolve with the row's shape.
- `twin.observe` and `rhow_true` honour an SRFModel's shape.  Before, they
  used only `.fwhm`, and a table SRF was silently observed as a Gaussian.
- `labcal`:
  - `lab_shape` (median shape over the clean lines with S/N ≥ 10);
  - `srf_shape_for` (Q7 class on the processor epoch, p at the physical
    epoch, L always Gaussian);
  - constants `CLEAN_LINES` and `SG_CHI2_RATIO = 0.5`.
- **New scripts:**
  - `hypernet/wiggles/lab_shape_fits.py`: fixed-p and free-p super-Gaussian
    E fits, on VEIT L1A and on the 129 Release 2 sequences of 120242,
    122302 and 122304;
  - `lab_l2_bias.py`: the bias re-estimate;
  - `lab_pre_lab_diff.py`: every product against its `_pre_lab` copy.
    Its outputs go to `$OS_COLOR/.../lab_data/pre_lab_diff/`.
- **Modified scripts:**
  - `phase1_instrument.py` (`--grids-only` / `--models-only`): grids for
    all seven instruments; `E_sg` models;
  - `phase0c_release2_srf.py` and `phase2_wavescale.py`: own E grids;
  - `phase2_srf_table.py`: shape and lab columns; the 122303 rule;
  - `phase2_veit_check.py`: an `srf_labshape` row;
  - `phase2_shape_check.py`: `lab_E_supergauss` cases and a `shaped`
    variant.

**Products**, each with its `*_pre_lab` copy kept beside it.  The
intermediates `release2_{template,models}_pre_lab.parquet` and
`wavescale_per_sequence_pre_lab.parquet` are kept under `$OS_COLOR`.

1. **`hypstar_grids.npz`.**  It now has `grid_E_<sn>` and `grid_L_<sn>`
   for all seven instruments, from `labcal.grid_for`.
   - E has 1539/1538/1538/1538/1537/1536/1535 px and L
     1540/1538/1538/1539/1538/1538/1536 px (order 120242, 121222, 121231,
     122302, 122303, 122304, 122305).
   - The old keys keep their meaning.  `grid_E_122304` changes by at most
     2.5e-5 nm, `grid_L_122304` by 9e-6 and `grid_L_122305` by 6e-6: the
     old values were the float32 L1A/L2B wavelengths.  So none is
     unchanged bit for bit; all agree to float32 precision.
   - The noise keys are identical.
   - `grid_L_122305` = 1536 px is correct.
2. **phase0c** (all 224 sequences OK), own E grid, `e_grid` = `own`
   everywhere:
   - Ld is unchanged (Δ = 0), and 122304 is unchanged in both periods.  It
     already used its own grid; the 2023-01 rows were only mislabelled
     `proxy`.
   - The E change in mean FWHM at 400 / 500 / 600 / 700 nm:

     | instrument | period | 400 | 500 | 600 | 700 |
     |---|---|---|---|---|---|
     | 120242 | | +0.055 | +0.026 | +0.005 | −0.007 |
     | 121222 | 2021 | +0.016 | +0.010 | +0.006 | +0.007 |
     | 121222 | 2025 | +0.017 | +0.012 | +0.010 | +0.010 |
     | 121231 | | +0.019 | +0.035 | +0.044 | +0.046 |
     | 122302 | | +0.005 | +0.013 | +0.012 | +0.001 |
     | 122303 | | −0.004 | −0.015 | −0.022 | −0.023 |
     | 122305 | | +0.041 | −0.007 | −0.033 | −0.038 |

   - ΔFWHM (Ld − E) at 450 / 600 nm:

     | instrument | period | 450 | 600 |
     |---|---|---|---|
     | 120242 | | 0.529 → 0.490 | 0.451 → 0.446 |
     | 121222 | 2021 | −0.027 → −0.039 | −0.038 → −0.044 |
     | 121222 | 2025 | +0.002 → −0.012 | −0.010 → −0.019 |
     | 121231 | | 0.040 → 0.012 | 0.070 → 0.025 |
     | 122302 | | 0.523 → 0.513 | 0.426 → 0.415 |
     | 122303 | | −0.013 → −0.003 | 0.033 → 0.055 |
     | 122305 | | 0.070 → 0.056 | −0.016 → 0.018 |

   - **Why so small:** with `via_grid` the model is evaluated at the grid's
     own (true) wavelengths and then interpolated, so a wrong E grid changes
     only the resampling phase, not the wavelength scale.  The proxy's
     0.15-1.3 nm same-pixel offsets were therefore never a width error of
     that size.
3. **`release2_srf_models.json`:**
   - `e_grid` = `own`; the E FWHMs move as in item 2;
   - new `E_sg` models (super-Gaussian, p fixed at the lab value) for
     120242@2024-11 (p 2.527), 122302@2023-01 (2.934), 122304@2023-01
     (2.757) and 122304@2024-11 (2.805);
   - `E` stays the Gaussian fit, which the twin and wavescale use.
4. **New `hypernet/data/veit_srf_model_lab_shape.json`** (VEIT L1A E,
   p 2.805, fixed).  FWHM at 400 / 550 / 850 nm is 2.614 / 2.501 / 2.863
   nm, against 2.346 / 2.212 / 2.552 for the Gaussian.  The free-p median
   is 2.544 (27 windows), and χ²_ν is the same (3.59 vs 3.60).
5. **`wavescale_window_pattern.csv`** changes by ≤ 0.0025 nm per window
   (rms 0.133 nm, correlation 0.79 as before).  The E scale medians
   (`phase2_wavescale.csv`) move ≤ 0.003 nm.  The proxy did not feed the
   common pattern.
6. **`srf_table.csv`**, E rows.  The Ld FWHMs are unchanged, and the scale
   columns move ≤ 0.01.

   | row | shape (p) | fwhm_550 old → new (Gaussian-equivalent) | ΔFWHM old → new (shaped) | method |
   |---|---|---|---|---|
   | 120242 2024-11 | sg 2.527 | 2.190 → 2.363 (2.226) | 0.390 → 0.353 (0.215) | srf |
   | 121222 2021-10 | gauss | 2.658 → 2.661 | −0.054 → −0.057 | linear |
   | 121222 2025-03 | gauss | 2.713 → 2.735 | −0.025 → −0.047 | linear |
   | 121231 2021-10 | gauss | 2.791 → 2.839 | 0.082 → 0.035 | linear |
   | 122302 2023-01 | sg 2.934 | 2.161 → 2.469 (2.188) | 0.464 → 0.438 (0.156) | srf |
   | 122303 2023-01 | gauss | 2.613 → 2.576 | 0.138 → **0.173** | **linear, from the lab** |
   | 122304 2023-01 R2 | sg 2.757 | 2.278 → 2.518 (2.278) | 0.491 → 0.491 (0.240) | srf |
   | 122304 2024-11 R2 | sg 2.805 | 2.342 → 2.609 (2.342) | 0.390 → 0.390 (0.123) | srf |
   | 122304 2024-11 VEIT | sg 2.805 | 2.212 → 2.501 (2.212) | 0.477 → 0.477 (0.189) | srf |
   | 122305 2023-01 | gauss | 2.691 → 2.674 | 0.011 → 0.026 | linear |

   - New columns: `shape`, `shape_p`, `shape_source` (= `lab`),
     `shape_p_field`, `fwhm_gauss_550`, `lab_epoch`,
     `lab_minus_table_{405,436,546}`, `lab_fwhm_def`,
     `delta_fwhm_shaped`, `lab_dfwhm`, `method_lab`, `method_from`.
   - `shape_p_field` (free-p median): 120242 2.514, 122302 2.856, 122304
     2.505 (2023-01), 2.430 (2024-11 R2), 2.544 (VEIT).  This is below the
     lab p on 122304.
   - Lab validation (like with like, 436 / 546 nm), for the super-Gaussian
     E rows:

     | row | 436 nm | 546 nm |
     |---|---|---|
     | 120242 | −0.026 | −0.043 |
     | 122302 | −0.014 | +0.018 |
     | 122304 2023 | −0.080 | −0.003 |
     | 122304 2024-11 R2 | −0.128 | −0.087 |
     | 122304 VEIT | −0.065 | +0.021 |

     Gaussian E rows: −0.05 to +0.13; Ld rows: −0.08 to +0.05.
- **122303 verdict: it does not come out `linear` from the data.**  With
  its own E grid, its Gaussian E narrows by 0.02-0.04 nm and ΔFWHM rises
  from 0.138 to 0.173, which is above the 0.15 threshold.  The lab ΔFWHM is
  −0.024, so the method is set from the lab, with a note and
  `method_from = lab`.  `method_sensitive` is True, because the L1A-bias
  estimate says `srf`.  This is the only row where data and lab disagree.
- **A decision JXP should review: G1 on a shaped row uses the
  Gaussian-equivalent ΔFWHM.**
  - The FWHM of a flat-topped E is 0.25-0.30 nm larger than its
    Gaussian-equivalent width.
  - Taken literally, ΔFWHM from the shaped FWHM (`delta_fwhm_shaped`)
    flips 122304 @ 2024-11 R2 to `linear` (0.123) and leaves 122302 at
    0.156.  The lab says 0.39-0.47 for both.
  - G1's 0.15 nm threshold was set with Gaussian SRFs, and the lab's
    quoted FWHM behaves like a Gaussian-fit width: lab − our Gaussian refit
    ≈ −0.05 nm, lab − our super-Gaussian refit ≈ −0.15 nm at 546 E.  So
    `delta_fwhm` (and `method`) use the Gaussian template E, while the
    convolution uses the shaped SRF.
  - For the same reason, the lab-validation column compares a
    super-Gaussian row with our super-Gaussian refit of the lamp lines
    (`fwhm_sg`), not the lab's quoted value.  Compared with the quoted
    value, the shaped rows read 0.13-0.25 nm too wide.
- **The "E from L2 reads 0.06-0.09 nm wide" bias, re-estimated**
  (`lab_l2_bias.csv`, Hg 436 / 546; bias = (field − lab)_E − (field − lab)_L):
  - The median over the 9 Release 2 models is −0.055 nm with own grids
    (−0.061 with the proxy).  This is the opposite sign to the claim.
  - The lab's FWHM definition alone gives −0.048 (lab − our Gaussian refit,
    E minus L).  The residual, −0.007 nm, is consistent with zero.
  - **The lamp shows no L2 inflation of E.**  The +0.064 nm of
    `phase0c_validation` is the single VEIT sequence's L2-vs-L1A
    difference, and VEIT's L1A Gaussian E reads narrow against the lamp
    (−0.129).
  - `delta_fwhm_l1a_est` is kept as that one-sequence estimate, with a note
    in the docstring.  Whether to drop it, and so `method_sensitive`, is
    JXP's call.
- **How much of srf's gain remains:**
  - **Twin** (`phase2_shape_check.csv`, 25 scenes, 122304 grids; truth E
    super-Gaussian p 2.805, L Gaussian):
    - a Gaussian correction keeps 63 % (fitted widths) or 49 % (true FWHM);
    - with ΔFWHM 0.2 that falls to 45 % / 27 %;
    - a correction that knows the E shape (`shaped`) keeps 100 % in both
      cases.
    - Lab p and an L that stays Gaussian are thus worse for a Gaussian
      correction than task 9b's 76-84 % / 66-78 %.  The old cases
      reproduce exactly.
  - **VEIT** (`phase2_veit_check.csv`; line-region residual over the
    away-from-line level, 1 = no excess):
    - linear 2.45, srf (Gaussian) 1.60 (59 % of the excess removed),
      `srf_labshape` 2.54 (−6 %).
    - By line, the lab shape improves the solar lines further: Ca K 0.0090
      → 0.0087, Ca H 0.0076 → 0.0071, G band 0.0088 → 0.0075, Mg b 0.0037
      → 0.0025, Hα 0.0033 → 0.0029.
    - It makes the O₂ bands much worse: O₂-A 0.0076 → 0.0198, O₂-B 0.0028
      → 0.0048.  Na D is also worse (0.0025 → 0.0031).
    - The free-p field fit prefers p ≈ 2.5, against the lab's 2.8.  The
      O₂-A result says the lab p is too flat for sunlight on this unit, or
      that the shape interacts with Emod's saturated band.  This needs
      looking into before `srf_labshape` is the default.
- **Tests:** `pytest -q` gives **151 passed** (144 before).
  - New tests:
    - `test_srf.py`: fixed and starting p; super-Gaussian p = 2 matrix =
      Gaussian (atol 1e-13); shaped convolution and cache keys; JSON
      backward compatibility;
    - `test_labcal.py`: npz = `grid_for` to 1e-9 nm; the shape classes;
    - `test_srftable.py`: the rows carry the lab shape.
  - Updated tests:
    - `test_srftable.py::test_table_contents`: G1 now holds except where
      `method_from = lab` (122303);
    - `test_srftable.py::test_method_for_matches_the_rule`:
      `choose_method` on a shaped E is no longer G1, so the test checks
      `delta_fwhm` against `method_for`;
    - `test_twin.py`: 122304 now has `E_sg`.
  - `test_snippets::test_phase1_scene_through_the_snippet` failed until
    `twin.observe` honoured the shape.  It passes without being changed.
- **Not re-run:** `phase0c --validate` (the VEIT grid is unchanged) and
  `phase0d_refine`.  Phase 0d reads `release2_template.parquet`; its
  committed CSVs are from the pre-lab parquet.  Nothing that is Phase 1
  only was re-run.
- **Learned:**
  - The proxy E grid hardly mattered: the E widths move ≤ 0.055 nm and the
    scales ≤ 0.003 nm.
  - The real consequences of the lab data are the 122303 class and the
    shape.
  - Holding p at the lab value is consistent with the lamp when like is
    compared with like.  But on VEIT it buys solar-line gains at the cost of
    O₂-A.

### 2026-10-18 -- Build #4 (Opus 5.5)

Task 4 (the common window pattern vs the lab polynomial residual), with
Setup Q13's default.  Task 5 was not started, and no Phase 3 metric was run.

**Method** (new `hypernet/wiggles/lab_window_pattern.py`; module docstring):

- **Sign.**  The processor grid is the lab polynomial, so at a lamp line
  `resid` = poly(mu) − λ_ref is the reported-minus-true error r of the
  pixel that sees the line.  The field pattern is −dlam, so the lab
  prediction of the pattern is −r.
- **Detrending.**  The field pattern is a residual about each sequence's
  shift + stretch line, so the lab curve is compared after removing a
  straight line over the same 27 windows.
- **Between lines**, only the Q13 transfer is used, in two variants:
  - **v1** (Q13's default): r = poly_std − poly_all for the standard-set
    instruments (13 lines), the mean over 121222 and 121231 2021_07.  It is
    0 between lines for the Kr instruments.
  - **v2** (the Q13 lead): truth = the 4th-order fit *without Xe*.  Then
    r = poly_std − poly_noXe for the standard-set instruments, and
    poly_all − poly_noXe for the Kr ones.
  - The curves agree between instruments and channels to ≤ 0.005 nm (v1)
    and ≤ 0.012 nm (v2), which justifies the transfer (Q13's ±0.01).

**Refits of the Kr epochs.**  Unweighted 4th order through our line
centres; the file polynomials equal the all-line fit to ≤ 0.015 nm.

| epoch | all lines | no Kr | no Xe |
|---|---|---|---|
| 121222 L / E | 0.068 / 0.064 (19) | 0.073 / 0.068 (13) | **0.028 / 0.029** (14) |
| 121231 L / E | 0.061 / 0.071 | 0.066 / 0.077 | **0.026 / 0.027** |
| 120242 2020_09 L / E (has no Xe) | 0.028 / 0.026 (14) | 0.019 / 0.012 (8) | n/a |

Values are rms residuals in nm, with the number of lines in brackets.
`sigma` = sqrt(Σr²/(n−5)) is 0.071-0.098 (all, no Kr) and 0.032-0.037 (no
Xe).

- Residuals at 405 / 546 nm: all lines +0.11 / −0.08 (both Kr epochs);
  no Kr +0.09 / −0.10; no Xe +0.045 / −0.04.
- The no-Xe fit cannot reach Xe 1083.8 (+0.57 to +0.69 nm when
  extrapolated).  It is valid over 297-1014 nm.

**The Q13 lead: confirmed in part.**

- Dropping Kr barely changes the fit.  Dropping Xe cuts the rms from
  0.061-0.071 to 0.026-0.029 nm.  So the ±0.09 nm "common residual" is a
  4th-order fit spanning 297-1084 nm with the Xe 823-1084 lines, not
  instrument physics.
- It is still a real error of the processor's grid, because the processor
  uses that polynomial.
- **120242 2020_09** (fitted without Xe, σ 0.03) against 2024_05 (Xe,
  σ 0.09): the same-pixel difference poly_2024_05 − poly_2020_09 is
  +0.019 at 405, −0.083 at 546, −0.107 at 600, −0.091 at 700 and +0.022
  at 882 nm (L).  It has the shape of the v2 curve: rms 0.059, falling to
  0.030 after subtracting v2 but only to 0.053 after v1.  The 2020→2022
  "polynomial drift" of 120242 (−0.09 nm at 550-700) is thus mostly the
  line-set change, not the instrument moving.
- **Not all of it.**  Even without Xe, every Kr epoch keeps about +0.045 nm
  at Hg 404.66 and −0.04 at 546.07.  The 404.66 part is plausibly its blend
  with Hg 407.78, which pulls the centroid red.  So no-Xe is a better
  "truth", but not the truth.

**Pattern rms before and after** (nm, over 27 windows; Ld = lab L):

| group | before | after v1 | after v2 | variance removed v1 / v2 | corr v2 |
|---|---|---|---|---|---|
| 120242 E / Ld | 0.102 / 0.134 | 0.100 / 0.132 | 0.082 / 0.119 | 5 / 37 %, 2 / 21 % | 0.64 / 0.47 |
| 121222 E / Ld | 0.135 / 0.147 | = | 0.120 / 0.137 | 0 / 20 %, 0 / 13 % | 0.46 / 0.36 |
| 121231 E / Ld | 0.153 / 0.128 | = | 0.148 / 0.126 | 0 / 7 %, 0 / 4 % | 0.27 / 0.24 |
| 122302 E / Ld | 0.104 / 0.136 | 0.102 / 0.135 | 0.085 / 0.124 | 4 / 34 %, 2 / 18 % | 0.60 / 0.43 |
| 122303 E / Ld | 0.153 / 0.131 | 0.152 / 0.129 | 0.146 / 0.119 | 2 / 9 %, 3 / 18 % | 0.31 / 0.43 |
| 122304 E / Ld | 0.100 / 0.127 | 0.098 / 0.125 | 0.087 / 0.118 | 4 / 23 %, 2 / 13 % | 0.49 / 0.36 |
| 122305 E / Ld | 0.139 / 0.154 | 0.137 / 0.153 | 0.124 / 0.138 | 2 / 20 %, 2 / 20 % | 0.46 / 0.47 |
| common | 0.119 | 0.117 | 0.106 | 2.5 / 21 % | 0.46 |

**Line by line**, common pattern: field / lab −r (raw line residual, no
interpolation) / v2 about its line, in nm:

| window | field | lab −r | v2 about its line | |
|---|---|---|---|---|
| 405 | +0.035 | −0.087 | −0.057 | **opposite sign** |
| 435 | −0.081 | −0.023 | −0.041 | **lab too small** |
| 545 | +0.108 | +0.093 | +0.037 | agree |
| 585 | +0.053 | +0.077 (Kr instruments) | +0.052 | agree |
| 855 | −0.084 | −0.033 (Kr) | −0.107 | |

- Every instrument has the same signs at 405 and 435; there field and lab
  disagree everywhere.
- At the windows where the pattern is largest, the lab is nowhere near:

  | window | field | v2 |
  |---|---|---|
  | 415 | −0.197 | −0.052 |
  | 505 | −0.226 | +0.012 |
  | 605 | +0.245 | +0.055 |
  | 665 | −0.184 | +0.045 |

  No lab line lies near these windows, and the lab residual is a difference
  of 4th-order polynomials.  It is smooth, and cannot make ±0.2 nm
  window-to-window structure.

**Verdict on Phase 2 open concern #1: partly resolved, and mostly not.**

- The lab polynomial explains a smooth, low-order part of the pattern.  It
  is real, and largest near 545-605 and 855 nm: about 20 % of the variance
  with the no-Xe variant (up to 37 % on 120242 E), and 2-5 % with Q13's
  default.  It also explains why the pattern is common to every instrument:
  all of them share one fitting recipe and line set.
- It does not explain the ±0.2 nm window-to-window structure (415, 505,
  605, 665 nm), which carries ~80 % of the variance.  Nor does it explain
  the sign at 405 nm or the size at 435 nm.
- That structure is still unexplained.  It is not the SRF shape (task 9b),
  not HSRS line positions, not the E grid (task 3: the pattern moved
  ≤ 0.0025 nm), and not the lab polynomial.

**srfcal option: not added.**

- Under Q13's default the lab residual removes 2.5 % of the common
  pattern's variance, which is not significant; under v2, 21 %.
- The empirical pattern already contains that component, so a
  `load_window_pattern(source='lab')` option would replace a correction
  that works with a weaker one.
- `phase2_wavescale.py` was not re-run, and `wavescale_window_pattern.csv`
  is unchanged.  It still has its task-3 `_pre_lab` copy only; no new copy
  was needed.
- The useful product of this task for the processor side is the lead
  itself: a fit without Xe, or a better pixel → λ model, would remove the
  smooth ~0.05-0.1 nm error from the processor grid.  That is for the note
  to Joel.

**Products:**

- new `hypernet/wiggles/lab_window_pattern.py`;
- `lab_window_pattern.csv` (per window × group: field, lab line, v1 and v2
  raw and detrended, after);
- `lab_window_pattern_summary.csv` (the rms table and the at-window values);
- `figs/lab/window_pattern_vs_lab.png`.  Panels: (a) field vs lab;
  (b) the refits; (c) the 120242 transfer test;
- intermediates in `$OS_COLOR/hypernet/wiggles/lab_data/window_pattern/`
  (`kr_refits.csv`, `transfer_curves.csv`, `transfer_test_120242.csv`,
  `run.log`).

**Tests:** no package code changed, so no tests were added.  `pytest -q`:
151 passed.

**Learned:**

- The same-pixel "drift" between lab epochs mixes instrument change with
  line-set change.  On 120242, most of the 2020→2022 polynomial change is
  the switch from a Kr to a Xe lamp.  This qualifies the "≤ 0.10 nm
  polynomial stability" of the Context; the instrument moved less than
  that.
- A line-set refit that matches Joel's files needs no weights:
  `file_minus_fit_max` is ≤ 0.015 nm for the all-line refits.

### 2026-10-18 -- Build #5 (Opus 5.5)

Task 5 (validate the L2 inversion), with Setup Q10's default.  ρw was not
recomputed and no Phase 3 metric was run.  Task 6 was not started.

**Code** (`hypernet/edinterp.py`):

- `interp_matrix(wav_irr, wav_rad)`: sparse M with M E = `np.interp`, built
  from `linear_weights`, the single definition.
- `recover_native_e(irradiance_L, wav_L, wav_E, ..., return_info=False,
  null_tol=1e-6)`: least squares (`lsqr`, scaled) on 1-D or
  (n_L, n_seq) input.  Edge rules:
  - L pixels beyond the E range are clamped rows and stay in the system;
  - NaN rows of E_L are dropped per column;
  - E pixels with a zero column are NaN.
  - **New:** so are E pixels in the null space of the remaining system,
    found by a dense SVD that runs only when rows were dropped or columns
    are empty.  HYPSTAR has n_L ≈ n_E (two of the L rows are the clamped
    ones), so a gap of 4 dropped rows leaves a stretch undetermined.  The
    minimum-norm solution put errors of up to 0.47 there, so these pixels
    are now NaN, and only those.
  - **New:** an E grid running past the last L pixel loses the clamped
    anchor, and ~90 red-end pixels become NaN.  This is correct.
  - `info` returns `n_dropped_rows`, `n_unconstrained` (NaN pixels) and
    `forward_rel_rms`.
- `second_difference_power(E, wav, band=(400, 900), nyquist=False)`: rms of
  the relative second difference.  With `nyquist=True`, d2 minus its
  3-pixel running mean (the "ringing index").

**VEIT recovery** (L1C irradiance → native E with `grid_for('122304','E')`,
against the L1A E, 400-900 nm):

- relative error rms **1.55e-5**, max **5.3e-5**, second-difference error
  0.19 %;
- the same with the L1A grid itself;
- `cond(M)` as before.

**The central finding: the L2 data cannot verify which E grid the processor
used.**

- Solving M_wrong x = M_true E makes the wrong grid's piecewise-linear
  curve match the true one at the L pixels.  So x is the same spectrum
  resampled at the wrong grid's own positions, correctly labelled.
- A 0.05 nm grid error therefore gives neither ringing nor a scale error;
  it changes only the double-interpolation smoothing.
- Against the L1A truth resampled at the wrong positions, it costs 7e-4 to
  1e-3 in E and 24-36 % in the second differences (2024_09 grid: 6.6e-4 and
  17 %).  Phase 3 would feel this, but the recovered spectrum shows nothing.
- Three detectors were tried:

  | detector | result |
  |---|---|
  | second-difference power, plain / Nyquist | moves 1-4 % on VEIT and < 1 % on synthetic spectra (`lab_inversion_synthetic.py`, 400-6000 lines, noise 0-3e-3) |
  | minimum of the ringing index over the grid shifted by δ (`delta_min`) | at the truth on VEIT (δ = 0; the ±0.05 injections found), but on synthetic spectra with a correct grid at −0.06 to −0.08 nm; the `own+0.05` "control" moves it by −0.05 only because the index depends on the absolute grid position, so it is not a control |
  | native minus `via_grid` wavelength scale (`dshift`) | ±0.02 nm for every grid, including the 0.05 nm control |

- What the ringing index does catch is a grossly wrong grid with a
  different pixel-phase pattern.  The 122304 proxy rings 1.4-2.6× more than
  the own grid on 5 of 6 instruments (not on 122303: 0.83×).  Every own grid
  rings at VEIT's level: 0.0067-0.0092, against 0.0077 on VEIT L1C with the
  true grid.
- My first single-sequence reading of `delta_min` (121222 at +0.06, 122302
  at −0.04) looked like evidence against the inferred E epochs.  The
  synthetic test shows that a correct grid produces such offsets, so it is
  not evidence.

**E-epoch verdict** (Release 2, 5 sequences per instrument × period, own
grid; all recovered E positive over 400-900 nm; `sdp` 0.0085-0.0139):

| instrument | period | ring (own) | proxy / own | delta_min (diagnostic) | dshift | verdict |
|---|---|---|---|---|---|---|
| 120242 | 2024-11 | 0.0079 | 1.39 | −0.08 | +0.006 | consistent (gross); ≤ 0.05 nm not testable |
| 121222 | 2021-10 | 0.0071 | 1.84 | +0.06 | −0.001 | consistent (gross); ≤ 0.05 nm not testable |
| 121222 | 2025-03 | 0.0068 | 2.03 | +0.06 | −0.002 | consistent (gross); ≤ 0.05 nm not testable |
| 121231 | 2021-10 | 0.0067 | 2.59 | 0.00 | −0.002 | consistent (gross); ≤ 0.05 nm not testable |
| 122302 | 2023-01 | 0.0084 | 1.65 | −0.04 | +0.005 | consistent (gross); ≤ 0.05 nm not testable |
| 122303 | 2023-01 | 0.0088 | 0.83 | 0.00 | −0.014 | consistent (gross); ≤ 0.05 nm not testable |
| 122304 | 2023-01 | 0.0091 | – | +0.01 | −0.000 | verified (VEIT L1A) |
| 122304 | 2024-11 | 0.0092 | – | +0.01 | +0.007 | verified (VEIT L1A) |
| 122305 | 2023-01 | 0.0086 | 2.19 | +0.01 | −0.004 | consistent (gross); ≤ 0.05 nm not testable |

- The other-epoch controls (120242 2022_06, 121222 `_tst`, 122302 and
  122304 2024_09) are indistinguishable from the own grid on every
  detector, as expected for 0.01-0.05 nm differences.
- The evidence for the E epoch on the six instruments other than 122304
  therefore stays the L-grid match (`inferred_L`).  Only L1A (or the
  processor's own E wavelengths) can settle it.  This belongs with the
  delivery request and the note to Joel/Kevin.

**Widths** (E FWHM at 546 nm, median of 5 sequences): native-grid fit of
the recovered E / `via_grid` fit of the L2 E / lab.  The lab value is the
lab FWHM for Gaussian fits and our super-Gaussian lamp refit for the
super-Gaussian ones.

| instrument | period | Gaussian: native / via / lab | super-Gaussian (lab p): native / via / lab |
|---|---|---|---|
| 120242 | 2024-11 | 2.181 / 2.183 / 2.225 | 2.370 / 2.375 / 2.334 |
| 121222 | 2021-10 | 2.757 / 2.759 / 2.741 | – |
| 121222 | 2025-03 | 2.730 / 2.729 / 2.741 | – |
| 121231 | 2021-10 | 2.871 / 2.871 / 2.807 | – |
| 122302 | 2023-01 | 2.217 / 2.217 / 2.258 | 2.513 / 2.515 / 2.490 |
| 122303 | 2023-01 | 2.715 / 2.726 / 2.727 | – |
| 122304 | 2023-01 | 2.321 / 2.320 / 2.368 | 2.577 / 2.578 / 2.508 |
| 122304 | 2024-11 | 2.353 / 2.353 / 2.373 | 2.621 / 2.623 / 2.513 |
| 122305 | 2023-01 | 2.711 / 2.707 / 2.759 | – |

- Native and `via_grid` agree to ≤ 0.011 nm everywhere.  They carry the
  same information, so the `via_grid` route of task 3 loses nothing.
- Field − lab is −0.05 to +0.06 for the Gaussian rows and −0.02 to +0.11
  for the super-Gaussian rows.  This is the 122304 over-width at the lab p
  seen in task 3.

**Products:**

- new `hypernet/wiggles/lab_inversion_check.py`, writing:
  - `lab_inversion_check.csv` (155 rows: per sequence × grid);
  - `lab_inversion_summary.csv` (per instrument × period × grid, with the
    verdict and the widths);
  - `figs/lab/inversion_check.png`.  Panels: (a) VEIT errors vs the L1A
    truth at 425-445 nm; (b) the Release 2 ringing scans, labelled as a
    diagnostic; (c) the widths;
  - intermediates in `$OS_COLOR/hypernet/wiggles/lab_data/inversion/`
    (`veit.csv`, `ring_scans.csv`, `widths.csv`, `run.log`).
- New helper scripts:
  - `hypernet/wiggles/lab_inversion_synthetic.py`: detector behaviour on
    synthetic spectra;
  - `lab_inversion_curves.py`: prints the median ringing curves.

**Tests** (`test_edinterp.py`, +5 tests; `pytest -q`: **156 passed**):

- Tier 1:
  - `interp_matrix` = `np.interp`, with clamped rows;
  - the round trip M·recover(M·E) = M·E, exactness, and 2-D input;
  - NaN rows and unconstrained / null-space pixels.  Only pixels near the
    gap or the red end become NaN, and the rest is exact to 1e-6;
  - `test_wrong_grid_is_detectable_only_against_the_truth`: a ±0.05 nm
    error gives a > 5 % second-difference error against the truth, but
    moves the plain and Nyquist powers by < 2 %.
  - The prompt's "detected as rising second-difference power" is not
    achievable, so the test pins the honest version.
- Tier 2: `test_recover_native_e_veit` (VEIT ≤ 1e-4 over 400-900 nm; it
  skips without the sample).

**Learned:**

- Inverting a linear resampling with the wrong nodes is self-consistent.
  The data cannot reject a grid that is off by ≤ 0.05 nm; they can reject
  only a grossly different pixel-phase pattern.
- The task's premise ("a wrong grid shows up as excess second-difference
  power") holds only against a truth.
- For Phase 3, the inversion route is as good as its grid, and the grid
  cannot be checked from L2.  This adds weight to waiting for the L1A.

### 2026-10-18 -- Gate #6 (Opus 5.5)

Task 6, the provisional gate G0.  No package code changed, and task 7 was
not started.

- **Read:**
  - `docs/wiggles_planning.md` §4 (Phase 0, Gate G0) and §9 in full;
  - `gate_G1.md` and `phase2_exit.md` for the structure;
  - the Phase 0 Logs 8c-8f;
  - Builds #2-#5 here;
  - the committed tables.
- **Wrote `claude_prompts/wiggles/gate_G0.md`** (not under `docs/`), headed
  "provisional: Release 2 L2 + lab data, no L1A".
- **New `hypernet/wiggles/lab_g0_numbers.py`** prints the numbers not
  already in the Logs:
  - the field FWHM std per group;
  - the SZA and sky-index slopes turned into changes over each group's
    range;
  - the recalibration and site pairs;
  - the lab epoch-to-epoch ranges;
  - the offsets and scale fits.

  It writes nothing.
- **Verdict: provisional pass, with conditions.**
  - **G0(a): pass per instrument.**  H2 is needed on 120242, 122302 and
    122304 (field 0.35-0.49, lab 0.39-0.47 nm) and not on the other four
    (lab −0.06 to +0.02; 122303 set from the lab).  The criterion's "most
    instruments" wording would drop H2: edit 1.
  - **G0(b): pass.**
    - Field max std 0.021-0.059 nm; SZA ≤ 0.08 nm; sky ≤ 0.09 nm, on Ld
      only (Ring-like, 2-4σ).
    - Recalibration ≤ 0.03 nm (122304 E +0.025-0.029, 5-7σ); site
      ≤ 0.012 nm.
    - Lab epochs: FWHM ≤ 0.12 nm, polynomial ≤ 0.10 nm (120242, mostly the
      line-set change).
  - **G0(c): triggered, and handled.**
    - Offsets ≤ 0.10 nm (122304 Ld, 0.09-0.10).
    - The Ld dispersion trend is 0.12-0.24 nm across 390-870 nm, plus the
      ±0.09 nm lab-polynomial residual.
    - Phase 2's shift + stretch correction is the §4 wavelength correction.
- **Settings:** `TABLE_IS_DEFAULT = True` and
  `APPLY_CENTROID_OFFSET = False` are supported by the evidence, and were
  not changed.
- **Open issues listed:**
  1. the lab E shape costs O₂-A on VEIT;
  2. the G1 threshold on Gaussian-equivalent widths;
  3. 122303 data vs lab;
  4. 80 % of the window pattern unexplained;
  5. the withdrawn L2 bias, still in `delta_fwhm_l1a_est`;
  6. the E epoch is untestable from L2 on six instruments;
  7. the items carried from Phase 2.

  The L1A checklist has eight items.
- **Q&A "Gate #6"**: five proposed §4 edits, each quoting the current
  text:
  1. G0(a) per instrument;
  2. G0(b) with lab epochs and a ~0.1 nm jump bound;
  3. G0(c) including the dispersion and the lab grids;
  4. the Phase 0 SRF shape;
  5. a G0 status line.

  `docs/wiggles_planning.md` was not edited.
- **Learned:**
  - The plan's G0 was written as one global switch.  The data make it
    per-instrument at every step: class, shape and epoch.
  - The only criterion that triggers is (c), through the dispersion, which
    Phase 2 already corrects.

### 2026-10-18 -- Report #7 (Opus 5.5)

**A. §4 edits applied.**  JXP approved all five Gate #6 proposals ("Ok use
your Proposed").  I applied them to `docs/wiggles_planning.md` exactly as
proposed, and left the rest of each paragraph unchanged:

1. G0(a) per instrument;
2. G0(b) with lab epochs and the ~0.1 nm recalibration bound;
3. G0(c) with shift/stretch and the `labcal` grids;
4. the super-Gaussian clause in Phase 0's first bullet;
5. the provisional-G0 status line.

The only unpublished reference added under `docs/` is the path
`claude_prompts/wiggles/gate_G0.md`, as proposed.  In `gate_G0.md`, the
last section now records the edits as applied (2026-10-18, JXP's approval),
and the Verdict's forward reference was updated.

**B. Task 7, the interim report.**

- **Kept the previous version** as
  `reports/wiggles_interim_report_2026-10-15.md` (an unchanged copy).
- **Edited `reports/wiggles_interim_report.md`** in place, same voice and
  readers; dated 18 October, "updated from the 15 October version".
  - **Summary.**  E ≈ L tolerance 0.07 → 0.06 nm; the 122303 "borderline"
    wording removed.  New bullets:
    - the lab confirmation of the two classes, and the air scale;
    - the flat-topped E, with the twin's 49-63 % and the VEIT O₂-A caveat;
    - the window-pattern verdict: Xe-driven ±0.09 nm residual, 2.5-21 % of
      the variance, ~80 % unexplained.
  - **§3.1.**  Stability 0.01-0.06 → 0.02-0.06 nm, with a note that the
    own E grids replaced the stand-in (≤ 0.055 nm effect).
    - **Table 1** now has a shape column and the old values in brackets.
    - 122303 moves to E ≈ L / `linear`, set from the lab.
    - The "E from L2 reads 0.06-0.09 nm too wide" note is withdrawn.
  - **§3.3.**
    - 137 → 156 tests;
    - the table carries the E shape;
    - the lab-shape twin row (49-63 %, 27-45 % with ΔFWHM 0.2, 100 % shape-
      aware);
    - window pattern 0.13 → 0.12 nm rms (common pattern; 0.13 was over the
      per-instrument patterns), and not the stand-in grid (≤ 0.0025 nm).
  - **New §3.4** (the lab data), with Figure 8 (`figs/lab/fwhm_lab_vs_table.png`)
    and Figure 9 (`figs/lab/window_pattern_vs_lab.png`).  Table 3 gives lab
    vs field ΔFWHM, the E p, and the processor epoch with its status.  It
    also covers:
    - the grids are the lab polynomials (≤ 4e-5 nm), one epoch per
      instrument, and the 120242 pixel window;
    - the shaped E and its VEIT trade-off;
    - the withdrawn bias;
    - the window pattern;
    - the L2 inversion: 1.5e-5 on VEIT; the E epoch is not testable from L2.
  - **§4.**  The provisional G0 verdict ((a) pass per instrument, (b) pass,
    (c) triggered and handled), "the SRF table is now validated against the
    lab", and the two pending decisions (shaped E default vs O₂-A; the G1
    threshold on Gaussian-equivalent widths).  The Phase 0b row adds the E
    epoch check; Open items adds the ~80 % pattern and O₂-A.
  - **§5.**
    - #1 now notes that the L1A is the only way to check the E epochs.
    - #2: received, with thanks to Joel, plus his five questions (the
      epoch after recalibration, a 2025 epoch for 121222, the nominal `wl`
      grid, the lamp illumination of the diffuser, which lines were fitted)
      and the public-tables note.
    - #4: the lamp floor ≤ 3e-4.
    - #7: answered (air).
    - New #8: the suggestion to fit without the Xe lines, or a better
      pixel → λ model, phrased as a suggestion.
  - **Footer:** adds `labcal.py` and `hypernet/data/lab_*.csv`.
- **Numbers that moved** (15 Oct → 18 Oct):

  | quantity | 15 Oct | 18 Oct |
  |---|---|---|
  | ΔFWHM 122302 | 0.46 | 0.44 |
  | ΔFWHM 120242 | 0.39 | 0.35 |
  | ΔFWHM 122303 | 0.14 (borderline) | 0.17 (E ≈ L / `linear` from the lab) |
  | ΔFWHM 121231 | 0.08 | 0.03 |
  | ΔFWHM 122305 | 0.01 | 0.03 |
  | ΔFWHM 121222 | −0.05 to −0.03 | −0.06 to −0.05 |
  | ΔFWHM 122304 | 0.39-0.49 | 0.39-0.49 (unchanged) |
  | E ≈ L agreement at 450 nm | 0.07 nm | 0.06 nm |
  | stability | 0.01-0.06 nm | 0.02-0.06 nm |
  | window pattern | 0.13 nm | 0.12 nm rms |
  | tests | 137 | 156 |
  | L2 E-wide bias | 0.06-0.09 nm | withdrawn |

- **New `hypernet/wiggles/lab_report_check.py`:**
  - it checks that every relative link resolves from `reports/` (9 targets,
    0 missing), and that Figures 1-9 are in sequence;
  - it checks Table 1 (field ΔFWHM, p) and Table 3 (lab ΔFWHM) against
    `srf_table.csv` and `lab_fwhm_compare_summary.csv`;
  - the first run caught two of my rounding errors (121231 0.04 → 0.03;
    122304 p 2.81 → 2.80) and a range-parsing bug in the checker; all are
    fixed;
  - result: all checks passed.
- **`correspondence/wiggles_note_lab_data.txt` was not edited.**  Its numbers
  still hold.  Its "about 2.2" for the VEIT field exponent is the Phase 0
  weighted mean; task 3's free-p fit gives 2.54, and the report says
  "~2.2-2.5".
- No package code changed, so pytest was not re-run (156 passed after task
  5).

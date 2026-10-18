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

## Logs

### 2026-10-18 -- Doc created

- Written from `wiggles_prompts.md` Q&A "New data #1" and JXP's answers:
  - proceed with items 1-5 in order, and add item 6 (update the interim
    report);
  - Phase 3 waits for the L1A;
  - grids are computed from committed coefficients.
- Items 1-6 map to tasks 2-7, with Setup as task 1.

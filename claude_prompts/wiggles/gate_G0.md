# Gate G0 -- SRFs, stability and wavelength scale

**provisional: Release 2 L2 + lab data, no L1A**

*2026-10-18, lab data prompts task 6 (Opus 5.5).  Not part of the published
docs.*

*Inputs:*
- *Builds #2-#5 of `wiggles_lab_data_prompts.md`;*
- *Phase 0 tasks 8c-8f (`wiggles_phase0_prompts.md`);*
- *`phase2_exit.md`;*
- *the committed tables `hypernet/data/srf_table.csv`,
  `hypernet/wiggles/lab_fwhm_compare_summary.csv`,
  `lab_window_pattern_summary.csv`, `lab_inversion_summary.csv` and
  `phase0c_stability*.csv`.*

*The numbers not in the Logs come from
`python -m hypernet.wiggles.lab_g0_numbers`.  The criteria are those of
`docs/wiggles_planning.md` §4.*

## Verdict

**Provisional pass, with conditions.  G0 stays open until Phase 0b fits the
L1A delivery.**

- **G0(a), the SRF mismatch (H2): pass, per instrument.**
  - H2 is real and large on three instruments (120242, 122302, 122304:
    ΔFWHM 0.35-0.49 nm field, 0.39-0.47 nm lab).  It is absent on four
    (121222, 121231, 122303, 122305).
  - Field and lab agree on every class.  For 122303 that holds only once the
    method is set from the lab.
  - The criterion as written ("on most instruments") would drop H2, because
    4 of 7 are E ≈ L.  It needs the per-instrument wording proposed below.
- **G0(b), stability: pass.**  The field spread is ≤ 0.06 nm.  Recalibration
  jumps are ≤ 0.03 nm and the lab epochs agree to ≤ 0.12 nm, all against
  0.2 nm.  **Phase 2 ships the per-instrument SRF table**, so
  `TABLE_IS_DEFAULT` stays True.
- **G0(c), the wavelength scale: triggered, and already handled.**
  - The centroid offsets are ≤ 0.10 nm, on the air scale.
  - The dispersion is not: shift + stretch reaches 0.1-0.24 nm at the ends
    of 390-870 nm, systematically on every Ld channel.  The lab polynomial
    also misses its own lines by ±0.09 nm.
  - The wavelength correction §4 calls for is Phase 2's per-sequence shift
    + stretch, with the common window pattern removed.  It stays in front
    of everything else.  `APPLY_CENTROID_OFFSET` stays False.
- **What would change it:** an L1A E SRF that differs from the table by
  > 0.07 nm, or a class flip on any instrument.  Also an E grid that is not
  the inferred lab epoch, which is testable only with the L1A: this would
  break the E grids and the inversion route.  Or a stability failure at the
  scan level (low sun, broken cloud), which L2 cannot see.

## Criteria

| # | Criterion (plan §4) | Result | Status |
|---|---|---|---|
| (a) | FWHM_E(λ) = FWHM_L(λ) within the fit uncertainty on most instruments → drop H2 | ΔFWHM = L − E, median 400-680 nm (Gaussian-equivalent): field 0.35 / 0.44 / 0.39-0.49 nm and lab 0.44 / 0.47 / 0.39 nm on 120242 / 122302 / 122304.  On the other four: field −0.06 to +0.04 (122303 +0.17), lab −0.06 to +0.02.  Two classes: in the lab there is a 0.37 nm gap (0.02 → 0.39 nm) around the G1 threshold of 0.15 nm | **Pass per instrument** (H2 needed on 3/7).  The wording is ambiguous |
| (b) | SRFs stable: spread across SZA, sky and season < ~0.2 nm, no jump at recalibration that the cal files do not explain → ship a per-instrument table | Field (224 L2B sequences, own E grids): largest FWHM std 0.021-0.059 nm in every group; over each group's SZA range ≤ 0.08 nm; over its sky-index range ≤ 0.09 nm (Ld only, 2-4σ, Ring-like; E ≤ 0.03, n.s.); season slopes ≤ 0.05 nm/unit.  Recalibration: 122304 E +0.025-0.029 nm (5-7σ), Ld ≤ 0.02; 121222 E −0.028 nm at 400, Ld n.s.  Site (122302 BEFR/THFR) ≤ 0.012 nm.  Lab epochs: FWHM ≤ 0.012 nm (122302, 122304, 2023→2024) and ≤ 0.12 nm (120242 E, 2020→2024); polynomial ≤ 0.04 nm, and ≤ 0.10 nm on 120242, which is mostly a lamp/line-set change (task 4) | **Pass** (no cal files to "explain" jumps, but every jump is < 0.03 nm) |
| (c) | Centroid offsets > ~0.1 nm systematically on some instrument → wavelength recalibration first | Table offsets vs HSRS (air): −0.05 to +0.10 nm (largest: 122304 Ld, 0.09-0.10).  Scale: shift +0.00 to +0.05 nm at 630 nm, stretch +0.03 to +0.11 nm per 240 nm (Ld 121231 0.11); the Phase 0d trend is 0.12-0.24 nm across 390-870 nm on Ld, 0.03-0.12 on E.  The lab polynomial misses its own lines by +0.09 (405) / −0.09 (546) / ~±0.1 nm (NIR) on every instrument | **Triggered** (the dispersion, not the centroid).  Handled by Phase 2's scale correction |

### G0(a) in detail: two classes, one definition issue

- **The lab confirms the field split.**  Lab minus table FWHM at 436/546
  nm, comparing like with like (task 3, `srf_table.csv`
  `lab_minus_table_*`):
  - Gaussian E rows: −0.05 to +0.13 nm;
  - Ld rows: −0.08 to +0.05 nm;
  - super-Gaussian E rows: −0.13 to +0.02 nm.

  Lab and table give the same G1 method on every row, once 122303 is set
  from the lab.
- **122303 is the exception in the data.**
  - With its own E grid, its Gaussian E narrows by 0.02-0.04 nm, and the
    field ΔFWHM rises from 0.138 to 0.173 nm, over the 0.15 nm threshold.
  - The lab ΔFWHM is −0.024, and E ≈ L on every lamp line.
  - The table sets `linear` from the lab (`method_from = lab`, with a
    note).  This is the only data-lab disagreement.
- **The flat-topped E.**
  - On 120242, 122302 and 122304, the lab E lines are super-Gaussian:
    p = 2.53 / 2.93 / 2.76-2.81, with χ²_sg / χ²_g = 0.49 / 0.09 / 0.50
    against Setup Q7's threshold of 0.5.  The L lines are Gaussian on every
    unit.
  - The FWHM of a flat-topped SRF is 0.25-0.30 nm larger than its
    Gaussian-equivalent width.  G1's threshold was set with Gaussian SRFs,
    and the lab's quoted FWHMs behave like Gaussian-fit widths.
  - So the table's ΔFWHM and `method` use the Gaussian-equivalent E width,
    while the convolution uses the shaped SRF.  On shaped widths,
    `delta_fwhm_shaped` would give 0.12-0.24 nm, and would flip 122304 @
    2024-11 to `linear`.
- **The "E from L2 reads 0.06-0.09 nm wide" bias is withdrawn.**  Against
  the lamp, with own E grids, the E-specific excess has a median of
  −0.055 nm over the 9 Release 2 models.  The lab's FWHM definition alone
  explains −0.048 of that.  The +0.064 nm L2-vs-L1A difference is one VEIT
  sequence.

## Per-instrument SRF summary

| instrument | class | ΔFWHM field (Gaussian-equivalent) | ΔFWHM lab | E shape (p) | method | method_from | processor epoch (L / E status) | E epoch testable from L2? |
|---|---|---|---|---|---|---|---|---|
| 120242 | narrow-E | 0.353 (2024-11) | +0.436 (2024_05) | super-Gaussian 2.527 | srf | data | 2024_05, pixel window of 2020_09 (`CROP_EPOCH`) (verified / inferred_L) | no: consistent at the gross level only |
| 121222 | E ≈ L | −0.057 (2021-10), −0.047 (2025-03) | −0.039 (2021_07) | Gaussian | linear | data | 2021_07 for both periods; no lab epoch for 2025-03 (verified / inferred_L) | no |
| 121231 | E ≈ L | 0.035 | +0.018 | Gaussian | linear | data | 2021_07 (verified / inferred_L) | no |
| 122302 | narrow-E | 0.438 | +0.469 | super-Gaussian 2.934 | srf | data | 2023_01; the 2024_09 lab epoch is unused (verified / inferred_L) | no |
| 122303 | E ≈ L | 0.173 | −0.024 | Gaussian | **linear** | **lab** | 2023_01 (verified / inferred_L) | no |
| 122304 | narrow-E | 0.491 (2023-01), 0.390 (2024-11 R2), 0.477 (VEIT L1A) | +0.395 / +0.389 | super-Gaussian 2.757 / 2.805 | srf | data | 2023_01 for both periods, after the 2024-11 recalibration too (verified / **verified**) | **yes** (VEIT L1A) |
| 122305 | E ≈ L | 0.026 | −0.062 | Gaussian | linear | data | 2023_01 (verified / inferred_L) | no |

Every row is `provisional`; `e_grid` = `own` everywhere.  The grids come from
`labcal.grid_for`.  They reproduce VEIT's L1A E and L to ≤ 2.5e-5 nm, and
every Release 2 L grid to ≤ 3.9e-5 nm.

"Testable from L2" (task 5): inverting the L2 irradiance with a grid that is
off by ≤ 0.05 nm is self-consistent.  It returns the same spectrum, resampled
at that grid's positions, with no ringing and no scale signature.  So only
gross grid errors are excluded: the 122304 proxy rings 1.4-2.6× more on 5 of
6 instruments.

## Supporting findings

- **The scale is air**, settled by the lab: all 19 lamp reference lines are
  NIST air wavelengths.
- **Processor grids.**  There is one lab polynomial per instrument, not one
  per calibration period, cropped to 350 ≤ λ < 1100 nm.  On 120242 the
  pixel window is the one its 2020_09 polynomial selects, with 2024_05
  wavelengths (`labcal.CROP_EPOCH`), so its E has 1539 px.
- **The old 122304 proxy E grid hardly mattered.**  Re-fitting with own E
  grids moved the E widths by ≤ 0.055 nm and the E scale by ≤ 0.003 nm, and
  the common window pattern by ≤ 0.0025 nm.  `via_grid` uses the grid only
  as sampling positions.
- **The common polynomial residual** (task 4) is a fitting artefact of
  4th-order polynomials over 297-1084 nm that include the Xe lines.
  - Without Xe, the rms falls from 0.061-0.071 to 0.026-0.029 nm.  120242
    2020_09, which has Kr and no Xe, has σ 0.03.
  - It is still a real error of the processor grid.
  - It explains a smooth part of the field window pattern: 21 % of its
    variance with the no-Xe variant, 2.5 % with Q13's default.
  - It does not explain the ±0.2 nm window-to-window structure, nor the
    405 nm sign or the 435 nm size.
- **Native-grid E from L2**, by inversion, matches VEIT's L1A E to 1.5e-5
  (5e-5 max) over 400-900 nm with the right grid.  Its widths equal the
  `via_grid` route's to ≤ 0.011 nm.

## What it means for Phase 2

- **SRF source: the table.**  G0(b) passes, so `TABLE_IS_DEFAULT = True`
  stays (current value in `hypernet/srftable.py`: True).
  `srfcal.calibrate_srf` stays the alternative.  All rows stay
  `provisional`; `SOURCE_RANK` puts Phase 0b L1A rows first when they
  exist.
- **Wavelength correction: keep it, first.**
  - The per-sequence shift + stretch, with the period median as fallback
    (`scale_or_fallback`) and the empirical common pattern removed
    (`load_window_pattern`).
  - Task 4 found no case for a lab-derived pattern option, so none was
    added.
  - The processor grids are now the lab polynomials (`grid_for`).  No grid
    recalibration beyond the scale fit is supported by the evidence.
- **Centroid offset: not applied.**  `APPLY_CENTROID_OFFSET = False` stays
  (current value: False).  The offsets (≤ 0.10 nm) are the shift term of the
  scale fit, which is applied instead.
- **The method** is G1's rule on the Gaussian-equivalent ΔFWHM:
  - `srf` on 120242, 122302, 122304;
  - `linear` on 121222, 121231, 122305, and on 122303 (from the lab).
- **The SRF shape** is in the table: super-Gaussian E with the lab p on the
  three narrow-E units.  `srf_for` returns it, and `interpolate_ed_to_l`
  convolves with it.  See open issue 1 before treating that as the default.

## Open issues G0 inherits

1. **The lab E shape costs O₂-A on VEIT** (task 3).
   - `srf` with the lab super-Gaussian E (p 2.805) improves the solar lines
     beyond the Gaussian `srf`: Ca K 0.0090→0.0087, G band 0.0088→0.0075,
     Mg b 0.0037→0.0025, Hα 0.0033→0.0029.
   - But it worsens O₂-A (0.0076→0.0198) and O₂-B (0.0028→0.0048).  The
     overall line/away ratio goes from 1.60 to 2.54, and linear is 2.45.
   - The free-p field fit prefers p ≈ 2.4-2.5, against the lab's 2.8.
   - In the twin, a shape-aware correction recovers 100 %, and a Gaussian
     one 49-63 %.
   - The table now carries the lab shape, so the snippet's default
     convolution is shaped.  **Recommend a decision** before Phase 3:
     Gaussian rows until L1A E shapes exist, or a per-band rule.
2. **The G1 threshold for shaped rows.**  It is applied to
   Gaussian-equivalent widths, a choice made in task 3 that JXP should
   confirm.  The threshold was never derived for a shaped E.
3. **122303: data vs lab.**  The field (own grid) gives ΔFWHM 0.173, the lab
   −0.024.  The table follows the lab.  The L1A E width decides.
4. **The common window pattern: about 80 % unexplained** (task 4: partly
   resolved).  The ±0.2 nm per-window structure at 415, 505, 605 and 665 nm
   is not the SRF shape, HSRS line positions, the E grid or the lab
   polynomial.  It stays removed empirically.
5. **The withdrawn L2 bias.**  `delta_fwhm_l1a_est` and so
   `method_sensitive` still use the one-sequence VEIT +0.064 nm.  The lamp
   says the bias is ≈ 0, so whether to drop that column is JXP's call.
6. **The E epoch on six instruments is inferred from L only** (task 5): the
   L2 inversion cannot check it at ≤ 0.05 nm.  This matters most for:
   - 121222 2025-03 (no lab epoch for that recalibration);
   - 122302, and 122304 after 2024 (lab epochs 2024_09 exist but are
     unused);
   - 120242's E crop (1539 px).
7. **From Phase 2 (carried):**
   - H₂O 936 nm with a fixed PWV;
   - Emod's 380-1000 nm range;
   - shape errors invisible to G3's trigger;
   - the processor's `mcsteps = 0`.

## What the L1A delivery must still confirm

- [ ] **E SRF widths on native grids** for all 7 instruments (Phase 0b,
      `srf_models.parquet`), replacing the L2-derived rows (`SOURCE_RANK`
      `phase0b_l1a`).  The L2-derived widths must agree to ≲ 0.07 nm, the
      G1 tolerance on FWHM_E.
- [ ] **The class of 122303** from its L1A E: ΔFWHM < 0.15 nm (lab) or
      0.17 nm (field, L2).
- [ ] **The E wavelength grids** the processor used, on the six instruments
      other than 122304:
  - [ ] 121222 after its 2025-03 recalibration;
  - [ ] 122302 and 122304 after 2024-09 (lab 2024_09 unused?);
  - [ ] 120242's E pixel window (1539 px under `CROP_EPOCH`).
- [ ] **The E SRF shape in sunlight**: free-p fits on L1A E against the lab
      p.  Then whether the table uses the lab p, the field p or a Gaussian,
      which bears on the O₂-A result.
- [ ] **Scan-level stability**: per-scan fits, `ed_cv_750` (broken cloud),
      and low-sun E in the blue (`calibrate_srf`'s low-blue thresholds,
      untested on real data).
- [ ] **Ld = Lu**, the internal check of the L channel, and the Lu/Ld scan
      split.
- [ ] **The cal `.nc` files** (`HYPERNETS_CAL_HYPSTAR_*_{RAD,IRR}_v2.3.nc`):
      what they hold (wavelength coefficients, `bandwidth`), and whether
      they explain the recalibration jumps (122304 E +0.03 nm).
- [ ] **Then** Phase 0 task 14 (the final G0), with Phase 0b's tables
      outranking these.

## Proposed edits to `docs/wiggles_planning.md` §4 (not applied; see Q&A)

See the Q&A entry "Gate #6" in `wiggles_lab_data_prompts.md`.

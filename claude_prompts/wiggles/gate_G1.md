# Gate G1 -- Twin experiment verdict

*2026-10-04, Phase 1 task 12 (Opus 5.5).  Inputs: tasks 9-11 of
`wiggles_phase1_prompts.md`; tables `hypernet/wiggles/phase1_twin_metrics.csv`,
`phase1_controls.csv`, `phase1_degradation.csv`, `phase1_prediction.csv`;
figures in `hypernet/wiggles/figs/phase1/`.  Not part of the published docs.*

## Verdict

**PASS for narrow-E instruments, with two knowledge requirements.  The
correction is not needed, and is risky, on E ≈ L instruments.**

- The `srf` method (generalised eq. 14, with the L SRF in the numerator)
  removes the SRF-mismatch wiggles.  `ruddick2023` (eq. 14) does not, so
  **H2 is a requirement, not a refinement**.
- **Phase 2 default method: `srf`, applied when ΔFWHM = FWHM_L − FWHM_E ≥
  0.15 nm; otherwise `linear`** (the processor's current method).
- It needs the E wavelength scale relative to the model irradiance to
  ≤ 0.05 nm, and the E and L SRF widths to ≲ 0.07 / 0.1 nm.  The SRF
  tolerance is met by Phase 0's per-instrument tables.  The wavelength
  tolerance calls for a per-sequence (or per-period) shift-and-stretch fit
  of E.

## Set-up in one paragraph

High-resolution Ed, Ld and Lw were built from TSIS-1 HSRS × HITRAN O₂/H₂O ×
Serdyuchenko O₃ (0.01 nm) and multiplied by smooth OSOAA fields.  There are 25
scenes: SZA 30/50/70 plus a VEIT-like 40°, AOT(550) 0.05/0.25, and four
water types, with fluorescence and Raman controls in Lw.  The scenes were
observed with Gaussian SRFs from the Release 2 template fits, on the real
HYPSTAR grids, for two instruments: 122304 (narrow-E, ΔFWHM ≈ 0.5 nm) and
122305 (E ≈ L).  The 16 cases (i)-(vii) and five methods were compared by the
rms of ρw'' − ρw''_true (h = 1 and 5 nm) within ±5 nm of the ten lines and
away from them.  Medians over scenes, h = 1 nm, unless stated.

## Criteria

| # | Criterion (plan §4) | Result | Status |
|---|---|---|---|
| 1 | `srf` reduces the line-region ρw'' error > 80 % in (ii)-(iii), to the away-from-lines level | 122304: **~100 %** in (ii) and (iii) (1.45e-4 → 7e-8, numerical floor); 122305: 99.5-99.8 %.  With noise (vii) the line/away ratio falls from 1.70 (linear) to 1.14 (srf) | **Pass** |
| 2 | Degrades gracefully under (iv)-(vi) | (vi) wrong Emod: ≥ 98 % for SZA ±10°, PWV × 0.5-3, low aerosol.  (iv)/(v) on 122304: monotone, never worse than linear.  80 % needs \|δλ_E\| ≤ 0.05 nm (shift; ±0.1 → 62 %, ±0.3 → 12 %), ≈ 0.07 nm stretch, \|δFWHM_E\| ≲ 0.07 nm, \|δFWHM_L\| ≲ 0.1 nm.  On 122305 **not graceful**: 0.3 nm SRF error → 5-7× worse than linear | **Pass (narrow-E), conditional on the tolerances; Fail (E ≈ L)** |
| 3 | Leaves the fluorescence and Raman controls unchanged within the (vii) noise floor | Every method, case and instrument: control ρw'' error ≤ 0.20 of the floor.  In (iii) `srf` leaves both controls exact (1e-13); linear distorts them by 23 % (fl) and 44 % (Raman) of their own ρw'', ruddick2023 by 31 % and 46 % | **Pass** (but see concern 6) |
| 4 | Does `ruddick2023` alone pass (iii) at the VEIT SRF difference? | **No**: 0 % in (ii), −6 % in (iii) on 122304; −128 % on 122305.  It corrects the grid offset (case i: 99 %) but not a width mismatch.  In the ΔFWHM sweep it is −1 to −56 % for every ΔFWHM > 0 | **No → H2 is required** |

## Supporting findings

- **H1 alone is small.**  With equal SRFs (case i), linear interpolation's
  line error is ~4 % of the 122304 mismatch error.  `ruddick2023`, `srf` and
  even `cubic` remove 98-99 % of it.
- **The E ≈ L floor is not H1.**  On 122305 the rms relative ρw error near
  the lines is 0.08 % from the grid offset alone and 0.46 % from offset plus
  the residual (< 0.07 nm) mismatch.  Phase 0 measured 0.5-2 % (Ca H/K
  1.8 %).  The twin reaches the low end, not the high end, so something
  outside the twin is involved: Ring in Ld, SRF shape, or the wavelength
  scale (the metrics differ, so this is indicative).
- **Prediction for G3** (`phase1_prediction.csv`): judge on the observable
  excess √(rms²_line − rms²_away) of ρw'' at h = 1 nm.
  - For ΔFWHM 0.5 nm, a single sequence: **60-75 %** reduction (realistic
    to noise-only knowledge), with a ceiling of ~92 %.
  - For ΔFWHM < 0.1 nm: 0-10 %.
  - Per line, the reduction follows SNR: Ca H/K ~77 %, G band 65 %, O₂-A
    24 %, Hα 13 %, H₂O ~0 at ΔFWHM 0.5.
  - h = 5 nm is dominated by real structure (a perfect correction shows
    28 %), so it is not a useful G3 metric.
- **Why 0.15 nm.**  With perfect knowledge `srf` is never worse than linear.
  With realistic knowledge (0.05 nm wavelength and FWHM errors), it is break-
  even at ΔFWHM 0.1 nm (+1-3 %) and clearly positive from 0.2 nm (+10-21 %).
  Linear's own error crosses the true line structure at ~0.1 nm.

## Open concerns (to carry into Phase 2 and 3)

1. **The E wavelength scale.**  This is the binding requirement:
   ≤ 0.05 nm relative to Emod.  Phase 0 found an air scale within 0.07 nm
   but a dispersion trend reaching > 0.1 nm at the ends.  Phase 2 should fit
   a shift and stretch of E against the model (HSRS) per sequence or per
   calibration period, and propagate its uncertainty.  This ties to G0's
   wavelength-recalibration decision.
2. **E SRFs from L1A.**  The tolerance on FWHM_E (≲ 0.07 nm) is within
   Phase 0's stability (0.01-0.06 nm).  But the Release 2 E widths come
   from L2 irradiance, which reads ~0.1 nm wide (Phase 0 8c), and use a proxy
   E grid except 122304@2024-11.  The L1A delivery is needed before the
   per-instrument E SRF table is trusted, and before ΔFWHM is compared with
   the 0.15 nm threshold.
3. **SRF shape.**  The twin's truth and correction are both Gaussian, so a
   shape error (Phase 0's pvoigt/supergauss fits) is untested.  Suggest a
   Phase 2 twin check with a non-Gaussian truth and a Gaussian correction.
4. **Physics outside the twin.**  It omits the Ring effect in Ld (seen in
   Phase 0), the multiplicative Lu NIR excess (ρ 0.038 vs 0.027), and
   angular effects.  These can leave a residual that correlates with
   Ed''/Ed after correction, which is G3's trigger to revisit the SRF model.
5. **E ≈ L instruments.**  `srf` gains nothing measurable and is fragile
   there.  The 0.15 nm threshold handles this, but the residual 0.5-2 % floor
   on these instruments stays unexplained (see above).
6. **The control criterion is weak as written.**  The controls' own ρw''
   is below the single-sequence noise floor (fluorescence 0.17, Raman 0.05
   of it at h = 1 nm), so any method passes.  Proposed rewording below.
   The sharper measure is how much a method distorts a control relative to
   its own size: `srf` leaves it exact (iii) or no worse than linear (iv-vi).
7. **Single-sequence noise.**  With case (vii) noise, `srf`'s twin
   reduction at the VEIT ΔFWHM is 45-51 %, not ~100 %.  G3 should use the
   noise row of the prediction table and allow for multi-sequence averages.

## Proposed edits to `docs/wiggles_planning.md` §4 (not applied; see Q&A)

See the Q&A entry "Gate #12" in `wiggles_phase1_prompts.md`.

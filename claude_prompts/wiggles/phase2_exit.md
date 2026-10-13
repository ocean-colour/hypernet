# Phase 2 exit check

*2026-10-12, Phase 2 task 10 (Opus 5.5).  Not part of the published docs.
Inputs: `wiggles_phase2_prompts.md` tasks 1-9b and their Logs; `gate_G1.md`.*

## Verdict

**Phase 2 is complete against plan §7 items 1, 2 and 4.**

- All deliverables exist and are tested: `pytest -q` gives **137 passed, 0
  skipped** (2026-10-12, with `$OS_COLOR` mounted).
- G1 is applied.  The G0 parts stay **provisional**, because G0 waits on the
  L1A delivery.
- Three open concerns (below) carry into Phase 3.

## Deliverables (plan §7)

| §7 item | Deliverable | Module | Tests |
|---|---|---|---|
| 1 | three-method interpolation `interpolate_ed_to_l`, `choose_method` (G1 rule), `delta_fwhm` | `hypernet/edinterp.py` | `test_edinterp.py` (23, with item 4) |
| 1 | pixel convolution with a cached sparse matrix (moved from `twin`) | `hypernet/srf.py` | `test_srf.py` (25), `test_edinterp.py` |
| 1 | SRF self-calibration `calibrate_srf` (HSRS template), low-blue flag | `hypernet/srfcal.py` | `test_srfcal.py` (12) |
| 1 | wavelength-scale fit `fit_wavelength_scale`, `correct_wavelengths`, `scale_or_fallback`, common window pattern | `hypernet/srfcal.py`, `hypernet/data/wavescale_window_pattern.csv` | `test_srfcal.py` |
| 1 | Emod builder and cache `emod_for`, `emod_on_grid`, `emod_perturbations` | `hypernet/emod.py` | `test_emod.py` (11) |
| 1 | OSOAA wrapper (Phase 1) | `hypernet/rt/osoaa.py` | `test_osoaa.py` (6) |
| 2 | drop-in `InterpolationWavSRF` (punpy `MeasurementFunction`), `add_srf_variables` | `hypernet/snippets/interpolate_wav_srf.py` | `test_snippets.py` (5) |
| 2 | per-instrument SRF table, `srf_for` / `scale_for` / `method_for` | `hypernet/data/srf_table.csv`, `hypernet/srftable.py`; built by `hypernet/wiggles/phase2_srf_table.py` | `test_srftable.py` (6) |
| 4 | punpy propagation into the three L1C components, `propagate_ed_to_l` | `hypernet/edinterp.py` | `test_edinterp.py` |

**The machine-precision special cases hold**
(`test_constant_emod_equals_linear`, `test_srf_with_equal_srfs_is_eq14`,
rtol 1e-12/1e-14):

- with a constant Emod, every method is `np.interp`;
- with SRF_L = SRF_E, `srf` is eq. 14, pixel by pixel.

## What the snippet ships with

- **Method:** gate G1's rule.  It uses `srf` where the median FWHM_Ld −
  FWHM_E over 400-680 nm is ≥ 0.15 nm, and `linear` otherwise.
  - `srf`: 120242, 122302, 122304.
  - `linear`: 121222, 121231, 122305, and 122303 (0.138 nm; 0.20 with the
    estimated L2 bias corrected, flagged `method_sensitive`).
- **SRF source:** the SRF table (`TABLE_IS_DEFAULT = True`).  The sources
  are the Release 2 clear-sky template fits, plus the VEIT L1A model,
  preferred for 122304 @ 2024-11.  `calibrate_srf` per sequence is the
  alternative.
- **Wavelength scale:** a shift + stretch correction of E and Ld.  It is
  per sequence where the fit passes 0.03 nm, else the calibration-period
  median from the table, with the common window pattern removed.  The SRF
  centroid offset is *not* applied (`APPLY_CENTROID_OFFSET = False`).
- **Emod:** `emod_for(mean SZA)`, in 5° SZA bins with PWV 15 mm.
- **Unwired** (the processor's 3-argument `setup`): it falls back to linear
  interpolation and reproduces `InterpolationWavLinear` exactly.

## Results that carry forward

- **VEIT (task 9).**
  - `linear` reproduces the L1C irradiance to a median of 1.3e-5 (plan
    §2.3: 2e-5).
  - On Ld/Ed_L, `srf` cuts the line-region excess from 2.45× to 1.60× the
    away-from-lines level, and Ca H/K about 4×.  `ruddick2023` makes it
    worse (3.13×).
  - New uncertainty terms for `srf`: median 0.036 %, up to 3.4 % in O₂-A,
    against 0.13 % measurement noise.
- **Wavelength scale (task 4b).**
  - The E error reaches ~0.1 nm at the ends, above G1's 0.05 nm, so the
    correction is needed.
  - It is measurable to 0.02-0.05 nm per sequence (≲ 0.01 nm per
    calibration period) once a **0.13 nm common per-window pattern** is
    removed.
- **Non-Gaussian SRFs (task 9b).**  `srf` keeps:
  - 97-98 % for a pseudo-Voigt;
  - 80-90 % for an asymmetric SRF;
  - 76-84 % (Phase 0 weighted mean) or 66-78 % (median) for a flat-topped
    super-Gaussian;
  - 88 % / 77 % for an L wavelength error of 0.05 / 0.1 nm.

  The fitted-FWHM pipeline does better than the true FWHM.

## Provisional while G0 is open

- Every SRF-table row is `provisional`: Release 2 L2 fits.  E from L2 reads
  ~0.06-0.09 nm wide, and the E grids are a 122304 proxy.
- The G0 switches are set from the preliminary findings:
  `TABLE_IS_DEFAULT` (G0(b)) and `APPLY_CENTROID_OFFSET` (G0(c)).
- 122303's and 121231's methods are to be decided from their L1A E widths.
- `calibrate_srf`'s low-blue thresholds (S/N 100, 0.05 nm) are untested on
  real low-sun E.
- When the delivery is fitted, rerun `phase2_srf_table.py` from Phase 0b's
  `srf_models.parquet` (source `phase0b_l1a`, ranked first).

## Open concerns (into Phase 3)

1. **The common window pattern** (0.13 nm rms, the same on every
   instrument and channel).
   - It is not SRF shape, and not HSRS vs Kurucz line positions.
   - Telluric lines explain part of it in the twin, but not on real data.
   - It is removed empirically.
   - Phase 3 should check whether the `srf` residuals concentrate at 415,
     505, 565, 605 and 665 nm.
2. **The H₂O band (936 nm)** gets worse with `srf` on VEIT.  Emod's water
   vapour is a fixed 15 mm; a per-sequence PWV fit, or masking the band, is
   needed.
3. **Emod covers 380-1000 nm.**  The L grid spans 350-1100 nm, so ~290 edge
   pixels fall back to linear.  Extending it needs a HITRAN fetch beyond
   10000 cm⁻¹.
4. **Shape errors are invisible to G3's trigger:** the srf residual hardly
   correlates with Ed''/Ed (|r| < 0.2).  Judge G3 against task 9b's
   super-Gaussian rows as well as the Gaussian prediction.
5. **The processor:**
   - its water default `mcsteps = 0` means none of the uncertainty terms
     propagate in operations;
   - its `InterpolationWavClearSky` (an eq.-14-style interpolator with a
     3 nm placeholder SRF) returns values on the irradiance pixels, not the
     radiance pixels;
   - both are for the next note to Kevin (Build #2 Q2).

## Proposed edits to `docs/wiggles_planning.md`

See the Q&A entry "Gate #10" in `wiggles_phase2_prompts.md`.  None have been
applied.

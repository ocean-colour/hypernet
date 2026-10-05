# Removing spectral wiggles in WATERHYPERNET water reflectance: a plan

**Draft — for discussion between J. X. Prochaska and K. Ruddick** · 28 September 2026

*This page is a working plan, not a result. It records what we have agreed so
far, what one sequence of data already says, what we propose to do, and what
we need from RBINS to do it. It is meant to seed the technical note that will
accompany the code.*

## 1. Summary

HYPSTAR measures downwelling irradiance (E) and radiance (L) through different
optical paths on one spectrometer. Their wavelength grids differ (1536 vs 1538
pixels, ~0.48 nm spacing, a sub-pixel offset that drifts through a full pixel
across the spectrum), so the processor interpolates E linearly onto the L grid
before forming ρw = π Lw / Ed. Kevin's hypothesis (Ruddick et al. 2023) is that
this linear interpolation, applied to irradiance with unresolved Fraunhofer and
telluric structure, injects wiggles into ρw with a large, spurious second
derivative — the quantity pigment-detection algorithms depend on.

Working on the one VEIT sequence Kevin sent (2026-06-04, HYPSTAR 122304), we
confirmed the symptom and found a second, apparently larger, cause: **the E and
L channels have different spectral response functions (SRFs)**. Gaussian fits to
the same solar lines give E about 0.5–0.9 nm sharper than L in the blue,
converging by 700 nm. A resolution mismatch cannot be removed by any
interpolation on the E grid, however clever, but it *can* be removed by a
one-symbol generalisation of Kevin's eq. (14): evaluate the model irradiance
with the **L** SRF in the numerator. That generalisation contains eq. (14)
and linear interpolation as special cases, so all three methods run on one
code path and are tested together.

We therefore propose two hypotheses, tested side by side:

- **H1 (Ruddick et al. 2023):** grid offset + unresolved structure + linear
  interpolation.
- **H2 (this work):** E-vs-L SRF-width mismatch.

The plan has four phases with explicit gates: (0) SRF and wavelength diagnosis
over ~150 sequences; (1) a twin experiment with TSIS-1 HSRS × HITRAN/HAPI gas
transmittance on a 0.01 nm grid and OSOAA for the smooth radiative-transfer
fields; (2) code in `hypernet/` with a drop-in snippet for `hypernets_processor`
and FIDUCEO-style uncertainty via `punpy`; (3) an in-situ test on the requested
data, judged by ρw'' at and away from the lines. Deliverable is code plus a
technical note (option (c) in our earlier exchange).

**What we need from Kevin:** L1A_IRR, L1A_RAD, L1C_ALL and L2A_REF for the
sequences listed in {download}`wiggles_data_request.csv <wiggles_data_request.csv>`
(150 primary + 74 spares over five sites and seven instruments), the RAD and IRR calibration files for those seven
instruments (both versions for the two recalibrated ones), any laboratory
line-spread or SRF data, and the processor version that produced Release 2.
Phases 0 (on the VEIT sample) and 1 (no data needed) start before the data
arrive.

## 2. The problem and two hypotheses

### 2.1 H1: grid offset and unresolved structure (Ruddick et al. 2023)

The SPIE paper explains the mechanism with a square absorption line narrower
than the SRF, equal 3 nm E and L SRFs on 1 nm grids, and centre offsets
λ'_E ≠ λ'_L of 0–0.75 nm. The measured E and L are trapezoidal smearings of the
line, sampled at different phases; linear interpolation of E to the L centres
then produces a positive and a negative lobe either side of the line. Eqs.
(14)–(15) propose weighting the interpolation by a modelled high-resolution
irradiance convolved with the E SRF. The paper does not test them on data:
§3.2 sets the model equal to the measurement envelope, which makes the result
exact by construction.

### 2.2 H2: E and L do not have the same SRF

If the E channel resolves a line more deeply than the L channel, the ratio
L/E retains a residual line profile whose amplitude scales roughly with line
depth × (FWHM_L² − FWHM_E²)/FWHM². This is not a sampling error and is
insensitive to how E is interpolated. Kevin's paper lists SRF differences
among the wiggle sources (§1.5, Table 1) but models them only as centre
offsets. JXP's reading is that E and L travel different paths through the
spectrometer, at different angles, so both the wavelength scale and the slit
illumination can differ; the second changes the line-spread function.

H1 and H2 are not exclusive. The method in §3 treats both.

### 2.3 What one VEIT sequence says

Sequence SEQ20260604T084543 at VEIT_H, HYPSTAR 122304 (post-recalibration),
files supplied by Kevin; scripts `wavecal/sanity_checks_veit.py` and
`wavecal/line_fits_veit.py`. **This is one sequence; read it as a lead, not a
result.**

- **The processor does what Kevin says.** L1C `irradiance` equals
  `np.interp(wL, wE, mean L1A_IRR)` to a median relative difference of
  2 × 10⁻⁵.
- **Grids.** E has 1536 pixels, L 1538; spacing 0.465–0.498 nm. The position
  of an L pixel within its bracketing E pixels spans the full 0.00–1.00 range
  across the spectrum, so the grid offset is not a constant phase. At ~6
  samples per FWHM the spectra are well sampled — the opposite of the
  undersampled case of Chance et al. (2005). The `bandwidth` variable is a
  flat 3.0 nm at every pixel, i.e. a placeholder, not a measured SRF.
- **Linear interpolation error is small.** Linear minus cubic-spline
  interpolation of Ed onto the L grid has rms 3–6 × 10⁻⁴ (relative) over
  400–700 nm, peaking at 0.6 % at Ca H/K and 1.8 % at O2-A.
- **Residual line structure in Ld/Ed is ten times larger.** With the same
  high-pass metric (rms of the departure from a 4-px Gaussian smooth), Ed
  itself has 1.6 × 10⁻², Ld/Ed with linear interpolation 4.7 × 10⁻³. So about
  **30 %** of the Fraunhofer structure survives the ratio, whereas the
  interpolation error accounts for roughly **3 %** of it. A rigid shift of the
  E grid does not help: the best shift is +0.1 nm and gains 3 %.
- **The E and L SRFs differ.** Gaussian line fits (FWHM, nm):

| line | λ (nm) | FWHM E | FWHM Lu | FWHM Ld | centroid Lu − E (nm) |
|---|---:|---:|---:|---:|---:|
| Ca K | 393.4 | 1.83 | 2.47 | 2.46 | +0.20 |
| Ca H | 396.9 | 1.91 | 2.81 | 2.73 | −0.14 |
| G band | 430.8 | 2.52 | 3.26 | 3.13 | +0.03 |
| Hβ | 486.1 | 2.97 | 3.44 | 3.47 | +0.05 |
| Mg b | 517.3 | 3.06 | 3.65 | 3.51 | −0.04 |
| Na D | 589.3 | 2.63 | 3.11 | 3.14 | 0.00 |
| Hα | 656.3 | 2.77 | 2.96 | 2.94 | +0.04 |
| O2-B | 687.0 | 3.05 | 3.05 | 3.14 | +0.08 |
| O2-A | 760.6 | 3.40 | 3.46 | 3.42 | +0.06 |

  E is sharper than L by 0.5–0.9 nm in the blue and the two converge by
  ~690 nm. Lu matches Ld throughout, as expected for one optical path. (Ca H/K
  is a blended pair and the 936 nm H₂O band is broad, so their centroids are
  less reliable; the 940 nm row is omitted.) *Correction (2026-09-28):* the
  first version of this table had the Ld and Lu columns swapped, because
  `viewing_zenith_angle` is measured from nadir (vza < 90 is the water view).
  The headers above are now right; the conclusions are unchanged. A joint
  Ca H/K fit (`hypernet/wiggles/phase0a_veit.py`) gives wider Ca H/K widths (E 2.7,
  L 3.0–3.4 nm), which shrinks the E–L difference there to 0.25–0.6 nm; the
  other entries reproduce to ≤ 0.05 nm except Hα in Lu (+0.12 nm) and the
  O2-A band (±0.3 nm), whose width depends on the fit weights.
- **The wavelength scales already agree.** At the clean lines the E–L
  centroid difference is 0.03–0.08 nm (0.1–0.2 px), and the absolute scale is
  within ~0.05 nm of the laboratory wavelengths at Hα and Na D. These offsets
  are much smaller than the 0.25–0.75 nm the paper simulates, so H1 alone
  predicts small wiggles here — consistent with the 3 % above.
- **Kevin's symptom is real.** Over 400–700 nm, corr(ρw'', Ed''/Ed) = −0.56,
  while corr(ρw'', Lu''/Lu) = 0.03. The rms of ρw'' (1.9 × 10⁻⁴) is comparable
  to rms(ρw · Ed''/Ed) (2.6 × 10⁻⁴): the irradiance structure is the right
  size to be the wiggle.

## 3. Method: one function, three methods

Notation: E_i is the calibrated irradiance at E pixel centre λ_E^i, ω_E^i(λ)
its SRF, λ_L^j an L pixel centre with λ_E^i ≤ λ_L^j < λ_E^{i+1}, and
w_j = (λ_L^j − λ_E^i)/(λ_E^{i+1} − λ_E^i) the linear weight. Emod(λ) is a
modelled irradiance at sub-SRF resolution.

**Linear interpolation** (paper eqs. 11–12; `interpolate_wav_linear.py`):

```{math}
E_{\rm LIN}(\lambda_L^j) = (1 - w_j)\,E_i + w_j\,E_{i+1}
```

**Model-adjusted interpolation** (paper eqs. 14–15), with the model convolved
with the **E** SRF on both sides:

```{math}
E_{\rm MINT}(\lambda_L^j) = (1 - w_j)\,\frac{E^{\rm mod}_E(\lambda_L^j)}{E^{\rm mod}_E(\lambda_E^i)}\,E_i
  + w_j\,\frac{E^{\rm mod}_E(\lambda_L^j)}{E^{\rm mod}_E(\lambda_E^{i+1})}\,E_{i+1},
\qquad
E^{\rm mod}_E(\lambda_c) = \frac{\int \omega_E(\lambda;\lambda_c)\,E^{\rm mod}(\lambda)\,d\lambda}{\int \omega_E(\lambda;\lambda_c)\,d\lambda}
```

**Generalised form** (this work): the numerator is evaluated with the **L**
SRF at the L pixel, so the output has the resolution of the channel it will be
divided into:

```{math}
E_{d,L}(\lambda_L^j) = E^{\rm mod}_L(\lambda_L^j)\left[(1 - w_j)\,\frac{E_i}{E^{\rm mod}_E(\lambda_E^i)}
  + w_j\,\frac{E_{i+1}}{E^{\rm mod}_E(\lambda_E^{i+1})}\right],
\qquad
E^{\rm mod}_L(\lambda_c) = \frac{\int \omega_L(\lambda;\lambda_c)\,E^{\rm mod}(\lambda)\,d\lambda}{\int \omega_L(\lambda;\lambda_c)\,d\lambda}
```

The bracket is the measured-over-model ratio, which is smooth if the model is
good (it carries calibration, aerosol and cloud effects, all broad), linearly
interpolated; the prefactor restores the line structure at L resolution.
Special cases:

- ω_L = ω_E → Kevin's eq. (14) exactly.
- Emod constant → linear interpolation exactly.

So one function, `method = {"linear", "ruddick2023", "srf"}`, covers all
three, and every comparison in Phases 1 and 3 is a like-for-like switch on one
code path. The properties listed after eq. (15) in the paper carry over: if the
measurements follow the model shape the output is α·Emod_L; if they equal it,
the output is Emod_L.

What the generalised form needs beyond eq. (14): the L SRF as well as the E
SRF (Phase 0), and a model irradiance whose *line shapes* are right, since it
now sets the sub-pixel structure of the output. Its broad-band accuracy matters
much less, because the ratio in the bracket absorbs it.

## 4. Phases and gates

### Phase 0 — Diagnosis (needs L1A; starts on the VEIT sample)

*Goal:* per-channel, per-instrument SRFs and wavelength offsets, and their
stability, from the data themselves.

- Parameterise each channel's SRF as a Gaussian whose FWHM is a quadratic in
  λ, fitted to ~10 clean solar/telluric lines per spectrum (the lines in the
  table above plus Ca H/K as a blended pair). Fit E from L1A_IRR and L from the
  Ld and Lu scans of L1A_RAD separately; Ld = Lu is the internal check.
- Line depths in Ld are affected by rotational-Raman filling-in (the Ring
  effect); widths much less so. We fit widths and centroids and treat depths
  as diagnostic only.
- Stability: FWHM(λ) and centroid offset vs SZA, sky condition, season and
  time, and across the recalibrations of 122304 and 121222.
  Same-instrument-two-sites (122302 at BEFR then THFR) and two-instruments-one-
  site (VEIT, GAIT, MAFR) separate site from instrument effects.
- Error budget: the ρw'' wiggle power at each line decomposed into the H1
  (interpolation, from linear-minus-template on the same E) and H2
  (SRF mismatch, from the FWHM difference) contributions, per sequence.
- Compare the fitted SRFs to whatever the calibration files and any lab
  line-spread data contain.

*Gate G0:* (a) If FWHM_E(λ) = FWHM_L(λ) within the fit uncertainty on most
instruments, H2 is dropped and eq. (14) suffices. (b) If the SRFs are stable
per instrument (spread across SZA, sky and season below ~0.2 nm, and no jump
at recalibration that the cal files do not explain), Phase 2 ships a per-
instrument SRF table. If not, Phase 2 must self-calibrate the SRF from each
sequence's own Fraunhofer lines, and the plan's cost rises. (c) If centroid
offsets exceed ~0.1 nm systematically on some instrument, a wavelength
recalibration goes into the function before anything else.

### Phase 1 — Twin experiment (no data dependency)

*Goal:* know, before touching real data, which method removes which wiggle,
and what it does to real features.

- **Task 0: build OSOAA** in JXP's fork (`cp gen/Makefile_OSOAA.gfortran
  Makefile; make all`; gfortran 16 may need `-std=legacy`), set `OSOAA_ROOT`,
  smoke-test with the fork's `tests/test_water_Ed.py`. Write
  `hypernet/rt/osoaa.py`: a merged wrapper (keyword dict →
  `OSOAA_MAIN.exe -KEY value`), the existing `Flux.txt` (direct/diffuse/total
  Ed per level) and `LUM_vsVZA.txt` (upward radiance) parsers, and a **new**
  `LUM_Advanced_Down.txt` parser for Ld. **Check `OSOAA.View.Level` in
  `src/OSOAA_MAIN.F` first**: the fork's docs say ±1, the working notebooks use
  1 = TOA and 4 = 0−, and we need 0+ (3 in the original manual).
- **High-resolution irradiance:** TSIS-1 HSRS F₀ × HAPI (HITRAN2020) O₂ and
  H₂O line-by-line transmittance × O₃ Chappuis/Huggins cross-sections, on a
  0.01 nm grid over 380–1000 nm, with separate air masses for the direct and
  diffuse paths.
- **Smooth fields from OSOAA on a 5 nm grid:** direct and diffuse Ed at 0+,
  Ld at the HYPSTAR sky geometry, Lw at 0+, for a set of SZA (30, 50, 70°),
  aerosol loads and a chlorophyll/sediment grid. These multiply the
  high-resolution spectrum; OSOAA has no gas absorption and no Ring effect, so
  the split is clean.
- **Ring and inelastic scattering:** as agreed, the twin experiment ignores
  the Ring effect in Ld (it is a risk item, §9). Chlorophyll fluorescence (a
  ~25 nm-wide peak near 683 nm) and water Raman are **controls**: they are
  added to the Lw library, and a method fails if it changes their second
  derivative.
- **Instrument model:** Gaussian SRFs with the VEIT FWHM(λ) curves for E and
  L, the real 1536/1538-pixel grids, and cases: (i) grid offset only, equal
  SRFs (Kevin's case); (ii) SRF-width mismatch, no offset; (iii) both; (iv)
  wavelength-calibration error of ±0.1 and ±0.3 nm; (v) wrong SRF width in
  the correction (±0.3 nm); (vi) wrong Emod (SZA off by 10°, water vapour ×2,
  no aerosol); (vii) photon noise at HYPSTAR levels.
- **Methods:** linear, cubic/sinc (as a null), `ruddick2023`, `srf`.
- **Metric:** rms of ρw'' − ρw''_true within ±5 nm of the lines and away
  from them, with h = 1 nm and h = 5 nm second differences as in the paper.

*Gate G1:* the `srf` method reduces the line-region ρw'' error to the
away-from-lines level (target > 80 % reduction in cases (ii)–(iii)), degrades
gracefully under (iv)–(vi), and changes the fluorescence and Raman controls'
ρw'' by no more than the linear method does, and by ≪ the controls themselves
in case (iii) (the controls' ρw'' is itself below the single-sequence noise).
If `ruddick2023` already
achieves this in case (iii) at the VEIT SRF difference, H2 is a refinement
rather than a requirement and Phase 2 defaults to eq. (14).

**G1 outcome (2026-10-04):** passed for narrow-E instruments. `srf` removes
~100 % of the mismatch error in (ii)–(iii); `ruddick2023` removes 0 %, so H2 is
required. Tolerances for 80 %: E wavelength scale ≤ 0.05 nm relative to Emod,
FWHM_E ≲ 0.07 nm, FWHM_L ≲ 0.1 nm; Emod errors are harmless. On E ≈ L
instruments the correction gains nothing and is fragile. See
`claude_prompts/wiggles/gate_G1.md`.

### Phase 2 — Code

- **Default method (G1):** `srf`, applied where the measured ΔFWHM ≥ 0.15 nm;
  `linear` otherwise.
- **E wavelength scale (G1):** fit a shift and linear stretch of E against the
  HSRS-based Emod (per sequence or per calibration period) to ≤ 0.05 nm, and
  carry its uncertainty.
- `hypernet/wavecal/` (module name to be settled): the SRF model, the
  Fraunhofer-line self-calibration (centroids and FWHM(λ) from a spectrum and
  a line list), the Emod builder (HSRS × HAPI × O₃ for a given SZA and
  pressure), and the interpolation function, roughly

  ```python
  def interpolate_ed_to_l(wav_irr, irradiance, wav_rad, *, emod,
                          srf_irr, srf_rad=None, method="srf",
                          u_rel_random=None, u_rel_systematic=None):
      """Return Ed on the L grid, with u_rel_random / u_rel_systematic."""
  ```

- **Uncertainty:** `punpy` Monte-Carlo propagation of the SRF parameters
  (centroid and FWHM, with their fit covariances), the Emod inputs (F₀,
  column amounts, SZA) and the wavelength grids into `u_rel_random` and
  `u_rel_systematic` on the output, in the `obsarray` conventions the
  processor uses.
- **Drop-in snippet:** a measurement-function class in the style of
  `hypernets_processor/interpolation/measurement_functions/interpolate_wav_linear.py`
  (same `function` / argument-name interface), with the SRF table and Emod
  cache as inputs. We will not run `hypernets_processor` end to end.
- Tests under `hypernet/tests/`: the two special cases must reproduce
  `np.interp` and eq. (14) to machine precision on synthetic data.

### Phase 3 — In-situ test

On every requested sequence: recompute ρw with the three methods from L1A
(using the L1C/L2A Lw and skyglint terms unchanged) and compare.

*Metrics:* rms ρw'' within ±5 nm of the ten lines vs away from them, before
and after; corr(ρw'', Ed''/Ed), which should go to ~0; the O2-A and Hα
regions in particular. *Not erasing real features:* the 683 nm fluorescence
peak at the clear and turbid sites and the pure-water absorption shoulders
above 600 nm must be unchanged away from the lines; the spectrum must not
become smoother than the away-from-lines floor. A cross-system control: at
VEIT the PANTHYR system (VEIT_P, separate TriOS RAMSES E and L radiometers, ~10 nm FWHM)
gives an L2 ρw'' with different wiggle physics; no data are requested for it.

*Gate G3:* the line-region ρw'' excess is reduced by the amount Phase 1
predicts for the measured SRF difference (the noise row of
`hypernet/wiggles/phase1_prediction.csv`, observable excess
√(rms²_line − rms²_away) of ρw'' at h = 1 nm; e.g. 60–75 % at ΔFWHM 0.5 nm,
≈ 0 below 0.1 nm), on most sequences and instruments,
without a detectable change in the control regions. If the residual after
correction still correlates with Ed''/Ed, the SRF model (Gaussian, quadratic
FWHM) is revisited before the tech note is written.

## 5. Data request

The candidate list is {download}`wiggles_data_request.csv <wiggles_data_request.csv>`
(written by `wavecal/select_sequences.py`; columns `row, priority, optional,
site, instrument, cal_dates_rad_irr, sequence_time, azimuth, sza, sky, file`).
Sequences were drawn from the Release 2 L2 archive within each
site × instrument × calibration period, spread over SZA tertiles (SZA ≤ 75°)
and month; the instrument id and SZA were read from each file. Each row has N
primary picks plus N/2 spares. **We ask Kevin to swap in spares where a
primary is unsuitable.** The `sky` column is filled by us, not Kevin, from
the sky index Ld(750)/Ed(750) that Kevin recommended (Ruddick et al. 2006,
eqs. 23–24): below 0.05 sr⁻¹ is clear sky, and higher values mean cloud in
the sky-viewing or sun direction. It can be computed from the Release 2 files
we already hold. Over the 224 candidates, the 10th, 50th and 90th
percentiles are 0.013, 0.024 and 0.114 sr⁻¹, and 78 % are clear. The
clear-sky value does not vary with SZA, and every SZA tertile includes
cloudy cases. Partly cloudy skies, which the 0.05 switch does not resolve,
are flagged from the scan-to-scan variability of Ed in L1A.

| row | site | instrument | cal (RAD/IRR) | primary + spare | SZA | months | rationale |
|---|---|---|---|---:|---:|---:|---|
| VEIT 122304 pre-recal | VEIT_H | 122304 | 2023-01-27/25 | 20 + 10 | 25–73° | 9 | Kevin's site, clear water, AERONET-OC alongside; before recalibration |
| VEIT 122305 | VEIT_H | 122305 | 2023-01-27/25 | 20 + 10 | 29–73° | 12 | second instrument at the same site, 2024-06 to 2025-07 |
| VEIT 122304 post-recal | VEIT_H | 122304 | 2024-11-08/12 | 20 + 10 | 23–72° | 12 | same instrument after recalibration; Kevin's sample sequence lives here |
| BEFR 122302 | BEFR_H | 122302 | 2023-01-27/25 | 20 + 10 | 21–71° | 7 | one instrument, two dark lagoons, one cal file: separates site from instrument |
| THFR 122302 | THFR_H | 122302 | 2023-01-27/25 | 20 + 10 | 25–74° | 12 | dark water, where Ed-borne wiggles dominate ρw |
| MAFR 121231 | MAFR_H | 121231 | 2021-10-05/04 | 22 + 11 | 23–71° | 12 | turbid, high-Lw regime; the oldest calibration in the network |
| MAFR 122303 | MAFR_H | 122303 | 2023-01-27/25 | 8 + 4 | 23–38° | 3 | the 2026-05 instrument swap at the same site |
| GAIT 121222 pre-recal *(optional)* | GAIT_H | 121222 | 2021-10-05/04 | 7 + 3 | 28–56° | 5 | second before/after-recalibration pair |
| GAIT 120242 *(optional)* | GAIT_H | 120242 | 2024-11-08/12 | 6 + 3 | 23–51° | 4 | loan instrument, 2025 |
| GAIT 121222 post-recal *(optional)* | GAIT_H | 121222 | 2025-03-25/24 | 7 + 3 | 27–57° | 5 | |
| **total** | | | | **150 + 74** | 21–74° | | 130 primary without GAIT |

If 150 is a ceiling, the GAIT rows go first; LPAR is the turbid alternative to
MAFR if MAFR is inconvenient.

**Per sequence:** L1A_IRR, L1A_RAD, L1C_ALL and L2A_REF (as for the VEIT
sample). L1A is the essential part: Release 2 as distributed holds only
L2A/L2B, where Ed is already on the L grid.

**Calibration:** `HYPERNETS_CAL_HYPSTAR_{122302, 122303, 122304, 122305,
121231, 121222, 120242}_{RAD,IRR}_v2.3.nc`, with **both** versions for 122304
(2023-01 and 2024-11) and 121222 (2021-10 and 2025-03); any laboratory
line-spread, slit-function or wavelength-calibration data behind them (the
files' flat 3.0 nm `bandwidth` suggests no per-pixel SRF is stored); and the
processor version that produced Release 2 L1C, since the public
`hypernets_processor` is not up to date with it.

**A note on attributes:** as Kevin flagged, `instrument_calibration_file_rad`
and `instrument_calibration_file_irr` both point to the IRR file
(`HYPERNETS_CAL_HYPSTAR_122304_IRR_v2.3.nc`) in the sample; we see the same in
L1A_RAD, L1C and L2A. We will pair cal files by instrument id and date rather
than by that attribute, and would welcome confirmation of the RAD file names.

## 6. Reference spectra and radiative-transfer tools

- **TSIS-1 HSRS** (Coddington et al. 2021; v2, Coddington et al. 2023) as F₀,
  at its highest native resolution.  The paper's Figure 4 used the same
  spectrum (v2), but at 1 nm, interpolated to 6SV's 2.5 nm spacing -- too
  coarse for eq. (15).
- **HITRAN2020 via HAPI** (Gordon et al. 2022; Kochanov et al. 2016) for O₂
  (A, B, γ bands) and H₂O line-by-line absorption on the 0.01 nm grid, Voigt
  profiles at surface pressure and temperature; column amounts from standard
  atmospheres, with water vapour as a sensitivity variable. HAPI is installed
  in `ocean14` from `~/Oceanography/python/HAPI`.
- **O₃** Chappuis–Huggins cross-sections from tabulated high-resolution
  measurements (Serdyuchenko et al. 2014), 300 DU default. Smooth at 3 nm but
  needed for the broad-band ratio.
- **OSOAA** (Chami et al. 2015; JXP's fork) for the smooth fields only.
  Limits that shape the plan: no gaseous absorption at all (stated in the
  fork's `docs/science/atmosphere_model.rst`), strictly monochromatic, 30 s to
  15 min per wavelength, so it runs on a 5 nm grid (~125 wavelengths per case)
  and never on the 0.01 nm grid. Its Rayleigh, aerosol and hydrosol terms are
  smooth, so this division of labour loses nothing at 3 nm resolution.
- **libRadtran** (REPTRAN fine, ~0.03–0.06 nm) is **deferred**: it would be
  the independent check on the HSRS × HAPI product, built only if the twin
  experiment and the data disagree. Py6S/6SV is too coarse (2.5 nm) except as
  a sanity check; MODTRAN and TAPAS are not open.

## 7. Deliverables and uncertainty

1. **Code in `hypernet/`:** SRF self-calibration, Emod builder, the
   three-method interpolation function, `hypernet/rt/osoaa.py`, tests.
2. **A drop-in snippet** for `hypernets_processor` in the
   `interpolate_wav_linear.py` style, with a per-instrument SRF table (a new
   calibration product: per-channel FWHM(λ) and centroid offset) and a cached
   Emod per SZA bin.
3. **A technical note** (this page, grown): the diagnosis over ~150 sequences,
   the twin-experiment results, the in-situ before/after, and the uncertainty
   budget. Authorship and any later paper to be agreed with Kevin.
4. **Uncertainty:** FIDUCEO-style, following De Vis et al. (2024). The new
   step propagates, with `punpy`, the SRF-parameter covariances, the Emod
   inputs and the grid uncertainties into `u_rel_random` and
   `u_rel_systematic` on Ed_L, so the downstream ρw uncertainties in the
   processor pick it up without change. `punpy` 1.1.0, `comet_maths` 1.0.10
   and `obsarray` 1.0.3 are installed and listed in `requirements.txt`.

## 8. Timeline (weeks from agreement on this plan)

| weeks | activity | depends on |
|---|---|---|
| 0 | send this plan and the CSV to Kevin; Kevin swaps spares; we compute the sky index | — |
| 0–1 | Phase 1 task 0: build OSOAA, `hypernet/rt/osoaa.py`, Level check | gfortran (present) |
| 0–2 | Phase 0 on the VEIT sample: SRF fitting code, FWHM(λ) model, unit tests | VEIT sample (in hand) |
| 1–4 | Phase 1 twin experiment; gate G1 | task 0 |
| 3–6 | Phase 0 on the full request; gate G0 | **data arrival** |
| 5–8 | Phase 2: function, self-calibration, punpy, snippet | G0, G1 |
| 8–11 | Phase 3 in-situ test; gate G3 | Phase 2, data |
| 11–13 | technical note and code hand-over | G3 |

Phases 0 (on the sample) and 1 need nothing from RBINS and start first. The
critical path is the data: every week the request slips moves weeks 3–13 by
the same amount. A separate prompt document will drive the Phase 0/1 work.

## 9. Open issues and risks

- **SRF not stable.** If the E SRF depends on SZA or the diffuse fraction
  (the diffuser filling the slit differently), a fixed per-instrument
  correction fails and the function must self-calibrate per sequence from the
  lines, which needs adequate signal in the blue at low sun. Phase 0 is the
  gate.
- **Ring effect in Ld.** Rotational Raman fills in Fraunhofer lines in sky
  radiance (a few per cent of line depth). It enters ρw through the skyglint
  term ρ_f·Ld, is absent from Ed and from OSOAA, and would look like a
  residual after an otherwise perfect E correction. Mitigation: use widths not
  depths in Phase 0; check whether Phase 3 residuals scale with ρ_f·Ld/Lu;
  model it (Chance & Spurr 1997) only if they do.
- **No per-pixel SRF in the files.** The 3.0 nm `bandwidth` is a placeholder,
  so the SRF table we derive is a new product with its own uncertainty and
  needs the lab data, if any, to anchor it.
- **Gaussian may be the wrong shape.** Line-spread functions of grating
  spectrometers are often asymmetric or have wings; the blue lines in the
  table (Ca H/K, G band) are the sensitive test. Non-Gaussian shapes are a
  Phase 2 extension if G3 fails.
- **Emod line shapes.** With the L SRF in the numerator, errors in the HSRS
  line profiles or in the HAPI Voigt parameters map directly into ρw at line
  centres. Case (vi) of Phase 1 bounds this; libRadtran is the fallback check.
- **Processor drift.** The public `hypernets_processor` is not current with
  Release 2 and unsupported; the snippet is written against the linear
  interpolation class's interface as published, and we rely on RBINS to
  integrate it.
- **Real high-ρw'' features near lines.** Kevin's doc notes that inelastic
  scattering and angular effects at low-transmittance wavelengths can produce
  true structure at the lines. The controls in Phase 1 cover the inelastic
  part; the angular part (direct vs diffuse Ed weighting inside a line) is
  untested and is one reason to keep the direct/diffuse split in the twin
  experiment.
- **Sample size.** Everything in §2.3 rests on one sequence from one
  instrument.

## 10. References

- Chami, M., Lafrance, B., Fougnie, B., Chowdhary, J., Harmel, T. & Waquet, F.
  (2015), "OSOAA: a vector radiative transfer model of coupled
  atmosphere-ocean system for a rough sea surface application to the
  estimates of the directional variations of the water leaving reflectance to
  better process multi-angular satellite sensors data over the ocean",
  *Opt. Express* 23(21), 27829–27852.
- Chance, K., Kurosu, T. P. & Sioris, C. E. (2005), "Undersampling correction
  for array detector-based satellite spectrometers", *Appl. Opt.* 44(7),
  1296–1304.
- Chance, K. & Spurr, R. J. D. (1997), "Ring effect studies: Rayleigh
  scattering, including molecular parameters for rotational Raman scattering,
  and the Fraunhofer spectrum", *Appl. Opt.* 36(21), 5224–5230.
- Coddington, O. M., Richard, E. C., Harber, D., Pilewskie, P., Woods, T. N.,
  Chance, K., Liu, X. & Sun, K. (2021), "The TSIS-1 Hybrid Solar Reference
  Spectrum", *Geophys. Res. Lett.* 48, e2020GL091709.
- Coddington, O. M., *et al.* (2023), "Version 2 of the TSIS-1 Hybrid Solar
  Reference Spectrum and Extension to the Full Spectrum", *Earth Space Sci.*
  10(3), e2022EA002637.
- De Vis, P., Goyens, C., Hunt, S., Vanhellemont, Q., *et al.* (2024), "Generating Hyperspectral Reference Measurements for Surface
  Reflectance from the LANDHYPERNET and WATERHYPERNET Networks",
  *Front. Remote Sens.* 5, 1347230.
- Gordon, I. E., *et al.* (2022), "The HITRAN2020 molecular spectroscopic
  database", *J. Quant. Spectrosc. Radiat. Transfer* 277, 107949.
- Kochanov, R. V., Gordon, I. E., Rothman, L. S., Wcisło, P., Hill, C. &
  Wilzewski, J. S. (2016), "HITRAN Application Programming Interface (HAPI): a
  comprehensive approach to working with spectroscopic data", *J. Quant.
  Spectrosc. Radiat. Transfer* 177, 15–30.
- Ruddick, K. G., De Cauwer, V., Park, Y.-J. & Moore, G. (2006), "Seaborne
  measurements of near infrared water-leaving reflectance: the similarity
  spectrum for turbid waters", *Limnol. Oceanogr.* 51(2), 1167–1179
  (`context/papers/ruddick2006.pdf`).
- Ruddick, K. G., De Vis, P., Goyens, C., Kuusk, J., Lavigne, H. &
  Vanhellemont, Q. (2023), "Second derivative water reflectance spectra for
  phytoplankton species detection — origin, impact and removal of spectral
  wiggles", *Proc. SPIE Remote Sensing*, Amsterdam, September 2023
  (`context/papers/ruddick2023.pdf`; volume and paper number to be confirmed).
- Ruddick, K. G., *et al.* (2024), "WATERHYPERNET: a prototype network of
  automated in situ measurements of hyperspectral water reflectance for
  satellite validation and water quality monitoring", *Front. Remote Sens.* 5.
- Serdyuchenko, A., Gorshelev, V., Weber, M., Chehade, W. & Burrows, J. P.
  (2014), "High spectral resolution ozone absorption cross-sections – Part 2:
  Temperature dependence", *Atmos. Meas. Tech.* 7, 625–636.
- Emde, C., *et al.* (2016), "The libRadtran software package for radiative
  transfer calculations (version 2.0.1)", *Geosci. Model Dev.* 9, 1647–1672.
- CoMet toolkit: `punpy`, `comet_maths`, `obsarray` (NPL),
  <https://www.comet-toolkit.org/> — software; no single paper reference
  confirmed here.

# WATERHYPERNET -- Wiggles Phase 2: Code

## Goal

Turn the Phase 0 SRF model and the Phase 1 interpolator into the deliverable
code in `hypernet/`: one function `interpolate_ed_to_l` with
`method = {"linear", "ruddick2023", "srf"}`, Fraunhofer-line SRF
self-calibration, a cached Emod builder, `punpy` uncertainty propagation into
`u_rel_random` / `u_rel_systematic` in `obsarray` conventions, a per-instrument
SRF table, and a drop-in measurement-function snippet in the style of
`hypernets_processor`'s `interpolate_wav_linear.py`.  We will not run
`hypernets_processor` end to end.

**Gate:** the plan sets no numbered gate for Phase 2.  The exit criterion is
plan §7 items 1, 2 and 4 delivered with tests, in particular: the two special
cases reproduce `np.interp` and eq. (14) to machine precision on synthetic
data, and the G0 and G1 verdicts are applied (table vs self-calibration,
wavelength recalibration or not, default method).

**G1 has closed (2026-10-04, `claude_prompts/wiggles/gate_G1.md`):**

- Default `srf` where the measured ΔFWHM = FWHM_L − FWHM_E ≥ 0.15 nm;
  `linear` otherwise.
- `srf` needs the E wavelength scale to ≤ 0.05 nm relative to Emod (a shift +
  stretch fit; task 4b) and FWHM_E/FWHM_L to ≲ 0.07/0.1 nm.
- A non-Gaussian SRF twin check is added (task 9b).

**G0 has not closed** (Phase 0 task 14 waits on Kevin's L1A delivery); its
preliminary findings are in Context.

## Conventions

- `ocean14`.  Run via `conda run -n ocean14 python ...`; `conda activate` fails
  non-interactively.
- **Model.** Use Opus 5.5 (`claude-opus-5-5`) for this work, including any
  subagents (pass `model: opus`).
- **JXP runs git.**  Claude does not run any state-changing git command.
- Run scripts **from the repository root** so `hypernet` imports resolve.
- Run `pytest -q` after each step where relevant.  92 tests on 2026-10-05; the
  archive-dependent ones skip themselves when `$OS_COLOR` is not mounted, the
  OSOAA ones when the exe is missing.
- **Tier 2** -- the VEIT check needs `$OS_COLOR`.  Do **not** unset it.
- Any calculation goes into a script on disk, not into the chat.
- Data intermediates go **outside** the repo under
  `$OS_COLOR/hypernet/wiggles/ref/` (Emod cache) and
  `$OS_COLOR/hypernet/wiggles/phase2/`.  Only figures and small tables are
  committed.
- Every `.md` in `docs/` is published.  Nothing unpublished goes there.
- Ask questions in the Q&A section below; log completed work under `## Logs`.

## Context

### Where things are

- **Plan:** `docs/wiggles_planning.md` §3 (the equations), §4 "Phase 2", §7
  (deliverables and uncertainty), §9 (Gaussian shape, Emod line shapes,
  processor drift).
- **Gate reports:** `claude_prompts/wiggles/gate_G1.md` (done; read it
  before task 2).  `gate_G0.md` does not exist yet (Phase 0 task 14).  Until
  it does, use the preliminary G0 findings from the Phase 0 Logs (Build #8c,
  #8d, #8f):
  - **G0(a):** E ≠ L on the narrow-E instruments (122302, 122304, 120242:
    Ld − E ≈ 0.5 nm); E ≈ L within 0.07 nm on 121222, 121231, 122303 and
    122305.  H2 stays.
  - **G0(b), preliminary pass:** SRFs are stable to 0.01-0.06 nm over SZA,
    season and recalibration (L channel from Release 2).
  - **G0(c):** air wavelength scale, mean offset +0.007 to +0.047 nm, so no
    relative (L − E) recalibration is needed.  But there is a dispersion trend
    in L (−0.03 to −0.05 nm per 100 nm) that exceeds 0.1 nm at the ends.
    Together with G1's 0.05 nm tolerance, this motivates task 4b.
- **Phase 0 code:** `hypernet/srf.py`:
  - the line fits: `LINES`, `fit_line`, `fit_lines`;
  - the FWHM model: `fit_fwhm_model`, `SRFModel` (JSON);
  - **the HSRS template fit**, the standard SRF measurement since Phase 0
    task 3b: `fit_srf_template` with `via_grid`, `shape` in `SHAPES`
    (gauss/pvoigt/supergauss) and `absorber`, and `fit_template_windows`.
    Gaussian line fits measure SRF ⊗ the line's intrinsic width, so they are
    diagnostics.

  Readers are in `hypernet/whn_l1a.py` (incl. `load_l2b`, `sky_index`),
  `hypernet/whn_srf.py` and `hypernet/refspec.py` (HSRS).
- **Phase 0 SRF products:**
  - `hypernet/data/veit_srf_model.json` (template; the empirical line-width
    version is `veit_srf_model_empirical.json`);
  - `hypernet/data/release2_srf_models.json`: E and Ld per instrument × cal
    period, clear-sky median template fits over the 224 Release 2 sequences
    (E from L2 reads ~0.1 nm wide; E grids are a 122304 proxy except
    122304@2024-11);
  - `$OS_COLOR/hypernet/wiggles/phase0/` (`srf_models.parquet` = VEIT only
    until the delivery; `release2_template.parquet`, `template_fits.parquet`).
- **Phase 1 code:**
  - the minimal `interpolate_ed_to_l` in the flat module
    **`hypernet/edinterp.py`** (Phase 1 Q&A Q1), `METHODS = (linear,
    ruddick2023, srf, cubic, sinc)`, where cubic/sinc are twin nulls;
  - `hypernet/emod.py` (`build_emod`, `ref/emod/` cache);
  - `hypernet/twin.py` (`observe`, `Case`, `case_table`,
    `load_instrument_srfs`, `load_grids`) for synthetic test data;
  - `hypernet/rt/osoaa.py`;
  - the twin scripts `hypernet/wiggles/phase1_twin.py` (`--task 9|10`) and
    `phase1_prediction.py`;
  - `hypernet/tests/test_edinterp.py`.
- **Phase 1 numbers** (in `hypernet/wiggles/`): `phase1_twin_metrics.csv`,
  `phase1_controls.csv`, `phase1_degradation.csv` (the tolerances), and
  `phase1_prediction.csv` (G3's yardstick).
- **hypernets_processor:** public repo `github.com/HYPERNETS/hypernets_processor`,
  file `hypernets_processor/interpolation/measurement_functions/interpolate_wav_linear.py`.
  Kevin's idea doc links it at commit
  `9a12819a3ffb224de9e5f761b1217d130fa69f83` (the version he points to; the
  public repo is not current with Release 2).  Its interface has not yet been
  re-read; Task 2 fetches it.
- **CoMet toolkit:** `punpy` 1.1.0, `comet_maths` 1.0.10, `obsarray` 1.0.3 in
  `ocean14` and in `requirements.txt`.  De Vis et al. (2024) describes the
  processor's use of them.
- **Data:** VEIT sample at `$OS_COLOR/WATERHYPERNET/Wavelengths/`; the full
  request (if delivered by now) alongside it, indexed by Phase 0b at
  `$OS_COLOR/hypernet/wiggles/phase0/request_index.parquet`.
- **Paper:** `context/papers/ruddick2023.pdf`.

### Dependencies

- Depends on G0 (SRF table vs self-calibration; wavelength recalibration) and
  G1 (default method).  If the full data have not arrived, G0 is provisional
  on the VEIT sample: build both the table loader and the self-calibration
  path, and leave the default to be set when G0 closes.
- G1 is closed; G0 is provisional (above).  Tasks 4-6 build both the table
  and the self-calibration paths so that G0 can set the default later.
- Module naming: the phase scripts live in the subpackage `hypernet/wiggles/`
  (run as `python -m hypernet.wiggles.<script>`), and the reusable code in
  flat modules.  `interpolate_ed_to_l` already lives in `hypernet/edinterp.py`
  (Phase 1 Q1), so the default is to harden it there and put the
  self-calibration, wavelength fit and SRF table next to it.  The plan's
  `hypernet/wavecal/` name is superseded unless JXP prefers a subpackage
  (Setup Q&A).

## Prompts

### Setup

1. Read the plan sections, gate reports and Phase 0/1 modules listed in
   Context, and this doc.  Put your questions in Q&A:
   - module layout (default: flat modules next to `hypernet/edinterp.py`);
   - the snippet's location (proposal `hypernet/snippets/`);
   - how to proceed while G0 is provisional;
   - how the method is chosen (where the ΔFWHM ≥ 0.15 nm rule lives, and
     ΔFWHM over which band);
   - the granularity of the E wavelength fit (per sequence or per
     calibration period);
   - whether to target Kevin's pinned commit or the current default branch
     of `hypernets_processor`;
   - the obsarray variable names to target.

   Do not write code yet.  Log your work.

### Build

2. **Read the processor's interface.**  Fetch
   `interpolate_wav_linear.py` from `github.com/HYPERNETS/hypernets_processor`
   at commit `9a12819a3ffb224de9e5f761b1217d130fa69f83` (and diff against the
   default branch if JXP asks) and the module that registers
   measurement functions; save a copy under `context/hypernets_processor/`
   with the commit hash in a `SOURCE.txt`.  Record in the Logs the class name,
   `function` signature, argument names, return convention, and how
   uncertainties and `obsarray` attributes flow through it.  Log your work.

3. **Harden `interpolate_ed_to_l`.**  Final signature per plan §4:
   `interpolate_ed_to_l(wav_irr, irradiance, wav_rad, *, emod, srf_irr,
   srf_rad=None, method="srf", u_rel_random=None, u_rel_systematic=None)`.
   Add input validation, NaN and out-of-range L pixels handled explicitly,
   2-D `irradiance` (wavelength, scan), a precomputed SRF-convolution matrix
   on the 0.01 nm grid cached per `SRFModel`, and full docstrings with the
   equations.
   - Keep `cubic`/`sinc` only as test nulls (or move them to the twin), as
     agreed in Q&A.
   - Add the G1 selection rule, e.g. `choose_method(srf_irr, srf_rad,
     threshold=0.15)` → `"srf"` if the median FWHM_L − FWHM_E over the
     agreed band is ≥ threshold, else `"linear"`.
   - Extend `test_edinterp.py`: the two machine-precision special cases
     remain, plus shape, NaN, edge, 2-D and threshold tests, and a timing test
     that a 1536 → 1538 call runs under a stated budget.

   Log your work.

4. **SRF self-calibration.**  `calibrate_srf(wav, spectrum, channel="E", ...)`
   → `SRFModel` (centroid offset and FWHM(λ) quadratic with covariance).
   - Build it on the HSRS template fit (`fit_template_windows`, then
     `SRFModel.from_lines`), the Phase 0 standard; `fit_lines` stays a
     diagnostic.
   - Add a quality flag for low blue signal (the low-sun risk in plan §9).
   - Test: an `observe()`d synthetic spectrum with a known FWHM(λ) is
     recovered within its reported uncertainty; a noisy one returns the
     flag.

   Log your work.

4b. **E wavelength-scale fit (G1).**  `fit_wavelength_scale(wav, spectrum,
   srf, emod=None)` → shift and linear stretch (and covariance) of a
   spectrum's wavelength scale relative to the HSRS-based model, fitted
   across the template windows.
   - Apply it as a wavelength correction to E before interpolation, at the
     granularity agreed in Q&A.
   - Requirement: ≤ 0.05 nm (G1); report the achieved uncertainty.
   - Tests: a synthetic E with a known shift/stretch (twin case iv) is
     recovered within its uncertainty and to < 0.02 nm at VEIT noise.
   - Run it on the VEIT sample and the Release 2 E spectra.  Compare with
     Phase 0 8d's dispersion trend, in a small table
     `hypernet/wiggles/phase2_wavescale.csv`.

   Log your work.

5. **Emod cache.**  `emod_for(sza, pwv_mm, ozone_du, pressure_hpa)` returning
   the nearest cached entry on a stated grid (SZA bins of 5°, a few water
   vapour columns), populated on demand under `ref/emod/`, plus cached
   E-SRF and L-SRF convolutions per `SRFModel` id.  Test: two calls in the
   same bin return the same array; a cache miss builds and stores.  Log your
   work.

6. **Per-instrument SRF table.**  A committed `hypernet/data/srf_table.csv`
   (instrument, channel, cal period, quadratic coefficients, centroid offset,
   uncertainties, n_sequences, ΔFWHM and the method the G1 rule selects)
   written by `hypernet/wiggles/phase2_srf_table.py`.
   - Source: Phase 0b's `srf_models.parquet` once the delivery has been
     fitted.  Until then use `hypernet/data/release2_srf_models.json`
     (Release 2 template fits, clear sky) plus the VEIT model, flagged as
     provisional.
   - Provide `load_srf_table()` and `srf_for(instrument, channel, date)`.  If G0(b) failed, the default path
   is `calibrate_srf` per sequence and the table is the fallback.  If G0(c)
   fired, apply the centroid offset as a wavelength recalibration before
   interpolation.  Tests on the CSV.  Log your work.

7. **Uncertainty with punpy.**  Propagate with `punpy.MCPropagation`: the
   SRF parameters (centroid and FWHM, with the fit covariance), the Emod
   inputs (F₀ scale, water vapour, ozone, SZA) and the wavelength grids,
   together with the input `u_rel_random` / `u_rel_systematic`, into
   `u_rel_random` and `u_rel_systematic` on Ed_L following the `obsarray`
   conventions found in task 2.  Tests: seeded reproducibility; for
   `method="linear"` with zero SRF/Emod uncertainty the output equals the
   analytic linear propagation of the input uncertainties; systematic and
   random components separate as expected for a common scale error.  Log
   your work.

8. **Drop-in snippet.**  In the location agreed in Q&A: a measurement-
   function class in the style of `interpolate_wav_linear.py` -- same
   `function` and argument-name interface -- taking the SRF table and Emod
   cache as inputs and calling `interpolate_ed_to_l`.  A test imports it and
   runs it on the Phase 1 synthetic scene.  A short module docstring says
   what RBINS must wire up; no end-to-end processor run.  Log your work.

9. **VEIT consistency check.**  `hypernet/wiggles/phase2_veit_check.py`: run the
   three methods on the VEIT L1A_IRR with the sample's SRF models; confirm
   `linear` reproduces the L1C `irradiance` to the 2 × 10⁻⁵ of plan §2.3, and
   show Ld/Ed_L for the three methods with their `u_rel_*` bands.  Figure
   `hypernet/wiggles/figs/phase2/veit_three_methods.png`; table
   `hypernet/wiggles/phase2_veit_check.csv`.  Log your work.

9b. **Non-Gaussian SRF check (G1 concern 3).**
    `hypernet/wiggles/phase2_shape_check.py`: rerun twin case (iii) (and the
    Phase 1 prediction at a few ΔFWHM) with the truth observed through the
    Phase 0 pvoigt/supergauss SRF shapes (`srf.SHAPES`, fitted parameters
    from the VEIT template fits) and the correction still Gaussian.
    - Report how much of the reduction survives, and whether the residual
      correlates with Ed''/Ed (G3's trigger).
    - Table `hypernet/wiggles/phase2_shape_check.csv` and a figure under
      `figs/phase2/`.

    Log your work.

### Gate

10. **Phase 2 exit check.**  Against plan §7 items 1, 2 and 4 and the G0/G1
    verdicts: list each deliverable, its module, its tests and the
    `pytest -q` count; confirm the machine-precision special cases; state the
    default `method` (and the ΔFWHM rule), the SRF source and wavelength-
    scale treatment the snippet ships with, the task 4b and 9b results, and
    what remains provisional while G0 is open.  Write a short report (proposed
    location `claude_prompts/wiggles/phase2_exit.md`; JXP to confirm -- not
    under `docs/`).  Propose the corresponding edits to
    `docs/wiggles_planning.md` §4 and §7 in Q&A and apply them only if JXP
    approves.  Log your work.

## Q&A

### Setup #1 -- 2026-10-07 (Opus 5.5)

**What I read:** plan §3, §4 (Phase 2, and the G1 outcome now in Phase 1),
§7 and §9; `gate_G1.md`; this doc; `hypernet/edinterp.py`, `srf.py`
(`SRFModel`), `emod.py` and `twin.py`; and the VEIT L1C file's irradiance and
uncertainty variables.  Nothing was fetched or run except two read-only
checks: the L1C variable list, and the time of one `srf` call.

**Findings that bear on the questions**

- **The L1C has three uncertainty components, not two.**
  `u_rel_random_irradiance` (err_corr random × random),
  `u_rel_systematic_indep_irradiance` and
  `u_rel_systematic_corr_rad_irr_irradiance` (both `err_corr_matrix` along
  wavelength, systematic along scan; the matrices are named in
  `err_corr_1_params`).  Units are **%**, dims (wavelength, scan); the
  radiances follow the same pattern.  The plan's `u_rel_random` /
  `u_rel_systematic` signature does not match this.  Phase 0 also found
  every L1C `u_rel_random` value to be 0, a placeholder.
- **`edinterp.py` imports `convolve_to_grid` from `twin.py`**, so the
  deliverable depends on the twin-experiment module.
- **Speed is not an issue.**  One 1536 → 1538 `srf` call with an SRFModel
  takes 0.02 s today.  A cached convolution matrix (task 3) is for punpy's
  Monte Carlo (hundreds of draws), not for a single call.
- **`SRFModel` already carries what punpy needs**: `coeffs` with `cov`
  (FWHM quadratic), and `offset` with `offset_err`.
- **G1 tested only E wavelength errors.**  `srf` also evaluates Emod at the
  L pixels, so an L wavelength error misplaces the model's lines relative to
  the measured L.  The twin has not measured this.

**Questions for JXP** (each has a default; say "defaults" to accept them all)

1. **Module layout.**  *Default: flat modules.*
   - *Harden `interpolate_ed_to_l` in `hypernet/edinterp.py`; add
     `choose_method` there.*
   - *Move `convolve_to_grid` (and the cached convolution matrix) into
     `hypernet/srf.py`, so `edinterp` no longer imports `twin`.  `twin`
     re-exports it, so the Phase 1 scripts are unchanged.*
   - *New `hypernet/srfcal.py`: `calibrate_srf` (task 4) and
     `fit_wavelength_scale` (task 4b).*
   - *New `hypernet/srftable.py`: `load_srf_table`, `srf_for` (task 6).*
   - *`emod_for` in `hypernet/emod.py` (task 5).*
   - *The punpy propagation in `edinterp.py`, behind the `u_rel_*`
     arguments (task 7).*

   *The plan's `hypernet/wavecal/` name is dropped.*

   >A. Use the default
2. **The nulls.**  *Default: public `METHODS = ("linear", "ruddick2023",
   "srf")`; `cubic` and `sinc` move to a private `NULL_METHODS` in the same
   module, still accepted by `interpolate_ed_to_l` so that
   `phase1_twin.py` reproduces its tables.*
   >A. Use the default
3. **The snippet.**  *Default: `hypernet/snippets/interpolate_wav_srf.py`,
   a class named after the processor's linear one (e.g.
   `InterpolateWavSRF`; the exact name follows task 2).*
   - *`hypernets_processor` is not installed, so the test uses a minimal
     stand-in for its base class.*
   - *`hypernet/snippets/` is not imported by `hypernet/__init__`.*
   >A. Use the default
4. **The processor version.**  *Default: target Kevin's pinned commit
   `9a12819a…`.  Diff it against the current default branch in task 2 for
   information only, and record any interface change in the Logs.*
   >A. Use the default
5. **Uncertainty components.**  *Default: match the L1C, not the plan's two
   names.*
   - *`interpolate_ed_to_l` takes and returns the three components,
     `u_rel_random`, `u_rel_systematic_indep` and
     `u_rel_systematic_corr_rad_irr`, in % on (wavelength[, scan]), with
     their err_corr matrices.*
   - *The input components are interpolated (linear weights) and
     propagated.*
   - *The new terms (SRF parameters, Emod inputs, the wavelength fit) are
     common to every scan of an instrument, so they are added to
     `u_rel_systematic_indep`, with a punpy-derived wavelength error-
     correlation matrix.*
   - *Plan §4/§7's signature gets updated at the Phase 2 exit (task 10).*
   >A. Use the default
6. **While G0 is open.**  *Default:*
   - *the snippet ships with the SRF table as its source (G0(b) passes
     preliminarily), and `calibrate_srf` per sequence is an option;*
   - *the table's Release 2 rows are flagged `provisional`;*
   - *no L − E centroid recalibration (G0(c) relative passes); E's absolute
     scale is handled by task 4b;*
   - *revisit all three when `gate_G0.md` exists.*
   >A. Use the default
7. **The method rule.**  *Default:*
   - *ΔFWHM = the median of FWHM_Ld(λ) − FWHM_E(λ) over 400-680 nm, where
     both template models are constrained (the VEIT model is unconstrained
     at 680-850 nm).*
   - *It is computed per instrument × calibration period and stored in the
     task 6 table, with the selected method.*
   - *`choose_method(srf_irr, srf_rad, threshold=0.15)` does the same for
     any pair of models; the snippet reads the table and allows an
     override.*
   >A. Use the default
8. **The wavelength-scale fit (task 4b).**  *Default:*
   - *Fit per sequence, on the scan-mean E, against the HSRS-based Emod
     through the E SRF, for the sequence's SZA.*
   - *Fall back to the calibration-period median when the fit is flagged
     (low blue signal) or its uncertainty exceeds 0.03 nm.*
   - *Fit L (Ld) the same way.*
   - *Add an **L** wavelength-error case to the task 9b twin check, since
     G1 never measured it.*
   >A. Use the default
9. **The Emod cache grid (task 5).**  *Default:*
   - *SZA bins of 5° (each entry builds in ~3.5 s);*
   - *water vapour snapped to {7.5, 15, 30} mm;*
   - *ozone 300 DU and pressure 1013 hPa fixed.*

   *G1 found the correction insensitive to SZA ± 10° and water vapour
   × 0.5-3, so this grid is generous.*
   >A. Use the default
10. **Edges.**  *Default: L pixels whose SRF reaches beyond the Emod range
    (380-1000 nm), or beyond the E grid, get `linear` and a flag; NaN
    inputs propagate as NaN at the affected pixels only.*
   >A. Use the default
11. **The exit report.**  *Default:
    `claude_prompts/wiggles/phase2_exit.md`, as proposed in task 10.*
   >A. Use the default

### Build #2 -- 2026-10-07 (Opus 5.5)

The interface is recorded in the Logs ("Build #2").  It raises two choices
for later tasks.

1. **How the new uncertainty terms enter (task 7).**  In the processor,
   punpy propagates only input *dataset variables* (their `u_rel_*` and
   err_corr), and the water default is `mcsteps = 0`.  *Default: both
   routes, with the first one shipped.*
   - *(a) **Shipped in the snippet:** the SRF parameters (FWHM_E, FWHM_L
     coefficients) and the E/L wavelength shift and stretch become extra
     measurement-function inputs.  Their uncertainties sit as
     `u_rel_systematic_indep_*` variables in the irradiance dataset, so the
     processor propagates them without code changes once `mcsteps > 0`.*
   - *(b) **In `interpolate_ed_to_l`:** its own punpy MC that returns the
     three components directly, for our Phase 3 use and for tests.*
   - *(a) and (b) must agree in a test.*
   >A. Go with the default
2. **Telling Kevin about `InterpolationWavClearSky`.**  The processor
   already has an eq.-14-style interpolator.  It uses a 3 nm placeholder SRF
   and, as written, returns values on the irradiance pixels rather than the
   radiance pixels.  *Default: no message now.  Mention it in the Phase 2
   exit report and in the next note to Kevin, as context for why the
   measured-SRF version matters, and ask whether RBINS ever ran it on water
   data.*
   >A. Go with the default

### Build #4b -- 2026-10-10 (Opus 5.5)

The E and Ld window offsets carry a **common per-window pattern**:
- 0.13 nm rms about each sequence's straight line;
- the same on all 7 instruments and both channels (correlation 0.54-0.95,
  median 0.79);
- large at 415 (−0.20), 445 (+0.14), 505 (−0.23), 535 (+0.12), 605 (+0.25)
  and 855-865 nm (−0.10 to −0.15).

Spectrometers with different pixel grids cannot share a pixel dispersion
error, so it comes from the reference or the method.  I now subtract it
before the shift + stretch fit (shipped as
`hypernet/data/wavescale_window_pattern.csv`).  That brings χ²ν from ~5 to
~1.5 and the per-sequence error to 0.02-0.05 nm.

1. **Where the pattern comes from.**  There are two candidates, and they
   matter differently for `srf`:
   - **(a) A method artifact.**  An asymmetric HYPSTAR SRF (all HYPSTARs
     share the design), or line-shape mismatch with a Gaussian SRF, shifts
     each window's apparent centroid by an amount set by its line mix.
     Then nothing is wrong with any wavelength.  The pattern only has to be
     removed from the scale fit, as now.
   - **(b) HSRS line positions or shapes off by ~0.1-0.3 nm in some
     windows.**  Then Emod carries them too, and `srf` places the model's
     lines wrongly there.  That costs gain in those windows (the G1
     tolerance is 0.05 nm).

   *Default: extend task 9b (the non-Gaussian SRF twin check) with an
   **asymmetric** SRF.  Observe the twin through it, run the template /
   wavelength-scale fit, and see whether the window pattern appears.  If it
   does, it is (a).  If not, compare HSRS with a second solar atlas (e.g.
   Kurucz 2005 or the Chance & Kurucz 2010 irradiance) in the worst
   windows.*
2. **The Q8 fallback threshold.**  Per-sequence errors with the pattern
   removed are 0.021/0.027 nm on VEIT L1A (E/Ld, with per-pixel errors).  On
   Release 2 L2B (no per-pixel errors) the medians are 0.023-0.053 nm, and
   0-100 % of sequences per instrument pass 0.03 nm.  *Default: keep
   0.03 nm.  The fallback is the calibration-period median (task 6 stores
   it); the period scatter is 0.01-0.03 nm, so the median is good to
   ≲ 0.01 nm.  Revisit with the L1A delivery, which has per-pixel errors.*
3. **Is the correction needed?**  *Default: yes, apply it.*  E's scale
   (pattern removed) has a shift of 0 to +0.05 nm and a stretch of +0.03 to
   +0.08 nm.  At the ends that is up to ~0.13 nm, above G1's 0.05 nm.  Ld's
   stretch is +0.04 to +0.11 nm.

>A. Go with all of your defaults

### Build #6 -- 2026-10-10 (Opus 5.5)

1. **122303 sits on the threshold.**  Its Release 2 ΔFWHM (median 400-680 nm)
   is **0.138 nm**, which the G1 rule maps to `linear`.  But E from L2
   reads too wide:
   - the VEIT validation (Phase 0c) puts the L2 ΔFWHM bias at +0.064 nm;
   - VEIT's own L1A ΔFWHM is 0.477 vs 0.390 from L2 (+0.087).

   With the bias corrected, 122303 is **0.20 nm → `srf`**, and 121231 is
   0.146 (just under).  Phase 0 put 122303 among the E ≈ L instruments
   ("within 0.07 nm") on a different statistic.
   *Default: ship the L2-measured method (`linear` for 122303), with
   `method_sensitive = True` in the table.  Decide 122303 and 121231 from
   their L1A E widths when the delivery arrives (task 6 reruns from Phase
   0b's `srf_models.parquet`).*

>A.  Your default is fine

## Logs

### 2026-10-05 -- Doc refreshed after G1 (Opus 5.5)

- Brought the doc up to date before Setup #1, from the Phase 1 session:
  - **Goal:** the G1 outcome, and G0 still open.
  - **Context:**
    - the preliminary G0(a)/(b)/(c) findings;
    - the HSRS template fit as the standard SRF measurement;
    - `release2_srf_models.json`;
    - `hypernet/edinterp.py` as the interpolator's home;
    - the Phase 1 scripts and tables.
  - The test count, now 92.
- **Tasks:**
  - Setup #1: new questions on the method rule and the wavelength-fit
    granularity.
  - Task 3: the `choose_method` rule (ΔFWHM ≥ 0.15 nm).
  - Task 4: built on the template fit.
  - **New task 4b:** the E wavelength shift + stretch fit (G1, ≤ 0.05 nm).
  - Task 6: the Release 2 source and the ΔFWHM/method columns.
  - **New task 9b:** the non-Gaussian SRF twin check.
  - Task 10: states the new items.
- Nothing else changed; no code.

### 2026-10-07 -- Setup #1 (Opus 5.5)

- Read the plan (§3, §4 Phase 2 and the G1 outcome, §7, §9), `gate_G1.md`,
  this doc and the Phase 0/1 modules.  Wrote 11 questions in Q&A "Setup #1",
  each with a default.  No code written.
- **Learned:**
  - The L1C carries three relative-uncertainty components in %
    (`u_rel_random_`, `u_rel_systematic_indep_`,
    `u_rel_systematic_corr_rad_irr_` + variable), with err_corr matrices
    named in `err_corr_1_params`.  The plan's two-name signature needs
    updating (Q5).
  - `edinterp` depends on `twin` through `convolve_to_grid` (Q1).
  - One `srf` call takes 0.02 s.
  - `SRFModel` already holds `cov` and `offset_err` for punpy.
  - G1 never tested an **L** wavelength error, which `srf` is exposed to
    through Emod_L (Q8).
- Read-only checks only: the VEIT L1C variable list and one timed call.


### 2026-10-07 -- Build #2: the processor's interface (Opus 5.5)

- JXP accepted all 11 Setup #1 defaults.
- **Fetched** `hypernets_processor` at Kevin's pinned commit `9a12819a…`
  (2026-03-20) into `context/hypernets_processor/` (same paths as the repo;
  `SOURCE.txt` lists them):
  - the `interpolation/` package (driver, factory, the wavelength and time
    measurement functions);
  - the water and land default configs;
  - the ATBD page `interpolate.rst`.
- **Diff against `main` (Q4):** none.  On 2026-10-07 `main` *is*
  `9a12819a…`; the public repo has not moved since March.
- **The class to mirror:** `InterpolationWavLinear(punpy.MeasurementFunction)`
  in `interpolation/measurement_functions/interpolate_wav_linear.py`.
  - `meas_function(self, rad_wavs, irr_wavs, irr)` returns the irradiance on
    `rad_wavs`, via `scipy.interpolate.interp1d(irr_wavs, irr, axis=0,
    fill_value="extrapolate")`.  It **extrapolates** at the ends (unlike
    `np.interp`, which holds the end values).  If `irr_wavs` is 2-D (punpy's
    broadcast), it uses column 0.
  - `get_name()` (staticmethod) returns `"InterpolationWavLinear"`, the key
    the config selects.
  - `get_argument_names()` returns `["radiance_wavelength", "wavelength",
    "irradiance"]`: the dataset variable names punpy reads, in argument
    order.
  - No `setup` (it inherits punpy's no-op).
- **How it is called** (`interpolation/interpolate.py`,
  `Interpolate.interpolate_irradiance`):
  - **Registration:** `interpolation_factory.py` holds a list
    `interpolators`.  `InterpolationFactory(*args, **kwargs)` instantiates
    each and keys it by `get_name()`.  The config key
    `measurement_function_interpolate_wav` picks one; water and land both
    default to `InterpolationWavLinear`.  **A new class needs one import and
    one list entry there, plus the config value.**
  - **Construction (water):** `factory(prop=punpy.MCPropagation(mcsteps,
    dtype="float32", parallel_cores=1), corr_dims="wavelength",
    yvariable="irradiance")`.
  - **`setup` is always called:** `setup(np.nanmean(SZA of the irradiance
    scans), network, irradiance wavelengths)`.  There is **no instrument ID
    or date**, which our snippet needs to pick its SRFs.
  - **Inputs:** the L1C dataset with `wavelength` renamed to
    `radiance_wavelength`, and the L1B irradiance dataset (`wavelength`,
    `irradiance` on (wavelength, scan)).
  - **If `mcsteps > 0`:** `propagate_ds_specific(["random",
    "systematic_indep", "systematic_corr_rad_irr"], ds_l1c, ds_l1b_irr,
    ds_out_pre=…, use_ds_out_pre_unmodified=True, store_unc_percent=True)`.
    punpy draws the inputs from their `u_rel_<comp>_<var>` and `err_corr_*`
    (obsarray) attributes, runs `meas_function` per draw, and writes
    `u_rel_<comp>_irradiance` in % with the error-correlation matrices.
    **The measurement function computes only the value; uncertainty flows
    through the datasets.**
  - **Otherwise:** `run(...)` gives only the value.
  - Then a separate time interpolation (`InterpolationTimeLinearCoscorrected`)
    takes Ed to the Lu scan times, with the same three components.
- **The water default config has `mcsteps= 0`** (land: 100).  So the
  operational water L1C propagates **no** uncertainty through interpolation.
  This most likely explains Phase 0's all-zero L1C `u_rel_random`.
- **A model-based interpolator already exists:** `InterpolationWavClearSky`.
  - Its `setup` loads a clear-sky high-resolution BOA irradiance for the
    nearest of SZA {0, 10, 20, 40, 60, 70, 80}
    (`data/quality_comparison_data/solar_irradiance_hypernets_sza*_highres_w.nc`,
    in the repo).
  - It uses a **3.0 nm Gaussian "bandwidth" for every irradiance pixel**
    (10 nm in the SWIR) and comet_maths `interpolate_1d_along_example`, then
    band-integrates with `matheo`.
  - It is in effect eq. 14 with a placeholder SRF, not the measured E or L
    SRFs.
  - As written, **its output is band-integrated on the *irradiance* pixels:
    `rad_wavs` is never used.**  Unless that is intended upstream, it would
    not return values on the L grid.
  - The ATBD says it "has been implemented, but is not currently
    operationally used … [pending whether the improvement justifies] their
    significantly slower runtime".
  - Worth telling Kevin: our snippet is the same idea with measured, per-
    channel SRFs.
- **Versions:** the processor requires `punpy>=1.0.0`, `comet_maths>=1.0.0`,
  `obsarray>=1.0.0` and `matheo`.  `ocean14` has punpy 1.1.0 and obsarray
  1.0.3, compatible.  `matheo` is not installed (only the clear-sky class
  needs it).
- **Consequences for later tasks:**
  - **Task 8 (snippet).**
    - Subclass `punpy.MeasurementFunction` with the same
      `meas_function(rad_wavs, irr_wavs, irr)` and argument names.
    - Do the SRF-table and Emod loading in `setup(szas, network, irr_wavs,
      system_id=None, date=None)`; the extra keywords are what RBINS must
      wire in (one call site in `interpolate.py`).
    - Match the linear class's extrapolation at the ends for `linear`.
    - Handle 2-D `irr_wavs` the same way.
  - **Task 7 (punpy).**
    - In the processor, punpy propagates only dataset variables.  So for the
      new SRF, Emod and wavelength-fit terms to flow through unchanged
      processor code, they must enter as extra **input variables** of the
      measurement function: e.g. `get_argument_names()` += SRF parameters
      carried in the irradiance dataset with their own
      `u_rel_systematic_indep_*` and err_corr.
    - The alternative is a post-hoc quadrature term, which the processor
      would not apply without changes.
    - I'll put this design choice in Q&A before task 7.
  - **Precision:** the MC runs in float32.
- `context/` is outside `docs/`, so nothing new is published.  No package
  code changed; `pytest -q` not rerun (92 at task 11).
- Slip: a first fetch loop relied on word-splitting that zsh doesn't do, and
  created empty stray directories (`hypernets_processor/…` at the repo root
  and under `context/`).  I removed them (empty, untracked) and refetched
  with a bash array.

### 2026-10-08 -- Build #3: `interpolate_ed_to_l` hardened (Opus 5.5)

- JXP accepted both Build #2 defaults (uncertainty routes (a) + (b); mention
  `InterpolationWavClearSky` in the exit report and the next note to Kevin).
- **`hypernet/srf.py` (Setup Q1).**  New pixel-convolution section:
  - `convolve_to_grid` moved here from `twin.py`; `twin` re-exports it, so
    the Phase 1 code is unchanged;
  - `convolution_matrix(lam_hr, grid, fwhm, nsig=5, cache=True)` returns a
    sparse CSR matrix, rows normalised;
  - an LRU cache (64 entries) keyed on the Emod grid, the pixel grid and
    `srf_key(fwhm)`.  SRFModels, constants and arrays are cached; plain
    callables are not (two lambdas share nothing hashable);
  - `fwhm_values(fwhm, grid)` evaluates any SRF spec.

  `edinterp` no longer imports `twin`.
- **`hypernet/edinterp.py` rewritten.**  Module docstring with the three
  equations, plus:
  - `METHODS = (linear, ruddick2023, srf)`; `NULL_METHODS = (cubic, sinc)`,
    still accepted (Q2).
  - `UNC_COMPONENTS` = the three L1C names (Q5); `SRF_THRESHOLD = 0.15`;
    `DFWHM_BAND = (400, 680)` (Q7).
  - `delta_fwhm(srf_irr, srf_rad, band)` and `choose_method(...)`: G1's
    rule, inclusive with a 1e-9 nm tolerance.  On the Release 2 models,
    122304 → `srf` and 122305 → `linear`.
  - `linear_weights(wav_irr, wav_rad)` → (i, w), with `np.interp`'s hold at
    the ends.
  - `interpolate_ed_to_l(wav_irr, irradiance, wav_rad, *, emod, srf_irr,
    srf_rad, method, u_rel_random, u_rel_systematic_indep,
    u_rel_systematic_corr_rad_irr, full_output)`:
    - **validation:** wav_irr 1-D, finite, strictly increasing; wav_rad
      finite; irradiance (n_E,) or (n_E, n_scan); method; required SRFs and
      Emod;
    - **2-D irradiance:** the model factors are computed once, and the ratio
      is interpolated per scan;
    - **NaN locality:** a NaN in E reaches only the L pixels that use it with
      non-zero weight;
    - **edges (Q10):** an L pixel beyond the E grid, or whose SRF (or its E
      neighbours', ±3σ) reaches past Emod, gets `linear` and
      `info['fallback'] = True`;
    - **uncertainty:** the input components are propagated to first order.
      Random adds in quadrature over the two neighbours, weighted by their
      share of the (model-ratio) sum.  Systematic is taken as fully
      correlated between neighbours; task 7 refines this with the err_corr
      matrices and adds the SRF/Emod/wavelength terms.
    - **Returns** `Ed`, or `(Ed, info)` when `full_output` (the default if any
      `u_rel_*` is given): `method`, `fallback`, `u_rel_<component>` in %.
- **Speed.**  On a full-size 1536 → 1538 `srf` call with SRFModels and a
  0.01 nm 380-1000 nm Emod: 25 ms cold, **4.5 ms** with the matrices cached.
- **Tests:** `test_edinterp.py` grew from 6 to 16 tests.
  - The machine-precision special cases are unchanged.
  - New: 2-D vs per-column, validation errors, NaN locality, edge fallback,
    the threshold rule (synthetic and on the 122304/122305 Release 2
    models), the analytic linear propagation (random and systematic), a
    constant systematic preserved by all three methods, cache identity, and
    timing (budget < 1 s cold, < 50 ms warm).
- **Test fix in `test_twin.py`:** `test_delta_line_observed_fwhm3` failed
  after the move.  The new convolution differs from the old loop by
  3e-15, but on exactly noise-free data `curve_fit`'s covariance is
  undefined, so the test had passed on rounding luck.  It now adds 1e-5
  seeded noise, as the other synthetic fit tests do.
- **Phase 1 reproduces:** `hypernet/wiggles/phase2_check_phase1_repro.py`
  reruns task 9's metrics for one scene.  The maximum relative difference
  from the stored parquet is 1.9e-10 (VEIT-like scene, 340 rows), on a value
  at the numerical floor.  It writes nothing.
- `pytest -q`: **102 passed** (was 92).

### 2026-10-08 -- Build #4: SRF self-calibration (Opus 5.5)

- New flat module `hypernet/srfcal.py` (Setup Q1):
  - `calibrate_srf(wav, spectrum, channel='E', err=None, ref=None,
    windows=None, veil=False, via_grid=None, shape='gauss', scale_cov=True,
    instrument='', return_fits=False, meta=None)` → `SRFModel`, and the
    window table with `return_fits`.
  - It runs `srf.fit_template_windows` against TSIS-1 HSRS in air (loaded
    once, `default_reference`), then `SRFModel.from_lines`: FWHM(λ)
    quadratic with covariance (chi²-inflated), and the weighted-mean centroid
    offset with its error.
  - **Veil fixed at 0** by default, as in the shipped VEIT models.
  - Gaussian `fit_lines` is not used (a diagnostic only).
- **Quality flags** in `model.meta['flags']`; `meta['ok']` only if there
  are none.
  - **`low_blue_signal`** (plan §9's low-sun risk) is raised by any of:
    - fewer than 3 good windows in 390-450 nm;
    - median blue S/N < 100 (when `err` is given);
    - FWHM error at 410 nm > 0.05 nm (inside G1's ≲ 0.07 nm FWHM_E
      tolerance).
  - Also `few_windows` (< 8 good windows) and `fit_failed`.
  - `meta` also records `n_windows`, `n_blue`, `blue_snr`,
    `fwhm_err_blue`, `method`, `veil` and `shape`.
- **Tests** (`hypernet/tests/test_srfcal.py`, 5):
  - **Synthetic:** a 1500-line "solar" spectrum observed by `twin.observe`
    with a known quadratic FWHM(λ) (2.72 nm at 410 → 2.27 at 700) and 0.1 %
    noise.
    - FWHM at 410-850 nm is recovered with |pull| < 3 (≈ ±1.2σ in
      practice) and errors < 0.03 nm; offset consistent with 0.
    - A +0.05 nm sampling shift is recovered as offset −0.05 nm (error
      < 0.01).
    - With the blue scaled to 2 % (S/N ≈ 17), the fit returns
      `low_blue_signal`; it is also raised without `err`, from the FWHM
      error.
    - Mismatched lengths raise; an all-NaN spectrum gives `fit_failed` +
      `few_windows`.
  - **VEIT (tier 2, ran):** the VEIT E (scan mean, flattened scan-scatter
    errors) reproduces `hypernet/data/veit_srf_model.json` E within
    0.01 nm at 420-660 nm.  It passes the flags with margin: blue S/N 455,
    FWHM error at 410 nm 0.029 nm, 27 windows (6 blue); offset
    +0.058 ± 0.018 nm.
  - A synthetic run takes ~0.1 s.
- `pytest -q`: **107 passed**.
- Not done here: the per-sequence fallback to the calibration-period median
  (Q8) belongs to task 4b/6; whether the 0.05 nm / S/N 100 thresholds catch
  real low-sun sequences is best checked on the delivery (or Release 2 E
  at SZA ≥ 60°) in task 6.

### 2026-10-10 -- Build #4b: wavelength-scale fit (Opus 5.5)

- **`hypernet/srfcal.py`** gained:
  - **`fit_wavelength_scale(wav, spectrum, srf=None, ref=None, err=None,
    windows=None, via_grid=None, pattern=None, return_fits=False)`.**
    - In each template window it fits the centroid offset dlam with the SRF
      width fixed by `srf` (`_fit_window_shift`: continuum × the reference
      through a Gaussian, evaluated on `via_grid` and resampled for L2
      irradiance), or free (`fit_srf_template`) if `srf` is None.
    - It subtracts the optional common `pattern`, then fits a weighted
      straight line to −dlam (chi²-inflated covariance).
    - Convention (twin case iv): a pixel reported at λ samples
      λ + shift + stretch (λ − 630)/240.
    - Returns shift, stretch, cov, their errors, `max_err` (largest 1σ of
      the correction over 390-870 nm), `slope_per100nm` (8d's dlam slope),
      chi2_nu, n_windows, and flags `uncertain` (max_err > 0.03, Q8),
      `few_windows` and `fit_failed`.
  - `correct_wavelengths(wav, scale)`: the true grid.  Apply it to E's grid
    before `interpolate_ed_to_l`.
  - `scale_or_fallback(scale, fallback)`: Q8's per-sequence-else-period
    rule.  The period medians come with task 6's table.
  - `load_window_pattern()`, and the constants `LAM0`, `HALF`,
    `MAX_SCALE_ERR = 0.03` and `G1_TOLERANCE = 0.05`.
- **Tests** (`test_srfcal.py`, 5 → 12):
  - **Known shift/stretch** (twin case iv: (0.1, 0), (−0.05, 0.08),
    (0.03, −0.05)) at VEIT-like noise (0.2 %/pixel): recovered within 3σ
    and to < 0.02 nm, with max_err < 0.03.  Correcting the grid removes the
    error.  SRF free and SRF fixed agree.
  - The `via_grid` (L2 resampling) case is recovered.
  - 5 % noise gives `uncertain`, and the fallback is used.
  - A constant pattern moves the shift by exactly that constant.
  - The shipped pattern file loads.
- **`hypernet/wiggles/phase2_wavescale.py`.**  VEIT L1A (E and Ld, with
  per-pixel errors, `veit_srf_model.json` widths) and the 224 Release 2 L2B
  (Ld native; E through the proxy 122304 E grid; `release2_srf_models.json`
  widths; unweighted), each with and without the pattern.  450 fits in 52 s.
  Writes:
  - `$OS_COLOR/.../phase2/wavescale_per_sequence.parquet`;
  - `hypernet/wiggles/phase2_wavescale.csv` (per instrument × channel);
  - `hypernet/wiggles/phase2_wavescale_windows.csv` (the per-window
    pattern per group);
  - `hypernet/data/wavescale_window_pattern.csv` (the common pattern; now
    in `package_data`, with `data/*.csv` added to `setup.py`).

  Diagnostic: `hypernet/wiggles/phase2_wavescale_veit_windows.py` prints
  VEIT's window offsets.
- **Results:**
  - **The slope agrees with Phase 0 8d.**  Raw dlam slope per 100 nm: Ld
    −0.019 to −0.055, E −0.006 to −0.026; 8d had −0.026 to −0.050 and
    −0.005 to −0.024.
  - **Without the pattern, no sequence reaches 0.03 nm.**  Per-sequence
    max_err is 0.05-0.09 nm and χ²ν ≈ 5 (VEIT 29).
  - **The common pattern:** 0.13 nm rms, correlated 0.54-0.95 across all
    7 instruments × 2 channels, with a standard error of 0.013 nm per
    window.  Not instrumental (see Q&A "Build #4b").
  - **With it removed:**
    - χ²ν ≈ 1-2.
    - Per-sequence max_err medians are 0.021-0.048 nm for E and
      0.027-0.053 for Ld.  All are inside G1's 0.05 nm but one (121231 Ld,
      0.053).
    - VEIT: E 0.021, Ld 0.027 nm; both pass Q8's 0.03.
  - **Scale errors** (pattern removed, medians):
    - E: shift −0.03 to +0.05 nm, stretch +0.01 to +0.08 nm;
    - Ld: shift −0.04 to +0.05, stretch +0.04 to +0.11.
    - At the spectrum ends the E correction reaches ~0.1 nm, which exceeds
      G1's tolerance, so the correction is needed.
    - Sequence-to-sequence scatter within an instrument is 0.01-0.03 nm, so
      the calibration-period median is a good fallback.
- `pytest -q`: **114 passed**.
- Open (Q&A "Build #4b"): the pattern's origin, a method artifact (e.g. an
  asymmetric SRF) or HSRS.  It decides whether `srf` inherits it.  The
  default is an asymmetric-SRF case in task 9b.  The pattern was estimated
  from the same Release 2 data it is applied to.

### 2026-10-10 -- Build #5: the Emod cache (Opus 5.5)

- JXP accepted all three Build #4b defaults.  The asymmetric-SRF case joins
  task 9b; the fallback threshold stays 0.03 nm; the scale correction is
  applied.
- **`hypernet/emod.py`, new section** (Setup Q9):
  - **The grid:** `SZA_STEP = 5` (clipped to 0-85°), `PWV_GRID = (7.5, 15,
    30)` mm, `OZONE_FIXED = 300` DU, `PRESSURE_FIXED = 1013.25` hPa, and
    `DIFFUSE_FRACTION = 0.15`.  `snap(sza, pwv)` gives the grid point:
    nearest 5°, and nearest PWV in log.
  - **`emod_for(sza, pwv_mm=15, ozone_du=None, pressure_hpa=None,
    diffuse_fraction=0.15, cache_dir=None, builder=None)`.**
    - It snaps, then returns one in-memory entry per bin: **the same object
      for every call in the bin**.
    - A miss calls `build_emod`, which reads or writes its npz under
      `ref/emod/` (~3.5 s to build).
    - The entry adds `E = F0[(1 − f_d) T_direct + f_d T_diffuse]`, a model
      *shape* (no cos SZA, no aerosol, which the measured/model ratio
      absorbs), plus `sza_bin`, `pwv_bin`, the fixed ozone/pressure actually
      used, and a `key`.
    - The ozone and pressure arguments are accepted but fixed by the grid,
      and the values used are reported.
  - **`emod_on_grid(entry, grid, srf)`**: Emod_X on a grid, cached per
    (entry, grid, SRF key), returned read-only, LRU of 64.  The sparse
    matrix underneath is already cached by `srf.convolution_matrix`.
  - `clear_memory_cache()`.
- **`edinterp.model_on_grid`** now uses `emod_on_grid` when given an
  `emod_for` entry, so `interpolate_ed_to_l(..., emod=emod_for(sza, pwv))`
  reuses both convolutions across calls and scans.
- **Tests** (`test_emod.py`, +4):
  - **Tier 1, with a fake builder:**
    - the snapping (including the log-nearest PWV and the clipping);
    - two calls in a bin return the same object and build once, with the
      fixed ozone and pressure passed; another bin is a miss;
    - `emod_on_grid` is cached and read-only, and `edinterp.model_on_grid`
      returns the cached array;
    - `interpolate_ed_to_l` with an entry reproduces Emod_L when
      E = Emod_E.
  - **Tier 2 (ran):** a miss builds into a tmp dir and stores
    `emod_sza35.00_pwv15.0…npz`.  After clearing memory, a call in the same
    bin (SZA 34) reads it back (`cached`), identical.
- `pytest -q`: **118 passed**.

### 2026-10-10 -- Build #6: the per-instrument SRF table (Opus 5.5)

- No new JXP answers since Build #5.
- **`hypernet/wiggles/phase2_srf_table.py`** writes
  **`hypernet/data/srf_table.csv`**: 20 rows, the 7 instruments × E/Ld ×
  9 calibration periods, plus the VEIT L1A pair.  All rows are
  `provisional`.
  - **Sources:**
    - `release2_l2b`: `release2_srf_models.json`, the clear-sky median
      template fits;
    - `veit_l1a`: `veit_srf_model.json`, preferred for 122304 @
      2024-11-08/12.
    - When the delivery is fitted, a `phase0b_l1a` source ranks first
      (`SOURCE_RANK`).
  - **Per row:**
    - the FWHM quadratic and its 6 covariance terms, lam_ref/scale/min/max,
      chi2, and the centroid offset ± error;
    - `fwhm_*` and `fwhm_err_*` at 400/550/700 nm;
    - `valid_from`, **the channel's own calibration date** (IRR for E, RAD
      for Ld), and `n_sequences`, `e_grid` and a note;
    - `delta_fwhm` (`edinterp.delta_fwhm`, median over 400-680 nm) and
      `method` (G1 rule);
    - for L2 rows, `delta_fwhm_l1a_est` = delta + 0.064 nm, the L2-vs-L1A
      width bias from `phase0c_validation.csv` (Ld − E bias, median
      400-680 nm, one VEIT sequence), with `method_l1a_est` and
      `method_sensitive`;
    - the **wavelength-scale correction**: calibration-period medians of
      task 4b's shift and stretch (pattern removed), with errors
      1.2533 sd/√n, plus n and frac_ok.  The VEIT rows use their single
      sequence and its fit errors.
- **`hypernet/srftable.py`:**
  - `load_srf_table`, `row_for(instrument, channel, date)` (latest
    `valid_from` ≤ date, then the preferred source; dates as
    'YYYY-MM-DD', datetime or 'YYYYMMDDTHHMM');
  - `model_from_row`, `srf_for(...)` → SRFModel, `scale_for(...)` → dict
    for `srfcal.correct_wavelengths` / `scale_or_fallback`, and
    `method_for(instrument, date)`;
  - **G0 switches, provisional:** `TABLE_IS_DEFAULT = True` (G0(b)
    preliminary pass) and `APPLY_CENTROID_OFFSET = False` (G0(c) relative
    pass; the task 4b scale correction is applied instead).  The docstring
    says how to flip them when G0 closes.
- **Results:**
  - **`srf`:** 120242 (ΔFWHM 0.39), 122302 (0.46), 122304 (0.49 @2023;
    0.39 L2 / **0.48 L1A** @2024-11).
  - **`linear`:** 121222 (−0.05, −0.03), 121231 (0.08), 122305 (0.01),
    and 122303 (0.138).
  - The L1A-estimated ΔFWHM flips **122303 to `srf` (0.20)**, and leaves
    121231 at 0.146, just under 0.15.  See Q&A "Build #6".
  - **Scale corrections** (period medians): E shift 0.00-0.05 nm, stretch
    0.03-0.08; Ld shift −0.05 to +0.05, stretch 0.03-0.11; errors ≤ 0.016.
- **Tests** (`hypernet/tests/test_srftable.py`, 6):
  - the CSV's columns, instruments and provisional flags;
  - the G1 rule row by row, and one method per pair;
  - valid_from = the channel's calibration date;
  - all covariances PSD, with `model_from_row` reproducing `fwhm_550`;
  - `srf_for` preference and dates, including RAD/IRR dates that differ
    (on 2024-11-10, Ld uses the new period and E the old), and KeyErrors;
  - `method_for` agrees with `choose_method` on the looked-up pair;
  - `scale_for` with `correct_wavelengths` and the fallback;
  - the G0 switches.
- `pytest -q`: **124 passed**.
- Not yet tested: whether `calibrate_srf`'s low-blue flag catches real
  low-sun sequences.  Release 2 E comes from L2 (no per-pixel errors, the
  proxy grid), so this waits for the L1A delivery.

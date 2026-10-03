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

## Conventions

- `ocean14`.  Run via `conda run -n ocean14 python ...`; `conda activate` fails
  non-interactively.
- **Model.** Use Opus 5.5 (`claude-opus-5-5`) for this work, including any
  subagents (pass `model: opus`).
- **JXP runs git.**  Claude does not run any state-changing git command.
- Run scripts **from the repository root** so `hypernet` imports resolve.
- Run `pytest -q` after each step where relevant.  21 tests today (more after
  Phases 0 and 1); the archive-dependent ones skip themselves when `$OS_COLOR`
  is not mounted.
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
- **Gate reports:** `claude_prompts/wiggles/gate_G0.md` and `gate_G1.md`
  (locations to be confirmed there).  Read both before task 2.
- **Phase 0 code:** the SRF module (default `hypernet/srf.py`: `LINES`,
  `fit_line`, `fit_lines`, `fit_fwhm_model`, `SRFModel`) and the readers
  (default `hypernet/whn_l1a.py`); `hypernet/data/veit_srf_model.json`;
  Phase 0b's `$OS_COLOR/hypernet/wiggles/phase0/srf_models.parquet`.
- **Phase 1 code:** the minimal `interpolate_ed_to_l` (module per Phase 1
  Q&A), `hypernet/emod.py` (`build_emod`, `ref/emod/` cache),
  `hypernet/twin.py` (`observe`, cases) for synthetic test data,
  `hypernet/tests/test_edinterp.py`.
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
- Module naming: since 2026-09-29 the phase scripts live in the subpackage
  `hypernet/wiggles/` (run as `python -m hypernet.wiggles.<script>`), and the
  reusable code in flat modules (`hypernet/srf.py`, `hypernet/whn_l1a.py`,
  `hypernet/whn_srf.py`, `hypernet/refspec.py`; Phase 1 adds e.g.
  `hypernet/emod.py`).  Whether the Phase 2 function (`interpolate_ed_to_l`)
  goes in a flat module or in `hypernet/wiggles/` is JXP's call in Setup Q&A.

## Prompts

### Setup

1. Read the plan sections, gate reports and Phase 0/1 modules listed in
   Context, and this doc.  Put your questions in Q&A: module layout (flat vs
   subpackage), the snippet's location (proposal `hypernet/snippets/`), which
   G0/G1 outcomes apply, whether to target Kevin's pinned commit or the
   current default branch of `hypernets_processor`, and the obsarray variable names to target.  Do not write code
   yet.  Log your work.

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
   equations.  Extend `test_edinterp.py`: the two machine-precision special
   cases remain, plus shape, NaN, edge and 2-D tests, and a timing test that
   a 1536 → 1538 call runs under a stated budget.  Log your work.

4. **SRF self-calibration.**  `calibrate_srf(wav, spectrum, lines=LINES,
   channel="E")` → `SRFModel` (centroid offset and FWHM(λ) quadratic with
   covariance), built on `fit_lines` and `fit_fwhm_model`, with a quality
   flag for low blue signal (the low-sun risk in plan §9).  Test: an
   `observe()`d synthetic spectrum with a known FWHM(λ) is recovered within
   its reported uncertainty; a noisy one returns the flag.  Log your work.

5. **Emod cache.**  `emod_for(sza, pwv_mm, ozone_du, pressure_hpa)` returning
   the nearest cached entry on a stated grid (SZA bins of 5°, a few water
   vapour columns), populated on demand under `ref/emod/`, plus cached
   E-SRF and L-SRF convolutions per `SRFModel` id.  Test: two calls in the
   same bin return the same array; a cache miss builds and stores.  Log your
   work.

6. **Per-instrument SRF table.**  A committed `hypernet/data/srf_table.csv`
   (instrument, channel, cal period, quadratic coefficients, centroid offset,
   uncertainties, n_sequences) written by `hypernet/wiggles/phase2_srf_table.py` from
   Phase 0b's `srf_models.parquet` (or from the VEIT sample alone if the data
   have not arrived, flagged as provisional), with `load_srf_table()` and
   `srf_for(instrument, channel, date)`.  If G0(b) failed, the default path
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

### Gate

10. **Phase 2 exit check.**  Against plan §7 items 1, 2 and 4 and the G0/G1
    verdicts: list each deliverable, its module, its tests and the
    `pytest -q` count; confirm the machine-precision special cases; state the
    default `method` and SRF source the snippet ships with, and what remains
    provisional if the data have not arrived.  Write a short report (proposed
    location `claude_prompts/wiggles/phase2_exit.md`; JXP to confirm -- not
    under `docs/`).  Propose the corresponding edits to
    `docs/wiggles_planning.md` §4 and §7 in Q&A and apply them only if JXP
    approves.  Log your work.

## Q&A

## Logs

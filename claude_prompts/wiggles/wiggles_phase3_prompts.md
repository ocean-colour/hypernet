# WATERHYPERNET -- Wiggles Phase 3: In-situ test

## Goal

On every requested sequence, recompute ρw from L1A with the three methods
(`linear`, `ruddick2023`, `srf`), keeping the L1C/L2A Lw and skyglint terms
unchanged, and judge the result by ρw'' at and away from the Fraunhofer and
telluric lines, by corr(ρw'', Ed''/Ed), and by whether real features (the
683 nm fluorescence peak, the pure-water shoulders above 600 nm) survive.
The PANTHYR system at VEIT is a cross-system L2 control with different wiggle
physics.  The phase ends with the technical note.

**Gate G3** (plan §4): the line-region ρw'' excess is reduced by the amount
Phase 1 predicts for the measured SRF difference, on most sequences and
instruments, without a detectable change in the control regions.  If the
residual after correction still correlates with Ed''/Ed, the SRF model
(Gaussian, quadratic FWHM) is revisited before the tech note is written.

## Conventions

- `ocean14`.  Run via `conda run -n ocean14 python ...`; `conda activate` fails
  non-interactively.
- **Model.** Use Opus 5.5 (`claude-opus-5-5`) for this work, including any
  subagents (pass `model: opus`).
- **JXP runs git.**  Claude does not run any state-changing git command.
- Run scripts **from the repository root** so `hypernet` imports resolve.
- Run `pytest -q` after each step where relevant.  21 tests today (more after
  Phases 0-2); the archive-dependent ones skip themselves when `$OS_COLOR` is
  not mounted.
- **Tier 2** -- every analysis step here reads `$OS_COLOR`.  Do **not** unset
  it.
- Any calculation goes into a script on disk, not into the chat.
- Data intermediates go **outside** the repo under
  `$OS_COLOR/hypernet/wiggles/phase3/`; only figures and small tables are
  committed.
- Every `.md` in `docs/` is published.  Gate reports stay out of `docs/`; the
  technical note is meant to be published, but its destination is JXP's
  decision (task 10).
- Ask questions in the Q&A section below; log completed work under `## Logs`.

## Context

### Where things are

- **Plan:** `docs/wiggles_planning.md` §2.3 (the VEIT baseline numbers:
  corr(ρw'', Ed''/Ed) = −0.56, rms ρw'' 1.9 × 10⁻⁴), §4 "Phase 3", §5
  (which sequences, water types: VEIT/GAIT clear, BEFR/THFR dark, MAFR
  turbid), §7 item 3 (the tech note), §9 (Ring effect, real high-ρw''
  features near lines).
- **Gate reports and exit check:** `claude_prompts/wiggles/gate_G0.md`,
  `gate_G1.md`, `phase2_exit.md` (locations as confirmed there).  Phase 1's
  `hypernet/wiggles/phase1_prediction.csv` is the yardstick for G3.
- **Code (Phases 0-2):** the readers (default `hypernet/whn_l1a.py`), the
  SRF module and `srf_table.csv` / `calibrate_srf`, `interpolate_ed_to_l`
  with `punpy` uncertainties, `emod_for`.  Metric helpers from
  `wavecal/sanity_checks_veit.py` (`hf_power`, `d2`) and
  `hypernet/wiggles/phase1_twin.py` (line/away masks, h = 1 and 5 nm second
  differences) should be promoted to a shared module (default
  `hypernet/wiggle_metrics.py`) in task 3 rather than copied.
- **Data:** the delivered request under `$OS_COLOR/WATERHYPERNET/Wavelengths/`,
  indexed by Phase 0b at `$OS_COLOR/hypernet/wiggles/phase0/request_index.parquet`
  (site, instrument, cal period, SZA, and the sky index `ld_ed_750` /
  `ed_cv_750` / `sky` of Phase 0 task 8e: Ld(750)/Ed(750), Ruddick et al.
  2006 eqs. 23-24, clear < 0.05).  Ld = vza ≥ 90.  L2A variables to reuse
  unchanged: `water_leaving_radiance`, `rhof`, `downwelling_radiance`,
  `upwelling_radiance`, `irradiance`, `reflectance`, `reflectance_nosc`,
  `epsilon`.  PANTHYR VEIT_P L2 files are in the local
  `$OS_COLOR/WATERHYPERNET/RELEASE_2` archive (`hypernet.whn_explore.build_index`,
  `load_spectrum`); no PANTHYR data were requested.
- **Paper:** `context/papers/ruddick2023.pdf` (Fig. 4-6 style plots of ρw'').

### Dependencies

- Needs the full L1A delivery (not yet arrived at the time of writing), Gate
  G0, and the Phase 2 code.  If only part of the request is delivered, run
  everything on what is there and say so in the gate report.
- The Ring check (task 5) may send work back to Phase 2 (non-Gaussian SRF or
  a Ring model); that is a decision for JXP, not a task here.

## Prompts

### Setup

1. Read the plan sections, gate reports and code listed in Context, and this
   doc.  Put your questions in Q&A: how many sequences arrived and whether
   any rows of the request are empty, the figure location (proposal
   `hypernet/wiggles/figs/phase3/`, copied into `docs/figs/` when the tech note takes
   them), the exact recomposition of ρw from L2A terms (task 2), and which
   VEIT_P sequences serve as the control.  Do not run anything yet.  Log
   your work.

### Analysis

2. **Recompute ρw.**  `hypernet/wiggles/phase3_recompute.py`: for every indexed
   sequence, Ed_L with the three methods from the mean L1A_IRR (SRF from the
   table or `calibrate_srf` per G0; `emod_for` at the sequence SZA), then
   ρw_method = π Lw / Ed_L with `water_leaving_radiance` from L2A unchanged
   (so the skyglint term ρ_f·Ld is untouched), and likewise `reflectance_nosc`.
   Check that `linear` reproduces the L2A `reflectance` to the processor's
   precision on every sequence and flag those that do not.  Write one parquet
   per sequence under `$OS_COLOR/hypernet/wiggles/phase3/rhow/` plus a
   manifest.  Log your work.

3. **Metrics.**  Promote the metric helpers into `hypernet/wiggle_metrics.py`
   (`second_difference(y, wav, h)`, `line_mask(wav, lines, half=5.0)`,
   `rms_in_out(y2, mask)`, `corr_with_ed(rhow2, ed2_over_ed)`) with tests on
   synthetic input.  `hypernet/wiggles/phase3_metrics.py`: per sequence and method,
   rms ρw'' within ±5 nm of the ten lines vs away, for h = 1 and 5 nm;
   corr(ρw'', Ed''/Ed) over 400-700 nm; the O2-A and Hα windows separately.
   Table `$OS_COLOR/hypernet/wiggles/phase3/metrics.parquet`, committed
   summary `hypernet/wiggles/phase3_metrics.csv` by instrument, water type and
   `sky` class (the method gains should not depend on the class; if they
   do, say how, and plot the gain against `ld_ed_750`).  Log your work.

4. **Feature preservation.**  `hypernet/wiggles/phase3_features.py`: the 683 nm
   fluorescence peak (height and ρw'' at 670-700 nm) at VEIT/GAIT and MAFR,
   the pure-water absorption shoulders above 600 nm, and the smoothness floor
   (rms ρw'' away from lines must not fall below the `linear` value).  Table
   `hypernet/wiggles/phase3_features.csv` and a figure.  Log your work.

5. **Ring check.**  `hypernet/wiggles/phase3_ring.py`: for the `srf` residual at each
   line (ρw'' excess after correction), regress against ρ_f·Ld/Lu per
   sequence and line, with `ld_ed_750` as a covariate (cloud raises Ld and
   changes the Ring filling-in).  A significant positive slope means the Ring effect is
   in the residual; report it, do not model it (plan §9 -- model only if it
   scales).  Table and figure.  Log your work.

6. **PANTHYR cross-system control.**  `hypernet/wiggles/phase3_panthyr.py`: from the
   local RELEASE_2 archive, take VEIT_P L2 spectra nearest in time to the
   VEIT_H requested sequences (same day, within one hour when possible),
   compute the same ρw'' metrics on the ~10 nm-FWHM TriOS data, and compare
   with VEIT_H before and after correction.  Table
   `hypernet/wiggles/phase3_panthyr.csv` and a figure.  Log your work.

7. **Measured vs predicted reduction.**  `hypernet/wiggles/phase3_vs_prediction.py`:
   per instrument, the measured reduction of the line-region ρw'' excess for
   `ruddick2023` and `srf` against the Phase 1 prediction at that
   instrument's ΔFWHM(λ) and line depths; the residual correlation with
   Ed''/Ed after correction.  Table `hypernet/wiggles/phase3_vs_prediction.csv` and
   a figure.  Log your work.

8. **Figure set for the tech note.**  `hypernet/wiggles/phase3_figures.py` writing
   under `hypernet/wiggles/figs/phase3/`: ρw and ρw'' before/after for one clear, one
   dark and one turbid sequence (paper Fig. 4-6 style); the by-instrument
   summary; the 683 nm control; the PANTHYR comparison.  Log your work.

### Gate

9. **Gate G3.**  Evaluate the criteria with the numbers from tasks 3-7: the
   reduction against the Phase 1 prediction, per instrument and for "most
   sequences"; any detectable change in the control regions; the residual
   correlation with Ed''/Ed and, if it persists, the recommendation to
   revisit the SRF model.  Write a short gate report (proposed location
   `claude_prompts/wiggles/gate_G3.md`; JXP to confirm -- not under `docs/`).
   Propose the corresponding edits to `docs/wiggles_planning.md` §4 in Q&A
   and apply them only if JXP approves.  Log your work.

### Technical note

10. **Outline.**  Propose in Q&A how `docs/wiggles_planning.md` grows into the
    technical note (plan §7 item 3): section order (diagnosis over the
    sequences, twin experiment, in-situ before/after, uncertainty budget,
    code), which figures and tables from Phases 0-3 go in, what stays in the
    plan page, and the destination file and toctree placement in `docs/` --
    **both are JXP's decisions**.  Do not write the note yet.  Log your work.

11. **Write the technical note.**  After JXP approves the outline and
    destination: write it, copying the chosen figures into `docs/figs/`,
    building the Sphinx site once to check for warnings, and leaving
    authorship and the "to be agreed with Kevin" items as placeholders.
    Update `docs/index.md` only if JXP asked for a toctree entry.  Log your
    work.

12. **Code hand-over.**  Check that every deliverable of plan §7 is in
    `hypernet/` with tests, that `pytest -q` passes, that the snippet's
    docstring tells RBINS what to wire up, and that the SRF table carries its
    provenance (which sequences, which cal period).  List the state in the
    Logs.  Log your work.

## Q&A

## Logs

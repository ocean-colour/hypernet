# CLAUDE.md

This file provides guidance to Claude Code when working in this repository.

## Project Overview

hypernet holds code for the WATER HYPERNET network of automated hyperspectral
radiometers: ingesting and analyzing HYPSTAR water-leaving radiance and
reflectance data, improving wavelength calibration, and producing new estimates
of the uncertainty in the radiance and reflectance measurements.  Metrics and
diagnostic plots are intended to be shared with the community.

Repository layout:

- `hypernet/` — the Python package source: the reusable WATERHYPERNET data
  layer.
  - `whn_explore.py` — archive root resolution, filename indexing, single-
    spectrum reads with the agreed product/units conventions (Rrs = ρw/π),
    pooling, shape clustering and the ~100-spectrum-per-site sample.
    `python -m hypernet.whn_explore 1|2`.
  - `whn_simspec_check.py` — the Similarity-Spectrum over-subtraction scan at
    the dark sites.  `python -m hypernet.whn_simspec_check`.
  - `hypernet/tests/` — pytest tests (`pytest.ini` points `testpaths` here).
- `docs/` — **the Sphinx source root** for the WATERHYPERNET user guide, built
  by Read the Docs (`.readthedocs.yaml` at the repo root, `docs/conf.py`,
  `docs/requirements.txt`).  The guide is Markdown parsed by MyST: `index.md`
  carries the byline and the toctree, and the eleven page files follow the
  reading order set in `claude_prompts/howto_prompts.md`.
  - ⚠ **Every `.md` file in `docs/` becomes a published page.**  One left out of
    the toctree is still built and still served at a guessable URL, so anything
    not meant for the public site must go in `exclude_patterns` in `conf.py` or
    live outside `docs/`.  This is why the review correspondence sits in
    `correspondence/`.
  - `whn_figures.py` writes `figs/*.png` + `summary_table.{csv,md}`
    (`summary_table.md` is excluded from the build).  `make_pdf.py` is
    **deprecated** — superseded by the RTD build, macOS/Chrome-only, and its
    default target no longer exists.  Generated PDFs are gitignored.
- `correspondence/` — `kevin_comments.md` and `respond_to_kevin.md`, the review
  exchange with Kevin Ruddick.  Deliberately outside the documentation build;
  the response has not been sent.
- `claude_prompts/` — prompts and task definitions that drive this work.
  `start_up_prompts.md` is the bootstrapping doc; read the relevant prompt doc
  before acting, and do the numbered task you were pointed at, not the whole file.
- `.claude/settings.json` — committed permission policy (copied from `IOPtics`).
  Read-only git is allowed; `git push`/`commit`/`reset`/`rebase` are denied.
- `.claude/skills/` — prompt-behavior skills copied from `IOPtics`:
  `critical-partner` (constructive disagreement) and `grill-me` (one-question-
  at-a-time design interview).

Packaging follows the house pattern in the sibling repos: `requirements.txt` +
`setup.py` (not `pyproject.toml`), `pytest.ini` at the root, and a snake_case
package directory named after the repo.

## Working Conventions

- **Git:** The user (J. Xavier Prochaska) performs all git commands (add,
  commit, push, etc.).  Do not run any git command that changes repository
  state.  Read-only git commands (e.g. `git status`, `git diff`, `git log`) are
  fine when helpful.
- **Calculations:** If you do any calculation, generate it as a Python script
  and write it to disk so that it can be added to the repository.  Do not
  perform one-off calculations only in memory or in the chat.
- **Python environment:** If you need to run Python, use the `ocean14` conda
  environment (e.g. `conda run -n ocean14 python script.py`).
- **Running the scripts:** from the repository root, so that `hypernet` imports
  resolve (`python docs/whn_figures.py`, `python -m hypernet.whn_explore 1`).
  `pip install -e .` removes the need for the `sys.path` bootstrap at the top of
  `docs/whn_figures.py`.
- **Data artifacts:** parquet/npz intermediates are written outside the repo,
  under `$OS_COLOR/hypernet/whn_explore`; only figures and tables are committed.
- **Logging:** Record completed work under the `## Logs` section of the prompt
  doc you were working from, dated, including what you learned about the repo.

## Related Repositories

- **IOPtics:** Related ocean optics code lives on this computer at
  `/Users/xavier/Oceanography/python/IOPtics`.  It is also the upstream source
  for this repo's `.claude/skills/` and `.claude/settings.json`.

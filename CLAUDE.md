# CLAUDE.md

This file provides guidance to Claude Code when working in this repository.

## Project Overview

hypernet holds code for the WATER HYPERNET network of automated hyperspectral
radiometers: ingesting and analyzing HYPSTAR water-leaving radiance and
reflectance data, improving wavelength calibration, and producing new estimates
of the uncertainty in the radiance and reflectance measurements.  Metrics and
diagnostic plots are intended to be shared with the community.

Repository layout:

- `hypernet/` — the Python package source.
  - `hypernet/tests/` — pytest tests.
- `context/WHN/` — the WATERHYPERNET Release-2 exploration, moved over wholesale
  from the `IOPtics` repo (see `claude_prompts/explore_prompts.md`, Prompts/Move).
  `whn_explore.py` (indexing/sampling), `whn_figures.py` (summary table +
  `figs/`), `whn_simspec_check.py` (similarity-spectrum over-subtraction scan)
  and `test_whn_explore.py` — standalone scripts that import each other flat, so
  they must stay in one directory; run them from inside it.  `WATERHYPERNET.md`
  and `respond_to_kevin.md` are outward-facing documents; `make_pdf.py` renders
  them and **shells out to a hard-coded macOS Google Chrome path**.  Generated
  PDFs are gitignored; parquet/npz intermediates stay outside the repo under
  `$OS_COLOR`.
- `pytest.ini` — `testpaths = hypernet/tests context/WHN`.
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
- **Logging:** Record completed work under the `## Logs` section of the prompt
  doc you were working from, dated, including what you learned about the repo.

## Related Repositories

- **IOPtics:** Related ocean optics code lives on this computer at
  `/Users/xavier/Oceanography/python/IOPtics`.  It is also the upstream source
  for this repo's `.claude/skills/` and `.claude/settings.json`.

# Getting started

## Goals

This repository will house code for the WATER HYPERNET network of automated
hyperspectral radiometers -- ingesting and analyzing HYPSTAR water-leaving
radiance and reflectance data, and related work.  
We are also likely to include code to improve wavelength calibration and
new estimates of the uncertainty in the radiance and reflectance measurements.  We will also generate metrics and plots
diagnostics to share with the community.

## Prompts

1. Read this file.  Execute the 1st task under "Claude/CLAUDE.md file"
2. Read this file.  Execute the 1st task under "Claude/Skills"
3. Read this file.  Execute the 1st task under "Claude/Settings"
4. Read this file.  Execute the 1st task under "Basic start up"

## Claude

### CLAUDE.md file

1. Please generate a basic CLAUDE.md file for this project.  Have it indicate:

- I will perform git commands
- Add to the CLAUDE.md file:  If you do any calculation, generate it as a python script and write it to disk so that I can add it to the Repository.
- Add to the CLAUDE.md file:  If you need to run Python, use the "ocean14" conda environment.
- Add to the CLAUDE.md file:  Related optics code lives in the sibling repository /Users/xavier/Oceanography/python/IOPtics

### Skills

1. Copy over the skills/ files from the `IOPtics` repository.

### Settings

1. Copy over the settings.json file from the `IOPtics` repository.  Copy the
   policy only -- prune the accumulated path-specific allow entries that refer
   to IOPtics' own scratchpad paths and test files.

## Basic start up

1. Generate the basic files that one needs for a Python GitHub repository, e.g. a file for dependencies.  Examine the other Repositories in Oceanography/python to see how I tend to organize things.  Also, make a suggestion of which of their settings.json files to copy into this one.

## Logging

The "Logs" section will record Claude's work.  Please use the following format:

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

...

## Logs

### 2026-09-18 (Generate the basic CLAUDE.md)

Executed prompt #1: generated `CLAUDE.md` at the repo root.

What was done:

- Read the sibling `CLAUDE.md` files in `IOPtics`, `OETHER`, and
  `cugn-climatology` first, to match house style rather than invent a format.
  Used `cugn-climatology`'s shape (Project Overview / Working Conventions /
  Related Repositories) since it is the only sibling that also carries a
  cross-repo pointer, which this repo needs.
- Recorded the four directives requested by the task: git is the user's alone
  (read-only git permitted), calculations become scripts on disk, Python runs in
  the `ocean14` conda environment, and related optics code lives at
  `/Users/xavier/Oceanography/python/IOPtics`.
- Wrote the Project Overview from the `## Goals` block above, including the
  wavelength-calibration and uncertainty-estimate aims.
- Added a Logging convention line so a session that never opens a prompt doc
  still knows work gets logged.

What was learned:

- The repo is currently near-empty: `LICENSE`, `README.md`, `.gitignore`,
  `claude_prompts/`, and now `CLAUDE.md`.  There is no Python package, no
  `requirements.txt`/`setup.py`, and no `.claude/` directory yet — those arrive
  with prompts #2-#4.  The layout section of `CLAUDE.md` therefore lists only
  `claude_prompts/` and carries a note to update it once the package exists,
  rather than documenting directories that do not exist.
- The `bin/claude_prompts/` referenced when bootstrapping this repo is
  `~/bin/claude_prompts/`, not a `bin/` inside this repo.
- The prompt doc here is named `start_up_prompts.md`, whereas the other seven
  repos use `start_up.md`.
- No git commands were run.

### 2026-09-18 (Copy skills from IOPtics)

Executed prompt #2: copied `.claude/skills/` from the `IOPtics` repository.

What was done:

- Created `.claude/skills/` and copied both skills from
  `/Users/xavier/Oceanography/python/IOPtics/.claude/skills/`:
  - `critical-partner/SKILL.md` (55 lines) — constructive-disagreement mode:
    surface an untested assumption before agreeing, argue the strongest
    opposing case, retreat only on new evidence, review weaknesses first, name
    emotional investment, do not manufacture flaws, end with a question rather
    than a summary.  No flattery openers.
  - `grill-me/SKILL.md` (10 lines) — interview one question at a time down each
    branch of the design tree, always offering a recommended answer, and
    explore the codebase instead of asking what the codebase can answer.
- Read both files before copying and confirmed they contain nothing
  IOPtics-specific, so a verbatim copy was correct — no adaptation needed.
- Verified with `diff -r` that the copies are identical to the IOPtics originals.
- Updated the layout section of `CLAUDE.md` to list `.claude/skills/`.

What was learned:

- IOPtics has exactly two skills, and both are repo-agnostic prompt-behavior
  skills rather than project tooling.  That is why the one-line "copy over the
  skills/ files" instruction works unchanged across repos.
- `.claude/` did not exist in this repo before this task; it now contains only
  `skills/`.  `settings.json` arrives with prompt #3.
- No git commands were run.

### 2026-09-18 (Copy settings.json from IOPtics, policy only)

Executed prompt #3: created `.claude/settings.json` from the `IOPtics` policy.

What was done:

- Copied the full permission policy: read-only shell (`ls`, `cat`, `find`,
  `grep`, `rg`, `head`, `tail`, `wc`, `sort`, `which`, `env`, `pwd`, `echo`),
  benign file ops (`mkdir`, `touch`, `cp`, `mv`), read-only git (`status`,
  `diff`, `log`, `show`, `branch`), conda/python/pip/pytest/jupyter including
  `conda run -n ocean14`, and the four publisher `WebFetch` domains (Optica,
  ASLO/Wiley, doi.org, Taylor & Francis).
- Kept the deny list verbatim — `sudo`, `rm -rf /`, `rm -rf ~`, and
  `git push` / `git commit` / `git reset` / `git rebase`, so the git rule is
  enforced by permissions and not only by prose in `CLAUDE.md`.  Kept
  `ask` = `rm:*`.
- Pruned 8 accumulated path-specific allow entries from the IOPtics tail:
  two `curl` calls with a hardcoded DOI and a hardcoded bibliographic query,
  two long `python3 -c` Crossref JSON parsers, two absolute-path pytest
  invocations (one `/home/xavier/miniconda3`, one `/Users/xavier/miniforge3`),
  an absolute-path `python -c "import bing"` check, and a `grep` of
  `IOPtics/reports/gloria_fits_report.md`.
- Added one entry, `WebFetch(domain:api.crossref.org)`, to preserve Crossref
  literature checking as a first-class activity without carrying over the
  hardcoded-DOI `curl` cruft.  This is the one deliberate departure from a
  straight copy.
- Validated the result parses as JSON: 38 allow, 7 deny, 1 ask (IOPtics has 44
  allow).

What was learned:

- The two pruned pytest entries point at *different machines*
  (`/home/xavier/miniconda3` = Linux, `/Users/xavier/miniforge3` = Mac), so the
  IOPtics allow-list has accreted cross-machine cruft.  The generic
  `Bash(pytest:*)` and `Bash(conda run -n ocean14:*)` entries already cover
  these, so nothing was lost by dropping them.
- `.claude/settings.local.json` is **not** ignored by any repo's own
  `.gitignore`; it is ignored globally by
  `~/.config/git/ignore`, which contains `**/.claude/settings.local.json`.
  Verified with `git check-ignore -v` in this repo and five siblings.  No
  `.gitignore` entry is needed here.
- `.claude/` in this repo now holds `settings.json` and `skills/`, matching
  IOPtics minus the machine-local `settings.local.json`.
- No git commands were run beyond read-only `git check-ignore`.

### 2026-09-18 (Basic Python repository files)

Executed prompt #4: generated the basic files for a Python GitHub repository,
after surveying the seven sibling repos in `Oceanography/python`.

What was done:

- Surveyed `IOPtics`, `OETHER`, `PAB`, `cugn-climatology`, `retrieve-or-bust`,
  `victor-class`, and `cugn` to copy the house layout rather than apply generic
  best practice.  Used `cugn-climatology` as the closest model: it is the
  smallest, most recent, and the only other repo whose `setup.py` documents a
  sibling-package dependency.
- Created the package: `hypernet/__init__.py` (empty, as in the siblings),
  `hypernet/tests/__init__.py`, and `hypernet/tests/test_import.py` — a smoke
  test mirroring `cugn_climatology/tests/test_import.py` verbatim in style.
- Created `pytest.ini` with `testpaths = hypernet/tests`,
  `python_files = test_*.py`, `addopts = -ra`, matching cugn-climatology.
- Created `requirements.txt` and `setup.py` with matching dependency lists, per
  the "keep this list in sync" comment the siblings carry.  `setup.py` follows
  the house `setup_keywords` dict idiom exactly, including the trailing
  `if os.path.isdir('bin')` scripts block: name `hypernet`, version
  `0.0.dev0`, BSD (confirmed against the repo's own BSD 3-Clause `LICENSE`),
  `python_requires >= 3.12`, url `https://github.com/ocean-colour/hypernet`
  (read from `git remote`, not guessed).
- Expanded `README.md`, keeping the original two lines verbatim at the top and
  adding What it does / Installation / Tests / Related sections in the style of
  the `IOPtics` and `cugn-climatology` READMEs.
- Updated the layout section of `CLAUDE.md` for the new package and removed the
  placeholder note left there by prompt #1.
- Ran `conda run -n ocean14 python -m pytest -q`: 1 passed.

Answer to the second half of the task ("suggest which settings.json to copy"):
`IOPtics` — already done under prompt #3.  It is the upstream source of truth
for `.claude/` across these repos, and its policy (read-only shell + git,
`conda run -n ocean14`, publisher `WebFetch` domains, deny
`git push`/`commit`/`reset`/`rebase`) is the one the others inherit.

What was learned:

- The house packaging convention is `requirements.txt` + `setup.py`, **never**
  `pyproject.toml` — confirmed across all seven repos.  `setup.py` is the same
  `setup_keywords` dict template in each, differing only in name/description/
  url/install_requires.
- The reliable core is `<package>/`, `claude_prompts/`, `CLAUDE.md`,
  `README.md`, `LICENSE`, `requirements.txt`, `setup.py`, `.gitignore`.
  `pytest.ini` appears in PAB, cugn-climatology, and victor-class but not
  IOPtics or OETHER.  Optional extras vary widely: `context/`, `docs/`,
  `reports/`, `papers/`, `notebooks/`, `data/`, `runs/`.
- The `hypernet` GitHub remote is `ocean-colour/hypernet` — the same org as
  IOPtics, ocpy, and BING.
- **Open question for the user:** IOPtics pulls the sibling packages `ocpy` and
  `BING` from GitHub via `requirements.txt`.  I did not add them here, since
  nothing in the goals requires an IOP retrieval engine yet.  If this repo will
  reuse `ocpy` utilities, add
  `git+https://github.com/ocean-colour/ocpy` to `requirements.txt`.
- Dependencies were chosen from the stated goals (netCDF HYPSTAR products,
  calibration, uncertainty) and are a starting point, not a researched set.
- No git commands were run beyond read-only `git remote -v`.

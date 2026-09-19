# WATERHYPERNET -- the User HOWTO

## Goal

Turn `docs/WATERHYPERNET.md` -- currently an independent data-user's summary
written for the WATERHYPERNET community -- into a published **User HOWTO** for
the archive, exposed via **readthedocs**.

## Conventions

- `ocean14`.  Run via `conda run -n ocean14 python ...`; `conda activate` fails
  non-interactively.
- **JXP runs git.**  Claude does not run any state-changing git command.
- Run scripts **from the repository root** so `hypernet` imports resolve
  (`python docs/whn_figures.py`, `python -m hypernet.whn_explore 1`).
  `pip install -e .` removes the need for the `sys.path` bootstrap at the top of
  `docs/whn_figures.py`.
- Run `pytest -q` after each step where relevant.  19 tests today; the archive-
  dependent ones skip themselves when `$OS_COLOR` is not mounted.
- **Tier 2** -- the figure and scan steps need the `$OS_COLOR` data tree mounted.
  Do **not** unset `$OS_COLOR`.
- Any calculation goes into a script on disk, not into the chat.
- Ask questions in the Q&A section below; log completed work under `## Logs`.

## Context

### Where things are (after the Move, 2026-09-19)

- `hypernet/` -- the reusable data layer.
  - `whn_explore.py` -- `whn_root()`, `out_root()`, filename indexing,
    `load_spectrum`, pooling, shape clustering, the ~100-per-site sample.
    Stages: `python -m hypernet.whn_explore 1` (index) and `... 2` (pool,
    cluster, sample).
  - `whn_simspec_check.py` -- the Similarity-Spectrum over-subtraction scan.
    `python -m hypernet.whn_simspec_check [--sites all]`.
  - `tests/test_whn_explore.py` -- 18 tests over both modules, two tiers.
- `docs/` -- the document and the code that exists only to produce it.
  - `WATERHYPERNET.md` -- **the draft HOWTO**, 509 lines, 15 numbered sections
    plus a "How this was produced" appendix.
  - `whn_figures.py` -- writes `figs/*.png` and `summary_table.{csv,md}` here.
  - `make_pdf.py` -- Markdown -> PDF.  **macOS + Google Chrome only** (hard-coded
    path at line 26); raises a clean `FileNotFoundError` elsewhere.  Generated
    PDFs are gitignored.
  - `kevin_comments.md`, `respond_to_kevin.md` -- the review correspondence.
  - `simspec_flagged.csv` -- the 30 flagged spectra, cited by §9.2.
- Intermediates (parquet/npz, ~18 MB) live **outside** the repo at
  `$OS_COLOR/hypernet/whn_explore`; the archive itself is
  `$OS_COLOR/WATERHYPERNET/RELEASE_2`.

### What the document currently is

Authored "J. Xavier Prochaska (University of California, Santa Cruz)", dated
18 September 2026, and framed as *"an independent exploration ... intended to
complement -- not replace -- the official release notes"*.  Sections 1-13 are
descriptive (contents, scope, sites, layout, file contents, units, grids,
uncertainty, QC, product choice, geometry, what the spectra look like, loader
hazards); §14 is "Implications for IOP retrieval work"; §15 is licence and
citation.  Every number in it was measured from the archive.

**The tension to resolve first:** it was written as a *findings* document for the
WATERHYPERNET community and reviewed as one by Kevin Ruddick.  A User HOWTO is a
different genre -- task-ordered rather than survey-ordered, and it would normally
open with "how do I load a spectrum correctly" rather than "what the release
contains".  Most of the raw material is right; the ordering, the voice and the
audience are the open questions.

### Live threads to pick up

- **Kevin's response has not been sent.**  `docs/respond_to_kevin.md` now points
  at `github.com/ocean-colour/hypernet`.  It promises him a stable repo link,
  offers right of reply on §9.2 and on the "Hazards for an automated loader"
  list, and takes up his offer of an "exploitation tools" link.
- **The deferred docstring sweep** (Q&A M6 in `explore_prompts.md`, held over to
  this effort).  Stale internal references, none of them functional:
  - `hypernet/whn_explore.py:8` -- cites `claude_prompts/waterhypernet_prompts.md`
    (this repo has `explore_prompts.md`).
  - `hypernet/whn_explore.py:28` -- "mirroring the `build_v1.py` driver pattern",
    an IOPtics convention with no counterpart here.
  - `hypernet/tests/test_whn_explore.py:3-4` -- points at
    `ioptics/tests/conftest.py` as the testing convention.
- **§14 "Implications for IOP retrieval work"** is the one section written for a
  different audience (it exists because this exploration started inside IOPtics).
  It may belong in the IOPtics docs rather than in a WATERHYPERNET user HOWTO.
- **Nothing readthedocs exists yet.**  `docs/` is flat Markdown: no Sphinx, no
  `conf.py`, no `.readthedocs.yaml`, and `docs/` is not a package.
  `requirements.txt` has no Sphinx/MyST/MkDocs entry.

### Numbers worth keeping straight

56,669 files, 9.8 GB, 11 sites (7 HYPSTAR `_H`, 4 PANTHYR `_P`; VEIT has both).
Products are **ρw, not Rrs** -- `Rrs = ρw/π`.  `reflectance` everywhere except
LPAR_H, MAFR_H and O1BE_P, which use `reflectance_nosc`.  Over-subtraction at the
dark sites: **30 of 8,947 (0.335 %)**, all `quality_flag = 0`, and that is a
**lower bound** (a sign test only).

## Prompts

### HOWTO

1. Read the Context above and `docs/WATERHYPERNET.md` in full.  Then propose how
   it should be restructured as a User HOWTO: the section order, what a
   task-oriented reader needs first, what to cut or move (§14 in particular),
   and how the "independent summary" voice should change now that it is
   documentation rather than a findings report.  **Do not edit the document
   yet** -- put the proposal and your questions in the Q&A/HOWTO section below. Use Opus 5.
   Log your work.

2. I have answered your questions.  Restructure `docs/WATERHYPERNET.md`
   accordingly and do the deferred docstring sweep listed in Context.  Keep every
   measured number intact, and tell me about any claim the restructure would
   strand.  Run `pytest -q`.  Log your work.

3. Set up readthedocs for `docs/`: the builder, the config, the theme, and the
   navigation.  Decide with me first whether the source stays Markdown (Sphinx +
   MyST) or moves to reStructuredText, and what happens to `make_pdf.py` once
   readthedocs can build a PDF itself.  Log your work.

4. Finalise and send-ready the Kevin correspondence: check
   `docs/respond_to_kevin.md` against the restructured document, regenerate both
   PDFs, and tell me exactly what is still outstanding on his side.  Log your
   work.

## Q&A

### HOWTO

**Defects I found while reading, which want fixing whatever we decide.** Not
questions -- flagging them so they are not lost in the restructure.

- **Line 75 is a broken cross-reference.** §2 ends "...that is discussed in
  [§11](#11-implications-for-iop-retrieval-work)", but IOP retrieval is **§14**;
  §11 is "Geometry and ancillary data". Both the label and the anchor are stale
  from an earlier draft, so the link is dead in the Markdown *and* in the PDF
  that is queued to go to Kevin.
- **"18 tests" is now wrong, and that one is mine.** The appendix run-block says
  `pytest -q  # 18 tests`. When I retargeted that block to the new repo last
  session I changed the command from `pytest -q test_whn_explore.py` to a bare
  `pytest -q`, which collects **19** (18 WHN + the package import test). Either
  restore the narrow path or change the number.
- **The two CSVs in `docs/` are never linked from the document.**
  `simspec_flagged.csv` is exactly the evidence behind §9.2's table of 30, and
  `summary_table.csv` is the machine-readable form of §3's site table. Both
  deserve links.
- **Five hard-coded `§N` cross-references** (lines 75, 265, 271, 367, 503) and 15
  manually numbered headings. Every one of them breaks on reorder, including the
  appendix's "Scope of each number" paragraph, which cites §9.2 by number.
- **Our own loader walks into hazard 5, and I only found this by testing it.**
  `_filled` is careful about the *crash* half (it casts to float before
  `filled`, so the integer `TypeError` never fires), but nothing recovers the
  masked value, so a PANTHYR `quality_flag` of 0 -- a legitimate **pass** --
  comes back as `nan`. Measured on a real file:

  ```
  O1BE_P  quality_flag via load_spectrum : nan     raw, masking off : 0
  THFR_H  quality_flag via load_spectrum : 0.0
  ```

  So anyone filtering `quality_flag == 0` through our loader silently discards
  **every PANTHYR spectrum -- 12,080 of them, 21% of the archive** -- while
  keeping all the HYPSTAR ones. The document itself is not wrong here: §9 says
  the 440-file check was "read with masking off so fill and zero are
  distinguishable", which is exactly right. It is `load_spectrum` that is lossy,
  and it is a one-line fix (`set_auto_mask(False)` for the flag, or special-case
  the PANTHYR variable). Worth doing before the HOWTO tells people to use it --
  and it makes an honest worked example for the Traps page.

**The structural read.** Three things stand out once you read it end to end as a
HOWTO rather than as a report:

1. **It has no code.** Eight fenced blocks, none of them Python: a directory
   tree, two filename grammars, the π identity and the bash reproduction block.
   The first question a user has is "how do I open one of these files and get
   Rrs", and the document never answers it in a form you can run -- even though
   `hypernet.whn_explore.load_spectrum(path, site)` does precisely that, and in
   doing so already handles hazards **1, 2 and 3** from §13, plus the crash half
   of hazard 5 (see the loader bug below for the half it does not).
2. **§13 "Hazards for an automated loader" is a recap, and it should be the
   spine.** Its ten items re-state points already made in §5 (two schemas), §6
   (ρw not Rrs), §7 (grids vary), §9.1 (negatives) and §11 (fill values). In a
   survey a closing checklist is good practice; in a HOWTO the warning belongs at
   the step where the reader would hit it, with the fix attached.
3. **The order is "what this archive is", not "what you are trying to do".**
   Everything a user needs to load data correctly -- units (§6), product choice
   (§10), schema differences (§5), grids (§7) -- is spread across four sections
   separated by the site table, the map and the coverage discussion.

**Proposed order.** Task-ordered, with the current sections in brackets:

| # | Page | From |
|---|---|---|
| 1 | **About the archive** -- what WATERHYPERNET is, Release 2 in numbers, the two systems, how to get it | §1, practical half of §15 |
| 2 | **Quickstart: read one spectrum** -- layout and filename grammar, then ~10 runnable lines to Rrs on the native grid | §4 + new code |
| 3 | **Getting Rrs right** -- ρw ÷ π (and σ ÷ π), which product at which site, the std-variable name split | §6, §10, part of §5 |
| 4 | **Reading many files** -- the two schemas, grids vary by serial, no CF-decodable time, site+time is not a unique key | §5, §7, hazards 1/3/4/7/8 |
| 5 | **Traps** -- PANTHYR `_FillValue = 0`, the integer-fill `TypeError`, the metadata bugs, the release-window overrun | remainder of §13 |
| 6 | **Quality: what to filter and what not to** -- `quality_flag`, negatives are deliberate, low-sun screening, what σ does and does not include | §9, §9.1, §8, geometry half of §11 |
| 7 | **Choosing your data** -- site table, map, coverage unevenness, what the spectra look like | §3, §12 |
| 8 | **Known issue: Similarity-Spectrum over-subtraction at dark sites** | §9.2 |
| 9 | **Licence, citation and acknowledgement** | §15 |
| -- | **Appendix: how these numbers were produced** | current appendix |

That moves §6 and §10 from the middle to the front (they are the two things a
user is most likely to get silently wrong), turns §13 from a recap into the
running commentary of pages 3-5, and leaves the descriptive material (§3, §12) as
a "help me choose" page rather than the opening.

**H1 -- Who is this for: the archive, the package, or both?** Today it documents
the *archive* and never mentions that this repo can read it. readthedocs for a
Python repo normally documents the package. Options: (a) archive guide only;
(b) two top-level sections -- "Using the archive" (this document) and "hypernet
API" (autodoc from the existing docstrings, which are already numpydoc);
(c) fully merged, so the HOWTO's examples *are* the package tutorial. **I
recommend (b)**, with the quickstart on page 2 using `load_spectrum` and also
showing the plain-`netCDF4` equivalent -- that makes the loader discoverable
without making a WATERHYPERNET user depend on our package to follow the guide.

**H2 -- How much of the "independent summary" voice survives?** The byline, the
date, "intended to complement -- not replace -- the official release notes",
"any remaining errors are mine" and the acknowledgement of Kevin's review. This
framing was negotiated with him and the document is authored under your name for
the WATERHYPERNET community, so I do not think it should simply vanish. Options:
(a) keep it at the top as now; (b) move it to an "About this guide" page and let
the task pages speak in plain instructional second person; (c) drop it.
**I recommend (b)** -- but it is your name on it, so it is your call.

**H3 -- §14 "Implications for IOP retrieval work": move, keep or cut?** It exists
because the exploration began inside IOPtics, and it is the one section addressed
to a different reader. **I recommend moving it to the IOPtics docs** with a
one-line pointer here -- but it does **not** split cleanly, and that is the part
worth your attention: its closing paragraph is Kevin's argument that flow
cytometry, turbidimeters and fluorimeters would be more useful additions to the
network than IOPs. That is a statement about **WATERHYPERNET**, not about
IOPtics, and I think it belongs here, in the scope discussion (§2), rather than
travelling to the other repo. Do you agree with splitting it that way?

**H4 -- §9.2 is a finding, and publishing it is a decision.** The document
promises Kevin right of reply on exactly this section and on the loader-hazards
list, and `respond_to_kevin.md` has not been sent. Putting the HOWTO on
readthedocs publishes both. Options: (a) publish in full as a "Known issues"
page; (b) reduce the published page to the user-facing caveat plus the sanity
check, keeping the full analysis in the correspondence; (c) hold the section out
of the first published build until he has replied. **I recommend (b) for the
shape and (c) for the timing** -- which makes prompt 4 (the Kevin thread) a
blocker on prompt 3 (going live), not a follow-up. Do you want the prompt order
in this file changed to reflect that?

**H5 -- One long page or nine?** Nine pages read better on readthedocs; one long
page renders better as the PDF you promised Kevin as a citable frozen artifact.
Sphinx can do both from one source. **I recommend multi-page with a PDF build**,
which would also retire `make_pdf.py` and its hard-coded Chrome path -- but that
is properly a prompt-3 decision, so treat this as a heads-up rather than a
question needing an answer now.

**H6 -- Manual section numbers.** I would drop `## N.` from the headings and use
named anchors, so reordering cannot produce another dead link like line 75.
Section numbers would then be generated. Any objection to losing the stable "§9.2"
style references, given Kevin's comments and your response both cite sections by
number?

## Logs

### 2026-09-19 — HOWTO prompt 1 (document read; restructure proposed, nothing edited)

Opus 5, no subagents. `docs/WATERHYPERNET.md` read in full (509 lines, 15
numbered sections + appendix). **No file was edited** — the prompt says propose,
not act — apart from this prompt doc.

- **Proposed a task-ordered structure** of nine pages plus the provenance
  appendix, mapped section by section onto the current ones, in Q&A/HOWTO. The
  two substantive moves: §6 (units) and §10 (product choice) go to the front,
  because they are what a user gets silently wrong; and §13 "Hazards for an
  automated loader" stops being a closing recap and becomes the running
  commentary of the three loading pages, since its ten items already re-state
  §5, §6, §7, §9.1 and §11.
- **Six questions posed (H1-H6)**, each with a recommendation: audience
  (archive guide vs package docs vs both), how much of the "independent summary"
  voice survives, what to do with §14, whether §9.2 can be published before Kevin
  replies, one page vs nine, and dropping the manual section numbers.
- **H3 and H4 are the two I would not decide alone.** §14 does *not* split
  cleanly — its closing paragraph is Kevin's argument about flow cytometry and
  fluorimeters, which is about the network, not about IOPtics, so it should stay
  here even if the rest of §14 leaves. And §9.2 plus the hazards list are exactly
  what the document promised Kevin right of reply on, so publishing to
  readthedocs before he answers would break that promise — which makes prompt 4
  a blocker on prompt 3 rather than a follow-up.
- **Four defects found by reading, one by testing.** Reading: the cross-reference
  at line 75 points at `#11-implications-for-iop-retrieval-work` when that
  section is §14 (dead in the Markdown and in the PDF queued for Kevin); the two
  CSVs sitting in `docs/` are never linked from the text; five hard-coded `§N`
  references will break on any reorder; and the appendix run-block says "18
  tests" when a bare `pytest -q` now collects 19 — that last one is mine, from
  retargeting the command during the Move.
- **The one found by testing is the real one.** I claimed in a first draft of the
  proposal that `load_spectrum` handles hazards 1, 2, 3 and 5, then checked it
  against a real file instead of trusting the docstring. It does not handle 5.
  `_filled` dodges the integer-fill `TypeError` by casting to float before
  `filled`, but nothing recovers the masked value, so PANTHYR's
  `quality_flag = 0` — a **pass** — returns `nan`, while HYPSTAR returns `0.0`.
  Anyone filtering `quality_flag == 0` through our loader drops all 12,080
  PANTHYR spectra (21% of the archive) and keeps every HYPSTAR one. The document
  is not wrong here — §9's 440-file census explicitly read with masking off — it
  is our loader that is lossy. Corrected the claim in the proposal and flagged
  the fix as something to land before the HOWTO recommends the function.
- **What the document lacks as a HOWTO, in one line:** no Python. Eight fenced
  blocks, none of them runnable — a directory tree, two filename grammars, the π
  identity and the reproduction block — while the repo contains the loader that
  answers the reader's first question.

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

2. I have answered your questions.  Read my responses and see if you need 
   any further information.  If so, ask another round of questions.
   Use Opus 5.  Log your work.

3. I have answered your Round 2 questions.  Read my responses and see if you need 
   any further information.  If so, ask another round of questions.
   Use Opus 5.  Log your work.

4. Ok, I have answered your Round 3 questions.  Please revise the prompts that  
   follow based on our Discussion and I will execute them in order.
   Use Opus 5.  Log your work.

*Everything settled in Q&A rounds 1-3 is digested below, so prompts 5-8 can be
executed in order without re-reading the discussion.*

**Settled — apply throughout.**

- **Scope.** The guide introduces WATERHYPERNET and gives the basics of
  **accessing, loading and viewing** the data (H1). The `hypernet` package is
  intentionally small for now and gets its own docs later, so leave room in the
  RTD tree for a package section without building one.
- **Voice.** Keep the independent-summary framing — byline, date, "complement,
  not replace, the official release notes", "any remaining errors are mine", and
  the acknowledgement of Kevin's review — **on the landing page only**, not
  repeated per page (H2).
- **Page order** (round 1, agreed): 1 About the archive *(including access)* ·
  2 Quickstart: read one spectrum · 3 Getting Rrs right · 4 Reading many files ·
  5 Traps · 6 Quality: what to filter and what not to · 7 Choosing your data ·
  8 Known issue: Similarity-Spectrum over-subtraction · 9 Licence, citation and
  acknowledgement · Appendix: how these numbers were produced.
- **Numbering.** Drop the manual `## N.` headings and use named anchors, so a
  reorder cannot produce another dead link (H6).
- **The IOP section.** Keep it, but generalised to IOP retrieval rather than
  IOPtics-specific, with a brief mention of the IOPtics effort linking both
  <https://github.com/ocean-colour/IOPtics> and
  <https://ioptics.readthedocs.io/en/stable/> (H3).
- **The SWIR note.** In the over-subtraction section, record that the behaviour
  is known and is being addressed with a SWIR-based approach, attributed as
  *(WATERHYPERNET team, pers. comm., September 2026)*. Scope it to **all HYPSTAR
  sites**, with the caveat that SWIR bands are not distributed in these files
  (H4, R2-Q1, R2-Q2, R3-Q1).
- **Access.** Request via the form at <https://waterhypernet.org/contact/> — name,
  e-mail and affiliation required, purpose "Request for Data access" — after
  which FTP credentials for <https://ftp.waterhypernet.org/> follow by reply. The
  site does not itself document the credential step, so do not imply that it
  does (R2-Q3, R3).
- **PDF.** RTD is the primary documentation source. `docs/make_pdf.py` and the two
  built PDFs stay on disk, marked deprecated (H5, R2-Q5).

**Defects to fix — found in rounds 1-2, none yet corrected.** D5 is handled in
prompt 5, D1-D4 in prompt 6.

- **D1 — the §5 product table is factually wrong** (line 148). It says `reflectance`
  is "corrected with the Similarity Spectrum test (SWIR-based for HYPSTAR in
  Release 2, NIR for PANTHYR)". Per the release notes the **correction is
  NIR-based for both systems**; what is SWIR-based for HYPSTAR is the **QC
  test**. Fix both halves and keep the two ideas distinct wherever they appear.
- **D2 — broken cross-reference** (line 75): `[§11](#11-implications-for-iop-retrieval-work)`
  points at a section that is §14.
- **D3 — "18 tests"** in the appendix run-block: a bare `pytest -q` collects 19.
- **D4 — the two CSVs are never linked** — `simspec_flagged.csv` is the evidence
  behind the over-subtraction table, `summary_table.csv` the machine-readable
  site table.
- **D5 — `load_spectrum` masks PANTHYR `quality_flag`**, returning `nan` where the
  real value is 0 (a pass), so filtering on it discards 21% of the archive.

---

5. **Code first, before the guide recommends the loader.** Fix the
   `quality_flag` masking bug in `hypernet/whn_explore.py` — narrow to that
   variable only, since PANTHYR's angles are genuinely fill and must stay `nan` —
   with a regression test covering one PANTHYR and one HYPSTAR file. Do the
   deferred docstring sweep listed in Context. Add the deprecation note to
   `docs/make_pdf.py`. Run `pytest -q`. Use Opus 5. Log your work.

6. **Restructure the document in place.** Rewrite `docs/WATERHYPERNET.md` into the
   agreed page order as a **single file** — do not split it yet. Apply the
   settled decisions above and fix defects D1-D4. **Keep every measured number
   intact**, and tell me about any claim the restructure would strand or leave
   load-bearing in a new way. Run `pytest -q`. Use Opus 5. Log your work.

7. **Split to readthedocs.** Split the restructured document into the nine pages
   plus appendix, write the toctree and `.readthedocs.yaml`, add the docs build
   requirements and a theme, and leave a slot for a future package-API section.
   Decide with me first whether the source stays Markdown (Sphinx + MyST — my
   recommendation, since every page is already Markdown) or moves to
   reStructuredText. Confirm the figures and both CSVs are carried into the
   build. Use Opus 5. Log your work.

8. **Kevin — deferred; do not act until I say so.** Check
   `docs/respond_to_kevin.md` against the restructured document: the
   over-subtraction numbers, the section references, and the line promising "A
   PDF of the current version accompanies this note", which the RTD decision
   makes stale. Tell me what is outstanding on his side and propose wording.
   Send nothing and regenerate nothing. Use Opus 5. Log your work.

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

>A. Ok, let's do (b) but the code base in this Repo is intentionally small at the
moment.  It will advance, and when it does there will be new and separate docs 
for the package.  This HOWTO is primarily to introduce the user to WATERHYPERNET
and provide a few basics on how to access, load, and view the data.

**H2 -- How much of the "independent summary" voice survives?** The byline, the
date, "intended to complement -- not replace -- the official release notes",
"any remaining errors are mine" and the acknowledgement of Kevin's review. This
framing was negotiated with him and the document is authored under your name for
the WATERHYPERNET community, so I do not think it should simply vanish. Options:
(a) keep it at the top as now; (b) move it to an "About this guide" page and let
the task pages speak in plain instructional second person; (c) drop it.
**I recommend (b)** -- but it is your name on it, so it is your call.

>A. (a)

**H3 -- §14 "Implications for IOP retrieval work": move, keep or cut?** It exists
because the exploration began inside IOPtics, and it is the one section addressed
to a different reader. **I recommend moving it to the IOPtics docs** with a
one-line pointer here -- but it does **not** split cleanly, and that is the part
worth your attention: its closing paragraph is Kevin's argument that flow
cytometry, turbidimeters and fluorimeters would be more useful additions to the
network than IOPs. That is a statement about **WATERHYPERNET**, not about
IOPtics, and I think it belongs here, in the scope discussion (§2), rather than
travelling to the other repo. Do you agree with splitting it that way?

>A. Keep - but make it generic to IOP retrievals and have a brief mention of
the `IOPtics` effort.  Link to both the RTD and GitHub pages for `IOPtics`:
   - https://github.com/ocean-colour/IOPtics
   - https://ioptics.readthedocs.io/en/stable/

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

>A. Let's not worry about responding to Kevin for now.  But do update the text to note that this is a known problem and the Team is addressing it by using a new
approach where there is SWIR data.

**H5 -- One long page or nine?** Nine pages read better on readthedocs; one long
page renders better as the PDF you promised Kevin as a citable frozen artifact.
Sphinx can do both from one source. **I recommend multi-page with a PDF build**,
which would also retire `make_pdf.py` and its hard-coded Chrome path -- but that
is properly a prompt-3 decision, so treat this as a heads-up rather than a
question needing an answer now.

>A. Ok, multi-page is fine.  We no longer need the PDF.  We will use RTD as the primary documentation source.

**H6 -- Manual section numbers.** I would drop `## N.` from the headings and use
named anchors, so reordering cannot produce another dead link like line 75.
Section numbers would then be generated. Any objection to losing the stable "§9.2"
style references, given Kevin's comments and your response both cite sections by
number?

>A. Yes, that is fine.

### HOWTO — round 2

All six answered, and five of them need nothing further: **H1(b)** with the
package docs deferred, **H2(a)** framing kept at the top, **H3** §14 kept but
generalised with the two IOPtics links, **H5** multi-page and no PDF, **H6** manual
numbers dropped. I am ready to build all of that.

**H4 sent me to the release notes, and they do not say quite what the answer
assumes — they say something better, plus they expose an error of ours.**

To write "the Team is addressing it with a new approach where there is SWIR data"
I needed a source. It is not in `kevin_comments.md` — he says only "Cluster 5 ...
is obviously not good. We need to see how that passed QC." So I read
`WATERHYPERNET_ReleaseNotes_2-0.pdf` directly. It distinguishes **two different
SWIR things** that our document has merged into one:

1. **The QC test** (Release-2 changes): *"Use of a new version of the Similarity
   Spectrum quality control (QC) test based on SWIR bands [6]. This replaces for
   all HYPSTAR sites the previous NIR-based SimSpec QC, which could not be
   applied adequately at the most turbid sites (LPAR, MAFR)..."*
2. **The reflectance correction** (Future Evolution): *"The 'reflectance'
   product, which is currently based on a **NIR-SimSpec correction for both
   PANTHYR and HYPSTAR** systems, **may be upgraded with a SWIR-SimSpec
   correction for HYPSTAR** data."*

**So our line 148 is wrong.** The §5 product table reads *"`reflectance` |
corrected with the Similarity Spectrum test (SWIR-based for HYPSTAR in Release 2,
NIR for PANTHYR)"*. Per the notes the **correction is NIR-based for both
systems**; only the **QC test** is SWIR-based for HYPSTAR. This matters more than
a typo: it is the table that defines the two products a user must choose between,
and a reader who believes HYPSTAR is already SWIR-corrected will misread §9.2 —
whose whole diagnosis is that an **NIR**-referenced offset over-subtracts in dark
water, which is coherent only because the correction really is NIR-based.

**R2-Q1 — Confirm the fix, and tell me how hard to lean on the SWIR sentence.**
I will correct line 148 regardless unless you object. For §9.2 I can write the
note two ways:

- *Sourced (what I can defend today):* "The release notes list a planned upgrade
  of the `reflectance` product to a SWIR-SimSpec correction for HYPSTAR, which
  would bear directly on this behaviour." — accurate, but the notes present that
  as a general improvement, **not** as a fix for dark-site over-subtraction.
- *Stronger (needs your source):* "the team is addressing this via a SWIR-based
  approach" — attributing intent to them. If that came from Kevin directly rather
  than from the notes, say so and I will write it that way and attribute it.

Which? I do not want to put a motive in the network's mouth that their own
documentation does not state.

>A. Let's go with the stronger version, which I learned yesterday at the Team meeting.

**R2-Q2 — "where there is SWIR data": which sites?** Your phrasing implies a
subset, but the notes say the SWIR QC replaced the NIR QC *"for all HYPSTAR
sites"*. Meanwhile our own measurement (§7) is that **no distributed HYPSTAR grid
extends past 1100 nm**, though `epsilon_SWIR` variables are present — so SWIR
radiance clearly existed at processing time but is not in the released files. Do
you mean (a) all HYPSTAR sites, as the notes say; (b) a subset with SWIR-capable
instruments; or (c) you are not sure either? I would rather write (a) with the
"not distributed in these files" caveat than invent a site scope.

>A. (a)

**R2-Q3 — "Access" is in the HOWTO's remit but not in the material.** H1 says the
guide should cover "access, load, and view". Load and view I have in hand.
*Access* exists today as a single clause in §15: distributed via the
password-protected `https://ftp.waterhypernet.org/`. Do you want a real access
page — how to request credentials, who to contact, what the data policy commits
the user to — and if so, can you give me the procedure? I cannot source it from
the archive or the notes.

>A. Yes, let's include that. The procedure is to fill the form out at this webpage: `https://waterhypernet.org/contact/`.  One is then provided FTP access.

**R2-Q4 — The loader bug: fix it in the next prompt, or separately?** The
`quality_flag` masking bug from round 1 still stands (PANTHYR reads `nan` instead
of `0`, so filtering on it drops 21% of the archive). Prompt 3 covers the
restructure and the docstring sweep, and this is neither. I would rather fix it
*before* the HOWTO tells people to call `load_spectrum` — shall I fold it into
prompt 3 with a regression test, or keep it out?

>A. Yes, fix it now.

**R2-Q5 — `make_pdf.py` and the two built PDFs.** "We no longer need the PDF."
Does that mean delete `docs/make_pdf.py` and the two generated files, or leave
the script in place unused? I will not delete anything without you saying so.
Note one knock-on: `respond_to_kevin.md` currently says *"A PDF of the current
version accompanies this note"*, which goes stale either way — flagging it, not
acting on it, since H4 says leave Kevin for now.

>A. You can leave the script and files.  Add a note saying it is deprecated.

**Assumptions I will proceed on unless you say otherwise** — no answer needed:
the H2(a) byline/framing goes on the **landing page only**, not repeated on every
page; H1(b) means I lay the RTD tree out so a package-API section can slot in
later, not that I generate autodoc for two modules now; and I keep the nine-page
split as proposed, though your "a few basics" could argue for merging pages 4 and
5 (**Reading many files** + **Traps**) into one.

### HOWTO — round 3

Short round. Four of the five answers are actionable as given and I have started
nothing that depends on the two below.

**I checked the access URL before agreeing to publish it**, since it will sit in
a public guide. `https://waterhypernet.org/contact/` resolves and carries a real
form: **name, e-mail, affiliation** (all required), **country** (optional), a
required **purpose** dropdown whose options are *"Request for Data access"* and
*"General questions/comments"*, a message box and a privacy-policy checkbox. So
your procedure is right and I can write the access page concretely rather than
vaguely. One nuance worth stating in the text: the page itself says nothing about
FTP or credentials — it only takes the request — so the guide should say that
credentials for `https://ftp.waterhypernet.org/` follow by reply, which is your
knowledge rather than something the site documents.

**R3-Q1 — How should the SWIR sentence be attributed?** You have a source for the
stronger wording (yesterday's Team meeting), so I will use it. The open question
is what the *published* sentence says, because an unminuted meeting is not
something a reader can check, and this sits in the one section the document
promises Kevin right of reply on. Three forms:

- **(a) Bare assertion** — "The team is addressing this with a SWIR-based
  approach." Cleanest to read; a reader has no way to verify it.
- **(b) Personal communication** *(my recommendation)* — "...is being addressed
  with a SWIR-based approach (WATERHYPERNET team, pers. comm., September 2026)",
  optionally naming Kevin. Standard scientific practice for exactly this, and it
  signals to Kevin on review that you are reporting the team's own statement
  rather than speculating.
- **(c) Notes-anchored plus confirmation** — lead with the published Future
  Evolution bullet (which is citable) and add that the team has since confirmed
  it is the intended route.

I recommend (b). Do you want Kevin named in it, or the team collectively?

>A. Ok, (b) is fine

**R3-Q2 — Does prompt 4 produce one file or nine?** This is the one that changes
what I hand you. H5 settled on multi-page, but prompt 4 says "restructure
`docs/WATERHYPERNET.md`" (one file) while prompt 5 covers the RTD "navigation",
which is arguably where page boundaries belong. Options:

- **(a) Restructure in place, split later** *(my recommendation)* — prompt 4
  yields a single reordered `WATERHYPERNET.md` with the new section order, the
  corrected line 148, the §14 rework, the access section and the SWIR note;
  prompt 5 splits it into pages and writes the toctree. The reason I prefer it:
  reordering and splitting in one step gives you nine new files and one deletion
  to review, which is effectively unreviewable — whereas a single-file reorder is
  a diff you can actually read and argue with.
- **(b) Split now** — prompt 4 emits `docs/howto/*.md` directly. Fewer total
  steps, but you review the prose and the page boundaries at the same time, and
  `WATERHYPERNET.md` disappears in the same commit that rewrites it.

>A. (a)

**Decisions I am taking without asking** (say if any is wrong):

- **The loader fix will be narrow — `quality_flag` only.** Hazard 5 lists three
  PANTHYR variables with `_FillValue = 0`: the flag, the angles and `bandwidth`.
  Un-masking all of them would be a regression, not a fix: §11 records that
  PANTHYR's `solar_azimuth_angle` is fill in every file sampled, so a blanket
  `set_auto_mask(False)` would turn genuinely-missing angles into confident
  `0.0`. Only `quality_flag` has the property that 0 is a real, meaningful value
  (no bits set = passed). Narrow fix, regression test covering both a PANTHYR and
  a HYPSTAR file.
- **Deprecation note for `make_pdf.py`** goes in its module docstring, pointing at
  the RTD build as the replacement. Script and both PDFs stay on disk.
- **Access becomes a section of page 1 ("About the archive")**, not its own page,
  in keeping with "a few basics".

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

### 2026-09-19 — HOWTO prompt 2 (answers reviewed; round 2 posed)

Opus 5, no subagents. **No file edited but this one** — the prompt asks whether
more information is needed, not for the restructure.

- Read all six answers. **Five need nothing further** (H1(b) with package docs
  deferred, H2(a), H3 keep-and-generalise §14 with the two IOPtics links, H5
  multi-page and no PDF, H6 drop manual numbers); I am ready to build those.
- **H4 could not be written as given, so I went looking for a source.** The
  answer says to note that the over-subtraction "is a known problem and the Team
  is addressing it by using a new approach where there is SWIR data". That claim
  is **not** in `kevin_comments.md` — his only remark on it is "Cluster 5 ... is
  obviously not good. We need to see how that passed QC." Rather than write an
  unsourced statement about a third party's plans into a document published under
  JXP's name, I extracted and read `WATERHYPERNET_ReleaseNotes_2-0.pdf`
  (`pdftotext -layout`; the archive's only non-NetCDF file).
- **The notes support the claim, but distinguish two SWIR things our document has
  merged.** Release-2 changes: *"a new version of the Similarity Spectrum
  **quality control (QC) test** based on SWIR bands ... replaces for all HYPSTAR
  sites the previous NIR-based SimSpec QC"*. Future Evolution: *"The
  'reflectance' product, which is currently based on a **NIR-SimSpec correction
  for both PANTHYR and HYPSTAR** systems, may be upgraded with a SWIR-SimSpec
  correction for HYPSTAR data."* So the **QC test** is SWIR-based today; the
  **correction** is NIR-based for both and the SWIR version is planned.
- **That makes line 148 of our document factually wrong.** The §5 product table
  says `reflectance` is "corrected with the Similarity Spectrum test (SWIR-based
  for HYPSTAR in Release 2, NIR for PANTHYR)". It conflates the QC test with the
  correction. Not cosmetic: it is the table defining the two products a user
  chooses between, and §9.2's whole diagnosis — an NIR-referenced offset
  over-subtracting in dark water — only coheres because the correction really is
  NIR-based. A reader who believed HYPSTAR was already SWIR-corrected would read
  §9.2 as incoherent. Queued for the restructure; quotes verified against the
  extracted text line by line before writing them down.
- **Five questions posed (R2-Q1..Q5).** The one that matters is Q1: I can write
  the SWIR note *sourced* ("the notes list a planned upgrade ... which would bear
  on this") or *stronger* ("the team is addressing this"), and the notes present
  the upgrade as a general improvement rather than as a fix for dark-site
  over-subtraction — so the stronger form needs a source I do not have. Also
  asked: which sites "where there is SWIR data" means, given the notes say *all*
  HYPSTAR sites while no distributed grid passes 1100 nm; whether the guide needs
  a real "access" page (I cannot source the credential procedure); whether to fold
  the `quality_flag` loader fix into the restructure; and whether "no longer need
  the PDF" means deleting `make_pdf.py` or just retiring it.
- Stated three assumptions I will proceed on without an answer: byline on the
  landing page only, RTD laid out so a package-API section can slot in later
  rather than autodoc now, and the nine-page split retained.

### 2026-09-19 — HOWTO prompt 3 (round-2 answers reviewed; short round 3 posed)

Opus 5, no subagents. **No file edited but this one** — prompt 3, like prompt 2,
asks whether more is needed rather than for the restructure.

- Read all five round-2 answers. Four are actionable as given: the stronger SWIR
  wording (JXP has it from a Team meeting), scope (a) all HYPSTAR sites, fix the
  loader bug now, and keep `make_pdf.py` with a deprecation note.
- **Verified the access URL before agreeing to publish it.** `waterhypernet.org/
  contact/` resolves and carries a real form — name, e-mail and affiliation
  required, country optional, a required purpose dropdown offering *"Request for
  Data access"* / *"General questions/comments"*, a message box and a privacy
  checkbox. JXP's procedure is correct and the access section can now be written
  concretely. Noted one nuance for the text: the page itself says nothing about
  FTP or credentials, so "credentials follow by reply" is JXP's knowledge rather
  than documented on the site, and the guide should not imply otherwise.
- **Two questions posed (R3-Q1, R3-Q2)**, both of which change the next
  deliverable rather than merely refining it:
  - *Attribution.* An unminuted Team meeting is not checkable by a reader, and
    the sentence lands in the section the document promises Kevin right of reply
    on. Offered bare assertion / pers.-comm. citation / notes-anchored-plus-
    confirmation; recommended the pers.-comm. form and asked whether to name
    Kevin or the team.
  - *One file or nine.* H5 chose multi-page, but prompt 4 says restructure
    `WATERHYPERNET.md` while prompt 5 owns RTD navigation. Recommended
    restructuring in place and splitting in prompt 5, on the grounds that
    reordering and splitting together produces nine new files plus a deletion —
    an unreviewable change — whereas a single-file reorder is a readable diff.
- **Three decisions taken without asking**, each recorded with its reason: the
  loader fix stays narrow to `quality_flag`, because §11 shows PANTHYR's
  `solar_azimuth_angle` is genuinely fill everywhere and a blanket
  `set_auto_mask(False)` would turn missing angles into confident zeros — only
  the flag has 0 as a real value; the `make_pdf.py` deprecation note goes in the
  module docstring; and access becomes a section of page 1 rather than a page of
  its own.

### 2026-09-19 — HOWTO prompt 4 (prompts revised against the settled discussion)

Opus 5, no subagents. **Only the Prompts section of this file changed** — no code,
no document, no config.

- Read the round-3 answers: pers.-comm. attribution (b), restructure-in-place
  (a), and no objection to the three decisions I had taken. With rounds 1-3 all
  closed, the old prompts 5-7 were written before most of those answers existed
  and no longer described the work.
- **Added a settled-decisions digest** ahead of the prompts, so each can be run
  standalone without re-reading ~200 lines of Q&A: scope, voice, the nine-page
  order, numbering, the IOP section with both IOPtics links, the SWIR note with
  its attribution and site scope, the access procedure, and the PDF decision —
  each tagged with the question it came from.
- **Listed the five defects as D1-D5** and routed them: D5 (the loader bug) to
  prompt 5, D1-D4 to prompt 6. Relabelled them from plain numbers because the
  defect list and the prompt list would otherwise both count from 5.
- **Split the old prompt 5 in two, which is the one structural change I made
  rather than was told to make.** It bundled the loader fix, the docstring sweep
  and the whole document restructure. The loader fix has to land *before* the
  guide tells readers to call `load_spectrum`, and mixing a code change plus a
  regression test into the same step as a 500-line prose rewrite makes both
  harder to review. So: prompt 5 is code only (loader fix, docstring sweep,
  deprecation note, `pytest -q`), prompt 6 is the single-file restructure.
- **Prompt 7 (RTD) absorbed the page split**, per R3-Q2(a), and lost the
  `make_pdf.py` clause, which R2-Q5 settled. The Markdown-vs-RST decision stays
  open as the prompt's one check-in, with MyST recommended since every page is
  already Markdown. Added: carry the figures and both CSVs into the build, and
  leave a slot for the deferred package-API section.
- **Prompt 8 (Kevin) rewritten as explicitly deferred** — "send nothing and
  regenerate nothing" — since H4 said to leave him for now, and its old
  "regenerate both PDFs" instruction directly contradicted the H5/R2-Q5 decision
  to retire the PDF. It now checks the letter against the restructured document
  and reports, including the stale "A PDF of the current version accompanies this
  note" line.

### 2026-09-19 — HOWTO prompt 5 (loader bug fixed, docstring sweep, make_pdf deprecated)

Opus 5, no subagents. Code only; `docs/WATERHYPERNET.md` untouched. `pytest -q`
→ **21 passed** (19 before, +2 from the new parametrised regression test).

- **Read the variable's metadata before writing the fix**, which turned a guess
  into a certainty. PANTHYR's `quality_flag` is `uint8` with `_FillValue = 0` on
  a bitmask (`flag_masks = 1, 2, 4 ...`), so 0 means *no bits set* = **passed**;
  the release has declared its meaningful value to be its fill value. HYPSTAR's
  is `uint32` with **no `_FillValue` at all**. That is what makes the narrow fix
  safe rather than merely plausible: bypassing the mask cannot change HYPSTAR,
  because HYPSTAR has no mask to bypass.
- **D5 fixed** with a dedicated `_quality_flag(ds)` in `hypernet/whn_explore.py`
  rather than a flag on `_scalar`, so the special case is confined and
  self-documenting. Its docstring states the trade-off honestly: a genuinely
  absent flag is now indistinguishable from a pass, which is unavoidable when the
  fill value *is* the meaningful value, and a pass is what 0 means here.
- **Kept narrow, as agreed.** PANTHYR's angles carry the same `_FillValue = 0`
  but are *genuinely* fill — the document records `solar_azimuth_angle` as empty
  in every file sampled — so they keep the mask and stay NaN. Only the flag is
  read raw.
- **Regression test, and I checked it actually regresses.** Parametrised over
  `O1BE_P` and `THFR_H`, asserting the loaded flag is finite and equals the raw
  unmasked value. Reverting the one-line fix makes the PANTHYR case fail
  (`isfinite(nan)`) while the HYPSTAR case still passes — which is the point:
  it proves the test catches the bug *and* that the fix is a no-op for HYPSTAR.
  Fix restored and verified afterwards.
- **Docstring sweep done, and it was wider than the three items in Context.**
  The listed three: `whn_explore.py:8` (`waterhypernet_prompts.md` →
  `explore_prompts.md`), `whn_explore.py:28` (the `build_v1.py` driver-pattern
  reference dropped), `test_whn_explore.py:3-4` (the `ioptics/tests/conftest.py`
  reference dropped). Three more found while in there, all stale since the Move:
  the run commands still read `python whn_explore.py 1` instead of
  `python -m hypernet.whn_explore 1`; the module cross-reference still pointed at
  `:mod:`whn_figures`` as if it were a package module rather than
  `docs/whn_figures.py`; and the product bullet described the correction as
  "NIR/SWIR", **the same conflation as document defect D1** — corrected here too,
  with an explicit note that the SWIR thing is the QC test. A grep for
  `waterhypernet_prompts|build_v1|ioptics` across the package now returns nothing.
- **`docs/make_pdf.py` marked deprecated** in its module docstring, pointing at
  the RTD build as the replacement and repeating the macOS-only caveat. Script
  and both PDFs left on disk, as agreed.
- **One thing found and deliberately not changed.** `whn_simspec_check.py:75`
  reads `quality_flag` with no mask control, which looks like the same bug but
  is not: `np.asarray` discards the mask, so it already yields 0 — verified
  against a real PANTHYR file, identical to `_quality_flag`. It is also
  unreachable from PANTHYR data (the scan is HYPSTAR-only, by default and under
  `--sites all`). So it is correct, but correct *by accident*, relying on
  `np.asarray`'s behaviour rather than on intent. Left alone as out of the
  prompt's scope and behaviourally a no-op; flagged for JXP, since it is a
  one-line change to route it through the shared helper if he wants the codebase
  consistent.


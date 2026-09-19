# WATERHYPERNET -- Exploring this dataset

## Goal

Examine at a high-level the WATERHYPERNET dataset, Release 2

## Conventions

- `ocean14`; JXP runs git; after each step run `pytest -q` where relevant, Q&A, Log.
- Run via the env interpreter directly
  (`/home/xavier/miniforge3/envs/ocean14/bin/python …`); `conda activate` fails
  non-interactively. This is a **Tier-2** operation — it needs the `$OS_COLOR`
  data tree (L23) mounted; do **not** unset `$OS_COLOR` here.
- `report`/`plotting` need matplotlib (`Agg` is fine headless) + Bokeh + Sphinx.

## Context

- **Config** (`run_v1.yaml`): `sweep_id: expb_giop_L23_v1`, `datasets: [L23]`,
  `noise_model: pace`, `seed: 1234`. **`expb_pow` = standard 5-parameter BING**
  (ExpBricaud a_nw + Pow bb_nw; `Adg, Sdg, Aph, Bnw, beta`, k=5) — **MCMC-fit on
  the `mcmc_subset` (1d)**, and χ²-fit over all records like every algorithm;
  **`giop`** = k=3 contrast (GIOP a_nw + Lee bb_nw; `Aexp, Aph, Bnw`) — **LM/χ²
  only (1b)**. `fit_method: chisq` (sweep default → giop); `mcmc_subset: 200`
  (D4; the expb_pow spectra also MCMC-fit).
- **Driver** (`build_v1.py`, `main(flg)`): `1` → `run.run_sweep(cfg)`; `2` →
  `metrics.compute(cfg.sweep_id)`; `3` → `report.standard.build(sweep_id,
  kind='cross_algorithm')` + `leaderboard.update()` +
  `rst.write_leaderboard_landing(...)`.
- **Artifact split (design Q22):** the `.rst` page + lightweight display assets
  (figure PNGs, `interactive_scatter.html`, CSV tables) are copied into
  `docs/source/reports/<sweep_id>/` and committed; parquet/chains stay in `runs/`.
- **Leaderboard** lands at `$OS_COLOR/IOPtics/leaderboard.parquet` (runs-root
  sibling) and is not committed.

### ⚠ Read before running — the run knobs

- **Error policy / fail-fast.** By default `run.run_sweep` is `strict=True` — a
  **single** failed fit aborts the whole sweep (the run-error-policy TODO). For
  the first full sweep use `strict=False` (D3): failures become
  `status='fit_failed'` rows + reduced coverage.
- **Cost (post-1d): χ² fast, MCMC is the long pole.** χ² over all 3320 × 2 is
  minutes (pooled at `n_cores`). But `expb_pow` (BING) is MCMC-fit on the
  **200-spectrum subset**, and those MCMC fits are **serial** (~40 000 steps
  each). The 5-fit smoke took ~5 min → the full **200 MCMC fits project to
  ≈3 hours**. Plan to run the full sweep in the **background** (step 3).
- **Getting the knobs in.** `build_v1.main` now accepts
  `n_cores`/`strict`/`obs_ids` (step-1 tweak) and threads them into
  `run.run_sweep`.

### WATERHYPERNET

The data files are located in `$OS_COLOR/WATERHYPERNET` under RELEASE_2

## Prompts

### Explore

1. Read the files in Context.  Then explore the data tree for WATERHYPERNET.
We are going to eventually write a context file in named 
`context/WHN/WATERHYPERNET.md`.
First, ask me a series of questions in the Q&A/Explore section below.  
Use Fable if you can.  Log your work.

2. I have answered your first round of Q&A.  Please see my answers and
ask another round if needed.  
Use Fable if you can.  Log your work.

3. I have answered your 2nd round of Q&A.  Please see my answers and
proceed to finish your exploration, generate the Python scripts and figures, and
write the `context/WHN/WATERHYPERNET.md` file.
Use Fable if you can.  Log your work.

4. I shared your exploration file with Kevin Ruddick and he provided some comments.  I put these in a file named `context/WHN/kevin_comments.md`. Please see his comments and discuss them with me in the Q&A/Explore section below.  We will modify the file after our discussion. Use Fable if you can.  Log your work.

5. I have answered your queestions and comments that came from reading Kevin's comments.  Please see my answers and proceed to finish your exploration, generate the Python scripts and figures, and update the `context/WHN/WATERHYPERNET.md` file.  Also produce the `context/WHN/respond_to_kevin.md` file. Use Fable if you can.  Log your work.

### Move

1. After running through all of the prompts in Explore from within the `IOPtics` repo, I have decided to generate a stand-alone Repo for this effort.  This file is a copy of the `waterhypernet_prompts.md` file.  Please copy over all the files from the `IOPtics` repo related to WATERHYPERNET to this new repo.  I will then remove the branch on the `IOPtics` repo.  Use Fable if you can.  Log your work.  If you have any questions, ask me in Q&A/Move section below.

2. Ok, thanks for that effort.  I have answered your questions below.  After reading them, generate a new prompt doc named `howto_prompts.md` in `claude_prompts` and I will carry on with the HOWTO there.  Use Opus 5.  Log your work.

## Q&A

### Move

**What I did before asking.** All 17 files under `context/WHN/` are copied
byte-for-byte into `hypernet/context/WHN/` (verified with `shasum`: 19/19 files
identical, counting the two gitignored PDFs). `pytest -q` at the repo root now
collects and passes **19 tests** (18 WHN + the existing import test) against the
live `$OS_COLOR/WATERHYPERNET/RELEASE_2`. The scripts import nothing from
`ioptics`, `ocpy` or `bing`, so they run unchanged in the new repo.

**Scope settled from git topology, not just grep.** Only **three** commits are
unique to the `waterhypernet` branch (`2843c25 explore`, `c1e1cd5 report`,
`ebee762 ok`); every older commit is also on `develop` and five other branches.
Those three touch exactly `context/WHN/*` (17 files),
`claude_prompts/waterhypernet_prompts.md`, and two `.gitignore` lines -- nothing
else. **So deleting the branch loses no non-WHN IOPtics work.** One caveat: that
`.gitignore` commit also added `papers/*.pdf`, which is IOPtics-specific and will
go with the branch; re-add it on `develop` if you want it.

**M1 -- Layout: keep `context/WHN/`, or restructure?** I did a faithful copy, so
the tree still reads `context/WHN/` even though this whole repo is now WHN. Three
options: (a) leave it -- zero risk, relative `figs/` links and the flat
`from whn_explore import ...` imports keep working; (b) promote the code into the
package -- `hypernet/whn_explore.py`, `whn_figures.py`, `whn_simspec_check.py`
with the tests under `hypernet/tests/`, documents and figures staying in
`context/WHN/`; (c) flatten to top-level `context/` + `scripts/`. I recommend (a)
for this prompt and (b) as its own follow-up, since (b) means rewriting the
cross-imports and the `fig_dir()`/`__file__` anchoring. Which do you want?

>A. We should treat the document we have created as a draft User HOWTO doc.  I have just created a `docs/` sub-folder.  Move doc-related code into that sub-folder.  Items that might be useful generically show go into `hypernet\`.  I am going to pick up the HOWTO effort in a separate prompt doc and expose via `readthedocs`.

**M2 -- Where should the intermediates live?** `out_root()` still writes to
`$OS_COLOR/IOPtics/whn_explore`, which currently holds 18 MB of cache
(`index.parquet`, `pool.parquet`, `pool_spectra.npz`, `sample.parquet`,
`simspec_check.parquet`). Since IOPtics no longer owns this work, I'd rename it to
`$OS_COLOR/hypernet/whn_explore` and `mv` the existing directory across, so
nothing has to be regenerated. Say the word and I'll do both (it is a one-line
change in `whn_explore.py` plus the move), or tell me to leave the path alone.

>A. Yes, rename as you suggest

**M3 -- Two outward-facing sentences now point at the wrong repo.** These are the
only stale references in prose, and both are in documents meant for Kevin:

- `WATERHYPERNET.md:477` -- "All code is in the public [IOPtics repository](...)
  under `context/WHN/`" (the "How this was produced" section).
- `respond_to_kevin.md:160` -- promises him a link to
  `github.com/ocean-colour/IOPtics` under `context/WHN/` "once it is on `main`".

`hypernet` is a separate public repo, so the honest replacement is
`github.com/ocean-colour/hypernet`. **Has the response to Kevin already gone out?**
If yes, the promised path needs a correction in your next note to him rather than a
silent edit here. I left both untouched pending your answer. (The §14 "Implications
for IOP retrieval work" reference to IOPtics at line 408 is still correct -- that
one is about the IOP package as a real sibling project, so I would keep it.)

>A. I don't think I've sent Kevin that response, but do update it to point at the right place(s).

**M4 -- The PDFs.** `WATERHYPERNET.pdf` and `respond_to_kevin.pdf` came across but
are gitignored, so they will not be committed. They also go stale the moment M3 is
actioned. Regenerate them with `make_pdf.py` after the wording is settled, or do
you want the PDFs committed in this repo so there is a citable frozen artifact?

>A. Ok, we will regnerate.  We won't add the PDF to this Repo.

**M5 -- Four small edits I made outside `context/WHN/`; tell me if any is
unwanted.** (i) `.gitignore` -- added `context/WHN/*.pdf`, carried over from
IOPtics. (ii) `requirements.txt` + `setup.py` -- added **cartopy**
(`whn_figures.py:171`, the site map) and **tabulate** (`whn_figures.py:123`,
`df.to_markdown`). Both were declared in *neither* repo; both happen to be in
`ocean14`, which is why nothing ever failed. (iii) `pytest.ini` -- `testpaths`
widened to `hypernet/tests context/WHN`, because otherwise the 18 WHN tests sat
in the repo and never ran; (iii) becomes moot under M1(b). (iv) `CLAUDE.md` --
the layout section described a repo without `context/WHN/` and stated
`testpaths` pointed only at `hypernet/tests`; I documented the new directory,
the flat-import constraint, and the Chrome dependency below.

>A. These mods are fine

**M6 -- Internal references inside the scripts, less urgent than M3.** These
mislead a reader but break nothing, and most are docstrings:
`whn_explore.py:8` cites `claude_prompts/waterhypernet_prompts.md` (here it is
`explore_prompts.md`); `whn_explore.py:28` says the stages mirror "the
`build_v1.py` driver pattern" and `:34` "the project's artifact split", both
IOPtics conventions with no counterpart in this repo;
`test_whn_explore.py:3-4` points at `ioptics/tests/conftest.py` as the testing
convention. Want these swept in one pass, together with whatever M1/M2/M3
decide? I'd rather do it once than touch the files repeatedly.

>A. I'll save the sweep for the HOWTO effort.

**Not a question, just a prerequisite worth recording:** `make_pdf.py:26`
hard-codes `/Applications/Google Chrome.app/.../Google Chrome` and shells out to
it headless. PDF regeneration is therefore macOS + Chrome only; it raises a clean
`FileNotFoundError` elsewhere. Noted in `CLAUDE.md`; say if you want it in
`README.md` too.

### Explore

**What I found first (so the questions make sense).** 56,669 NetCDF files, 9.8 GB,
11 sites, `SITE_X/YYYY/MM/DD/*.nc`, one above-water hyperspectral measurement per
file. `_H` = HYPSTAR (7 sites), `_P` = PANTHYR (4 sites); VEIT has both. **No IOPs,
no Chl, no in-water measurements of any kind** — the release is reflectance +
geometry + QC only. Both systems ship `reflectance` (NIR/SWIR-SimSpec-corrected)
and `reflectance_nosc` (uncorrected); I verified `reflectance_nosc == pi*Lw/Ed`
exactly, so these are **rho_w, not Rrs** (Rrs = rho_w/pi).

**Q1 — What is this dataset *for* in IOPtics?** WATERHYPERNET carries no truth
IOPs and no Chl, so none of the existing truth-scored metrics apply. Is the intent
(a) Rrs-space closure only — fit, score the Rrs residual, and compare algorithm
*self-consistency* / retrieved-IOP spread with no truth; (b) a match-up target to
be joined later against an external IOP/Chl source (PANGAEA, satellite, PI data);
or (c) purely a spectral-shape exploration for now, no sweep?
>A. I believe there are IOPs and other in-situ measurements.  I will find out, but assume so for now.  Add to the report this glaring absence, however.

**Q2 — Which reflectance product?** The release notes call `reflectance_nosc`
"definitely recommended" at the turbid sites (LPAR, MAFR, O1BE) because the
SimSpec-corrected `reflectance` is "known to be poor" there. Do we (a) always use
`reflectance_nosc`, (b) always use `reflectance`, or (c) use `reflectance` except
at the three turbid sites? And do you want the stored quantity to be Rrs
(= rho_w/pi) or rho_w?
>A. Let's use (c) for now.  I'll ask the experts for further guidance too.

**Q3 — Wavelength trim and grid.** README says <400 nm, >900 nm and ~762 nm
(O2-A) are unreliable and not recommended. Trim to 400-900 nm and notch 762?
Separately, the **HYPSTAR grid changes with instrument serial number** — GAIT,
LPAR, MAFR and VEIT_H each contain two different grids (1536-1541 points,
~0.49 nm), while PANTHYR is a fixed 237-point 355-945 @ 2.5 nm grid. Keep each
record on its native grid (the adapter contract allows it) or resample everything
to one common grid — and if so, which (PANTHYR 2.5 nm? a 5 nm hyperspectral grid?
PACE bands?)?
>A. Right now we are just exploring, so no need to trim or resample yet. 

**Q4 — Scope.** All 11 sites / 56,669 spectra, or a subset? Note the per-site
imbalance: VEIT_H alone is 17,668 files (31%) while WRUK is 889. If we subset, by
what — per-site cap, one spectrum per day, a date window, specific sites?
>A. How about we aim for a ~100 repreesentative spectra per site?

**Q5 — HYPSTAR and PANTHYR: one dataset or two?** They have genuinely different
schemas (dims, dtypes, variable names: `std_reflectance` vs `reflectance_std`,
`irradiance` vs `downwelling_irradiance_mean`) and very different spectral
resolution. Register as one `WATERHYPERNET` adapter that normalises both, or two
(`WHN_HYPSTAR` / `WHN_PANTHYR`) so results are never silently pooled? VEIT is
co-located for both, which would make a nice cross-system check.
>A. I like keeping them separate with `WHN_HYPSTAR` and `WHN_PANTHYR` adapters.

**Q6 — Uncertainty / fit weights.** The release reports **no** uncertainty budget
(`unc_comps = []`, and the README calls uncertainties "not yet mature"). The
`std_*` variables are scan-to-scan radiance variability only — and in the HYPSTAR
files I sampled `std_reflectance` is entirely NaN. Options: (a) use `std_*` as
`Rrs_err` where finite and fall back to the flat 5% `insitu` model otherwise
(what PANGAEA does today), (b) always flat 5%, (c) a fixed floor you'd prefer.
>A. Let's use (a) but with a fall back of 10%

**Q7 — QC and negative reflectance.** `quality_flag` was 0 in every one of the
3,300 files sampled (HYPSTAR has a 30-bit flag definition; PANTHYR's is
placeholder `flag1..flag8` and, with `_FillValue=0`, xarray masks it). So the
release looks pre-filtered. But negative reflectance is common and *deliberately*
retained — 100% of sampled BEFR/GAIT/THFR/VEIT_H/WRUK files contain some negative
band (min -0.039 at LPAR). Do we (a) pass them through untouched, (b) drop
records with negatives in the fit window, (c) drop records where the blue/green
peak itself is negative? This materially changes how many spectra are fittable.
>A.  Do not drop negative Rrs.  I will ask about QC algorithms

**Q8 — Where does the loader live?** Follow the GLORIA/PANGAEA pattern —
`ocpy.insitu.waterhypernet` doing the I/O, thin `WATERHYPERNETAdapter` in
`ioptics/datasets.py` — or keep everything in IOPtics? (`ocpy/insitu/` currently
has only `gloria.py` and `pangaea.py`; nothing WATERHYPERNET exists anywhere yet.)
>A. Let's put the loader in ocpy.  I have moved onto a `whn` branch

**Q9 — obs_id and an index.** Scanning 56k files per sweep is not viable. I'd
propose building a manifest parquet once (site, system, datetime, rel. azimuth,
sza/vza, path, grid id) and keying `obs_id` off it. Filename stem as the id? Note
WRUK has 345 acquisition-timestamp collisions (two azimuths, same `sequence_id`),
so site+timestamp alone is not unique.
>A. That sounds fine, but we aren't doing a sweep yet.  We are just exploring.

**Q10 — How much do you want written down?** Should `context/WATERHYPERNET.md` be
(a) a data-description doc only, or (b) data description + the adapter/loader
design decisions (the answers above), so it doubles as the implementation spec?
>A. Let's do (a) for now. 

### Explore — round 2

**Two of your answers rest on facts I had wrong or unverified, so I checked first.**

**R2-Q1 — The IOPs are definitively not in this release; where should I look?**
You said to assume in-situ IOPs exist. They are not in RELEASE_2. I read every
variable and global-attribute name in 440 files spanning all 11 sites: the sets
are perfectly constant (29 variables HYPSTAR, 22 PANTHYR) and contain nothing
optical-in-water — no absorption, backscatter, chlorophyll, TSS, turbidity,
salinity (the only name that even pattern-matched was `system_temperature`, an
instrument housekeeping attribute, and it is NaN). The release notes agree:
"The core business of HYPERNETS is to provide uniquely valuable in situ
measurements of hyperspectral water reflectance" (p.4), and the word
absorption/chlorophyll/IOP appears nowhere in the document. So: is the companion
data (a) a separate WATERHYPERNET/HYPERNETS product or portal I should be
pointed at, (b) PI-held data per site (the notes name a PI per site) you will
request, or (c) an entirely different archive we match up to? **I have recorded
this absence prominently, as you asked.** One lead from the notes: AERONET-OC is
co-located at VEIT, CBUS, TBBE and LPAR — though that is also radiometry, not
IOPs.
>A. I will find out, but for now don't worry about it.

**R2-Q2 — Your Q6 answer was based on a wrong fact of mine; please re-decide.**
I reported that HYPSTAR `std_reflectance` was all-NaN. That was wrong. It is
present and fully finite in **all 220 files I sampled across all 11 sites**, both
`std_reflectance` and `std_reflectance_nosc` (PANTHYR likewise). Typical
magnitude, median over 450-650 nm of std/|reflectance|: VEIT_P 0.8%, CBUS 1.1%,
MAFR 1.5%, WRUK 1.6%, VEIT_H 1.9%, BEFR 2.1%, GAIT 2.5%, O1BE/LPAR/TBBE 2.5-2.7%,
THFR 2.9% (the `_nosc` values run slightly higher, 1.0-4.0%). So the measured
sigma is always available and your 10% fallback would essentially never fire.
Given that — and given this std is **scan-to-scan variability only**, excluding
calibration and glint-correction systematics, so it understates true uncertainty
— do you want (a) sigma = std as-is, (b) sigma = max(std, X% of Rrs) with a floor
you choose, or (c) keep the flat 10% regardless and ignore the std?
>A. Go with (a)

**R2-Q3 — You did not answer the rho_w vs Rrs half of Q2.** Both products are
rho_w. I recommend the loader store **Rrs = rho_w/pi** (and sigma = std/pi) so
records are directly comparable with L23/PANGAEA/GLORIA and need no special-
casing downstream. Confirm? Also, on your (c) choice: the notes name exactly
LPAR, MAFR and O1BE as sites where the corrected `reflectance` is poor. Is the
`nosc` switch exactly those three, or do you want it per-system (all PANTHYR on
one product) or per-site-configurable?
>A. I confirm.  Let's go with exactly those 3 for now.

**R2-Q4 — What does "representative" mean for the ~100 per site?** The sites are
not comparable in coverage: BEFR is 2023-06 to 2024-01 only, WRUK is a single
year (2024) over 7 months, CBUS 7 months, GAIT 9 — only 6 of 11 sites see a full
seasonal cycle. So pick the axis: (a) uniform in time (every Nth spectrum,
date-sorted — simple, reproducible, but VEIT_H's 17,668 and WRUK's 889 sample
very differently in density); (b) stratified by year x month; (c) stratified by
**optical diversity** — cluster the normalised spectra and sample across
clusters, which is what actually matters for exercising IOP algorithms; (d)
stratified by solar zenith. My preference is (c) with (a) as a tiebreak. Also:
11 site-directories x 100 = ~1,100 spectra, and VEIT counts twice (VEIT_H and
VEIT_P) — is that what you want? And should the sample be forced to cover both
azimuths where they exist (LPAR 4670x270 vs 1312x090; WRUK 476x270 vs 413x225)
and both instrument serials at the four HYPSTAR sites that changed instrument
(each has its own wavelength grid)?
>A. Yes, go with (c)

**R2-Q5 — What should exploration actually produce?** You said no trim, no
resample, no sweep — just explore. So what do you want to look at: a per-site
summary table (N, date span, sza range, Rrs magnitude/shape stats), per-site
median spectrum with percentile envelopes, a spectral-shape/OWT clustering view,
band time series, or something else? And where do the figures go — IOPtics
`papers/`, ocpy `nb/`, or a scratch location?
>A. Let's do:
   - A per site summary Table
   - A figure showing all the WHN locations on the globe
   - A figure indicating the amount of data at each site
   - Per-site median spectrum with percentile envelopes
   - A spectral-shape/OWT clustering view
   - Band time series

Put the Python scripts and figures in the `context/WHN` folder


**R2-Q6 — Code placement, now that the loader is going in ocpy.** Confirmed ocpy
is on the `whn` branch. Should the *exploration* script/notebook also live in
ocpy (`nb/` or `papers/`), or stay in IOPtics? And do you want a skip-guarded
`ocpy/tests/test_waterhypernet.py` following the `test_pangaea.py` pattern, given
the data are not bundled?
>A. Stay in IOPtics, in `context/WHN`.  Yes, have the skip-guarded test.

### Explore — round 3: Kevin Ruddick's comments

He accepts the summary, praises four things specifically (the HYPSTAR/PANTHYR
harmonisation gaps, the loader hazards, the temporal-coverage plot, and the
document's value alongside their own documentation), and raises four things that
need a decision from you. Nothing in `WATERHYPERNET.md` has been changed yet.

**K1 — Shared y-axis on the per-site spectra. Agreed in principle; the question
is how.** His reason is good: a common axis shows instantly which sites are dark
and which are bright. The catch is the dynamic range — median Rrs(560) runs from
0.0030 sr^-1 (WRUK) to 0.036 (MAFR), a factor of 12, and the peaks span 0.003 to
~0.05. On one **linear** axis the four dark sites (WRUK, THFR, BEFR, and to a
degree VEIT/GAIT/CBUS) collapse to near-flat lines and all spectral shape is
lost — which is most of what the figure currently conveys. Three ways:

  (a) **Shared linear** — literally what he asked for; magnitude ordering is
      instant, shape detail at the dark sites is gone.
  (b) **Shared symlog** — both magnitude and shape stay legible, and it handles
      the negative NIR bands honestly; slightly more effort for the reader.
  (c) **Shared linear as the primary figure, keeping the present autoscaled
      version as a second figure.**

I recommend **(c)**: it gives Kevin exactly the comparison he wants without
losing the per-site detail that makes the figure useful for algorithm work. If
you would rather ship one figure, I would pick (b). Which?

>A. Yes, let's use (c)

**K2 — Cluster 5 "obviously not good ... how did that pass QC?" I have the
answer, and I think it is a finding worth sending back to him.** All 9 spectra:

  - **all HYPSTAR, all at the three darkest sites** — THFR_H 5, BEFR_H 2,
    WRUK_H 2. No PANTHYR, no turbid site.
  - **all have `quality_flag = 0` and all 6/6 scans valid.** They passed QC
    cleanly; nothing in the flag bitmask fired.
  - solar zenith is *not* the explanation: five are high-sun-angle (79-86 deg)
    but four are ordinary (23, 26, 52, 60 deg).
  - **the diagnosis is the Similarity Spectrum correction over-subtracting.** In
    every one of the 9, `reflectance_nosc` is positive and entirely plausible
    (median over 440-600 nm: 0.0003 to 0.0071), while the corrected
    `reflectance` is at or below zero or close to it. The correction removes
    more than the entire signal.

That is the mirror image of the caveat in the release notes: they warn that
`reflectance` is poor at the **turbid** sites (LPAR, MAFR, O1BE), whereas this
says it also fails at the **darkest** sites, where the subtracted offset exceeds
the water signal. Rarity, over the 4,400-spectrum pool: 9 in cluster 5 (0.2%),
and 4 spectra whose visible-band (440-600 nm) median Rrs is actually <= 0
(THFR_H 3, BEFR_H 1).

Three questions: do you want this (i) added to `WATERHYPERNET.md` as a short
subsection, (ii) sent to Kevin as a separate note with the 9 filenames, or
(iii) both? And should I **extend the check across the full archive** for the
three dark sites (THFR + BEFR + WRUK = 8,947 files, roughly 10 minutes) so he
gets a real rate rather than a pool estimate? I would do (iii) plus the full
scan — it is the most useful thing in his reply and it is cheap.

>A. Let's do both.  We will create a file named `context/WHN/respond_to_kevin.md`.  Yes, do the full scan too.

**K3 — He pushes back on the IOP framing, and I think he is substantially
right.** My section is titled *"The glaring absence: there are no IOPs, and no
biogeochemistry at all"*. His points: HYPERNETS is deliberately a single-
parameter network designed as a core for expansion; reflectance is in his view
the most important parameter; a reflectance-only archive supports four real
applications (satellite validation, water-quality monitoring, spectral-
variability research, AERONET-OC cross-comparison); and IOPs are "an
intermediate, explanatory parameter, not an end-user parameter".

"Glaring absence" frames a deliberate design choice as a deficiency, and he has
the standing to correct that. But you asked in round 1 for the absence to be
recorded prominently, and the *fact* matters for IOPtics regardless — so I do
not want to quietly drop it. My proposal is to keep every verified fact and the
IOPtics consequence, and change only the judgement:

  - retitle to something like **"Scope: reflectance only, by design"**;
  - keep the exhaustive verification (440 files, constant variable set) and keep
    the plain statement that no truth-scored IOPtics metric applies;
  - add his four use cases and the design rationale, attributed to him;
  - add that some PIs hold extra measurements from short-duration deployments,
    and that he will ask them and suggest HYPERNETS advertise this.

Confirm the reframing, or tell me to hold the original tone?

>A. Yes, reframe.

**K4 — Dissemination. This one is yours, not mine.** He offers two options, not
exclusive: link from a WATERHYPERNET "exploitation tools" page to the GitHub
section *"if it is stable"*, and/or freeze a PDF with author (XP) and a date into
the WATERHYPERNET documentation. Things to settle before you answer him:

  - **Authorship and credit.** He proposes "(XP)". Your call entirely — including
    whether and how this conversation's tooling is acknowledged.
  - **Stability.** `github.com/ocean-colour/IOPtics` is public, but `context/WHN/`
    is currently **uncommitted on the `waterhypernet` branch** — nothing for him
    to link to yet. A link implies you will keep that path stable; a frozen PDF
    carries no such commitment. If you want both, the PDF is the safer first
    step.
  - **Audience.** The document is written for IOPtics. If it becomes a
    WATERHYPERNET-facing artifact, a few things should change: the
    "Consequence for IOPtics" framing, the `$OS_COLOR` paths, the artifact-split
    note, and the internal references to `claude_prompts/`. I would produce an
    external variant rather than repoint the internal one.
  - He should probably also get right of reply on anything that reads as
    criticism of their processing (K2 above, and the hazards list) before it is
    published under their documentation.

>A. Answers:
   - Use JXP for the Authorship and credit.
   - We will get this report onto `main` soon.  But we should also create a PDF
   - Let us make it a WATERHYPERNET-facing artifact.  There can be a final section discussing the IOPtics implications.
   - Agreed, put comemnts on this in the `context/WHN/respond_to_kevin.md` file.


**K5 — A forward-looking note, no action needed yet.** What he would most like
to see added to the sites is automated flow cytometry, turbidimeters and
chlorophyll fluorimeters — explicitly *not* IOPs. If companion data ever
arrives, it will most likely be **Chl, turbidity and phytoplankton composition**.
That is better news for IOPtics than it sounds: `Chl` and `tss` are already
scalar truth keys in `ioptics/datasets.py`, so that data would map onto existing
machinery, whereas spectral `a_ph`/`a_dg`/`bb_p` truth would not arrive at all.
Worth a sentence in the document? I think yes.

>A. Sounds good.  Indeed, Chl, turbidity and phytoplankton composition would be terrific.  And, maybe WATERHYPERNET will be more valuable for my retrieve-or-bust (RoB) effort. 

## Logs

### 2026-09-15 — Explore prompt 1 (survey, no code written)

- Read the Context files: `ioptics/runs/prototypes/expb_giop/run_v1.yaml` and
  `build_v1.py` (L23-only today; datasets resolved through the
  `ioptics.datasets` registry).
- Reviewed the adapter contract in `ioptics/datasets.py` —
  `obs_ids(**opts)` + `load_obs(obs_id) -> RawObs(wave, Rrs, truth, Rrs_err,
  meta)`, `register_dataset(...)`. Three precedents: `L23Adapter` (synthetic,
  full spectral truth), `PANGAEAAdapter` (per-family native grids, no
  `Rrs_err` -> flat 5%), `GLORIAAdapter` (hyperspectral + measured per-band
  std -> genuine `insitu` weights, scalar-only truth).
- Confirmed no WATERHYPERNET loader exists in `ocpy.insitu` (only `gloria.py`,
  `pangaea.py`) or in IOPtics. `context/` was empty.
- Surveyed `$OS_COLOR/WATERHYPERNET/RELEASE_2` with a Fable agent (read-only):
  file counts/sizes per site, filename grammar, full NetCDF structure for one
  HYPSTAR and one PANTHYR file, QC flag definitions, and the Release-2 notes PDF
  (`0_README/WATERHYPERNET_ReleaseNotes_2-0.pdf`).
- Spot-verified by hand: 56,669 `.nc` / 9.8 GB; BEFR HYPSTAR file = 1539 bands,
  350.14-1099.75 nm, `quality_flag=0`, 75 negative bands; O1BE PANTHYR file =
  237 bands, 355-945 nm, `quality_flag` masked by `_FillValue=0`.
- Established the product definition numerically: `reflectance_nosc` /
  (`water_leaving_radiance` / Ed) = 3.1416 = pi exactly for both systems, while
  the SimSpec-corrected `reflectance` gives 3.20 (HYPSTAR) / 2.96 (PANTHYR).
  So both products are **rho_w**, `water_leaving_radiance` is the *uncorrected*
  Lw, and **Rrs = rho_w / pi**.
- Key blocker for a sweep: the release contains **no truth** (no IOPs, no Chl,
  no TSS) and **no uncertainty budget**. Questions posed above; no code written
  and nothing modified under `$OS_COLOR`.

### 2026-09-15 — Explore prompt 2 (answers reviewed, round 2 posed)

- Read the round-1 answers. Two of them rested on claims that needed checking
  before a second round was worth asking, so I verified both.
- **Correction to my own round-1 report:** I stated HYPSTAR `std_reflectance`
  was "entirely NaN". It is not. Sampling 20 files per site x 11 sites (220
  files), `std_reflectance` and `std_reflectance_nosc` are **fully finite in
  every file**, read both with netCDF4 auto-masking on and off with `_FillValue`
  excluded explicitly. Median std/|reflectance| over 450-650 nm ranges 0.8%
  (VEIT_P) to 2.9% (THFR). The Q6 answer (10% fallback) was made on the wrong
  premise and is re-asked as R2-Q2.
- **Verified the absence of IOPs.** Union of all variable and global-attribute
  names over 440 files spanning all 11 sites: constant at 29 vars / 57 attrs
  (HYPSTAR) and 22 vars / 58 attrs (PANTHYR). No absorption, backscatter,
  chlorophyll, TSS, turbidity or salinity variable exists anywhere; the only
  regex hit was the housekeeping attribute `system_temperature` (NaN).
- Extracted the full text of `0_README/WATERHYPERNET_ReleaseNotes_2-0.pdf`
  (6 pages, 14,605 chars) and searched it: no mention of IOPs, chlorophyll, TSS
  or water sampling anywhere. p.4 states the network's "core business" is
  hyperspectral water reflectance. The only in-situ lead is AERONET-OC, noted as
  co-located at VEIT, CBUS, TBBE and LPAR. Negative reflectance is retained
  deliberately (p.2: "metrologically entirely valid", and filtering biases the
  statistics).
- Filename-only coverage census (no file opens): per-site year/month/hour spans
  and azimuth-token counts. Key result for the ~100-per-site sample — only 6 of
  11 sites span a full seasonal cycle; BEFR (2023-06..2024-01), WRUK (2024 only),
  CBUS (7 months) and GAIT (9 months) do not. Azimuth counts re-derived
  independently and they match round 1.
- Confirmed ocpy is on branch `whn`; reviewed the `ocpy.insitu` house pattern
  (`gloria.py` / `pangaea.py`: `$OS_COLOR/<DATASET>` path resolution, module-level
  loaders, skip-guarded tests).
- Two Fable agents launched for the census work died without reporting after
  ~30 min (over-specified prompts: per-file double reads plus per-band stats over
  ~1,650 files each). Redone directly with narrower scans, as described above.
- No code written; nothing modified under `$OS_COLOR`; `context/` still empty.

### 2026-09-15 — Explore prompt 3 (exploration finished; scripts, figures, context doc)

Built to the round-2 answers: `reflectance` except LPAR_H/MAFR_H/O1BE_P
(`reflectance_nosc`); Rrs = rho_w/pi and sigma = sigma(rho_w)/pi; archive sigma
used as-is, no floor; optical-diversity sampling; no trim, no resample; negatives
kept; code and figures in `context/WHN`; skip-guarded test.

- **`context/WHN/whn_explore.py`** — data layer. Stage 1 indexes all 56,669 files
  from their filenames alone (both grammars parse; zero unparsable). Stage 2
  reads 400 spectra per site spread uniformly through each site's record (4,400
  total, **zero read failures**), clusters the L2-normalised 400-800 nm shapes
  across all sites at once (KMeans k=8, seed 1234), and apportions each site's
  100 spectra across the clusters it occupies by largest remainder, min one per
  occupied cluster, spread in time within a cluster. L2 rather than area
  normalisation because retained negatives can drive an integral to zero.
- **`context/WHN/whn_figures.py`** — `summary_table.{csv,md}` plus the five
  requested figures in `figs/`: site map (world + Europe inset), data volume
  (bar + month x site coverage heatmap), per-site median spectra with 25-75 and
  10-90 envelopes, clustering view (cluster shapes, per-site composition,
  absolute-unit cluster medians, blue-green vs magnitude scatter), and band time
  series at 490/560/665 nm.
- **`context/WHN/test_whn_explore.py`** — 14 tests, Tier 1 logic + Tier 2
  `@needs_whn`. All pass; full suite **282 passed**.
- **`context/WHN/WATERHYPERNET.md`** — the context document (data description
  only, per Q10), leading with the absence of IOPs.
- Intermediates (index/pool/sample parquet, spectra npz) go to
  `$OS_COLOR/IOPtics/whn_explore/`, not the repo, per the artifact split.

Findings added beyond rounds 1-2:

- **Optical range is wide.** Median Rrs(560) spans 0.0030 sr^-1 (WRUK) to 0.036
  (MAFR), a factor of 12; blue-green Rrs(490)/Rrs(560) spans 0.52 (BEFR) to 1.03
  (VEIT_H). Sites occupy distinct clusters rather than one continuum.
- **Grid count corrected.** Counting distinct (n_bands, start) pairs rather than
  band count alone: GAIT, LPAR, MAFR and VEIT_H each carry 2 HYPSTAR grids
  (MAFR's two instruments both give 1538 bands but start 0.18 nm apart). PANTHYR
  sites have used 2-3 sensor serials each but the 237-band grid never changes.
- **quality_flag verified myself** this time: 0 in all 440 files checked (40 per
  site, raw values with masking off). The earlier 3,300-file figure came from the
  agent survey and is not repeated in the doc.
- **A negative-spectrum cluster exists** (n=9 of 4,400): spectra negative across
  most of the range. Retained, not dropped, per the release policy.
- Three bugs found and fixed while building: netCDF4 cannot fill an integer
  variable with NaN (PANTHYR `quality_flag` is uint8 with `_FillValue=0`), which
  killed every PANTHYR read until the cast-to-float-first fix; cartopy silently
  auto-selects a finer Natural Earth scale for a zoomed extent and blocks on the
  download (scales now pinned); and cartopy does not clip annotation text, so
  out-of-extent site labels distorted the saved map.

### 2026-09-18 — Explore prompt 4 (Kevin Ruddick's comments reviewed, discussion posed)

- Read `context/WHN/kevin_comments.md`. He accepts the summary and singles out
  four things as useful (HYPSTAR/PANTHYR harmonisation gaps, the loader-hazards
  list, the temporal-coverage figure, and its value beside their own docs).
  Four items need decisions; they are posed in Q&A round 3. **No change made to
  `WATERHYPERNET.md` yet** -- the prompt says we modify after discussing.
- **Investigated his Cluster-5 question** ("obviously not good... how did that
  pass QC?") rather than just relaying it. All 9 members are HYPSTAR at the three
  darkest sites (THFR_H 5, BEFR_H 2, WRUK_H 2); every one has `quality_flag = 0`
  and 6/6 valid scans; solar zenith does not explain them (5 at 79-86 deg but 4
  at 23-60 deg). The cause is the Similarity Spectrum correction over-
  subtracting: in all 9, `reflectance_nosc` is positive and plausible (median
  440-600 nm, 0.0003-0.0071) while the corrected `reflectance` sits at or below
  zero. This is the mirror image of the release notes' warning, which flags
  `reflectance` as poor at the *turbid* sites -- it also fails at the *darkest*
  ones, where the offset exceeds the water signal.
- Quantified over the 4,400-spectrum pool: cluster 5 is 9 spectra (0.2%); 4
  spectra have a visible-band (440-600 nm) median Rrs <= 0 (THFR_H 3, BEFR_H 1).
  Offered to extend the scan to all 8,947 files at the three dark sites for a
  true rate.
- On the shared-y-axis request: measured the obstacle rather than assuming one.
  Median Rrs(560) spans 0.0030 (WRUK) to 0.036 sr^-1 (MAFR), a factor of 12, so a
  shared *linear* axis flattens the dark sites. Three options put to JXP, with
  shared-linear-plus-autoscaled-companion recommended.
- On his push-back over the IOP framing: recommended keeping every verified fact
  and the IOPtics consequence but retitling away from "glaring absence", since
  that judges a deliberate single-parameter network design. Held for JXP's
  confirmation rather than acting, because JXP asked in round 1 for the absence
  to be prominent.
- Checked one fact relevant to his dissemination offer: the IOPtics GitHub repo
  is public, but `context/WHN/` is uncommitted on the `waterhypernet` branch, so
  there is nothing for him to link to yet. Authorship, stability commitment and
  the audience shift are flagged as JXP's decisions.
- No code or figures changed this round.

### 2026-09-18 — Explore prompt 5 (Kevin's comments actioned; doc reframed, response drafted)

Built to the round-3 answers: figure option (c); cluster-5 finding both in the
document and in a response file, plus the full scan; IOP framing reframed;
WATERHYPERNET-facing artifact authored "J. Xavier Prochaska" with a final
IOPtics section; PDF produced.

- **Full over-subtraction scan** (`whn_simspec_check.py`, new). Scanned all
  **8,947 files** at the three dark sites. Flag = `median(reflectance) <= 0 <
  median(reflectance_nosc)` over 440-600 nm. Result **30 of 8,947 (0.335%)** --
  THFR_H 26/6,260 (0.42%), BEFR_H 4/1,798 (0.22%), WRUK_H 0/889. All 30 have
  `quality_flag = 0` and 6/6 valid scans. Solar zenith spans 21.1-87.6 deg
  (median 54.4; only a third above 70), so low sun is **not** the explanation --
  BEFR's four are all Dec 2023 at sza 85-87, but THFR's 26 span all months and
  all sun angles. The common factor is **water darkness**: flagged spectra have
  uncorrected rho_w(440-600) median 0.00043 against site medians of 0.0106 /
  0.0097 / 0.0072, i.e. ~20x darker than typical. The subtracted offset is
  0.00063, about **145% of the entire signal**. Spread over 2023/2025/2026, so
  persistent. Exported as `simspec_flagged.csv`.
- Noted carefully that this rate is a **lower bound**: it counts only sign
  changes, whereas the shape clustering also caught WRUK spectra distorted
  without the median going negative (which is why WRUK reads 0.00% here but
  contributed 2 of the 9 cluster-5 members). Both statements are in the doc and
  the response.
- **Figures**: `fig_median_spectra_autoscale.png` added; `fig_median_spectra.png`
  is now the shared-axis version, and both are in the document as a pair. The
  shared axis does exactly what Kevin predicted -- LPAR/MAFR/O1BE dominate while
  BEFR/THFR/WRUK (and, less obviously, VEIT/GAIT) are nearly flat.
- **`WATERHYPERNET.md` rewritten as an external artifact**: authored J. Xavier
  Prochaska, dated, numbered sections, reframed section 2 ("Scope: reflectance
  only, by design") carrying Kevin's rationale and his four use cases, new
  section 9.2 on the over-subtraction, and a final section 14 confining the
  IOPtics implications. Internal references (`$OS_COLOR`, `claude_prompts/`, the
  artifact-split note) removed in favour of the public repo path.
- **`respond_to_kevin.md`** written to Kevin: what changed, the cluster-5
  investigation in full with the 10 worst cases tabulated, the reframing and
  thanks for the correction, and the dissemination answer -- PDF first (no
  maintenance commitment), repo link once merged to `main`, plus an explicit
  offer of right of reply on section 9.2 and the hazards list, including pulling
  9.2 from the public version if he would rather raise it internally first.
- **PDF**: none of pandoc/weasyprint/wkhtmltopdf/python-markdown are installed in
  `ocean14`, so `make_pdf.py` converts the Markdown subset used here to styled
  HTML and prints it with the headless Chrome already on the machine.
  `WATERHYPERNET.pdf` (1.8 MB, ~14 pp) and `respond_to_kevin.pdf` verified by
  reading the rendered pages.
- Tests extended to 18 (4 new for the over-subtraction check, including a
  regression guard on a specific flagged spectrum). Full suite **286 passed**.
- Nothing committed; `__pycache__` is already gitignored.

### 2026-09-18 — Move prompt 1 (WHN files copied from IOPtics into hypernet)

- **Defined the WHN set from git topology, not just grep.** `git log --oneline
  main..waterhypernet` looks alarming (~50 commits, 19k lines, all of stages 3-6,
  the reporting subsystem, the docs site), but `git branch --contains` shows only
  **three** commits are unique to the branch: `2843c25 explore`, `c1e1cd5 report`,
  `ebee762 ok`. Everything older is also on `develop`, `first-go`,
  `improve-reporting`, `ls2` and `rt-tests`. `git diff --name-status
  26a9801..waterhypernet` gives the true WHN footprint: the 17 files under
  `context/WHN/`, `claude_prompts/waterhypernet_prompts.md`, and `.gitignore`.
  A `grep -ril` for waterhypernet/hypstar/panthyr/WHN over the whole IOPtics tree
  returned the same set, so the two methods agree. **Deleting the branch is safe**
  (noted in Q&A: it does take the IOPtics-only `papers/*.pdf` gitignore line with
  it).
- **Copied** `context/WHN/` → `hypernet/context/WHN/` with `rsync --exclude
  __pycache__`, then verified with `shasum`: **19/19 files byte-identical**,
  including the two gitignored PDFs. Nothing was modified in IOPtics.
- `claude_prompts/waterhypernet_prompts.md` needed no copying — `diff` confirms
  `explore_prompts.md` is already it, plus the Move section.
- **The scripts are self-contained.** No import of `ioptics`, `ocpy` or `bing`
  anywhere in `context/WHN/*.py`; the only cross-references are flat
  (`from whn_explore import ...`), which survive the move because the directory
  moved whole.
- **Verified it actually runs in the new repo**, not just that files landed:
  `conda run -n ocean14 python -m pytest -q` → **19 passed** (18 WHN + the
  pre-existing import test) against the live
  `$OS_COLOR/WATERHYPERNET/RELEASE_2`. All six `figs/*.png` links in
  `WATERHYPERNET.md` resolve.
- **Three repo-level fixes** (all listed in Q&A M5 for approval): `.gitignore`
  gained `context/WHN/*.pdf`; `requirements.txt` + `setup.py` gained **cartopy**,
  a hard import at `whn_figures.py:171` for the site map that was declared in
  *neither* repo (latent because `ocean14` already has it, v0.25.0); `pytest.ini`
  `testpaths` widened to `hypernet/tests context/WHN`, since the default
  `hypernet/tests` collected 1 test and silently ignored the 18 WHN ones.
- **What I deliberately did not change:** the two outward-facing sentences that
  still point readers at `github.com/ocean-colour/IOPtics` under `context/WHN/`
  (`WATERHYPERNET.md:477`, `respond_to_kevin.md:160`). The second is a written
  promise to Kevin Ruddick about where the code will appear, so whether it gets
  edited or corrected in a follow-up note is JXP's call (Q&A M3). Likewise
  `out_root()` still points at `$OS_COLOR/IOPtics/whn_explore`, where the 18 MB
  parquet/npz cache lives (Q&A M2).
- **About the repo:** `hypernet` follows the house pattern — `requirements.txt` +
  `setup.py` (no `pyproject.toml`), `pytest.ini` at root, package dir
  `hypernet/`. CLAUDE.md puts package source under `hypernet/` and tests under
  `hypernet/tests/`, so the copied `context/WHN/` layout is a deliberate
  faithful-copy choice pending the M1 answer, not the repo's own convention.
- **Fable audit (read-only, run in parallel) confirmed the set is complete.**
  Independently swept the working tree, the filenames, and the branch history
  (including `git log -G` pickaxe over non-WHN paths) for the site codes and
  every WHN term; found nothing outside the copied set. The only other hits were
  genuine literature citations -- "Park & Ruddick (2005)" in `docs/context.md`
  and Doxaran/Ruddick/Nechad in `reports/gloria_fits_report.md`. It also
  confirmed `papers/*.pdf` is IOPtics-only (`IOPtics/papers/` holds
  `lange2020.pdf`, `moana_atbd.pdf`) and that the `.claude/settings.json` files
  carry no WHN entries.
- **Three things the audit caught that I had missed**, all verified by hand
  before acting:
  - **`tabulate` was undeclared** -- `whn_figures.py:123` calls
    `df.to_markdown()`, which needs it. Same latent-failure shape as cartopy
    (present in `ocean14` at 0.10.0, absent from both repos' requirements).
    Added to `requirements.txt` + `setup.py`.
  - **`make_pdf.py:26` hard-codes the macOS Chrome path** and drives it headless
    to render the PDFs, so regeneration is macOS-only. Recorded in `CLAUDE.md`.
  - **`CLAUDE.md` had gone stale** the moment I widened `testpaths` -- it still
    claimed `pytest.ini` pointed only at `hypernet/tests`, and its layout section
    had no `context/WHN/`. Rewritten, including the flat-import constraint
    (the scripts must stay in one directory) and the Chrome prerequisite.
  - Also flagged: stale *internal* references in the script docstrings
    (`whn_explore.py:8,28,34`, `test_whn_explore.py:3-4`) pointing at
    `waterhypernet_prompts.md`, `build_v1.py` and `ioptics/tests/conftest.py`.
    Left alone, posed as Q&A M6 so they can be swept in one pass with M1-M3.
- One audit claim I checked and did not act on: it reported
  `context/WHN/__pycache__/` as having been copied across. It was not -- the
  `rsync` excluded it, and the `.pyc` timestamps match my own pytest run. It is
  gitignored regardless.
- Re-ran `pytest -q` after the dependency and config edits: **19 passed**.
- No git state changed in either repo (JXP runs git).

### 2026-09-19 — Move prompt 2 (answers actioned; repo restructured; HOWTO doc seeded)

Built to the Q&A/Move answers. Opus 5 throughout, as asked; no subagents.

- **M1 restructure.** Split `context/WHN/` along the line JXP drew — doc-related
  code into `docs/`, generically useful code into `hypernet/` — and removed
  `context/` entirely. To `hypernet/`: `whn_explore.py` (the data layer) and
  `whn_simspec_check.py` (an archive-wide QC diagnostic, reusable independently
  of the write-up), with `test_whn_explore.py` to `hypernet/tests/`. To `docs/`:
  `WATERHYPERNET.md`, `respond_to_kevin.md`, `kevin_comments.md`,
  `whn_figures.py`, `make_pdf.py`, `figs/`, `summary_table.{csv,md}`,
  `simspec_flagged.csv`. Plain `mv`, not `git mv` (JXP runs git).
- **The split forced four functional edits**, none of them optional:
  - `fig_dir()` lived in `whn_explore.py` but anchors on `__file__`, so after the
    move it would have written figures into `hypernet/figs`. Only `whn_figures.py`
    ever used it, so it moved there and now resolves `docs/figs`.
  - `whn_figures.py` and `whn_simspec_check.py` imported `whn_explore` flat; both
    now import `hypernet.whn_explore`.
  - `test_whn_explore.py` dropped its `sys.path.insert` hack for
    `from hypernet import whn_explore as wx` / `whn_simspec_check as sc`.
  - `whn_figures.py` sits outside the package, so it needs the repo root on
    `sys.path`; added a three-line bootstrap that `pip install -e .` makes a
    no-op. `pytest.ini` went back to `testpaths = hypernet/tests`.
- **M2 done.** `out_root()` now resolves `$OS_COLOR/hypernet/whn_explore`, and
  the existing 18 MB cache was moved there rather than regenerated — all five
  artifacts (`index.parquet`, `pool.parquet`, `pool_spectra.npz`,
  `sample.parquet`, `simspec_check.parquet`) verified in place afterwards.
  `$OS_COLOR/IOPtics/` keeps its own `leaderboard.parquet` and `runs/`.
- **M3 done**, since the response had not gone out. `WATERHYPERNET.md` "How this
  was produced" now names the hypernet repo and the new run commands; the §14
  IOPtics link stays but reads "the companion IOPtics package". In
  `respond_to_kevin.md` the whole "once it is on `main`" promise was rewritten —
  the old text pledged a merge to IOPtics `main` that is now never going to
  happen, so it would have been wrong even with the URL swapped.
- **Verified the restructure rather than assuming it.** `pytest -q` → **19
  passed**. Both `python -m hypernet.whn_explore --help` and
  `python -m hypernet.whn_simspec_check --help` resolve. `python
  docs/whn_figures.py` regenerated all six figures into `docs/figs/` and rewrote
  the summary table — the per-site numbers came back identical to the committed
  ones, so the move changed nothing about the science.
- **M4 done.** Both PDFs regenerated from the corrected Markdown
  (`WATERHYPERNET.pdf` 1.78 MB, `respond_to_kevin.pdf` 0.15 MB) and left
  gitignored, per "we won't add the PDF to this Repo". `.gitignore` entry
  retargeted `context/WHN/*.pdf` → `docs/*.pdf`.
- **M6 deliberately not done** — JXP is saving the docstring sweep for the HOWTO
  effort, so the three stale references (`whn_explore.py:8,28`,
  `test_whn_explore.py:3-4`) are carried into `howto_prompts.md` as a listed
  task. I did update the docstring lines that the code changes made factually
  wrong (the `out_root` path, the figure output directory, the `-m` invocations),
  which is M2/M1 fallout rather than the deferred sweep.
- **`CLAUDE.md` and `README.md` rewritten** for the new layout: the package/docs
  split, run-from-the-root and `pip install -e .`, the `$OS_COLOR/hypernet`
  artifact location, and Chrome as a macOS-only prerequisite for PDFs.
- **Wrote `claude_prompts/howto_prompts.md`** — Goal, Conventions (corrected for
  this repo; the ones inherited here still describe IOPtics' Linux paths and
  Bokeh/Sphinx needs), Context, a proposed four-prompt sequence, Q&A and Logs.
  The Context section is the part worth reading: current layout, what the
  document actually is, the live threads (Kevin unsent, §14 written for the wrong
  audience, nothing readthedocs exists yet), and the numbers to keep straight.
  Flagged the genre problem explicitly — the document is a findings survey and a
  HOWTO is task-ordered, which is the first thing to settle.
- No git commands run; `$OS_COLOR/WATERHYPERNET` untouched.

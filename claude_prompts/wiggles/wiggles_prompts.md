# WATERHYPERNET -- Wiggles

## Goal

Attempt to improve the wiggles in the WATERHYPERNET dataset

## Conventions

- `ocean14`.  Run via `conda run -n ocean14 python ...`; `conda activate` fails
  non-interactively.
- **JXP runs git.**  Claude does not run any state-changing git command.
- Run scripts **from the repository root** so `hypernet` imports resolve
  (`python docs/whn_figures.py`, `python -m hypernet.whn_explore 1`).
  `pip install -e .` removes the need for the `sys.path` bootstrap at the top of
  `docs/whn_figures.py`.
- Run `pytest -q` after each step where relevant.  21 tests today; the archive-
  dependent ones skip themselves when `$OS_COLOR` is not mounted.
- **Tier 2** -- the figure and scan steps need the `$OS_COLOR` data tree mounted.
  Do **not** unset `$OS_COLOR`.
- Any calculation goes into a script on disk, not into the chat.
- Ask questions in the Q&A section below; log completed work under `## Logs`.

## Context

### Where things are 

- `hypernet/` -- the reusable data layer.
  - `whn_explore.py` -- `whn_root()`, `out_root()`, filename indexing,
    `load_spectrum`, pooling, shape clustering, the ~100-per-site sample.
    Stages: `python -m hypernet.whn_explore 1` (index) and `... 2` (pool,
    cluster, sample).
  - `whn_simspec_check.py` -- the Similarity-Spectrum over-subtraction scan.
    `python -m hypernet.whn_simspec_check [--sites all]`.
  - `tests/` -- `test_whn_explore.py` (19 tests over both modules, two tiers)
    and `test_import.py`; `pytest -q` runs 21 today.
- `docs/` -- the Sphinx source for the Read the Docs user guide (MyST
  Markdown): `index.md` plus the page files (`quickstart.md`, `rrs.md`,
  `over-subtraction.md`, ...).  The old single-file `WATERHYPERNET.md` has been
  split into these pages.  Every `.md` in `docs/` is published.
  - `whn_figures.py` -- writes `figs/*.png` and `summary_table.{csv,md}` here.
  - `make_pdf.py` -- **deprecated** (superseded by the RTD build).
  - `simspec_flagged.csv` -- the 30 flagged spectra from the SimSpec scan.
- `correspondence/` -- `kevin_comments.md`, `respond_to_kevin.md` (outside the
  docs build).
- `wavecal/` -- exploratory scripts for this effort (`sanity_checks_veit.py`,
  `line_fits_veit.py`, `instrument_timeline.py`, `select_sequences.py`).  Reusable code will go into
  `hypernet/`.  HAPI clone: `~/Oceanography/python/HAPI`.
- `context/` -- `kevin_wave.txt` (Kevin's covering email) and
  `ruddick2023.pdf` (the "SPIEwiggles" paper, Ruddick et al. 2023).
- Intermediates (parquet/npz, ~18 MB) live **outside** the repo at
  `$OS_COLOR/hypernet/whn_explore`; the archive itself is
  `$OS_COLOR/WATERHYPERNET/RELEASE_2` (L2A/L2B only -- no L1A).
- Kevin Ruddick's files for this task are on Google Drive at
  `AIOcean:data/Color/WATERHYPERNET/Wavelengths/` (read with `rclone`), mirrored
  locally to `$OS_COLOR/WATERHYPERNET/Wavelengths/`: his idea doc
  `XVwiggles_idea_2026-09-21.docx` and one VEIT sequence (L1A_IRR, L1A_RAD,
  L1C_ALL, L2A_REF netCDF + PNGs).
- **Lab wavelength-calibration data** (from Joel, 2026-10-17):
  `$OS_COLOR/hypernet/wiggles/wiggles_lab_data.zip`, unzipped to
  `.../wiggles/lab_data/hypstar_<sn>/calibration/wavelength/<YYYY_MM>/`
  (`fwhm_vnir.dat`, `deltawl_vnir.dat`, `*_wl_coefs_*.dat`, `lamp_vnir.dat`;
  7 instruments, 13 epochs, VNIR only).  Read by `wavecal/lab_data_examine.py`
  (outputs in `.../lab_data/derived/`).
- Radiative transfer: JXP's OSOAA fork, cloned at
  `/Users/xavier/Oceanography/python/RadiativeTransferCode-OSOAA`
  (`https://github.com/profxj/RadiativeTransferCode-OSOAA`).

## Prompts

### Setup

1. Read the Context above and `AIOcean:data/Color/WATERHYPERNET/Wavelengths/XVwiggles_idea_2026-09-21.docx` in full.  Then start a conversation with me about how we will proceed.  Put your questions in the Q&A section below.  Use Fable if you can.  Log your work.

2. I have answered your questions.  Read my responses and see if you need to 
ask any further questions.  If so, ask another round of questions.  Use Fable if you can.  Log your work.

3. I have answered your secound round of questions.  Read my responses and see if you need to 
ask any further questions.  If so, ask another round of questions.  Use Fable if you can.  Log your work.

4. I have answered your third round of questions.  Read my responses and move on to generate the planning doc.  Call it `docs/wiggles_planning.md`.  
Use Fable if you can.  Log your work.

5. The planning doc looks great.  Please generate a series of prompt docs in this `claude_prompts/wiggles` folder.  I think one per phase will work fine, with names `wiggles_phase<#>_prompts.md`.  
Have them use Opus 5.5.  Use Fable if you can for this task.  
Log your work.  And remind me where the note to Kevin is located.  
Or write it if you haven't yet.

6. Kevin responded to our note and will provide the data.  He was confused by one request -- the sky column in the data request.  Kevin wrote:

```
I’m not sure what you want for the “sky” in the csv, but I recommend to use Ld_750/Ed_750 as a good indicator of blue sky (low) vs clouds in the sky-viewing or sun-direction (high).

 

[Ruddick et al, 2006] equation (24) suggested to use the threshold of 0.05 for this parameter (but the Rhof(wind) formulae used in that equation is obsolete and should not be used).
```

I have added the Ruddick et al, 2006 paper to the `context/papers` directory.  Please read it and see if you agree with Kevin's suggestion.
If so, modify the prompt docs accordingly.  
Use Fable if you can.  Log your work.

### New data

1. Joel shared a zip file with new data.  It is the "wiggles_lab_data.zip" in `$OS_COLOR/hypernet/wiggles`.  Please examine it and discuss with me whether you can proceed with the undone phases of the plan.
Use Fable if you can.  Log your work.

## Q&A

### Setup #1 -- 2026-09-26 (drafted with Fable)

**What one VEIT sequence already says** (SEQ20260604T084543, HYPSTAR 122304;
`wavecal/sanity_checks_veit.py`, `wavecal/line_fits_veit.py`).  One sequence,
so read it as a lead, not a result.

- L1C `irradiance` is exactly `np.interp(wL, wE, mean L1A_IRR)`.  The pipeline
  does what Kevin says.
- E grid: 1536 px, L grid: 1538 px, spacing 0.465-0.498 nm.  The sub-pixel
  offset between them is **not** constant; it drifts through a full pixel
  across the spectrum.  `bandwidth` = 3.0 nm for every pixel, which looks like a
  placeholder.  About 6 samples per FWHM means the spectra are well sampled,
  the opposite of the Chance et al. (2005) undersampling case.
- The error from linear interpolation of Ed (cubic minus linear) is small:
  rms 3-6e-4 relative over 400-700 nm, peaking at 0.6 % at Ca H/K and 1.8 % at
  O2-A.
- The Fraunhofer structure left in Ld/Ed is about **10x larger** than that.
  Linear interpolation accounts for only ~3 % of the line structure; ~30 % of it
  survives the ratio.  A rigid shift of the E grid doesn't help: the best
  shift is +0.1 nm, with negligible gain.
- **The spectral response functions (SRFs) of E and L differ.**  Gaussian FWHM
  of the same solar lines, Ed vs Ld (nm): Ca K 1.8 vs 2.5, G band 2.5 vs 3.3,
  Hβ 3.0 vs 3.4, Na D 2.6 vs 3.1, Hα 2.8 vs 3.0, O2-A 3.4 vs 3.5.  E is sharper
  than L in the blue, and Lu matches Ld.
- The E-vs-L wavelength scales already agree to about 0.03-0.08 nm
  (0.1-0.2 px) at the clean lines, and the absolute scale is within ~0.05 nm
  (Hα, Na D).
- Kevin's symptom is real: corr(ρw'', Ed''/Ed) = -0.56 over 400-700 nm.

**My reading (critical-partner mode):** on this sequence the wiggles come
mainly from the **E/L SRF mismatch**, not from interpolating Ed.  Kevin's idea
still works if we generalise it.  We would replace "interpolate Ed to the L
grid" with "fit a high-resolution solar × atmosphere template to Ed, then
evaluate it with the L SRF on the L grid."  That fixes both the interpolation
and the SRF mismatch, and it is an SRF/wavelength calibration exercise, which
fits this branch's goal.

**Proposed phases.**

- **0: Diagnosis.**  Line centroids and widths for E/L, plus a wiggle budget,
  over ~100 sequences, which needs L1A data.
- **1: Twin experiment.**  TSIS-1 HSRS × atmospheric transmittance, smooth ρw
  library with fluorescence/Raman controls.  Cases: grid offset only, SRF
  mismatch, wavelength-cal error.  Methods: linear, cubic/sinc, template
  forward model.
- **2: Code.**  An `Ed_E → Ed_L` function plus self-calibration of the SRFs from
  the Fraunhofer lines, in `hypernet/` and as a drop-in snippet.
- **3: In-situ test.**  Compare ρw'' before and after.

**Questions for JXP**

1. **Reframe now or later?**  Do we tell Kevin now that the SRF mismatch looks
   dominant?  Or do we first run his twin experiment as written, as a control
   that tests his hypothesis?
>A. Yes, I will send him an email.  But I think we should also run the twin experiment.  No harm in that.

2. **What does "improve the wavelength calibration" mean here?**  On this
   instrument E and L agree to ~0.05 nm.  Was the goal prompted by something
   seen elsewhere (other instruments/sites, the SimSpec scan)?  Or is SRF
   characterisation the real target?  The answer decides whether Phase 0 is a
   side check or the main event.
>A. Fair point.  If the wavelength calibration is spot on then the focus will become the SRF

3. **L1A data.**  RELEASE_2 locally holds only L2A/L2B, where Ed is already on
   the L grid.  L2A can still support the SRF diagnosis from the Lu/Ld lines,
   but the interpolation test needs L1A_IRR + L1A_RAD.  Do we ask Kevin for a
   bulk L1A pull?  Which sites/instruments, and how many sequences?
>A. Yes, we are going to need to ask for those data.  We will develop a plan doc that includes a table of data needed and the sites/instruments/sequences.

4. **Calibration files and measured SRFs.**  Can Kevin supply
   `HYPERNETS_CAL_HYPSTAR_*_{RAD,IRR}_v2.3.nc` and any lab line-spread data?
   Also, why do E and L differ in width at all?  It could be the diffuser vs
   fore-optics filling the slit differently, or two spectrometer configurations.
   If no lab SRFs exist, we estimate them from the solar lines.  That estimate
   becomes part of the deliverable.
>A. Yes, he can provide these.  I believe the E and L data traveled through different paths through the spectrometer.  And it's all about angles.  It's why I also suspected the wavelength calibration could differ.

5. **The SPIEwiggles paper** (eqs. 14-15) is not in the Drive folder.  Can you
   get it from Kevin?  We need it to test *his* equations rather than our
   reconstruction of them.
>A. I have downloaded a copy and put it in `context/papers/ruddick2023.pdf`

6. **High-resolution reference.**  I propose TSIS-1 HSRS (Coddington et al.
   2021) as the solar spectrum.  For telluric O2/H2O/O3, is a simplified
   transmittance model OK at 3 nm, or do you want libRadtran/MODTRAN (neither
   is known to be installed)?
>A. Let's do both, if there is an open source option.   I have a fork of a rad transfer code here: `https://github.com/profxj/RadiativeTransferCode-OSOAA`.  Let's use that if we can. I just cloned a copy to `Oceanography/python/RadiativeTransferCode-OSOAA`.

7. **Ring effect and real features.**  Ld (and therefore ρf·Ld) carries
   rotational-Raman filling-in of the Fraunhofer lines that Ed lacks.  Options:
   (a) ignore it, (b) model it, or (c) diagnose with Lw/Ed only.  Separately,
   the chlorophyll-fluorescence and water-Raman controls, which the method must
   not erase, would be built into Phase 1.  Agreed?
>A. Agreed.
 
8. **Deliverable and authorship.**  (a) A Kevin-led paper, (b) a JXP-led methods
   paper, or (c) code + tech note?  And for the code: `hypernet/` plus a
   standalone snippet in the style of `interpolate_wav_linear.py`, with no
   attempt to run `hypernets_processor` end-to-end.  OK?
>A. (c).  And yes, we will put the code in the `hypernet/` directory.

9. **Uncertainty.**  Should the new step carry a FIDUCEO-style uncertainty
   term (De Vis et al. 2024) in the snippet, or is that left to the processor
   team?
>A. Yes

10. **Housekeeping.**  The "Where things are" context above is stale:
    `kevin_comments.md` and `respond_to_kevin.md` now live in
    `correspondence/`, `WATERHYPERNET.md` has been split into the Sphinx pages,
    and there are 19 tests, not 18.  Shall I update it?  I also mirrored Kevin's
    Drive folder to `$OS_COLOR/WATERHYPERNET/Wavelengths/` and added a
    top-level `wavecal/` directory for the scripts.  Say if you want them
    elsewhere, e.g. `hypernet/`.
>A. Sure, update it

### Setup #2 -- 2026-09-27 (drafted with Fable)

**What I learned since round 1**

- **Ruddick et al. 2023, eqs. (14)-(15).**
  - Eq. 15 convolves a high-resolution model irradiance Emod(λ) with the
    **E** SRF at each E pixel.
  - Eq. 14 linearly interpolates the ratio E_i / Emod_E(λ_E^i) (measured over
    model) and multiplies back by the model at the L centre.
  - Inputs: Emod(λ) at sub-SRF resolution (TSIS-1 HSRS × atmospheric
    transmittance for the SZA), the per-pixel E SRF, and both grids.
  - The paper never tests it: §3.2 sets Emod = the measurement envelope, which
    makes the result trivially perfect.
  - **Because the E SRF appears on both sides, eq. 14 fixes a grid offset but
    cannot fix an E/L SRF-width mismatch.**  Every simulation in the paper
    uses equal 3-nm widths with centre offsets of 0-0.75 nm.  The paper does
    mention SRF differences (§1.2, §1.5, Table 1), but only ever models them
    as offsets.
- **Our generalisation is a one-symbol change.**  Put the **L**-SRF-convolved
  model in the numerator:
  `Ed_L(λ_L^j) = Emod_L(λ_L^j) · [(1-w)·E_i/Emod_E(λ_E^i) + w·E_{i+1}/Emod_E(λ_E^{i+1})]`.
  - Setting ω_L = ω_E gives exactly Kevin's eq. 14.
  - A flat Emod gives plain linear interpolation.
  - So one function covers all three methods, and the twin experiment can
    compare them on one code path.
- **This squares with round 1.**  Our measured E-L offsets (0.03-0.08 nm) are
  far smaller than the offsets the paper simulates.  His mechanism therefore
  predicts small wiggles here.  That matches our finding that linear
  interpolation explains only ~3 % of the residual line structure.
- **OSOAA** (JXP fork):
  - **No gaseous absorption at all**: `docs/science/atmosphere_model.rst`
    says "OSOAA V2.0 does not include explicit gaseous absorption modeling".
  - Strictly monochromatic, 299-1000 nm, taking ~30 s to 15 min per
    wavelength.  At 1e4-1e5 wavelengths that is days to weeks, which rules it
    out for sub-nm work.
  - It doesn't need to do that job: Rayleigh, aerosol and hydrosol terms are
    smooth.
  - Not built: no `exe/OSOAA_MAIN.exe`, and no `gfortran` on this machine
    (conda-forge has one).
  - The Python wrapper (`OSOAASimulation`) lives inside a notebook and
    `tests/test_water_Ed.py`, not in a module.
  - No Ring effect either, which fits the Q7 decision.
- **Open-source telluric options:**
  - **HAPI** (`hitran-api`, pure Python; HITRAN O2/H2O lines on any grid).
    Add O3 Chappuis cross-sections from tables.
  - **libRadtran** (build with conda gfortran; REPTRAN fine ≈ 0.03-0.06 nm)
    as an independent check.
  - Py6S/6SV is too coarse (2.5 nm) except as a sanity check.
  - MODTRAN and TAPAS are not open.
- **Sample data:** L1A carries no uncertainty variables and no real SRF
  information (only the flat 3.0 nm `bandwidth`).  So a per-pixel E/L FWHM
  would be a **new calibration product**, not a fix to an existing one.
- Housekeeping: the paper is at `context/ruddick2023.pdf`, not
  `context/papers/`.  The Context section has been updated accordingly.

**Round-2 questions for JXP** (each has a default; say "defaults" to accept
them all)

1. **Division of RT labour.**  OSOAA would run on a coarse (5 nm) grid for the
   smooth direct, diffuse, Ld and Lw fields.  Its output gets multiplied by
   TSIS-1 HSRS F0 × HAPI gas transmittance on a 0.01 nm grid, with separate
   air masses for the direct and diffuse light.  The open question is whether
   libRadtran is built now, as the "both" cross-check.
   *Default: OSOAA + HAPI now; libRadtran later, and only if the two disagree
   with the data.*
>A. Yes, let's use OSOAA+HAPI.  Put the HAPI package in `Oceanography/python/HAPI`.

2. **Build OSOAA?**  That means installing conda-forge `gfortran` into
   `ocean14` (or a new env), running `make all`, and moving the notebook's
   `OSOAASimulation` class into a module.
   *Default: build in `ocean14` and put the wrapper in `hypernet/rt/`.*  Or
   would you rather the wrapper live in your OSOAA fork?
>A. I have added `gcc` and `gfortran` to this Mac.  I am pretty confident that
I have code out of a Notebook that works on OSOAA.  Check my other ocean colour Repos,
e.g. BING

3. **One method function.**  Implement the generalised eq. 14 with Kevin's
   eq. 14 (ω_L = ω_E) and linear interpolation as special cases.
   *Default: yes.*
>A. Yes

4. **SRF parameterisation.**  Options: a Gaussian whose FWHM is a smooth
   (quadratic) function of λ, fit per channel (E, L) per instrument to about
   10 clean solar lines; or a per-pixel / non-Gaussian shape.
   *Default: Gaussian + quadratic FWHM(λ); revisit when the cal files and any
   lab line-spread data arrive.*
>A. Let's go with your Default

5. **Is the SRF stable?**  You suspect E and L differ because of path and
   angles.  If the diffuser's angular fill changes with SZA or diffuse
   fraction, the E SRF could change from sequence to sequence, and a fixed
   correction would fail.  That makes Phase 0 (FWHM vs SZA / cloudiness /
   temperature / time) the gate for everything else.
   *Default: yes, design Phase 0 to test stability, and ask for clear and
   overcast sequences across SZA.*
>A. Let's go with your Default

6. **Size of the data request** (for the plan-doc table).
   *Default: 3 sites (VEIT + one turbid + one dark/clear), each with L1A_IRR +
   L1A_RAD + L1C_ALL + L2A_REF for ~50 sequences spanning SZA and sky
   condition (~150 in total).  Plus `HYPERNETS_CAL_HYPSTAR_*_{RAD,IRR}_v2.3.nc`
   for each instrument, any lab SRF data, and the processor version.*  Do you
   have preferred sites/instruments, e.g. more than one instrument at the same
   site, to separate site from instrument effects?
>A. That sounds fine, and I don't have any preferred sites/instruments.  Pick a diversity

7. **Uncertainty framework.**  The processor uses NPL's `punpy` /
   `comet_maths` / `obsarray` (De Vis et al. 2024), none of which is
   installed.
   *Default: adopt `punpy` so the snippet fits their chain, propagating SRF
   parameters, Emod and grid uncertainties into `u_rel_random` /
   `u_rel_systematic`.*
>A. Install all of those and add them to the `requirements.txt` file for this Repo

8. **Order and home of the plan doc.**
   *Default: start Phase 1 (twin experiment, no data dependency) and Phase 0
   on the VEIT sample now, while the data request is out.  Write the plan doc
   + data table in `correspondence/` (unpublished) and move it to `docs/`
   when you and Kevin agree on it, where it seeds the tech note.*
>A. Yes, it is fine to start those experiments before we have the data.  But wait until
we are done conversing and the plan is written.  I will ask for a separate prompt doc
to execute the Phase 1 and Phase 0 experiments.

### Setup #3 -- 2026-09-27 (drafted with Fable)

**Done since round 2**

- **HAPI** cloned to `~/Oceanography/python/HAPI`
  (github.com/hitranonline/hapi) and `pip install -e`'d into `ocean14`
  (hitran-api 1.3.0.0).
- **`punpy` 1.1.0, `comet_maths` 1.0.10, `obsarray` 1.0.3** installed in
  `ocean14`.  They are added to `requirements.txt` and `setup.py`, together
  with `hitran-api`.  `pytest -q`: 21 passed.

**What I learned**

- **The working OSOAA wrapper is not in BING.**  BING's only "osoaa" hits are
  bytes inside embedded PNGs in notebooks.  The code is in your OSOAA fork:
  - `notebooks/OSOAA_Downwelling_Irradiance_Tutorial.ipynb` and
    `..._Water_Leaving_Radiance_Tutorial.ipynb`, both executed on the Linux
    box with outputs saved.
  - The same `OSOAASimulation` class, pasted into `tests/test_water_Ed.py` and
    `tests/test_water_rt.py`.
  - A third copy, never executed, in
    `correct-atmosphere/nb/05_osoaa_toa_simulation.ipynb`.
- **What it covers.**  It turns a dict of keywords into an
  `OSOAA_MAIN.exe -KEY value` call via `subprocess`, one wavelength per call,
  in a temp work dir.  It parses two outputs:
  - `Flux.txt`: direct, diffuse and total Ed at every level, which is exactly
    the Ed split we need at 0+.
  - `LUM_vsVZA.txt`: upward radiance at one level, giving Lw at 0-.
  - **Nothing parses `LUM_Advanced_Down.txt`**, which is where Ld (sky
    radiance) comes from.  That parser is new but small (~30 lines).
- **Reuse means copy, not import.**  There are three near-identical copies,
  none in a module.  A merged `hypernet/rt/osoaa.py` (Flux + vsVZA + Adv.Down
  parsers, `OSOAA_ROOT` from the environment) would be ~250 lines.
- **Build looks routine.**  There is no `OSOAA_MAIN.exe` on this Mac yet.  The
  build is `cp gen/Makefile_OSOAA.gfortran Makefile; make all`: 12 fixed-form
  files.  gfortran 16.2 may need `-std=legacy`.
- **Level numbering conflict.**  The fork's docs give `OSOAA.View.Level` as
  ±1, but the working wrappers use 1 = TOA and 4 = 0-.  We need 0+ (3 in the
  original manual), so check `src/OSOAA_MAIN.F` first.
- **Instrument structure of the archive** comes from the netCDF attributes of
  every 20th RELEASE_2 file (`wavecal/instrument_timeline.py`; output in
  `$OS_COLOR/hypernet/wavecal/`).
  - **One instrument at two sites:** HYPSTAR_122302 at BEFR (2023-06 to
    2024-01), then THFR (2025-05 on), with the same cal file.
  - **Two instruments at one site:**
    - VEIT_H: 122304, then 122305 (2024-06 to 2025-07), then 122304 again.
    - GAIT_H: 121222, with a 120242 loan in 2025.
    - MAFR_H: 121231, then 122303 from 2026-05.
    - LPAR_H: 122307, then 120251.
  - **Same instrument before and after recalibration:** VEIT 122304
    (recalibrated 2024-11; Kevin's sample sequence is after it) and GAIT
    121222 (recalibrated 2025-03).
  - Water types: LPAR, MAFR and O1BE are turbid; THFR, BEFR and WRUK are dark;
    VEIT and GAIT are clear.  Every site reaches SZA 60-77°.

**Proposed data request.**  About 150 sequences in total.  Each sequence comes
with L1A_IRR, L1A_RAD, L1C_ALL and L2A_REF, and the set spreads over SZA,
season and sky condition.

| Site | Instrument(s) | N | Why |
|---|---|---|---|
| VEIT_H | 122304 (pre-recal), 122305, 122304 (post-recal) | 60 | two instruments + a recalibration at one clear site; Kevin's site, AERONET-OC alongside |
| BEFR_H + THFR_H | 122302 | 20 + 20 | one instrument, two dark lagoons, one cal file: separates site from instrument; dark water is where Ed-borne wiggles dominate ρw |
| MAFR_H | 121231 (+ 122303) | 30 | turbid, high-Lw regime; oldest cal (2021-10); 2026 instrument swap |
| GAIT_H (optional) | 121222 pre/post recal, 120242 | 20 | a second before/after-recalibration pair; dropped first if 150 is a ceiling |

We would also ask for:
- the cal files
  `HYPERNETS_CAL_HYPSTAR_{122302,122303,122304,122305,121231,121222,120242}_{RAD,IRR}_v2.3.nc`,
  with both cal versions for 122304 and 121222;
- any lab line-spread data;
- the processor version.

**Round-3 questions** (each has a default; say "defaults" to accept them all)

1. **When to build OSOAA.**  Building is infrastructure rather than an
   experiment, but you asked for nothing to run until the plan is written.
   *Default: make it task 0 of the Phase 1 prompt doc.  Build in place in the
   fork, set `OSOAA_ROOT`, and use the fork's `tests/test_water_Ed.py` as the
   smoke test.*
>A Use the default
2. **Where the wrapper lives.**
   *Default: a new `hypernet/rt/osoaa.py` that merges the two parsers and adds
   the Adv.Down parser.  Nothing is imported from the fork's `tests/` or
   notebooks, and we offer it back to the fork later.*
>A Use the default
3. **Sites.**  *Default: the four rows above.  If Kevin prefers three sites,
   drop GAIT; LPAR is the turbid alternative to MAFR.*
>A Use the default
4. **Who picks the sequences.**  We can choose sequence IDs from the L2
   archive by instrument period, SZA and season.  L2 does not tell us sky
   condition, though.
   *Default: we send Kevin a candidate list of IDs per row and ask him to tag
   them clear/overcast or swap some out.  The data table then lists IDs, not
   just counts.*
>A Use the default
5. **Shape of the plan doc.**  *Default: `correspondence/wiggles_plan.md`,
   written for Kevin, with these sections:*
   - *Problem & two hypotheses: H1 grid offset (Kevin) and H2 SRF mismatch
     (ours), both tested on one code path.*
   - *Phases 0-3 with gates.*
   - *Data request table, plus `wiggles_data_request.csv`.*
   - *Reference spectra & RT: TSIS-1 HSRS, HAPI, OSOAA on a 5 nm grid.*
   - *Method: the generalised eq. 14.*
   - *Deliverables & uncertainty (punpy).*
   - *Timeline.*

   *VEIT_P (PANTHYR) gets one sentence as a cross-system L2 control, with no
   data requested.*
>A Use the default

### Setup #6 -- 2026-10-03 (drafted with Fable): Kevin's sky index

**Verdict: agree, with two modifications.**

- **What eq. 24 is.**  Ruddick et al. (2006, L&O 51:1167; p. 1173) switch
  the sky-glint coefficient on L0+sky(750)/E0+d(750).  Below 0.05 the sky is
  clear and ρ_sky follows the wind (eq. 23).  At or above 0.05 it is cloudy
  and ρ_sky = 0.0256 (eq. 24).  The text quotes ~0.02 for clear skies
  (Mobley 1999) and ~0.3 for full overcast; an isotropic sky gives 1/π.
  Kevin's caveat holds: the obsolete ρ(wind) formula sits only in eq. 23's
  glint term, and the index uses measured Ld and Ed with no ρ.
- **We can compute it ourselves.**  Ld and Ed are in L1A, L1C_ALL, L2A_REF
  and the Release 2 L2B we already hold, so Kevin need not fill anything.
  VEIT sample: 0.0270.  Over the 224 requested sequences (L2B): 10/50/90 %
  quantiles 0.013 / 0.024 / 0.114, max 0.313, and 78 % below 0.05
  (`wavecal/ld_ed_750_request.py`).
- **The 0.05 threshold transfers to HYPSTAR.**  Every requested sequence is
  at relative azimuth 90°, while the 2006 paper used 135°.  Even so, the
  clear floor (10th percentile, 0.011-0.015) is flat with SZA, so one
  threshold suffices.  Every SZA bin holds both clear and cloudy cases.
- **Modification 1: keep the number, not just a tag.**  Phase 0 regresses
  FWHM and offsets on the continuous `ld_ed_750`.  The classes are for plots:
  clear < 0.05, cloudy 0.05-0.25, overcast ≥ 0.25.
- **Modification 2: add a broken-cloud flag.**  The paper says neither
  branch fits partly cloudy skies, and it dropped those stations.  The index
  also mixes "cloud in the sky view, sun clear" with "sun covered, blue
  view", and for the diffuser/SRF question it is the sun direction that
  matters.  So we add the within-sequence Ed scan CV at 750 nm from L1A
  (`broken` if > 2 %; VEIT has 0.46 %).  We do not use the difference
  between the two sky series as a flag: it was already 15 % on the clear
  VEIT sample.
- **Edits made:**
  - Phase 0: Context → Sky index; new task 8e; tasks 9 and 11; G0(b).
  - Phase 1 task 5: RAA 90° and a clear-floor check.
  - Phase 3: Context, metrics by sky class, and a Ring covariate.
  - `docs/wiggles_planning.md`: §3, §5, the timeline and the references.
  - `docs/wiggles_data_request.csv` is unchanged (as sent).

**Questions for JXP** (defaults stand unless you say otherwise)

1. Thresholds: clear < 0.05 (Kevin's), overcast ≥ 0.25, broken if Ed CV > 2 %?
>A. Yes, those thresholds are fine.
2. Shall I draft a one-paragraph reply to Kevin (in `correspondence/`) saying
   we adopt his index and compute it ourselves, so he only needs to swap
   spares?
>A. Yes, draft a short reply.  Then hold tight.

### New data #1 -- 2026-10-17 (drafted with Fable): Joel's lab data

**What the zip is.**  It holds no field data.  For each of our seven
instruments it has the **lab wavelength calibration** from 1 to 4 epochs,
13 epochs in all.  Each epoch has the same four files:
- per-line FWHM of the L and E channels at Hg/Ar/Xe (+Kr) lamp lines;
- the 4th-order pixel→λ polynomial for each channel;
- the lamp spectra seen by both channels (one L/E column pair per lamp,
  peak-normalised);
- `deltawl`.

This is our ask #2 (lab slit-function and wavelength-calibration data).
There are no L1A sequences and no `HYPERNETS_CAL_*.nc` files, so Phase 0b
tasks 9-12 and 14 and all of Phase 3 are still waiting on the delivery.

**What it says** (`wavecal/lab_data_examine.py`; I re-ran it and checked the
key numbers and the FWHM figure)
- **The lab confirms H2, instrument by instrument.**  Median ΔFWHM = L − E
  over the 400-680 nm lines (lab | our field table):

  | instrument | lab | field |
  |---|---|---|
  | 120242 | +0.44 | 0.39 |
  | 122302 | +0.47 | 0.46 |
  | 122304 | +0.39 | 0.39-0.49 |
  | 121222 | −0.04 | −0.05 |
  | 121231 | +0.02 | 0.08 |
  | 122305 | −0.06 | 0.01 |
  | **122303** | **−0.02** | **0.14** |

  The srf/linear split matches exactly, and 122303 is clearly `linear`
  (it is currently flagged `method_sensitive`).  Our field FWHM(λ) models
  sit within ~0.1 nm of the lamp widths, channel by channel.  We had thought
  E from L2 read 0.06-0.09 nm too wide; the lamp does not support that.
- **The narrow-E instruments have flat-topped E SRFs.**  Super-Gaussian
  exponent p ≈ 2.5-2.9 for E on 120242, 122302, 122303, 122304 and 122305,
  against p ≈ 2.1 (Gaussian) for L everywhere and for E on 121222/121231.
  χ²_ν drops 2-5× for E with a super-Gaussian.  There are no Lorentzian
  wings (pseudo-Voigt η = 0), and the profiles are slightly wider on the
  blue side.  Phase 2 task 9b found that `srf` keeps 76-84 % of its gain
  for a flat-topped SRF, so that is the relevant row.
  - Caveat: our Fraunhofer fit on VEIT gave p_E = 2.23.  Lamp light filling
    the diffuser need not mimic sky light, so the field shape still matters.
- **We now have every instrument's own E wavelength grid.**  The processor's
  grids are the lab polynomials cropped to 350 ≤ λ < 1100 nm.  This
  reproduces the VEIT L1A E and L grids to 2e-5 nm, and the Release 2 L grid
  of every instrument (one red-end pixel is uncertain on 120242).
  - The 122304 proxy we used for the other instruments' E grids is off by
    0.15-1.3 nm (up to more than 2 pixels).  It also has the wrong E-pixel
    phase relative to L.
  - Own-grid E − L offsets at 400/550/700 nm: about −0.5/−0.4/−0.3 nm on
    the narrow-E units, and −0.1 to +0.01 nm on 121222/121231.
- **Release 2 uses one wavelength calibration per instrument, not per
  period.**  122304 after its 2024-11 recalibration still runs on the
  2023_01 polynomial; the 2024_09 lab polynomial differs by 0.02-0.05 nm
  rms, up to 0.10 nm.  122302's 2024_09 lab epoch is not used either.
  121222 has no lab epoch for its 2025-03 recalibration.
- **The scale is air.**  All 19 reference lines are NIST air wavelengths.
  This settles interim-report ask #7.
- **A common lab polynomial error.**  The 4th-order fit misses the same
  lines by the same amounts on all seven instruments and both channels:
  +0.09 nm at Hg 404.7, −0.09 at Hg 546.1, ~±0.1 nm in the NIR.  It matches
  our 0.13 nm "common window pattern" (Phase 2 open concern #1) near 545 and
  585 nm, but not at 405/435.  Between 440 and 690 nm only one lamp line
  constrains it.
- **Stability on the lab side.**  Between epochs, FWHM moves ≤ 0.12 nm and
  the polynomial ≤ 0.09 nm (120242, over four years).  122302 and 122304
  change by ≤ 0.04 nm from 2023 to 2024.  All of this is well under G0(b)'s
  0.2 nm.
- **Stray light.**  The Xe floor is ≤ 3e-4 of a bright line.  The near
  wings at ±3 nm are 1 % on the red side and 3-5 % on the blue side.
- **`deltawl`** is reference λ − nominal-grid λ.  The nominal grid is the
  `wl` column of `lamp_vnir.dat`, a fixed per-instrument (probably factory)
  polynomial.  The model fits to 0.002-0.013 nm rms.  It is not a fit
  residual.

**A route around the missing L1A** (`wavecal/invert_l2_irradiance.py`).
L1C/L2 `irradiance` is linear interpolation of the mean E onto the L grid:
E_L = M·E, where M depends only on the two grids.  We now know both grids,
so E on its native grid can be recovered by least squares.
- On VEIT, recovering E from L1C with the lab 2023_01 grid matches the true
  L1A E to 1.5e-5 rms (5e-5 max) over 400-900 nm.  Its second differences
  match to 0.2 % of the E signal.  cond(M) = 42.
- **It is only as good as the grid.**  The 2024_09 polynomial (0.05 nm off)
  gives an 18 % error in the second differences.  So we must use the
  polynomial the processor actually used.  That is verified for 122304
  only; for the others we infer it from the matching L grid.
- If this holds, we could do several things now on the Release 2 L2B files,
  which carry E_L, Ld, Lu, Lw and ρw:
  - fit E widths on the native grid with the correct grid, not the 122304
    proxy, so no L2 E-width bias is left to estimate;
  - run Phase 0b tasks 10-12 in L2 form on the 224 requested sequences
    (beyond 8c);
  - run Phase 3's in-situ test (recompute ρw with `linear` /
    `ruddick2023` / `srf`).

  Things that still need L1A: per-scan data (`ed_cv_750`, scan-level fits),
  the Lu/Ld scan splits, the cal `.nc` files, and a check of the grid
  inference on instruments other than 122304.

**Can we proceed?**  Partly.  What I propose, in order:
1. **Close the lab half of task 13.**  Promote the reader to a module
   (default `hypernet/labcal.py`).  Commit the lab table
   (`hypernet/data/lab_srf.csv`: per instrument × epoch × channel × line,
   FWHM, super-Gaussian p, the polynomial) and per-instrument E/L grids.
   Add tests.
2. **Revise the Phase 2 SRF table.**
   - Use each instrument's own E grid instead of the proxy, and refit the
     Release 2 E widths with it.
   - Set 122303 to `linear`.
   - Add an SRF-shape column (super-Gaussian p for E on the srf
     instruments), then re-run the VEIT check and task 9b with the lab
     shape.
   - Keep the field values as the source of truth for widths, with the lab
     as validation.
3. **Revisit the wavelength-scale correction.**  Test whether the lab's
   common polynomial error explains the common window pattern, and if so
   replace the empirical pattern with it.
4. **Test the L2 inversion on the archive.**  Check that the recovered E is
   physically sensible across instruments, e.g. that E-line fits on
   recovered E agree with the lab FWHM.  If so, run Phase 3 on the 224
   requested sequences from L2B as a *provisional* G3, and re-run it when
   the L1A arrives.
5. **A provisional G0 report** (`gate_G0.md`) from Release 2 + lab: (a)
   confirmed, (b) passes, (c) air scale, offsets ≤ 0.1 nm.  It would be
   marked provisional until the L1A is fitted.

**Questions for JXP**
1. OK to proceed with 1-5 in that order?  Should 2-4 reopen Phase 2 / start
   Phase 3 in their own prompt docs, or go here?
>A. Yes, proceed with 1-5 in that order.  And add a step 6 which is to udpate the `wiggles_interim_report.md` with the new analysis.  Also, please generate a new prompt doc named `wiggles_lab_data_prompts.md` which will be used to prompt the agent to perform 1-5 + step 6.
2. Item 4 runs Phase 3 before the L1A delivery, on reconstructed E.  Do you
   want that, or would you rather wait for the L1A?
>A. Let's wait for the L1A.
3. Questions back to Joel/Kevin:
   - Why does the processor keep the 2023_01 wavelengths after the
     2024-11 recalibration of 122304, and likewise for 122302?
   - Is there a 2025 lab epoch for 121222?
   - What is the nominal `wl` grid in `lamp_vnir.dat`?
   - Is the lamp illuminating the diffuser the same way the sky does?
   - The L1A delivery and cal files are still outstanding.

   Shall I draft that note?
>A. Yes, draft it.  I will send by email, so use plain text, not markdown.
4. Where should the per-instrument grids live: `hypernet/data/` (≈10 small
   CSVs) or computed on the fly from the committed coefficients?  Default:
   commit the coefficients and compute the grids.
>A. Use your default.

## Logs

### 2026-09-26 -- Setup #1

- Read the Context and Kevin's `XVwiggles_idea_2026-09-21.docx` (fetched with
  `rclone` from `AIOcean:`) plus `context/kevin_wave.txt`.  Kevin's
  hypothesis: linear interpolation of Ed onto the L grid
  (`hypernets_processor/.../interpolate_wav_linear.py`) creates spurious ρw''.
  His proposal is a twin experiment followed by an in-situ test.
- Mirrored the Drive folder (L1A_IRR, L1A_RAD, L1C_ALL, L2A_REF netCDF + PNGs
  for VEIT 2026-06-04) to `$OS_COLOR/WATERHYPERNET/Wavelengths/`.
- A Fable subagent drafted the plan and questions and wrote two diagnostic
  scripts.  I re-ran them and got the same numbers, then saved them to the new
  top-level `wavecal/` directory (data path taken from `$OS_COLOR`):
  - `wavecal/sanity_checks_veit.py`: L1C = linear interp check, interpolation
    error, rigid-shift scan, ρw''-vs-Ed'' correlation.
  - `wavecal/line_fits_veit.py`: Gaussian centroids, widths and depths of 10
    solar/telluric lines in Ed, Ld and Lu.
- Learned: the E (1536 px) and L (1538 px) grids differ, with a drifting
  sub-pixel offset.  `bandwidth` is a flat 3.0 nm placeholder.  Local RELEASE_2
  has no L1A, so Ed there is already interpolated.  Kevin's cal-file attribute
  bug is present in L1A_RAD, L1C and L2A.
- Nothing run through pytest (no package code changed).
### 2026-09-27 -- Setup #2

- Read JXP's round-1 answers.  Decisions: run Kevin's twin experiment as a
  control as well; if the wavelength calibration is fine, the focus moves to
  the SRF; Kevin will supply L1A data and cal files; deliverable is code in
  `hypernet/` + a tech note (option c), with FIDUCEO-style uncertainty; use
  JXP's OSOAA fork for the RT.
- Updated the stale Context "Where things are" section (Sphinx split,
  `correspondence/`, `wavecal/`, `context/`, local mirror of Kevin's folder,
  OSOAA clone) and the test count (21 pass today).
- A Fable subagent read Ruddick et al. 2023 in full and surveyed the OSOAA
  fork; I spot-checked its key claims (the no-gas-absorption note in
  `docs/science/atmosphere_model.rst`, no `gfortran`, the eq. 14/15 text).
  Learned: eq. 14 uses the E SRF on both sides, so it only fixes grid offsets;
  our L-SRF generalisation contains it as a special case.  OSOAA is
  monochromatic, gas-free, unbuilt, and too slow for sub-nm grids.  It is
  suited to the smooth fields only; HAPI (or libRadtran) is needed for
  telluric lines.
- Wrote the Setup #2 findings and 8 round-2 questions (each with a default)
  to the Q&A section.  No code changed; `pytest -q`: 21 passed.

### 2026-09-27 -- Setup #3

- Read JXP's round-2 answers: all defaults accepted.  Use OSOAA + HAPI, with
  HAPI in `Oceanography/python/HAPI`; install the NPL CoMet packages and add
  them to requirements; pick a diverse set of sites; no experiments until the
  plan is written, since a separate prompt doc will cover Phases 0/1.
- Cloned HAPI and installed it editable into `ocean14`.  Installed `punpy`,
  `comet_maths` and `obsarray`.  Added all four to `requirements.txt` and
  `setup.py`.  `pytest -q`: 21 passed.
- A Fable subagent searched for your OSOAA code and surveyed instrument/site
  structure in RELEASE_2.  I spot-checked two claims: BING has no OSOAA code,
  and HYPSTAR_122302 appears at BEFR.  I saved its survey script as
  `wavecal/instrument_timeline.py`, with output redirected to
  `$OS_COLOR/hypernet/wavecal/instrument_timeline.csv`.
- Learned: the OSOAA wrapper lives in the fork's notebooks/tests, and an Ld
  parser is missing.  Instruments moved between sites and were recalibrated,
  which gives natural site-vs-instrument and before/after-cal controls.  Wrote
  Setup #3 findings, a draft data-request table and 5 round-3 questions.

### 2026-09-28 -- Setup #4

- JXP accepted every round-3 default.  The prompt doc has been renamed
  `wiggles_prompts.md`, and the paper now lives at
  `context/papers/ruddick2023.pdf`.
- Wrote `wavecal/select_sequences.py`, which builds the data-request table
  `docs/wiggles_data_request.csv`: 224 rows, 150 primary + 74 spare.
  - It derives the site x instrument x calibration periods from
    `instrument_timeline.csv`.
  - Within each period it samples from the RELEASE_2 index, spread over SZA
    tertiles and month, with SZA <= 75 deg (an earlier pass picked sequences
    up to 88 deg).
  - It reads SZA and `system_id` from each file, so every pick is checked
    against the intended instrument.  There are no duplicates.
  - Most picks are L2B files, because the archive holds more L2B than L2A.
    Kevin can identify sequences by site + `sequence_time` + azimuth either
    way.
  - The `sky` column is blank, for Kevin to fill.
- A Fable subagent wrote `docs/wiggles_planning.md` from the Q&A, Kevin's doc,
  the paper and the scripts.  It re-ran the VEIT scripts to quote their
  numbers.  I reviewed the whole draft and fixed three things:
  - Figure 4 of the paper used TSIS v2 at 1 nm, interpolated to 6SV's 2.5 nm,
    not at native resolution.
  - The De Vis et al. 2024 author list had two authors added from memory; it
    is now trimmed to "et al.".
  - The PANTHYR description now names its TriOS RAMSES radiometers.
- Equations use MyST ```` ```{math} ```` blocks, because `dollarmath` is not
  enabled.  The Sphinx build succeeds.  Its only new warning is that
  `wiggles_planning.md` is not in any toctree.  The page is still built and
  served (`wiggles_planning.html`), and the CSV is published as a download.
- Open items in the doc: the SPIE volume/paper number of Ruddick et al. 2023;
  a citation for CoMet; the exact interface of `interpolate_wav_linear.py`
  (not re-read from GitHub).

### 2026-09-28 -- Setup #5

- A Fable subagent wrote four phase prompt docs in `claude_prompts/wiggles/`,
  following the house structure (Goal / Conventions / Context / Prompts /
  Q&A / Logs) and taking scope from `docs/wiggles_planning.md`:
  - `wiggles_phase0_prompts.md`: 14 prompts.  0a covers SRF fitting on the
    VEIT sample now; 0b covers the full request once L1A arrives; then gate
    G0.
  - `wiggles_phase1_prompts.md`: 12 prompts.  OSOAA build and
    `hypernet/rt/osoaa.py`, Emod builder, twin experiment cases (i)-(vii),
    gate G1.
  - `wiggles_phase2_prompts.md`: 10 prompts.  `interpolate_ed_to_l`,
    self-calibration, punpy uncertainty, drop-in snippet, exit check.
  - `wiggles_phase3_prompts.md`: 12 prompts.  In-situ test, Ring check, the
    PANTHYR control, gate G3, the tech note.
- Each doc carries a **Model** bullet: use Opus 5.5 (`claude-opus-5-5`),
  including for subagents.  Gate reports go to
  `claude_prompts/wiggles/gate_G*.md`, never `docs/`.  Open decisions (module
  names, figure locations and so on) are left for JXP in each Setup prompt.
- I corrected the Phase 2 doc.  It claimed no `hypernets_processor` commit was
  pinned, but Kevin's idea doc links commit `9a12819a3ffb...`, so Task 2 now
  fetches that commit.
- Wrote the draft note to Kevin, `correspondence/wiggles_note_to_kevin.md`
  (not sent).  It gives the one-sequence finding, the generalised eq. 14, the
  data/cal/SRF requests and the attribute-bug confirmation.  It attaches the
  plan doc and `docs/wiggles_data_request.csv`.

### 2026-10-03 -- Setup #6

- Kevin proposed Ld(750)/Ed(750) with a 0.05 threshold (Ruddick et al. 2006,
  eq. 24) in place of the `sky` column.  A Fable subagent read the paper in
  full and computed the index.  I checked the eq. 23/24 text, the
  0.02/0.3 quote and the partly-cloudy caveat in the PDF.  I also re-ran
  both scripts and got the same numbers.
- Saved the scripts as `wavecal/ld_ed_750_veit.py` (VEIT sample, variable
  inventory) and `wavecal/ld_ed_750_request.py` (all 224 request rows →
  `$OS_COLOR/hypernet/wavecal/ld_ed_750_request.csv`).
- Edited the Phase 0, 1 and 3 prompt docs and `docs/wiggles_planning.md`
  (see Q&A Setup #6).  Phase 2 needed no change.  Updating the code
  (`whn_l1a.sky_index`, `phase0b_index.py`, `phase0b_stability.py`) is now
  Phase 0 task 8e; no package code changed today.
- Learned:
  - All requested HYPSTAR sequences are at relative azimuth 90°, not the
    135° of the protocols.
  - L2B carries `downwelling_radiance` and `irradiance`, so the sky index
    is available archive-wide without L1A.
  - HYPERNETS L2 still uses the Mobley 1999 ρ (rhof 0.026-0.038).
- JXP's answers to Setup #6: the thresholds are accepted, and I am to draft a
  short reply to Kevin, then hold.  The reply draft is given in the chat for
  JXP and was not saved to `correspondence/`.  A task-8e run was started and
  then stopped before it changed any file, so 8e is still open.

### 2026-10-17 -- New data #1

- Read Joel's `wiggles_lab_data.zip` (52 files, 7 instruments, 13 lab
  wavelength-calibration epochs) and unzipped it to
  `$OS_COLOR/hypernet/wiggles/lab_data/`.  It holds no L1A and no cal
  `.nc` files.
- A Fable subagent wrote `wavecal/lab_data_examine.py` (965 lines; outputs
  in `.../lab_data/derived/`): the inventory, the file semantics, the grid
  crops, the FWHM comparison with `hypernet/data/srf_table.csv`, line shapes,
  stability and wings.  I re-ran it and checked the ΔFWHM table, the VEIT
  grid crop and the line-shape summary against its output, plus the FWHM
  comparison figure.
- Wrote `wavecal/invert_l2_irradiance.py`.  On VEIT, recovering native E
  from the L1C irradiance with the lab grid reproduces L1A E to 1.5e-5.
  The 2024_09 grid gives 18 % error in the second differences, so the grid
  must be the one the processor used.
- Wrote New data #1 in the Q&A: findings, what can proceed (the lab half of
  task 13, the SRF-table revision, the window pattern, the L2 inversion
  route to Phase 3, a provisional G0), and 4 questions.  No package code
  changed.
- Learned:
  - the lab data confirm the field E/L split on all seven instruments and
    put 122303 firmly in the `linear` group;
  - the narrow-E units have flat-topped E SRFs (p ≈ 2.5-2.9);
  - Release 2 wavelengths come from one lab polynomial per instrument;
  - the HYPSTAR scale is air;
  - a ~0.1 nm common polynomial error sits on every instrument.
  - The `egrid_122304_2024-11*.csv` written by the script uses the 2024_09
    polynomial, which the processor does not use.  Use the 2023_01 grid for
    that period.

### 2026-10-18 -- New data #1, answers

- JXP's answers:
  - proceed with items 1-5 in order, plus an item 6 (update
    `reports/wiggles_interim_report.md`), all driven by a new prompt doc;
  - **Phase 3 waits for the L1A**, so nothing runs on reconstructed E;
  - draft the note to Joel/Kevin in plain text for email;
  - grids are computed from committed coefficients (the default).
- Wrote `claude_prompts/wiggles/wiggles_lab_data_prompts.md` in the house
  structure (Opus 5.5).
  - Tasks:
    1. Setup;
    2. `hypernet/labcal.py` with `lab_wl_coefs.csv`, `lab_srf.csv`,
       `grid_for` and `PROCESSOR_EPOCH`;
    3. the SRF-table revision: own E grids replace the 122304 proxy in
       `hypstar_grids.npz`, phase0c, phase2_srf_table and phase2_wavescale;
       122303 → `linear`; `shape`/`shape_p` columns; re-runs of the VEIT
       check and task 9b;
    4. the window pattern vs the lab polynomial residual;
    5. validation of the L2 inversion only;
    6. a provisional `gate_G0.md`;
    7. the interim-report update, keeping the 2026-10-15 version.
  - Replaced committed products are kept as `*_pre_lab`.
- Wrote `correspondence/wiggles_note_lab_data.txt` (plain text, not sent).
  It thanks Joel, gives six findings and asks five questions: the
  wavelength epoch after recalibration, a 2025 lab epoch for 121222, the
  nominal `wl` grid, how the lamp illuminates the E diffuser, and the L1A /
  cal-file status.
- `pytest -q`: 137 passed.  No package code changed.

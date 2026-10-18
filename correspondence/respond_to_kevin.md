# Response to Kevin Ruddick's comments on the WATERHYPERNET Release 2 summary

**From:** J. Xavier Prochaska (UC Santa Cruz) · **Date:** 18 September 2026

Kevin — thank you for reading it so carefully, and for the context on the
network's design. Every one of your points has been acted on. Below: what
changed, what I found when I followed up on Cluster 5, and answers to your
dissemination question.

---

## 1. Common y axis on the per-site spectra — done

You were right that this is the fastest way to see which sites are dark and which
are bright, and the effect is striking: on one axis, Río de la Plata, the Gironde
and Oostende dominate, while Berre, Thau, Wraysbury — and, less obviously,
Acqua Alta and Lake Garda — are nearly flat. The network spans a factor of 12 in
median Rrs(560), 0.0030 sr⁻¹ at Wraysbury to 0.036 sr⁻¹ at the Gironde.

The one cost is that a common axis suppresses all spectral shape at the dark
sites, which is much of what the original figure conveyed. So the document now
carries **both**: the common-axis version as the primary figure, and the
per-panel autoscaled version immediately after it. Site order and colours are
unchanged between the two, so they read as a pair.

---

## 2. Cluster 5 — I followed this up, and I think it is a real finding

You wrote that Cluster 5 is "obviously not good" and asked how it passed QC.
Short answer: **it passed cleanly — every one of these spectra has
`quality_flag = 0` and 6/6 valid scans** — and the cause appears to be the
Similarity Spectrum correction over-subtracting in very dark water.

### What I did

Rather than rely on the 4,400-spectrum sample the clustering used, I scanned
**every file at the three darkest sites — all 8,947 of them**. The test, per
file, is a simple sign check over 440–600 nm:

> flag when `median(reflectance) ≤ 0 < median(reflectance_nosc)`

i.e. the corrected product has gone non-physical across the visible while the
uncorrected product is still positive. That can only mean the correction removed
more than the entire water signal.

### Result

| site | files | flagged | rate |
|---|---:|---:|---:|
| THFR_H (Thau) | 6,260 | 26 | 0.42% |
| BEFR_H (Berre) | 1,798 | 4 | 0.22% |
| WRUK_H (Wraysbury) | 889 | 0 | 0.00% |
| **total** | **8,947** | **30** | **0.335%** |

Characteristics of the 30:

- **All have `quality_flag = 0`; all have 6/6 valid scans.** Nothing in the
  30-bit mask fires, including `simil_fail`.
- **Low sun does not explain it.** Solar zenith spans 21.1–87.6°, median 54.4°,
  and only a third exceed 70°. Berre's four are all December 2023 at
  sza 85–87°, which does look like a low-sun winter effect — but Thau's 26 are
  spread across all months and all sun angles, including sza 21–31°.
- **Water darkness is the common factor.** For the flagged spectra the
  *uncorrected* ρw(440–600) median is 0.00043, against site medians of 0.0106
  (Berre), 0.0097 (Thau) and 0.0072 (Wraysbury) — about 20× darker than typical
  for the same site.
- **The subtracted offset exceeds the signal.** Median ρw removed is 0.00063
  against an uncorrected signal of 0.00043 — roughly **145% of the whole
  signal** — leaving a corrected median of −0.00019.
- They are spread over 2023 (4), 2025 (9) and 2026 (17), so this is persistent
  behaviour, not a single bad episode.

### Interpretation, offered tentatively

This looks like the mirror image of the caveat already in your release notes.
Those warn that `reflectance` is poor at the **turbid** sites (LPAR, MAFR, O1BE),
where the NIR signal is high. What I see is the opposite extreme: in the
**darkest** water, the offset the similarity correction subtracts becomes
comparable to or larger than the water signal itself, and the corrected spectrum
inverts. `reflectance_nosc` is positive and entirely plausible in all 30 cases,
so this concerns the correction rather than the acquisition.

Two caveats on my numbers. The 0.335% is a **lower bound** on affected spectra,
because it counts only cases where the visible median actually changes sign; the
shape-based clustering also flagged Wraysbury spectra that are badly distorted
without the median going negative, and those do not appear in the table above.
And I have only scanned the three dark sites — the bright sites show nothing of
this kind in the 400-spectrum-per-site sample, but I have not scanned them
exhaustively. I am happy to run the full HYPSTAR archive if that would be useful.

### The 30 files

The complete list, with solar zenith and the two medians, is in
`simspec_flagged.csv` alongside this note. The 10 worst by corrected median:

| site | acquisition | sza | median ρw corrected | median ρw uncorrected |
|---|---|---:|---:|---:|
| THFR_H | 2025-11-28 07:30 | 84.9 | −0.00238 | 0.00144 |
| THFR_H | 2025-08-21 08:15 | 55.1 | −0.00208 | 0.00198 |
| THFR_H | 2025-11-17 07:00 | 87.2 | −0.00120 | 0.00061 |
| THFR_H | 2026-07-20 10:15 | 30.2 | −0.00100 | 0.00078 |
| THFR_H | 2025-12-01 07:15 | 87.6 | −0.00087 | 0.00045 |
| BEFR_H | 2023-12-30 07:45 | 85.4 | −0.00070 | 0.00052 |
| THFR_H | 2026-05-28 09:30 | 34.8 | −0.00070 | 0.00116 |
| BEFR_H | 2023-12-14 07:30 | 86.4 | −0.00069 | 0.00402 |
| BEFR_H | 2023-12-26 07:45 | 85.3 | −0.00050 | 0.00024 |
| THFR_H | 2026-03-12 09:00 | 60.7 | −0.00045 | 0.00109 |

If it would help, the check is a single short script
(`whn_simspec_check.py`, ~150 lines, numpy + netCDF4) and you are welcome to
take it as-is for your own processing chain.

---

## 3. The single-parameter network — reframed, and thank you for the correction

My draft headed that section "The glaring absence", which framed a deliberate
design decision as a deficiency. That was not fair, and it is fixed. The section
is now **"Scope: reflectance only, by design"** and says so explicitly: HYPERNETS
is a single-parameter network conceived as a core that site PIs can expand as
resources allow, that a common multi-measurand protocol at every site is not
logistically or financially reachable, and that some PIs hold extra measurements
from short deployments.

I have also added your four applications for a reflectance-only archive —
satellite validation, local water-quality monitoring, aquatic-optics research
into spectral variability, and AERONET-OC cross-comparison — because a reader
coming to the archive cold deserves to see them stated.

What I have kept is the factual finding, because it is genuinely useful to a data
user: the variable and attribute sets are constant across all 11 sites (verified
over 440 files), so anyone hoping for co-located in-water measurements knows
immediately where they stand. The consequences specific to my own use case are
now confined to a clearly-labelled final section rather than framing the whole
document.

On your closing point: I agree, and it changes my own priorities. Automated flow
cytometry, turbidimeters and chlorophyll fluorimeters would be more directly
useful to me than IOPs — IOPs are, as you say, an intermediate explanatory
quantity rather than an end-user product, and chlorophyll and turbidity map
straight onto quantities a retrieval-evaluation framework already scores. If the
PIs do have such data from short deployments, even episodic coverage at one or
two sites would be valuable, and I would be glad to hear what turns up when you
ask them.

---

## 4. Dissemination — yes, with a preference for the order

Both options are appropriate, and I am glad you find it useful. My preferences:

**A frozen PDF first.** Authored "J. Xavier Prochaska" and dated, for the
WATERHYPERNET documentation. That carries no maintenance commitment on either
side and gives you something stable to cite immediately. A PDF of the current
version accompanies this note.

**A link to the repository second.** The code and document now live in their own
public repository,
[github.com/ocean-colour/hypernet](https://github.com/ocean-colour/hypernet) —
the document under `docs/`, the analysis code in the `hypernet` package. I will
send you the stable path once it settles, at which point an "exploitation tools"
link would be very welcome. Everything is re-runnable: indexing, sampling, the
figures and the over-subtraction scan, with a test suite that skips cleanly when
the archive is not mounted.

**On right of reply.** The document now contains two things that read as
criticism of your processing — the over-subtraction finding in §9.2 and the
"Hazards for an automated loader" list — and I would not want either published
under WATERHYPERNET's own documentation without your sign-off on the wording. You
have the final say on both. If you would rather the over-subtraction material was
raised privately with the processing team first and left out of the public
version until you have looked into it, say so and I will pull it from the
external release; it is your data and your reputation, and the finding is a small
one.

I would also be glad to add a line noting that the summary was reviewed by you,
if you are comfortable with that — it would tell readers the numbers have had a
knowledgeable eye on them. And naturally, if anything in the document is simply
wrong, tell me and I will correct it.

---

## 5. Small things

- The harmonisation gaps you mention are collected in §5 (schema differences),
  §4 (filename grammars) and §15 (the licence attribute present on HYPSTAR and
  absent on PANTHYR), so they sit together if that is useful for planning.
- The temporal-coverage figure you liked is now §3; it is generated directly
  from the filenames, so it will regenerate against any future release with no
  changes.

Thanks again for the careful reading — the Cluster 5 question in particular sent
me somewhere I would not otherwise have looked.

With best wishes,
J. Xavier Prochaska

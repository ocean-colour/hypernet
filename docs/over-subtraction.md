# Known issue: Similarity-Spectrum over-subtraction at the darkest sites

Clustering the spectra by shape isolated a small group that is negative across
most of the visible — not near-zero-with-noise, but systematically negative.
Following it up produced a finding worth reporting back to the network, since
**every one of these spectra passed QC with `quality_flag = 0`**.

The test: compare the median of `reflectance` and of `reflectance_nosc` over
440–600 nm. Where the corrected product is ≤ 0 while the uncorrected one is still
positive, the Similarity-Spectrum correction has removed more than the entire
water signal. Scanning **all 8,947 files** at the three darkest sites:

| site | files | flagged | rate |
|---|---:|---:|---:|
| THFR_H | 6,260 | 26 | 0.42% |
| BEFR_H | 1,798 | 4 | 0.22% |
| WRUK_H | 889 | 0 | 0.00% |
| **total** | **8,947** | **30** | **0.335%** |

The 30 flagged files, with their measured values, are listed in
{download}`simspec_flagged.csv <simspec_flagged.csv>`.

Characteristics of the 30:

- **All have `quality_flag = 0` and all have 6/6 valid scans.** Nothing in the
  flag bitmask fires.
- **Low sun is not the explanation.** Solar zenith spans 21.1–87.6°, median
  54.4°; only a third exceed 70°. (BEFR's four are all December 2023 at
  sza 85–87°, but THFR's 26 are spread across the whole year and all sun angles.)
- **Water darkness is the common factor.** For flagged spectra the *uncorrected*
  ρw(440–600) median is 0.00043, against site medians of 0.0106 (BEFR), 0.0097
  (THFR) and 0.0072 (WRUK) — roughly 20× darker than typical for the same site.
- **The offset subtracted exceeds the signal.** Median ρw removed is 0.00063
  against an uncorrected signal of 0.00043 — about **145% of the entire
  signal** — leaving a corrected median of −0.00019.
- They are spread across months and years (2023: 4, 2025: 9, 2026: 17), so this
  is a persistent behaviour rather than one bad episode.

This is the mirror image of the caveat already in the release notes: those warn
that `reflectance` is poor at the **turbid** sites (LPAR, MAFR, O1BE), where the
NIR signal is high. The behaviour here appears at the opposite extreme — the
**darkest** water, where the subtracted offset is comparable to or larger than
the water signal itself.

**This is known to the network and is being addressed.** The team is moving the
`reflectance` correction to a SWIR-based approach for all HYPSTAR sites
(WATERHYPERNET team, pers. comm., September 2026); the release notes likewise
list a SWIR-SimSpec upgrade of the `reflectance` product among the planned
evolutions. Note that the SWIR bands this relies on are not themselves
distributed in these files — no HYPSTAR grid in Release 2 extends past 1100 nm.

Two notes on scope. First, the rate is small (0.3%) and is a *lower* bound on
affected spectra, since it counts only cases where the visible median actually
changes sign; the shape-based clustering also picked up WRUK spectra that are
badly distorted without the median going negative. Second, `reflectance_nosc` is
positive and plausible in every one of the 30, so the underlying measurements
look sound — this concerns the correction, not the data acquisition.

Until the corrected product is reissued, a user working at the darkest sites is
well advised to apply the sign test above as a sanity check, or to work from
`reflectance_nosc`.

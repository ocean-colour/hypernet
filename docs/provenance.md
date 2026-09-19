# How these numbers were produced

All code is in the public [hypernet repository](https://github.com/ocean-colour/hypernet)
and runs end to end against a local copy of the archive, from the repository
root:

```bash
python -m hypernet.whn_explore 1       # index all 56,669 files from their names
python -m hypernet.whn_explore 2       # read 400/site, cluster, sample 100/site
python docs/whn_figures.py             # summary table + figures
python -m hypernet.whn_simspec_check   # over-subtraction scan at the dark sites
pytest -q                              # data-dependent tests self-skip
```

**Conventions.** `reflectance` everywhere except LPAR_H, MAFR_H and O1BE_P, which
use `reflectance_nosc`; Rrs = ρw/π and σ = σ(ρw)/π; the archive's own σ used as-is
with no floor; no wavelength trimming and no resampling of the data itself;
negative values retained.

**Sampling.** The ~100 spectra per site were drawn by optical diversity: 400
spectra per site spread uniformly through that site's record were read, their
L2-normalised 400–800 nm shapes clustered across all sites at once (k-means,
k=8, seed 1234), and each site's 100 apportioned across the clusters it occupies
by largest remainder — at least one per occupied cluster — spreading the picks
evenly in time within each cluster. L2 rather than area normalisation, because
retained negative values can drive an integral towards zero.

**Scope of each number.** File counts, date ranges, azimuth counts and the site
table come from all 56,669 files. The over-subtraction rates come from all 8,947
files at the three dark sites. Optical statistics (σ/Rrs, negative fractions, Rrs
magnitudes, clusters) come from the 4,400-spectrum pool — 400 per site. The
variable/attribute census covers 440 files across all 11 sites; the σ coverage
check covers 220; the `quality_flag` check covers 440. Figures are drawn on a
display grid (350–900 nm at 2.5 nm) used only for plotting and clustering, never
as a record grid.

# Getting Rrs right

## The products are ρw, not Rrs

Both `reflectance` and `reflectance_nosc` are **water-leaving reflectance ρw**
(dimensionless; `preferred_symbol = rhow`). This was confirmed numerically rather
than taken from the attribute: for both systems,

```
reflectance_nosc / (water_leaving_radiance / Ed) = 3.14159…  = π   (exactly)
```

so `water_leaving_radiance` is the *uncorrected* Lw, and

> **Rrs [sr⁻¹] = ρw / π**, and likewise **σ(Rrs) = σ(ρw) / π**.

The same ratio computed with the similarity-corrected `reflectance` gives 3.20
(HYPSTAR) and 2.96 (PANTHYR) — i.e. it differs from π exactly by the similarity
correction, which is applied to the reflectance but not to the stored radiance.

**All figures and tables in this guide are Rrs = ρw/π.**

## Which of the two products to use

Both systems ship two reflectance products:

| variable | meaning |
|---|---|
| `reflectance` | corrected with the **NIR** Similarity-Spectrum correction |
| `reflectance_nosc` | **not** similarity-corrected |

A point worth keeping straight, because the release notes describe two different
SWIR-related things: the **correction** applied to the `reflectance` product is
NIR-based for **both** systems. What changed in Release 2 is the separate
Similarity-Spectrum **quality-control test**, which is now SWIR-based for all
HYPSTAR sites, replacing an NIR-based QC that could not be applied adequately at
the most turbid sites. Correction and QC test are different things.

The release notes single out **LPAR, MAFR and O1BE** as sites where the corrected
`reflectance` is "known to be poor" — high NIR reflectance in turbid water — and
state `reflectance_nosc` is "definitely recommended" there. Everywhere else the
corrected `reflectance` is the standard product, and that is the convention used
throughout this guide.

If you work at the darkest sites, read
[the over-subtraction issue](over-subtraction.md)
as well: there is a case for sanity-checking the corrected product there too.

## The standard deviation is spelled differently per system

`std_reflectance` on HYPSTAR, `reflectance_std` on PANTHYR — and correspondingly
`std_reflectance_nosc` and `reflectance_nosc_std`. Any loader covering both
systems needs the mapping; see [Reading many files](reading-many-files.md) for the
rest of the schema differences.

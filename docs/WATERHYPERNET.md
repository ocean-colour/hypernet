# WATERHYPERNET Data Release 2 — an independent data user's guide

**J. Xavier Prochaska** (University of California, Santa Cruz) · 18 September 2026

*An independent guide to the WATERHYPERNET Release 2 archive, written from the
perspective of a data user preparing the dataset for ocean-colour algorithm
work. It is intended to complement — not replace — the official release notes.
A draft was reviewed by Kevin Ruddick, whose comments are incorporated
throughout; any remaining errors are mine.*

Every number here was measured from the archive itself or read from the release
notes. Scope of each measurement is stated in
[How these numbers were produced](#how-these-numbers-were-produced), and all of
the code is available (see the same section) so any figure can be reproduced or
extended.

**Where to start.** If you just want to read a spectrum, go to
[Quickstart](#quickstart-read-one-spectrum) and then
[Getting Rrs right](#getting-rrs-right) — between them they cover the two things
most likely to go silently wrong. If you are choosing which sites to work with,
start at [Choosing your data](#choosing-your-data).

---

## About the archive

WATERHYPERNET is a network of automated above-water hyperspectral radiometers
providing water-leaving reflectance for satellite validation and water-quality
monitoring. Release 2 (dated 2026-09-11) holds **56,669 measurements** —
**9.8 GB**, one NetCDF file per measurement — from **11 site–instrument
combinations** at 10 physical locations, spanning **2023-01-04 to 2026-09-09**.

Two instrument systems contribute:

| | HYPSTAR (`*_H`) | PANTHYR (`*_P`) |
|---|---|---|
| Sites | 7 | 4 |
| Spectra | 44,589 (79%) | 12,080 (21%) |
| Wavelengths | ~1,538 bands, 350–1100 nm, ~0.49 nm | 237 bands, 355–945 nm, 2.5 nm |
| Grid stability | **varies with instrument serial** | fixed |
| Processor | `hypernets_processor` 2.1.0 | `v20240912` |
| Licence attribute | CC BY-NC-ND | *(none)* |

Acqua Alta (VEIT) is the only location carrying both systems, which makes
`VEIT_H` vs `VEIT_P` a ready-made cross-system comparison.

### Getting access

The archive is not public. Request access through the contact form at
<https://waterhypernet.org/contact/>, which asks for your **name**, **e-mail**
and **affiliation** (all required), your country (optional), and a **purpose** —
choose *"Request for Data access"* — plus a free-text message and agreement to
the privacy policy.

Credentials for the password-protected distribution server,
<https://ftp.waterhypernet.org/>, follow by reply. The contact page itself does
not describe that step, so treat the form as the way in and expect the FTP
details to come back to you rather than to be published anywhere.

Before you use the data, read [Licence and citation](#licence-and-citation) —
the network asks for specific things in return, including PI priority use and,
in some cases, co-authorship.

### Scope: reflectance only, by design

The archive contains reflectance, viewing/solar geometry and QC — and nothing
else. There is no absorption, backscatter or attenuation, and no chlorophyll,
TSS/SPM, turbidity, salinity or water temperature.

This was checked exhaustively rather than assumed. The union of every variable
and global-attribute name over **440 files spanning all 11 sites** is perfectly
constant: **29 variables / 57 attributes** (HYPSTAR) and **22 variables / 58
attributes** (PANTHYR). No file anywhere carries a variable the others do not.
The only name matching any in-water-property pattern is `system_temperature`, an
instrument-housekeeping attribute, and it is `NaN` in every file examined.

**This is a deliberate design, not an omission.** As Kevin Ruddick puts it,
HYPERNETS is a single-parameter network — water reflectance — conceived as a
*core* that individual site PIs can expand according to their resources and
interests. Adding a common protocol for further measurands at every site is
logistically and financially out of reach for the network as a whole. Some PIs
do hold additional measurements, typically from short-duration deployments of
extra instruments; these are not part of this release, and the network is
considering how to advertise their existence to data users.

It is worth stating plainly what a reflectance-only archive is good for, because
it is a great deal:

1. **Satellite validation** — the network's primary purpose.
2. **Water-quality monitoring** for local users: time series of chlorophyll-a and
   phytoplankton type for water managers and aquatic biologists.
3. **Aquatic-optics research** into spectral variability itself.
4. **Cross-comparison with AERONET-OC**, co-located at VEIT, CBUS, TBBE and
   LPAR, which informs data quality and processing improvements for both.

What this scope *does* constrain is retrieval validation, discussed in
[Using these data for IOP retrieval](#using-these-data-for-iop-retrieval).

---

## Quickstart: read one spectrum

### Finding a file

```
RELEASE_2/
  0_README/WATERHYPERNET_ReleaseNotes_2-0.pdf      # the only non-NetCDF file
  <SITE>/<YYYY>/<MM>/<DD>/<one file per measurement>.nc
```

The two systems use different filename grammars. Note that the azimuth and
processing-time fields are **swapped** between them, and that acquisition time is
minute-resolution for HYPSTAR but second-resolution for PANTHYR:

```
HYPSTAR:  HYPERNETS_W_{SITE}_L2B_REF_{acqYYYYMMDDThhmm}_{procYYYYMMDDThhmm}_{RAA}_v2.1.nc
PANTHYR:  PANTHYR_W_{SITE}_L2A_REF_{acqYYYYMMDDThhmmss}_{AZ}_{procYYYYMMDDThhmmss}_v20240912_QA.nc
```

All 56,669 filenames parse cleanly under these two patterns.

### Reading it

With the `hypernet` package, which applies the per-site product choice and the
ρw → Rrs conversion for you:

```python
import glob, os
from hypernet.whn_explore import whn_root, load_spectrum

path = sorted(glob.glob(os.path.join(whn_root(), 'VEIT_H', '*', '*', '*', '*.nc')))[0]
spec = load_spectrum(path, 'VEIT_H')

spec['wave']          # native wavelength grid [nm]
spec['Rrs']           # remote-sensing reflectance [1/sr]
spec['sigma']         # its per-band standard deviation
spec['quality_flag']  # 0 = passed
```

Or with nothing but `netCDF4`:

```python
import netCDF4
import numpy as np

ds = netCDF4.Dataset(path)
wave = np.asarray(ds.variables['wavelength'][:], dtype=float).ravel()
rho_w = np.ma.filled(ds.variables['reflectance'][:].astype(float), np.nan).ravel()
ds.close()

Rrs = rho_w / np.pi          # the archive stores rho_w, NOT Rrs
```

Both give, for that file, 1,538 bands over 350.1–1099.9 nm and
Rrs(560) = 0.00938 sr⁻¹.

Two things the second version does not do for you, and which the rest of this
guide is largely about: it hard-codes `reflectance`, which is the wrong product
at three of the eleven sites, and it divides by π, which you must remember to do
to the standard deviation as well. Both are covered next.

---

## Getting Rrs right

### The products are ρw, not Rrs

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

### Which of the two products to use

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
[the over-subtraction issue](#known-issue-similarity-spectrum-over-subtraction-at-the-darkest-sites)
as well: there is a case for sanity-checking the corrected product there too.

### The standard deviation is spelled differently per system

`std_reflectance` on HYPSTAR, `reflectance_std` on PANTHYR — and correspondingly
`std_reflectance_nosc` and `reflectance_nosc_std`. Any loader covering both
systems needs the mapping; see [Reading many files](#reading-many-files) for the
rest of the schema differences.

---

## Reading many files

One spectrum per file in both systems (HYPSTAR `series=1`; PANTHYR
`sequence=1, time=1`). NETCDF4, no groups. Besides the two reflectance products
each file carries `water_leaving_radiance`, the downwelling/upwelling radiances
and irradiance, per-band standard deviations, the similarity-spectrum epsilons,
the air–water interface reflectance factor `rhof` (Mobley 1999) and its wind
input, scan counts, geometry and a quality flag.

### The two schemas differ

These are the harmonisation gaps most likely to trip up a data user writing one
loader for both systems:

| quantity | HYPSTAR | PANTHYR |
|---|---|---|
| dims | `wavelength, series` | `wavelength, sequence, time` |
| dtype | float32 | float64 |
| reflectance std | `std_reflectance` | `reflectance_std` |
| downwelling irradiance | `irradiance` | `downwelling_irradiance_mean` |
| downwelling radiance | `downwelling_radiance` | `downwelling_radiance_mean` |
| time | `acquisition_time`, uint32 epoch seconds, `units = "s"` | scalar float + an ISO string coordinate |
| instrument id | `system_id` (e.g. `HYPSTAR_122302`) | `l_sensor_sn` / `e_sensor_sn` (TriOS SAM serials) |

### The wavelength grid is not fixed

- **PANTHYR** is a fixed 237-band grid, 355–945 nm at exactly 2.5 nm, identical
  at all four sites and across sensor swaps (each PANTHYR site has used 2–3
  radiometer serials; the grid never changes).
- **HYPSTAR** grids change with the instrument serial number. Four sites contain
  two distinct grids each — GAIT, LPAR, MAFR and VEIT_H — differing in length
  (1536–1541 bands) and/or start wavelength (MAFR's two instruments both give
  1538 bands but start 0.18 nm apart). Spacing is ~0.463–0.498 nm; `bandwidth`
  is 3 nm FWHM throughout.

Spectra therefore **cannot be stacked without interpolation**, and a loader must
not assume a per-site grid — not even that the grid is the one the first file at
that site happened to use.

The release notes flag **<400 nm, >900 nm, and the ~762 nm O₂-A band** as
possibly unreliable and not recommended for satellite validation. Note that
despite the presence of `epsilon_SWIR` variables there is no SWIR coverage in the
distributed files — no HYPSTAR grid extends past 1100 nm.

### Time and identity

There is **no CF-decodable time**. `wavelength` has no `units` attribute in
either system; HYPSTAR `acquisition_time` has `units = "s"` (not
`seconds since …`) and PANTHYR's has none, so xarray will not decode either as a
datetime. The filename is the reliable source of acquisition time.

**Site + time is not a unique key.** WRUK has 345 pairs of files sharing an
acquisition minute and `sequence_id`, differing only in azimuth. More generally
the azimuth token is the relative azimuth for HYPSTAR (90 / 270, plus 225 at
WRUK) and a constant absolute pointing direction for PANTHYR (90 at CBUS, 270
elsewhere). LPAR is the only site with a substantially mixed azimuth population
(4,670×270 vs 1,312×090); WRUK splits 476×270 / 413×225.

Finally, the filename says `L2B` for HYPSTAR while `product_level` inside the
file says `W_L2A`. Trust neither over the other for anything load-bearing.

---

## Traps

The items below bite in code rather than in interpretation. The first is the one
that silently changes your results.

### PANTHYR's fill value collides with a real value

`_FillValue = 0` is declared on PANTHYR's `quality_flag`, its angles and
`bandwidth`. For the angles that is fine — they really are absent. For
`quality_flag` it is not: the variable is a bitmask whose zero value means *no
flags set*, i.e. the measurement **passed**, so an ordinary masked read turns
every passing PANTHYR measurement into a missing value. Filter on
`quality_flag == 0` after such a read and you discard all 12,080 PANTHYR spectra
while keeping every HYPSTAR one.

Read that one variable with masking disabled. HYPSTAR declares no `_FillValue`
here at all, so doing so changes nothing for it. `hypernet`'s `load_spectrum`
handles this for you.

`bandwidth` is entirely fill on PANTHYR.

### Integer variables cannot hold NaN

Cast to float **before** filling, or netCDF4 raises
`TypeError: Cannot convert fill_value nan to dtype uint8`.

### Negative reflectance is expected

Any log-space handling or area-normalisation must tolerate it. See
[Quality](#quality-what-to-filter-and-what-not-to) for why it is there and why
you should think twice before removing it.

### Two metadata oddities

- **Release window overrun.** The notes say data run to 2026-07-31, but O1BE_P
  and TBBE_P contain files to 2026-09-09.
- `instrument_calibration_file_rad` points at an `_IRR_` file, same as the
  irradiance attribute — apparently a metadata bug.

### Hazards for an automated loader — the checklist

1. **Two schemas.** Same variable names for the core products, different names
   for the std and irradiance variables, different dims and dtypes.
   → [Reading many files](#the-two-schemas-differ)
2. **The products are ρw, not Rrs.** Divide by π — and divide the std too.
   → [Getting Rrs right](#getting-rrs-right)
3. **HYPSTAR grids vary within a site**, by instrument serial.
   → [The wavelength grid is not fixed](#the-wavelength-grid-is-not-fixed)
4. **No CF-decodable time.** → [Time and identity](#time-and-identity)
5. **PANTHYR `_FillValue = 0`** masks a legitimate `quality_flag` of 0, and
   integer variables cannot hold NaN.
   → [above](#panthyrs-fill-value-collides-with-a-real-value)
6. **Negative reflectance is expected.**
   → [Quality](#quality-what-to-filter-and-what-not-to)
7. **Filename says `L2B` for HYPSTAR, but `product_level` inside says `W_L2A`.**
8. **WRUK timestamp collisions** — site+time is not a unique key.
9. **Release window overrun** — files past the documented end date.
10. `instrument_calibration_file_rad` points at an `_IRR_` file.

---

## Quality: what to filter and what not to

### The quality flag

`quality_flag` was **0 in all 440 files checked** — 40 per site, spread through
each site's record, read with masking off so fill and zero are distinguishable —
so the release appears to ship only measurements that already passed QC. HYPSTAR
defines a 30-bit flag (`lon_default`, `bad_pointing`, `outliers`, `L0_threshold`,
`dark_masked`, `not_enough_dark_scans`, `no_clear_sky_sequence`, `simil_fail`,
…); PANTHYR's is a uint8 whose `flag_meanings` are still placeholders
(`flag1 … flag8`).

### Negative reflectance is expected, and deliberate

The release notes give two reasons for retaining it: it is "metrologically
entirely valid to record a negative reflectance" when the true value is near zero
and zero lies within the uncertainty range, and excluding negatives biases the
statistics of any comparison. Frequency of spectra with at least one negative
band in 400–900 nm, over the pool sampled here:

| VEIT_H | GAIT_H | BEFR_H | THFR_H | WRUK_H | TBBE_P | VEIT_P | MAFR_H | CBUS_P | LPAR_H | O1BE_P |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 57% | 56% | 32% | 30% | 25% | 2.8% | 2.2% | 1.0% | 0% | 0% | 0% |

The three sites reading 0% are exactly the turbid ones read here via
`reflectance_nosc`, plus Chesapeake: bright water and no similarity correction
leave no near-zero bands to go negative. **Do not filter these out** without a
specific reason to.

### Uncertainty, and what it does not include

There is no uncertainty budget: `unc_comps` is empty on the reflectance and
radiance variables, and the release notes state uncertainties are "not yet
mature" and "currently not reported".

What *is* present is a per-band standard deviation (`std_reflectance` /
`reflectance_std` and their `_nosc` counterparts), and it is populated
everywhere — finite in **all 220 files sampled across all 11 sites**. Typical
magnitude, as median σ/|Rrs| over 450–650 nm:

| site | σ/Rrs | site | σ/Rrs |
|---|---:|---|---:|
| VEIT_P | 0.9% | VEIT_H | 2.0% |
| CBUS_P | 1.3% | LPAR_H | 2.2% |
| MAFR_H | 1.6% | O1BE_P | 2.4% |
| TBBE_P | 1.8% | GAIT_H | 2.5% |
| WRUK_H | 2.0% | THFR_H | 2.6% |
| BEFR_H | 2.5% | | |

**Caveat.** This is scan-to-scan variability within a measurement sequence only.
It excludes calibration uncertainty (only pre-deployment calibration was used)
and any error in the glint/sky-reflectance correction, so it is a lower bound on
the true uncertainty rather than an estimate of it.

### Geometry, and what you may want to screen on

Viewing zenith is ~40° throughout (HYPSTAR 39.8–40.3°, PANTHYR exactly 40.0°).
Solar zenith in the sampled pool spans 15.7–87.7°, i.e. the archive includes very
low-sun measurements that a user may wish to screen. Wind speed is present as the
GDAS-derived `rhof_wind` variable (HYPSTAR) or the `fresnel_wind` global
attribute (PANTHYR). PANTHYR applies a clear-sky test (Lsky/Ed at 750 nm < 0.05);
HYPSTAR does not, and may therefore pass lower-quality data.

Two gotchas: PANTHYR's `solar_azimuth_angle` is fill (0) in every file sampled —
use the `timestamp`/geometry attributes instead — and HYPSTAR's meteorological
attributes (`system_temperature`, pressure, humidity, illuminance) are all NaN.

---

## Choosing your data

### The sites

![Site locations](figs/fig_sites_map.png)

| site | location | system | spectra | first | last | months | lat | lon | PI |
|---|---|---|---:|---|---|---:|---:|---:|---|
| BEFR_H | Berre lagoon, FR | HYPSTAR | 1,798 | 2023-06-19 | 2024-01-04 | 7 | 43.44231 | 5.09718 | D. Doxaran (LOV) |
| CBUS_P | Chesapeake Bay, US | PANTHYR | 1,950 | 2023-07-11 | 2024-12-23 | 7 | 39.12 | −76.35 | — |
| GAIT_H | Lake Garda, IT | HYPSTAR | 6,373 | 2023-07-27 | 2026-07-31 | 9 | 45.57700 | 10.57942 | V. Brando (CNR) |
| LPAR_H | Río de la Plata, AR | HYPSTAR | 5,987 | 2024-03-24 | 2026-07-31 | 12 | −34.81797 | −57.89595 | A. Dogliotti (CONICET) |
| MAFR_H | Gironde / MAGEST, FR | HYPSTAR | 5,614 | 2024-06-21 | 2026-07-31 | 12 | 45.54765 | −1.04050 | D. Doxaran (LOV) |
| O1BE_P | Oostende RT1, BE | PANTHYR | 4,607 | 2023-04-25 | 2026-09-09 | 12 | 51.25 | 2.92 | — |
| TBBE_P | Thornton Bank, BE | PANTHYR | 2,021 | 2023-05-11 | 2026-09-06 | 12 | 51.53 | 2.96 | — |
| THFR_H | Thau lagoon, FR | HYPSTAR | 6,260 | 2025-05-28 | 2026-07-29 | 12 | 43.43486 | 3.66416 | D. Doxaran (LOV) |
| VEIT_H | Acqua Alta AAOT, IT | HYPSTAR | 17,668 | 2023-04-24 | 2026-07-31 | 12 | 45.31420 | 12.50830 | V. Brando (CNR) |
| VEIT_P | Acqua Alta AAOT, IT | PANTHYR | 3,502 | 2023-01-04 | 2026-05-31 | 12 | 45.31 | 12.51 | — |
| WRUK_H | Wraysbury reservoir, UK | HYPSTAR | 889 | 2024-03-05 | 2024-09-15 | 7 | 51.45758 | −0.53219 | A. Bialek (NPL) |

The same table in machine-readable form: [`summary_table.csv`](summary_table.csv),
which carries additional per-site columns (band counts, solar-zenith range,
median Rrs(560), blue–green ratio, relative σ, negative fraction).

HYPSTAR coordinates are fixed per-site attributes. PANTHYR files instead carry
per-file GPS averages that jitter at the ~1e-5° level, so the values above are
medians. The coordinate is the radiometer's position; the water target is 3–20 m
away. GAIT and WRUK coordinates were wrong in Release 1 and are corrected here.

### Coverage is very uneven

![Data volume and temporal coverage](figs/fig_data_volume.png)

VEIT_H alone is 31% of the archive; WRUK_H is 1.6%. More importantly for anyone
designing a sampling scheme, **only 6 of the 11 sites span a full seasonal
cycle**: BEFR (2023-06 → 2024-01), WRUK (2024 only), CBUS (7 distinct months)
and GAIT (9 months) do not. Several sites have multi-month gaps mid-record.

### What the spectra look like

![Per-site median spectra, common y axis](figs/fig_median_spectra.png)

On a common y axis the range of the network is immediate: the Río de la Plata,
Gironde and Oostende are bright, while Berre, Thau, Wraysbury and — perhaps
surprisingly — Acqua Alta and Lake Garda are dark enough to be nearly flat at
this scale. Median Rrs(560) runs from 0.0030 sr⁻¹ at Wraysbury to 0.036 sr⁻¹ at
the Gironde, a factor of 12, and the blue–green ratio Rrs(490)/Rrs(560) from 0.52
(Berre) to 1.03 (Acqua Alta) — from strongly green/turbid to near-marine.

The same data with each panel autoscaled recovers the spectral shape at the dark
sites, which the common axis necessarily suppresses:

![Per-site median spectra, per-panel autoscale](figs/fig_median_spectra_autoscale.png)

![Spectral-shape clustering](figs/fig_clusters.png)

Clustering the L2-normalised 400–800 nm shapes across all sites (k=8) shows the
sites occupy genuinely different optical regimes rather than one continuum: LPAR
is almost entirely one cluster, CBUS and O1BE another, while VEIT and GAIT split
across several. Cluster 5 (n=9) is the group discussed in
[the over-subtraction issue](#known-issue-similarity-spectrum-over-subtraction-at-the-darkest-sites).

![Band time series](figs/fig_band_timeseries.png)

Seasonal cycles are visible at GAIT and VEIT; LPAR and MAFR are persistently
turbid with little seasonal structure.

---

## Known issue: Similarity-Spectrum over-subtraction at the darkest sites

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
[`simspec_flagged.csv`](simspec_flagged.csv).

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

---

## Using these data for IOP retrieval

A reflectance-only archive can exercise inherent-optical-property (IOP) retrieval
algorithms but cannot, on its own, validate what they retrieve. An algorithm can
be run on these spectra and the Rrs residual scored; algorithms can be compared
against one another for consistency and for the spread of their retrieved IOPs.
But nothing in the release constrains retrieved `a_ph`, `a_dg`, `bb_p`,
chlorophyll or the spectral slopes, so no truth-referenced error metric is
available. That is a statement about the intended scope of the network, not a
defect.

Three things nonetheless make the archive valuable for this work:

- **Optical diversity with a stable instrument.** Eleven sites spanning a factor
  of 12 in brightness and 0.52–1.03 in blue–green ratio, measured with two
  well-characterised systems on a fixed geometry, is a good stress test for
  algorithm robustness — particularly the very turbid sites, where many
  ocean-colour IOP algorithms degrade.
- **Hyperspectral resolution.** The HYPSTAR grid (~0.49 nm) is far finer than
  the satellite bands most IOP algorithms were built for, so the same spectrum
  can be convolved to PACE, MODIS, SeaWiFS or SBG bands and the retrieval
  compared across spectral configurations.
- **Repeat measurements at fixed stations.** Thousands of spectra from one
  location allow the *stability* of a retrieval to be assessed — how much of the
  retrieved variance is water and how much is algorithm noise — which
  single-visit field campaigns cannot support.

One such effort is [IOPtics](https://github.com/ocean-colour/IOPtics)
([documentation](https://ioptics.readthedocs.io/en/stable/)), a framework for
testing and comparing IOP retrieval algorithms, which is what motivated this
exploration.

Looking ahead, the measurands Kevin Ruddick would most like to see added to the
sites — automated flow cytometry (phytoplankton composition), turbidimeters and
chlorophyll fluorimeters — would be **more** directly useful for validating
in-water retrievals than IOPs would, and he notes that IOPs are usually an
intermediate, explanatory parameter rather than an end-user product. Chlorophyll
and turbidity/TSS also map onto truth quantities that retrieval-evaluation
frameworks already handle, so such data could be folded in without new machinery.

---

## Licence and citation

HYPSTAR files carry `licence = "Attribution-NonCommercial-NoDerivs CC BY-NC-ND"`
and a long `acknowledgement` attribute asking users to respect PI priority use,
cite the key HYPERNETS papers, and offer PI co-authorship where the data are a
principal component of a publication. **PANTHYR files carry no licence
attribute** — one of the harmonisation gaps noted above. The release notes carry
no explicit licence clause; they refer to the WATERHYPERNET data policy and ask
users to cite Ruddick et al. (2024) and De Vis et al. (2024). Data are
distributed via the password-protected `https://ftp.waterhypernet.org/`.

Key references from the release notes:

1. Ruddick, K. G., *et al.* (2024), "WATERHYPERNET: a prototype network of
   automated in situ measurements of hyperspectral water reflectance for
   satellite validation and water quality monitoring", *Front. Remote Sens.* 5.
2. De Vis, P., *et al.* (2024), "Generating Hyperspectral Reference Measurements
   for Surface Reflectance from the LANDHYPERNET and WATERHYPERNET Networks",
   *Front. Remote Sens.* 5, 1347230.
3. Kuusk, J., *et al.* (2024), "HYPSTAR: a Hyperspectral Pointable System for
   Terrestrial and Aquatic Radiometry", *Front. Remote Sens.*
4. Vansteenwegen, D., *et al.* (2019), "The Pan-and-Tilt Hyperspectral
   Radiometer System (PANTHYR)…", *Remote Sensing* 11, 1360.
5. Mobley, C. D. (1999), "Estimation of the remote-sensing reflectance from
   above-surface measurements", *Appl. Opt.* 38(36), 7442–7455.
6. Ruddick, K., *et al.* (2006), "Seaborne measurements of near infrared
   water-leaving reflectance: the similarity spectrum for turbid waters",
   *Limnol. Oceanogr.* 51(2), 1167–1179.

---

## How these numbers were produced

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

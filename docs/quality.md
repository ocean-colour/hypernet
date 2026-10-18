# Quality: what to filter and what not to

## The quality flag

`quality_flag` was **0 in all 440 files checked** — 40 per site, spread through
each site's record, read with masking off so fill and zero are distinguishable —
so the release appears to ship only measurements that already passed QC. HYPSTAR
defines a 30-bit flag (`lon_default`, `bad_pointing`, `outliers`, `L0_threshold`,
`dark_masked`, `not_enough_dark_scans`, `no_clear_sky_sequence`, `simil_fail`,
…); PANTHYR's is a uint8 whose `flag_meanings` are still placeholders
(`flag1 … flag8`).

## Negative reflectance is expected, and deliberate

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

## Uncertainty, and what it does not include

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

## Geometry, and what you may want to screen on

Viewing zenith is ~40° throughout (HYPSTAR 39.8–40.3°, PANTHYR exactly 40.0°).
Solar zenith in the sampled pool spans 15.7–87.7°, i.e. the archive includes very
low-sun measurements that a user may wish to screen. Wind speed is present as the
GDAS-derived `rhof_wind` variable (HYPSTAR) or the `fresnel_wind` global
attribute (PANTHYR). PANTHYR applies a clear-sky test (Lsky/Ed at 750 nm < 0.05);
HYPSTAR does not, and may therefore pass lower-quality data.

Two gotchas: PANTHYR's `solar_azimuth_angle` is fill (0) in every file sampled —
use the `timestamp`/geometry attributes instead — and HYPSTAR's meteorological
attributes (`system_temperature`, pressure, humidity, illuminance) are all NaN.

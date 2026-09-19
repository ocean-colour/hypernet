# Reading many files

One spectrum per file in both systems (HYPSTAR `series=1`; PANTHYR
`sequence=1, time=1`). NETCDF4, no groups. Besides the two reflectance products
each file carries `water_leaving_radiance`, the downwelling/upwelling radiances
and irradiance, per-band standard deviations, the similarity-spectrum epsilons,
the air–water interface reflectance factor `rhof` (Mobley 1999) and its wind
input, scan counts, geometry and a quality flag.

## The two schemas differ

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

## The wavelength grid is not fixed

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

## Time and identity

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

# Traps

The items below bite in code rather than in interpretation. The first is the one
that silently changes your results.

## PANTHYR's fill value collides with a real value

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

## Integer variables cannot hold NaN

Cast to float **before** filling, or netCDF4 raises
`TypeError: Cannot convert fill_value nan to dtype uint8`.

## Negative reflectance is expected

Any log-space handling or area-normalisation must tolerate it. See
[Quality](quality.md) for why it is there and why
you should think twice before removing it.

## Two metadata oddities

- **Release window overrun.** The notes say data run to 2026-07-31, but O1BE_P
  and TBBE_P contain files to 2026-09-09.
- `instrument_calibration_file_rad` points at an `_IRR_` file, same as the
  irradiance attribute — apparently a metadata bug.

## Hazards for an automated loader — the checklist

1. **Two schemas.** Same variable names for the core products, different names
   for the std and irradiance variables, different dims and dtypes.
   → [Reading many files](reading-many-files.md#the-two-schemas-differ)
2. **The products are ρw, not Rrs.** Divide by π — and divide the std too.
   → [Getting Rrs right](rrs.md)
3. **HYPSTAR grids vary within a site**, by instrument serial.
   → [The wavelength grid is not fixed](reading-many-files.md#the-wavelength-grid-is-not-fixed)
4. **No CF-decodable time.** → [Time and identity](reading-many-files.md#time-and-identity)
5. **PANTHYR `_FillValue = 0`** masks a legitimate `quality_flag` of 0, and
   integer variables cannot hold NaN.
   → [above](#panthyrs-fill-value-collides-with-a-real-value)
6. **Negative reflectance is expected.**
   → [Quality](quality.md)
7. **Filename says `L2B` for HYPSTAR, but `product_level` inside says `W_L2A`.**
8. **WRUK timestamp collisions** — site+time is not a unique key.
9. **Release window overrun** — files past the documented end date.
10. `instrument_calibration_file_rad` points at an `_IRR_` file.

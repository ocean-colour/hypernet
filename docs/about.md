# About the archive

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

## Getting access

The archive is not public. Request access through the contact form at
<https://waterhypernet.org/contact/>, which asks for your **name**, **e-mail**
and **affiliation** (all required), your country (optional), and a **purpose** —
choose *"Request for Data access"* — plus a free-text message and agreement to
the privacy policy.

Credentials for the password-protected distribution server,
<https://ftp.waterhypernet.org/>, follow by reply. The contact page itself does
not describe that step, so treat the form as the way in and expect the FTP
details to come back to you rather than to be published anywhere.

Before you use the data, read [Licence and citation](licence.md) —
the network asks for specific things in return, including PI priority use and,
in some cases, co-authorship.

## Scope: reflectance only, by design

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
[Using these data for IOP retrieval](iop-retrieval.md).

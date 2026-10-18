# Quickstart: read one spectrum

## Finding a file

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

## Reading it

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

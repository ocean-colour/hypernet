"""Which L1A_RAD scans are sky (Ld) and which are water (Lu)?

line_fits_veit.py assigns Ld = viewing_zenith_angle < 90.  L1A_IRR has
vza = 180 (looking up), which suggests the opposite convention (vza measured
from nadir, so 40 deg = water, 140 deg = sky).  Settle it by comparing the
L1A_RAD means of each vza group with L1C `downwelling_radiance` and
`upwelling_radiance` (all on the same L grid).
"""
import glob
import os
import numpy as np
import xarray as xr

d = os.path.join(os.environ['OS_COLOR'], 'WATERHYPERNET', 'Wavelengths', '')
irr = xr.open_dataset(glob.glob(d + '*L1A_IRR*.nc')[0])
rad = xr.open_dataset(glob.glob(d + '*L1A_RAD*.nc')[0])
l1c = xr.open_dataset(glob.glob(d + '*L1C_ALL*.nc')[0])

print('L1A_IRR vza:', irr.viewing_zenith_angle.values)
vza = rad.viewing_zenith_angle.values
print('L1A_RAD vza:', vza)
R = rad.radiance.values
lo = np.nanmean(R[:, vza < 90], axis=1)
hi = np.nanmean(R[:, vza >= 90], axis=1)
Ldc = np.nanmean(l1c.downwelling_radiance.values, axis=1)
Luc = np.nanmean(l1c.upwelling_radiance.values, axis=1)
w = rad.wavelength.values
band = (w > 400) & (w < 700)


def medrel(a, b):
    return np.nanmedian(np.abs(a[band] - b[band]) / np.abs(b[band]))


print('median |rel diff| 400-700 nm:')
print('  vza<90  vs L1C downwelling_radiance: %.3e' % medrel(lo, Ldc))
print('  vza<90  vs L1C upwelling_radiance:   %.3e' % medrel(lo, Luc))
print('  vza>=90 vs L1C downwelling_radiance: %.3e' % medrel(hi, Ldc))
print('  vza>=90 vs L1C upwelling_radiance:   %.3e' % medrel(hi, Luc))
for lam in (450, 550, 650):
    j = np.argmin(np.abs(w - lam))
    print('  %d nm: mean(vza<90) = %.4g, mean(vza>=90) = %.4g' % (lam, lo[j], hi[j]))
for v in ('downwelling_radiance', 'upwelling_radiance'):
    print(v, dict(l1c[v].attrs))

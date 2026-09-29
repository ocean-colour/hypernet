"""How do the L1C quantities of the VEIT sample relate to L1A?

Checks, per scan and for the means:

1. L1C `irradiance` vs `np.interp` of L1A_IRR (mean and per scan).
2. L1C `upwelling_radiance` vs the L1A_RAD water-view scans (vza < 90).
3. L1C `downwelling_radiance` vs the L1A_RAD sky scans (vza >= 90), mean
   and per scan.
4. Per-scan metadata: series_id, quality_flag, acquisition_time, SZA.
5. L1C `u_rel_random_*` vs the relative scan-to-scan scatter of L1A.

Run from the repository root: ``python wiggles/phase0a_l1c_consistency.py``.
"""
import glob
import os

import numpy as np
import xarray as xr

d = os.path.join(os.environ['OS_COLOR'], 'WATERHYPERNET', 'Wavelengths', '')
irr = xr.open_dataset(glob.glob(d + '*L1A_IRR*.nc')[0])
rad = xr.open_dataset(glob.glob(d + '*L1A_RAD*.nc')[0])
l1c = xr.open_dataset(glob.glob(d + '*L1C_ALL*.nc')[0])

wE = irr.wavelength.values
wL = rad.wavelength.values
band = (wL > 400) & (wL < 900)


def rel(a, b):
    """median and max |a-b|/|b| over 400-900 nm."""
    r = np.abs(a[band] - b[band]) / np.abs(b[band])
    return np.nanmedian(r), np.nanmax(r)


for name, ds in [('L1A_IRR', irr), ('L1A_RAD', rad), ('L1C', l1c)]:
    print('==', name)
    for v in ('series_id', 'quality_flag', 'acquisition_time',
              'viewing_zenith_angle', 'viewing_azimuth_angle',
              'solar_zenith_angle'):
        print('   %-22s %s' % (v, ds[v].values))
print('   bandwidth unique (IRR, RAD):', np.unique(irr.bandwidth.values),
      np.unique(rad.bandwidth.values))

# 1. irradiance
E = irr.irradiance.values
Ec = l1c.irradiance.values
print('\n[1] L1C irradiance vs np.interp(L1A mean): median %.2e max %.2e'
      % rel(Ec.mean(axis=1), np.interp(wL, wE, E.mean(axis=1))))
for k in range(Ec.shape[1]):
    print('    L1C scan %d vs interp(L1A mean) %.2e / %.2e; vs interp(L1A scan %d) %.2e / %.2e'
          % ((k,) + rel(Ec[:, k], np.interp(wL, wE, E.mean(axis=1))) + (k,)
             + rel(Ec[:, k], np.interp(wL, wE, E[:, k]))))

# 2./3. radiances
vza = rad.viewing_zenith_angle.values
R = rad.radiance.values
Lu_s, Ld_s = R[:, vza < 90], R[:, vza >= 90]
Luc, Ldc = l1c.upwelling_radiance.values, l1c.downwelling_radiance.values
print('\n[2] L1C Lu scan k vs L1A water scan k:')
for k in range(Luc.shape[1]):
    print('    %d: %.2e / %.2e' % ((k,) + rel(Luc[:, k], Lu_s[:, k])))
print('\n[3] L1C Ld: scans identical to each other? max rel spread %.2e'
      % np.nanmax(np.abs(Ldc[band] - Ldc[band, :1]) / Ldc[band, :1]))
print('    L1C Ld scan 0 vs L1A sky mean (6): %.2e / %.2e' % rel(Ldc[:, 0], Ld_s.mean(axis=1)))
print('    vs mean of first 3 sky: %.2e / %.2e' % rel(Ldc[:, 0], Ld_s[:, :3].mean(axis=1)))
print('    vs mean of last 3 sky:  %.2e / %.2e' % rel(Ldc[:, 0], Ld_s[:, 3:].mean(axis=1)))
print('    vs median of 6 sky:     %.2e / %.2e' % rel(Ldc[:, 0], np.median(Ld_s, axis=1)))
for k in range(Ld_s.shape[1]):
    print('    vs sky scan %d:          %.2e / %.2e' % ((k,) + rel(Ldc[:, 0], Ld_s[:, k])))
print('    L1C std_downwelling_radiance / L1C Ld (median 400-900): %.2e'
      % np.nanmedian(l1c.std_downwelling_radiance.values[band] / Ldc[band, 0]))

# 5. uncertainties
print('\n[5] relative uncertainty, median over 400-900 nm:')
for lab, scans, u in [('E', E, 'u_rel_random_irradiance'),
                      ('Lu', Lu_s, 'u_rel_random_upwelling_radiance'),
                      ('Ld', Ld_s, 'u_rel_random_downwelling_radiance')]:
    s = np.nanstd(scans, axis=1, ddof=1) / np.nanmean(scans, axis=1)
    if lab == 'E':
        s = np.interp(wL, wE, s)
    urr = l1c[u].values
    print('    %-2s scan-to-scan std/mean %.2e   L1C %s (units %s) median %.2e'
          % (lab, np.nanmedian(s[band]), u, l1c[u].attrs.get('units', '?'),
             np.nanmedian(urr[band])))

# 6. Is L1C Ld the sky mean interpolated in time to the water-view time?
t = rad.acquisition_time.values.astype(float)
t_sky = np.unique(t[vza >= 90])
t_lu = np.unique(t[vza < 90])[0]
m1 = Ld_s[:, :3].mean(axis=1)
m2 = Ld_s[:, 3:].mean(axis=1)
w1 = (t_sky[1] - t_lu) / (t_sky[1] - t_sky[0])
print('\n[6] sky series at t = %s, water at t = %d -> weight on first sky series %.3f'
      % (t_sky, t_lu, w1))
print('    L1C Ld vs time-interpolated sky mean: %.2e / %.2e'
      % rel(Ldc[:, 0], w1 * m1 + (1 - w1) * m2))
tE = np.unique(irr.acquisition_time.values.astype(float))
wE1 = (tE[1] - t_lu) / (tE[1] - tE[0])
print('    E series at t = %s -> weight on first E series %.3f (0.5 = plain mean)'
      % (tE, wE1))

# 7. How much of the scan-to-scan scatter is broadband (a smooth scale factor
#    per scan) rather than per-pixel noise?  Divide each scan by a smooth fit of
#    its ratio to the mean, then recompute std/mean.
from scipy.ndimage import gaussian_filter1d  # noqa: E402

print('\n[7] relative scan-to-scan scatter, median 400-900 nm: raw vs after removing'
      ' a smooth (20-px Gaussian) scale per scan')
for lab, scans, w in [('E', E, wE), ('Lu', Lu_s, wL), ('Ld (all 6)', Ld_s, wL),
                      ('Ld (1st 3)', Ld_s[:, :3], wL), ('Ld (2nd 3)', Ld_s[:, 3:], wL)]:
    mean = np.nanmean(scans, axis=1)
    raw = np.nanstd(scans, axis=1, ddof=1) / mean
    ratio = scans / mean[:, None]
    smooth = gaussian_filter1d(ratio, 20, axis=0, mode='nearest')
    flat = scans / smooth
    fl = np.nanstd(flat, axis=1, ddof=1) / np.nanmean(flat, axis=1)
    b = (w > 400) & (w < 900)
    print('    %-11s raw %.2e   flattened %.2e   (x%.0f)'
          % (lab, np.nanmedian(raw[b]), np.nanmedian(fl[b]),
             np.nanmedian(raw[b]) / np.nanmedian(fl[b])))

# 8. Where does L1C irradiance depart from np.interp of the L1A mean?
Elin = np.interp(wL, wE, E.mean(axis=1))
r = np.abs(Ec[:, 0] - Elin) / Elin
print('\n[8] |L1C E - interp(L1A mean E)| / interp, max by range:')
for lo, hi in [(350, 380), (380, 400), (400, 900), (900, 1000), (1000, 1100)]:
    m = (wL >= lo) & (wL < hi)
    print('    %4d-%4d nm: max %.2e (at %.1f nm), median %.2e, median E %.3g'
          % (lo, hi, np.nanmax(r[m]), wL[m][np.nanargmax(r[m])], np.nanmedian(r[m]),
             np.nanmedian(Elin[m])))

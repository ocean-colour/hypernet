"""Quick read-only sanity checks on the VEIT sample sequence.

1. Is L1C irradiance exactly the linear interpolant of mean L1A_IRR onto the L grid?
2. How big is the linear-vs-cubic interpolation difference of Ed (a proxy for the
   linear interpolation error)?  Where is it largest?
3. Is Ld/Ed (sky radiance / irradiance) smooth?  Scan a rigid wavelength shift of
   the E grid and see which shift minimises the high-frequency power of Ld/Ed.
4. Does rho_w'' correlate with Ed''?
"""
import glob
import os
import numpy as np
import xarray as xr
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter1d

d = os.path.join(os.environ['OS_COLOR'], 'WATERHYPERNET', 'Wavelengths', '')
irr = xr.open_dataset(glob.glob(d + '*L1A_IRR*.nc')[0])
rad = xr.open_dataset(glob.glob(d + '*L1A_RAD*.nc')[0])
l1c = xr.open_dataset(glob.glob(d + '*L1C_ALL*.nc')[0])
l2a = xr.open_dataset(glob.glob(d + '*L2A_REF*.nc')[0])

wE = irr.wavelength.values
wL = rad.wavelength.values
Ed_scans = irr.irradiance.values  # (wavelength, scan)
print('irr shape', Ed_scans.shape, 'rad shape', rad.radiance.values.shape)
Ed = np.nanmean(Ed_scans, axis=1)

# ---- 1. L1C irradiance vs linear interpolant of L1A mean
Ed_lin = np.interp(wL, wE, Ed)
Ed_l1c = np.nanmean(l1c.irradiance.values, axis=1)
ok = np.isfinite(Ed_l1c) & (Ed_l1c > 0)
rel = (Ed_l1c[ok] - Ed_lin[ok]) / Ed_l1c[ok]
print('\n[1] L1C irradiance vs np.interp(L1A mean): median |rel| = %.2e, max |rel| = %.2e'
      % (np.median(np.abs(rel)), np.max(np.abs(rel))))
print('    (L1C is a 6-scan mean; L1A_IRR also 6 scans -> should match if same scans kept)')

# ---- 2. linear vs cubic-spline interpolation of Ed
cs = CubicSpline(wE, Ed)
Ed_cub = cs(wL)
inside = (wL > wE[0]) & (wL < wE[-1]) & (Ed_lin > 0)
diff = (Ed_lin - Ed_cub) / Ed_cub
print('\n[2] (linear - cubic)/cubic for Ed onto the L grid:')
for lo, hi in [(380, 400), (400, 500), (500, 600), (600, 700), (700, 800), (800, 950), (950, 1050)]:
    m = inside & (wL >= lo) & (wL < hi)
    print('    %4d-%4d nm: rms %.2e  max|.| %.2e at %.1f nm' % (
        lo, hi, np.sqrt(np.mean(diff[m] ** 2)), np.max(np.abs(diff[m])),
        wL[m][np.argmax(np.abs(diff[m]))]))
# the classic lines
for name, lam in [('Ca K 393', 393.4), ('Ca H 397', 396.8), ('G 430', 430.8), ('Hb 486', 486.1),
                  ('Mg b 517', 517.5), ('Na D 589', 589.3), ('Ha 656', 656.3), ('O2 B 687', 687.0),
                  ('O2 A 760', 760.5), ('H2O 940', 940.0)]:
    j = np.argmin(np.abs(wL - lam))
    sl = slice(max(j - 6, 0), j + 7)
    print('    %-10s max|lin-cub|/cub within +-3 nm: %.2e' % (name, np.max(np.abs(diff[sl]))))

# offset between the two grids (in pixels)
j = np.searchsorted(wE, wL)
j = np.clip(j, 1, wE.size - 1)
frac = (wL - wE[j - 1]) / (wE[j] - wE[j - 1])
print('    fractional-pixel position of L pixels on the E grid: min %.2f median %.2f max %.2f'
      % (frac[inside].min(), np.median(frac[inside]), frac[inside].max()))

# ---- 3. Ld / Ed smoothness and rigid wavelength-shift scan
Ld = np.nanmean(l1c.downwelling_radiance.values, axis=1)  # on L grid
Lu = np.nanmean(l1c.upwelling_radiance.values, axis=1)


def hf_power(y, sigma_px=4):
    """rms of (y - smooth(y))/smooth(y): high-frequency wiggle power."""
    s = gaussian_filter1d(y, sigma_px)
    return np.sqrt(np.nanmean(((y - s) / s) ** 2))


band = (wL > 400) & (wL < 700)
print('\n[3] Ld/Ed high-frequency wiggle power (400-700 nm) vs rigid shift of the E grid:')
best = None
for dl in np.arange(-1.0, 1.01, 0.1):
    Ed_s = cs(wL - dl)   # cubic interpolant of Ed evaluated on a shifted grid
    r = Ld[band] / Ed_s[band]
    p = hf_power(r)
    if best is None or p < best[1]:
        best = (dl, p)
    print('    shift %+.1f nm: %.3e' % (dl, p))
print('    best rigid shift: %+.1f nm (wiggle power %.3e)' % best)
print('    same metric with LINEAR interp, shift 0: %.3e' % hf_power(Ld[band] / Ed_lin[band]))
print('    intrinsic wiggle power of Ed itself (same metric): %.3e' % hf_power(Ed_cub[band]))
print('    intrinsic wiggle power of Ld itself: %.3e' % hf_power(Ld[band]))

# ---- 4. rho_w'' vs Ed''
rho = l2a.reflectance.values.ravel()
rho_nosc = l2a.reflectance_nosc.values.ravel()
h = np.gradient(wL)


def d2(y):
    return np.gradient(np.gradient(y, wL), wL)


m = band & np.isfinite(rho) & (rho > 0)
r2 = d2(rho)
e2 = d2(Ed_lin) / Ed_lin
l2 = d2(Lu) / Lu
print('\n[4] correlation of rho_w\'\' with (Ed\'\'/Ed) and (Lu\'\'/Lu) over 400-700 nm:')
print('    corr(rho\'\', Ed\'\'/Ed) = %.2f' % np.corrcoef(r2[m], e2[m])[0, 1])
print('    corr(rho\'\', Lu\'\'/Lu) = %.2f' % np.corrcoef(r2[m], l2[m])[0, 1])
print('    rms rho\'\' = %.2e, rms(rho * Ed\'\'/Ed) = %.2e, rms(rho*Lu\'\'/Lu) = %.2e'
      % (np.sqrt(np.mean(r2[m] ** 2)), np.sqrt(np.mean((rho[m] * e2[m]) ** 2)),
         np.sqrt(np.mean((rho[m] * l2[m]) ** 2))))
# relative wiggle amplitude of rho_w itself
print('    rho_w hf wiggle power (400-700): %.3e; rho_nosc: %.3e'
      % (hf_power(rho[m]), hf_power(rho_nosc[m])))
print('    rho_w median 400-700: %.4f' % np.median(rho[m]))

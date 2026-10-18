"""Fit Gaussian absorption lines in Ed (E grid) and Ld / Lu (L grid) of the VEIT
sample sequence.  Centroid differences -> relative wavelength calibration of E vs
L; width differences -> whether the E and L spectral response functions differ.
Line depths are confounded by the Ring effect in Ld, widths much less so.
"""
import glob
import os
import numpy as np
import xarray as xr
from scipy.optimize import curve_fit

d = os.path.join(os.environ['OS_COLOR'], 'WATERHYPERNET', 'Wavelengths', '')
irr = xr.open_dataset(glob.glob(d + '*L1A_IRR*.nc')[0])
rad = xr.open_dataset(glob.glob(d + '*L1A_RAD*.nc')[0])
wE = irr.wavelength.values
wL = rad.wavelength.values
Ed = np.nanmean(irr.irradiance.values, axis=1)
vza = rad.viewing_zenith_angle.values
print('RAD scan vza:', vza)
R = rad.radiance.values
# NB: labels are swapped -- vza is from nadir, so vza < 90 is the water view (Lu)
# and vza >= 90 the sky (Ld).  See wavecal/check_vza_convention_veit.py.
Ld = np.nanmean(R[:, vza < 90], axis=1)
Lu = np.nanmean(R[:, vza >= 90], axis=1)
print('n Ld scans', (vza < 90).sum(), 'n Lu scans', (vza >= 90).sum())


def model(x, c0, c1, a, mu, sig):
    return (c0 + c1 * (x - mu)) * (1 - a * np.exp(-0.5 * ((x - mu) / sig) ** 2))


def fit(w, y, lam, half=4.0):
    m = (w > lam - half) & (w < lam + half) & np.isfinite(y)
    x, yy = w[m], y[m]
    p0 = [np.median(yy), 0.0, 0.2, lam, 1.3]
    try:
        p, cov = curve_fit(model, x, yy, p0=p0, maxfev=20000)
        e = np.sqrt(np.diag(cov))
        return p, e
    except Exception as ex:
        return None, None


lines = [('Ca K', 393.37), ('Ca H', 396.85), ('G band', 430.8), ('H beta', 486.13),
         ('Mg b', 517.3), ('Na D', 589.3), ('H alpha', 656.28), ('O2 B', 687.0),
         ('O2 A', 760.6), ('H2O', 936.5)]
print('\n%-8s %22s %22s %22s   dmu(L-E)  FWHM_E FWHM_Ld FWHM_Lu' % ('line', 'Ed mu,sig,depth', 'Ld mu,sig,depth', 'Lu mu,sig,depth'))
for name, lam in lines:
    pe, ee = fit(wE, Ed, lam)
    pl, el = fit(wL, Ld, lam)
    pu, eu = fit(wL, Lu, lam)
    if pe is None or pl is None:
        print(name, 'fit failed')
        continue
    fw = 2.3548
    s_u = '%7.2f %5.2f %5.3f' % (pu[3], abs(pu[4]), pu[2]) if pu is not None else ' ' * 20
    print('%-8s %7.2f %5.2f %5.3f   %7.2f %5.2f %5.3f   %s   %+6.3f    %5.2f  %5.2f  %5s' % (
        name, pe[3], abs(pe[4]), pe[2], pl[3], abs(pl[4]), pl[2], s_u,
        pl[3] - pe[3], fw * abs(pe[4]), fw * abs(pl[4]),
        '%5.2f' % (fw * abs(pu[4])) if pu is not None else '-'))

"""Can native-grid E be recovered from the irradiance stored on the L grid?

The processor builds L1C/L2 ``irradiance`` as ``np.interp(wL, wE, E)`` -- a
linear map E_L = M @ E with M (n_L x n_E) fixed by the two wavelength grids.
With the lab wavelength polynomials (``wavecal/lab_data_examine.py``) every
instrument's E grid is known, so E can be estimated by least squares from
E_L.  If that is accurate, the Phase 3 test can run on the Release 2 L2B
archive without the L1A delivery.

Test on the VEIT sample (HYPSTAR 122304), where both the L1A E (truth) and
the L1C irradiance exist: solve with the true E grid (from L1A) and with the
lab-polynomial E grid, and report the recovery error, total and in the
second-difference (wiggle) band that Phase 3 measures.

Run from the repository root:  python wavecal/invert_l2_irradiance.py
"""
import glob
import os

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import lsqr

from hypernet import whn_l1a as wl

WAV = os.path.join(os.environ['OS_COLOR'], 'WATERHYPERNET', 'Wavelengths')
LAB = os.path.join(os.environ['OS_COLOR'], 'hypernet', 'wiggles', 'lab_data', 'derived')


def interp_matrix(wL, wE):
    """Sparse M with M @ E == np.interp(wL, wE, E) (edges clamped)."""
    j = np.clip(np.searchsorted(wE, wL) - 1, 0, len(wE) - 2)
    w = np.clip((wL - wE[j]) / (wE[j + 1] - wE[j]), 0.0, 1.0)
    rows = np.repeat(np.arange(len(wL)), 2)
    cols = np.stack([j, j + 1], axis=1).ravel()
    vals = np.stack([1 - w, w], axis=1).ravel()
    return sparse.csr_matrix((vals, (rows, cols)), shape=(len(wL), len(wE)))


def d2(y, wave):
    """Second difference on adjacent pixels."""
    return y[2:] - 2 * y[1:-1] + y[:-2]


def report(label, E_true, E_hat, wE, lo=400., hi=900.):
    m = (wE >= lo) & (wE <= hi)
    rel = (E_hat - E_true) / E_true
    md = m[1:-1]
    r2 = d2(E_hat, wE)[md] - d2(E_true, wE)[md]
    s2 = d2(E_true, wE)[md]
    print('%-34s rel err rms %.2e  max %.2e | d2 err rms / d2(E) rms = %.3f'
          % (label, np.sqrt(np.mean(rel[m] ** 2)), np.abs(rel[m]).max(),
             np.sqrt(np.mean(r2 ** 2)) / np.sqrt(np.mean(s2 ** 2))))


def main():
    irr = wl.load_l1a_irr(glob.glob(os.path.join(WAV, '*L1A_IRR*.nc'))[0])
    l1c = wl.load_l1c(glob.glob(os.path.join(WAV, '*L1C_ALL*.nc'))[0])
    wE, E = irr['wave'], irr['mean']
    wL = l1c['wave']
    EL = np.nanmean(l1c['irradiance'], axis=1)
    print('n_E = %d, n_L = %d' % (len(wE), len(wL)))

    M = interp_matrix(wL, wE)
    fr = np.abs(M @ E - EL) / EL
    mid = (wL >= 400) & (wL <= 900)
    print('forward check: max |M E - L1C| / L1C = %.2e (400-900 nm), %.2e (all;'
          ' worst at %.1f nm)' % (fr[mid].max(), fr.max(), wL[np.argmax(fr)]))
    sv = np.linalg.svd(M.toarray(), compute_uv=False)
    print('singular values of M: max %.3f  min %.3e  condition %.2e'
          % (sv.max(), sv.min(), sv.max() / sv.min()))

    # 1. true E grid, noiseless (round-off of the stored float32 E_L only)
    E_hat = lsqr(M, EL, atol=1e-14, btol=1e-14, iter_lim=20000)[0]
    report('true grid, from L1C', E, E_hat, wE)

    # 2. lab-polynomial grids (2023_01 = what Release 2 uses; 2024_09 for contrast)
    for f in sorted(glob.glob(os.path.join(LAB, 'egrid_122304_*.csv'))):
        g = pd.read_csv(f)['wavelength_E'].values
        if len(g) != len(wE):
            print('%s: %d px, skipped' % (os.path.basename(f), len(g)))
            continue
        print('%s: max |grid - L1A grid| = %.2e nm'
              % (os.path.basename(f), np.abs(g - wE).max()))
        Mg = interp_matrix(wL, g)
        E_hat = lsqr(Mg, EL, atol=1e-14, btol=1e-14, iter_lim=20000)[0]
        report('lab grid ' + os.path.basename(f)[13:23], E, E_hat, wE)

    # 3. sensitivity to errors NOT of the form M @ E (e.g. if L2 averaged
    #    scans after interpolating with per-scan filtering, or rounding).
    #    Noise on E itself passes through M and back unamplified.
    rng = np.random.default_rng(1)
    for s in (1e-4, 1e-3):
        ELn = EL * (1 + s * rng.standard_normal(len(EL)))
        E_hat = lsqr(M, ELn, atol=1e-14, btol=1e-14, iter_lim=20000)[0]
        report('true grid, %.1f%% noise on E_L' % (100 * s), E, E_hat, wE)


if __name__ == '__main__':
    main()

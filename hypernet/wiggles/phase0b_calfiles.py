"""Phase 0b, task 13: what the HYPSTAR calibration files contain.

For every ``HYPERNETS_CAL_HYPSTAR_<id>_{RAD,IRR}_v*.nc`` under
``$OS_COLOR/WATERHYPERNET/Wavelengths`` (all subfolders):

1. list every variable (name, dims, shape, dtype, units, long_name, and a
   value summary) and every global attribute;
2. pick out anything wavelength- or SRF-like -- names matching wavelength,
   wave_coef, coefficient, bandwidth, fwhm, srf, slit, lsf, spectral_response,
   resolution, sigma -- and summarise it;
3. compare with the fits of task 10 for the same instrument: the cal-file
   wavelength grid (if any) against the L1A ``wavelength`` of the matching
   sequences (max |difference|), and any bandwidth / FWHM variable against
   the fitted template FWHM(lambda) at 400-850 nm (E for IRR, Ld for RAD).

Also lists any other file in the tree that looks like lab line-spread data
(names containing lsf, slit, srf, line_spread, fwhm).

Outputs: ``hypernet/wiggles/phase0_calfiles.csv`` (one row per cal file),
``hypernet/wiggles/phase0_calfiles_variables.csv`` (one row per variable) and, if fits
exist, ``hypernet/wiggles/phase0_calfiles_compare.csv``.  With no cal files present it
says so and writes nothing.

Run from the repository root: ``python -m hypernet.wiggles.phase0b_calfiles``.
"""
import glob
import os
import re

import numpy as np
import pandas as pd

from hypernet import whn_l1a as wl  # noqa: E402
from hypernet.wiggles import DATA_DIR, FIGDIR, OUT, REPO, WIGGLES_DIR  # noqa: F401

CAL_RE = re.compile(r'HYPERNETS_CAL_HYPSTAR_(?P<inst>\d+)_(?P<kind>RAD|IRR)_(?P<ver>v[\d.]+)\.nc$')
SRF_LIKE = re.compile(r'wavelength|wave_?coef|coefficient|bandwidth|fwhm|srf|slit|lsf|'
                      r'spectral_response|resolution|sigma', re.I)
LSF_FILE = re.compile(r'lsf|slit|srf|line_spread|fwhm', re.I)


def find_calfiles(root=None):
    root = wl.wavelengths_root(root)
    out = []
    for p in sorted(glob.glob(os.path.join(root, '**', 'HYPERNETS_CAL_*.nc'), recursive=True)):
        m = CAL_RE.search(os.path.basename(p))
        out.append(dict(path=p, instrument='HYPSTAR_' + m['inst'] if m else None,
                        kind=m['kind'] if m else None, version=m['ver'] if m else None))
    lsf = [p for p in glob.glob(os.path.join(root, '**', '*'), recursive=True)
           if os.path.isfile(p) and LSF_FILE.search(os.path.basename(p))
           and 'HYPERNETS_CAL_' not in os.path.basename(p)]
    return pd.DataFrame(out), lsf


def _summary(a):
    a = np.asarray(a)
    if a.dtype.kind in 'fiu' and a.size:
        f = a[np.isfinite(a)] if a.dtype.kind == 'f' else a
        if not f.size:
            return 'all NaN'
        u = np.unique(f)
        if u.size <= 3:
            return 'values %s' % u.tolist()
        return 'min %.6g, median %.6g, max %.6g, n_unique %d' % (f.min(), np.median(f), f.max(),
                                                                  u.size)
    return str(a.ravel()[:3].tolist())[:80]


def inspect(path):
    import xarray as xr
    ds = xr.open_dataset(path)
    try:
        variables = []
        for v in ds.variables:
            x = ds[v]
            variables.append(dict(variable=v, dims=','.join(x.dims), shape=str(x.shape),
                                  dtype=str(x.dtype), units=x.attrs.get('units'),
                                  long_name=x.attrs.get('long_name'),
                                  srf_like=bool(SRF_LIKE.search(v) or SRF_LIKE.search(
                                      str(x.attrs.get('long_name', '')))),
                                  summary=_summary(x.values)))
        attrs = {k: str(v)[:200] for k, v in ds.attrs.items()}
        wave = ds['wavelength'].values.astype(float) if 'wavelength' in ds else None
        bw = None
        for name in ('bandwidth', 'fwhm', 'FWHM'):
            if name in ds:
                bw = ds[name].values.astype(float)
                break
    finally:
        ds.close()
    return variables, attrs, wave, bw


def compare(cal, wave, bw, index, models):
    """Cal-file wavelengths vs L1A grids; bandwidth vs fitted FWHM."""
    rows = []
    prod = 'L1A_IRR' if cal['kind'] == 'IRR' else 'L1A_RAD'
    if index is not None and wave is not None:
        seq = index[(index.get('system_id') == cal['instrument']) & index[prod].notna()]
        for _, r in seq.iterrows():
            d = wl.load_l1a_irr(r[prod]) if prod == 'L1A_IRR' else wl.load_l1a_rad(r[prod])
            w = d['wave']
            dmax = float(np.max(np.abs(w - wave))) if w.size == wave.size else np.nan
            rows.append(dict(cal=os.path.basename(cal['path']), sequence=r['sequence_time'],
                             quantity='wavelength grid vs L1A', n_cal=wave.size, n_l1a=w.size,
                             max_abs_diff_nm=dmax, cal_period=r.get('cal_period')))
    if models is not None and bw is not None and wave is not None:
        ch = 'E' if cal['kind'] == 'IRR' else 'Ld'
        mm = models[(models['system_id'] == cal['instrument']) & (models['method'] == 'template')
                    & (models['channel'] == ch)]
        for lam in (400, 500, 600, 700, 850):
            j = np.argmin(np.abs(wave - lam))
            rows.append(dict(cal=os.path.basename(cal['path']), quantity='bandwidth vs fitted FWHM',
                             lam=lam, cal_value=float(bw[j]), fitted_median=float(
                                 mm['fwhm_%d' % lam].median()) if len(mm) else np.nan,
                             n_sequences=len(mm)))
    return rows


def main():
    cals, lsf = find_calfiles()
    print('lab line-spread-like files: %s' % (lsf or 'none'))
    if not len(cals):
        print('no HYPERNETS_CAL_HYPSTAR_*_{RAD,IRR}_v*.nc files under %s -- nothing to do '
              '(waiting on the delivery)' % wl.wavelengths_root())
        return
    ip = os.path.join(OUT, 'request_index.parquet')
    mp = os.path.join(OUT, 'srf_models.parquet')
    index = pd.read_parquet(ip) if os.path.exists(ip) else None
    models = pd.read_parquet(mp) if os.path.exists(mp) else None
    files, variables, comp = [], [], []
    for _, c in cals.iterrows():
        v, attrs, wave, bw = inspect(c['path'])
        for r in v:
            variables.append(dict(cal=os.path.basename(c['path']), **r))
        files.append(dict(c, n_variables=len(v), n_srf_like=sum(r['srf_like'] for r in v),
                          srf_like=';'.join(r['variable'] for r in v if r['srf_like']),
                          has_wavelength=wave is not None, has_bandwidth=bw is not None,
                          bandwidth_summary=_summary(bw) if bw is not None else None,
                          attrs='; '.join('%s=%s' % kv for kv in attrs.items())))
        comp += compare(c, wave, bw, index, models)
    pd.DataFrame(files).to_csv(os.path.join(WIGGLES_DIR, 'phase0_calfiles.csv'), index=False)
    pd.DataFrame(variables).to_csv(os.path.join(WIGGLES_DIR, 'phase0_calfiles_variables.csv'),
                                   index=False)
    if comp:
        pd.DataFrame(comp).to_csv(os.path.join(WIGGLES_DIR, 'phase0_calfiles_compare.csv'),
                                  index=False, float_format='%.5g')
    pd.set_option('display.width', 220)
    print(pd.DataFrame(files)[['instrument', 'kind', 'version', 'n_variables', 'srf_like',
                               'bandwidth_summary']].to_string(index=False))
    print('wrote hypernet/wiggles/phase0_calfiles*.csv')


if __name__ == '__main__':
    main()

"""Phase 0b, task 10: SRF fits of every delivered sequence.

For every sequence in ``request_index.parquet`` (task 9) with both L1A files,
runs :func:`hypernet.whn_srf.fit_sequence`: empirical line fits
(``srf.fit_lines``) and HSRS template fits (``srf.fit_template_windows``, veil
fixed at 0 as in the shipped VEIT models) on the mean E, Ld and Lu, and the
FWHM(lambda) model per method and channel.  A sequence that fails is logged
in ``status`` and skipped; it does not stop the run.

Outputs, under ``$OS_COLOR/hypernet/wiggles/phase0/``, keyed by site,
instrument (``system_id``), cal period and sequence:

- ``line_fits.parquet`` -- empirical, one row per sequence x channel x line;
- ``template_fits.parquet`` -- one row per sequence x channel x window;
- ``srf_models.parquet`` -- one row per sequence x method x channel
  (coefficients, errors, flattened covariance, FWHM at 400-850 nm, offset);
- ``scan_fits.parquet`` -- per-scan empirical fits, with ``--per-scan``;
- ``fit_status.csv`` -- per sequence: ok / error message, runtime.

Run from the repository root::

    python -m hypernet.wiggles.phase0b_fit_all [--per-scan] [--limit N] [--no-template]
"""
import argparse
import os
import time
import traceback

import pandas as pd

from hypernet import whn_srf  # noqa: E402
from hypernet.wiggles import DATA_DIR, FIGDIR, OUT, REPO, WIGGLES_DIR  # noqa: F401

KEYS = ['site_code', 'sequence_time', 'system_id', 'cal_period', 'status', 'row',
        'priority', 'sza_l1a', 'sky', 'month', 'water_type']


def sequences(index):
    """Delivered (requested or extra) sequences with both L1A files."""
    ok = index['status'].isin(['delivered', 'extra']) & \
        index['L1A_IRR'].notna() & index['L1A_RAD'].notna()
    return index[ok].reset_index(drop=True)


def run(index, per_scan=False, template=True, limit=None, verbose=True):
    seqs = sequences(index)
    if limit:
        seqs = seqs.head(limit)
    ref = whn_srf.load_reference('air') if template else None
    lines, temps, mods, scans, status = [], [], [], [], []
    for i, r in seqs.iterrows():
        keys = {k: r.get(k) for k in KEYS}
        keys['request_status'] = keys.pop('status')
        t0 = time.time()
        try:
            res = whn_srf.fit_files({'L1A_IRR': r['L1A_IRR'], 'L1A_RAD': r['L1A_RAD']},
                                    ref=ref, template=template, veil=False,
                                    per_scan=per_scan)
            for key, lst in (('lines', lines), ('template', temps), ('scans', scans)):
                if res[key] is not None:
                    lst.append(res[key].assign(**keys))
            mods.append(whn_srf.models_table(res['models']).assign(**keys))
            st = 'ok'
        except Exception as ex:  # keep going; the message goes in fit_status
            st = '%s: %s' % (type(ex).__name__, ex)
            if verbose:
                traceback.print_exc()
        dt = time.time() - t0
        status.append(dict(keys, fit_status=st, seconds=round(dt, 1)))
        if verbose:
            print('[%d/%d] %s %s %s: %s (%.1f s)' % (i + 1, len(seqs), r['site_code'],
                                                      r['sequence_time'], r.get('system_id'),
                                                      st, dt))
    cat = lambda x: pd.concat(x, ignore_index=True) if x else pd.DataFrame()  # noqa: E731
    return dict(lines=cat(lines), template=cat(temps), models=cat(mods), scans=cat(scans),
                status=pd.DataFrame(status))


def main(argv=None):
    ap = argparse.ArgumentParser(description='Phase 0b batch SRF fits')
    ap.add_argument('--per-scan', action='store_true', help='also fit every scan')
    ap.add_argument('--no-template', action='store_true', help='skip the HSRS template fits')
    ap.add_argument('--limit', type=int, default=None, help='first N sequences only')
    a = ap.parse_args(argv)
    index = pd.read_parquet(os.path.join(OUT, 'request_index.parquet'))
    res = run(index, per_scan=a.per_scan, template=not a.no_template, limit=a.limit)
    for key, name in (('lines', 'line_fits'), ('template', 'template_fits'),
                      ('models', 'srf_models'), ('scans', 'scan_fits')):
        if len(res[key]):
            res[key].to_parquet(os.path.join(OUT, name + '.parquet'))
    res['status'].to_csv(os.path.join(OUT, 'fit_status.csv'), index=False)
    st = res['status']
    print('\n%d sequences: %d ok, %d failed' % (len(st), (st['fit_status'] == 'ok').sum(),
                                                 (st['fit_status'] != 'ok').sum()))
    m = res['models']
    if len(m):
        cols = ['site_code', 'sequence_time', 'system_id', 'method', 'channel', 'npts',
                'chi2_nu', 'fwhm_400', 'fwhm_600', 'fwhm_850', 'offset', 'offset_err']
        print(m[cols].to_string(index=False, float_format='%.3f'))
    print('wrote', OUT)


if __name__ == '__main__':
    main()

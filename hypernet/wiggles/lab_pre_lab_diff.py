"""Differences between the committed products and their ``*_pre_lab`` copies
(lab data prompts, task 3: "keep what you replace ... and report the
differences row by row").

Parts (``--part``, default all that exist):

- ``grids``: ``hypernet/data/hypstar_grids.npz`` vs ``hypstar_grids_pre_lab.npz``,
  key by key;
- ``phase0c``: ``hypernet/wiggles/phase0c_summary.csv`` and
  ``phase0c_stability.csv`` vs their ``_pre_lab`` copies: E FWHM at
  400/500/600/700 nm and Delta FWHM (Ld - E at 450 and 600 nm) per
  instrument x calibration period;
- ``models``: ``release2_srf_models.json`` vs ``_pre_lab``: FWHM at 450/550/700
  nm per key x channel;
- ``table``: ``srf_table.csv`` vs ``srf_table_pre_lab.csv``, per row;
- ``pattern``: ``wavescale_window_pattern.csv`` vs ``_pre_lab``, per window;
- ``wavescale``: ``hypernet/wiggles/phase2_wavescale.csv`` vs ``_pre_lab``: the
  median shift / stretch per instrument x channel.

Each part writes ``$OS_COLOR/hypernet/wiggles/lab_data/pre_lab_diff/<part>.csv``
and prints it.

Run from the repository root:
``python -m hypernet.wiggles.lab_pre_lab_diff [--part grids|phase0c|models|table|pattern]``
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

from hypernet import labcal, srf
from hypernet.wiggles import DATA_DIR, WIGGLES_DIR

OUT = os.path.join(labcal.lab_root(), 'pre_lab_diff')
LAMS = (400, 500, 600, 700)


def _pre(path):
    b, e = os.path.splitext(path)
    return b + '_pre_lab' + e


def grids():
    new = np.load(os.path.join(DATA_DIR, 'hypstar_grids.npz'))
    old = np.load(os.path.join(DATA_DIR, 'hypstar_grids_pre_lab.npz'))
    rows = []
    for k in sorted(set(new.files) | set(old.files)):
        d = dict(key=k, in_old=k in old.files, in_new=k in new.files)
        if d['in_old'] and d['in_new']:
            a, b = old[k], new[k]
            d.update(n_old=a.size, n_new=b.size, same_shape=a.shape == b.shape,
                     identical=bool(a.shape == b.shape and np.array_equal(a, b)),
                     max_abs_diff=float(np.abs(a - b).max()) if a.shape == b.shape else np.nan)
        elif d['in_new']:
            d.update(n_new=new[k].size)
        rows.append(d)
    return pd.DataFrame(rows)


def phase0c():
    p = os.path.join(WIGGLES_DIR, 'phase0c_stability.csv')
    new, old = pd.read_csv(p), pd.read_csv(_pre(p))
    keys = ['system_id', 'cal_period', 'channel']
    cols = ['fwhm_%d_mean' % l for l in LAMS]
    m = old[keys + cols + ['e_grid']].merge(new[keys + cols + ['e_grid']], on=keys,
                                            suffixes=('_old', '_new'))
    for c in cols:
        m['d_' + c] = m[c + '_new'] - m[c + '_old']
    s = os.path.join(WIGGLES_DIR, 'phase0c_summary.csv')
    sn, so = pd.read_csv(s), pd.read_csv(_pre(s))
    k2 = ['system_id', 'cal_period']
    q = ['dfwhm_450_Ld_E', 'dfwhm_600_Ld_E', 'fwhm_450_E', 'fwhm_600_E']
    d = so[k2 + q].merge(sn[k2 + q], on=k2, suffixes=('_old', '_new'))
    for c in q:
        d['d_' + c] = d[c + '_new'] - d[c + '_old']
    return m, d


def _models(path):
    doc = json.load(open(path))
    out = []
    for key, chans in doc['models'].items():
        if '@' not in key:
            continue
        for ch, dd in chans.items():
            m = srf.SRFModel.from_dict(dd)
            out.append(dict(key=key, channel=ch, shape=m.shape, shape_par=m.shape_par,
                            e_grid=m.meta.get('e_grid'), n_sequences=m.meta.get('n_sequences'),
                            **{'fwhm_%d' % l: m.fwhm(l) for l in (450, 550, 700)}))
    return pd.DataFrame(out)


def models():
    p = os.path.join(DATA_DIR, 'release2_srf_models.json')
    new, old = _models(p), _models(_pre(p))
    m = old.merge(new, on=['key', 'channel'], how='outer', suffixes=('_old', '_new'))
    for l in (450, 550, 700):
        m['d_fwhm_%d' % l] = m['fwhm_%d_new' % l] - m['fwhm_%d_old' % l]
    return m


def table():
    p = os.path.join(DATA_DIR, 'srf_table.csv')
    new, old = pd.read_csv(p), pd.read_csv(_pre(p))
    keys = ['instrument', 'channel', 'cal_period', 'source']
    cols = ['e_grid', 'fwhm_400', 'fwhm_550', 'fwhm_700', 'offset', 'delta_fwhm', 'method',
            'delta_fwhm_l1a_est', 'method_l1a_est', 'method_sensitive', 'scale_shift',
            'scale_stretch']
    extra = [c for c in ('shape', 'shape_p', 'shape_source', 'shape_p_field',
                         'lab_minus_table_405', 'lab_minus_table_436', 'lab_minus_table_546')
             if c in new]
    m = old[keys + cols].merge(new[keys + cols + extra], on=keys, how='outer',
                               suffixes=('_old', '_new'))
    for c in ('fwhm_400', 'fwhm_550', 'fwhm_700', 'offset', 'delta_fwhm', 'delta_fwhm_l1a_est',
              'scale_shift', 'scale_stretch'):
        m['d_' + c] = m[c + '_new'] - m[c + '_old']
    m['method_changed'] = m['method_old'] != m['method_new']
    return m


def pattern():
    p = os.path.join(DATA_DIR, 'wavescale_window_pattern.csv')
    new, old = pd.read_csv(p), pd.read_csv(_pre(p))
    m = old.merge(new, on='lam_center', how='outer', suffixes=('_old', '_new'))
    m['d_pattern'] = m['pattern_new'] - m['pattern_old']
    return m


def wavescale():
    p = os.path.join(WIGGLES_DIR, 'phase2_wavescale.csv')
    new, old = pd.read_csv(p), pd.read_csv(_pre(p))
    keys = ['source', 'instrument', 'channel']
    cols = ['shift_med', 'stretch_med', 'max_err_med', 'frac_ok', 'chi2_med']
    m = old[keys + cols].merge(new[keys + cols], on=keys, suffixes=('_old', '_new'))
    for c in cols:
        m['d_' + c] = m[c + '_new'] - m[c + '_old']
    return m


def main(argv=None):
    ap = argparse.ArgumentParser(description='committed products vs their _pre_lab copies')
    ap.add_argument('--part', choices=('grids', 'phase0c', 'models', 'table', 'pattern', 'wavescale'),
                    action='append')
    a = ap.parse_args(argv)
    os.makedirs(OUT, exist_ok=True)
    pd.set_option('display.width', 260)
    pd.set_option('display.max_columns', 60)
    for part in a.part or ('grids', 'phase0c', 'models', 'table', 'pattern', 'wavescale'):
        try:
            res = globals()[part]()
        except FileNotFoundError as ex:
            print('[%s] skipped: %s' % (part, ex))
            continue
        res = res if isinstance(res, tuple) else (res,)
        for i, r in enumerate(res):
            name = part if i == 0 else '%s_%d' % (part, i)
            r.to_csv(os.path.join(OUT, name + '.csv'), index=False, float_format='%.5g')
            print('\n[%s]' % name)
            print(r.to_string(index=False, float_format='%.4f'))


if __name__ == '__main__':
    main()

"""Collect the numbers for the Phase 0a report slides (task 8b) into JSON.

Everything comes from the committed Phase 0a products (tasks 3b, 5-7):
``hypernet/wiggles/phase0_veit_{lines,scans,budget}.csv`` and
``hypernet/data/veit_srf_model.json``.  The deck itself is built by
``build_wiggles_phase0_report.js`` from the JSON written here.

Run from the repository root::

    python docs/slides/wiggles_phase0_report_data.py
    node docs/slides/build_wiggles_phase0_report.js
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from PIL import Image

# as in docs/whn_figures.py: make ``hypernet`` importable without pip install -e
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from hypernet import srf  # noqa: E402
from hypernet.wiggles import DATA_DIR, FIGDIR, WIGGLES_DIR  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'wiggles_phase0_report_data.json')

# lines shown in the tables (a readable subset at >= 20 pt)
SHOW = ['Ca K', 'G band', 'H beta', 'Mg b', 'Na D', 'H alpha', 'Ca II 854.2', 'Ca II 866.2']
LABEL = {'Ca K': 'Ca K 393', 'G band': 'G band 431', 'H beta': 'Hβ 486', 'Mg b': 'Mg b 517',
         'Na D': 'Na D 589', 'H alpha': 'Hα 656', 'Ca II 854.2': 'Ca II 854',
         'Ca II 866.2': 'Ca II 866', 'Ca H/K': 'Ca H/K 395',
         'Ca II 849.8/854.2': 'Ca II 850/854'}


def task5():
    f = pd.read_csv(os.path.join(WIGGLES_DIR, 'phase0_veit_lines.csv'))
    f = f[f['err_kind'] == 'flat'].set_index(['name', 'channel'])
    rows = []
    for n in SHOW:
        E, Ld, Lu = (f.loc[(n, c)] for c in ('E', 'Ld', 'Lu'))
        d = Ld['fwhm'] - E['fwhm']
        e = np.hypot(Ld['fwhm_err'], E['fwhm_err'])
        rows.append([LABEL[n], '%.2f' % E['fwhm'], '%.2f' % Ld['fwhm'], '%.2f' % Lu['fwhm'],
                     '%+.2f ± %.2f' % (d, e), '%.0fσ' % (d / e)])
    return rows


def task6():
    s = pd.read_csv(os.path.join(WIGGLES_DIR, 'phase0_veit_scans.csv')).set_index(
        ['name', 'channel'])
    rows = []
    for n in SHOW:
        a, b, e = (s.loc[(n, c)] for c in ('Ld', 'Lu', 'E'))
        err = max(np.hypot(a['fwhm_mean_err_fit'], b['fwhm_mean_err_fit']),
                  np.hypot(a['fwhm_mean_err_scan'], b['fwhm_mean_err_scan']))
        d = a['fwhm_mean'] - b['fwhm_mean']
        ew = lambda r: r['depth_mean'] * r['sigma_mean']  # noqa: E731  (sqrt(2 pi) cancels)
        rows.append([LABEL[n], '%.1f' % a['fwhm_ratio'], '%+.2f ± %.2f' % (d, err),
                     '%.2f' % (a['depth_mean'] / e['depth_mean']), '%.2f' % (ew(a) / ew(e))])
    return rows


def task7():
    b = pd.read_csv(os.path.join(WIGGLES_DIR, 'phase0_veit_budget.csv')).set_index('line')
    rows = []
    for n in ['Ca H/K', 'G band', 'H beta', 'Mg b', 'Na D', 'H alpha']:
        r = b.loc[n]
        rows.append([LABEL[n], '%.1e' % r['rms_meas_LdEd'], '%.1e' % r['rms_H2'],
                     '%.1e' % r['rms_H1'],
                     '%.2f ± %.2f' % (r['alpha2_LdEd'], r['alpha2_err_LdEd']),
                     '%.2f ± %.2f' % (r['alpha2_rho'], r['alpha2_err_rho'])])
    c = b[b['use_for_srf'] & b['fit_ok']]
    stats = {}
    for pre in ('LdEd', 'LuEd', 'rho'):
        w = 1 / c['alpha2_err_' + pre] ** 2
        stats[pre] = ['%.2f' % (np.sum(w * c['alpha2_' + pre]) / w.sum()),
                      '%.2f' % (1 / np.sqrt(w.sum()))]
    stats['h2_over_h1'] = '%.0f' % (c['rms_H2'] / c['rms_H1']).median()
    return rows, stats


def task3b():
    mods = srf.load_srf_models(os.path.join(DATA_DIR, 'veit_srf_model.json'))
    rows = []
    for lam in (400, 500, 600, 700, 850):
        E, L = mods['E'], mods['Ld']
        d = L.fwhm(lam) - E.fwhm(lam)
        e = np.hypot(L.fwhm_err(lam), E.fwhm_err(lam))
        rows.append(['%d' % lam, '%.2f' % E.fwhm(lam), '%.2f' % L.fwhm(lam),
                     '%+.2f ± %.2f' % (d, e)])
    return rows


def images():
    out = {}
    for k in ('veit_fwhm_vs_lambda', 'veit_scan_scatter', 'veit_budget', 'veit_template_fwhm'):
        p = os.path.join(FIGDIR, k + '.png')
        w, h = Image.open(p).size
        out[k] = dict(path=p, aspect=w / h)
    return out


def main():
    t7, st7 = task7()
    data = dict(task5=task5(), task6=task6(), task7=t7, task7_stats=st7, task3b=task3b(),
                images=images())
    with open(OUT, 'w') as fh:
        json.dump(data, fh, indent=1, ensure_ascii=False)
    print(json.dumps({k: v for k, v in data.items() if k != 'images'}, indent=1,
                     ensure_ascii=False))
    print('wrote', OUT)


if __name__ == '__main__':
    main()

"""Numbers and figure sizes for the Phase 1 slides (Phase 1 task 13).

Reads the committed Phase 1 tables in ``hypernet/wiggles/`` and the figures in
``hypernet/wiggles/figs/phase1/``, and writes ``slides/wiggles_phase1_data.json``
for ``slides/build_wiggles_phase1.js``.  Medians over the 25 scenes, h = 1 nm.

Run from the repository root: ``python slides/wiggles_phase1_data.py``.
"""
import json
import os

import pandas as pd
from PIL import Image

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WIG = os.path.join(REPO, 'hypernet', 'wiggles')
FIG = os.path.join(WIG, 'figs', 'phase1')
I4, I5 = 'HYPSTAR_122304', 'HYPSTAR_122305'


def sci(x):
    """1.45e-4 -> '1.5×10⁻⁴'."""
    m, e = ('%.1e' % x).split('e')
    sup = str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
    return '%s×10%s' % (m, str(int(e)).translate(sup))


def pct(x):
    return ('%.0f %%' % (100 * x)).replace('-', '−')


PRETTY = {'Ca K': 'Ca K', 'G band': 'G band', 'H beta': 'Hβ', 'Na D': 'Na D',
          'O2 B': 'O₂-B', 'O2 A': 'O₂-A', 'H alpha': 'Hα', 'H2O': 'H₂O'}


def main():
    out = {'images': {}}
    for f in sorted(os.listdir(FIG)):
        if f.endswith('.png'):
            w, h = Image.open(os.path.join(FIG, f)).size
            out['images'][f[:-4]] = dict(path=os.path.join(FIG, f), aspect=w / h)

    # task 9: cases x methods
    t = pd.read_csv(os.path.join(WIG, 'phase1_twin_metrics.csv'))
    t = t[t['h'] == 1.0].set_index(['instrument', 'case', 'method'])
    rows = []
    for case, lab in (('i_offset', '(i) grid offset only'),
                      ('ii_mismatch', '(ii) SRF mismatch only'),
                      ('iii_both', '(iii) offset + mismatch'),
                      ('iv_shift+0.1', '(iv) E λ shift 0.1 nm'),
                      ('iv_shift+0.3', '(iv) E λ shift 0.3 nm'),
                      ('v_wrongsrf_E-0.3', '(v) FWHM_E off 0.3 nm'),
                      ('vi_wrongemod_pwv_factor', '(vi) Emod PWV × 2'),
                      ('vii_noise', '(vii) VEIT noise')):
        g = lambda m, c='rms_line': t.loc[(I4, case, m), c]  # noqa: E731
        rows.append([lab, sci(g('linear')), pct(g('ruddick2023', 'reduction')),
                     pct(g('srf', 'reduction'))])
    out['t9_rows'] = rows
    out['t9'] = dict(
        true_line=sci(t.loc[(I4, 'iii_both', 'srf'), 'rms_true_line']),
        lin_iii=sci(t.loc[(I4, 'iii_both', 'linear'), 'rms_line']),
        lin_i=sci(t.loc[(I4, 'i_offset', 'linear'), 'rms_line']),
        offset_frac=pct(t.loc[(I4, 'i_offset', 'linear'), 'rms_line'] /
                        t.loc[(I4, 'iii_both', 'linear'), 'rms_line']),
        ratio_lin_true='%.0f×' % (t.loc[(I4, 'iii_both', 'linear'), 'rms_line'] /
                                  t.loc[(I4, 'iii_both', 'srf'), 'rms_true_line']),
        lin5_iii=sci(t.loc[(I5, 'iii_both', 'linear'), 'rms_line']),
        rud5_iii=pct(t.loc[(I5, 'iii_both', 'ruddick2023'), 'reduction']),
        srf5_iii=pct(t.loc[(I5, 'iii_both', 'srf'), 'reduction']),
        srf5_v='%.0f–%.0f×' % (lambda v: (min(v), max(v)))([
            t.loc[(I5, c, 'srf'), 'rms_line'] / t.loc[(I5, c, 'linear'), 'rms_line']
            for c in ('v_wrongsrf_E-0.3', 'v_wrongsrf_E+0.3', 'v_wrongsrf_L-0.3',
                      'v_wrongsrf_L+0.3')]),
        rho5_i='%.2f %%' % (100 * t.loc[(I5, 'i_offset', 'linear'), 'rms_rho_line']),
        rho5_iii='%.2f %%' % (100 * t.loc[(I5, 'iii_both', 'linear'), 'rms_rho_line']),
        lo_lin_vii='%.2f' % t.loc[(I4, 'vii_noise', 'linear'), 'line_over_away'],
        lo_srf_vii='%.2f' % t.loc[(I4, 'vii_noise', 'srf'), 'line_over_away'],
    )

    # task 10: controls
    c = pd.read_csv(os.path.join(WIG, 'phase1_controls.csv'))
    c1 = c[c['h'] == 1.0]
    ci = c1[(c1['instrument'] == I4) & (c1['case'] == 'iii_both')].set_index(
        ['control', 'method'])
    out['t10'] = dict(
        max_over_noise='%.2f' % c['ctl_err_over_noise'].max(),
        amp_fl='%.2f' % ci.loc[('fl', 'srf'), 'ctl_amp_over_noise'],
        amp_raman='%.2f' % ci.loc[('raman', 'srf'), 'ctl_amp_over_noise'],
        rows=[[m, pct(ci.loc[('fl', m), 'ctl_err_over_amp']),
               pct(ci.loc[('raman', m), 'ctl_err_over_amp'])]
              for m in ('linear', 'ruddick2023', 'srf')])

    # task 10: degradation
    d = pd.read_csv(os.path.join(WIG, 'phase1_degradation.csv'))
    d = d[(d['h'] == 1.0) & (d['instrument'] == I4) & (d['method'] == 'srf')]
    red = lambda k, v: d[(d['kind'] == k) & (abs(d['value'] - v) < 1e-9)]['reduction'].mean()  # noqa: E731,E501
    out['t10d'] = dict(
        shift05=pct(red('shift', 0.05)), shift10=pct(red('shift', 0.1)),
        shift30=pct(red('shift', 0.3)), stretch05=pct(red('stretch', 0.05)),
        e05=pct((red('dfwhm_E', 0.05) + red('dfwhm_E', -0.05)) / 2),
        e10=pct((red('dfwhm_E', 0.1) + red('dfwhm_E', -0.1)) / 2),
        l10=pct((red('dfwhm_L', 0.1) + red('dfwhm_L', -0.1)) / 2),
        emod_min=pct(d[d['kind'].isin(['sza_offset', 'pwv_factor'])]['reduction'].min()))

    # task 11: prediction
    p = pd.read_csv(os.path.join(WIG, 'phase1_prediction.csv'))
    p = p[(p['h'] == 1.0) & (p['method'] == 'srf')]
    dfs = (0.1, 0.2, 0.3, 0.5, 0.7, 1.0)
    pa = p[p['line'] == 'all']
    get = lambda sc, col, x: pa[(pa['scenario'] == sc) & (abs(pa['dfwhm'] - x) < 1e-9)][col].iloc[0]  # noqa: E731,E501
    out['t11'] = dict(
        dfwhm=['%g' % x for x in dfs],
        rows=[['noise only'] + [pct(get('noise', 'reduction_obs', x)) for x in dfs],
              ['+ 0.05 nm λ and FWHM errors'] +
              [pct(get('realistic', 'reduction_obs', x)) for x in dfs],
              ['perfect correction'] + [pct(get('noise', 'reduction_obs_max', x))
                                        for x in dfs]])
    pl = p[(p['scenario'] == 'noise') & (abs(p['dfwhm'] - 0.5) < 1e-9)].set_index('line')
    out['t11']['lines'] = [[PRETTY[n], '%.2f' % pl.loc[n, 'depth'], pct(pl.loc[n, 'reduction'])]
                           for n in ('Ca K', 'G band', 'H beta', 'Na D', 'O2 B', 'O2 A',
                                     'H alpha', 'H2O')]
    path = os.path.join(REPO, 'slides', 'wiggles_phase1_data.json')
    with open(path, 'w') as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
    print(json.dumps({k: v for k, v in out.items() if k != 'images'}, indent=1,
                     ensure_ascii=False))
    print('wrote', path)


if __name__ == '__main__':
    main()

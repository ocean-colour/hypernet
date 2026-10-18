"""Numbers for the provisional gate G0 (lab data prompts, task 6).

Collects, from the committed tables, what ``claude_prompts/wiggles/gate_G0.md``
quotes for the plan's G0 criteria (``docs/wiggles_planning.md`` §4):

- G0(b) field stability (``phase0c_stability*.csv``, post-lab, own E grids):
  the largest FWHM std per instrument x period x channel against 0.2 nm; the
  SZA, sky-index and season slopes turned into a change over each group's
  range; the recalibration and site pairs;
- G0(b) lab stability across epochs (``lab_srf.csv``, ``lab_wl_coefs.csv``):
  FWHM at 405/436/546 nm and the polynomial at 400/550/700/850 nm, epoch to
  epoch, per instrument and channel;
- G0(c) centroid offsets (``srf_table.csv`` ``offset``) and the
  wavelength-scale shift/stretch (``phase2_wavescale.csv``), and the
  Phase 0d dispersion trend (``phase0d_trend.csv``).

It prints the tables; nothing is written.

Run from the repository root:  ``python -m hypernet.wiggles.lab_g0_numbers``
"""
import os

import numpy as np
import pandas as pd

from hypernet import labcal
from hypernet.wiggles import DATA_DIR, WIGGLES_DIR

LAMS = (400, 500, 600, 700)


def field_stability():
    st = pd.read_csv(os.path.join(WIGGLES_DIR, 'phase0c_stability.csv'))
    a = st[['system_id', 'cal_period', 'channel', 'n', 'max_std', 'passes_0p2']]
    tr = pd.read_csv(os.path.join(WIGGLES_DIR, 'phase0c_stability_trends.csv'))
    # ranges of the covariates in the 224 requested sequences, per group
    sky = pd.read_csv(os.path.join(WIGGLES_DIR, 'phase0_sky_index.csv'), dtype={'sequence_time': str})
    m = pd.read_parquet(os.path.join(os.environ['OS_COLOR'], 'hypernet', 'wiggles', 'phase0',
                                     'release2_models.parquet'))
    g = m.drop_duplicates(['system_id', 'cal_period', 'sequence_time']).groupby(
        ['system_id', 'cal_period'])
    rng = g.agg(sza_lo=('sza_l1a', 'min'), sza_hi=('sza_l1a', 'max'),
                month_n=('month', 'nunique')).reset_index()
    sky['site_code'] = sky['site'].str[:-2]
    mm = m.drop_duplicates(['system_id', 'cal_period', 'sequence_time']).merge(
        sky[['site_code', 'sequence_time', 'ld_ed_750']], on=['site_code', 'sequence_time'],
        how='left')
    rs = mm.groupby(['system_id', 'cal_period'])['ld_ed_750'].agg(['min', 'max']).reset_index()
    tr = tr.merge(rng, on=['system_id', 'cal_period']).merge(rs, on=['system_id', 'cal_period'])
    tr['d_sza'] = tr['slope_sza'] * (tr['sza_hi'] - tr['sza_lo'])
    tr['d_sky'] = tr['slope_sky'] * (tr['max'] - tr['min'])
    tr['sig_sza'] = (tr['slope_sza'] / tr['slope_sza_err']).abs()
    tr['sig_sky'] = (tr['slope_sky'] / tr['slope_sky_err']).abs()
    pr = pd.read_csv(os.path.join(WIGGLES_DIR, 'phase0c_stability_pairs.csv'))
    return a, tr, pr


def lab_stability():
    s = labcal.load_lab_srf()
    s = s[~s['epoch'].str.contains('tst')]
    rows = []
    for (sn, ch), g in s.groupby(['sn', 'channel']):
        eps = sorted(g['epoch'].unique())
        if len(eps) < 2:
            continue
        p = g[g['lam_ref'].isin([404.6563, 435.8328, 546.0735])].pivot_table(
            index='lam_ref', columns='epoch', values='lab_fwhm')
        dmax = float((p.max(axis=1) - p.min(axis=1)).max())
        cf = {e: labcal.coefs(sn, ch, e) for e in eps}
        px = labcal.px_of(cf[eps[0]], np.array([400.0, 550.0, 700.0, 850.0]))
        lam = np.array([labcal.poly(cf[e], px) for e in eps])
        rows.append(dict(sn=sn, channel=ch, epochs=' '.join(eps),
                         fwhm_range_max=dmax,
                         poly_range_max=float((lam.max(axis=0) - lam.min(axis=0)).max())))
    return pd.DataFrame(rows)


def offsets():
    t = pd.read_csv(os.path.join(DATA_DIR, 'srf_table.csv'))
    ws = pd.read_csv(os.path.join(WIGGLES_DIR, 'phase2_wavescale.csv'))
    tr = pd.read_csv(os.path.join(WIGGLES_DIR, 'phase0d_trend.csv'))
    return (t[['instrument', 'channel', 'cal_period', 'source', 'offset', 'offset_err']],
            ws[['source', 'instrument', 'channel', 'shift_med', 'stretch_med', 'max_err_med',
                'frac_ok']], tr)


def main():
    pd.set_option('display.width', 220)
    pd.set_option('display.max_columns', 30)
    a, tr, pr = field_stability()
    print('G0(b) field: largest FWHM std over 400-700 nm (nm):')
    print(a.to_string(index=False, float_format='%.3f'))
    print('  overall max std %.3f nm; all pass 0.2: %s' % (a['max_std'].max(), a['passes_0p2'].all()))
    q = tr[tr['quantity'].str.startswith('fwhm_600')]
    print('\nG0(b) field trends of FWHM at 600 nm: change over each group\'s SZA and sky range (nm):')
    print(q[['system_id', 'cal_period', 'quantity', 'n', 'sza_lo', 'sza_hi', 'd_sza', 'sig_sza',
             'd_sky', 'sig_sky', 'slope_season', 'slope_season_err']]
          .to_string(index=False, float_format='%.3f'))
    print('  max |d_sza| %.3f, max |d_sky| %.3f nm' % (q['d_sza'].abs().max(), q['d_sky'].abs().max()))
    print('\nG0(b) pairs (recalibration, site):')
    print(pr[pr['quantity'].str.startswith('fwhm_')][['pair', 'quantity', 'diff', 'err', 'z']]
          .to_string(index=False, float_format='%.3f'))
    print('  max |diff| over FWHM pairs %.3f nm'
          % pr[pr['quantity'].str.startswith('fwhm_')]['diff'].abs().max())
    ls = lab_stability()
    print('\nG0(b) lab epochs: max range of FWHM (405/436/546) and of the polynomial (nm):')
    print(ls.to_string(index=False, float_format='%.3f'))
    o, ws, t8d = offsets()
    print('\nG0(c) centroid offsets in srf_table (nm, vs HSRS air):')
    print(o.to_string(index=False, float_format='%+.3f'))
    print('  max |offset| %.3f nm' % o['offset'].abs().max())
    print('\nG0(c) wavelength scale (phase2_wavescale.csv, shift at 630 nm, stretch per 240 nm):')
    print(ws.to_string(index=False, float_format='%+.3f'))
    print('\nPhase 0d dispersion trend:')
    print(t8d.to_string(index=False, float_format='%+.3f'))


if __name__ == '__main__':
    main()

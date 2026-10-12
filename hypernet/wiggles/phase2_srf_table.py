"""Phase 2, task 6: write the per-instrument SRF table ``hypernet/data/srf_table.csv``.

Sources (until the L1A delivery is fitted, every row is ``provisional``):

- ``release2_l2b``: ``hypernet/data/release2_srf_models.json``, the clear-sky
  median template fits of the 224 Release 2 L2B sequences, per instrument x
  IRR calibration (Phase 1 task 7).  E comes from L2 irradiance (resampled
  onto the L grid; a proxy E grid except 122304 @ 2024-11-12).
- ``veit_l1a``: ``hypernet/data/veit_srf_model.json``, the L1A template fit of
  Kevin's VEIT sequence (122304, IRR 2024-11-12 / RAD 2024-11-08), preferred
  over the Release 2 row of that period.

Per row:

- the FWHM(lambda) quadratic, its covariance, the centroid offset and error;
- ``valid_from``: the channel's own calibration date (IRR for E, RAD for Ld);
- ``delta_fwhm`` = median FWHM_Ld - FWHM_E over 400-680 nm
  (``edinterp.delta_fwhm``) and ``method`` = gate G1's rule (>= 0.15 nm:
  srf);
- for L2-derived pairs, ``delta_fwhm_l1a_est``: delta_fwhm corrected by the
  L2-vs-L1A width bias measured on VEIT (``phase0c_validation.csv``, the "E
  via E grid" route), and ``method_l1a_est``.  ``method_sensitive`` marks rows
  where the two methods differ;
- the wavelength-scale correction: the calibration-period median shift and
  stretch (common window pattern removed) of
  ``$OS_COLOR/hypernet/wiggles/phase2/wavescale_per_sequence.parquet`` (task
  4b), with errors 1.2533 sd / sqrt(n); the VEIT rows use its one sequence.

Run from the repository root: ``python -m hypernet.wiggles.phase2_srf_table``.
"""
import json
import os

import numpy as np
import pandas as pd

from hypernet import edinterp, srf
from hypernet.wiggles import DATA_DIR, WIGGLES_DIR

ROOT = os.path.join(os.getenv('OS_COLOR', '.'), 'hypernet', 'wiggles', 'phase2')
VEIT_CAL = '2024-11-08/2024-11-12'
LAMS = (400.0, 550.0, 700.0)


def l2_bias():
    """Median over 400-680 nm of (L2 - L1A) width bias, Ld minus E: the amount
    by which an L2-derived delta_fwhm is low (one VEIT sequence)."""
    v = pd.read_csv(os.path.join(WIGGLES_DIR, 'phase0c_validation.csv'))
    v = v[v['route'] == 'L2 route, E via E grid']
    lam = np.arange(400.0, 681.0, 1.0)
    d = {ch: np.interp(lam, g['lam'], g['diff']) for ch, g in v.groupby('channel')}
    return float(np.median(d['E'] - d['Ld']))


def _row(m, source, cal_period, provisional, e_grid, n_seq, note=''):
    rad, irr = cal_period.split('/')
    cov = np.asarray(m.cov, float)
    r = dict(instrument=m.instrument, channel=m.channel, cal_period=cal_period,
             cal_date_rad=rad, cal_date_irr=irr, valid_from=irr if m.channel == 'E' else rad,
             source=source, provisional=provisional, e_grid=e_grid if m.channel == 'E' else 'own',
             n_sequences=int(n_seq), c0=m.coeffs[0], c1=m.coeffs[1], c2=m.coeffs[2],
             cov00=cov[0, 0], cov01=cov[0, 1], cov02=cov[0, 2], cov11=cov[1, 1],
             cov12=cov[1, 2], cov22=cov[2, 2], lam_ref=m.lam_ref, lam_scale=m.lam_scale,
             lam_min=m.lam_min, lam_max=m.lam_max, chi2_nu=m.chi2_nu, offset=m.offset,
             offset_err=m.offset_err, note=note)
    for lam in LAMS:
        r['fwhm_%d' % lam] = m.fwhm(lam)
        r['fwhm_err_%d' % lam] = m.fwhm_err(lam)
    return r


def srf_rows():
    rows = []
    doc = json.load(open(os.path.join(DATA_DIR, 'release2_srf_models.json')))
    for key, chans in doc['models'].items():
        if '@' not in key:
            continue                                  # the plain keys alias the latest period
        ms = {ch: srf.SRFModel.from_dict(d) for ch, d in chans.items()}
        meta = chans['E']['meta']
        pair = []
        for ch in ('E', 'Ld'):
            mm = chans[ch]['meta']
            pair.append(_row(ms[ch], 'release2_l2b', mm['cal_period'], True, meta['e_grid'],
                             mm['n_sequences'], note='clear-sky median of L2B template fits'))
        rows.append((pair, ms))
    veit = srf.load_srf_models(os.path.join(DATA_DIR, 'veit_srf_model.json'))
    ms = {'E': veit['E'], 'Ld': veit['Ld']}
    for m in ms.values():
        m.instrument = m.instrument or 'HYPSTAR_122304'
    pair = [_row(ms[ch], 'veit_l1a', VEIT_CAL, True, 'own', 1,
                 note='one L1A sequence (VEIT 2026-06-04 08:45)') for ch in ('E', 'Ld')]
    rows.append((pair, ms))
    return rows


def add_method(rows, bias):
    out = []
    for pair, ms in rows:
        d = edinterp.delta_fwhm(ms['E'], ms['Ld'])
        meth = edinterp.choose_method(ms['E'], ms['Ld'])
        l2 = pair[0]['source'] == 'release2_l2b'
        d_est = d + bias if l2 else d
        meth_est = ('srf' if d_est >= edinterp.SRF_THRESHOLD - 1e-9 else 'linear') if l2 else meth
        for r in pair:
            r.update(delta_fwhm=d, method=meth, delta_fwhm_l1a_est=d_est,
                     method_l1a_est=meth_est, method_sensitive=meth != meth_est)
            out.append(r)
    return pd.DataFrame(out)


def add_scale(t):
    p = pd.read_parquet(os.path.join(ROOT, 'wavescale_per_sequence.parquet'))
    p = p[np.isfinite(p['shift'].astype(float))].copy()
    p.loc[p['source'] == 'VEIT L1A', 'cal_period'] = VEIT_CAL
    p['src'] = np.where(p['source'] == 'VEIT L1A', 'veit_l1a', 'release2_l2b')
    sem = lambda v: 1.2533 * v.std() / np.sqrt(len(v)) if len(v) > 1 else np.nan  # noqa: E731
    g = p.groupby(['instrument', 'channel', 'cal_period', 'src'])
    s = g.agg(scale_shift=('shift', 'median'), scale_stretch=('stretch', 'median'),
              scale_shift_err=('shift', sem), scale_stretch_err=('stretch', sem),
              scale_n=('shift', 'size'), scale_frac_ok=('ok', 'mean'),
              _se1=('shift_err', 'median'), _se2=('stretch_err', 'median')).reset_index()
    # one sequence: its own fit errors
    one = s['scale_n'] == 1
    s.loc[one, 'scale_shift_err'] = s.loc[one, '_se1']
    s.loc[one, 'scale_stretch_err'] = s.loc[one, '_se2']
    s = s.drop(columns=['_se1', '_se2']).rename(columns={'src': 'source'})
    return t.merge(s, on=['instrument', 'channel', 'cal_period', 'source'], how='left')


def main():
    bias = l2_bias()
    t = add_scale(add_method(srf_rows(), bias))
    t = t.sort_values(['instrument', 'valid_from', 'channel', 'source']).reset_index(drop=True)
    path = os.path.join(DATA_DIR, 'srf_table.csv')
    t.to_csv(path, index=False, float_format='%.6g')
    pd.set_option('display.width', 250)
    print('L2-vs-L1A delta_fwhm bias (VEIT, 400-680 nm): +%.3f nm' % bias)
    cols = ['instrument', 'channel', 'cal_period', 'source', 'n_sequences', 'fwhm_550',
            'delta_fwhm', 'method', 'delta_fwhm_l1a_est', 'method_l1a_est',
            'method_sensitive', 'scale_shift', 'scale_stretch', 'scale_shift_err',
            'scale_stretch_err', 'scale_n']
    print(t[cols].to_string(index=False, float_format='%.3f'))
    print('wrote', path)
    return t


if __name__ == '__main__':
    main()

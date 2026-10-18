"""Phase 2, task 6: write the per-instrument SRF table ``hypernet/data/srf_table.csv``.

Sources (until the L1A delivery is fitted, every row is ``provisional``):

- ``release2_l2b``: ``hypernet/data/release2_srf_models.json``, the clear-sky
  median template fits of the 224 Release 2 L2B sequences, per instrument x
  IRR calibration (Phase 1 task 7).  E comes from L2 irradiance (resampled
  onto the L grid), fitted through each instrument's own E grid
  (``labcal.grid_for``; lab data prompts task 3).
- ``veit_l1a``: ``hypernet/data/veit_srf_model.json``, the L1A template fit of
  Kevin's VEIT sequence (122304, IRR 2024-11-12 / RAD 2024-11-08), preferred
  over the Release 2 row of that period.

**SRF shape** (lab data prompts task 3; Setup Q6-Q8).  ``shape``, ``shape_p``,
``shape_source``: E is ``supergauss`` on the instruments whose lab E lines
prefer it (``labcal.srf_shape_for``: chi2_sg < 0.5 chi2_g on the processor
epoch), with p the lab median at the period's lab epoch; every other row is
``gauss``.  A super-Gaussian row's FWHM(lambda) is the fit *with that shape*
(p fixed at the lab value): the ``E_sg`` models of ``release2_srf_models.json``
and ``veit_srf_model_lab_shape.json`` (``hypernet.wiggles.lab_shape_fits``).
``shape_p_field`` is the free-p field fit's median (a diagnostic) and
``fwhm_gauss_550`` the Gaussian template FWHM at 550 nm for comparison.
``srftable.model_from_row`` / ``srf_for`` return the shaped SRFModel, so
``edinterp.interpolate_ed_to_l`` convolves with it.

**Gate G1 on a shaped row.**  ``delta_fwhm`` (and so ``method``) uses the
*Gaussian-equivalent* E width, the Gaussian template fit (``fwhm_gauss_550``):
G1's 0.15 nm threshold was set with Gaussian SRFs, and the lab's own FWHMs
behave like Gaussian-fit widths too.  The FWHM of a flat-topped SRF is
0.25-0.30 nm larger than its Gaussian-equivalent width, so a Delta FWHM from
it (``delta_fwhm_shaped``, a diagnostic) is not comparable with the
threshold: it would flip 122304 @ 2024-11 to ``linear``.

**Lab validation.**  ``lab_epoch`` (``labcal.physical_epoch`` at the
channel's calibration date) and ``lab_minus_table_{405,436,546}``: the lab
width of that channel minus the row's FWHM(lambda) at Hg 404.66 (a blend),
435.83 and 546.07 nm (``hypernet/data/lab_srf.csv``), like with like
(``lab_fwhm_def``): for a Gaussian row the lab's FWHM (``lab``), for a
super-Gaussian row our super-Gaussian refit of the lamp line (``lab_sg_refit``,
``fwhm_sg``), since the lab's FWHM reads 0.1-0.15 nm below that on the
flat-topped E lines.  ``method_lab`` is gate
G1's rule on the lab Delta FWHM; where it disagrees with ``method`` the row is
reported, and for 122303 (lab data prompts task 3) the lab method is used with
a note.

Per row:

- the FWHM(lambda) quadratic, its covariance, the centroid offset and error;
- ``valid_from``: the channel's own calibration date (IRR for E, RAD for Ld);
- ``delta_fwhm`` = median FWHM_Ld - FWHM_E over 400-680 nm
  (``edinterp.delta_fwhm``) and ``method`` = gate G1's rule (>= 0.15 nm:
  srf);
- for L2-derived pairs, ``delta_fwhm_l1a_est``: delta_fwhm corrected by the
  L2-vs-L1A width bias measured on VEIT (``phase0c_validation.csv``, the "E
  via E grid" route), and ``method_l1a_est``.  ``method_sensitive`` marks rows
  where the two methods differ.  The lamp does not support that bias: against
  the lab FWHMs, L2 E (own grid) is no wider than Ld once the lab's FWHM
  definition is allowed for (``hypernet/wiggles/lab_l2_bias.py``,
  ``lab_l2_bias.csv``); the column is kept as the one-sequence L1A estimate;
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

from hypernet import edinterp, labcal, srf
from hypernet.wiggles import DATA_DIR, WIGGLES_DIR

ROOT = os.path.join(os.getenv('OS_COLOR', '.'), 'hypernet', 'wiggles', 'phase2')
VEIT_CAL = '2024-11-08/2024-11-12'
LAMS = (400.0, 550.0, 700.0)
LAB_LINES = {'405': 404.6563, '436': 435.8328, '546': 546.0735}
#: Instruments whose method is taken from the lab if the data disagree
#: (lab data prompts task 3).
LAB_METHOD = ('HYPSTAR_122303',)


def l2_bias():
    """Median over 400-680 nm of (L2 - L1A) width bias, Ld minus E: the amount
    by which an L2-derived delta_fwhm is low (one VEIT sequence)."""
    v = pd.read_csv(os.path.join(WIGGLES_DIR, 'phase0c_validation.csv'))
    v = v[v['route'] == 'L2 route, E via E grid']
    lam = np.arange(400.0, 681.0, 1.0)
    d = {ch: np.interp(lam, g['lam'], g['diff']) for ch, g in v.groupby('channel')}
    return float(np.median(d['E'] - d['Ld']))


def _row(m, source, cal_period, provisional, e_grid, n_seq, note='', gauss=None):
    rad, irr = cal_period.split('/')
    cov = np.asarray(m.cov, float)
    sh = labcal.srf_shape_for(m.instrument, m.channel, irr if m.channel == 'E' else rad)
    if sh['shape'] != m.shape or (m.shape != 'gauss' and abs(sh['shape_p'] - m.shape_par) > 1e-9):
        raise ValueError('%s %s %s: model shape %s %s, lab says %s %s' % (
            m.instrument, m.channel, cal_period, m.shape, m.shape_par, sh['shape'],
            sh['shape_p']))
    r = dict(instrument=m.instrument, channel=m.channel, cal_period=cal_period,
             cal_date_rad=rad, cal_date_irr=irr, valid_from=irr if m.channel == 'E' else rad,
             source=source, provisional=provisional, e_grid=e_grid if m.channel == 'E' else 'own',
             n_sequences=int(n_seq), c0=m.coeffs[0], c1=m.coeffs[1], c2=m.coeffs[2],
             cov00=cov[0, 0], cov01=cov[0, 1], cov02=cov[0, 2], cov11=cov[1, 1],
             cov12=cov[1, 2], cov22=cov[2, 2], lam_ref=m.lam_ref, lam_scale=m.lam_scale,
             lam_min=m.lam_min, lam_max=m.lam_max, chi2_nu=m.chi2_nu, offset=m.offset,
             offset_err=m.offset_err, note=note, shape=m.shape,
             shape_p=m.shape_par, shape_source=sh['shape_source'],
             shape_p_field=m.meta.get('shape_p_field', np.nan) if m.shape != 'gauss' else np.nan,
             fwhm_gauss_550=(gauss or m).fwhm(550.0), lab_epoch=sh['lab_epoch'])
    for lam in LAMS:
        r['fwhm_%d' % lam] = m.fwhm(lam)
        r['fwhm_err_%d' % lam] = m.fwhm_err(lam)
    lab = labcal.load_lab_srf()
    g = lab[(lab['sn'] == m.instrument[-6:]) & (lab['epoch'] == sh['lab_epoch']) &
            (lab['channel'] == ('E' if m.channel == 'E' else 'L'))].set_index('lam_ref')
    col = 'fwhm_sg' if m.shape == 'supergauss' else 'lab_fwhm'
    r['lab_fwhm_def'] = 'lab_sg_refit' if m.shape == 'supergauss' else 'lab'
    for k, lam in LAB_LINES.items():
        r['lab_minus_table_' + k] = float(g.loc[lam, col]) - float(m.fwhm(lam))
    return r


def _e_model(chans):
    """The E model of a Release 2 key: the lab-shape fit where there is one."""
    if 'E_sg' in chans:
        return srf.SRFModel.from_dict(chans['E_sg']), srf.SRFModel.from_dict(chans['E'])
    return srf.SRFModel.from_dict(chans['E']), None


def srf_rows():
    rows = []
    doc = json.load(open(os.path.join(DATA_DIR, 'release2_srf_models.json')))
    for key, chans in doc['models'].items():
        if '@' not in key:
            continue                                  # the plain keys alias the latest period
        e, eg = _e_model(chans)
        ms = {'E': e, 'Ld': srf.SRFModel.from_dict(chans['Ld']), 'E_gauss': eg or e}
        pair = []
        for ch in ('E', 'Ld'):
            mm = ms[ch].meta
            note = 'clear-sky median of L2B template fits'
            if ch == 'E' and eg is not None:
                note += ', super-Gaussian with p fixed at the lab value'
            pair.append(_row(ms[ch], 'release2_l2b', mm['cal_period'], True,
                             chans['E']['meta']['e_grid'], mm['n_sequences'], note=note,
                             gauss=eg if ch == 'E' else None))
        rows.append((pair, ms))
    veit = srf.load_srf_models(os.path.join(DATA_DIR, 'veit_srf_model.json'))
    ms = {'E': veit['E'], 'Ld': veit['Ld']}
    eg = None
    if labcal.srf_shape_for('122304', 'E', VEIT_CAL.split('/')[1])['shape'] != 'gauss':
        eg = ms['E']
        ms['E'] = srf.load_srf_models(os.path.join(DATA_DIR, 'veit_srf_model_lab_shape.json'))['E']
    ms['E_gauss'] = eg or ms['E']
    for m in ms.values():
        m.instrument = m.instrument or 'HYPSTAR_122304'
    note = 'one L1A sequence (VEIT 2026-06-04 08:45)'
    pair = [_row(ms[ch], 'veit_l1a', VEIT_CAL, True, 'own', 1,
                 note=note + (', super-Gaussian with p fixed at the lab value'
                              if ch == 'E' and eg is not None else ''),
                 gauss=eg if ch == 'E' else None) for ch in ('E', 'Ld')]
    rows.append((pair, ms))
    return rows


def lab_dfwhm(sn, epoch):
    """Lab Delta FWHM (median L - E over the 400-680 nm lines with S/N >= 10)."""
    s = labcal.load_lab_srf()
    g = s[(s['sn'] == sn) & (s['epoch'] == epoch)]
    p = g.pivot_table(index='lam_ref', columns='channel', values='lab_fwhm')
    weak = g.groupby('lam_ref')['low_snr'].any().reindex(p.index).values
    m = (p.index >= edinterp.DFWHM_BAND[0]) & (p.index <= edinterp.DFWHM_BAND[1]) & ~weak
    return float((p['L'] - p['E'])[m].median())


def add_lab_method(t):
    """Gate G1's rule on the lab Delta FWHM; for LAB_METHOD instruments the
    lab method replaces a disagreeing data method, with a note."""
    t = t.copy()
    t['lab_dfwhm'] = [lab_dfwhm(i[-6:], e) for i, e in zip(t['instrument'], t['lab_epoch'])]
    t['method_lab'] = np.where(t['lab_dfwhm'] >= edinterp.SRF_THRESHOLD - 1e-9, 'srf', 'linear')
    t['method_from'] = 'data'
    bad = t['method'] != t['method_lab']
    for k in t.index[bad]:
        print('METHOD DISAGREES: %s %s %s %s: data %s (delta_fwhm %.3f), lab %s (%.3f)'
              % tuple(t.loc[k, ['instrument', 'channel', 'cal_period', 'source', 'method',
                                'delta_fwhm', 'method_lab', 'lab_dfwhm']]))
    fix = bad & t['instrument'].isin(LAB_METHOD)
    if fix.any():
        t.loc[fix, 'note'] = t.loc[fix, 'note'] + '; method set from the lab (data gave ' + \
            t.loc[fix, 'method'] + ')'
        t.loc[fix, 'method'] = t.loc[fix, 'method_lab']
        t.loc[fix, 'method_from'] = 'lab'
        t['method_sensitive'] = t['method'] != t['method_l1a_est']
    return t


def add_method(rows, bias):
    out = []
    for pair, ms in rows:
        d = edinterp.delta_fwhm(ms['E_gauss'], ms['Ld'])          # Gaussian-equivalent
        meth = edinterp.choose_method(ms['E_gauss'], ms['Ld'])
        d_shaped = edinterp.delta_fwhm(ms['E'], ms['Ld'])
        l2 = pair[0]['source'] == 'release2_l2b'
        d_est = d + bias if l2 else d
        meth_est = ('srf' if d_est >= edinterp.SRF_THRESHOLD - 1e-9 else 'linear') if l2 else meth
        for r in pair:
            r.update(delta_fwhm=d, delta_fwhm_shaped=d_shaped, method=meth,
                     delta_fwhm_l1a_est=d_est,
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
    t = add_lab_method(add_scale(add_method(srf_rows(), bias)))
    t = t.sort_values(['instrument', 'valid_from', 'channel', 'source']).reset_index(drop=True)
    path = os.path.join(DATA_DIR, 'srf_table.csv')
    t.to_csv(path, index=False, float_format='%.6g')
    pd.set_option('display.width', 250)
    print('L2-vs-L1A delta_fwhm bias (VEIT, 400-680 nm): +%.3f nm' % bias)
    cols = ['instrument', 'channel', 'cal_period', 'source', 'n_sequences', 'shape', 'shape_p',
            'shape_p_field', 'fwhm_550', 'fwhm_gauss_550', 'delta_fwhm', 'delta_fwhm_shaped',
            'method', 'method_lab',
            'method_from', 'delta_fwhm_l1a_est', 'method_l1a_est', 'method_sensitive',
            'lab_minus_table_436', 'lab_minus_table_546', 'scale_shift', 'scale_stretch',
            'scale_n']
    print(t[cols].to_string(index=False, float_format='%.3f'))
    print('wrote', path)
    return t


if __name__ == '__main__':
    main()

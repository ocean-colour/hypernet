"""Phase 2, task 4b: the wavelength-scale error (shift + stretch) of E and Ld.

``hypernet.srfcal.fit_wavelength_scale`` on:

- **VEIT** (Kevin's L1A sample, HYPSTAR 122304): the scan-mean E and Ld on
  their native grids, with the flattened scan-scatter errors, and the SRF
  widths of ``hypernet/data/veit_srf_model.json``;
- **Release 2** (the 224 L2B files of ``docs/wiggles_data_request.csv``):
  Ld on its native grid, and E through the E-grid resampling (``via_grid`` =
  the instrument's own E grid, ``labcal.grid_for``, as in Phase 0c since lab
  data prompts task 3; before, the 122304 E grid was a proxy).  The SRF
  widths are the Gaussian ones of ``release2_srf_models.json`` (channel
  ``E``, not ``E_sg``) for the instrument and IRR calibration: the window
  shifts are fitted with a Gaussian kernel.  L2B has no per-pixel errors, so the
  window fits are unweighted (their errors come from the residuals).

**The common window pattern.**  The per-window offsets scatter ~0.13 nm
about each straight line, far beyond their errors, and the scatter is the
*same* on all seven instruments and both channels (``window_pattern``:
correlation ~0.8).  Spectrometers with different pixel grids cannot share a
pixel dispersion error, so this is a property of the reference or the method
(HSRS line positions/shapes in those windows, or a common SRF asymmetry).
Its median over the 14 instrument x channel groups is written to
``hypernet/data/wavescale_window_pattern.csv`` and subtracted before the
straight-line fit (``pattern=``); every fit is also done without it
(``*_raw`` columns).

The error convention is twin case iv's: a pixel reported at lambda really
samples lambda + shift + stretch (lambda - 630) / 240.

Outputs:

- ``$OS_COLOR/hypernet/wiggles/phase2/wavescale_per_sequence.parquet``;
- ``hypernet/wiggles/phase2_wavescale.csv``: per instrument x channel
  (Release 2) and VEIT, the medians of shift, stretch and their errors, the
  fraction that passes Setup Q8's 0.03 nm, the sequence-to-sequence scatter,
  and Phase 0 8d's slope for comparison;
- ``hypernet/data/wavescale_window_pattern.csv`` (window centre, the common
  pattern, its scatter between groups);
- ``hypernet/wiggles/phase2_wavescale_windows.csv``: the per-window pattern
  beyond the straight line (``window_pattern``), from the Phase 0c Release 2
  template fits, and how it correlates between instruments and channels.

Run from the repository root:
``python -m hypernet.wiggles.phase2_wavescale [--workers 12] [--limit N]``.
"""
import argparse
import os
import time

import numpy as np
import pandas as pd

from hypernet import labcal, srf, srfcal, twin, whn_srf
from hypernet import whn_l1a as wl
from hypernet.wiggles import DATA_DIR, OUT, REPO, WIGGLES_DIR

ROOT = os.path.join(os.getenv('OS_COLOR', '.'), 'hypernet', 'wiggles', 'phase2')
REQUEST = os.path.join(REPO, 'docs', 'wiggles_data_request.csv')
VEIT = dict(site='VEIT', seq_time='20260604T0845')
KEYS = ('shift', 'stretch', 'shift_err', 'stretch_err', 'max_err', 'slope_per100nm',
        'chi2_nu', 'n_windows', 'ok')


def _row(r, raw, **meta):
    d = {k: r[k] for k in KEYS}
    d['flags'] = ','.join(r['flags'])
    d.update({k + '_raw': raw[k] for k in KEYS})
    return dict(meta, **d)


def _fit2(wave, spec, mod, pattern, **kw):
    """The fit with the common pattern removed, and without (raw)."""
    return (srfcal.fit_wavelength_scale(wave, spec, srf=mod, pattern=pattern, **kw),
            srfcal.fit_wavelength_scale(wave, spec, srf=mod, **kw))


def veit(pattern):
    f = wl.sequence_files(**VEIT)
    irr, rad = wl.load_l1a_irr(f['L1A_IRR']), wl.load_l1a_rad(f['L1A_RAD'])
    models = srf.load_srf_models(os.path.join(DATA_DIR, 'veit_srf_model.json'))
    rows = []
    for ch, wave, scans in (('E', irr['wave'], irr['scans']),
                            ('Ld', rad['wave'], rad['Ld']['scans'])):
        mean, e_mean, _ = srf.scan_errors(scans, flatten_px=whn_srf.FLATTEN_PX)
        r, raw = _fit2(wave, mean, models[ch], pattern, err=e_mean)
        rows.append(_row(r, raw, source='VEIT L1A', instrument='HYPSTAR_122304', channel=ch,
                         sequence_time=VEIT['seq_time'], site_code='VEIT', cal_period='',
                         e_grid='own'))
    return rows


def _release2_one(args):
    r, pattern = args
    srfs = twin.load_instrument_srfs()
    out = []
    try:
        d = wl.load_l2b(wl.release2_path(r['site'], r['sequence_time'], r['file']))
        m = d['meta']
        sysid = m['system_id']
        irr_date = m['instrument_calibration_date_irr']
        cal = '%s/%s' % (m['instrument_calibration_date_rad'], irr_date)
        mod = srfs.get('%s@%s' % (sysid, irr_date), srfs.get(sysid))
        gE = labcal.grid_for(sysid, 'E')
        for ch, spec, via in (('E', d['irradiance'][:, 0], gE),
                              ('Ld', d['downwelling_radiance'][:, 0], None)):
            res, raw = _fit2(d['wave'], spec, mod[ch], pattern, via_grid=via)
            out.append(_row(res, raw, source='Release 2 L2B', instrument=sysid, channel=ch,
                            sequence_time=r['sequence_time'], site_code=r['site'][:-2],
                            cal_period=cal, sza=float(np.nanmean(d['sza'])),
                            e_grid='own'))
    except Exception as ex:
        out.append(dict(source='Release 2 L2B', sequence_time=r['sequence_time'],
                        site_code=r['site'][:-2], flags='error: %s' % ex, ok=False))
    return out


def release2(workers, pattern, limit=None):
    from multiprocessing import Pool
    req = pd.read_csv(REQUEST, dtype=str)
    if limit:
        req = req.head(limit)
    rows = []
    with Pool(workers) as pool:
        for o in pool.imap_unordered(_release2_one, [(r, pattern) for _, r in req.iterrows()]):
            rows += o
    return rows


def summarise(df):
    d = df[np.isfinite(df['shift'].astype(float))]
    g = d.groupby(['source', 'instrument', 'channel'])
    s = g.agg(n_seq=('sequence_time', 'nunique'), shift_med=('shift', 'median'),
              shift_scatter=('shift', 'std'), stretch_med=('stretch', 'median'),
              stretch_scatter=('stretch', 'std'), shift_err_med=('shift_err', 'median'),
              stretch_err_med=('stretch_err', 'median'), max_err_med=('max_err', 'median'),
              frac_ok=('ok', 'mean'), slope_per100nm_med=('slope_per100nm', 'median'),
              shift_raw_med=('shift_raw', 'median'), stretch_raw_med=('stretch_raw', 'median'),
              max_err_raw_med=('max_err_raw', 'median'), frac_ok_raw=('ok_raw', 'mean'),
              chi2_med=('chi2_nu', 'median'), chi2_raw_med=('chi2_nu_raw', 'median'),
              slope_raw_med=('slope_per100nm_raw', 'median'),
              e_grid=('e_grid', 'first')).reset_index()
    p8d = pd.read_csv(os.path.join(WIGGLES_DIR, 'phase0d_trend.csv'))
    p8d = p8d.rename(columns={'slope_mean': 'slope_per100nm_8d', 'slope_sem': 'slope_sem_8d'})
    s = s.merge(p8d[['instrument', 'channel', 'slope_per100nm_8d', 'slope_sem_8d']],
                on=['instrument', 'channel'], how='left')
    # 8d's dlam slope as this convention's stretch: stretch = -slope * HALF / 100
    s['stretch_8d'] = -s['slope_per100nm_8d'] * srfcal.HALF / 100.0
    return s


def write_pattern(pat):
    """The median over instrument x channel groups of the per-window pattern."""
    g = pat.groupby('lam_center')['resid_med']
    out = pd.DataFrame(dict(lam_center=g.median().index, pattern=g.median().values,
                            scatter=g.std().values, n_groups=g.size().values))
    out.to_csv(os.path.join(DATA_DIR, 'wavescale_window_pattern.csv'), index=False,
               float_format='%.4f')
    return dict(zip(out['lam_center'], out['pattern']))


def window_pattern():
    """Per-window residual of -dlam about each sequence's straight line,
    median over the sequences of each instrument x channel (Phase 0c Release 2
    template fits, SRF width free).  Returns the long table and the matrix of
    correlations between the instrument x channel patterns."""
    t = pd.read_parquet(os.path.join(OUT, 'release2_template.parquet'))
    t = t[t['ok'] & np.isfinite(t['dlam_err']) & (t['dlam_err'] > 0)]
    rows = []
    for (inst, ch, seq), g in t.groupby(['system_id', 'channel', 'sequence_time']):
        if len(g) < 8:
            continue
        x = (g['lam_center'].values - srfcal.LAM0) / srfcal.HALF
        y, w = -g['dlam'].values, 1.0 / g['dlam_err'].values ** 2
        A = np.column_stack([np.ones_like(x), x])
        coef = np.linalg.solve(A.T @ (A * w[:, None]), A.T @ (w * y))
        rows.append(pd.DataFrame(dict(instrument=inst, channel=ch, sequence_time=seq,
                                      lam_center=g['lam_center'].values, resid=y - A @ coef)))
    r = pd.concat(rows, ignore_index=True)
    pat = (r.groupby(['instrument', 'channel', 'lam_center'])['resid']
           .agg(resid_med='median', resid_sem=lambda v: 1.2533 * v.std() / np.sqrt(len(v)),
                n='size').reset_index())
    wide = pat.pivot_table(index='lam_center', columns=['instrument', 'channel'],
                           values='resid_med')
    return pat, wide.corr()


def main(argv=None):
    ap = argparse.ArgumentParser(description='Phase 2 task 4b: wavelength-scale fits')
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--limit', type=int, default=None, help='first N Release 2 rows')
    a = ap.parse_args(argv)
    t0 = time.time()
    pat, corr = window_pattern()
    pattern = write_pattern(pat)
    rows = veit(pattern) + release2(a.workers, pattern, a.limit)
    df = pd.DataFrame(rows)
    os.makedirs(ROOT, exist_ok=True)
    df.to_parquet(os.path.join(ROOT, 'wavescale_per_sequence.parquet'))
    s = summarise(df)
    s.to_csv(os.path.join(WIGGLES_DIR, 'phase2_wavescale.csv'), index=False,
             float_format='%.4f')
    pat.to_csv(os.path.join(WIGGLES_DIR, 'phase2_wavescale_windows.csv'), index=False,
               float_format='%.4f')
    pd.set_option('display.width', 250)
    wide = pat.pivot_table(index='lam_center', columns=['instrument', 'channel'],
                           values='resid_med')
    print('per-window residual about the straight line (median over sequences), nm:')
    print(wide.to_string(float_format='%+.3f'))
    c = corr.values[np.triu_indices_from(corr.values, 1)]
    print('pattern correlation between the %d instrument x channel groups: median %.2f, '
          'range %.2f to %.2f; rms of the pattern %.3f nm'
          % (corr.shape[0], np.median(c), c.min(), c.max(),
             float(np.sqrt(np.nanmean(wide.values ** 2)))))
    print('%d fits in %.0f s; errors: %d' % (len(df), time.time() - t0,
                                             df['flags'].astype(str).str.startswith('error').sum()))
    cols = ['instrument', 'channel', 'n_seq', 'shift_med', 'stretch_med', 'max_err_med',
            'frac_ok', 'chi2_med', 'shift_raw_med', 'stretch_raw_med', 'max_err_raw_med',
            'frac_ok_raw', 'chi2_raw_med', 'slope_raw_med', 'slope_per100nm_8d',
            'shift_scatter', 'stretch_scatter']
    print(s[['source'] + cols].to_string(index=False, float_format='%.3f'))
    return df, s


if __name__ == '__main__':
    main()

"""Lab wavelength-calibration tables (lab data prompts, task 2).

Reads Joel's lab tree (:func:`hypernet.labcal.lab_root`) and writes the
committed tables

- ``hypernet/data/lab_wl_coefs.csv`` (26 rows: 13 epochs x L/E);
- ``hypernet/data/lab_srf.csv`` (every line of every epoch, both channels);

then checks the grids (Setup Q9) against the processor's:

- VEIT L1A E and L (122304) and one Release 2 L2B L grid per instrument x
  calibration period of ``docs/wiggles_data_request.csv``, with
  :func:`hypernet.labcal.grid_for`, and with the plain crop rule on the
  processor epoch (no :data:`~hypernet.labcal.CROP_EPOCH` exception);
- for 120242, the pixel window selected by each of its epochs' polynomials
  with the 2024_05 wavelengths (the 2020_09-window hypothesis).

Intermediates: ``$OS_COLOR/hypernet/wiggles/lab_data/lab_tables/``
(``grid_check.csv``, ``crop_hypothesis_120242.csv``).

Run from the repository root:  ``python -m hypernet.wiggles.lab_tables``
"""
import os

import numpy as np
import pandas as pd

from hypernet import labcal
from hypernet import whn_l1a as wl

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REQUEST = os.path.join(REPO, 'docs', 'wiggles_data_request.csv')
OUT = os.path.join(labcal.lab_root(), 'lab_tables')


def build_tables():
    coef_rows, srf_rows = [], []
    for sn, ep, d in labcal.epochs():
        e = labcal.read_epoch(d)
        coef_rows += labcal.coef_rows(sn, ep, e)
        srf_rows.append(labcal.fit_epoch_lines(sn, ep, e))
    coefs = pd.DataFrame(coef_rows)
    srf = pd.concat(srf_rows, ignore_index=True)
    print('wrote', labcal.write_wl_coefs(coefs), len(coefs), 'rows')
    print('wrote', labcal.write_lab_srf(srf), len(srf), 'rows')
    return coefs, srf


def _cmp(label, sn, ch, target, grid):
    d = dict(target=label, sn=sn, channel=ch, n_target=len(target), n_grid=len(grid))
    if len(grid) == len(target):
        d.update(max_abs=float(np.abs(grid - target).max()),
                 rms=float(np.sqrt(np.mean((grid - target) ** 2))))
    else:
        n = min(len(grid), len(target))
        d.update(max_abs_common_start=float(np.abs(grid[:n] - target[:n]).max()),
                 first_grid=grid[0], first_target=target[0], last_grid=grid[-1],
                 last_target=target[-1])
    return d


def targets():
    """(label, sn, channel, wavelengths) of the processor grids we hold."""
    out = []
    f = wl.sequence_files(site='VEIT', seq_time='20260604T0845')
    out.append(('VEIT L1A E', '122304', 'E', wl.load_l1a_irr(f['L1A_IRR'])['wave']))
    out.append(('VEIT L1A L', '122304', 'L', wl.load_l1a_rad(f['L1A_RAD'])['wave']))
    req = pd.read_csv(REQUEST, dtype=str)
    for (inst, cal), g in req.groupby(['instrument', 'cal_dates_rad_irr']):
        for _, r in g.iterrows():
            p = wl.release2_path(r['site'], r['sequence_time'], r['file'])
            if os.path.exists(p):
                out.append(('R2 L2B L %s %s' % (cal, r['site']), inst.split('_')[1], 'L',
                            wl.load_l2b(p)['wave']))
                break
    return out


def check_grids():
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for label, sn, ch, w in targets():
        rows.append(dict(_cmp(label, sn, ch, w, labcal.grid_for(sn, ch)), method='grid_for'))
        c = labcal.coefs(sn, ch)
        s, e = labcal.crop_window(c)
        rows.append(dict(_cmp(label, sn, ch, w, labcal.poly(c, np.arange(s, e + 1))),
                         method='plain rule on processor epoch'))
        if sn == '122304':
            g = labcal.grid_for(sn, ch, epoch='2024_09')
            rows.append(dict(_cmp(label, sn, ch, w, g), method='2024_09 (guard)'))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT, 'grid_check.csv'), index=False)
    pd.set_option('display.width', 250)
    print(t.round(6).to_string(index=False))

    # Q9: which epoch's window reproduces 120242's Release 2 L grid?
    w = [x for x in targets() if x[1] == '120242'][0][3]
    cw = labcal.coefs('120242', 'L')
    hyp = []
    for ep in ('2020_09', '2022_05', '2022_06', '2024_05'):
        for ch in labcal.CHANNELS:
            s, e = labcal.crop_window(labcal.coefs('120242', ch, ep))
            d = dict(window_epoch=ep, channel=ch, start=s, end=e, n=e - s + 1)
            if ch == 'L':
                g = labcal.poly(cw, np.arange(s, e + 1))
                d['matches_R2_L'] = len(g) == len(w) and bool(np.abs(g - w).max() < 1e-4)
            hyp.append(d)
    hyp = pd.DataFrame(hyp)
    hyp.to_csv(os.path.join(OUT, 'crop_hypothesis_120242.csv'), index=False)
    print('\n120242: pixel window of each epoch polynomial, 2024_05 wavelengths on it:')
    print(hyp.to_string(index=False))
    print('R2 L grid: %d px, %.4f-%.4f nm' % (len(w), w[0], w[-1]))
    return t, hyp


def main():
    build_tables()
    check_grids()


if __name__ == '__main__':
    main()

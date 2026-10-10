"""Phase 2 task 4b diagnostic: VEIT per-window centroid offsets behind the
shift + stretch fit (why real data scatter more than their errors).

Prints, for E and Ld, each window's -dlam with its error, the straight-line
fit, and the residual / error; writes nothing.

Run from the repository root:
``python -m hypernet.wiggles.phase2_wavescale_veit_windows``.
"""
import os

import numpy as np

from hypernet import srf, srfcal, whn_srf
from hypernet import whn_l1a as wl
from hypernet.wiggles import DATA_DIR


def main():
    f = wl.sequence_files(site='VEIT', seq_time='20260604T0845')
    irr, rad = wl.load_l1a_irr(f['L1A_IRR']), wl.load_l1a_rad(f['L1A_RAD'])
    models = srf.load_srf_models(os.path.join(DATA_DIR, 'veit_srf_model.json'))
    for ch, wave, scans in (('E', irr['wave'], irr['scans']),
                            ('Ld', rad['wave'], rad['Ld']['scans'])):
        mean, e_mean, _ = srf.scan_errors(scans, flatten_px=whn_srf.FLATTEN_PX)
        r, fits = srfcal.fit_wavelength_scale(wave, mean, srf=models[ch], err=e_mean,
                                              return_fits=True)
        print('\n%s: shift %+.4f ± %.4f, stretch %+.4f ± %.4f, max_err %.4f, chi2_nu %.1f'
              % (ch, r['shift'], r['shift_err'], r['stretch'], r['stretch_err'],
                 r['max_err'], r['chi2_nu']))
        g = fits[fits['ok']]
        x = (g['lam_center'] - srfcal.LAM0) / srfcal.HALF
        model = r['shift'] + r['stretch'] * x
        for (_, w), mm in zip(g.iterrows(), model):
            print('  %5.0f  -dlam %+.4f ± %.4f  resid/err %+5.1f  win chi2 %.1f'
                  % (w['lam_center'], -w['dlam'], w['dlam_err'],
                     (-w['dlam'] - mm) / w['dlam_err'], w['chi2_nu']))


if __name__ == '__main__':
    main()

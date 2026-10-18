"""Does E from L2 read wide?  The bias re-estimated with the right E grids
(lab data prompts, task 3).

Phase 0c's validation (``phase0c_validation.csv``, the VEIT sequence: L2
route vs L1A) found E from L2 ~0.06-0.09 nm wider than from L1A, and
``phase2_srf_table.l2_bias`` adds that to the L2-derived Delta FWHM
(``delta_fwhm_l1a_est``).  The lamp gives an independent check on every
instrument: if L2 inflates E, then (field - lab) is larger for E than for L.

Per Release 2 model (instrument x IRR calibration) and the VEIT L1A model,
at Hg 435.83 and 546.07 nm (the clean lines; 404.66 is a blend):

- ``d_E`` = field Gaussian FWHM_E - lab FWHM_E; ``d_L`` likewise for Ld;
- ``bias`` = d_E - d_L, the E-specific excess (the old claim: +0.06-0.09);
- for three sets of field models: ``pre_lab`` (the 122304 proxy E grid,
  ``release2_srf_models_pre_lab.json``), ``own`` (own E grid, Gaussian E,
  ``release2_srf_models.json``), ``labshape`` (own grid, super-Gaussian E
  with the lab p where the lab prefers it; ``E_sg``);
- ``lab_def`` = (lab - our Gaussian refit of the lamp line) for E minus L
  (``lab_srf.csv``): how far the lab's FWHM definition alone separates the
  channels, since the template FWHM is a Gaussian-equivalent width.

The lab epoch is ``labcal.physical_epoch`` at each channel's calibration date.

Output: ``hypernet/wiggles/lab_l2_bias.csv``.

Run from the repository root:  ``python -m hypernet.wiggles.lab_l2_bias``
"""
import json
import os

import numpy as np
import pandas as pd

from hypernet import labcal, srf
from hypernet.wiggles import DATA_DIR, WIGGLES_DIR

LINES = (435.8328, 546.0735)


def _models(path):
    doc = json.load(open(path))
    for key, chans in doc['models'].items():
        if '@' in key:
            yield key, {ch: srf.SRFModel.from_dict(d) for ch, d in chans.items()}


def _lab(sn, epoch, ch, lam, col='lab_fwhm'):
    s = labcal.load_lab_srf()
    g = s[(s['sn'] == sn) & (s['epoch'] == epoch) & (s['channel'] == ch) & (s['lam_ref'] == lam)]
    return float(g[col].iloc[0])


def rows_for(label, key, ms, cal):
    sn = key.split('@')[0][-6:]
    rad, irr = cal.split('/')
    epE, _ = labcal.physical_epoch(sn, irr)
    epL, _ = labcal.physical_epoch(sn, rad)
    out = []
    E = ms.get('E_sg', ms['E']) if label == 'labshape' else ms['E']
    for lam in LINES:
        dE = float(E.fwhm(lam)) - _lab(sn, epE, 'E', lam)
        dL = float(ms['Ld'].fwhm(lam)) - _lab(sn, epL, 'L', lam)
        defE = _lab(sn, epE, 'E', lam) - _lab(sn, epE, 'E', lam, 'fwhm_gauss')
        defL = _lab(sn, epL, 'L', lam) - _lab(sn, epL, 'L', lam, 'fwhm_gauss')
        out.append(dict(models=label, key=key, instrument='HYPSTAR_' + sn, cal_period=cal,
                        lab_epoch_E=epE, lab_epoch_L=epL, e_shape=E.shape, lam=lam, d_E=dE,
                        d_L=dL, bias=dE - dL, lab_def=defE - defL))
    return out


def main():
    rows = []
    for label, path in (('pre_lab', 'release2_srf_models_pre_lab.json'),
                        ('own', 'release2_srf_models.json'),
                        ('labshape', 'release2_srf_models.json')):
        for key, ms in _models(os.path.join(DATA_DIR, path)):
            rows += rows_for(label, key, ms, ms['E'].meta['cal_period'])
    veit = srf.load_srf_models(os.path.join(DATA_DIR, 'veit_srf_model.json'))
    veit_sg = srf.load_srf_models(os.path.join(DATA_DIR, 'veit_srf_model_lab_shape.json'))
    for label, ms in (('veit_l1a', veit), ('veit_l1a_labshape', dict(veit, E_sg=veit_sg['E']))):
        rows += rows_for('labshape' if label.endswith('labshape') else 'own',
                         'HYPSTAR_122304@2024-11-12', ms, '2024-11-08/2024-11-12')
        for r in rows[-len(LINES):]:
            r['models'] = label
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(WIGGLES_DIR, 'lab_l2_bias.csv'), index=False, float_format='%.4f')
    s = t.groupby(['models', 'key']).agg(e_shape=('e_shape', 'first'), d_E=('d_E', 'mean'),
                                         d_L=('d_L', 'mean'), bias=('bias', 'mean'),
                                         lab_def=('lab_def', 'mean')).reset_index()
    pd.set_option('display.width', 200)
    print('per model, mean over Hg 435.83 / 546.07 nm (nm):')
    print(s.to_string(index=False, float_format='%+.3f'))
    print('\nmedian over the Release 2 models:')
    r2 = s[~s['models'].str.startswith('veit')]
    print(r2.groupby('models')[['d_E', 'd_L', 'bias', 'lab_def']].median()
          .to_string(float_format='%+.3f'))
    v = pd.read_csv(os.path.join(WIGGLES_DIR, 'phase0c_validation.csv'))
    v = v[v['route'] == 'L2 route, E via E grid']
    print('\nVEIT L2-vs-L1A (phase0c_validation.csv, E via E grid), E diff - Ld diff at 450/500/600:',
          ' '.join('%+.3f' % (v[(v.channel == 'E') & (v.lam == l)]['diff'].iloc[0] -
                              v[(v.channel == 'Ld') & (v.lam == l)]['diff'].iloc[0])
                   for l in (450, 500, 600)))


if __name__ == '__main__':
    main()

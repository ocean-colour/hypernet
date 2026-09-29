"""Phase 0a, task 8: the numbers behind the 0a interim report (Q&A Build #8).

Reads the committed task 5-7 products and prints:

1. the FWHM(lambda) models of ``hypernet/data/veit_srf_model.json`` at
   400-850 nm with their propagated errors (cov inflated by chi2_nu), and
   the model differences Ld - E and Lu - E with their significance (the
   three models are fitted independently, so the errors add in quadrature);
2. per clean line, the E - L FWHM difference and its significance, and the
   same for the blend lines;
3. relative centroid offsets L - E per line and their weighted means (clean
   lines), against the 0.1 nm threshold of G0(c); the absolute offsets in air
   and vacuum are in ``phase0a_veit.py``'s stdout;
4. the Gaussian-adequacy summary: per-scan chi2_nu by line class, and the
   FWHM(lambda) chi2_nu.

Run from the repository root: ``python wiggles/phase0a_interim.py``.
"""
import os
import sys

import numpy as np
import pandas as pd

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, REPO)
from hypernet import srf  # noqa: E402

mods = srf.load_srf_models(os.path.join(REPO, 'hypernet', 'data', 'veit_srf_model.json'))
lines = pd.read_csv(os.path.join(REPO, 'wiggles', 'phase0_veit_lines.csv'))
scans = pd.read_csv(os.path.join(REPO, 'wiggles', 'phase0_veit_scans.csv'))
lines = lines[lines['err_kind'] == 'flat']

print('[1] FWHM(lambda) models, nm (cov x chi2_nu)')
lam = np.array([400., 450., 500., 550., 600., 650., 700., 750., 800., 850.])
tab = pd.DataFrame({'lam': lam})
for ch, m in mods.items():
    tab[ch] = m.fwhm(lam)
    tab[ch + '_err'] = m.fwhm_err(lam)
for ch in ('Ld', 'Lu'):
    tab[ch + '-E'] = tab[ch] - tab['E']
    tab[ch + '-E_err'] = np.hypot(tab[ch + '_err'], tab['E_err'])
    tab[ch + '-E_z'] = tab[ch + '-E'] / tab[ch + '-E_err']
print(tab.to_string(index=False, float_format='%.2f'))
for ch, m in mods.items():
    e = np.sqrt(np.diag(m.cov))
    print('  %-2s coeffs %s +- %s  chi2_nu %.1f  npts %d  range %.0f-%.0f nm' % (
        ch, np.round(m.coeffs, 3), np.round(e, 3), m.chi2_nu, m.npts, m.lam_min, m.lam_max))

print('\n[2] per line: FWHM_L - FWHM_E and significance')
piv = lines.pivot(index='name', columns='channel', values='fwhm')
err = lines.pivot(index='name', columns='channel', values='fwhm_err')
meta = lines[lines['channel'] == 'E'].set_index('name')[['lam_air', 'blend', 'use_for_srf']]
out = meta.copy()
for ch in ('Ld', 'Lu'):
    out['d_' + ch] = piv[ch] - piv['E']
    out['z_' + ch] = out['d_' + ch] / np.hypot(err[ch], err['E'])
out = out.loc[list(srf.LINES['name'])]
print(out.to_string(float_format='%.2f'))
use = out['use_for_srf']
for ch in ('Ld', 'Lu'):
    for lab, sel in (('clean', use & ~out['blend']), ('blend', use & out['blend'])):
        z = out.loc[sel, 'z_' + ch]
        print('  %-2s - E, %s SRF lines: %d of %d at z > 3; range %.2f to %.2f nm'
              % (ch, lab, (z > 3).sum(), sel.sum(), out.loc[sel, 'd_' + ch].min(),
                 out.loc[sel, 'd_' + ch].max()))

print('\n[3] relative centroid offsets L - E (nm)')
mu = lines.pivot(index='name', columns='channel', values='mu')
mue = lines.pivot(index='name', columns='channel', values='mu_err')
cen = meta.copy()
for ch in ('Ld', 'Lu'):
    cen['d_' + ch] = mu[ch] - mu['E']
    cen['e_' + ch] = np.hypot(mue[ch], mue['E'])
cen = cen.loc[list(srf.LINES['name'])]
print(cen.to_string(float_format='%.3f'))
for ch in ('Ld', 'Lu'):
    for lab, sel in (('clean', cen['use_for_srf'] & ~cen['blend']),
                     ('all SRF', cen['use_for_srf'])):
        d, e = cen.loc[sel, 'd_' + ch], cen.loc[sel, 'e_' + ch]
        w = 1 / e ** 2
        mean = np.sum(w * d) / w.sum()
        chi2 = np.sum(w * (d - mean) ** 2) / max(len(d) - 1, 1)
        emean = np.sqrt(max(chi2, 1.0) / w.sum())
        print('  %-2s - E, %-7s lines: weighted mean %+.3f +- %.3f nm (chi2_nu %.1f); '
              'max |d| %.3f; n |d| > 0.1 nm: %d of %d'
              % (ch, lab, mean, emean, chi2, np.abs(d).max(), (np.abs(d) > 0.1).sum(), len(d)))

print('\n[4] Gaussian adequacy: median per-scan chi2_nu by line class')
s = scans.merge(meta.reset_index()[['name']], on='name')
s['cls'] = np.where(~s['use_for_srf'], 'band (diag.)',
                    np.where(s['blend'], 'blend', 'clean'))
print(s.groupby(['cls', 'channel'])['chi2_nu_scan_med'].median().unstack().round(2))

"""Phase 2, task 9b (Build #4b Q1): is the common window pattern in HSRS?

The per-window centroid offsets of every HYPSTAR (E and Ld) share a ~0.13 nm
rms pattern (task 4b).  Part B of ``phase2_shape_check.py`` shows that no
plausible SRF shape produces it.  Here a second solar atlas stands in for the
real Sun:

- Kurucz (2005) high-resolution irradiance (``irradthu.dat``, 300-1000 nm,
  Thuillier calibration; fetched from kurucz.harvard.edu, SHA-256
  ``6013b1398718ed9488cdcbfd44ef93e8adcd7fd536a329d35252f5fdf7f330d0``,
  cached under ``$OS_COLOR/hypernet/ref``).
- Its wavelength frame (air or vacuum) is found empirically: the median
  centroid offset against HSRS in air should be ~0 for air and ~+0.16 nm at
  600 nm for vacuum.
- Kurucz observed through a Gaussian (the 122304 E FWHM) on the 122304 E grid
  is fitted against HSRS (air) with the same windows and method as the real
  data (``srfcal.fit_wavelength_scale``, SRF width free).  If HSRS's line
  positions or shapes are off in some windows, the Kurucz "instrument" shows
  the real pattern.

A second candidate is tested the same way: **weak telluric lines inside the
"clean" windows**.  The real E carries O2/H2O/O3 absorption that bare HSRS
does not, so the twin's Ed (HSRS x the HAPI/O3 transmittances, VEIT-like
scene) observed through a Gaussian is fitted against bare HSRS.

Third, on the real VEIT E and Ld: the wavelength-scale fit against bare
HSRS (as task 4b) and against the Emod for the sequence (HSRS x T, which
carries the telluric lines; Ld with the diffuse air mass), raw (no pattern
removed): does chi2_nu and the per-window scatter shrink?

Output: ``hypernet/wiggles/phase2_reference_check.csv`` (per window: Kurucz-vs-
HSRS -dlam about its line, the real pattern) and a printed summary.

Run from the repository root: ``python -m hypernet.wiggles.phase2_reference_check``.
"""
import os

import numpy as np
import pandas as pd

from hypernet import srf, srfcal, twin
from hypernet.wiggles import WIGGLES_DIR

KURUCZ = os.path.join(os.getenv('OS_COLOR', '.'), 'hypernet', 'ref', 'kurucz2005_irradthu.dat')


def load_kurucz(lo=375.0, hi=905.0):
    rows = []
    with open(KURUCZ) as fh:
        for line in fh:
            p = line.split()
            if len(p) == 2:
                try:
                    rows.append((float(p[0]), float(p[1])))
                except ValueError:
                    pass
    a = np.array(rows)
    a = a[(a[:, 0] >= lo) & (a[:, 0] <= hi)]
    return a[:, 0], a[:, 1]


def instrument(lam, flux, grid, fwhm_fn):
    return srf.convolve_to_grid(lam, flux, grid, fwhm_fn, cache=False)


def main():
    hsrs = srfcal.default_reference('air')
    lam_k, f_k = load_kurucz()
    grids = twin.load_grids()
    gE = grids['grid_E_122304']
    gE = gE[(gE > 385) & (gE < 895)]
    fE = twin.load_instrument_srfs()['HYPSTAR_122304']['E'].fwhm
    real = srfcal.load_window_pattern()
    out = {}
    for frame in ('as_given', 'vac_to_air'):
        lam = lam_k if frame == 'as_given' else srf.vac_to_air(lam_k)
        spec = instrument(lam, f_k, gE, fE)
        res, fits = srfcal.fit_wavelength_scale(gE, spec, ref=hsrs, return_fits=True)
        out[frame] = (res, fits)
        print('Kurucz %-10s vs HSRS(air): shift %+.3f nm, stretch %+.3f nm (dlam at 630 nm '
              '%+.3f), %d windows' % (frame, res['shift'], res['stretch'], -res['shift'],
                                     res['n_windows']))
    # air-vacuum at 600 nm is ~0.166 nm: the frame with the smaller |shift| is right
    frame = min(out, key=lambda k: abs(out[k][0]['shift']))
    res, fits = out[frame]
    g = fits[fits['ok']].copy()
    x = (g['lam_center'] - srfcal.LAM0) / srfcal.HALF
    g['kurucz_resid'] = -g['dlam'] - (res['shift'] + res['stretch'] * x)
    g['real_pattern'] = [real.get(float(c), np.nan) for c in g['lam_center']]
    ok = np.isfinite(g['real_pattern'])
    r = float(np.corrcoef(g.loc[ok, 'kurucz_resid'], g.loc[ok, 'real_pattern'])[0, 1])
    slope = float(np.polyfit(g.loc[ok, 'kurucz_resid'], g.loc[ok, 'real_pattern'], 1)[0])
    rms_k = float(np.sqrt(np.mean(g.loc[ok, 'kurucz_resid'] ** 2)))
    rms_r = float(np.sqrt(np.mean(g.loc[ok, 'real_pattern'] ** 2)))
    resid = g.loc[ok, 'real_pattern'] - g.loc[ok, 'kurucz_resid']
    print('frame used: %s.  Per-window pattern, Kurucz-through-a-Gaussian vs HSRS: rms '
          '%.3f nm; real HYPSTAR pattern rms %.3f nm; correlation %.2f, slope (real on '
          'Kurucz) %.2f; real minus Kurucz rms %.3f nm'
          % (frame, rms_k, rms_r, r, slope, float(np.sqrt(np.mean(resid ** 2)))))
    # candidate 2: weak telluric lines (twin Ed = HSRS x T, fitted against bare HSRS)
    from hypernet import emod
    fields = twin.load_fields()
    for sza, pwv in ((40.0, 15.0), (40.0, 30.0), (70.0, 30.0)):
        e = emod.build_emod(sza, pwv_mm=pwv)
        i = 0 if sza == 40.0 else list(np.asarray(fields['name']).astype(str)).index(
            'sza70_aot0.05_chl1')
        sc = twin.compose_scene(e, fields, i)
        spec = instrument(sc['lam'], sc['Ed'], gE, fE)
        res2, fits2 = srfcal.fit_wavelength_scale(gE, spec, ref=hsrs, return_fits=True)
        g2 = fits2[fits2['ok']].set_index('lam_center')
        x2 = (g2.index.values - srfcal.LAM0) / srfcal.HALF
        tr = -g2['dlam'] - (res2['shift'] + res2['stretch'] * x2)
        col = 'telluric_sza%d_pwv%d' % (sza, pwv)
        g[col] = [tr.get(c, np.nan) for c in g['lam_center']]
        m2 = ok & np.isfinite(g[col])
        print('telluric (SZA %d, PWV %d mm) vs bare HSRS: pattern rms %.3f nm, correlation with '
              'real %.2f, shift %+.3f stretch %+.3f'
              % (sza, pwv, float(np.sqrt(np.mean(g.loc[m2, col] ** 2))),
                 float(np.corrcoef(g.loc[m2, col], g.loc[m2, 'real_pattern'])[0, 1]),
                 res2['shift'], res2['stretch']))
    # candidate 2 on real data: fit VEIT against Emod (with telluric lines)
    from hypernet import whn_srf
    from hypernet import whn_l1a as wl
    from hypernet.wiggles import DATA_DIR
    f = wl.sequence_files(site='VEIT', seq_time='20260604T0845')
    irr, rad = wl.load_l1a_irr(f['L1A_IRR']), wl.load_l1a_rad(f['L1A_RAD'])
    models = srf.load_srf_models(os.path.join(DATA_DIR, 'veit_srf_model.json'))
    for pwv in (7.5, 15.0, 30.0):
        e = emod.build_emod(36.7, pwv_mm=pwv)
        refs = {'E': (e['lam'], e['F0'] * (0.85 * e['T_direct'] + 0.15 * e['T_diffuse'])),
                'Ld': (e['lam'], e['F0'] * e['T_diffuse'])}
        for ch, wave, scans in (('E', irr['wave'], irr['scans']),
                                ('Ld', rad['wave'], rad['Ld']['scans'])):
            mean, err, _ = srf.scan_errors(scans, flatten_px=whn_srf.FLATTEN_PX)
            m = (wave > 385) & (wave < 990)
            for name, ref in (('bare HSRS', hsrs), ('Emod PWV %g' % pwv, refs[ch])):
                if name == 'bare HSRS' and pwv != 15.0:
                    continue
                rr, ff = srfcal.fit_wavelength_scale(wave[m], mean[m], srf=models[ch],
                                                     ref=ref, err=err[m], return_fits=True)
                gg = ff[ff['ok']]
                xx = (gg['lam_center'] - srfcal.LAM0) / srfcal.HALF
                res_ = -gg['dlam'] - (rr['shift'] + rr['stretch'] * xx)
                print('VEIT %-2s vs %-12s: chi2_nu %5.1f, per-window rms about the line '
                      '%.3f nm, max_err %.3f, shift %+.3f stretch %+.3f'
                      % (ch, name, rr['chi2_nu'], float(np.sqrt(np.mean(res_ ** 2))),
                         rr['max_err'], rr['shift'], rr['stretch']))
    t = g[['lam_center', 'dlam', 'dlam_err', 'kurucz_resid', 'real_pattern'] +
          [c for c in g.columns if c.startswith('telluric_')]].assign(
        frame=frame, corr=r, slope=slope)
    t.to_csv(os.path.join(WIGGLES_DIR, 'phase2_reference_check.csv'), index=False,
             float_format='%.4f')
    print(t[['lam_center', 'kurucz_resid', 'real_pattern'] +
            [c for c in t.columns if c.startswith('telluric_')]].to_string(
        index=False, float_format='%+.3f'))
    return t


if __name__ == '__main__':
    main()

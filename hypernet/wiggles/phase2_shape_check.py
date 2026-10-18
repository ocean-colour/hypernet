"""Phase 2, task 9b: non-Gaussian SRFs in the twin (G1 concern 3), the origin
of the common window pattern (Build #4b Q1), and an L wavelength error
(Setup Q8).

**Part A -- does the srf gain survive a non-Gaussian SRF?**  Twin case (iii)
on the 25 Phase 1 scenes, 122304 grids.  The truth and every observation (E
on the E grid; Lu, Ld, Lw and Ed on the L grid) go through a non-Gaussian
SRF with the 122304 FWHM(lambda) (``release2_srf_models.json``); the
correction is still Gaussian, with either

- ``true_fwhm``: the true FWHM(lambda), or
- ``fitted``: the Gaussian FWHM a template fit would measure
  (``srfcal.calibrate_srf`` on the noise-free E and Ld, against the scene's
  own Ed and Ld): the realistic pipeline.

Shape cases (``CASES``): the Gaussian control; the Phase 0 fits
(``phase0d_shape_summary.csv``: pseudo-Voigt eta 0.06, super-Gaussian p =
2.23 / 2.0 (E / Ld, weighted means) and 2.54 / 2.27 (medians)); an
asymmetric Gaussian (sigma (1 +- a) either side, FWHM kept) with a = 0.05,
0.10; a Gaussian with the L pixels really sampling +0.05 / +0.10 nm (Q8);
and, for the prediction, E = L - dFWHM (0.2, 0.8 nm) for the Gaussian, the
median super-Gaussian and asym 0.10.

**The lab shape** (lab data prompts task 3): ``lab_E_supergauss`` has E
super-Gaussian with the lab p of 122304 (``labcal.srf_shape_for`` at its
2024-11 calibration) and L Gaussian, as the lamp lines show, at the Gaussian
FWHMs (``_dfwhm0.2``: E = L - 0.2 nm).  Besides ``true_fwhm`` and ``fitted``
(Gaussian corrections) it has the variant ``shaped``: the correction knows
the E shape (``srf.ShapedFWHM`` with the true FWHM and p), as the SRF table
now does.

Metric as Phase 1 task 9: rms(rho_w'' - rho_w''_true), h = 1 nm, within +-5 nm
of the ten lines; ``reduction`` = 1 - rms / rms_linear (same case).  Also the
correlation of the srf residual with Ed''/Ed near the lines (G3's trigger).

**Part B -- where does the common window pattern come from?**  On the
VEIT-like scene, ``srfcal.fit_wavelength_scale`` (SRF width free, as the
Phase 0c fits) on the noise-free E and Ld observed through each shape, with
the scene's own spectra as the reference.  A Gaussian gives dlam = 0 in every
window; a shape that distorts line profiles gives a per-window pattern.  It
is compared with the real one (``hypernet/data/wavescale_window_pattern.csv``):
correlation, regression slope, rms.

Outputs: ``$OS_COLOR/hypernet/wiggles/phase2/shape_check.parquet``;
``hypernet/wiggles/phase2_shape_check.csv`` (part A medians over scenes) and
``phase2_shape_pattern.csv`` (part B); ``figs/phase2/shape_check.png``.

Run from the repository root:
``python -m hypernet.wiggles.phase2_shape_check [--workers 12] [--scenes N]``.
"""
import argparse
import os
import time

import numpy as np
import pandas as pd

from hypernet import edinterp, emod, labcal, srf, srfcal, twin
from hypernet.wiggles import FIGDIR, WIGGLES_DIR
from hypernet.wiggles.phase1_twin import L_RANGE, masks, second_diff

ROOT = os.path.join(os.getenv('OS_COLOR', '.'), 'hypernet', 'wiggles', 'phase2')
FIG2 = FIGDIR.replace('phase0', 'phase2')
INST = 'HYPSTAR_122304'
METHODS = ('linear', 'ruddick2023', 'srf')
H = 1.0


def _case(name, shape='gauss', par_E=None, par_L=None, dfwhm=None, l_shift=0.0,
          shape_L=None, shaped=False):
    return dict(name=name, shape=shape, par_E=par_E, par_L=par_L, dfwhm=dfwhm,
                l_shift=l_shift, shape_L=shape if shape_L is None else shape_L,
                shaped=shaped)


#: The lab E exponent of the twin instrument (122304 at its 2024-11 calibration).
P_LAB = labcal.srf_shape_for(INST, 'E', '2024-11-12')['shape_p']


CASES = [
    _case('gauss'),
    _case('pvoigt_0.06', 'pvoigt', 0.06, 0.06),
    _case('supergauss_wmean', 'supergauss', 2.23, 2.00),
    _case('supergauss_median', 'supergauss', 2.54, 2.27),
    _case('asym_0.05', 'asym', 0.05, 0.05),
    _case('asym_0.10', 'asym', 0.10, 0.10),
    _case('Lshift_+0.05', l_shift=0.05),
    _case('Lshift_+0.10', l_shift=0.10),
]
CASES += [_case('lab_E_supergauss', 'supergauss', P_LAB, None, shape_L='gauss', shaped=True),
          _case('lab_E_supergauss_dfwhm0.2', 'supergauss', P_LAB, None, dfwhm=0.2,
                shape_L='gauss', shaped=True)]
for _d in (0.2, 0.8):
    CASES += [_case('gauss_dfwhm%.1f' % _d, dfwhm=_d),
              _case('supergauss_median_dfwhm%.1f' % _d, 'supergauss', 2.54, 2.27, dfwhm=_d),
              _case('asym_0.10_dfwhm%.1f' % _d, 'asym', 0.10, 0.10, dfwhm=_d)]


def kernel(dx, sigma, shape, par):
    """SRF weights at offsets dx = pixel - lambda; 'asym' is a split Gaussian
    with sigma (1 + a) for dx > 0 and sigma (1 - a) for dx < 0 (same FWHM)."""
    if shape == 'asym':
        s = np.where(dx > 0, sigma * (1 + par), sigma * (1 - par))
        return np.exp(-0.5 * (dx / s) ** 2)
    return srf.srf_kernel(dx, sigma, shape, par)


def observe_shape(lam, specs, grid, fwhm_fn, shape, par, chunk=96):
    """Observe several 0.01 nm spectra (rows of ``specs``) through a
    (non-Gaussian) SRF of FWHM ``fwhm_fn(grid)`` on ``grid``; normalised
    weights, truncated at +-6 sigma (+-8 nm for the pseudo-Voigt's wings)."""
    specs = np.atleast_2d(specs)
    sig = np.asarray(fwhm_fn(grid), float) / srf.FWHM_PER_SIGMA
    dw = np.gradient(lam)
    out = np.empty((specs.shape[0], grid.size))
    for a in range(0, grid.size, chunk):
        g, s = grid[a:a + chunk], sig[a:a + chunk]
        half = 8.0 if shape == 'pvoigt' else 6.0 * s.max() * (1.2 if shape == 'asym' else 1)
        j0, j1 = np.searchsorted(lam, g.min() - half), np.searchsorted(lam, g.max() + half)
        dx = g[:, None] - lam[None, j0:j1]
        k = kernel(dx, s[:, None], shape, par) * (np.abs(dx) <= half) * dw[None, j0:j1]
        k /= k.sum(axis=1, keepdims=True)
        out[:, a:a + chunk] = specs[:, j0:j1] @ k.T
    return out


def run_scene(args):
    i_case, = args
    fields = twin.load_fields()
    lib = twin.rhow_library(fields)
    grids = twin.load_grids()
    srfs = twin.load_instrument_srfs()[INST]
    sza = float(fields['sza'][i_case])
    scene = twin.compose_scene(emod.build_emod(sza), fields, i_case, lib=lib)
    lam = scene['lam']
    name = str(np.asarray(fields['name'])[i_case])
    gE = grids['grid_E_122304']
    gE = gE[(gE >= L_RANGE[0] - 3) & (gE <= L_RANGE[1] + 3)]
    gL = grids['grid_L_122304']
    gL = gL[(gL >= L_RANGE[0]) & (gL <= L_RANGE[1])]
    near, away, _ = masks(gL)
    rho_eff = np.interp(gL, lam, scene['rho_eff'])
    mdl = (lam, scene['Ed'])
    rows, pattern_rows = [], []
    for c in CASES:
        fL = srfs['Ld'].fwhm
        fE = srfs['E'].fwhm if c['dfwhm'] is None else (lambda l, d=c['dfwhm']: fL(l) - d)
        tL = gL + c['l_shift']
        Eobs = observe_shape(lam, scene['Ed'], gE, fE, c['shape'], c['par_E'])[0]
        Lw_, Ed_, Lu, Ld = observe_shape(lam, np.vstack([scene['Lw'], scene['Ed'], scene['Lu'],
                                                         scene['Ld']]), tL, fL, c['shape_L'],
                                         c['par_L'])
        truth = np.pi * Lw_ / Ed_
        t2 = second_diff(gL, truth, H)
        edpp = second_diff(gL, Ed_, H) / Ed_
        corr = {'true_fwhm': (fE, fL)}
        if c['shape'] != 'gauss' or c['name'] == 'gauss':
            mE = srfcal.calibrate_srf(gE, Eobs, ref=(lam, scene['Ed']))
            mL = srfcal.calibrate_srf(gL, Ld, channel='Ld', ref=(lam, scene['Ld']))
            corr['fitted'] = (mE, mL)
        if c['shaped']:
            corr['shaped'] = (srf.ShapedFWHM(fE, c['shape'], c['par_E']), fL)
        lin = None
        for variant, (a, b) in corr.items():
            for m in METHODS:
                EdL = edinterp.interpolate_ed_to_l(gE, Eobs, gL, emod=mdl, srf_irr=a,
                                                   srf_rad=b, method=m)
                rho = np.pi * (Lu - rho_eff * Ld) / EdL
                d = second_diff(gL, rho, H) - t2
                rms = float(np.sqrt(np.mean(d[near] ** 2)))
                if m == 'linear':
                    lin = rms
                r = dict(scene=name, case=c['name'], shape=c['shape'], variant=variant,
                         method=m, rms_line=rms, rms_away=float(np.sqrt(np.mean(d[away] ** 2))),
                         rms_true_line=float(np.sqrt(np.mean(t2[near] ** 2))),
                         reduction=1 - rms / lin,
                         corr_resid_edpp=float(np.corrcoef(d[near], edpp[near])[0, 1]))
                if variant == 'fitted':
                    r.update(fitted_dfwhm_E=float(np.median(mE.fwhm(np.arange(400, 681)) -
                                                            fE(np.arange(400.0, 681)))),
                             fitted_dfwhm_L=float(np.median(mL.fwhm(np.arange(400, 681)) -
                                                            fL(np.arange(400.0, 681)))))
                rows.append(r)
        # part B, on the VEIT-like scene only
        if i_case == 0 and c['l_shift'] == 0 and c['dfwhm'] is None:
            for ch, wav, spec, ref in (('E', gE, Eobs, scene['Ed']), ('Ld', gL, Ld, scene['Ld'])):
                res, fits = srfcal.fit_wavelength_scale(wav, spec, ref=(lam, ref),
                                                        return_fits=True)
                g = fits[fits['ok']]
                x = (g['lam_center'] - srfcal.LAM0) / srfcal.HALF
                resid = -g['dlam'] - (res['shift'] + res['stretch'] * x)
                for lc, v, dv in zip(g['lam_center'], -g['dlam'], resid):
                    pattern_rows.append(dict(case=c['name'], channel=ch, lam_center=float(lc),
                                             minus_dlam=float(v), resid=float(dv),
                                             shift=res['shift'], stretch=res['stretch']))
    return rows, pattern_rows


def summarise(df):
    keys = ['case', 'shape', 'variant', 'method']
    s = df.groupby(keys).agg(rms_line=('rms_line', 'median'), rms_away=('rms_away', 'median'),
                             rms_true_line=('rms_true_line', 'median'),
                             reduction=('reduction', 'median'),
                             corr_resid_edpp=('corr_resid_edpp', 'median'),
                             n_scenes=('scene', 'nunique')).reset_index()
    if 'fitted_dfwhm_E' in df:
        f = df.groupby(keys)[['fitted_dfwhm_E', 'fitted_dfwhm_L']].median()
        s = s.join(f, on=keys)
    order = {c['name']: i for i, c in enumerate(CASES)}
    return s.sort_values(['case', 'variant', 'method'], key=lambda col: col.map(order)
                         if col.name == 'case' else col).reset_index(drop=True)


def pattern_compare(p):
    real = srfcal.load_window_pattern()
    out = []
    for (case, ch), g in p.groupby(['case', 'channel']):
        r = np.array([real.get(float(c), np.nan) for c in g['lam_center']])
        v = g['resid'].values
        ok = np.isfinite(r)
        has = np.sqrt(np.mean(v[ok] ** 2)) > 1e-3        # the Gaussian control has none
        rr = float(np.corrcoef(v[ok], r[ok])[0, 1]) if has else np.nan
        slope = float(np.sum(v[ok] * r[ok]) / np.sum(v[ok] ** 2)) if has else np.nan
        out.append(dict(case=case, channel=ch, n_windows=int(ok.sum()),
                        rms_pattern_twin=float(np.sqrt(np.mean(v[ok] ** 2))),
                        rms_pattern_real=float(np.sqrt(np.mean(r[ok] ** 2))),
                        corr_with_real=rr, real_over_twin_slope=slope,
                        shift=float(g['shift'].iloc[0]), stretch=float(g['stretch'].iloc[0])))
    return pd.DataFrame(out)


def figure(s, p, pc):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    os.makedirs(FIG2, exist_ok=True)
    fig, axs = plt.subplots(1, 2, figsize=(15, 5.6), gridspec_kw=dict(width_ratios=[1.2, 1]))
    ax = axs[0]
    main = [c['name'] for c in CASES]
    q = s[(s['method'] == 'srf')]
    x = np.arange(len(main))
    for k, (variant, col) in enumerate((('true_fwhm', '#1baf7a'), ('fitted', '#2a78d6'),
                                        ('shaped', '#4a3aa7'))):
        g = q[q['variant'] == variant].set_index('case')
        y = [100 * g.loc[c, 'reduction'] if c in g.index else np.nan for c in main]
        ax.plot(x + (k - 1) * 0.22, y, 'o', ms=6, color=col,
                label='srf, correction with the true E shape' if variant == 'shaped' else
                'srf, Gaussian correction with %s' % variant.replace('_', ' '))
    g = s[(s['method'] == 'ruddick2023') & (s['variant'] == 'true_fwhm')].set_index('case')
    ax.plot(x, [100 * g.loc[c, 'reduction'] if c in g.index else np.nan for c in main], 's',
            ms=5, mfc='none', color='#888', label='ruddick2023')
    ax.axhline(80, color='#888', lw=1, ls='-.')
    ax.axhline(0, color='#333', lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(main, rotation=50, ha='right', fontsize=8)
    ax.set_ylabel("reduction of the line-region ρw'' error vs linear (%)")
    ax.set_ylim(-60, 105)
    ax.set_title('(a) truth through a non-Gaussian SRF, 122304, 25 scenes (median)',
                 fontsize=10, loc='left')
    ax.legend(fontsize=8, frameon=False, loc='lower left')
    ax = axs[1]
    real = srfcal.load_window_pattern()
    lc = np.array(sorted(real))
    ax.plot(lc, [real[c] for c in lc], 'k-o', ms=4, lw=1.5, label='real (Release 2, common)')
    for case, col in (('asym_0.10', '#e34948'), ('asym_0.05', '#eda100'),
                      ('supergauss_median', '#2a78d6'), ('pvoigt_0.06', '#1baf7a')):
        g = p[(p['case'] == case) & (p['channel'] == 'E')].sort_values('lam_center')
        r = pc[(pc['case'] == case) & (pc['channel'] == 'E')]
        lab = '%s, E (r = %.2f)' % (case, r['corr_with_real'].iloc[0]) if len(r) else case
        ax.plot(g['lam_center'], g['resid'], '-o', ms=3, lw=1, color=col, label=lab)
    ax.axhline(0, color='#333', lw=0.8)
    ax.set_xlabel('window centre (nm)')
    ax.set_ylabel('−dlam about the shift + stretch line (nm)')
    ax.set_title('(b) per-window pattern: twin shapes vs the real common pattern', fontsize=10,
                 loc='left')
    ax.legend(fontsize=8, frameon=False)
    for a in axs:
        a.grid(alpha=0.25, lw=0.6)
        a.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    path = os.path.join(FIG2, 'shape_check.png')
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(description='Phase 2 task 9b: non-Gaussian SRF check')
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--scenes', type=int, default=None)
    a = ap.parse_args(argv)
    fields = twin.load_fields()
    n = len(fields['name']) if a.scenes is None else a.scenes
    for sza in sorted(set(float(x) for x in fields['sza'][:n])):
        emod.build_emod(sza)
    t0 = time.time()
    from multiprocessing import Pool
    rows, prow = [], []
    with Pool(a.workers) as pool:
        for r, pr in pool.imap_unordered(run_scene, [(i,) for i in range(n)]):
            rows += r
            prow += pr
    df, p = pd.DataFrame(rows), pd.DataFrame(prow)
    os.makedirs(ROOT, exist_ok=True)
    df.to_parquet(os.path.join(ROOT, 'shape_check.parquet'))
    s = summarise(df)
    s.to_csv(os.path.join(WIGGLES_DIR, 'phase2_shape_check.csv'), index=False,
             float_format='%.4g')
    pc = pattern_compare(p)
    p.to_csv(os.path.join(WIGGLES_DIR, 'phase2_shape_pattern.csv'), index=False,
             float_format='%.4f')
    path = figure(s, p, pc)
    pd.set_option('display.width', 220)
    print('%d rows in %.0f s' % (len(df), time.time() - t0))
    print(s.pivot_table(index='case', columns=['variant', 'method'], values='reduction')
          .reindex([c['name'] for c in CASES]).to_string(float_format='%.2f'))
    print('\nsrf: corr(residual, Ed\'\'/Ed) near lines and fitted FWHM bias:')
    q = s[s['method'] == 'srf']
    print(q[['case', 'variant', 'rms_line', 'rms_true_line', 'corr_resid_edpp'] +
            [c for c in ('fitted_dfwhm_E', 'fitted_dfwhm_L') if c in q]]
          .to_string(index=False, float_format='%.3g'))
    print('\npart B, per-window pattern vs the real one:')
    print(pc.to_string(index=False, float_format='%.3f'))
    print('wrote', path)


if __name__ == '__main__':
    main()

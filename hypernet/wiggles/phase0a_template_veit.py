"""Phase 0a, task 3b: HSRS template SRF fit of the VEIT sample.

Fits E, Ld and Lu against TSIS-1 HSRS (in air) convolved with a Gaussian SRF,
in contiguous 10 nm windows over 390-880 nm with the telluric bands excluded
(``srf.template_windows``), and fits FWHM(lambda) per channel.  Each channel
is run twice, with the additive veil free and with it fixed at 0.

Outputs:

- ``$OS_COLOR/hypernet/wiggles/phase0/veit_template.parquet`` and the
  committed copy ``hypernet/wiggles/phase0_veit_template.csv`` (every window, both
  veil settings);
- ``hypernet/data/veit_srf_model.json`` -- the template SRF models handed to
  Phase 1 (Q&A Q10), with the veil setting chosen below;
- ``hypernet/wiggles/figs/phase0/veit_template_fwhm.png``;
- on stdout: template vs empirical FWHM, E vs L, the absolute offsets
  against air and vacuum, the veil, and the chi2.

Run from the repository root: ``python -m hypernet.wiggles.phase0a_template_veit``.
"""
import os

import numpy as np
import pandas as pd

from hypernet import srf, whn_srf  # noqa: E402
from hypernet import whn_l1a as wl  # noqa: E402
from hypernet.wiggles import DATA_DIR, FIGDIR, OUT, REPO, WIGGLES_DIR  # noqa: F401

COLORS = {'E': '#2a78d6', 'Ld': '#eb6834', 'Lu': '#1baf7a'}
MARKERS = {'E': 'o', 'Ld': 's', 'Lu': '^'}
#: Veil setting of the shipped models.  The free veil is consistent with 0 in
#: every channel (median -0.03 +- 0.05) and makes a third of the windows fail,
#: so the shipped models fix it at 0; the veil-free fits are kept as a check.
SHIP_VEIL = {'E': False, 'Ld': False, 'Lu': False}


def run():
    files = wl.sequence_files(site='VEIT', seq_time='20260604T0845')
    irr = wl.load_l1a_irr(files['L1A_IRR'])
    rad = wl.load_l1a_rad(files['L1A_RAD'])
    ref = whn_srf.load_reference('air')
    res = {}
    for veil in (True, False):
        res[veil] = whn_srf.fit_sequence(irr, rad, ref=ref, template=True,
                                         empirical=not veil, veil=veil)
    return irr, rad, res


def model_of(res, method, ch, veil=None):
    for m in res.models if hasattr(res, 'models') else res['models']:
        if m.channel == ch and m.meta['method'] == method and \
                (veil is None or m.meta.get('veil') == veil):
            return m
    raise KeyError((method, ch, veil))


def report(res):
    tv, tf = res[True]['template'], res[False]['template']
    tv = tv.assign(veil_fit=True)
    tf = tf.assign(veil_fit=False)
    tab = pd.concat([tv, tf], ignore_index=True)
    lam = np.array([400, 450, 500, 550, 600, 650, 700, 800, 850])
    emp = {m.channel: m for m in res[False]['models'] if m.meta['method'] == 'empirical'}
    print('\nFWHM(lambda), nm: template (veil free | veil 0) vs empirical')
    rows = []
    for ch in whn_srf.CHANNELS:
        mv = model_of(res[True], 'template', ch)
        mf = model_of(res[False], 'template', ch)
        for l in lam:
            rows.append(dict(channel=ch, lam=l, tmpl_veil=mv.fwhm(l), err_veil=mv.fwhm_err(l),
                             tmpl_noveil=mf.fwhm(l), err_noveil=mf.fwhm_err(l),
                             empirical=emp[ch].fwhm(l)))
        print('  %-2s veil free: c = %s, chi2_nu %.2f, npts %d; veil 0: c = %s, chi2_nu %.2f'
              % (ch, np.round(mv.coeffs, 3), mv.chi2_nu, mv.npts, np.round(mf.coeffs, 3),
                 mf.chi2_nu))
    t = pd.DataFrame(rows)
    print(t.pivot(index='lam', columns='channel',
                  values=['tmpl_veil', 'tmpl_noveil', 'empirical']).round(2).to_string())
    # E vs L at the model level, shipped veil settings
    print('\nshipped models (veil: %s): L - E, nm' % SHIP_VEIL)
    mE = model_of(res[SHIP_VEIL['E']], 'template', 'E')
    for ch in ('Ld', 'Lu'):
        mL = model_of(res[SHIP_VEIL[ch]], 'template', ch)
        d = mL.fwhm(lam) - mE.fwhm(lam)
        e = np.hypot(mL.fwhm_err(lam), mE.fwhm_err(lam))
        print('  %-2s - E: ' % ch + '  '.join('%d: %+.2f±%.2f' % (l, a, b)
                                              for l, a, b in zip(lam, d, e)))
    # per window: E vs L in the same window (shipped settings)
    ship = pd.concat([tab[(tab['channel'] == ch) & (tab['veil_fit'] == SHIP_VEIL[ch])]
                      for ch in whn_srf.CHANNELS])
    ship = ship[ship['ok']]
    piv = ship.pivot(index='name', columns='channel', values='fwhm')
    pe = ship.pivot(index='name', columns='channel', values='fwhm_err')
    for ch in ('Ld', 'Lu'):
        z = (piv[ch] - piv['E']) / np.hypot(pe[ch], pe['E'])
        print('  per window %s - E: median %+.2f nm, %d of %d windows at z > 3'
              % (ch, np.nanmedian(piv[ch] - piv['E']), (z > 3).sum(), z.notna().sum()))
    # absolute offsets
    print('\nabsolute offset dlam (measured - HSRS), weighted over windows; air and vacuum:')
    for ch in whn_srf.CHANNELS:
        g = ship[ship['channel'] == ch]
        w = 1 / g['dlam_err'] ** 2
        d_air = np.sum(w * g['dlam']) / w.sum()
        va = srf.air_to_vac(g['lam_center']) - g['lam_center']
        d_vac = np.sum(w * (g['dlam'] - va)) / w.sum()
        chi2 = np.sum(w * (g['dlam'] - d_air) ** 2) / (len(g) - 1)
        e = np.sqrt(max(chi2, 1) / w.sum())
        # linear trend of dlam with wavelength (air)
        A = np.column_stack([np.ones(len(g)), (g['lam_center'] - 600) / 100])
        cf = np.linalg.lstsq(A * np.sqrt(w.values)[:, None], g['dlam'] * np.sqrt(w), rcond=None)[0]
        med_air = np.median(g['dlam'])
        med_vac = np.median(g['dlam'] - va)
        mad = 1.4826 * np.median(np.abs(g['dlam'] - med_air))
        print('  %-2s air %+.3f ± %.3f nm (scatter chi2_nu %.1f); vacuum %+.3f; '
              'trend %+.3f nm per 100 nm; range %+.2f to %+.2f; vac-air %.2f-%.2f'
              % (ch, d_air, e, chi2, d_vac, cf[1], g['dlam'].min(), g['dlam'].max(),
                 va.min(), va.max()))
        print('      robust: median air %+.3f, median vacuum %+.3f, MAD-sigma %.3f, '
              'error of median ~%.3f nm (n = %d)'
              % (med_air, med_vac, mad, 1.25 * mad / np.sqrt(len(g)), len(g)))
    mE_ = ship[ship['channel'] == 'E'].set_index('name')
    for ch in ('Ld', 'Lu'):
        mL_ = ship[ship['channel'] == ch].set_index('name')
        d = mL_['dlam'] - mE_['dlam']
        e = np.hypot(mL_['dlam_err'], mE_['dlam_err'])
        w = 1 / e ** 2
        mean = np.nansum(w * d) / np.nansum(w)
        chi2 = np.nansum(w * (d - mean) ** 2) / (d.notna().sum() - 1)
        print('  relative %s - E: %+.3f ± %.3f nm (chi2_nu %.1f)' % (
            ch, mean, np.sqrt(max(chi2, 1) / np.nansum(w)), chi2))
    print('\nveil (free fits), median by channel:',
          tv[tv['ok']].groupby('channel')['veil'].median().round(3).to_dict(),
          ' median error:', tv[tv['ok']].groupby('channel')['veil_err'].median().round(3).to_dict())
    print('window chi2_nu, median by channel (veil free):',
          tv[tv['ok']].groupby('channel')['chi2_nu'].median().round(2).to_dict(),
          '(veil 0):', tf[tf['ok']].groupby('channel')['chi2_nu'].median().round(2).to_dict())
    print('failed windows:', tab.loc[~tab['ok'], ['channel', 'name', 'veil_fit']].values.tolist())
    return tab, t


def figure(res, tab, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, (ax, bx, cx) = plt.subplots(3, 1, figsize=(8.5, 10), sharex=True,
                                     gridspec_kw=dict(height_ratios=[3, 2, 2], hspace=0.1))
    lam = np.linspace(390, 880, 300)
    offs = {'E': -1.2, 'Ld': 0.0, 'Lu': 1.2}
    emp = {m.channel: m for m in res[False]['models'] if m.meta['method'] == 'empirical'}
    mods = {ch: model_of(res[SHIP_VEIL[ch]], 'template', ch) for ch in whn_srf.CHANNELS}
    for ch in whn_srf.CHANNELS:
        c, m = COLORS[ch], mods[ch]
        g = tab[(tab['channel'] == ch) & (tab['veil_fit'] == SHIP_VEIL[ch]) & tab['ok']]
        y, e = m.fwhm(lam), m.fwhm_err(lam)
        ax.fill_between(lam, y - e, y + e, color=c, alpha=0.15, lw=0)
        ax.plot(lam, y, color=c, lw=2, label='%s template (χ²ν %.1f%s)'
                % (ch, m.chi2_nu, ', veil free' if SHIP_VEIL[ch] else ''))
        ax.plot(lam, emp[ch].fwhm(lam), color=c, lw=1.2, ls='--')
        ax.errorbar(g['lam_center'] + offs[ch], g['fwhm'], yerr=g['fwhm_err'],
                    fmt=MARKERS[ch], ms=5, color=c, lw=1, alpha=0.8)
        cx.errorbar(g['lam_center'] + offs[ch], g['dlam'], yerr=g['dlam_err'],
                    fmt=MARKERS[ch], ms=5, color=c, lw=1, label=ch)
        if ch != 'E':
            gE = tab[(tab['channel'] == 'E') & (tab['veil_fit'] == SHIP_VEIL['E'])
                     & tab['ok']].set_index('name')
            h = g.set_index('name')
            gE = gE.reindex(h.index)
            d = h['fwhm'] - gE['fwhm']
            de = np.hypot(h['fwhm_err'], gE['fwhm_err'])
            bx.errorbar(h['lam_center'] + offs[ch], d, yerr=de, fmt=MARKERS[ch], ms=5,
                        color=c, lw=1, label='%s − E' % ch)
            bx.plot(lam, m.fwhm(lam) - mods['E'].fwhm(lam), color=c, lw=2)
    ax.plot([], [], color='#555', lw=1.2, ls='--', label='empirical (task 5), for reference')
    ax.set_ylabel('SRF FWHM (nm)')
    ax.set_title('VEIT 2026-06-04, HYPSTAR 122304: Gaussian SRF from the TSIS-1 HSRS '
                 'template fit (10 nm windows)', fontsize=10, loc='left')
    ax.legend(fontsize=8, frameon=False, ncol=2, loc='lower left')
    bx.axhline(0, color='#888', lw=0.8)
    bx.set_ylabel('FWHM_L − FWHM_E (nm)')
    bx.legend(fontsize=8, frameon=False, loc='upper right')
    va = srf.air_to_vac(lam) - lam
    cx.plot(lam, va, color='#888', lw=1.2, ls=':', label='expected if the scale were vacuum')
    cx.axhline(0, color='#888', lw=0.8)
    cx.set_ylabel('offset vs HSRS in air,\nmeasured − reference (nm)')
    cx.set_xlabel('wavelength, air (nm)')
    cx.legend(fontsize=8, frameon=False, ncol=4, loc='upper left')
    for a in (ax, bx, cx):
        a.grid(alpha=0.25, lw=0.6)
        a.spines[['top', 'right']].set_visible(False)
        a.tick_params(labelsize=9)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)


def main():
    pd.set_option('display.width', 220)
    os.makedirs(OUT, exist_ok=True)
    irr, rad, res = run()
    tab, _ = report(res)
    tab.to_parquet(os.path.join(OUT, 'veit_template.parquet'))
    cols = ['channel', 'veil_fit', 'name', 'lo', 'hi', 'lam_center', 'ok', 'npix', 'sigma',
            'sigma_err', 'fwhm', 'fwhm_err', 'dlam', 'dlam_err', 'veil', 'veil_err', 'c0',
            'c1', 'chi2_nu']
    tab[cols].to_csv(os.path.join(WIGGLES_DIR, 'phase0_veit_template.csv'), index=False,
                     float_format='%.7g')
    mods = [model_of(res[SHIP_VEIL[ch]], 'template', ch) for ch in whn_srf.CHANNELS]
    meta = dict(res[True]['meta'],
                method='TSIS-1 HSRS v2 (air) x Gaussian SRF, 10 nm windows 390-880 nm, '
                       'telluric bands excluded (srf.fit_template_windows)',
                veil={k: bool(v) for k, v in SHIP_VEIL.items()},
                errors='scan_errors flatten_px=%d' % whn_srf.FLATTEN_PX,
                offset='measured - HSRS(air); add (air_to_vac(lam) - lam) for vacuum',
                empirical_models='hypernet/data/veit_srf_model_empirical.json',
                script='hypernet/wiggles/phase0a_template_veit.py')
    srf.save_srf_models(os.path.join(DATA_DIR, 'veit_srf_model.json'),
                        mods, meta=meta)
    figure(res, tab, os.path.join(FIGDIR, 'veit_template_fwhm.png'))
    print('\nwrote veit_template.parquet, hypernet/wiggles/phase0_veit_template.csv, '
          'hypernet/data/veit_srf_model.json, hypernet/wiggles/figs/phase0/veit_template_fwhm.png')


if __name__ == '__main__':
    main()

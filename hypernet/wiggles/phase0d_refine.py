"""Phase 0, task 8d: template-fit refinements.

Parts (``--part`` selects one; default all):

- ``telluric`` (i): fill the 680-850 nm gap.  In 10 nm windows across the
  O2-B, H2O and O2-A bands, fit VEIT E and Ld (L1A means, flattened errors)
  against HSRS x exp(-k tau), where tau is the dominant absorber's vertical
  optical depth (``hypernet.emod``, HITRAN via HAPI) with k fitted.  The minor
  absorber is applied at the geometric air mass.  Then compare FWHM(lambda),
  refitted over all windows, with the shipped model.
- ``shape`` (ii): in the non-telluric windows, fit Gaussian, pseudo-Voigt and
  super-Gaussian SRFs (``srf.fit_srf_template(shape=...)``).  Compare chi2_nu,
  the shape parameters and the FWHMs.
- ``lured`` (iii): the additive red component in Lu (task 6).
  - VEIT: template fits of Lu and Ld with the veil free, giving the additive
    level in radiance units against wavelength.
  - VEIT and Release 2: Lu / Ld at 900-1000 nm, where Lw ~ 0 in clear water,
    against the Fresnel rho_F ~ 0.026.
- ``trend`` (iv): the wavelength-dependent centroid offset, the per-sequence
  slope of dlam against lambda for E and Ld over the 224 Release 2
  sequences (``release2_template.parquet``, task 8c), per instrument.
- ``ring``: the Ld FWHM-vs-sky slope of task 8e.  Refit Release 2 Ld with the
  veil free, then test whether the veil tracks the sky index and absorbs the
  slope.

Outputs: ``hypernet/wiggles/phase0d_*.csv`` and
``figs/phase0/phase0d_*.png``; parquet under ``$OS_COLOR``.

Run from the repository root:
``python -m hypernet.wiggles.phase0d_refine [--part P]``.
"""
import argparse
import os

import numpy as np
import pandas as pd

from hypernet import emod, srf, whn_srf
from hypernet import whn_l1a as wl
from hypernet.wiggles import DATA_DIR, FIGDIR, OUT, REPO, WIGGLES_DIR

VEIT = dict(site='VEIT', seq_time='20260604T0845')
REQUEST = os.path.join(REPO, 'docs', 'wiggles_data_request.csv')
SZA_VEIT = 36.7
#: Telluric windows (10 nm tiles across the bands excluded in task 3b).
TELLURIC_WINDOWS = [(680.0, 690.0), (690.0, 700.0), (700.0, 710.0), (710.0, 720.0),
                    (720.0, 730.0), (730.0, 740.0), (740.0, 750.0), (750.0, 760.0),
                    (760.0, 770.0), (770.0, 780.0), (780.0, 790.0), (790.0, 800.0),
                    (800.0, 810.0), (810.0, 820.0), (820.0, 830.0), (830.0, 840.0),
                    (840.0, 850.0)]


def veit_channels():
    f = wl.sequence_files(**VEIT)
    irr, rad = wl.load_l1a_irr(f['L1A_IRR']), wl.load_l1a_rad(f['L1A_RAD'])
    out = {}
    for ch, (w, sc) in whn_srf.channel_spectra(irr, rad).items():
        mean, e, _ = srf.scan_errors(sc, flatten_px=whn_srf.FLATTEN_PX)
        out[ch] = (w, mean, e)
    return out


# --- (i) telluric -----------------------------------------------------------------

def telluric(ref):
    rw, rf = ref
    m = (rw > 670) & (rw < 865)
    rw, rf = rw[m], rf[m]
    tau = {mol: emod.optical_depth(rw, mol) for mol in ('O2', 'H2O')}
    am = 1.0 / np.cos(np.radians(SZA_VEIT))
    ch = veit_channels()
    rows = []
    for lo, hi in TELLURIC_WINDOWS:
        wm = (rw >= lo) & (rw <= hi)
        dom = max(tau, key=lambda k: tau[k][wm].sum())
        minor = [k for k in tau if k != dom][0]
        rf_eff = rf * np.exp(-am * tau[minor])
        for c in ('E', 'Ld', 'Lu'):
            w, y, e = ch[c]
            r = srf.fit_srf_template(w, y, rw, rf_eff, lo, hi, err=e, veil=False,
                                     absorber=tau[dom])
            r.update(channel=c, absorber=dom, airmass_geom=am,
                     tau_peak=float(tau[dom][wm].max()))
            rows.append(r)
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(WIGGLES_DIR, 'phase0d_telluric_windows.csv'), index=False,
             float_format='%.5g')
    print('\n(i) telluric windows (VEIT, k = fitted air mass x column scale):')
    print(t[['channel', 'lo', 'hi', 'absorber', 'ok', 'fwhm', 'fwhm_err', 'dlam', 'tau_scale',
             'tau_scale_err', 'chi2_nu']].to_string(index=False, float_format='%.3f'))
    # refit FWHM(lambda) over 3b windows + good telluric windows
    base = pd.read_csv(os.path.join(WIGGLES_DIR, 'phase0_veit_template.csv'))
    base = base[~base['veil_fit'] & base['ok']]
    old = srf.load_srf_models(os.path.join(DATA_DIR, 'veit_srf_model.json'))
    cmp_rows, new_models = [], {}
    for c in ('E', 'Ld', 'Lu'):
        b = base[base['channel'] == c][['lam_center', 'fwhm', 'fwhm_err', 'dlam', 'dlam_err']]
        tt = t[(t['channel'] == c) & t['ok'] & (t['chi2_nu'] < 20)]
        allw = pd.concat([b, tt[['lam_center', 'fwhm', 'fwhm_err', 'dlam', 'dlam_err']]])
        df = allw.rename(columns={'lam_center': 'lam_air', 'dlam': 'dmu',
                                  'dlam_err': 'mu_err'}).assign(ok=True, use_for_srf=True,
                                                                blend=False)
        mdl = srf.SRFModel.from_lines(c, df.reset_index(drop=True), scale_cov=True,
                                      instrument=old[c].instrument,
                                      meta=dict(old[c].meta, method='template + telluric '
                                                '(8d)', telluric_windows=int(len(tt))))
        new_models[c] = mdl
        for lam in (650, 700, 750, 800, 850):
            d = mdl.fwhm(lam) - old[c].fwhm(lam)
            e = np.hypot(mdl.fwhm_err(lam), old[c].fwhm_err(lam))
            cmp_rows.append(dict(channel=c, lam=lam, fwhm_old=old[c].fwhm(lam),
                                 err_old=old[c].fwhm_err(lam), fwhm_new=mdl.fwhm(lam),
                                 err_new=mdl.fwhm_err(lam), diff=d, z=d / e,
                                 n_tell=len(tt)))
    cmp = pd.DataFrame(cmp_rows)
    cmp.to_csv(os.path.join(WIGGLES_DIR, 'phase0d_telluric_model_compare.csv'), index=False,
               float_format='%.4f')
    print('\nFWHM(lambda) with the telluric windows vs the shipped model:')
    print(cmp.to_string(index=False, float_format='%.3f'))
    return t, cmp, new_models


# --- (ii) shape -------------------------------------------------------------------

def shape(ref):
    rw, rf = ref
    ch = veit_channels()
    rows = []
    for c in ('E', 'Ld'):
        w, y, e = ch[c]
        for lo, hi in srf.template_windows():
            for sh in srf.SHAPES:
                r = srf.fit_srf_template(w, y, rw, rf, lo, hi, err=e, veil=False, shape=sh)
                r.update(channel=c)
                rows.append(r)
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(WIGGLES_DIR, 'phase0d_shape_windows.csv'), index=False,
             float_format='%.5g')
    piv = t[t['ok']].pivot_table(index=['channel', 'lo'], columns='shape',
                                 values=['chi2_nu', 'fwhm', 'shape_par', 'shape_par_err'])
    summ = []
    for c in ('E', 'Ld'):
        p = piv.loc[c]
        g = p['chi2_nu'].dropna()
        row = dict(channel=c, n=len(g))
        for sh in ('pvoigt', 'supergauss'):
            ratio = g[sh] / g['gauss']
            par = p['shape_par'][sh].dropna()
            perr = p['shape_par_err'][sh].dropna()
            w = 1 / perr.loc[par.index] ** 2
            row['%s_chi2_ratio_median' % sh] = ratio.median()
            row['%s_par_median' % sh] = par.median()
            row['%s_par_wmean' % sh] = float(np.sum(w * par) / w.sum())
            row['%s_par_wmean_err' % sh] = float(1 / np.sqrt(w.sum()))
            row['%s_dfwhm_median' % sh] = (p['fwhm'][sh] - p['fwhm']['gauss']).median()
        row['gauss_chi2_median'] = g['gauss'].median()
        summ.append(row)
    summ = pd.DataFrame(summ)
    summ.to_csv(os.path.join(WIGGLES_DIR, 'phase0d_shape_summary.csv'), index=False,
                float_format='%.4f')
    print('\n(ii) SRF shape (VEIT, %d windows/channel):' % len(srf.template_windows()))
    print(summ.set_index('channel').T.to_string(float_format=lambda v: '%.3f' % v))
    return t, summ


# --- (iii) Lu red -----------------------------------------------------------------

def lured(ref):
    rw, rf = ref
    ch = veit_channels()
    rows = []
    wins = [wn for wn in srf.template_windows() if wn[0] >= 550]
    for c in ('Lu', 'Ld'):
        w, y, e = ch[c]
        for lo, hi in wins:
            r = srf.fit_srf_template(w, y, rw, rf, lo, hi, err=e, veil=True)
            mm = (w >= lo) & (w <= hi)
            # additive level in radiance units: veil * <R> * continuum ~ veil * level / (1 + veil)
            lev = float(np.nanmean(y[mm]))
            r.update(channel=c, level=lev,
                     additive=lev * r['veil'] / (1 + r['veil']) if r['ok'] else np.nan)
            rows.append(r)
    t = pd.DataFrame(rows)
    # NIR ratio Lu / Ld (VEIT)
    wL, Lu, _ = ch['Lu']
    _, Ld, _ = ch['Ld']
    nir = []
    for lo, hi in ((900, 925), (925, 950), (950, 975), (975, 1000), (1000, 1050)):
        mm = (wL >= lo) & (wL < hi)
        nir.append(dict(source='VEIT L1A', lo=lo, hi=hi, lu=np.nanmean(Lu[mm]),
                        ld=np.nanmean(Ld[mm]), ratio=np.nanmean(Lu[mm]) / np.nanmean(Ld[mm])))
    # Release 2: Lu/Ld at 975-1000 nm per instrument
    req = pd.read_csv(REQUEST, dtype=str)
    r2 = []
    for _, r in req.iterrows():
        try:
            d = wl.load_l2b(wl.release2_path(r['site'], r['sequence_time'], r['file']))
        except Exception:
            continue
        w = d['wave']
        mm = (w >= 975) & (w < 1000)
        mm7 = (w >= 745) & (w < 755)
        lu, ld = d['upwelling_radiance'][:, 0], d['downwelling_radiance'][:, 0]
        r2.append(dict(site=r['site'][:-2], instrument=d['meta']['system_id'],
                       ratio_985=np.nanmean(lu[mm]) / np.nanmean(ld[mm]),
                       ratio_750=np.nanmean(lu[mm7]) / np.nanmean(ld[mm7]),
                       rhof=float(np.nanmean(d.get('rhof', [np.nan])))))
    r2 = pd.DataFrame(r2)
    by = r2.groupby(['instrument', 'site']).agg(n=('ratio_985', 'size'),
                                                ratio_985=('ratio_985', 'median'),
                                                ratio_750=('ratio_750', 'median'),
                                                rhof=('rhof', 'median')).reset_index()
    t.to_csv(os.path.join(WIGGLES_DIR, 'phase0d_lured_windows.csv'), index=False,
             float_format='%.5g')
    pd.DataFrame(nir).to_csv(os.path.join(WIGGLES_DIR, 'phase0d_lured_nir_veit.csv'),
                             index=False, float_format='%.5g')
    by.to_csv(os.path.join(WIGGLES_DIR, 'phase0d_lured_release2.csv'), index=False,
              float_format='%.5g')
    print('\n(iii) Lu vs Ld with the veil free (VEIT, >= 550 nm):')
    print(t[['channel', 'lo', 'ok', 'veil', 'veil_err', 'level', 'additive', 'fwhm']]
          .to_string(index=False, float_format='%.4f'))
    print('\nVEIT Lu/Ld in the NIR (Lw ~ 0; expect ~ rho_F):')
    print(pd.DataFrame(nir).to_string(index=False, float_format='%.4f'))
    print('\nRelease 2: median Lu/Ld at 975-1000 and 745-755 nm, and the processor rho_F:')
    print(by.to_string(index=False, float_format='%.4f'))
    return t, nir, by


# --- (iv) offset trend --------------------------------------------------------------

def trend():
    t = pd.read_parquet(os.path.join(OUT, 'release2_template.parquet'))
    t = t[t['ok'] & np.isfinite(t['dlam_err']) & (t['dlam_err'] > 0)]
    rows = []
    for (seq, inst, c), g in t.groupby(['sequence_time', 'system_id', 'channel']):
        if len(g) < 8:
            continue
        x = (g['lam_center'].values - 600) / 100
        w = 1 / g['dlam_err'].values ** 2
        A = np.column_stack([np.ones_like(x), x])
        cov = np.linalg.inv(A.T @ (A * w[:, None]))
        b = cov @ (A.T @ (w * g['dlam'].values))
        res = g['dlam'].values - A @ b
        chi2 = np.sum(w * res ** 2) / (len(g) - 2)
        rows.append(dict(sequence_time=seq, instrument=inst, channel=c, n=len(g),
                         offset600=b[0], slope_per100nm=b[1],
                         slope_err=np.sqrt(cov[1, 1] * max(chi2, 1)), chi2_nu=chi2))
    s = pd.DataFrame(rows)
    by = s.groupby(['instrument', 'channel']).agg(
        n_seq=('slope_per100nm', 'size'), slope_mean=('slope_per100nm', 'mean'),
        slope_std=('slope_per100nm', 'std'), slope_err_median=('slope_err', 'median'),
        offset600_mean=('offset600', 'mean')).reset_index()
    by['slope_sem'] = by['slope_std'] / np.sqrt(by['n_seq'])
    by['t'] = by['slope_mean'] / by['slope_sem']
    by['shift_390_to_870_nm'] = by['slope_mean'] * 4.8
    s.to_csv(os.path.join(OUT, 'phase0d_trend_per_sequence.csv'), index=False)
    by.to_csv(os.path.join(WIGGLES_DIR, 'phase0d_trend.csv'), index=False, float_format='%.4f')
    print('\n(iv) slope of the centroid offset vs wavelength (nm per 100 nm), Release 2:')
    print(by.to_string(index=False, float_format='%.4f'))
    return by


# --- Ring: Ld veil vs sky index ------------------------------------------------------

def ring(ref):
    req = pd.read_csv(REQUEST, dtype=str)
    sky = pd.read_csv(os.path.join(WIGGLES_DIR, 'phase0_sky_index.csv'),
                      dtype={'sequence_time': str})
    rows = []
    for i, r in req.iterrows():
        try:
            d = wl.load_l2b(wl.release2_path(r['site'], r['sequence_time'], r['file']))
        except Exception:
            continue
        w, Ld = d['wave'], d['downwelling_radiance'][:, 0]
        for veil in (False, True):
            t = srf.fit_template_windows(w, Ld, ref[0], ref[1], veil=veil)
            t = t[t['ok']]
            m = srf.SRFModel.from_lines('Ld', t, scale_cov=True)
            blue = t[t['lam_center'] < 560]
            rows.append(dict(site=r['site'], sequence_time=r['sequence_time'],
                             instrument=d['meta']['system_id'], veil_fit=veil,
                             fwhm_450=m.fwhm(450), fwhm_600=m.fwhm(600),
                             veil_blue_median=blue['veil'].median() if veil else 0.0,
                             veil_all_median=t['veil'].median() if veil else 0.0))
        if (i + 1) % 50 == 0:
            print('ring: %d/%d' % (i + 1, len(req)))
    s = pd.DataFrame(rows).merge(sky[['site', 'sequence_time', 'ld_ed_750']],
                                 on=['site', 'sequence_time'], how='left')
    s.to_csv(os.path.join(OUT, 'phase0d_ring_per_sequence.csv'), index=False)
    out = []
    for (inst, veil), g in s.groupby(['instrument', 'veil_fit']):
        g = g.dropna(subset=['ld_ed_750'])
        if len(g) < 6:
            continue
        r = dict(instrument=inst, veil_fit=veil, n=len(g))
        for q in ('fwhm_600', 'fwhm_450', 'veil_blue_median'):
            if not veil and q.startswith('veil'):
                continue
            x, y = g['ld_ed_750'].values, g[q].values
            A = np.column_stack([np.ones_like(x), x])
            c, *_ = np.linalg.lstsq(A, y, rcond=None)
            res = y - A @ c
            cov = np.sum(res ** 2) / (len(g) - 2) * np.linalg.inv(A.T @ A)
            r['slope_%s' % q] = c[1]
            r['slope_%s_err' % q] = np.sqrt(cov[1, 1])
        out.append(r)
    out = pd.DataFrame(out)
    out.to_csv(os.path.join(WIGGLES_DIR, 'phase0d_ring.csv'), index=False, float_format='%.4f')
    print('\nRing test: slopes vs the sky index ld_ed_750 (per unit index), veil fixed / free:')
    print(out.to_string(index=False, float_format='%.4f'))
    return s, out


def figure(tell, shp, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    base = pd.read_csv(os.path.join(WIGGLES_DIR, 'phase0_veit_template.csv'))
    base = base[~base['veil_fit'] & base['ok']]
    cols = {'E': '#2a78d6', 'Ld': '#eb6834', 'Lu': '#1baf7a'}
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(9, 7.5), sharex=True,
                                 gridspec_kw=dict(hspace=0.1))
    for c in ('E', 'Ld'):
        b = base[base['channel'] == c]
        ax.errorbar(b['lam_center'], b['fwhm'], yerr=b['fwhm_err'], fmt='o', ms=4,
                    color=cols[c], alpha=0.6, label='%s, task 3b windows' % c)
        t = tell[(tell['channel'] == c) & tell['ok'] & (tell['chi2_nu'] < 20)]
        ax.errorbar(t['lam_center'], t['fwhm'], yerr=t['fwhm_err'], fmt='D', ms=6,
                    mfc='white', color=cols[c], label='%s, telluric windows (HAPI)' % c)
        s = shp[(shp['channel'] == c) & shp['ok'] & (shp['shape'] == 'pvoigt')]
        bx.errorbar(s['lam_center'], s['shape_par'], yerr=s['shape_par_err'], fmt='o', ms=4,
                    color=cols[c], label='%s: pseudo-Voigt η' % c)
    ax.set_ylabel('SRF FWHM (nm)')
    ax.legend(fontsize=8, frameon=False, ncol=2)
    ax.set_title('VEIT: template SRF with the telluric windows filled, and SRF shape',
                 fontsize=10, loc='left')
    bx.axhline(0, color='#888', lw=0.8)
    bx.set_ylabel('Lorentzian fraction η\n(0 = Gaussian)')
    bx.set_xlabel('wavelength, air (nm)')
    bx.legend(fontsize=8, frameon=False)
    for a in (ax, bx):
        a.grid(alpha=0.25, lw=0.6)
        a.spines[['top', 'right']].set_visible(False)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description='Phase 0 task 8d')
    ap.add_argument('--part', default='all',
                    choices=['all', 'telluric', 'shape', 'lured', 'trend', 'ring'])
    a = ap.parse_args(argv)
    pd.set_option('display.width', 220)
    ref = whn_srf.load_reference('air')
    tell = shp = None
    if a.part in ('all', 'telluric'):
        tell, cmp, new = telluric(ref)
        os.makedirs(OUT, exist_ok=True)
        srf.save_srf_models(os.path.join(OUT, 'veit_srf_model_with_telluric.json'),
                            list(new.values()), meta={'note': 'candidate, task 8d(i)'})
    if a.part in ('all', 'shape'):
        shp, _ = shape(ref)
    if a.part in ('all', 'lured'):
        lured(ref)
    if a.part in ('all', 'trend'):
        trend()
    if a.part in ('all', 'ring'):
        ring(ref)
    ft = os.path.join(WIGGLES_DIR, 'phase0d_telluric_windows.csv')
    fs = os.path.join(WIGGLES_DIR, 'phase0d_shape_windows.csv')
    if a.part in ('all', 'telluric', 'shape') and os.path.exists(ft) and os.path.exists(fs):
        tell = pd.read_csv(ft) if tell is None else tell
        shp = pd.read_csv(fs) if shp is None else shp
        os.makedirs(FIGDIR, exist_ok=True)
        figure(tell, shp[shp['ok']], os.path.join(FIGDIR, 'phase0d_telluric_shape.png'))


if __name__ == '__main__':
    main()

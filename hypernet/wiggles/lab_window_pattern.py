"""The common window pattern vs the lab polynomial residual (lab data prompts,
task 4; Setup Q13).

**What is compared.**  The processor's wavelengths are the lab 4th-order
polynomials (``labcal.grid_for``), so wherever that polynomial misses the true
pixel -> lambda relation, the processor grid is in error.  At a lamp line,
``resid`` = poly(mu) - lambda_ref (``hypernet/data/lab_srf.csv``) is the
reported minus the true wavelength of the pixel that sees the line.  The field
pattern (``wavescale_window_pattern.csv``, Phase 2 task 4b) is in -dlam
(dlam = measured - reference), which is the same thing with the opposite sign
convention: a pixel reported at lambda samples lambda + e with e = -dlam.  A
reported-minus-true error r gives e = -r.  So the lab prediction of the
pattern is **-r**.  The field pattern is a residual about each sequence's
straight line (shift + stretch), so the prediction is also taken about a
straight line fitted over the same windows (``*_det``).

**The lab residual between lines** (Setup Q13).  Only one interpolation is
used, the transfer of a refit difference curve, in two variants:

- ``v1`` (Q13's default): "truth" = the 4th-order fit to *all* lines of the
  two Kr processor epochs (121222 and 121231 2021_07, 19 lines; their file
  polynomials are this fit to <= 0.014 nm).  A standard-set instrument (13
  lines, Hg/Ar + Xe) then has r(lambda) = poly_std - poly_all, the mean over
  the two instruments, per channel.  For the Kr instruments themselves r = 0
  between lines.
- ``v2`` (the Q13 lead): "truth" = the fit without Xe (Hg/Ar + Kr, 14 lines),
  which has a much smaller residual.  r = poly_std - poly_noXe for the
  standard instruments, and poly_all - poly_noXe for the Kr ones.

The transfer is justified only if the curves agree between instruments and
channels (``spread``), and is tested on 120242, which has a Kr-without-Xe
epoch (2020_09) and a standard one (2024_05).  Each Kr epoch, and 120242
2020_09, is refitted with all lines, without Kr, and without Xe; we report
the rms of the residual (and sqrt(sum r^2 / (n - 5)), the convention of a
"residual standard deviation") and the residual at each line.

Inputs: ``lab_srf.csv`` (our Gaussian line centres ``mu_px`` and ``resid``),
``lab_wl_coefs.csv``, ``wavescale_window_pattern.csv`` and the per-instrument
patterns ``hypernet/wiggles/phase2_wavescale_windows.csv``.

Outputs:

- ``hypernet/wiggles/lab_window_pattern.csv``: per window x instrument x
  channel, the empirical pattern, the lab predictions (raw and about the line)
  and the empirical pattern minus each;
- ``hypernet/wiggles/lab_window_pattern_summary.csv``: per instrument x
  channel (and the common pattern), rms before and after, and at the lines;
- ``hypernet/wiggles/figs/lab/window_pattern_vs_lab.png``;
- intermediates in ``$OS_COLOR/hypernet/wiggles/lab_data/window_pattern/``:
  ``kr_refits.csv`` (sigma and per-line residual of each refit),
  ``transfer_curves.csv``, ``transfer_test_120242.csv``, ``lines.csv``.

Run from the repository root:  ``python -m hypernet.wiggles.lab_window_pattern``
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from hypernet import labcal, srfcal
from hypernet.wiggles import DATA_DIR, WIGGLES_DIR

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, 'figs', 'lab')
OUT = os.path.join(labcal.lab_root(), 'window_pattern')
KR_EPOCHS = (('121222', '2021_07'), ('121231', '2021_07'))
TEST = ('120242', '2020_09', '2024_05')
CH = {'E': 'E', 'Ld': 'L'}
#: windows holding a lab line (window centre -> line)
LINE_WINDOWS = {405.0: 404.6563, 435.0: 435.8328, 545.0: 546.0735, 585.0: 587.0916,
                855.0: 850.88728}
#: the windows where the field pattern is largest (|pattern| ~ 0.2 nm)
BIG_WINDOWS = (415.0, 505.0, 605.0, 665.0)
SETS = {'all': lambda sp: np.ones(len(sp), bool),
        'no_Kr': lambda sp: sp != 'Kr I',
        'no_Xe': lambda sp: sp != 'Xe I'}


def _lines(sn, epoch, ch):
    s = labcal.load_lab_srf()
    g = s[(s['sn'] == sn) & (s['epoch'] == epoch) & (s['channel'] == ch)]
    return g[np.isfinite(g['mu_px'])].sort_values('lam_ref')


def refit(sn, epoch, ch, which):
    """4th-order fit of lambda_ref(mu_px) over a line subset; unweighted, as
    the lab files appear to be."""
    g = _lines(sn, epoch, ch)
    use = SETS[which](g['species'].values)
    c = np.polyfit(g['mu_px'][use], g['lam_ref'][use], 4)[::-1]
    r = labcal.poly(c, g['mu_px'].values) - g['lam_ref'].values
    return c, g.assign(resid_fit=r, used=use)


def _on_lambda(c_a, c_b, c_map, lam):
    """poly_a - poly_b at the pixels where c_map puts lam."""
    px = labcal.px_of(c_map, np.asarray(lam, float))
    return labcal.poly(c_a, px) - labcal.poly(c_b, px)


def kr_refits():
    rows, curves = [], []
    for sn, ep in KR_EPOCHS + ((TEST[0], TEST[1]),):
        for ch in labcal.CHANNELS:
            cf = labcal.coefs(sn, ch, ep)
            fits = {w: refit(sn, ep, ch, w) for w in SETS
                    if not (w == 'no_Xe' and sn == TEST[0])}       # 2020_09 has no Xe
            for w, (c, g) in fits.items():
                u = g[g['used']]
                d = dict(sn=sn, epoch=ep, channel=ch, lines=w, n=int(u.shape[0]),
                         rms=float(np.sqrt(np.mean(u['resid_fit'] ** 2))),
                         sigma=float(np.sqrt(np.sum(u['resid_fit'] ** 2) / (len(u) - 5))),
                         file_minus_fit_max=float(np.abs(
                             _on_lambda(cf, c, cf, np.arange(400, 901, 5))).max()))
                for _, l in g.iterrows():
                    d['r_%.0f' % l['lam_ref']] = float(l['resid_fit'])
                rows.append(d)
            # difference curves on the window grid and a fine grid
            lam = np.arange(390.0, 1001.0, 1.0)
            ref_all = fits['all'][0]
            if 'no_Kr' in fits:
                curves.append(pd.DataFrame(dict(sn=sn, epoch=ep, channel=ch, curve='std_minus_all',
                                                lam=lam, r=_on_lambda(fits['no_Kr'][0], ref_all,
                                                                      cf, lam))))
            if 'no_Xe' in fits:
                ref = fits['no_Xe'][0]
                curves.append(pd.DataFrame(dict(sn=sn, epoch=ep, channel=ch, curve='std_minus_noXe',
                                                lam=lam, r=_on_lambda(fits['no_Kr'][0], ref,
                                                                      cf, lam))))
                curves.append(pd.DataFrame(dict(sn=sn, epoch=ep, channel=ch, curve='all_minus_noXe',
                                                lam=lam, r=_on_lambda(ref_all, ref, cf, lam))))
    return pd.DataFrame(rows), pd.concat(curves, ignore_index=True)


def transfer(curves):
    """Mean (and spread) over the two Kr processor epochs, per channel."""
    k = curves[curves['sn'].isin([s for s, _ in KR_EPOCHS])]
    m = k.groupby(['curve', 'channel', 'lam'])['r'].agg(['mean', 'std', 'min', 'max']).reset_index()
    return m


def _curve(tr, name, ch, lam):
    g = tr[(tr['curve'] == name) & (tr['channel'] == ch)]
    return np.interp(lam, g['lam'], g['mean'])


def lab_r(sn, ch, lam, tr, variant):
    """The lab reported-minus-true error r(lambda) of sn's processor grid."""
    kr = sn in [s for s, _ in KR_EPOCHS]
    if variant == 'v1':
        return np.zeros_like(lam) if kr else _curve(tr, 'std_minus_all', ch, lam)
    return _curve(tr, 'all_minus_noXe' if kr else 'std_minus_noXe', ch, lam)


def _detrend(lam, y):
    x = (lam - srfcal.LAM0) / srfcal.HALF
    A = np.column_stack([np.ones_like(x), x])
    coef = np.linalg.lstsq(A, y, rcond=None)[0]
    return y - A @ coef


def transfer_test(tr):
    """120242: the transferred curves vs its own epochs.  (a) the 2024_05
    (standard) line residuals vs the v1 / v2 prediction at the lines; (b) the
    'drift' poly_2024_05 - poly_2020_09 at the same pixels vs std_minus_noXe."""
    sn, e0, e1 = TEST
    rows = []
    for ch in labcal.CHANNELS:
        g = _lines(sn, e1, ch)
        g = g[~g['low_snr'] & ~g['outside_grid']]
        for v in ('v1', 'v2'):
            p = lab_r(sn, ch, g['lam_ref'].values, tr, v)
            for (_, l), pv in zip(g.iterrows(), p):
                rows.append(dict(test='resid_2024_05_at_lines', channel=ch, variant=v,
                                 lam=l['lam_ref'], observed=l['resid'], predicted=pv))
        c0, c1 = labcal.coefs(sn, ch, e0), labcal.coefs(sn, ch, e1)
        lam = np.array([404.66, 435.83, 500.0, 546.07, 600.0, 700.0, 763.51, 823.16, 881.94])
        obs = _on_lambda(c1, c0, c1, lam)
        for v in ('v1', 'v2'):
            pr = _curve(tr, 'std_minus_all' if v == 'v1' else 'std_minus_noXe', ch, lam)
            for l, o, pv in zip(lam, obs, pr):
                rows.append(dict(test='poly_2024_05_minus_2020_09', channel=ch, variant=v,
                                 lam=l, observed=o, predicted=pv))
    return pd.DataFrame(rows)


def compare(tr):
    common = pd.read_csv(os.path.join(DATA_DIR, 'wavescale_window_pattern.csv'))
    per = pd.read_csv(os.path.join(WIGGLES_DIR, 'phase2_wavescale_windows.csv'))
    per = pd.concat([per, common.rename(columns={'pattern': 'resid_med'})
                     .assign(instrument='common', channel='both')[per.columns.intersection(
                         ['instrument', 'channel', 'lam_center', 'resid_med'])]],
                    ignore_index=True)
    rows, summ = [], []
    for (inst, ch), g in per.groupby(['instrument', 'channel']):
        g = g.sort_values('lam_center')
        lam = g['lam_center'].values.astype(float)
        emp = g['resid_med'].values.astype(float)
        sns = [inst[-6:]] if inst != 'common' else [s for s in labcal.INSTRUMENTS]
        chans = [CH[ch]] if ch in CH else list(labcal.CHANNELS)
        pred = {}
        for v in ('v1', 'v2'):
            # the common pattern is the median of the groups: compare with the median
            # of the group predictions
            p = np.median([-lab_r(s, c, lam, tr, v) for s in sns for c in chans], axis=0)
            pred[v] = (p, _detrend(lam, p))
        # lab line residuals (raw, no interpolation) at the line windows
        line_pred = np.full(lam.size, np.nan)
        for k, l in enumerate(lam):
            if l in LINE_WINDOWS:
                vals = []
                for s in sns:
                    for c in chans:
                        gg = _lines(s, labcal.PROCESSOR_EPOCH[s]['epoch'], c)
                        gg = gg[np.isclose(gg['lam_ref'], LINE_WINDOWS[l], atol=1e-3) &
                                ~gg['low_snr']]
                        vals += list(-gg['resid'])
                line_pred[k] = np.median(vals) if vals else np.nan
        d = pd.DataFrame(dict(instrument=inst, channel=ch, lam_center=lam, empirical=emp,
                              lab_line=line_pred,
                              lab_v1=pred['v1'][0], lab_v1_det=pred['v1'][1],
                              lab_v2=pred['v2'][0], lab_v2_det=pred['v2'][1]))
        d['after_v1'] = d['empirical'] - d['lab_v1_det']
        d['after_v2'] = d['empirical'] - d['lab_v2_det']
        rows.append(d)
        rms = lambda v: float(np.sqrt(np.nanmean(np.asarray(v) ** 2)))  # noqa: E731
        lw = d[np.isfinite(d['lab_line'])]
        s = dict(instrument=inst, channel=ch, n_windows=len(d), rms_before=rms(emp),
                 rms_after_v1=rms(d['after_v1']), rms_after_v2=rms(d['after_v2']),
                 corr_v1=float(np.corrcoef(emp, d['lab_v1_det'])[0, 1])
                 if np.std(d['lab_v1_det']) > 0 else np.nan,
                 corr_v2=float(np.corrcoef(emp, d['lab_v2_det'])[0, 1]),
                 n_line_windows=len(lw), rms_lines_before=rms(lw['empirical']),
                 rms_lines_after_raw=rms(lw['empirical'] - lw['lab_line']))
        for l in tuple(LINE_WINDOWS) + BIG_WINDOWS:
            r = d[d['lam_center'] == l]
            if len(r):
                s['emp_%d' % l] = float(r['empirical'].iloc[0])
                s['lab_%d' % l] = float(r['lab_line'].iloc[0])
                s['v2det_%d' % l] = float(r['lab_v2_det'].iloc[0])
        s['frac_removed_v1'] = 1 - s['rms_after_v1'] ** 2 / s['rms_before'] ** 2
        s['frac_removed_v2'] = 1 - s['rms_after_v2'] ** 2 / s['rms_before'] ** 2
        summ.append(s)
    return pd.concat(rows, ignore_index=True), pd.DataFrame(summ)


def figure(long, tr, refits, test):
    os.makedirs(FIG, exist_ok=True)
    fig, axs = plt.subplots(1, 3, figsize=(17, 5.2), gridspec_kw=dict(width_ratios=[1.3, 1, 1]))
    a = axs[0]
    c = long[long['instrument'] == 'common']
    a.plot(c['lam_center'], c['empirical'], 'k-o', ms=4, lw=1.6, label='field pattern (common)')
    a.plot(c['lam_center'], c['lab_v1_det'], '-', color='#2a78d6', lw=1.6,
           label='lab v1: std − all-line fit (Kr epochs)')
    a.plot(c['lam_center'], c['lab_v2_det'], '-', color='#eb6834', lw=1.6,
           label='lab v2: fit − no-Xe fit (Kr epochs)')
    lw = c[np.isfinite(c['lab_line'])]
    a.plot(lw['lam_center'], lw['lab_line'], 's', ms=8, mfc='none', mec='#1baf7a', mew=1.8,
           label='lab line residual, −r (median, no interpolation)')
    a.axhline(0, color='0.5', lw=0.8)
    a.set_xlabel('window centre (nm)')
    a.set_ylabel('−dlam = reported-wavelength error (nm)')
    a.set_title('(a) field pattern vs lab polynomial residual', fontsize=10, loc='left')
    a.legend(fontsize=7.5, frameon=False, loc='lower left')
    a = axs[1]
    cols = {'all': '#2a78d6', 'no_Kr': '#e34948', 'no_Xe': '#1baf7a'}
    for (sn, ep, w), g in refits[refits['channel'] == 'L'].groupby(['sn', 'epoch', 'lines']):
        rc = [k for k in g.columns if k.startswith('r_')]
        lam = np.array([float(k[2:]) for k in rc])
        o = np.argsort(lam)
        v = g[rc].values[0][o]
        ok = np.isfinite(v)
        mk = {'121222': 'o', '121231': 's', '120242': '^'}[sn]
        a.plot(lam[o][ok], v[ok], '-' + mk, ms=4, lw=1, color=cols[w], alpha=0.85,
               label='%s %s %s (rms %.3f)' % (sn, ep, w, g['rms'].iloc[0]))
    a.axhline(0, color='0.5', lw=0.8)
    a.set_xlim(280, 1100)
    a.set_ylim(-0.15, 0.15)
    a.text(1095, 0.14, 'no-Xe fits at Xe 1083.8\n(+0.57 to +0.69 nm, extrapolated)\noff scale', ha='right',
           va='top', fontsize=6.5, color='#1baf7a')
    a.set_xlabel('line (nm)')
    a.set_ylabel('4th-order fit residual, L (nm)')
    a.set_title('(b) refits of the Kr epochs: all / no Kr / no Xe', fontsize=10, loc='left')
    a.legend(fontsize=6, frameon=False, ncol=1)
    a = axs[2]
    t = test[(test['test'] == 'poly_2024_05_minus_2020_09') & (test['channel'] == 'L')]
    for v, col in (('v1', '#2a78d6'), ('v2', '#eb6834')):
        g = t[t['variant'] == v]
        a.plot(g['lam'], g['predicted'], '-', color=col, lw=1.6,
               label='%s transfer: %s' % (v, 'std − all' if v == 'v1' else 'std − no-Xe'))
    g = t[t['variant'] == 'v1']
    a.plot(g['lam'], g['observed'], 'ko', ms=5, label='120242: poly 2024_05 − poly 2020_09')
    a.axhline(0, color='0.5', lw=0.8)
    a.set_xlabel('wavelength (nm)')
    a.set_ylabel('nm (L, same pixel)')
    a.set_title('(c) transfer test on 120242', fontsize=10, loc='left')
    a.legend(fontsize=7.5, frameon=False)
    for a in axs:
        a.grid(alpha=0.25, lw=0.5)
        a.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    p = os.path.join(FIG, 'window_pattern_vs_lab.png')
    fig.savefig(p, dpi=140)
    plt.close(fig)
    return p


def main():
    os.makedirs(OUT, exist_ok=True)
    pd.set_option('display.width', 250)
    pd.set_option('display.max_columns', 40)
    refits, curves = kr_refits()
    tr = transfer(curves)
    test = transfer_test(tr)
    long, summ = compare(tr)
    refits.to_csv(os.path.join(OUT, 'kr_refits.csv'), index=False, float_format='%.4f')
    tr.to_csv(os.path.join(OUT, 'transfer_curves.csv'), index=False, float_format='%.4f')
    test.to_csv(os.path.join(OUT, 'transfer_test_120242.csv'), index=False, float_format='%.4f')
    long.to_csv(os.path.join(WIGGLES_DIR, 'lab_window_pattern.csv'), index=False,
                float_format='%.4f')
    summ.to_csv(os.path.join(WIGGLES_DIR, 'lab_window_pattern_summary.csv'), index=False,
                float_format='%.4f')
    print('Kr-epoch refits (rms / sigma of the residual over the lines used, nm):')
    print(refits[['sn', 'epoch', 'channel', 'lines', 'n', 'rms', 'sigma',
                  'file_minus_fit_max']].to_string(index=False, float_format='%.3f'))
    rc = ['r_405', 'r_436', 'r_546', 'r_587', 'r_764', 'r_785', 'r_823', 'r_851', 'r_882',
          'r_975', 'r_992', 'r_1014']
    print(refits[['sn', 'channel', 'lines'] + [c for c in rc if c in refits]]
          .to_string(index=False, float_format='%+.3f'))
    w = tr[tr['lam'].isin([405.0, 435.0, 505.0, 545.0, 565.0, 585.0, 605.0, 665.0, 855.0])]
    print('\ntransfer curves at window centres (mean over 121222/121231, std = spread):')
    print(w.pivot_table(index='lam', columns=['curve', 'channel'], values=['mean', 'std'])
          .round(3).to_string())
    print('\ntransfer test, 120242:')
    print(test.pivot_table(index=['test', 'channel', 'lam'], columns='variant',
                           values=['observed', 'predicted']).round(3).to_string())
    for (tt, v), g in test.groupby(['test', 'variant']):
        print('  %-28s %s: rms(observed) %.3f, rms(observed - predicted) %.3f'
              % (tt, v, np.sqrt(np.mean(g['observed'] ** 2)),
                 np.sqrt(np.mean((g['observed'] - g['predicted']) ** 2))))
    print('\npattern rms before / after (nm), per instrument x channel:')
    print(summ[['instrument', 'channel', 'n_windows', 'rms_before', 'rms_after_v1', 'rms_after_v2',
                'frac_removed_v1', 'frac_removed_v2', 'corr_v1', 'corr_v2', 'rms_lines_before',
                'rms_lines_after_raw']].to_string(index=False, float_format='%.3f'))
    print('\nat the line windows (field vs lab -r, raw; and the v2 curve about its line):')
    cols = [c for c in summ.columns if c[:4] in ('emp_', 'lab_', 'v2de')]
    print(summ[['instrument', 'channel'] + sorted(cols, key=lambda c: (c.split('_')[-1], c))]
          .to_string(index=False, float_format='%+.3f'))
    print('wrote', figure(long, tr, refits, test))


if __name__ == '__main__':
    main()

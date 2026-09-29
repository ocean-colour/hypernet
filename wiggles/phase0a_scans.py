"""Phase 0a, task 6: per-scan line fits, the Ld = Lu check, and line depths.

For the VEIT sample:

1. Fit every line in every single scan (6 E, 6 Ld, 6 Lu), each scan weighted
   by the flattened single-scan error (``srf.scan_errors(..., flatten_px=20)``,
   Q&A Q9).  Per line and channel, compare the scan-to-scan standard
   deviation of the centroid and FWHM with the median ``curve_fit`` error:
   a ratio near 1 means the errors are right.
2. Ld vs its two sky series (before / after the water view): do the widths
   change over the ~50 s between them while the sky brightness changes ~5 %?
3. Ld = Lu per line, from the mean-spectrum fits of task 5: the difference
   in FWHM and centroid over the combined error.  Two errors are used: the
   ``curve_fit`` error of the mean fit, and the scan scatter / sqrt(N); the
   larger is the conservative one.
4. Depths and equivalent widths (EW = depth x sigma x sqrt(2 pi), the
   weak-line approximation) in E, Ld and Lu, as a diagnostic only.  A wider
   SRF lowers the depth but conserves EW; Ring filling-in lowers both.  So
   EW_L / EW_E < 1 at the Fraunhofer lines points at filling-in (or at
   continuum placement), not at the SRF.

Outputs: ``$OS_COLOR/hypernet/wiggles/phase0/veit_scan_fits.parquet`` (every
per-scan fit), ``wiggles/phase0_veit_scans.csv`` (per line and channel
summary), and ``wiggles/figs/phase0/veit_scan_scatter.png``.

Run from the repository root: ``python wiggles/phase0a_scans.py``.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
from hypernet import srf  # noqa: E402
import phase0a_veit as p5  # noqa: E402  (load(), paths, colours)

CHANNELS = ('E', 'Ld', 'Lu')


def per_scan_fits(chans, rad):
    """Fit every line in every scan; also the two Ld series means."""
    rows = []
    for ch, (w, scans) in chans.items():
        _, _, e_scan = srf.scan_errors(scans, flatten_px=p5.FLATTEN_PX)
        for k in range(scans.shape[1]):
            df = srf.fit_lines(w, scans[:, k], err=e_scan)
            df.insert(0, 'scan', k)
            df.insert(0, 'channel', ch)
            rows.append(df)
    # the two sky series, each a mean of 3 scans, error = e_scan / sqrt(3)
    Ld = rad['Ld']
    _, _, e_scan = srf.scan_errors(Ld['scans'], flatten_px=p5.FLATTEN_PX)
    for sid in np.unique(Ld['series_id']):
        sel = Ld['series_id'] == sid
        df = srf.fit_lines(rad['wave'], Ld['scans'][:, sel].mean(axis=1),
                           err=e_scan / np.sqrt(sel.sum()))
        df.insert(0, 'scan', -int(sid))       # negative = a series mean
        df.insert(0, 'channel', 'Ld_series')
        rows.append(df)
    return pd.concat(rows, ignore_index=True)


def scatter_summary(scans, mean_fits):
    """Per line and channel: scan spread vs curve_fit errors."""
    s = scans[scans['channel'].isin(CHANNELS)]
    out = []
    for (ch, name), g in s.groupby(['channel', 'name'], sort=False):
        ok = g[g['ok']]
        n = len(ok)
        m = mean_fits[(mean_fits['channel'] == ch) & (mean_fits['name'] == name)].iloc[0]
        row = dict(channel=ch, name=name, lam_air=m['lam_air'], blend=m['blend'],
                   use_for_srf=m['use_for_srf'], n_ok=n)
        for q in ('mu', 'fwhm', 'depth'):
            std = ok[q].std(ddof=1) if n > 1 else np.nan
            med_err = ok[q + '_err'].median() if n else np.nan
            row[q + '_scan_std'] = std
            row[q + '_scan_err_med'] = med_err
            row[q + '_ratio'] = std / med_err if n > 1 else np.nan
            # the mean-spectrum fit (task 5) and the error of the mean from
            # the scan scatter
            row[q + '_mean'] = m[q]
            row[q + '_mean_err_fit'] = m[q + '_err']
            row[q + '_mean_err_scan'] = std / np.sqrt(n) if n > 1 else np.nan
        row['chi2_nu_scan_med'] = ok['chi2_nu'].median() if n else np.nan
        row['sigma_mean'] = m['sigma']
        row['sigma_mean_err'] = m['sigma_err']
        out.append(row)
    return pd.DataFrame(out)


def ld_lu_test(summ):
    """Ld - Lu per line, over the fit and the scan-scatter errors."""
    a = summ[summ['channel'] == 'Ld'].set_index('name')
    b = summ[summ['channel'] == 'Lu'].set_index('name')
    out = pd.DataFrame(index=a.index)
    out['lam_air'] = a['lam_air']
    out['blend'] = a['blend']
    out['use_for_srf'] = a['use_for_srf']
    for q in ('fwhm', 'mu'):
        d = a[q + '_mean'] - b[q + '_mean']
        e_fit = np.hypot(a[q + '_mean_err_fit'], b[q + '_mean_err_fit'])
        e_scn = np.hypot(a[q + '_mean_err_scan'], b[q + '_mean_err_scan'])
        e_max = np.fmax(e_fit, e_scn)
        out['d_' + q] = d
        out['err_fit_' + q] = e_fit
        out['err_scan_' + q] = e_scn
        out['z_' + q] = d / e_max        # conservative
    return out.reset_index()


def depths(summ):
    """Depth and weak-line equivalent width per line and channel."""
    t = summ[['channel', 'name', 'lam_air', 'blend', 'depth_mean', 'sigma_mean']].copy()
    t['ew'] = t['depth_mean'] * t['sigma_mean'] * np.sqrt(2 * np.pi)
    piv_d = t.pivot(index='name', columns='channel', values='depth_mean')
    piv_w = t.pivot(index='name', columns='channel', values='ew')
    out = pd.DataFrame({'lam_air': t.drop_duplicates('name').set_index('name')['lam_air']})
    for ch in CHANNELS:
        out['depth_' + ch] = piv_d[ch]
        out['ew_' + ch] = piv_w[ch]
    for ch in ('Ld', 'Lu'):
        out['depth_%s_over_E' % ch] = piv_d[ch] / piv_d['E']
        out['ew_%s_over_E' % ch] = piv_w[ch] / piv_w['E']
    return out.loc[list(srf.LINES['name'])].reset_index()


def series_check(scans, summ):
    """Ld series A vs B (FWHM, centroid) vs the all-scan Ld fit."""
    s = scans[scans['channel'] == 'Ld_series']
    ids = sorted(s['scan'].unique(), reverse=True)     # -4, -10 -> series 4, 10
    a = s[s['scan'] == ids[0]].set_index('name')
    b = s[s['scan'] == ids[1]].set_index('name')
    out = pd.DataFrame({'lam_air': a['lam_air'],
                        'fwhm_A': a['fwhm'], 'fwhm_B': b['fwhm'],
                        'd_fwhm': b['fwhm'] - a['fwhm'],
                        'err_d_fwhm': np.hypot(a['fwhm_err'], b['fwhm_err']),
                        'd_mu': b['mu'] - a['mu'],
                        'err_d_mu': np.hypot(a['mu_err'], b['mu_err']),
                        'depth_A': a['depth'], 'depth_B': b['depth']})
    out['z_fwhm'] = out['d_fwhm'] / out['err_d_fwhm']
    out['z_mu'] = out['d_mu'] / out['err_d_mu']
    return out.loc[list(srf.LINES['name'])].reset_index()


def figure(summ, ldlu, dep, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    names = list(srf.LINES['name'])
    x = np.arange(len(names))
    offs = {'E': -0.22, 'Ld': 0.0, 'Lu': 0.22}
    fig, axs = plt.subplots(4, 1, figsize=(9, 11), sharex=True,
                            gridspec_kw=dict(hspace=0.12))
    # (a), (b): scan std / median fit error, FWHM and centroid
    for ax, q, lab in ((axs[0], 'fwhm', 'FWHM'), (axs[1], 'mu', 'centroid')):
        for ch in CHANNELS:
            g = summ[summ['channel'] == ch].set_index('name').loc[names]
            ax.plot(x + offs[ch], g[q + '_ratio'], ls='', marker=p5.MARKERS[ch],
                    ms=7, color=p5.COLORS[ch], label=ch)
        ax.axhline(1, color='#888', lw=0.8)
        ax.set_yscale('log')
        ax.set_yticks([0.25, 0.5, 1, 2, 4])
        ax.set_yticklabels(['0.25', '0.5', '1', '2', '4'])
        ax.minorticks_off()
        ax.set_ylabel('%s: scan std /\nmedian fit error' % lab)
    axs[0].legend(fontsize=8.5, frameon=False, ncol=3, loc='upper left')
    axs[0].set_title('VEIT 2026-06-04, HYPSTAR 122304: per-scan fits (6 scans per channel)',
                     fontsize=10.5, loc='left')
    # (c): Ld - Lu FWHM with the conservative error
    t = ldlu.set_index('name').loc[names]
    err = np.fmax(t['err_fit_fwhm'], t['err_scan_fwhm'])
    fill = np.where(t['use_for_srf'], p5.COLORS['Ld'], 'white')
    for i in range(len(names)):
        axs[2].errorbar(x[i], t['d_fwhm'].iloc[i], yerr=err.iloc[i], fmt='D', ms=6,
                        mfc=fill[i], mec='#444', ecolor='#444', lw=1.2)
    axs[2].axhline(0, color='#888', lw=0.8)
    axs[2].set_ylabel('FWHM_Ld − FWHM_Lu\n(nm)')
    axs[2].text(0.995, 0.95, 'open = not in the FWHM(λ) fit; error = max(fit, scan)',
                transform=axs[2].transAxes, ha='right', va='top', fontsize=7.5,
                color='#666')
    # (d): EW ratio L / E
    d = dep.set_index('name').loc[names]
    for ch in ('Ld', 'Lu'):
        axs[3].plot(x + offs[ch], d['ew_%s_over_E' % ch], ls='', marker=p5.MARKERS[ch],
                    ms=7, color=p5.COLORS[ch], label='%s / E' % ch)
        axs[3].plot(x + offs[ch], d['depth_%s_over_E' % ch], ls='',
                    marker=p5.MARKERS[ch], ms=7, mfc='none', mec=p5.COLORS[ch])
    axs[3].axhline(1, color='#888', lw=0.8)
    axs[3].set_ylabel('L / E\n(diagnostic)')
    axs[3].legend(fontsize=8.5, frameon=False, ncol=2, loc='lower left')
    axs[3].text(0.995, 0.05, 'filled = equivalent width, open = depth',
                transform=axs[3].transAxes, ha='right', va='bottom', fontsize=7.5,
                color='#666')
    axs[3].set_xticks(x)
    axs[3].set_xticklabels(names, rotation=40, ha='right', fontsize=8.5)
    for ax in axs:
        ax.grid(alpha=0.25, lw=0.6)
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(labelsize=9)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)


def main():
    pd.set_option('display.width', 200)
    irr, rad, chans = p5.load()
    lines = p5.fit_all(chans)
    mean_fits = lines[lines['err_kind'] == 'flat']
    scans = per_scan_fits(chans, rad)
    scans.to_parquet(os.path.join(p5.OUT, 'veit_scan_fits.parquet'))
    summ = scatter_summary(scans, mean_fits)
    summ.to_csv(os.path.join(p5.REPO, 'wiggles', 'phase0_veit_scans.csv'),
                index=False, float_format='%.7g')
    ff = '%.3f'
    print('\n[1] scan std / median curve_fit error (1 = errors right), and '
          'per-scan chi2_nu:')
    for q in ('fwhm', 'mu'):
        print('  ' + q)
        print(summ.pivot(index='name', columns='channel', values=q + '_ratio')
              .loc[list(srf.LINES['name'])].to_string(float_format='%.2f'))
    print('  median per-scan chi2_nu')
    print(summ.pivot(index='name', columns='channel', values='chi2_nu_scan_med')
          .loc[list(srf.LINES['name'])].to_string(float_format='%.2f'))
    print('  failed per-scan fits (n_ok < 6):')
    print(summ.loc[summ['n_ok'] < 6, ['channel', 'name', 'n_ok']].to_string(index=False))
    print('\n  FWHM of the mean spectrum: fit error vs scan scatter / sqrt(N):')
    print(summ[['channel', 'name', 'fwhm_mean', 'fwhm_mean_err_fit',
                'fwhm_mean_err_scan']].to_string(index=False, float_format=ff))

    ser = series_check(scans, summ)
    print('\n[2] Ld series A (before water) vs B (after): B - A')
    print(ser[['name', 'fwhm_A', 'fwhm_B', 'd_fwhm', 'err_d_fwhm', 'z_fwhm', 'd_mu',
               'err_d_mu', 'z_mu']].to_string(index=False, float_format=ff))

    ldlu = ld_lu_test(summ)
    print('\n[3] Ld - Lu per line (z uses the larger of the fit and scan errors):')
    print(ldlu.to_string(index=False, float_format=ff))
    use = ldlu['use_for_srf']
    for q in ('fwhm', 'mu'):
        z = ldlu.loc[use, 'z_' + q]
        e = np.fmax(ldlu.loc[use, 'err_fit_' + q], ldlu.loc[use, 'err_scan_' + q])
        w = 1 / e ** 2
        mean = np.sum(w * ldlu.loc[use, 'd_' + q]) / w.sum()
        print('  %-4s over the %d SRF lines: chi2 = %.1f for %d dof; |z| > 2 at %d; '
              'weighted mean Ld - Lu = %+.3f +- %.3f nm'
              % (q, use.sum(), np.sum(z ** 2), use.sum(), np.sum(np.abs(z) > 2),
                 mean, 1 / np.sqrt(w.sum())))

    dep = depths(summ)
    print('\n[4] depths and weak-line equivalent widths (diagnostic only):')
    print(dep.to_string(index=False, float_format=ff))

    figure(summ, ldlu, dep, os.path.join(p5.FIGDIR, 'veit_scan_scatter.png'))
    print('\nwrote veit_scan_fits.parquet, wiggles/phase0_veit_scans.csv, '
          'wiggles/figs/phase0/veit_scan_scatter.png')


if __name__ == '__main__':
    main()

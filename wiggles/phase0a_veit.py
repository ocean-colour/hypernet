"""Phase 0a, task 5: line fits and FWHM(lambda) SRF models for the VEIT sample.

Fits every line of :data:`hypernet.srf.LINES` in

- E  -- the mean of the 6 L1A_IRR scans (E grid),
- Ld -- the mean of the 6 sky scans of L1A_RAD (vza >= 90, L grid),
- Lu -- the mean of the 6 water-view scans (vza < 90, L grid),

weighted by the flattened scan-to-scan error of the mean
(``srf.scan_errors(..., flatten_px=20)``, Q&A Q9).  The same fits with the raw
scatter as error are kept alongside (``err_kind = 'raw'``).  Then fits a
quadratic FWHM(lambda) per channel (``SRFModel.from_lines``, covariance
inflated by chi2_nu when > 1).

Outputs:

- ``$OS_COLOR/hypernet/wiggles/phase0/veit_lines.parquet`` and the committed
  copy ``wiggles/phase0_veit_lines.csv``;
- ``hypernet/data/veit_srf_model.json`` (E, Ld, Lu; the empirical models --
  to be superseded by the template fit of task 3b);
- ``wiggles/figs/phase0/veit_fwhm_vs_lambda.png``;
- on stdout, the plan §2.3 table re-derived and compared.

Run from the repository root: ``python wiggles/phase0a_veit.py``.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from hypernet import srf  # noqa: E402
from hypernet import whn_l1a as wl  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
OUT = os.path.join(os.environ['OS_COLOR'], 'hypernet', 'wiggles', 'phase0')
FIGDIR = os.path.join(REPO, 'wiggles', 'figs', 'phase0')
FLATTEN_PX = 20

# The plan §2.3 table as printed.  Its "Ld" column came from the vza < 90
# scans, i.e. it is really Lu, and its "Lu" column is really Ld (Q&A Setup #1);
# the centroid column was (vza < 90) - E, i.e. Lu - E.
PLAN = pd.DataFrame(
    [('Ca K', 1.83, 2.47, 2.46, +0.20), ('Ca H', 1.91, 2.81, 2.73, -0.14),
     ('G band', 2.52, 3.26, 3.13, +0.03), ('H beta', 2.97, 3.44, 3.47, +0.05),
     ('Mg b', 3.06, 3.65, 3.51, -0.04), ('Na D', 2.63, 3.11, 3.14, 0.00),
     ('H alpha', 2.77, 2.96, 2.94, +0.04), ('O2 B', 3.05, 3.05, 3.14, +0.08),
     ('O2 A', 3.40, 3.46, 3.42, +0.06)],
    columns=['name', 'E', 'Lu', 'Ld', 'dmu_Lu_E'])  # relabelled columns

# Palette slots 1-3 (dataviz reference palette, light), fixed order.
COLORS = {'E': '#2a78d6', 'Ld': '#eb6834', 'Lu': '#1baf7a'}
MARKERS = {'E': 'o', 'Ld': 's', 'Lu': '^'}


def load():
    files = wl.sequence_files(site='VEIT', seq_time='20260604T0845')
    irr = wl.load_l1a_irr(files['L1A_IRR'])
    rad = wl.load_l1a_rad(files['L1A_RAD'])
    chans = {'E': (irr['wave'], irr['scans']),
             'Ld': (rad['wave'], rad['Ld']['scans']),
             'Lu': (rad['wave'], rad['Lu']['scans'])}
    return irr, rad, chans


def fit_all(chans):
    rows = []
    for ch, (w, scans) in chans.items():
        mean, e_flat, _ = srf.scan_errors(scans, flatten_px=FLATTEN_PX)
        _, e_raw, _ = srf.scan_errors(scans)
        for kind, e in (('flat', e_flat), ('raw', e_raw)):
            df = srf.fit_lines(w, mean, err=e)
            df.insert(0, 'err_kind', kind)
            df.insert(0, 'channel', ch)
            rows.append(df)
        b = (w > 400) & (w < 900)
        print('%-2s median relative error of the mean, 400-900 nm: raw %.2e, '
              'flattened %.2e (x%.1f)' % (ch, np.nanmedian(e_raw[b] / mean[b]),
                                          np.nanmedian(e_flat[b] / mean[b]),
                                          np.nanmedian(e_raw[b] / e_flat[b])))
    return pd.concat(rows, ignore_index=True)


def models(lines, meta):
    out = []
    for ch in ('E', 'Ld', 'Lu'):
        f = lines[(lines['channel'] == ch) & (lines['err_kind'] == 'flat')]
        out.append(srf.SRFModel.from_lines(ch, f.reset_index(drop=True),
                                           scale_cov=True,
                                           instrument=meta['system_id'],
                                           meta=dict(meta, channel=ch)))
    return out


def compare_plan(lines):
    f = lines[lines['err_kind'] == 'flat']
    piv = f.pivot(index='name', columns='channel', values='fwhm')
    perr = f.pivot(index='name', columns='channel', values='fwhm_err')
    mu = f.pivot(index='name', columns='channel', values='mu')
    muerr = f.pivot(index='name', columns='channel', values='mu_err')
    order = list(srf.LINES['name'])
    print('\nRe-derived §2.3 table (FWHM, nm; flattened-error weights).  '
          'Columns correctly labelled; plan values (relabelled) in brackets.')
    print('%-12s %7s  %-19s %-19s %-19s %-16s %-16s' % (
        'line', 'lam', 'FWHM E', 'FWHM Ld (sky)', 'FWHM Lu (water)',
        'Lu - E (nm)', 'Ld - E (nm)'))
    diffs = []
    for name in order:
        lam = srf.LINES.set_index('name').loc[name, 'lam_air']
        p = PLAN.set_index('name').loc[name] if name in set(PLAN['name']) else None

        def cell(ch):
            v = piv.loc[name, ch]
            s = '%5.2f±%.2f' % (v, perr.loc[name, ch]) if np.isfinite(v) else '  -  '
            if p is not None:
                s += ' [%4.2f]' % p[ch]
                diffs.append((name, ch, v - p[ch]))
            return s

        def dcell(ch):
            d = mu.loc[name, ch] - mu.loc[name, 'E']
            e = np.hypot(muerr.loc[name, ch], muerr.loc[name, 'E'])
            s = '%+5.2f±%.2f' % (d, e) if np.isfinite(d) else '  -  '
            if p is not None and ch == 'Lu':
                s += ' [%+4.2f]' % p['dmu_Lu_E']
                diffs.append((name, 'dmu_Lu_E', d - p['dmu_Lu_E']))
            return s

        print('%-12s %7.2f  %-19s %-19s %-19s %-16s %-16s' % (
            name, lam, cell('E'), cell('Ld'), cell('Lu'), dcell('Lu'), dcell('Ld')))
    d = pd.DataFrame(diffs, columns=['name', 'quantity', 'new_minus_plan'])
    print('\nDifferences from the plan table (new - plan), |d| > 0.02 nm:')
    big = d[np.abs(d['new_minus_plan']) > 0.02]
    print(big.to_string(index=False, float_format='%+.2f') if len(big) else '  none')
    return d


def absolute_offsets(lines):
    """Measured centroid minus laboratory wavelength, air and vacuum, per
    clean (unblended) line: which scale does the HYPSTAR calibration follow?"""
    f = lines[(lines['err_kind'] == 'flat') & ~lines['blend'] & lines['ok']].copy()
    f['dmu_vac'] = f['mu'] - f['lam_vac']
    f['vac_minus_air'] = f['lam_vac'] - f['lam_air']
    print('\nAbsolute centroid offsets of the clean lines (nm):')
    print(f[['channel', 'name', 'mu', 'mu_err', 'dmu', 'dmu_vac', 'vac_minus_air']]
          .rename(columns={'dmu': 'mu-air', 'dmu_vac': 'mu-vac'})
          .to_string(index=False, float_format='%.3f'))
    for ch, g in f.groupby('channel'):
        w = 1 / g['mu_err'] ** 2
        print('  %-2s weighted mean: mu-air %+.3f, mu-vac %+.3f' % (
            ch, np.sum(w * g['dmu']) / w.sum(), np.sum(w * g['dmu_vac']) / w.sum()))


def figure(lines, mods, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    f = lines[lines['err_kind'] == 'flat']
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(8, 7), sharex=True,
                                 gridspec_kw=dict(height_ratios=[3, 2], hspace=0.08))
    lam = np.linspace(385, 875, 300)
    offs = {'E': -2.0, 'Ld': 0.0, 'Lu': 2.0}   # nm, to separate error bars
    for m in mods:
        ch, c = m.channel, COLORS[m.channel]
        g = f[f['channel'] == ch]
        y, e = m.fwhm(lam), m.fwhm_err(lam)
        ax.fill_between(lam, y - e, y + e, color=c, alpha=0.15, lw=0)
        ax.plot(lam, y, color=c, lw=2, label='%s  (χ²ν = %.1f)' % (ch, m.chi2_nu))
        for sel, face, lab in (
                (g['use_for_srf'] & ~g['blend'], c, 'clean'),
                (g['use_for_srf'] & g['blend'], 'white', 'blend'),
                (~g['use_for_srf'], 'none', 'diagnostic')):
            h = g[sel & g['ok']]
            ax.errorbar(h['lam_air'] + offs[ch], h['fwhm'], yerr=h['fwhm_err'],
                        fmt=MARKERS[ch], ms=7, mfc=face, mec=c, ecolor=c,
                        alpha=0.45 if lab == 'diagnostic' else 1.0, lw=1.2,
                        capsize=0)
        if ch != 'E':
            # FWHM_L - FWHM_E per line, errors in quadrature
            gE = f[f['channel'] == 'E'].set_index('name')
            h = g.set_index('name')
            d = h['fwhm'] - gE['fwhm']
            de = np.hypot(h['fwhm_err'], gE['fwhm_err'])
            use = h['ok'] & gE['ok'] & h['use_for_srf']
            bx.errorbar(h.loc[use, 'lam_air'] + offs[ch], d[use], yerr=de[use],
                        fmt=MARKERS[ch], ms=7, mfc=c, mec=c, ecolor=c, lw=1.2,
                        label='%s − E' % ch)
            mE = [x for x in mods if x.channel == 'E'][0]
            bx.plot(lam, m.fwhm(lam) - mE.fwhm(lam), color=c, lw=2)
    ax.set_ylabel('Gaussian FWHM (nm)')
    fig.suptitle('VEIT 2026-06-04 08:45, HYPSTAR 122304: line widths and quadratic FWHM(λ)',
                 fontsize=10.5, x=0.08, ha='left', y=0.995)
    # legend: channels, then marker fills
    from matplotlib.lines import Line2D
    h1, l1 = ax.get_legend_handles_labels()
    extra = [Line2D([], [], ls='', marker='o', mfc='#555', mec='#555', ms=7),
             Line2D([], [], ls='', marker='o', mfc='white', mec='#555', ms=7),
             Line2D([], [], ls='', marker='o', mfc='none', mec='#555', ms=7, alpha=0.45)]
    ax.legend(h1 + extra, l1 + ['clean line', 'blend / band (in fit)',
                                'diagnostic (not in fit)'],
              fontsize=8.5, frameon=False, ncol=3, loc='lower left',
              bbox_to_anchor=(0.0, 1.0, 1.0, 0.1), mode='expand')
    ax.set_xlim(375, 885)
    ax.set_ylim(1.6, 4.0)
    ax.text(0.995, 0.97, 'H₂O 936 nm (diagnostic) off-plot', transform=ax.transAxes,
            ha='right', va='top', fontsize=7.5, color='#666')
    bx.axhline(0, color='#888', lw=0.8)
    bx.set_ylabel('FWHM_L − FWHM_E (nm)')
    bx.set_xlabel('wavelength, air (nm)')
    bx.legend(fontsize=8.5, frameon=False, loc='upper right')
    for a in (ax, bx):
        a.grid(alpha=0.25, lw=0.6)
        a.spines[['top', 'right']].set_visible(False)
        a.tick_params(labelsize=9)
    for name, lam0 in zip(srf.LINES['name'], srf.LINES['lam_air']):
        if name in ('Ca H', 'Ca II 854.2'):
            continue
        if name == 'H2O':
            continue
        ax.annotate(name.replace('Ca II 849.8', 'Ca II IR'), (lam0, 0.0),
                    xycoords=('data', 'axes fraction'), xytext=(0, 3),
                    textcoords='offset points', ha='center', va='bottom',
                    fontsize=7, color='#666', rotation=90)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(FIGDIR, exist_ok=True)
    irr, rad, chans = load()
    meta = dict(site='VEIT', sequence_id=irr['meta']['sequence_id'],
                system_id=irr['meta']['system_id'],
                cal_date_irr=irr['meta']['instrument_calibration_date_irr'],
                cal_date_rad=rad['meta']['instrument_calibration_date_rad'],
                method='empirical Gaussian line fits (srf.fit_lines); '
                       'to be superseded by the HSRS template fit (task 3b)',
                errors='scan_errors flatten_px=%d' % FLATTEN_PX,
                script='wiggles/phase0a_veit.py')
    lines = fit_all(chans)
    lines.to_parquet(os.path.join(OUT, 'veit_lines.parquet'))
    csv_cols = ['channel', 'err_kind', 'name', 'lam_air', 'lam_vac', 'group',
                'blend', 'use_for_srf', 'ok', 'npix', 'mu', 'mu_err', 'dmu',
                'sigma', 'sigma_err', 'fwhm', 'fwhm_err', 'depth', 'depth_err',
                'c0', 'c1', 'chi2_nu']
    lines[csv_cols].to_csv(os.path.join(REPO, 'wiggles', 'phase0_veit_lines.csv'),
                           index=False, float_format='%.7g')
    mods = models(lines, meta)
    srf.save_srf_models(os.path.join(REPO, 'hypernet', 'data', 'veit_srf_model.json'),
                        mods, meta=meta)
    print('\nSRF models (FWHM = c0 + c1 x + c2 x², x = (λ - 600)/100 nm):')
    for m in mods:
        e = np.sqrt(np.diag(m.cov))
        print('  %-2s c = %s ± %s  chi2_nu %.1f  npts %d  offset %+.3f ± %.3f nm  '
              'FWHM(400/500/600/700/850) = %s' % (
                  m.channel, np.round(m.coeffs, 3), np.round(e, 3), m.chi2_nu,
                  m.npts, m.offset, m.offset_err,
                  ' '.join('%.2f' % m.fwhm(x) for x in (400, 500, 600, 700, 850))))
    fl = lines[lines['err_kind'] == 'flat']
    print('\nPer-line chi2_nu (flattened errors), median by channel:',
          fl.groupby('channel')['chi2_nu'].median().round(2).to_dict())
    rw = lines[lines['err_kind'] == 'raw']
    print('Per-line chi2_nu (raw errors), median by channel:',
          rw.groupby('channel')['chi2_nu'].median().round(2).to_dict())
    print('Failed fits:', fl.loc[~fl['ok'], ['channel', 'name']].values.tolist())
    compare_plan(lines)
    absolute_offsets(lines)
    figure(lines, mods, os.path.join(FIGDIR, 'veit_fwhm_vs_lambda.png'))
    print('\nwrote', OUT, 'veit_lines.parquet; wiggles/phase0_veit_lines.csv; '
          'hypernet/data/veit_srf_model.json; wiggles/figs/phase0/veit_fwhm_vs_lambda.png')


if __name__ == '__main__':
    main()

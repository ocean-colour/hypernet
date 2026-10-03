"""Phase 0, task 8e: do the rho_w wiggles follow the E/L SRF mismatch?

Task 8c found that FWHM_Ld - FWHM_E is ~0.5 nm on HYPSTAR 122302, 122304
and 120242 ("narrow-E") but ~0 on 121222, 121231, 122303 and 122305
("E ~ L").  If H2 drives the wiggles, rho_w at the Fraunhofer lines should
be much rougher on the narrow-E instruments, at the same site and water type
too.

For each of the 224 requested Release 2 L2B files:

- rho_w is the site's agreed product (``whn_explore.product_for``:
  ``reflectance``, or ``reflectance_nosc`` at the turbid sites);
- for each line group of ``srf.LINES`` from Ca H/K to H alpha, three
  quantities within +-5 nm (task 7's windows and continuum,
  ``phase0a_budget_veit``):
  - the relative residual rho_w / continuum - 1;
  - the rho_w'' excess over the quiet level (task 12's metric);
  - the H2 prediction from the sequence's own E (``irradiance``) and Ld line
    fits (``srf.fit_lines``; E depth, centroid and sigma; Ld sigma).
- alpha2, the projection of the rho_w residual onto the H2 profile.

Per sequence, the "blue" summary is the median over Ca H/K, G band, H beta,
Mg b and Na D, where the wiggles are largest and rho_w is well above 0 at
every site.  Per instrument and class: medians and IQRs.  For the same-site
pairs: the ratio of the median blue rms residuals, with a bootstrap 68 %
interval.

Outputs: ``$OS_COLOR/hypernet/wiggles/phase0/release2_wiggles_lines.parquet``;
committed ``hypernet/wiggles/phase0e_{wiggles,summary,pairs,pairs_by_line,lines_by_class}.csv`` and
``figs/phase0/phase0e_wiggles.png``.

Run from the repository root: ``python -m hypernet.wiggles.phase0e_release2_wiggles``.
"""
import os

import numpy as np
import pandas as pd

from hypernet import srf
from hypernet import whn_l1a as wl
from hypernet.whn_explore import product_for
from hypernet.wiggles import FIGDIR, OUT, REPO, WIGGLES_DIR
from hypernet.wiggles import phase0a_budget_veit as b7
from hypernet.wiggles.phase0b_budget import rho2_baseline

REQUEST = os.path.join(REPO, 'docs', 'wiggles_data_request.csv')
NARROW_E = ('HYPSTAR_122302', 'HYPSTAR_122304', 'HYPSTAR_120242')
GROUPS = ['CaHK', 'G', 'Hb', 'Mgb', 'NaD', 'Ha']
BLUE = ['CaHK', 'G', 'Hb', 'Mgb', 'NaD']
PAIRS = [('VEIT: 122304 vs 122305', 'VEIT', 'HYPSTAR_122304', 'HYPSTAR_122305'),
         ('GAIT: 120242 vs 121222', 'GAIT', 'HYPSTAR_120242', 'HYPSTAR_121222'),
         ('MAFR: 121231 vs 122303 (both E ~ L)', 'MAFR', 'HYPSTAR_121231', 'HYPSTAR_122303')]
CLASS_COLORS = {'narrow-E': '#eb6834', 'E ~ L': '#2a78d6'}


def sequence_lines(d, rho, lines=srf.LINES):
    """Per line group: rho_w residual, rho'' excess, H2 prediction, alpha2."""
    w = d['wave']
    E = d['irradiance'][:, 0]
    Ld = d['downwelling_radiance'][:, 0]
    fE = srf.fit_lines(w, E).set_index('name')
    fL = srf.fit_lines(w, Ld).set_index('name')
    d2 = np.gradient(np.gradient(rho, w), w)
    base = rho2_baseline(w, rho)
    rows = []
    for grp in GROUPS:
        g = lines[lines['group'] == grp]
        names = list(g['name'])
        lo, hi = g['lam_air'].min() - b7.HALF, g['lam_air'].max() + b7.HALF
        ok = all(fE.loc[n, 'ok'] and fL.loc[n, 'ok'] for n in names)
        r_rho, m = b7._residual(w, rho, lo, hi)
        r = dict(group=grp, lam_air=float(g['lam_air'].mean()), fit_ok=ok,
                 rho_median=float(np.nanmedian(rho[m])),
                 rms_rho_rel=b7._rms(r_rho, m), rms_rho2=b7._rms(d2, m),
                 rho2_excess=b7._rms(d2, m) / base)
        if ok:
            comps = [(fE.loc[n, 'depth'], fE.loc[n, 'mu'], fE.loc[n, 'sigma'],
                      fL.loc[n, 'sigma']) for n in names]
            r_h2 = b7._residual(w, 1.0 + b7._h2_profile(w, comps), lo, hi)[0]
            a, ea = b7._project(r_rho[m], np.ones(m.sum()), np.column_stack([r_h2])[m])
            k = int(np.argmax([c[0] for c in comps]))
            r.update(rms_H2=b7._rms(r_h2, m), alpha2=a[0], alpha2_err=ea[0],
                     fwhm_E=srf.FWHM_PER_SIGMA * comps[k][2],
                     fwhm_L=srf.FWHM_PER_SIGMA * comps[k][3])
        rows.append(r)
    return pd.DataFrame(rows)


def run():
    req = pd.read_csv(REQUEST, dtype=str)
    out = []
    for i, r in req.iterrows():
        try:
            d = wl.load_l2b(wl.release2_path(r['site'], r['sequence_time'], r['file']))
            var, _ = product_for(r['site'])
            df = sequence_lines(d, d[var][:, 0])
            sysid = d['meta']['system_id']
            out.append(df.assign(row=r['row'], site=r['site'][:-2], system_id=sysid,
                                 sequence_time=r['sequence_time'], product=var,
                                 sza=float(np.nanmean(d['sza'])),
                                 cls='narrow-E' if sysid in NARROW_E else 'E ~ L'))
        except Exception as ex:
            print('%s %s: %s' % (r['site'], r['sequence_time'], ex))
        if (i + 1) % 25 == 0:
            print('%d/%d' % (i + 1, len(req)))
    return pd.concat(out, ignore_index=True)


def per_sequence(lines):
    b = lines[lines['group'].isin(BLUE)]
    keys = ['row', 'site', 'system_id', 'cls', 'sequence_time', 'product', 'sza']
    s = b.groupby(keys).agg(rms_rho_rel=('rms_rho_rel', 'median'),
                            rms_H2=('rms_H2', 'median'),
                            rho2_excess=('rho2_excess', 'median'),
                            alpha2=('alpha2', 'median'),
                            fwhm_L_minus_E=('fwhm_L', 'median'),
                            n_fit=('fit_ok', 'sum')).reset_index()
    fe = b.groupby(keys)['fwhm_E'].median().values
    s['fwhm_L_minus_E'] = s['fwhm_L_minus_E'] - fe
    return s


def summarise(s):
    rows = []
    for by in ('cls', 'system_id'):
        for k, g in s.groupby(by):
            r = dict(group_by=by, group=k, n=len(g),
                     sites=','.join(sorted(g['site'].unique())))
            for q in ('rms_rho_rel', 'rms_H2', 'rho2_excess', 'alpha2', 'fwhm_L_minus_E'):
                r[q + '_median'] = g[q].median()
                r[q + '_q25'] = g[q].quantile(0.25)
                r[q + '_q75'] = g[q].quantile(0.75)
            rows.append(r)
    return pd.DataFrame(rows)


def pairs(s, nboot=2000, seed=1):
    rng = np.random.default_rng(seed)
    rows = []
    for name, site, a, b in PAIRS:
        for q in ('rms_rho_rel', 'rho2_excess'):
            va = s[(s['site'] == site) & (s['system_id'] == a)][q].dropna().values
            vb = s[(s['site'] == site) & (s['system_id'] == b)][q].dropna().values
            if len(va) < 3 or len(vb) < 3:
                continue
            ratio = np.median(va) / np.median(vb)
            boot = [np.median(rng.choice(va, len(va))) / np.median(rng.choice(vb, len(vb)))
                    for _ in range(nboot)]
            lo, hi = np.percentile(boot, [16, 84])
            rows.append(dict(pair=name, quantity=q, n_a=len(va), n_b=len(vb),
                             median_a=np.median(va), median_b=np.median(vb), ratio=ratio,
                             ratio_lo68=lo, ratio_hi68=hi))
    return pd.DataFrame(rows)


def pair_lines(lines):
    """Per line group, for each same-site pair: median rms_rho_rel of a over b,
    and the H2 predictions of both."""
    rows = []
    for name, site, a, b in PAIRS:
        for grp in GROUPS:
            ga = lines[(lines['site'] == site) & (lines['system_id'] == a) &
                       (lines['group'] == grp)]
            gb = lines[(lines['site'] == site) & (lines['system_id'] == b) &
                       (lines['group'] == grp)]
            ma, mb = ga['rms_rho_rel'].median(), gb['rms_rho_rel'].median()
            rows.append(dict(pair=name, group=grp, rms_a=ma, rms_b=mb, ratio=ma / mb,
                             H2_a=ga['rms_H2'].median(), H2_b=gb['rms_H2'].median(),
                             # excess over the E ~ L floor, in quadrature
                             excess_a=np.sqrt(max(ma ** 2 - mb ** 2, 0.0))))
    return pd.DataFrame(rows)


def class_lines(lines):
    """Per line group and class: median measured rms, H2 prediction, and the
    narrow-E excess over the E ~ L floor (in quadrature)."""
    t = lines.groupby(['group', 'cls']).agg(rms=('rms_rho_rel', 'median'),
                                             H2=('rms_H2', 'median')).unstack()
    t = t.loc[GROUPS]
    out = pd.DataFrame({'rms_E~L': t[('rms', 'E ~ L')], 'rms_narrowE': t[('rms', 'narrow-E')],
                        'ratio': t[('rms', 'narrow-E')] / t[('rms', 'E ~ L')],
                        'H2_E~L': t[('H2', 'E ~ L')], 'H2_narrowE': t[('H2', 'narrow-E')]})
    out['excess_narrowE'] = np.sqrt(np.clip(out['rms_narrowE'] ** 2 - out['rms_E~L'] ** 2,
                                            0, None))
    return out


def figure(s, lines, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(12, 4.9), gridspec_kw=dict(wspace=0.28))
    order = (s.groupby('system_id')['fwhm_L_minus_E'].median().sort_values(ascending=False)
             .index.tolist())
    for i, inst in enumerate(order):
        g = s[s['system_id'] == inst]
        c = CLASS_COLORS[g['cls'].iloc[0]]
        x = i + np.random.default_rng(i).uniform(-0.15, 0.15, len(g))
        ax.semilogy(x, g['rms_rho_rel'], ls='', marker='o', ms=5, alpha=0.75, color=c)
        ax.plot([i - 0.3, i + 0.3], [g['rms_rho_rel'].median()] * 2, color='#333333', lw=2)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(['%s\n%s' % (k.replace('HYPSTAR_', ''),
                                    ','.join(sorted(s[s['system_id'] == k]['site'].unique())))
                        for k in order], fontsize=8.5)
    ax.set_ylabel('ρw relative residual at the lines\n(rms ±5 nm, median Ca H/K–Na D)')
    ax.set_title('per sequence, instruments ordered by FWHM_L − FWHM_E', fontsize=10,
                 loc='left')
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], ls='', marker='o', color=c, label=k)
                       for k, c in CLASS_COLORS.items()] +
              [Line2D([], [], color='#333333', lw=2, label='median')],
              fontsize=8.5, frameon=False, loc='upper right')
    for cls, c in CLASS_COLORS.items():
        g = s[s['cls'] == cls]
        bx.loglog(g['rms_H2'], g['rms_rho_rel'], ls='', marker='o', ms=5, alpha=0.7, color=c,
                  label=cls)
    lim = [np.nanmin(s[['rms_H2', 'rms_rho_rel']].values) / 1.5,
           np.nanmax(s[['rms_H2', 'rms_rho_rel']].values) * 1.5]
    bx.plot(lim, lim, color='#888', lw=0.8)
    bx.set_xlim(lim)
    bx.set_ylim(lim)
    bx.set_xlabel('H2 prediction (rms, from the sequence\'s own E and Ld fits)')
    bx.set_ylabel('measured ρw relative residual (rms)')
    bx.set_title('per sequence (median over Ca H/K–Na D); line = 1:1', fontsize=10, loc='left')
    bx.legend(fontsize=8.5, frameon=False, loc='upper left')
    for a in (ax, bx):
        a.grid(alpha=0.25, lw=0.6)
        a.spines[['top', 'right']].set_visible(False)
        a.tick_params(labelsize=9)
    fig.suptitle('Release 2: ρw wiggles at the Fraunhofer lines by instrument (%d sequences)'
                 % len(s), fontsize=10.5, x=0.06, ha='left')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)


def main():
    pd.set_option('display.width', 220)
    lines = run()
    os.makedirs(OUT, exist_ok=True)
    lines.to_parquet(os.path.join(OUT, 'release2_wiggles_lines.parquet'))
    s = per_sequence(lines)
    summ = summarise(s)
    pr = pairs(s)
    s.to_csv(os.path.join(WIGGLES_DIR, 'phase0e_wiggles.csv'), index=False, float_format='%.5g')
    summ.to_csv(os.path.join(WIGGLES_DIR, 'phase0e_summary.csv'), index=False,
                float_format='%.4g')
    pr.to_csv(os.path.join(WIGGLES_DIR, 'phase0e_pairs.csv'), index=False, float_format='%.4g')
    os.makedirs(FIGDIR, exist_ok=True)
    figure(s, lines, os.path.join(FIGDIR, 'phase0e_wiggles.png'))
    cols = ['group', 'n', 'sites', 'fwhm_L_minus_E_median', 'rms_rho_rel_median',
            'rms_rho_rel_q25', 'rms_rho_rel_q75', 'rms_H2_median', 'rho2_excess_median',
            'alpha2_median']
    print(summ[cols].to_string(index=False, float_format='%.4f'))
    print('\nsame-site pairs (ratio of medians, a / b, bootstrap 68 %):')
    print(pr.to_string(index=False, float_format='%.4f'))
    cl = class_lines(lines)
    pl = pair_lines(lines)
    cl.to_csv(os.path.join(WIGGLES_DIR, 'phase0e_lines_by_class.csv'), float_format='%.5g')
    pl.to_csv(os.path.join(WIGGLES_DIR, 'phase0e_pairs_by_line.csv'), index=False,
              float_format='%.5g')
    print('\nper line group, by class (median rms; excess = narrow-E over the E ~ L floor, '
          'in quadrature, vs the H2 prediction):')
    print(cl.to_string(float_format='%.4f'))
    print('\nsame-site pairs per line group:')
    print(pl.to_string(index=False, float_format='%.4f'))
    print('\nwrote phase0e_*.csv, figs/phase0/phase0e_wiggles.png')


if __name__ == '__main__':
    main()

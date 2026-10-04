"""Phase 0b, task 11: stability of the SRFs per instrument (gate G0(b)).

Reads ``srf_models.parquet`` (task 10) and, per instrument x calibration
period x channel (template method by default; Ld is the reference L channel,
Q&A Q11):

- FWHM at 400, 500, 600 and 700 nm, and the centroid offset L - E (the
  difference of the per-sequence model offsets): mean, std and range, and
  whether the std is below the 0.2 nm criterion of G0(b);
- the dependence of FWHM(600) and of FWHM_Ld - FWHM_E at 450 nm on SZA, sky
  condition, season and time (linear slopes with errors when n >= 4, means
  per sky class);
- the pair comparisons of the plan: 122304 and 121222 before/after
  recalibration, 122302 at BEFR vs THFR (site vs instrument), and the
  two-instruments-one-site pairs at VEIT, GAIT and MAFR (difference of means
  with its error and z).

Outputs: ``hypernet/wiggles/phase0_stability.csv``, ``hypernet/wiggles/phase0_stability_pairs.csv``,
``hypernet/wiggles/phase0_stability_trends.csv`` and
``hypernet/wiggles/figs/phase0/stability_vs_sza.png``, ``stability_vs_time.png``.

Run from the repository root::

    python -m hypernet.wiggles.phase0b_stability [--method template|empirical]
"""
import argparse
import os

import numpy as np
import pandas as pd
from hypernet.wiggles import DATA_DIR, FIGDIR, OUT, REPO, WIGGLES_DIR  # noqa: F401

LAMS = (400, 500, 600, 700)
CRITERION = 0.2   # nm, G0(b)
#: Categorical palette slots in fixed order (dataviz reference palette).
SLOTS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7',
         '#e34948']
#: One marker per sky class (whn_l1a.SKY_CLASSES; task 8e).
SKY_MARKERS = {'clear': 'o', 'cloudy': 's', 'overcast': '^', 'broken': 'D'}

#: The plan's pair comparisons.  Each side selects sequences by any of
#: system_id, site_code, cal_period; 'first'/'last' pick the earliest/latest
#: calibration period present for that instrument.
PAIRS = [
    ('122304 recalibration', dict(system_id='HYPSTAR_122304', cal_period='first'),
     dict(system_id='HYPSTAR_122304', cal_period='last')),
    ('121222 recalibration', dict(system_id='HYPSTAR_121222', cal_period='first'),
     dict(system_id='HYPSTAR_121222', cal_period='last')),
    ('122302 BEFR vs THFR (site)', dict(system_id='HYPSTAR_122302', site_code='BEFR'),
     dict(system_id='HYPSTAR_122302', site_code='THFR')),
    # 122305 ran between the two 122304 calibration periods: compare with each
    ('VEIT 122304 pre-recal vs 122305',
     dict(system_id='HYPSTAR_122304', site_code='VEIT', cal_period='first'),
     dict(system_id='HYPSTAR_122305', site_code='VEIT')),
    ('VEIT 122304 post-recal vs 122305',
     dict(system_id='HYPSTAR_122304', site_code='VEIT', cal_period='last'),
     dict(system_id='HYPSTAR_122305', site_code='VEIT')),
    ('GAIT 121222 vs 120242', dict(system_id='HYPSTAR_121222', site_code='GAIT'),
     dict(system_id='HYPSTAR_120242', site_code='GAIT')),
    ('MAFR 121231 vs 122303', dict(system_id='HYPSTAR_121231', site_code='MAFR'),
     dict(system_id='HYPSTAR_122303', site_code='MAFR')),
]


def per_sequence(models, method):
    """One row per sequence: FWHM per channel at LAMS, Ld - E, offsets."""
    m = models[models['method'] == method]
    keys = ['site_code', 'sequence_time', 'system_id', 'cal_period', 'sza_l1a', 'sky',
            'month', 'water_type']
    keys += [k for k in ('ld_ed_750', 'ed_cv_750') if k in m]
    rows = []
    for k, g in m.groupby(['site_code', 'sequence_time', 'system_id', 'cal_period'],
                          dropna=False):
        g = g.set_index('channel')
        r = {c: g[c].iloc[0] for c in keys}
        for ch in g.index:
            for lam in LAMS + (450,):
                r['fwhm_%d_%s' % (lam, ch)] = g.loc[ch, 'fwhm_%d' % lam] \
                    if 'fwhm_%d' % lam in g else np.nan
                r['fwhm_%d_err_%s' % (lam, ch)] = g.loc[ch, 'fwhm_%d_err' % lam] \
                    if 'fwhm_%d_err' % lam in g else np.nan
            r['offset_' + ch] = g.loc[ch, 'offset']
            r['offset_err_' + ch] = g.loc[ch, 'offset_err']
        for ch in ('Ld', 'Lu'):
            if ch in g.index and 'E' in g.index:
                r['dfwhm_450_%s_E' % ch] = r['fwhm_450_' + ch] - r['fwhm_450_E']
                r['doffset_%s_E' % ch] = r['offset_' + ch] - r['offset_E']
        rows.append(r)
    s = pd.DataFrame(rows)
    t = pd.to_datetime(s['sequence_time'], format='%Y%m%dT%H%M')
    s['date'] = t
    s['year'] = t.dt.year + (t.dt.dayofyear - 1) / 365.25
    s['sky'] = s['sky'].fillna('untagged')
    return s


def group_stats(s):
    """Per instrument x cal period x channel: mean, std, range, pass/fail."""
    rows = []
    for (inst, cal), g in s.groupby(['system_id', 'cal_period'], dropna=False):
        for ch in ('E', 'Ld', 'Lu'):
            r = dict(system_id=inst, cal_period=cal, channel=ch, n=len(g),
                     sites=','.join(sorted(g['site_code'].unique())))
            stds = []
            for lam in LAMS:
                v = g.get('fwhm_%d_%s' % (lam, ch), pd.Series(dtype=float)).dropna()
                r['fwhm_%d_mean' % lam] = v.mean()
                r['fwhm_%d_std' % lam] = v.std(ddof=1) if len(v) > 1 else np.nan
                r['fwhm_%d_range' % lam] = v.max() - v.min() if len(v) else np.nan
                stds.append(r['fwhm_%d_std' % lam])
            r['max_std'] = np.nanmax(stds) if np.any(np.isfinite(stds)) else np.nan
            r['passes_0p2'] = bool(r['max_std'] < CRITERION) if np.isfinite(r['max_std']) \
                else None
            if ch != 'E':
                d = g.get('doffset_%s_E' % ch, pd.Series(dtype=float)).dropna()
                r['offset_L_E_mean'] = d.mean()
                r['offset_L_E_std'] = d.std(ddof=1) if len(d) > 1 else np.nan
                dd = g.get('dfwhm_450_%s_E' % ch, pd.Series(dtype=float)).dropna()
                r['dfwhm_450_L_E_mean'] = dd.mean()
                r['dfwhm_450_L_E_std'] = dd.std(ddof=1) if len(dd) > 1 else np.nan
            rows.append(r)
    return pd.DataFrame(rows)


def _slope(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4 or np.ptp(x[m]) == 0:
        return np.nan, np.nan
    A = np.column_stack([np.ones(m.sum()), x[m]])
    c, res, *_ = np.linalg.lstsq(A, y[m], rcond=None)
    dof = m.sum() - 2
    s2 = np.sum((y[m] - A @ c) ** 2) / dof
    cov = s2 * np.linalg.inv(A.T @ A)
    return c[1], np.sqrt(cov[1, 1])


def trends(s):
    """Slopes of FWHM(600) (E, Ld) and Ld - E at 450 nm vs SZA, season, time."""
    rows = []
    season = np.cos(2 * np.pi * (s['month'].astype(float) - 6.5) / 12)  # +1 in June/July
    for (inst, cal), g in s.groupby(['system_id', 'cal_period'], dropna=False):
        gi = g.index
        for q in ('fwhm_600_E', 'fwhm_600_Ld', 'dfwhm_450_Ld_E', 'doffset_Ld_E'):
            if q not in g:
                continue
            y = g[q].values.astype(float)
            r = dict(system_id=inst, cal_period=cal, quantity=q, n=int(np.isfinite(y).sum()))
            xs = [('sza', g['sza_l1a'].values.astype(float)),
                  ('season', season.loc[gi].values), ('year', g['year'].values)]
            if 'ld_ed_750' in g:                       # task 8e: the sky index
                xs.append(('sky', g['ld_ed_750'].values.astype(float)))
            for name, x in xs:
                r['slope_' + name], r['slope_%s_err' % name] = _slope(x, y)
            for sky, h in g.groupby('sky'):
                r['mean_sky_' + sky] = h[q].mean()
            rows.append(r)
    return pd.DataFrame(rows)


def _select(s, sel):
    g = s.copy()
    for k in ('system_id', 'site_code'):
        if k in sel:
            g = g[g[k] == sel[k]]
    if 'cal_period' in sel and len(g):
        cals = sorted(g['cal_period'].dropna().unique(),
                      key=lambda c: c.split('/')[-1])
        if len(cals) < 2 and sel['cal_period'] in ('first', 'last'):
            return g.iloc[0:0]
        pick = cals[0] if sel['cal_period'] == 'first' else (
            cals[-1] if sel['cal_period'] == 'last' else sel['cal_period'])
        g = g[g['cal_period'] == pick]
    return g


def pairs(s):
    rows = []
    for name, a, b in PAIRS:
        ga, gb = _select(s, a), _select(s, b)
        for q in ['fwhm_%d_%s' % (lam, ch) for ch in ('E', 'Ld') for lam in LAMS] + \
                ['dfwhm_450_Ld_E', 'doffset_Ld_E']:
            va = ga.get(q, pd.Series(dtype=float)).dropna()
            vb = gb.get(q, pd.Series(dtype=float)).dropna()
            r = dict(pair=name, quantity=q, n_a=len(va), n_b=len(vb),
                     cal_a=','.join(ga['cal_period'].dropna().unique()),
                     cal_b=','.join(gb['cal_period'].dropna().unique()))
            if len(va) >= 2 and len(vb) >= 2:
                d = vb.mean() - va.mean()
                e = np.hypot(va.std(ddof=1) / np.sqrt(len(va)), vb.std(ddof=1) / np.sqrt(len(vb)))
                r.update(mean_a=va.mean(), mean_b=vb.mean(), diff=d, err=e, z=d / e)
            rows.append(r)
    return pd.DataFrame(rows)


def figures(s, prefix='stability', label='template fits'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    insts = sorted(s['system_id'].dropna().unique())
    col = {k: SLOTS[i % len(SLOTS)] for i, k in enumerate(insts)}
    panels = [('fwhm_600_E', 'FWHM_E at 600 nm (nm)'), ('fwhm_600_Ld', 'FWHM_Ld at 600 nm (nm)'),
              ('dfwhm_450_Ld_E', 'FWHM_Ld − FWHM_E at 450 nm (nm)'),
              ('doffset_Ld_E', 'centroid offset Ld − E (nm)')]
    for xname, xlab, fname in (('sza_l1a', 'SZA (deg)', prefix + '_vs_sza.png'),
                               ('date', 'date', prefix + '_vs_time.png')):
        fig, axs = plt.subplots(2, 2, figsize=(10, 7.5), sharex=True)
        for ax, (q, lab) in zip(axs.ravel(), panels):
            for inst in insts:
                g = s[s['system_id'] == inst]
                for sky, h in g.groupby('sky'):
                    ax.plot(h[xname], h[q], ls='', marker=SKY_MARKERS.get(sky, 'D'), ms=7,
                            color=col[inst], mfc=col[inst] if sky != 'untagged' else 'white',
                            label=None)
            ax.set_ylabel(lab, fontsize=9)
            ax.grid(alpha=0.25, lw=0.6)
            ax.spines[['top', 'right']].set_visible(False)
            ax.tick_params(labelsize=8.5)
        for ax in axs[1]:
            ax.set_xlabel(xlab)
        from matplotlib.lines import Line2D
        h = [Line2D([], [], ls='', marker='o', ms=7, color=col[i], label=i.replace('HYPSTAR_', ''))
             for i in insts]
        h += [Line2D([], [], ls='', marker=m, ms=7, color='#555', label=k)
              for k, m in SKY_MARKERS.items()]
        h += [Line2D([], [], ls='', marker='o', ms=7, mfc='white', mec='#555', label='sky untagged')]
        fig.suptitle('SRF stability per instrument (%s; n = %d sequences)' % (label, len(s)),
                     fontsize=10.5, y=0.995)
        fig.legend(handles=h, loc='upper center', bbox_to_anchor=(0.5, 0.965),
                   ncol=min(len(h), 7), frameon=False, fontsize=8.5)
        fig.tight_layout(rect=(0, 0, 1, 0.89))
        fig.savefig(os.path.join(FIGDIR, fname), dpi=150)
        plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description='Phase 0b SRF stability')
    ap.add_argument('--method', default='template', choices=['template', 'empirical'])
    a = ap.parse_args(argv)
    models = pd.read_parquet(os.path.join(OUT, 'srf_models.parquet'))
    s = per_sequence(models, a.method)
    st, tr, pr = group_stats(s), trends(s), pairs(s)
    st.to_csv(os.path.join(WIGGLES_DIR, 'phase0_stability.csv'), index=False,
              float_format='%.4f')
    pr.to_csv(os.path.join(WIGGLES_DIR, 'phase0_stability_pairs.csv'), index=False,
              float_format='%.4f')
    tr.to_csv(os.path.join(WIGGLES_DIR, 'phase0_stability_trends.csv'), index=False,
              float_format='%.4f')
    os.makedirs(FIGDIR, exist_ok=True)
    figures(s)
    pd.set_option('display.width', 220)
    print('%d sequences, %d instrument x cal-period groups (%s method)'
          % (len(s), st[['system_id', 'cal_period']].drop_duplicates().shape[0], a.method))
    cols = ['system_id', 'cal_period', 'channel', 'n'] + ['fwhm_%d_mean' % l for l in LAMS] + \
        ['max_std', 'passes_0p2', 'offset_L_E_mean', 'dfwhm_450_L_E_mean']
    print(st[[c for c in cols if c in st]].to_string(index=False, float_format='%.3f'))
    done = pr.dropna(subset=['z']) if 'z' in pr else pr.iloc[0:0]
    print('\npair comparisons with n >= 2 on both sides: %d of %d'
          % (len(done), len(pr)))
    if len(done):
        print(done[['pair', 'quantity', 'n_a', 'n_b', 'diff', 'err', 'z']]
              .to_string(index=False, float_format='%.3f'))
    print('wrote hypernet/wiggles/phase0_stability*.csv and hypernet/wiggles/figs/phase0/stability_vs_*.png')


if __name__ == '__main__':
    main()

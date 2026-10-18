"""Phase 0b, task 12: the H1/H2 wiggle budget of every delivered sequence.

Generalises ``hypernet/wiggles/phase0a_budget_veit.py`` (task 7): for every sequence in
``request_index.parquet`` with L1A_IRR, L1A_RAD and L2A_REF, runs its
:func:`budget` with the E and Ld empirical line fits of ``line_fits.parquet``
(task 10), and adds the rho_w'' excess at each line -- rms(rho_w'') within
+-5 nm over a baseline, the 25th percentile of rms(rho_w'') in 10 nm tiles
across 400-700 nm.  Most visible tiles contain Fraunhofer lines, so the
median tile is not line-free; the quietest quarter is closer to it.

Summaries, per sequence (weighted over the SRF lines) and then by instrument
and water type: alpha2 (fraction of the H2 profile present) in Ld/Ed, Lu/Ed
and rho_w; alpha1; rms_H2 / rms_H1; rho_w measured vs rho_w x H2.

Outputs: ``$OS_COLOR/hypernet/wiggles/phase0/budget_lines.parquet`` (every
sequence x line), ``hypernet/wiggles/phase0_budget.csv`` (per sequence),
``hypernet/wiggles/phase0_budget_summary.csv`` (by instrument and water type) and
``hypernet/wiggles/figs/phase0/budget_all.png``.

Run from the repository root: ``python -m hypernet.wiggles.phase0b_budget``.
"""
import os

import numpy as np
import pandas as pd

from hypernet import whn_l1a as wl  # noqa: E402
from hypernet.wiggles import phase0a_budget_veit as b7  # noqa: E402
from hypernet.wiggles import DATA_DIR, FIGDIR, OUT, REPO, WIGGLES_DIR  # noqa: F401

WATER_COLORS = {'clear': '#2a78d6', 'dark': '#eb6834', 'turbid': '#1baf7a'}
KEYS = ['site_code', 'sequence_time', 'system_id', 'cal_period', 'water_type', 'sza_l1a',
        'sky']


def rho2_baseline(wave, rho, lo=400.0, hi=700.0, width=10.0):
    """25th percentile over 10 nm tiles of rms(rho_w'') -- the quiet level."""
    d2 = np.gradient(np.gradient(rho, wave), wave)
    vals = []
    for a in np.arange(lo, hi, width):
        m = (wave >= a) & (wave < a + width) & np.isfinite(d2)
        if m.sum() > 3:
            vals.append(np.sqrt(np.mean(d2[m] ** 2)))
    return float(np.percentile(vals, 25)) if vals else np.nan


def _wmean(v, e):
    m = np.isfinite(v) & np.isfinite(e) & (e > 0)
    if not m.any():
        return np.nan, np.nan
    w = 1 / e[m] ** 2
    return float(np.sum(w * v[m]) / w.sum()), float(1 / np.sqrt(w.sum()))


def run(index, line_fits):
    ok = index['status'].isin(['delivered', 'extra']) & index['L1A_IRR'].notna() & \
        index['L1A_RAD'].notna() & index['L2A_REF'].notna()
    seqs = index[ok].reset_index(drop=True)
    per_line, per_seq = [], []
    for i, r in seqs.iterrows():
        keys = {k: r.get(k) for k in KEYS}
        f = line_fits[(line_fits['site_code'] == r['site_code']) &
                      (line_fits['sequence_time'] == r['sequence_time'])]
        fE = f[f['channel'] == 'E'].reset_index(drop=True)
        fL = f[f['channel'] == 'Ld'].reset_index(drop=True)
        if not len(fE) or not len(fL):
            print('[%d/%d] %s %s: no line fits, skipped' % (i + 1, len(seqs), r['site_code'],
                                                           r['sequence_time']))
            continue
        try:
            irr = wl.load_l1a_irr(r['L1A_IRR'])
            rad = wl.load_l1a_rad(r['L1A_RAD'])
            l2a = wl.load_l2a(r['L2A_REF'])
            df = b7.budget(irr, rad, l2a, fE, fL)
        except Exception as ex:
            print('[%d/%d] %s %s: %s' % (i + 1, len(seqs), r['site_code'], r['sequence_time'], ex))
            continue
        base = rho2_baseline(l2a['wave'], l2a['reflectance'][:, 0])
        df['rho2_baseline'] = base
        df['rho2_excess'] = df['rms_rho2'] / base
        df = df.assign(**keys)
        per_line.append(df)
        c = df[df['use_for_srf'] & df['fit_ok']]
        s = dict(keys, n_lines=len(c), rho2_baseline=base,
                 rho2_excess_median=float(c['rho2_excess'].median()),
                 h2_over_h1_median=float((c['rms_H2'] / c['rms_H1']).median()),
                 rho_meas_over_pred_H2_median=float((c['rms_rho_abs'] / c['rms_rho_pred_H2'])
                                                    .median()))
        for pre in ('LdEd', 'LuEd', 'rho'):
            for a in ('alpha2', 'alpha1'):
                s['%s_%s' % (a, pre)], s['%s_%s_err' % (a, pre)] = _wmean(
                    c['%s_%s' % (a, pre)].values, c['%s_err_%s' % (a, pre)].values)
        per_seq.append(s)
        print('[%d/%d] %s %s %s: alpha2 Ld/Ed %.2f, rho %.2f' % (
            i + 1, len(seqs), r['site_code'], r['sequence_time'], r.get('system_id'),
            s['alpha2_LdEd'], s['alpha2_rho']))
    lines = pd.concat(per_line, ignore_index=True) if per_line else pd.DataFrame()
    return lines, pd.DataFrame(per_seq)


def summarise(per_seq):
    cols = ['alpha2_LdEd', 'alpha2_LuEd', 'alpha2_rho', 'alpha1_LdEd', 'h2_over_h1_median',
            'rho_meas_over_pred_H2_median', 'rho2_excess_median']
    rows = []
    for by in ('system_id', 'water_type'):
        for k, g in per_seq.groupby(by):
            r = dict(group_by=by, group=k, n=len(g))
            for c in cols:
                r[c + '_median'] = g[c].median()
                r[c + '_iqr'] = g[c].quantile(0.75) - g[c].quantile(0.25) if len(g) > 1 else np.nan
            rows.append(r)
    return pd.DataFrame(rows)


def figure(lines, per_seq, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw=dict(wspace=0.3))
    insts = sorted(per_seq['system_id'].dropna().unique())
    for i, inst in enumerate(insts):
        g = per_seq[per_seq['system_id'] == inst]
        for j, (pre, mk) in enumerate((('LdEd', 'o'), ('rho', 'D'))):
            x = i + (j - 0.5) * 0.25 + np.linspace(-0.05, 0.05, len(g))
            ax.errorbar(x, g['alpha2_' + pre], yerr=g['alpha2_%s_err' % pre], fmt=mk, ms=6,
                        lw=1, color=['#333333', '#1baf7a'][j],
                        label=['Ld/Ed', 'ρw'][j] if i == 0 else None)
    ax.axhline(1, color='#eb6834', lw=0.8, ls='--')
    ax.axhline(0, color='#888', lw=0.8)
    ax.set_xticks(range(len(insts)))
    ax.set_xticklabels([s.replace('HYPSTAR_', '') for s in insts], fontsize=9)
    ax.set_ylabel('α₂: fraction of the H2 profile present\n(weighted over the SRF lines)')
    ax.set_title('per sequence, by instrument', fontsize=10, loc='left')
    ax.legend(fontsize=8.5, frameon=False)
    c = lines[lines['use_for_srf'] & lines['fit_ok']]
    for wt, g in c.groupby('water_type'):
        bx.loglog(g['rms_rho_pred_H2'], g['rms_rho_abs'], ls='', marker='o', ms=5, alpha=0.8,
                  color=WATER_COLORS.get(wt, '#888'), label=wt)
    lim = [np.nanmin(c[['rms_rho_pred_H2', 'rms_rho_abs']].values) / 2,
           np.nanmax(c[['rms_rho_pred_H2', 'rms_rho_abs']].values) * 2]
    bx.plot(lim, lim, color='#888', lw=0.8)
    bx.set_xlabel('predicted ρw residual, ρw × H2 (rms, ±5 nm)')
    bx.set_ylabel('measured ρw residual (rms, ±5 nm)')
    bx.set_title('per line (SRF lines), by water type', fontsize=10, loc='left')
    bx.legend(fontsize=8.5, frameon=False)
    for a in (ax, bx):
        a.grid(alpha=0.25, lw=0.6)
        a.spines[['top', 'right']].set_visible(False)
        a.tick_params(labelsize=9)
    fig.suptitle('H1/H2 wiggle budget, %d sequences' % len(per_seq), fontsize=10.5, x=0.06,
                 ha='left')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)


def main():
    pd.set_option('display.width', 220)
    index = pd.read_parquet(os.path.join(OUT, 'request_index.parquet'))
    line_fits = pd.read_parquet(os.path.join(OUT, 'line_fits.parquet'))
    lines, per_seq = run(index, line_fits)
    if not len(per_seq):
        print('no sequences with L1A + L2A and line fits')
        return
    lines.to_parquet(os.path.join(OUT, 'budget_lines.parquet'))
    per_seq.to_csv(os.path.join(WIGGLES_DIR, 'phase0_budget.csv'), index=False,
                   float_format='%.5g')
    summ = summarise(per_seq)
    summ.to_csv(os.path.join(WIGGLES_DIR, 'phase0_budget_summary.csv'), index=False,
                float_format='%.4g')
    os.makedirs(FIGDIR, exist_ok=True)
    figure(lines, per_seq, os.path.join(FIGDIR, 'budget_all.png'))
    print(summ[['group_by', 'group', 'n', 'alpha2_LdEd_median', 'alpha2_rho_median',
                'h2_over_h1_median_median', 'rho2_excess_median_median']]
          .to_string(index=False, float_format='%.2f'))
    print('wrote budget_lines.parquet, hypernet/wiggles/phase0_budget*.csv, '
          'hypernet/wiggles/figs/phase0/budget_all.png')


if __name__ == '__main__':
    main()

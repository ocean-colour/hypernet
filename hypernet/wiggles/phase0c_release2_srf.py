"""Phase 0, task 8c: SRF stability from the Release 2 L2B archive, before the
L1A delivery.

Release 2 L2B files hold the sequence-mean sky radiance Ld and the
irradiance on the native L grid (``whn_l1a.load_l2b``).

- **L channel:** Ld is never resampled, so the HSRS template fit
  (``srf.fit_template_windows``, veil 0) gives the L SRF directly.
- **E channel:** ``irradiance`` is E resampled onto the L grid with
  ``np.interp``.  The fit models that resampling (``via_grid`` = the E grid).
  Only 122304's E grid is known (from the VEIT L1A, IRR cal 2024-11-12); it
  is used as a proxy for every other instrument and calibration, and those E
  results are flagged ``e_grid = 'proxy'``.

**Errors:** L2B carries no per-scan spectra and no per-pixel error for Ld or
the irradiance (``std_*`` exist only for the reflectances and Lw).  The fits
are therefore unweighted, and ``curve_fit`` scales the covariance by the
residuals, so a parameter's error reflects the actual misfit in its window.
``--validate`` checks this route against the L1A fits of the VEIT sample.

Steps:

1. ``--validate``: VEIT sample, L2A (the same products as L2B) vs L1A
   template fits, with and without the E-grid resampling model.
2. Every row of ``docs/wiggles_data_request.csv`` (224 L2B files): template
   fits of Ld and E, SRF models per sequence.
3. The task 11 stability analysis (``phase0b_stability``) on those models.

Outputs: ``$OS_COLOR/hypernet/wiggles/phase0/release2_{template,models}.parquet``;
committed ``hypernet/wiggles/phase0c_validation.csv``,
``phase0c_stability{,_pairs,_trends}.csv`` and
``figs/phase0/phase0c_stability_vs_{sza,time}.png``,
``phase0c_fwhm_by_instrument.png``.

Run from the repository root::

    python -m hypernet.wiggles.phase0c_release2_srf [--validate] [--limit N]
"""
import argparse
import os
import time

import numpy as np
import pandas as pd

from hypernet import srf, whn_srf
from hypernet import whn_l1a as wl
from hypernet.wiggles import FIGDIR, OUT, REPO, WIGGLES_DIR
from hypernet.wiggles import phase0b_stability as stab

REQUEST = os.path.join(REPO, 'docs', 'wiggles_data_request.csv')
VEIT = dict(site='VEIT', seq_time='20260604T0845')
LAMS = (400, 450, 500, 600, 700, 850)
WATER_TYPE = {'VEIT': 'clear', 'GAIT': 'clear', 'BEFR': 'dark', 'THFR': 'dark',
              'WRUK': 'dark', 'MAFR': 'turbid', 'LPAR': 'turbid', 'O1BE': 'turbid'}
#: Palette slots in fixed order (dataviz reference palette).
SLOTS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7',
         '#e34948']


def e_grid():
    """The E grid of HYPSTAR 122304 (IRR cal 2024-11-12), from the VEIT L1A."""
    f = wl.sequence_files(**VEIT)
    return wl.load_l1a_irr(f['L1A_IRR'])['wave']


def fit_l2(wave, Ld, E, ref, gE, meta):
    """Template fits of Ld and E (E through the E-grid resampling)."""
    rows, mods = [], []
    for ch, spec, via in (('Ld', Ld, None), ('E', E, gE)):
        t = srf.fit_template_windows(wave, spec, ref[0], ref[1], veil=False, via_grid=via)
        t.insert(0, 'channel', ch)
        rows.append(t)
        mods.append(srf.SRFModel.from_lines(ch, t, scale_cov=True,
                                            instrument=meta.get('system_id', ''),
                                            meta=dict(meta, channel=ch, method='template')))
    return pd.concat(rows, ignore_index=True), mods


def validate(ref, gE):
    """VEIT sample: the L2A route (as L2B) vs the L1A template fits."""
    f = wl.sequence_files(**VEIT)
    l1 = whn_srf.fit_files(f, ref=ref, veil=False, empirical=False)
    l1m = {m.channel: m for m in l1['models']}
    l2 = wl.load_l2a(f['L2A_REF'])
    w = l2['wave']
    Ld, E = l2['downwelling_radiance'][:, 0], l2['irradiance'][:, 0]
    rows = []
    for label, via in (('L2 route, E via E grid', gE), ('L2 route, E fitted on L grid', None)):
        mods = {}
        for ch, spec, v in (('Ld', Ld, None), ('E', E, via)):
            t = srf.fit_template_windows(w, spec, ref[0], ref[1], veil=False, via_grid=v)
            mods[ch] = srf.SRFModel.from_lines(ch, t, scale_cov=True)
            mods[ch].meta['n_ok'] = int(t['ok'].sum())
        for ch in ('E', 'Ld'):
            for lam in LAMS:
                rows.append(dict(route=label, channel=ch, lam=lam,
                                 fwhm_l2=mods[ch].fwhm(lam), err_l2=mods[ch].fwhm_err(lam),
                                 fwhm_l1a=l1m[ch].fwhm(lam), err_l1a=l1m[ch].fwhm_err(lam),
                                 offset_l2=mods[ch].offset, offset_l1a=l1m[ch].offset))
    v = pd.DataFrame(rows)
    v['diff'] = v['fwhm_l2'] - v['fwhm_l1a']
    v.to_csv(os.path.join(WIGGLES_DIR, 'phase0c_validation.csv'), index=False, float_format='%.4f')
    pd.set_option('display.width', 200)
    print('[validate] VEIT sample: template FWHM via the L2 products vs via L1A (nm)')
    print(v.pivot_table(index=['channel', 'lam'], columns='route', values='diff')
          .round(3).to_string())
    print(v[v['route'] == 'L2 route, E via E grid'][['channel', 'lam', 'fwhm_l2', 'err_l2',
                                                     'fwhm_l1a', 'err_l1a']]
          .to_string(index=False, float_format='%.3f'))
    for ch in ('E', 'Ld'):
        g = v[(v['route'] == 'L2 route, E via E grid') & (v['channel'] == ch)].iloc[0]
        print('  %s offset: L2 %+.3f, L1A %+.3f nm' % (ch, g['offset_l2'], g['offset_l1a']))
    return v


def run(ref, gE, limit=None, verbose=True):
    req = pd.read_csv(REQUEST, dtype=str)
    if limit:
        req = req.groupby('row', group_keys=False).head(limit)
    temps, mods, status = [], [], []
    for i, r in req.reset_index(drop=True).iterrows():
        t0 = time.time()
        keys = dict(row=r['row'], priority=r['priority'], site_code=r['site'][:-2],
                    sequence_time=r['sequence_time'], sky=r['sky'] if isinstance(r['sky'], str)
                    else None, water_type=WATER_TYPE.get(r['site'][:-2]),
                    month=float(r['sequence_time'][4:6]))
        try:
            d = wl.load_l2b(wl.release2_path(r['site'], r['sequence_time'], r['file']))
            m = d['meta']
            sysid = m['system_id']
            cal = '%s/%s' % (m['instrument_calibration_date_rad'],
                             m['instrument_calibration_date_irr'])
            own = sysid == 'HYPSTAR_122304' and \
                m['instrument_calibration_date_irr'] == '2024-11-12'
            keys.update(system_id=sysid, cal_period=cal, sza_l1a=float(np.nanmean(d['sza'])),
                        e_grid='own' if own else 'proxy',
                        instrument_ok=sysid == r['instrument'],
                        cal_ok=cal == r['cal_dates_rad_irr'])
            t, ms = fit_l2(d['wave'], d['downwelling_radiance'][:, 0], d['irradiance'][:, 0],
                           ref, gE, dict(sequence_id=m['sequence_id'], system_id=sysid))
            temps.append(t.assign(**keys))
            tab = whn_srf.models_table(ms).drop(columns=['instrument'], errors='ignore')
            mods.append(tab.assign(**keys))
            st = 'ok'
        except Exception as ex:
            st = '%s: %s' % (type(ex).__name__, ex)
        status.append(dict(keys, fit_status=st, seconds=round(time.time() - t0, 1)))
        if verbose:
            print('[%d/%d] %s %s %s: %s' % (i + 1, len(req), r['site'], r['sequence_time'],
                                            keys.get('system_id'), st))
    cat = lambda x: pd.concat(x, ignore_index=True) if x else pd.DataFrame()  # noqa: E731
    return cat(temps), cat(mods), pd.DataFrame(status)


def add_sky(models):
    """Attach the task 8e sky index (``phase0_sky_index.csv``) if it exists."""
    p = os.path.join(WIGGLES_DIR, 'phase0_sky_index.csv')
    if not os.path.exists(p):
        return models
    sk = pd.read_csv(p, dtype={'sequence_time': str})
    sk['site_code'] = sk['site'].str[:-2]
    m = models.drop(columns=[c for c in ('sky', 'ld_ed_750') if c in models])
    return m.merge(sk[['site_code', 'sequence_time', 'ld_ed_750', 'sky']],
                   on=['site_code', 'sequence_time'], how='left')


def stability(models):
    """The task 11 analysis on the Release 2 models (Ld and E)."""
    models = add_sky(models)
    s = stab.per_sequence(models, 'template')
    st, tr, pr = stab.group_stats(s), stab.trends(s), stab.pairs(s)
    st = st[st['channel'].isin(['E', 'Ld'])]
    # E-grid provenance per group
    eg = models.groupby(['system_id', 'cal_period'])['e_grid'].first().reset_index()
    st = st.merge(eg, on=['system_id', 'cal_period'], how='left')
    st.to_csv(os.path.join(WIGGLES_DIR, 'phase0c_stability.csv'), index=False,
              float_format='%.4f')
    pr.to_csv(os.path.join(WIGGLES_DIR, 'phase0c_stability_pairs.csv'), index=False,
              float_format='%.4f')
    tr.to_csv(os.path.join(WIGGLES_DIR, 'phase0c_stability_trends.csv'), index=False,
              float_format='%.4f')
    os.makedirs(FIGDIR, exist_ok=True)
    stab.figures(s, prefix='phase0c_stability', label='Release 2 L2B, template fits')
    return s, st, tr, pr


def summary(s):
    """Per instrument x cal period: FWHM_E, FWHM_Ld and Ld - E at 450 and
    600 nm, absolute offsets (vs HSRS in air) and the L - E offset; mean
    and std over the sequences."""
    rows = []
    for (inst, cal), g in s.groupby(['system_id', 'cal_period']):
        r = dict(system_id=inst, cal_period=cal, n=len(g),
                 sites=','.join(sorted(g['site_code'].unique())))
        for q in ('fwhm_450_E', 'fwhm_450_Ld', 'dfwhm_450_Ld_E', 'fwhm_600_E', 'fwhm_600_Ld',
                  'offset_E', 'offset_Ld', 'doffset_Ld_E'):
            if q == 'dfwhm_450_Ld_E' or q in g:
                v = g[q].astype(float)
                r[q] = v.mean()
                r[q + '_std'] = v.std(ddof=1) if len(v) > 1 else np.nan
        r['dfwhm_600_Ld_E'] = (g['fwhm_600_Ld'] - g['fwhm_600_E']).mean()
        rows.append(r)
    return pd.DataFrame(rows)


def figure_by_instrument(s, path):
    """FWHM_Ld and FWHM_E at 450 and 600 nm per instrument x cal period."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    s = s.copy()
    s['group'] = s['system_id'].str.replace('HYPSTAR_', '') + '\n' + \
        s['cal_period'].str.split('/').str[-1].str[:7]
    groups = sorted(s['group'].unique())
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for ax, lam in zip(axs, (450, 600)):
        for j, (ch, mk, col) in enumerate((('E', 'o', SLOTS[0]), ('Ld', 's', SLOTS[1]))):
            for i, g in enumerate(groups):
                v = s.loc[s['group'] == g, 'fwhm_%d_%s' % (lam, ch)].dropna()
                x = i + (j - 0.5) * 0.3 + np.random.default_rng(i).uniform(-0.06, 0.06, len(v))
                ax.plot(x, v, ls='', marker=mk, ms=5, alpha=0.75, color=col,
                        label=ch if i == 0 else None)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(groups, fontsize=8)
        ax.set_title('FWHM at %d nm' % lam, fontsize=10, loc='left')
        ax.grid(alpha=0.25, lw=0.6)
        ax.spines[['top', 'right']].set_visible(False)
    axs[0].set_ylabel('template SRF FWHM (nm)')
    axs[0].legend(fontsize=9, frameon=False)
    fig.suptitle('Release 2 L2B: SRF per instrument and IRR calibration date '
                 '(E for instruments other than 122304 post-recal uses a proxy E grid)',
                 fontsize=10, x=0.06, ha='left')
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description='Phase 0 task 8c: Release 2 SRF stability')
    ap.add_argument('--validate', action='store_true', help='VEIT L2 vs L1A check only')
    ap.add_argument('--limit', type=int, default=None, help='first N per request row')
    ap.add_argument('--from-parquet', action='store_true',
                    help='skip the fits; re-run the analysis on release2_models.parquet')
    a = ap.parse_args(argv)
    ref = whn_srf.load_reference('air')
    gE = e_grid()
    if a.validate:
        validate(ref, gE)
        return
    if a.from_parquet:
        models = pd.read_parquet(os.path.join(OUT, 'release2_models.parquet'))
    else:
        temps, models, status = run(ref, gE, limit=a.limit)
        os.makedirs(OUT, exist_ok=True)
        temps.to_parquet(os.path.join(OUT, 'release2_template.parquet'))
        models.to_parquet(os.path.join(OUT, 'release2_models.parquet'))
        status.to_csv(os.path.join(OUT, 'release2_status.csv'), index=False)
        print('\n%d sequences: %d ok' % (len(status), (status['fit_status'] == 'ok').sum()))
    s, st, tr, pr = stability(models)
    summ = summary(s)
    summ.to_csv(os.path.join(WIGGLES_DIR, 'phase0c_summary.csv'), index=False,
                float_format='%.4f')
    print('\nper instrument x calibration (mean, std over sequences):')
    print(summ.to_string(index=False, float_format='%.3f'))
    figure_by_instrument(s, os.path.join(FIGDIR, 'phase0c_fwhm_by_instrument.png'))
    pd.set_option('display.width', 220)
    cols = ['system_id', 'cal_period', 'channel', 'n', 'e_grid'] + \
        ['fwhm_%d_mean' % l for l in (400, 500, 600, 700)] + \
        ['fwhm_%d_std' % l for l in (400, 500, 600, 700)] + ['passes_0p2']
    print(st[cols].to_string(index=False, float_format='%.3f'))
    done = pr.dropna(subset=['z']) if 'z' in pr else pr.iloc[0:0]
    print('\npairs:')
    print(done[['pair', 'quantity', 'n_a', 'n_b', 'mean_a', 'mean_b', 'diff', 'err', 'z']]
          .to_string(index=False, float_format='%.3f'))
    print('\ntrends (slopes per deg SZA / per season unit / per year):')
    print(tr[['system_id', 'cal_period', 'quantity', 'n', 'slope_sza', 'slope_sza_err',
              'slope_season', 'slope_season_err', 'slope_year', 'slope_year_err'] +
             [c for c in ('slope_sky', 'slope_sky_err') if c in tr]]
          .to_string(index=False, float_format='%.4f'))


if __name__ == '__main__':
    main()

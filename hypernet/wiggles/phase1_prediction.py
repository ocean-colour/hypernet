"""Phase 1, task 11: the prediction for Phase 3 (Gate G3).

Case (iii) -- real E and L grids (122304), the 122304 Ld SRF model as the L
SRF -- with a synthetic E SRF narrower by ``DFWHM`` (FWHM_E = FWHM_L - dFWHM
at every lambda, 0 to 1 nm), over the 25 OSOAA scenes, for ``linear``,
``ruddick2023`` and ``srf``.  Three scenarios of what the correction knows:

- ``ideal``: noise-free, the correction told the true SRFs and wavelengths;
- ``noise``: the case (vii) VEIT per-scan noise / sqrt(6) on E, Lu and Ld;
- ``realistic``: ``noise`` plus a +0.05 nm E wavelength shift and the
  correction told FWHM_E + 0.05 nm (the edges of Phase 0's wavelength-scale
  and SRF-stability findings).

Per scene, method and line (the ten lines of ``phase1_twin.LINES10``, +-5 nm):

- the twin metric: rms(rho_w'' - rho_w''_true) near the line, and its
  reduction relative to ``linear`` (``reduction``);
- the **observable** metric that Phase 3 can measure without a truth: the
  line-region excess, sqrt(rms_line(rho_w'')^2 - rms_away(rho_w'')^2) (0 if
  negative), and its reduction relative to linear (``reduction_obs``).  The
  truth has an excess of its own (real line-filling and gas-band peaks), so
  ``reduction_obs_max`` = 1 - excess_true / excess_linear is the most a
  perfect correction can show;
- the line depth: 1 - min(Ed_L) within +-1 nm of the line / max(Ed_L) within
  +-5 nm, of the true Ed through the L SRF.  Gas bands deepen with SZA.

Also ``line = 'all'``: the ten windows together (400-900 nm, as task 9).

Outputs:

- ``$OS_COLOR/hypernet/wiggles/phase1/twin_prediction.parquet`` (every row);
- ``hypernet/wiggles/phase1_prediction.csv``: medians over scenes per
  scenario x dFWHM x method x line (with the median depth), and per depth
  bin (``line = 'depth[a,b)'``);
- ``hypernet/wiggles/figs/phase1/twin_prediction.png``.

Run from the repository root:
``python -m hypernet.wiggles.phase1_prediction [--workers 12] [--scenes N]``.
"""
import argparse
import os
import time

import numpy as np
import pandas as pd

from hypernet import edinterp, emod, twin
from hypernet.wiggles import WIGGLES_DIR
from hypernet.wiggles.phase1_twin import (FIG1, HS, L_RANGE, LINES10, ROOT, WIN, masks,
                                          second_diff)

INSTRUMENT = 'HYPSTAR_122304'
DFWHM = (0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0)
METHODS = ('linear', 'ruddick2023', 'srf')
SCENARIOS = {'ideal': dict(noise=False, shift=0.0, told=0.0),
             'noise': dict(noise=True, shift=0.0, told=0.0),
             'realistic': dict(noise=True, shift=0.05, told=0.05)}
DEPTH_BINS = (0.0, 0.1, 0.2, 0.4, 0.6, 1.0)


def _rms(a):
    return float(np.sqrt(np.mean(a ** 2))) if a.size else np.nan


def _excess(near, away):
    return float(np.sqrt(max(near ** 2 - away ** 2, 0.0)))


def line_depths(grid, ed):
    """1 - min(ed, +-1 nm) / max(ed, +-5 nm) for each of ``LINES10``."""
    out = {}
    for n, l0 in LINES10.items():
        core = np.abs(grid - l0) <= 1.0
        win = np.abs(grid - l0) <= WIN
        out[n] = float(1 - ed[core].min() / ed[win].max()) if core.any() else np.nan
    return out


def run_scene(args):
    i_case, seed = args
    fields = twin.load_fields()
    lib = twin.rhow_library(fields)
    grids = twin.load_grids()
    sL = twin.load_instrument_srfs()[INSTRUMENT]['Ld']
    sza = float(fields['sza'][i_case])
    scene = twin.compose_scene(emod.build_emod(sza), fields, i_case, lib=lib)
    lam = scene['lam']
    name = str(np.asarray(fields['name'])[i_case])
    gL = grids['grid_L_122304']
    gL = gL[(gL >= L_RANGE[0]) & (gL <= L_RANGE[1])]
    gE = grids['grid_E_122304']
    mE = (gE >= L_RANGE[0] - 3) & (gE <= L_RANGE[1] + 3)
    gE = gE[mE]
    nE = (grids['noise_E_scan'] / np.sqrt(6))[mE]
    nLu = np.interp(gL, grids['grid_L_122304'], grids['noise_Lu_scan']) / np.sqrt(6)
    nLd = np.interp(gL, grids['grid_L_122304'], grids['noise_Ld_scan']) / np.sqrt(6)
    truth = twin.rhow_true(scene, sL, gL)
    rho_eff = np.interp(gL, lam, scene['rho_eff'])
    near, away, per = masks(gL)
    depth = line_depths(gL, twin.observe(lam, scene['Ed'], gL, sL))
    Lu0 = twin.observe(lam, scene['Lu'], gL, sL)
    Ld0 = twin.observe(lam, scene['Ld'], gL, sL)
    mdl = (lam, scene['Ed'])
    meta = dict(scene=name, sza=sza, aot550=float(fields['aot550'][i_case]),
                water=str(np.asarray(fields['water'])[i_case]))
    rows = []
    for k_s, (sc, p) in enumerate(SCENARIOS.items()):
        s0 = seed * 1000 + 10 * k_s
        Lu = twin.observe(lam, scene['Lu'], gL, sL, noise=nLu, seed=s0 + 1) if p['noise'] else Lu0
        Ld = twin.observe(lam, scene['Ld'], gL, sL, noise=nLd, seed=s0 + 2) if p['noise'] else Ld0
        for d in DFWHM:
            sE = (lambda l, d=d: sL.fwhm(l) - d)
            sE_told = (lambda l, d=d, t=p['told']: sL.fwhm(l) - d + t)
            Eobs = twin.observe(lam, scene['Ed'], gE, sE, true_grid=gE + p['shift'],
                                noise=nE if p['noise'] else None, seed=s0)
            for h in HS:
                t2 = second_diff(gL, truth, h)
                tr_away = _rms(t2[away])
                base = dict(**meta, scenario=sc, dfwhm=d, h=h)
                lin = {}
                for m in METHODS:
                    EdL = edinterp.interpolate_ed_to_l(gE, Eobs, gL, emod=mdl, srf_irr=sE_told,
                                                       srf_rad=sL.fwhm, method=m)
                    r2 = second_diff(gL, np.pi * (Lu - rho_eff * Ld) / EdL, h)
                    err, r_away = r2 - t2, _rms(r2[away])
                    for n, msk in [('all', near)] + list(per.items()):
                        r = dict(base, method=m, line=n,
                                 depth=depth.get(n, np.nan),
                                 rms_err=_rms(err[msk]),
                                 excess_obs=_excess(_rms(r2[msk]), r_away),
                                 excess_true=_excess(_rms(t2[msk]), tr_away))
                        if m == 'linear':
                            lin[n] = r
                        L = lin[n]
                        r['reduction'] = 1 - r['rms_err'] / L['rms_err']
                        r['reduction_obs'] = (1 - r['excess_obs'] / L['excess_obs']
                                              if L['excess_obs'] > 0 else np.nan)
                        r['reduction_obs_max'] = (1 - r['excess_true'] / L['excess_obs']
                                                  if L['excess_obs'] > 0 else np.nan)
                        rows.append(r)
    return rows


def summarise(df):
    keys = ['scenario', 'dfwhm', 'h', 'method', 'line']
    cols = ['depth', 'rms_err', 'excess_obs', 'excess_true', 'reduction', 'reduction_obs',
            'reduction_obs_max']
    s = df.groupby(keys)[cols].median().reset_index()
    s['n'] = df.groupby(keys).size().values
    # per depth bin, over individual lines in individual scenes
    d = df[df['line'] != 'all'].copy()
    d['line'] = pd.cut(d['depth'], DEPTH_BINS, right=False).map(
        lambda iv: 'depth[%.1f,%.1f)' % (iv.left, iv.right)).astype(str)
    b = d.groupby(keys)[cols].median().reset_index()
    b['n'] = d.groupby(keys).size().values
    return pd.concat([s, b], ignore_index=True)


def figure(summ):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    cols = {'linear': '#e34948', 'ruddick2023': '#2a78d6', 'srf': '#1baf7a'}
    lss = {'ideal': '-', 'noise': '--', 'realistic': ':'}
    s = summ[summ['h'] == 1.0]
    fig, axs = plt.subplots(1, 3, figsize=(16, 5.2))
    ax = axs[0]
    for sc, ls in lss.items():
        for m in ('ruddick2023', 'srf'):
            q = s[(s['scenario'] == sc) & (s['method'] == m) & (s['line'] == 'all')]
            ax.plot(q['dfwhm'], 100 * q['reduction'], ls=ls, marker='o', ms=4, lw=1.6,
                    color=cols[m], label='%s, %s' % (m, sc))
    ax.axhline(80, color='#888', lw=1, ls='-.')
    ax.set_ylim(-150, 105)
    ax.set_xlabel('ΔFWHM = FWHM_L − FWHM_E (nm)')
    ax.set_ylabel('reduction of the line-region ρw\'\' error vs linear (%)')
    ax.set_title('(a) twin metric, ten lines together', fontsize=10, loc='left')
    ax.legend(fontsize=7.5, frameon=False, ncol=2, loc='lower right')

    ax = axs[1]
    for sc, ls in lss.items():
        q = s[(s['scenario'] == sc) & (s['method'] == 'srf') & (s['line'] == 'all')]
        ax.plot(q['dfwhm'], 100 * q['reduction_obs'], ls=ls, marker='o', ms=4, lw=1.6,
                color=cols['srf'], label='srf, %s' % sc)
        ax.plot(q['dfwhm'], 100 * q['reduction_obs_max'], ls=ls, lw=1, color='#333',
                label='perfect correction, %s' % sc)
    ax.axhline(80, color='#888', lw=1, ls='-.')
    ax.set_ylim(-10, 105)
    ax.set_xlabel('ΔFWHM (nm)')
    ax.set_ylabel('reduction of the observable excess vs linear (%)')
    ax.set_title('(b) observable: √(rms²_line − rms²_away) of ρw\'\'', fontsize=10, loc='left')
    ax.legend(fontsize=7.5, frameon=False, loc='lower right')

    ax = axs[2]
    # per line, not per depth bin: with noise, depth is confounded with
    # wavelength (the deep H2O and O2-A bands sit where Lu is noisiest)
    q = s[(s['scenario'] == 'noise') & (s['method'] == 'srf') & s['line'].isin(list(LINES10))]
    grid = q.pivot_table(index='line', columns='dfwhm', values='reduction').reindex(
        sorted(LINES10, key=LINES10.get))
    dep = q[q['dfwhm'] == 0.5].set_index('line')['depth']
    im = ax.imshow(100 * grid.values, aspect='auto', origin='lower', cmap='RdYlGn',
                   vmin=-100, vmax=100)
    for (i, j), v in np.ndenumerate(grid.values):
        if np.isfinite(v):
            ax.text(j, i, '%.0f' % (100 * v), ha='center', va='center', fontsize=7.5)
    ax.set_xticks(range(grid.shape[1]))
    ax.set_xticklabels(['%g' % c for c in grid.columns], fontsize=8)
    ax.set_yticks(range(grid.shape[0]))
    ax.set_yticklabels(['%s (depth %.2f)' % (n, dep[n]) for n in grid.index], fontsize=8)
    ax.set_xlabel('ΔFWHM (nm)')
    ax.set_title('(c) srf, with noise: reduction (%) per line', fontsize=10, loc='left')
    fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
    for a in axs[:2]:
        a.grid(alpha=0.25, lw=0.6)
        a.spines[['top', 'right']].set_visible(False)
    fig.suptitle('Phase 1 prediction for Phase 3: case (iii), 122304 grids and L SRF, '
                 'FWHM_E = FWHM_L − ΔFWHM; medians over 25 scenes, h = 1 nm',
                 fontsize=10.5, x=0.05, ha='left')
    fig.tight_layout()
    fig.savefig(os.path.join(FIG1, 'twin_prediction.png'), dpi=150)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description='Phase 1 task 11: prediction for Phase 3')
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--scenes', type=int, default=None, help='first N scenes only')
    a = ap.parse_args(argv)
    fields = twin.load_fields()
    n = len(fields['name']) if a.scenes is None else a.scenes
    for sza in sorted(set(float(x) for x in fields['sza'][:n])):
        emod.build_emod(sza)
    t0 = time.time()
    from multiprocessing import Pool
    rows = []
    with Pool(a.workers) as pool:
        for r in pool.imap_unordered(run_scene, [(i, i) for i in range(n)]):
            rows += r
    df = pd.DataFrame(rows)
    os.makedirs(ROOT, exist_ok=True)
    df.to_parquet(os.path.join(ROOT, 'twin_prediction.parquet'))
    summ = summarise(df)
    summ.to_csv(os.path.join(WIGGLES_DIR, 'phase1_prediction.csv'), index=False,
                float_format='%.4g')
    figure(summ)
    pd.set_option('display.width', 220)
    print('%d rows in %.0f s' % (len(df), time.time() - t0))
    s = summ[(summ['h'] == 1.0) & (summ['method'] != 'linear')]
    for col in ('reduction', 'reduction_obs', 'reduction_obs_max'):
        print('\n%s, ten lines together (h = 1 nm):' % col)
        print(s[s['line'] == 'all'].pivot_table(index=['scenario', 'method'], columns='dfwhm',
                                                values=col).to_string(float_format='%.2f'))
    print('\nlinear rms error and observable excess, all lines (h = 1 nm, ideal):')
    q = summ[(summ['h'] == 1.0) & (summ['method'] == 'linear') & (summ['line'] == 'all') &
             (summ['scenario'] == 'ideal')]
    print(q[['dfwhm', 'rms_err', 'excess_obs', 'excess_true']].to_string(
        index=False, float_format='%.3g'))
    print('\nsrf reduction by line (h = 1 nm), with median depth:')
    q = s[(s['method'] == 'srf') & ~s['line'].str.startswith('depth')]
    print(q.pivot_table(index=['scenario', 'line'], columns='dfwhm', values='reduction')
          .join(q[q['dfwhm'] == 0.5].groupby(['scenario', 'line'])['depth'].first())
          .to_string(float_format='%.2f'))
    print('\nsrf reduction by depth bin (h = 1 nm):')
    q = s[(s['method'] == 'srf') & s['line'].str.startswith('depth')]
    print(q.pivot_table(index=['scenario', 'line'], columns='dfwhm', values='reduction')
          .to_string(float_format='%.2f'))


if __name__ == '__main__':
    main()

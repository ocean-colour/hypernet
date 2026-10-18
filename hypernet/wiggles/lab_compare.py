"""Lab vs field FWHM (lab data prompts, task 2).

Compares the lab lamp-line FWHMs (``hypernet/data/lab_srf.csv``) with the
per-instrument SRF table (``hypernet/data/srf_table.csv``) at the line
wavelengths inside each table row's fitted range (395-865 nm).

The lab epoch of a table row (column ``lab_epoch``, Setup Q8) is the latest
non-``_tst`` epoch with coefficient date <= that channel's calibration date
(``cal_date_irr`` for E, ``cal_date_rad`` for Ld): the lab state of the
instrument when it was calibrated.  It differs from the processor's
wavelength epoch (``processor_epoch``) only for 122304 at 2024-11, where it
is 2024_09.

Lab Delta FWHM = median over the 400-680 nm lines with S/N >= 10 of
FWHM_L - FWHM_E at the same epoch (the band of the G1 rule,
``edinterp.DFWHM_BAND``); ``lab_dfwhm_noblend`` drops Hg 404.66.

Outputs (committed):

- ``hypernet/wiggles/lab_fwhm_compare.csv`` -- one row per table row x line;
- ``hypernet/wiggles/lab_fwhm_compare_summary.csv`` -- one row per table row:
  median lab - table FWHM, lab and table Delta FWHM, both G1 methods;
- ``hypernet/wiggles/figs/lab/fwhm_lab_vs_table.png``.

Intermediates: ``$OS_COLOR/hypernet/wiggles/lab_data/lab_tables/``
(``lab_dfwhm_epochs.csv``: lab Delta FWHM of every epoch;
``lab_shape_clean.csv``: median shape parameters on the clean lines).

Run from the repository root:  ``python -m hypernet.wiggles.lab_compare``
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from hypernet import edinterp, labcal, srftable
from hypernet import srf as srfmod

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, 'figs', 'lab')
OUT = os.path.join(labcal.lab_root(), 'lab_tables')
CLEAN = (404.6563, 435.8328, 546.0735, 763.5106, 826.4522, 881.9411)
BAND = edinterp.DFWHM_BAND
COL = {'E': '#2a78d6', 'Ld': '#eb6834'}


def _usable(g):
    return g[~g['low_snr'] & ~g['outside_grid']]


def lab_dfwhm(srf, sn, epoch, noblend=False):
    """Median FWHM_L - FWHM_E over the usable lines in BAND."""
    g = srf[(srf['sn'] == sn) & (srf['epoch'] == epoch)]
    p = g.pivot_table(index='lam_ref', columns='channel', values='lab_fwhm')
    ok = g.groupby('lam_ref')['low_snr'].any()
    m = (p.index >= BAND[0]) & (p.index <= BAND[1]) & ~ok.reindex(p.index).values
    if noblend:
        m &= ~np.isin(p.index, labcal.BLENDS)
    d = (p['L'] - p['E'])[m]
    return float(d.median()), int(m.sum())


def compare():
    tab = srftable.load_srf_table()
    srf = labcal.load_lab_srf()
    rows, summ = [], []
    for _, r in tab.iterrows():
        sn = r['instrument'].split('_')[1]
        ch = 'L' if r['channel'] == 'Ld' else 'E'
        cal = r['cal_date_irr'] if ch == 'E' else r['cal_date_rad']
        ep, before = labcal.physical_epoch(sn, cal)
        g = srf[(srf['sn'] == sn) & (srf['epoch'] == ep) & (srf['channel'] == ch)]
        g = g[(g['lam_ref'] >= r['lam_min']) & (g['lam_ref'] <= r['lam_max'])]
        model = srfmod.fwhm_at(g['lam_ref'].values, [r['c0'], r['c1'], r['c2']],
                               r['lam_ref'], r['lam_scale'])
        keys = dict(instrument=r['instrument'], channel=r['channel'], cal_period=r['cal_period'],
                    source=r['source'], lab_epoch=ep, lab_epoch_before_cal=before,
                    processor_epoch=labcal.PROCESSOR_EPOCH[sn]['epoch'])
        d = pd.DataFrame(dict(keys, lam=g['lam_ref'].values, species=g['species'].values,
                              lab_fwhm=g['lab_fwhm'].values, table_fwhm=model,
                              lab_minus_table=g['lab_fwhm'].values - model,
                              our_gauss_fwhm=g['fwhm_gauss'].values,
                              our_sg_fwhm=g['fwhm_sg'].values, sg_p=g['sg_p'].values,
                              snr=g['snr'].values, blend=g['blend'].values,
                              low_snr=g['low_snr'].values))
        rows.append(d)
        u = d[~d['low_snr']]
        ub = u[~u['blend']]
        dl, n_band = lab_dfwhm(srf, sn, ep)
        dlb, _ = lab_dfwhm(srf, sn, ep, noblend=True)
        summ.append(dict(keys, n_lines=len(u), lab_minus_table_med=u['lab_minus_table'].median(),
                         lab_minus_table_min=u['lab_minus_table'].min(),
                         lab_minus_table_max=u['lab_minus_table'].max(),
                         lab_minus_table_med_noblend=ub['lab_minus_table'].median(),
                         n_band_lines=n_band, lab_dfwhm=dl, lab_dfwhm_noblend=dlb,
                         table_dfwhm=r['delta_fwhm'], dfwhm_lab_minus_table=dl - r['delta_fwhm'],
                         table_method=r['method'],
                         lab_method=('srf' if dl >= edinterp.SRF_THRESHOLD else 'linear')))
    per_line = pd.concat(rows, ignore_index=True)
    summ = pd.DataFrame(summ)
    per_line.to_csv(os.path.join(HERE, 'lab_fwhm_compare.csv'), index=False, float_format='%.4f')
    summ.to_csv(os.path.join(HERE, 'lab_fwhm_compare_summary.csv'), index=False,
                float_format='%.4f')
    return tab, srf, per_line, summ


def epoch_tables(srf):
    os.makedirs(OUT, exist_ok=True)
    ep = []
    for (sn, e), _ in srf.groupby(['sn', 'epoch']):
        dl, n = lab_dfwhm(srf, sn, e)
        dlb, nb = lab_dfwhm(srf, sn, e, noblend=True)
        ep.append(dict(sn=sn, epoch=e, processor=labcal.PROCESSOR_EPOCH[sn]['epoch'] == e,
                       n_band_lines=n, lab_dfwhm=dl, lab_dfwhm_noblend=dlb))
    ep = pd.DataFrame(ep)
    ep.to_csv(os.path.join(OUT, 'lab_dfwhm_epochs.csv'), index=False, float_format='%.4f')
    cl = srf[srf['lam_ref'].isin(CLEAN) & ~srf['low_snr']]
    sh = cl.groupby(['sn', 'epoch', 'channel']).agg(
        n=('lam_ref', 'size'), sg_p=('sg_p', 'median'), sg_p_min=('sg_p', 'min'),
        sg_p_max=('sg_p', 'max'), asym=('asym', 'median'), chi2_gauss=('chi2_gauss', 'median'),
        chi2_sg=('chi2_sg', 'median'), chi2_asym=('chi2_asym', 'median')).reset_index()
    sh['sg_over_gauss'] = sh['chi2_sg'] / sh['chi2_gauss']
    sh.to_csv(os.path.join(OUT, 'lab_shape_clean.csv'), index=False, float_format='%.4f')
    return ep, sh


def figure(tab, per_line, summ):
    os.makedirs(FIG, exist_ok=True)
    insts = sorted(set(tab['instrument']))
    fig, axs = plt.subplots(2, 4, figsize=(15, 7.2))
    lam = np.linspace(395, 865, 200)
    for a, inst in zip(axs.ravel(), insts):
        g = per_line[per_line['instrument'] == inst]
        for ch in ('E', 'Ld'):
            t = tab[(tab['instrument'] == inst) & (tab['channel'] == ch)]
            for k, (_, r) in enumerate(t.iterrows()):
                ls = ':' if r['source'] == 'veit_l1a' else ('-' if k == 0 else '--')
                a.plot(lam, srfmod.fwhm_at(lam, [r['c0'], r['c1'], r['c2']], r['lam_ref'],
                                           r['lam_scale']),
                       color=COL[ch], lw=1.6, ls=ls, alpha=0.9,
                       label='table %s %s%s' % (ch, r['cal_period'][:7],
                                                ' (VEIT L1A)' if r['source'] == 'veit_l1a' else ''))
            gg = g[(g['channel'] == ch)].drop_duplicates(['lam', 'lab_epoch'])
            ok = gg[~gg['low_snr']]
            weak = gg[gg['low_snr']]
            a.scatter(ok['lam'], ok['lab_fwhm'], s=36, color=COL[ch], zorder=3,
                      edgecolor='white', linewidth=0.8,
                      label='lab %s (%s)' % (ch, '/'.join(sorted(set(gg['lab_epoch'])))))
            a.scatter(weak['lam'], weak['lab_fwhm'], s=36, facecolor='none', edgecolor=COL[ch],
                      zorder=3, linewidth=1.2)
        a.set_title(inst.replace('HYPSTAR_', 'HYPSTAR '), fontsize=10)
        a.set_ylim(1.9, 3.9)
        a.grid(alpha=0.25, lw=0.5)
        a.legend(fontsize=6, loc='upper left', frameon=False)
        a.set_xlabel('wavelength (nm)', fontsize=8)
        a.set_ylabel('FWHM (nm)', fontsize=8)
        a.tick_params(labelsize=8)
    # last panel: Delta FWHM, lab vs table, one point per E row
    a = axs.ravel()[-1]
    s = summ[summ['channel'] == 'E']
    a.plot([-0.2, 0.6], [-0.2, 0.6], color='0.6', lw=1)
    a.axhline(edinterp.SRF_THRESHOLD, color='0.6', lw=0.8, ls='--')
    a.axvline(edinterp.SRF_THRESHOLD, color='0.6', lw=0.8, ls='--')
    a.scatter(s['table_dfwhm'], s['lab_dfwhm'], s=40, color='#4a3aa7', edgecolor='white', zorder=3)
    for _, r in s.drop_duplicates(['instrument', 'cal_period']).iterrows():
        a.annotate(r['instrument'][-6:], (r['table_dfwhm'], r['lab_dfwhm']), fontsize=6,
                   xytext=(3, 3), textcoords='offset points', color='0.25')
    a.set_xlabel('table ΔFWHM = L − E (nm)', fontsize=8)
    a.set_ylabel('lab ΔFWHM, 400-680 nm (nm)', fontsize=8)
    a.set_title('ΔFWHM: lab vs table (dashed: G1 threshold)', fontsize=10)
    a.grid(alpha=0.25, lw=0.5)
    a.tick_params(labelsize=8)
    fig.suptitle('Lab lamp-line FWHM (points; open = S/N < 10) vs srf_table FWHM(λ) (lines)',
                 fontsize=11)
    fig.tight_layout()
    p = os.path.join(FIG, 'fwhm_lab_vs_table.png')
    fig.savefig(p, dpi=130)
    plt.close(fig)
    return p


def main():
    tab, srf, per_line, summ = compare()
    ep, sh = epoch_tables(srf)
    pd.set_option('display.width', 250)
    print(summ[['instrument', 'channel', 'cal_period', 'source', 'lab_epoch', 'n_lines',
                'lab_minus_table_med', 'lab_minus_table_med_noblend', 'lab_dfwhm',
                'table_dfwhm', 'table_method', 'lab_method']].round(3).to_string(index=False))
    print('\nlab Delta FWHM per epoch:')
    print(ep.round(3).to_string(index=False))
    print('\nclean-line shapes:')
    print(sh.round(3).to_string(index=False))
    print('wrote', figure(tab, per_line, summ))


if __name__ == '__main__':
    main()

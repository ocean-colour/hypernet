"""Phase 2, task 9: the three methods on the VEIT sample.

Kevin's VEIT sequence (HYPSTAR 122304, 2026-06-04 08:45, SZA 36.7):

1. **L1C reproduction.**  The L1A E series means, interpolated in time to the
   mean water-view time (``whn_l1a.interp_in_time``, as Phase 0's
   ``check_against_l1c``), then ``linear`` on the
   L grid, against the L1C ``irradiance`` over 400-900 nm.  Plan §2.3 found
   a median relative difference of 2 x 10^-5; ``linear`` must stay there.  Also the snippet's
   ``InterpolationWavSRF`` unwired (linear, extrapolating).
2. **Three methods.**  ``linear`` (as the processor: no scale correction);
   ``ruddick2023`` and ``srf`` with the VEIT template SRFs
   (``hypernet/data/veit_srf_model.json``), Emod from ``emod.emod_for(36.7)``,
   and the per-sequence wavelength-scale corrections of E and Ld (task 4b,
   common window pattern removed; both pass Setup Q8's 0.03 nm).
   ``srf_labshape`` (lab data prompts task 3) is ``srf`` with the E SRF of
   the SRF table: super-Gaussian with the lab p and the FWHM(lambda) fitted
   with that shape (``hypernet/data/veit_srf_model_lab_shape.json``); Ld
   stays Gaussian.
3. **Uncertainties** (``edinterp.propagate_ed_to_l``, 200 draws, on the
   time-mean E).  The L1C ``u_rel_*`` are all 0 (the water default is
   ``mcsteps = 0``) and the calibration uncertainties are not in the files,
   so the inputs are:
   - random: E's scan scatter (flattened, error of the mean, ``srf.scan_errors``);
   - systematic_indep: no input part, only the new terms: the SRF
     covariances, the scale-fit covariances, and the Emod inputs
     (``emod.emod_perturbations``);
   - systematic_corr_rad_irr: none.

   So the bands show what this work adds, on top of a noise floor.
4. **Line residual of Ld / Ed_L** (Ld at the water-view time): the ratio
   over its 10 nm running median, minus 1; rms within +-5 nm of the ten
   lines (``phase1_twin.LINES10``) and away from them, 400-900 nm.  Ld has
   its own Ring filling-in (Phase 0 8d), which no E method removes.

Outputs: ``hypernet/wiggles/phase2_veit_check.csv`` and
``hypernet/wiggles/figs/phase2/veit_three_methods.png``.

Run from the repository root: ``python -m hypernet.wiggles.phase2_veit_check``.
"""
import os

import numpy as np
import pandas as pd
from scipy.ndimage import median_filter

from hypernet import edinterp, emod, srf, srfcal, whn_srf
from hypernet import whn_l1a as wl
from hypernet.snippets.interpolate_wav_srf import InterpolationWavSRF
from hypernet.wiggles import DATA_DIR, FIGDIR, WIGGLES_DIR
from hypernet.wiggles.phase1_twin import LINES10, masks

FIG2 = FIGDIR.replace('phase0', 'phase2')
VEIT = dict(site='VEIT', seq_time='20260604T0845')
METHODS = ('linear', 'ruddick2023', 'srf', 'srf_labshape')
#: label -> (edinterp method, E SRF: 'gauss' = veit_srf_model.json, 'lab' = lab shape)
VARIANTS = {'linear': ('linear', None), 'ruddick2023': ('ruddick2023', 'gauss'),
            'srf': ('srf', 'gauss'), 'srf_labshape': ('srf', 'lab')}
MCSTEPS = 200


def load():
    f = wl.sequence_files(**VEIT)
    irr, rad, l1c = (wl.load_l1a_irr(f['L1A_IRR']), wl.load_l1a_rad(f['L1A_RAD']),
                     wl.load_l1c(f['L1C_ALL']))
    E_t = wl.interp_in_time(irr['scans'], irr['series_id'], irr['time'], rad['Lu']['t_mean'])
    return irr, rad, l1c, E_t


def l1c_check(irr, rad, l1c, E_t):
    w = rad['wave']
    band = (w >= 400) & (w <= 900)
    lin = edinterp.interpolate_ed_to_l(irr['wave'], E_t, w, method='linear')
    mf = InterpolationWavSRF()
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        mf.setup(float(np.mean(irr['sza'])), 'w', irr['wave'])
    snip = mf.meas_function(w, irr['wave'], E_t)
    def rel(a):
        a, b = np.broadcast_arrays(np.asarray(a)[:, None] if np.ndim(a) == 1 else a,
                                   l1c['irradiance'])
        d = np.abs(a[band] / b[band] - 1)
        return float(np.nanmax(d)), float(np.nanmedian(d))
    return dict(edinterp_linear=rel(lin), snippet_unwired=rel(snip))


def run():
    irr, rad, l1c, E_t = load()
    chk = l1c_check(irr, rad, l1c, E_t)
    models = srf.load_srf_models(os.path.join(DATA_DIR, 'veit_srf_model.json'))
    sE, sL = models['E'], models['Ld']
    sE_lab = srf.load_srf_models(os.path.join(DATA_DIR, 'veit_srf_model_lab_shape.json'))['E']
    pattern = srfcal.load_window_pattern()
    mean, e_mean, _ = srf.scan_errors(irr['scans'], flatten_px=whn_srf.FLATTEN_PX)
    scale_E = srfcal.fit_wavelength_scale(irr['wave'], mean, srf=sE, err=e_mean, pattern=pattern)
    ld_mean, ld_err, _ = srf.scan_errors(rad['Ld']['scans'], flatten_px=whn_srf.FLATTEN_PX)
    scale_L = srfcal.fit_wavelength_scale(rad['wave'], ld_mean, srf=sL, err=ld_err,
                                          pattern=pattern)
    sza = float(np.mean(irr['sza']))
    entry = emod.emod_for(sza)
    pert = emod.emod_perturbations(entry)
    w = rad['wave']
    E1 = E_t if E_t.ndim == 1 else E_t.mean(axis=1)
    u_rand = 100.0 * e_mean / np.abs(mean)
    Ld = rad['Ld']['at_lu_time']
    near, away, per = masks(w)
    res, rows = {}, []
    for m in METHODS:
        meth, which = VARIANTS[m]
        model = meth != 'linear'
        kw = dict(emod=entry, srf_irr=sE_lab if which == 'lab' else sE, srf_rad=sL) \
            if model else {}
        sc = dict(scale_irr=scale_E, scale_rad=scale_L) if model else {}
        ed, info = edinterp.propagate_ed_to_l(
            irr['wave'], E1, w, method=meth, u_rel_random=u_rand, mcsteps=MCSTEPS, seed=42,
            emod_perturbations=pert if model else None, **kw, **sc)
        ratio = Ld / ed
        hp = ratio / median_filter(ratio, size=21, mode='nearest') - 1
        u_tot = np.sqrt(info['u_rel_random'] ** 2 + info['u_rel_systematic_indep'] ** 2)
        res[m] = dict(ed=ed, ratio=ratio, hp=hp, info=info, u_tot=u_tot)
        b = (w >= 400) & (w <= 900)
        r = dict(method=m,
                 max_rel_diff_vs_l1c=chk['edinterp_linear'][0] if m == 'linear' else np.nan,
                 median_rel_diff_vs_l1c=chk['edinterp_linear'][1] if m == 'linear' else np.nan,
                 rms_hp_near_lines=float(np.sqrt(np.mean(hp[near] ** 2))),
                 rms_hp_away=float(np.sqrt(np.mean(hp[away] ** 2))),
                 u_random_med=float(np.median(info['u_rel_random'][b])),
                 u_sys_indep_med=float(np.median(info['u_rel_systematic_indep'][b])),
                 u_sys_indep_near_lines_med=float(np.median(info['u_rel_systematic_indep'][near])),
                 u_sys_indep_max=float(np.max(info['u_rel_systematic_indep'][b])),
                 u_emod_med=float(np.median(info['u_rel_emod'][b])),
                 u_emod_max=float(np.max(info['u_rel_emod'][b])),
                 n_fallback=int(info['fallback'].sum()))
        r['line_over_away'] = r['rms_hp_near_lines'] / r['rms_hp_away']
        for n, msk in per.items():
            r['rms_hp_' + n.replace(' ', '_')] = (float(np.sqrt(np.mean(hp[msk] ** 2)))
                                                  if msk.any() else np.nan)
        rows.append(r)
    t = pd.DataFrame(rows)
    meta = dict(snippet_unwired_max_rel_diff=chk['snippet_unwired'],
                scale_E=(scale_E['shift'], scale_E['stretch'], scale_E['max_err']),
                scale_L=(scale_L['shift'], scale_L['stretch'], scale_L['max_err']),
                delta_fwhm=edinterp.delta_fwhm(sE, sL),
                delta_fwhm_labshape=edinterp.delta_fwhm(sE_lab, sL), p_lab=sE_lab.shape_par,
                sza=sza, sza_bin=entry['sza_bin'])
    return w, Ld, res, t, meta


def figure(w, Ld, res, meta):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    os.makedirs(FIG2, exist_ok=True)
    cols = {'linear': '#e34948', 'ruddick2023': '#2a78d6', 'srf': '#1baf7a',
            'srf_labshape': '#4a3aa7'}
    fig = plt.figure(figsize=(14, 9.5))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.1, 1, 0.9])
    ax0 = fig.add_subplot(gs[0, :])
    for m in METHODS:
        ax0.plot(w, res[m]['ratio'], color=cols[m], lw=0.8, label=m)
    ax0.set_xlim(385, 900)
    ax0.set_ylabel('Ld / Ed_L (sr⁻¹)')
    ax0.set_title('VEIT 2026-06-04 08:45, HYPSTAR 122304 (ΔFWHM %.2f nm): sky radiance over '
                  'irradiance on the L grid' % meta['delta_fwhm'], fontsize=10, loc='left')
    ax0.legend(fontsize=9, frameon=False, ncol=3)
    for k, (lab, lo, hi) in enumerate((('Ca H/K', 388, 402), ('Hα', 650, 663),
                                       ('Na D', 583, 596))):
        ax = fig.add_subplot(gs[1, k])
        m_ = (w > lo) & (w < hi)
        for m in METHODS:
            r, u = res[m]['ratio'], res[m]['u_tot'] / 100.0
            ax.fill_between(w[m_], (r * (1 - u))[m_], (r * (1 + u))[m_], color=cols[m],
                            alpha=0.15, lw=0)
            ax.plot(w[m_], r[m_], color=cols[m], lw=1.3)
        ax.set_title('%s (bands: ±1σ of Ed_L)' % lab, fontsize=9.5, loc='left')
        ax.set_xlabel('wavelength (nm)')
    fig.axes[1].set_ylabel('Ld / Ed_L')
    ax3 = fig.add_subplot(gs[2, :])
    info = res['srf']['info']
    for key, c, lab in (('u_rel_random', '#888888', 'random (E scan scatter)'),
                        ('u_rel_systematic_indep', '#1baf7a',
                         'systematic_indep: SRF + wavelength scale + Emod'),
                        ('u_rel_emod', '#eda100', 'of which Emod (SZA, PWV, O₃)')):
        ax3.semilogy(w, np.clip(info[key], 1e-5, None), color=c, lw=0.8, label=lab)
    ax3.set_xlim(385, 900)
    ax3.set_ylabel('u_rel of Ed_L, srf (%)')
    ax3.set_xlabel('wavelength (nm)')
    ax3.legend(fontsize=8.5, frameon=False, ncol=3, loc='upper right')
    for ax in fig.axes:
        ax.grid(alpha=0.25, lw=0.6)
        ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    path = os.path.join(FIG2, 'veit_three_methods.png')
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def main():
    w, Ld, res, t, meta = run()
    t.to_csv(os.path.join(WIGGLES_DIR, 'phase2_veit_check.csv'), index=False,
             float_format='%.4g')
    path = figure(w, Ld, res, meta)
    pd.set_option('display.width', 220)
    print('L1C reproduction, 400-900 nm (plan §2.3: median 2e-5): edinterp linear median '
          '%.2e, max %.2e; snippet unwired median %.2e, max %.2e'
          % (t.loc[0, 'median_rel_diff_vs_l1c'], t.loc[0, 'max_rel_diff_vs_l1c'],
             meta['snippet_unwired_max_rel_diff'][1], meta['snippet_unwired_max_rel_diff'][0]))
    print('scale E shift %+.3f stretch %+.3f (max_err %.3f); L shift %+.3f stretch %+.3f '
          '(max_err %.3f) nm' % (*meta['scale_E'], *meta['scale_L']))
    print('delta_fwhm: Gaussian E %.3f, lab-shape E %.3f (p %.3f)'
          % (meta['delta_fwhm'], meta['delta_fwhm_labshape'], meta['p_lab']))
    lin = t.set_index('method')['line_over_away']
    for m in METHODS[1:]:
        print('  %-13s line/away %.3f: %.0f %% of the linear excess removed'
              % (m, lin[m], 100 * (lin['linear'] - lin[m]) / (lin['linear'] - 1)))
    print(t[['method', 'rms_hp_near_lines', 'rms_hp_away', 'line_over_away', 'u_random_med',
             'u_sys_indep_med', 'u_sys_indep_near_lines_med', 'u_sys_indep_max', 'u_emod_med',
             'u_emod_max', 'n_fallback']].to_string(index=False, float_format='%.4g'))
    print(t[['method'] + [c for c in t.columns if c.startswith('rms_hp_') and
                          c not in ('rms_hp_near_lines', 'rms_hp_away')]]
          .to_string(index=False, float_format='%.4f'))
    print('wrote', path)


if __name__ == '__main__':
    main()

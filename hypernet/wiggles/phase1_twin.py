"""Phase 1, task 9: the twin experiment, cases x methods.

For each instrument (122304 narrow-E, 122305 E ~ L; Q&A Q6), each instrument
case of ``hypernet.twin.case_table`` (i-vii, 16 in all), each OSOAA scene (25:
SZA, aerosol, water) and each method of ``hypernet.edinterp`` (linear,
cubic, sinc, ruddick2023, srf):

1. **Scene.**  ``twin.compose_scene`` at 0.01 nm (Emod for the scene's SZA),
   with the fluorescence and Raman controls.
2. **Truth.**  ``twin.rhow_true``: Lw and Ed through the L SRF on the L grid.
3. **Measurements.**
   - Lu and Ld through the L SRF on the L grid;
   - E through the case's E SRF at the case's *true* E grid, reported on the
     nominal E grid (case iv);
   - case vii adds the VEIT per-scan noise / sqrt(6) to E, Lu and Ld.
4. **Correction.**  ``interpolate_ed_to_l`` with the SRFs the correction is
   told (case v: +-0.3 nm) and its Emod: the scene's own high-resolution Ed,
   i.e. a perfect model, or the case vi variant (SZA + 10, water vapour x 2,
   or the low-aerosol fields).
5. **rho_w** = pi (Lu - rho_eff Ld) / Ed_L, with rho_eff the scene's.
6. **rho_w''** by second differences, [rho(l + h) - 2 rho(l) + rho(l - h)] /
   h^2, for h = 1 and 5 nm (as in Ruddick et al. 2023), on the L grid.
7. **Metrics.**  rms(rho_w'' - rho_w''_true) within +-5 nm of the plan's
   ten lines (``LINES10``) and away from them, over 400-900 nm; also per
   line.

Outputs:

- ``$OS_COLOR/hypernet/wiggles/phase1/twin_metrics.parquet`` (every
  combination);
- ``hypernet/wiggles/phase1_twin_metrics.csv`` (medians over scenes per
  instrument x case x method x h, with the reduction relative to linear);
- ``hypernet/wiggles/figs/phase1/twin_rho2_lines.png`` and
  ``twin_reduction.png``.

**Task 10** (``--task 10``): controls and sensitivity.

- **Controls.**  Each scene is composed three times: with both controls,
  without fluorescence, and without Raman.  Ed (and so every method's Ed_L)
  is the same for all three, so a method's view of a control is
  rho_m(all) - rho_m(without it), and the truth's likewise.  Per instrument,
  case (all 16) and method: the rms of (control'' - control''_true), the
  control's own rms'', and the total rho_w'' error, in the control's region:
  fluorescence over 665-705 nm (``FL_WINDOW``; it peaks in O2-B), Raman
  within +-5 nm of the ten lines (it fills them in).
- **Noise floor.**  rms[(rho_m(vii) - rho_m(iii))''] in the same region: the
  case (vii) noise alone, per method.
- **Degradation curves.**  Case (iii) with one perturbation swept
  (``SWEEP``): rigid E wavelength shift, stretch, the correction's FWHM_E
  or FWHM_L error, and the Emod's SZA offset or water-vapour factor; for
  linear, ruddick2023 and srf.

Outputs: ``$OS_COLOR/hypernet/wiggles/phase1/twin_controls.parquet`` and
``twin_degradation.parquet``; committed ``hypernet/wiggles/phase1_controls.csv``
and ``phase1_degradation.csv`` (medians over scenes); figures
``figs/phase1/twin_controls.png`` and ``twin_degradation.png``.

Run from the repository root:
``python -m hypernet.wiggles.phase1_twin [--task 9|10] [--workers 12] [--scenes N]``.
"""
import argparse
import dataclasses
import os
import time

import numpy as np
import pandas as pd

from hypernet import edinterp, emod, twin
from hypernet.wiggles import FIGDIR, WIGGLES_DIR

ROOT = os.path.join(os.getenv('OS_COLOR', '.'), 'hypernet', 'wiggles', 'phase1')
FIG1 = FIGDIR.replace('phase0', 'phase1')
INSTRUMENTS = ('HYPSTAR_122304', 'HYPSTAR_122305')
METHODS = ('linear', 'cubic', 'sinc', 'ruddick2023', 'srf')
HS = (1.0, 5.0)
#: The ten lines of plan §2.3 (nm, air).
LINES10 = {'Ca K': 393.37, 'Ca H': 396.85, 'G band': 430.79, 'H beta': 486.13, 'Mg b': 517.27,
           'Na D': 589.29, 'H alpha': 656.28, 'O2 B': 687.0, 'O2 A': 760.6, 'H2O': 936.5}
WIN = 5.0
L_RANGE = (392.0, 988.0)          # L pixels used (SRF fully inside the 380-1000 scene)


def second_diff(grid, y, h):
    """[y(l + h) - 2 y(l) + y(l - h)] / h^2 on ``grid`` (linear interpolation)."""
    return (np.interp(grid + h, grid, y) - 2 * y + np.interp(grid - h, grid, y)) / h ** 2


def masks(grid):
    near = np.zeros(grid.size, bool)
    per = {}
    for n, l0 in LINES10.items():
        m = np.abs(grid - l0) <= WIN
        per[n] = m
        near |= m
    band = (grid >= 400) & (grid <= 900)
    return near & band, ~near & band, per


def _shifted(srf_model, d):
    if not d:
        return srf_model
    return lambda l: srf_model.fwhm(l) + d


def model_ed(e, fields, i_case, variant, scene, emod_cache):
    """The correction's high-resolution Ed model (perfect, or a case vi variant)."""
    if not variant:
        return (scene['lam'], scene['Ed'])
    sza = float(fields['sza'][i_case])
    names = list(np.asarray(fields['name']).astype(str))
    j = i_case
    e2 = e
    if 'sza_offset' in variant:
        key = ('sza', sza + variant['sza_offset'])
        if key not in emod_cache:
            emod_cache[key] = emod.build_emod(sza + variant['sza_offset'])
        e2 = emod_cache[key]
    if 'pwv_factor' in variant:
        key = ('pwv', sza, variant['pwv_factor'])
        if key not in emod_cache:
            emod_cache[key] = emod.build_emod(sza, pwv_mm=15.0 * variant['pwv_factor'])
        e2 = emod_cache[key]
    if variant.get('aerosol') == 'low':
        nm = names[i_case]
        alt = nm.replace('aot0.25', 'aot0.05').replace('aot0.10', 'aot0.05')
        if nm.startswith('veit'):
            alt = 'sza50_aot0.05_chl1'               # nearest low-aerosol grid case
        j = names.index(alt) if alt in names else i_case
    s = twin._spline_log
    lam5 = fields['lam']
    F0 = e2['F0']
    ed = F0 * (s(lam5, fields['ed_dir'][j], e2['lam']) * e2['T_direct'] +
               s(lam5, fields['ed_dif'][j], e2['lam']) * e2['T_diffuse'])
    return (e2['lam'], ed)


def run_scene(args):
    """All instruments x cases x methods for one OSOAA scene."""
    i_case, seed = args
    fields = twin.load_fields()
    lib = twin.rhow_library(fields)
    grids = twin.load_grids()
    srfs = twin.load_instrument_srfs()
    sza = float(fields['sza'][i_case])
    e = emod.build_emod(sza)
    scene = twin.compose_scene(e, fields, i_case, lib=lib)
    lam = scene['lam']
    rows, emod_cache, keep = [], {}, {}
    for inst in INSTRUMENTS:
        cases = twin.case_table(inst, grids=grids, srfs=srfs)
        sL = srfs[inst]['Ld']
        gL_full = cases[0].grid_L
        mL = (gL_full >= L_RANGE[0]) & (gL_full <= L_RANGE[1])
        gL = gL_full[mL]
        truth = twin.rhow_true(scene, sL, gL)
        Lu0 = twin.observe(lam, scene['Lu'], gL, sL)
        Ld0 = twin.observe(lam, scene['Ld'], gL, sL)
        rho_eff = np.interp(gL, lam, scene['rho_eff'])
        near, away, per = masks(gL)
        t2 = {h: second_diff(gL, truth, h) for h in HS}
        nE = grids['noise_E_scan'] / np.sqrt(6)
        nLu = np.interp(gL, grids['grid_L_122304'], grids['noise_Lu_scan']) / np.sqrt(6)
        nLd = np.interp(gL, grids['grid_L_122304'], grids['noise_Ld_scan']) / np.sqrt(6)
        for k, c in enumerate(cases):
            gE = c.grid_E
            mE = (gE >= L_RANGE[0] - 3) & (gE <= L_RANGE[1] + 3)
            gE = gE[mE]
            tE = c.true_grid_E[mE] if c.true_grid_E is not None else None
            noisy = c.noise
            s0 = seed * 1000 + k
            Eobs = twin.observe(lam, scene['Ed'], gE, c.srf_E, true_grid=tE,
                                noise=(np.interp(gE, grids['grid_E_122304'], nE) if noisy
                                       else None), seed=s0)
            Lu = twin.observe(lam, scene['Lu'], gL, sL, noise=nLu, seed=s0 + 1) if noisy else Lu0
            Ld = twin.observe(lam, scene['Ld'], gL, sL, noise=nLd, seed=s0 + 2) if noisy else Ld0
            mdl = model_ed(e, fields, i_case, c.emod_variant, scene, emod_cache)
            sE_c = _shifted(c.srf_E, c.corr_dfwhm_E)
            sL_c = _shifted(c.srf_L, c.corr_dfwhm_L)
            for meth in METHODS:
                EdL = edinterp.interpolate_ed_to_l(gE, Eobs, gL, emod=mdl, srf_irr=sE_c,
                                                   srf_rad=sL_c, method=meth)
                rho = np.pi * (Lu - rho_eff * Ld) / EdL
                if inst == INSTRUMENTS[0] and c.name in ('iii_both', 'i_offset') and \
                        str(np.asarray(fields['name'])[i_case]) == 'veit_sza40_aot0.10_chl1':
                    keep[(inst, c.name, meth)] = (gL, rho, truth)
                for h in HS:
                    d = second_diff(gL, rho, h) - t2[h]
                    r = dict(instrument=inst, case=c.name, label=c.label,
                             scene=str(np.asarray(fields['name'])[i_case]), sza=sza,
                             aot550=float(fields['aot550'][i_case]),
                             water=str(np.asarray(fields['water'])[i_case]), method=meth, h=h,
                             rms_line=float(np.sqrt(np.mean(d[near] ** 2))),
                             rms_away=float(np.sqrt(np.mean(d[away] ** 2))),
                             rms_true_line=float(np.sqrt(np.mean(t2[h][near] ** 2))),
                             rms_rho_line=float(np.sqrt(np.mean(((rho - truth) / truth)[near]
                                                                ** 2))))
                    for n, m in per.items():
                        r['rms_' + n.replace(' ', '_')] = float(np.sqrt(np.mean(d[m] ** 2))) \
                            if m.any() else np.nan
                    rows.append(r)
    return rows, keep


def summarise(df):
    g = df.groupby(['instrument', 'case', 'label', 'method', 'h'])
    s = g.agg(rms_line=('rms_line', 'median'), rms_away=('rms_away', 'median'),
              rms_true_line=('rms_true_line', 'median'),
              rms_rho_line=('rms_rho_line', 'median'), n_scenes=('scene', 'nunique'),
              **{('rms_' + n.replace(' ', '_')): ('rms_' + n.replace(' ', '_'), 'median')
                 for n in LINES10}).reset_index()
    lin = df[df['method'] == 'linear'].set_index(['instrument', 'case', 'scene', 'h'])['rms_line']
    d = df.join(lin.rename('rms_line_linear'), on=['instrument', 'case', 'scene', 'h'])
    d['reduction'] = 1 - d['rms_line'] / d['rms_line_linear']
    red = d.groupby(['instrument', 'case', 'method', 'h'])['reduction'].median()
    s = s.join(red, on=['instrument', 'case', 'method', 'h'])
    s['line_over_away'] = s['rms_line'] / s['rms_away']
    return s


def figures(keep, summ):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    os.makedirs(FIG1, exist_ok=True)
    cols = {'linear': '#e34948', 'cubic': '#eda100', 'sinc': '#e87ba4',
            'ruddick2023': '#2a78d6', 'srf': '#1baf7a'}
    fig, axs = plt.subplots(2, 3, figsize=(13, 7), sharey=False)
    for row, cname in enumerate(('iii_both', 'i_offset')):
        for ax, (lab, lo, hi) in zip(axs[row], (('Ca H/K', 388, 402), ('Hα', 650, 663),
                                                ('O₂-A', 754, 772))):
            gL, _, truth = keep[(INSTRUMENTS[0], cname, 'linear')]
            m = (gL > lo) & (gL < hi)
            ax.plot(gL[m], second_diff(gL, truth, 1.0)[m], color='#333333', lw=2.2,
                    label='truth')
            for meth in ('linear', 'ruddick2023', 'srf'):
                g2, rho, _ = keep[(INSTRUMENTS[0], cname, meth)]
                ax.plot(g2[m], second_diff(g2, rho, 1.0)[m], color=cols[meth], lw=1.3,
                        label=meth)
            ax.set_title('%s, case %s' % (lab, cname), fontsize=9.5, loc='left')
            ax.grid(alpha=0.25, lw=0.6)
            ax.spines[['top', 'right']].set_visible(False)
    axs[0, 0].set_ylabel("ρw'' (h = 1 nm), nm⁻²")
    axs[1, 0].set_ylabel("ρw'' (h = 1 nm), nm⁻²")
    axs[0, 0].legend(fontsize=8, frameon=False)
    for ax in axs[1]:
        ax.set_xlabel('wavelength (nm)')
    fig.suptitle('Twin experiment, HYPSTAR 122304 SRFs, VEIT-like scene (SZA 40°, AOT 0.1, '
                 'Chl 1)', fontsize=10.5, x=0.06, ha='left')
    fig.tight_layout()
    fig.savefig(os.path.join(FIG1, 'twin_rho2_lines.png'), dpi=150)
    plt.close(fig)

    s = summ[summ['h'] == 1.0]
    order = list(dict.fromkeys(s.sort_values('case')['case']))
    fig, axs = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
    for ax, inst in zip(axs, INSTRUMENTS):
        g = s[s['instrument'] == inst].set_index(['case', 'method'])
        x = np.arange(len(order))
        for k, meth in enumerate(METHODS):
            y = [g.loc[(c, meth), 'rms_line'] if (c, meth) in g.index else np.nan for c in order]
            ax.semilogy(x + (k - 2) * 0.12, y, ls='', marker='o', ms=6, color=cols[meth],
                        label=meth)
        away = [g.loc[(c, 'srf'), 'rms_away'] for c in order]
        ax.semilogy(x, away, ls='', marker='_', ms=16, mew=2, color='#888',
                    label='srf, away from lines')
        ax.set_ylabel("median rms ρw'' error\nnear lines (h = 1 nm)")
        ax.set_title(inst.replace('HYPSTAR_', 'HYPSTAR ') +
                     (' (narrow-E)' if inst.endswith('304') else ' (E ≈ L)'), fontsize=10,
                     loc='left')
        ax.grid(alpha=0.25, lw=0.6)
        ax.spines[['top', 'right']].set_visible(False)
    h_, l_ = axs[0].get_legend_handles_labels()
    fig.legend(h_, l_, fontsize=8.5, frameon=False, ncol=6, loc='lower center')
    axs[1].set_xticks(np.arange(len(order)))
    axs[1].set_xticklabels(order, rotation=40, ha='right', fontsize=8)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(os.path.join(FIG1, 'twin_reduction.png'), dpi=150)
    plt.close(fig)


# --- Task 10: controls and degradation curves ------------------------------------

#: The fluorescence control's region (it peaks in O2-B); Raman is judged near
#: the ten lines (``masks``), where it fills them in.
FL_WINDOW = (665.0, 705.0)
CTL_REGION = {'fl': 'fl', 'raman': 'lines'}
CTL_VARIANTS = {'all': ('fl', 'raman'), 'no_fl': ('raman',), 'no_raman': ('fl',)}
#: Case (iii) with one perturbation swept; 0 (or factor 1) is case (iii).
SWEEP = {'shift': (-0.3, -0.2, -0.1, -0.05, -0.02, 0.0, 0.02, 0.05, 0.1, 0.2, 0.3),
         'stretch': (-0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2),
         'dfwhm_E': (-0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3),
         'dfwhm_L': (-0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3),
         'sza_offset': (-10.0, -5.0, -2.0, 0.0, 2.0, 5.0, 10.0),
         'pwv_factor': (0.5, 0.75, 1.0, 1.5, 2.0, 3.0)}
SWEEP_METHODS = ('linear', 'ruddick2023', 'srf')
YMIN_CTL = 1e-6                   # floor of the controls figure's ratio axis


def _rms(a):
    return float(np.sqrt(np.mean(a ** 2))) if a.size else np.nan


def sweep_cases(instrument, grids=None, srfs=None):
    """(kind, value, Case) for every ``SWEEP`` point: case (iii) with one
    perturbation (a stretch is +-value at 390 / 870 nm, as in case iv)."""
    base = next(c for c in twin.case_table(instrument, grids=grids, srfs=srfs)
                if c.name == 'iii_both')
    gE = base.grid_E
    out = []
    for kind, vals in SWEEP.items():
        for v in vals:
            if kind == 'shift':
                kw = dict(true_grid_E=gE + v)
            elif kind == 'stretch':
                kw = dict(true_grid_E=gE + v * (gE - 630.0) / 240.0)
            elif kind == 'dfwhm_E':
                kw = dict(corr_dfwhm_E=v)
            elif kind == 'dfwhm_L':
                kw = dict(corr_dfwhm_L=v)
            elif kind == 'sza_offset':
                kw = dict(emod_variant={'sza_offset': v} if v else {})
            else:
                kw = dict(emod_variant={'pwv_factor': v} if v != 1 else {})
            out.append((kind, v, dataclasses.replace(base, name='%s%+g' % (kind, v),
                                                     label='(iii) %s = %g' % (kind, v), **kw)))
    return out


def run_scene_t10(args):
    """Controls (16 cases x 5 methods) and the sweeps (x 3 methods), one scene."""
    i_case, seed = args
    fields = twin.load_fields()
    lib = twin.rhow_library(fields)
    grids = twin.load_grids()
    srfs = twin.load_instrument_srfs()
    sza = float(fields['sza'][i_case])
    e = emod.build_emod(sza)
    scenes = {v: twin.compose_scene(e, fields, i_case, controls=ctl, lib=lib)
              for v, ctl in CTL_VARIANTS.items()}
    lam, Ed = scenes['all']['lam'], scenes['all']['Ed']
    name = str(np.asarray(fields['name'])[i_case])
    meta = dict(scene=name, sza=sza, aot550=float(fields['aot550'][i_case]),
                water=str(np.asarray(fields['water'])[i_case]))
    ctl_rows, swp_rows, keep, emod_cache = [], [], {}, {}
    nE = grids['noise_E_scan'] / np.sqrt(6)
    for inst in INSTRUMENTS:
        cases = twin.case_table(inst, grids=grids, srfs=srfs)
        sL = srfs[inst]['Ld']
        gL = cases[0].grid_L
        gL = gL[(gL >= L_RANGE[0]) & (gL <= L_RANGE[1])]
        rho_eff = np.interp(gL, lam, scenes['all']['rho_eff'])
        near, away, _ = masks(gL)
        reg = {'lines': near, 'fl': (gL >= FL_WINDOW[0]) & (gL <= FL_WINDOW[1])}
        truth = {v: twin.rhow_true(s, sL, gL) for v, s in scenes.items()}
        Lu0 = {v: twin.observe(lam, s['Lu'], gL, sL) for v, s in scenes.items()}
        Ld0 = twin.observe(lam, scenes['all']['Ld'], gL, sL)
        nLu = np.interp(gL, grids['grid_L_122304'], grids['noise_Lu_scan']) / np.sqrt(6)
        nLd = np.interp(gL, grids['grid_L_122304'], grids['noise_Ld_scan']) / np.sqrt(6)

        def retrieve(c, s0, methods, variants):
            gE = c.grid_E
            mE = (gE >= L_RANGE[0] - 3) & (gE <= L_RANGE[1] + 3)
            gE = gE[mE]
            tE = c.true_grid_E[mE] if c.true_grid_E is not None else None
            Eobs = twin.observe(lam, Ed, gE, c.srf_E, true_grid=tE, seed=s0,
                                noise=(np.interp(gE, grids['grid_E_122304'], nE) if c.noise
                                       else None))
            if c.noise:      # the same seed, so the same noise draw for every variant
                Lu = {v: twin.observe(lam, scenes[v]['Lu'], gL, sL, noise=nLu, seed=s0 + 1)
                      for v in variants}
                Ld = twin.observe(lam, scenes['all']['Ld'], gL, sL, noise=nLd, seed=s0 + 2)
            else:
                Lu, Ld = Lu0, Ld0
            mdl = model_ed(e, fields, i_case, c.emod_variant, scenes['all'], emod_cache)
            sE_c, sL_c = _shifted(c.srf_E, c.corr_dfwhm_E), _shifted(c.srf_L, c.corr_dfwhm_L)
            out = {}
            for m in methods:
                EdL = edinterp.interpolate_ed_to_l(gE, Eobs, gL, emod=mdl, srf_irr=sE_c,
                                                   srf_rad=sL_c, method=m)
                out[m] = {v: np.pi * (Lu[v] - rho_eff * Ld) / EdL for v in variants}
            return out

        # controls
        t2 = {(v, h): second_diff(gL, truth[v], h) for v in CTL_VARIANTS for h in HS}
        rho_iii = None
        for k, c in enumerate(cases):
            R = retrieve(c, seed * 1000 + k, METHODS, CTL_VARIANTS)
            if c.name == 'iii_both':
                rho_iii = R
                if inst == INSTRUMENTS[0] and name == 'veit_sza40_aot0.10_chl1':
                    keep = dict(grid=gL, truth=truth, rho=R)
            for m in METHODS:
                for h in HS:
                    d_all = second_diff(gL, R[m]['all'], h) - t2[('all', h)]
                    nz = (second_diff(gL, R[m]['all'] - rho_iii[m]['all'], h)
                          if c.noise else None)
                    for ctl, rname in CTL_REGION.items():
                        off, msk = 'no_' + ctl, reg[rname]
                        cm = R[m]['all'] - R[m][off]
                        ct = truth['all'] - truth[off]
                        ct2 = t2[('all', h)] - t2[(off, h)]
                        ctl_rows.append(dict(
                            instrument=inst, case=c.name, label=c.label, **meta, method=m, h=h,
                            control=ctl, region=rname, ctl_amp=_rms(ct2[msk]),
                            ctl_err=_rms((second_diff(gL, cm, h) - ct2)[msk]),
                            ctl_level_err=_rms((cm - ct)[msk]) / _rms(ct[msk]),
                            tot_err=_rms(d_all[msk]),
                            noise=_rms(nz[msk]) if nz is not None else np.nan))
        # degradation sweeps
        for k, (kind, v, c) in enumerate(sweep_cases(inst, grids=grids, srfs=srfs)):
            R = retrieve(c, seed * 1000 + 100 + k, SWEEP_METHODS, ('all',))
            for m in SWEEP_METHODS:
                for h in HS:
                    d = second_diff(gL, R[m]['all'], h) - t2[('all', h)]
                    swp_rows.append(dict(instrument=inst, kind=kind, value=v, **meta, method=m,
                                         h=h, rms_line=_rms(d[near]), rms_away=_rms(d[away]),
                                         rms_true_line=_rms(t2[('all', h)][near])))
    return ctl_rows, swp_rows, keep


def summarise_controls(df):
    keys = ['instrument', 'case', 'label', 'control', 'region', 'method', 'h']
    s = df.groupby(keys).agg(ctl_amp=('ctl_amp', 'median'), ctl_err=('ctl_err', 'median'),
                             ctl_level_err=('ctl_level_err', 'median'),
                             tot_err=('tot_err', 'median'),
                             n_scenes=('scene', 'nunique')).reset_index()
    floor = (df[df['case'] == 'vii_noise'].groupby(['instrument', 'region', 'method', 'h'])
             ['noise'].median().rename('noise_floor'))
    s = s.join(floor, on=['instrument', 'region', 'method', 'h'])
    s['ctl_err_over_noise'] = s['ctl_err'] / s['noise_floor']
    s['ctl_err_over_amp'] = s['ctl_err'] / s['ctl_amp']
    s['ctl_amp_over_noise'] = s['ctl_amp'] / s['noise_floor']
    s['tot_err_over_noise'] = s['tot_err'] / s['noise_floor']
    return s


def summarise_degradation(df):
    keys = ['instrument', 'kind', 'value', 'method', 'h']
    s = df.groupby(keys).agg(rms_line=('rms_line', 'median'), rms_away=('rms_away', 'median'),
                             rms_true_line=('rms_true_line', 'median'),
                             n_scenes=('scene', 'nunique')).reset_index()
    lin = df[df['method'] == 'linear'].set_index(['instrument', 'kind', 'value', 'scene', 'h'])
    d = df.join(lin['rms_line'].rename('lin'), on=['instrument', 'kind', 'value', 'scene', 'h'])
    lin0 = df[(df['method'] == 'linear') & (df['kind'] == 'shift') & (df['value'] == 0)]
    d = d.join(lin0.set_index(['instrument', 'scene', 'h'])['rms_line'].rename('lin0'),
               on=['instrument', 'scene', 'h'])
    d['reduction'] = 1 - d['rms_line'] / d['lin']
    d['reduction_vs_iii_linear'] = 1 - d['rms_line'] / d['lin0']
    s = s.join(d.groupby(keys)[['reduction', 'reduction_vs_iii_linear']].median(), on=keys)
    return s


def figures_t10(keep, csum, dsum):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    cols = {'linear': '#e34948', 'cubic': '#eda100', 'sinc': '#e87ba4',
            'ruddick2023': '#2a78d6', 'srf': '#1baf7a'}
    fig, axs = plt.subplots(2, 2, figsize=(13, 8.5))
    gL, truth, R = keep['grid'], keep['truth'], keep['rho']
    for ax, (ctl, lo, hi, lab) in zip(axs[0], (('fl', 660, 712, 'fluorescence control'),
                                               ('raman', 750, 775, 'Raman control at O₂-A'))):
        m = (gL > lo) & (gL < hi)
        off = 'no_' + ctl
        ax.plot(gL[m], second_diff(gL, truth['all'] - truth[off], 1.0)[m], color='#333333',
                lw=2.4, label='truth')
        for meth in ('linear', 'ruddick2023', 'srf'):
            ax.plot(gL[m], second_diff(gL, R[meth]['all'] - R[meth][off], 1.0)[m],
                    color=cols[meth], lw=1.2, label=meth)
        ax.set_title("%s: Δρw'' (h = 1 nm), 122304, case (iii), VEIT-like scene" % lab,
                     fontsize=9.5, loc='left')
        ax.set_xlabel('wavelength (nm)')
        ax.set_ylabel("Δρw'', nm⁻²")
        ax.legend(fontsize=8, frameon=False)
    s = csum[csum['h'] == 1.0]
    order = list(dict.fromkeys(s.sort_values('case')['case']))
    x = np.arange(len(order))
    for ax, inst in zip(axs[1], INSTRUMENTS):
        g = s[s['instrument'] == inst].set_index(['case', 'control', 'method'])
        for k, meth in enumerate(METHODS):
            for ctl, mk in (('fl', 'o'), ('raman', 's')):
                y = [max(g.loc[(c, ctl, meth), 'ctl_err_over_noise'], YMIN_CTL) for c in order]
                ax.semilogy(x + (k - 2) * 0.13, y, ls='', marker=mk, ms=5.5, color=cols[meth],
                            mfc=cols[meth] if ctl == 'fl' else 'none',
                            label='%s, %s' % (meth, 'fluorescence' if ctl == 'fl' else 'Raman'))
        ax.axhline(1.0, color='#888', lw=1, ls='--')
        ax.set_ylim(YMIN_CTL / 2, 3)
        ax.text(0.01, 0.93, 'points on the bottom edge: < %.0e (exact to rounding)' % YMIN_CTL,
                transform=ax.transAxes, fontsize=7.5, color='#555')
        ax.set_xticks(x)
        ax.set_xticklabels(order, rotation=40, ha='right', fontsize=7.5)
        ax.set_title('%s: control error / case (vii) noise floor (h = 1 nm)'
                     % inst.replace('HYPSTAR_', 'HYPSTAR '), fontsize=9.5, loc='left')
    axs[1, 0].set_ylabel("rms Δρw'' error / noise")
    for ax in axs.ravel():
        ax.grid(alpha=0.25, lw=0.6)
        ax.spines[['top', 'right']].set_visible(False)
    h_, l_ = axs[1, 0].get_legend_handles_labels()
    fig.legend(h_, l_, fontsize=7.5, frameon=False, ncol=5, loc='lower center')
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(os.path.join(FIG1, 'twin_controls.png'), dpi=150)
    plt.close(fig)

    s = dsum[dsum['h'] == 1.0]
    floor = csum[(csum['h'] == 1.0) & (csum['case'] == 'vii_noise') & (csum['method'] == 'srf')
                 & (csum['control'] == 'raman')].set_index('instrument')['noise_floor']
    panels = (('shift', 'E wavelength shift (nm)'), ('stretch', 'E stretch at the ends (nm)'),
              ('dfwhm', 'correction FWHM error (nm)\nsolid: FWHM_E, dashed: FWHM_L'), ('sza_offset', 'Emod SZA offset (deg)'),
              ('pwv_factor', 'Emod water-vapour factor'))
    fig, axs = plt.subplots(2, 5, figsize=(17, 7.5), sharey='row')
    for row, inst in enumerate(INSTRUMENTS):
        g = s[s['instrument'] == inst]
        for ax, (kind, xl) in zip(axs[row], panels):
            for meth in SWEEP_METHODS:
                for kk, ls in ((('dfwhm_E', '-'), ('dfwhm_L', '--')) if kind == 'dfwhm'
                               else ((kind, '-'),)):
                    q = g[(g['kind'] == kk) & (g['method'] == meth)].sort_values('value')
                    lab = meth if kind != 'dfwhm' else '%s, %s' % (meth, kk[-1])
                    ax.semilogy(q['value'], q['rms_line'], ls=ls, marker='o', ms=4, lw=1.6,
                                color=cols[meth], label=lab)
            ax.axhline(g['rms_true_line'].median(), color='#333', lw=1, ls=':',
                       label="true line structure")
            ax.axhline(floor[inst], color='#888', lw=1.2, ls='-.', label='(vii) noise floor')
            ax.set_xlabel(xl)
            ax.grid(alpha=0.25, lw=0.6)
            ax.spines[['top', 'right']].set_visible(False)
        axs[row, 0].set_ylabel("%s\nmedian rms ρw'' error\nnear lines (h = 1 nm)"
                               % inst.replace('HYPSTAR_', 'HYPSTAR '))
    h_, l_ = axs[0, 0].get_legend_handles_labels()
    fig.legend(h_, l_, fontsize=8, frameon=False, ncol=5, loc='lower center')

    fig.suptitle('Case (iii) degradation curves (medians over 25 scenes)', fontsize=10.5,
                 x=0.05, ha='left')
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(os.path.join(FIG1, 'twin_degradation.png'), dpi=150)
    plt.close(fig)


def main_t10(n, workers):
    fields = twin.load_fields()
    # build every Emod once (cached), so workers only read them
    for sza in sorted(set(float(x) for x in fields['sza'][:n])):
        emod.build_emod(sza)
        for v in SWEEP['sza_offset']:
            emod.build_emod(sza + v)
        for v in set(SWEEP['pwv_factor']) | {2.0}:
            emod.build_emod(sza, pwv_mm=15.0 * v)
    t0 = time.time()
    from multiprocessing import Pool
    crows, srows, keep = [], [], {}
    with Pool(workers) as pool:
        for c, s_, k in pool.imap_unordered(run_scene_t10, [(i, i) for i in range(n)]):
            crows += c
            srows += s_
            keep = keep or k
    cdf, sdf = pd.DataFrame(crows), pd.DataFrame(srows)
    os.makedirs(ROOT, exist_ok=True)
    cdf.to_parquet(os.path.join(ROOT, 'twin_controls.parquet'))
    sdf.to_parquet(os.path.join(ROOT, 'twin_degradation.parquet'))
    csum, dsum = summarise_controls(cdf), summarise_degradation(sdf)
    csum.to_csv(os.path.join(WIGGLES_DIR, 'phase1_controls.csv'), index=False,
                float_format='%.4g')
    dsum.to_csv(os.path.join(WIGGLES_DIR, 'phase1_degradation.csv'), index=False,
                float_format='%.4g')
    if keep:
        figures_t10(keep, csum, dsum)
    pd.set_option('display.width', 220)
    print('%d control rows, %d sweep rows in %.0f s' % (len(cdf), len(sdf), time.time() - t0))
    s = csum[csum['h'] == 1.0]
    for col in ('ctl_err_over_noise', 'tot_err_over_noise'):
        print('\n%s (h = 1 nm):' % col)
        print(s.pivot_table(index=['instrument', 'control', 'case'], columns='method',
                            values=col)[list(METHODS)].to_string(float_format='%.3g'))
    print('\ncontrol amplitude / noise floor (srf, case iii):')
    print(csum[(csum['method'] == 'srf') & (csum['case'] == 'iii_both')].pivot_table(
        index=['instrument', 'control'], columns='h', values='ctl_amp_over_noise')
        .to_string(float_format='%.3g'))
    print('\nworst case over (i)-(vii): control error / noise floor and / control amplitude')
    w = csum.groupby(['instrument', 'control', 'h', 'method'])[
        ['ctl_err_over_noise', 'ctl_err_over_amp']].max().unstack('method')
    print(w.to_string(float_format='%.3g'))
    print('\ncase iii: control error / control amplitude')
    print(csum[csum['case'] == 'iii_both'].pivot_table(
        index=['instrument', 'control', 'h'], columns='method', values='ctl_err_over_amp')
        [list(METHODS)].to_string(float_format='%.3g'))
    print('\ncontrol level error (rms / rms control), case iii, h = 1:')
    print(s[s['case'] == 'iii_both'].pivot_table(index=['instrument', 'control'],
                                                 columns='method', values='ctl_level_err')
          [list(METHODS)].to_string(float_format='%.2e'))
    d = dsum[dsum['h'] == 1.0]
    print('\ndegradation: rms_line (h = 1 nm)')
    print(d.pivot_table(index=['instrument', 'kind', 'value'], columns='method',
                        values='rms_line')[list(SWEEP_METHODS)].to_string(float_format='%.3g'))
    print('\ndegradation: srf reduction vs linear (same perturbation) / vs case iii linear')
    print(d[d['method'] == 'srf'].pivot_table(
        index=['kind', 'value'], columns='instrument',
        values=['reduction', 'reduction_vs_iii_linear']).to_string(float_format='%.2f'))


def main(argv=None):
    ap = argparse.ArgumentParser(description='Phase 1 tasks 9-10: twin experiment')
    ap.add_argument('--task', type=int, default=9, choices=(9, 10))
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--scenes', type=int, default=None, help='first N scenes only')
    a = ap.parse_args(argv)
    fields = twin.load_fields()
    n = len(fields['name']) if a.scenes is None else a.scenes
    if a.task == 10:
        return main_t10(n, a.workers)
    # build the Emods once (cached), so workers only read them
    for sza in sorted(set(float(x) for x in fields['sza'][:n])):
        emod.build_emod(sza)
        emod.build_emod(sza + 10.0)
        emod.build_emod(sza, pwv_mm=30.0)
    t0 = time.time()
    from multiprocessing import Pool
    rows, keep = [], {}
    with Pool(a.workers) as pool:
        for r, k in pool.imap_unordered(run_scene, [(i, i) for i in range(n)]):
            rows += r
            keep.update(k)
    df = pd.DataFrame(rows)
    os.makedirs(ROOT, exist_ok=True)
    df.to_parquet(os.path.join(ROOT, 'twin_metrics.parquet'))
    summ = summarise(df)
    summ.to_csv(os.path.join(WIGGLES_DIR, 'phase1_twin_metrics.csv'), index=False,
                float_format='%.4g')
    if keep:
        figures(keep, summ)
    pd.set_option('display.width', 220)
    print('%d rows in %.0f s' % (len(df), time.time() - t0))
    s = summ[summ['h'] == 1.0]
    print(s.pivot_table(index=['instrument', 'case'], columns='method', values='reduction')
          [list(METHODS)].to_string(float_format='%.2f'))
    print('\nline/away ratio of the rms error (h = 1 nm), srf vs linear:')
    print(s[s['method'].isin(['linear', 'srf'])].pivot_table(
        index=['instrument', 'case'], columns='method', values='line_over_away')
        .to_string(float_format='%.2f'))


if __name__ == '__main__':
    main()

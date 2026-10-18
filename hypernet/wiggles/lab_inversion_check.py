"""Validate the L2 inversion (lab data prompts, task 5; no Phase 3).

L1C/L2 ``irradiance`` is E resampled onto the L grid, E_L = M E
(``edinterp.interp_matrix``).  ``edinterp.recover_native_e`` inverts it with
the E grid the processor used (``labcal.grid_for``).  With the right grid the
recovered E is the native E; with a wrong one, the least-squares solution
picks up pixel-alternating errors.  This script measures that on VEIT (where
the L1A E is the truth) and on Release 2.

**Result: the L2 data cannot tell which E grid the processor used.**
Solving M_wrong x = M_true E makes the wrong grid's piecewise-linear curve
match the true one at the L pixels, so x is the *same spectrum resampled at
the wrong grid's own positions*, correctly labelled.  A grid error of 0.05 nm
therefore produces neither a wavelength-scale error nor ringing; it only
changes the double-interpolation smoothing.  On VEIT, against the L1A truth
resampled at the wrong grid's positions, a +-0.05 nm error gives 7e-4 to 1e-3
in E and 24-36 % in the second differences, and the 2024_09 grid 6.6e-4 and
17 % (with the right grid: 1.5e-5 and 0.2 %; ``veit.csv``).  Phase 3 would
feel that, but nothing in the recovered spectrum itself shows it.  Three candidate detectors were tried:

- ``edinterp.second_difference_power`` (plain and ``nyquist=True``, the
  "ringing index" ``ring``): it changes by 1-4 % on VEIT and < 1 % on
  synthetic spectra (``lab_inversion_synthetic.py``);
- the minimum of ring(delta) over the grid shifted by delta (``delta_min``):
  on VEIT it sits at the truth, but on synthetic spectra with a *correct*
  grid it sits at -0.06 to -0.08 nm, and ``own+0.05`` moves it by -0.05
  only because ring depends on the absolute grid position;
- native minus via_grid wavelength scale (``dshift``, ``dstretch``,
  ``srfcal.fit_wavelength_scale``): +-0.02 nm for every grid, including the
  0.05 nm control, as the resampling argument predicts.

All three are kept as diagnostics.  What the ringing index *does* catch is a
grossly wrong grid with a different pixel phase pattern: the 122304 proxy
raises it 1.4-2.6x on 5 of 6 instruments (not on 122303, 0.83x) (``ring_proxy_over_own``), while
every own grid rings at VEIT's level.  The verdict on the E epoch
(:func:`verdict`) is therefore ``verified`` for 122304 (the VEIT L1A E grid),
and elsewhere at most ``consistent (gross); <=0.05 nm not testable``: the
evidence for the epoch stays the L-grid match (``inferred_L``).

Steps:

1. **VEIT** (122304, L1A E = truth): recovery error with ``grid_for``, the L1A
   grid, the 2024_09 grid and grids shifted by +-0.05 nm; the metrics' response.
2. **Release 2**: up to 5 sequences per instrument x calibration period of
   ``docs/wiggles_data_request.csv``.  Recover with the own grid; check
   positivity, smoothness and ringing; scan delta; repeat with wrong-grid
   controls (another lab epoch where one exists; 122304's E grid as the old
   proxy; the own grid + 0.05 nm, a positive control of the scan).  Verdict per
   instrument x period (:func:`verdict`).
3. **Widths**: HSRS template fits on the recovered native-grid E (Gaussian;
   super-Gaussian with the lab p where the SRF table assigns one), against
   the ``via_grid`` route of task 3 (the same sequences' E through
   ``fit_template_windows(via_grid=...)``) and the lab FWHM.

Outputs:

- ``hypernet/wiggles/lab_inversion_check.csv`` (per sequence x grid);
- ``hypernet/wiggles/lab_inversion_summary.csv`` (per instrument x period:
  the E-epoch verdict, the control, the widths);
- ``hypernet/wiggles/figs/lab/inversion_check.png``;
- intermediates in ``$OS_COLOR/hypernet/wiggles/lab_data/inversion/``.

Run from the repository root:
``python -m hypernet.wiggles.lab_inversion_check [--workers 12] [--per-group 5]``
"""
import argparse
import os
import time

import numpy as np
import pandas as pd

from hypernet import edinterp, labcal, srf, srfcal, twin, whn_srf
from hypernet import whn_l1a as wl
from hypernet.wiggles import REPO, WIGGLES_DIR

REQUEST = os.path.join(REPO, 'docs', 'wiggles_data_request.csv')
OUT = os.path.join(labcal.lab_root(), 'inversion')
FIG = os.path.join(WIGGLES_DIR, 'figs', 'lab')
VEIT = dict(site='VEIT', seq_time='20260604T0845')
BAND = (400.0, 900.0)
DELTAS = np.round(np.arange(-0.10, 0.1001, 0.01), 3)


def ringing(E, wav, band=BAND):
    """The ringing index, ``edinterp.second_difference_power(nyquist=True)``."""
    return float(edinterp.second_difference_power(E, wav, band, nyquist=True))


def scan(EL, wL, wE):
    """ring(delta) for the grid wE shifted by DELTAS."""
    return np.array([ringing(edinterp.recover_native_e(EL, wL, wE + d), wE + d) for d in DELTAS])


def veit():
    f = wl.sequence_files(**VEIT)
    irr = wl.load_l1a_irr(f['L1A_IRR'])
    l1c = wl.load_l1c(f['L1C_ALL'])
    E, wT = irr['mean'], irr['wave']
    wL = l1c['wave']
    EL = np.nanmean(l1c['irradiance'], axis=1)
    m = (wT >= BAND[0]) & (wT <= BAND[1])
    vmodel = srf.load_srf_models(os.path.join(os.path.dirname(WIGGLES_DIR), 'data',
                                              'veit_srf_model.json'))['E']
    rows = []
    # the expected grid error of each test grid relative to the truth, as a
    # shift/stretch fitted to (grid - L1A grid) in task 4b's convention
    x = (wT - srfcal.LAM0) / srfcal.HALF
    grids = {'grid_for (2023_01)': labcal.grid_for('122304', 'E'), 'L1A grid': wT,
             '2024_09 grid': labcal.grid_for('122304', 'E', epoch='2024_09'),
             'grid_for +0.05 nm': labcal.grid_for('122304', 'E') + 0.05,
             'grid_for -0.05 nm': labcal.grid_for('122304', 'E') - 0.05}
    for label, g in grids.items():
        if g.size != wT.size:
            rows.append(dict(grid=label, note='%d px, not comparable pixel by pixel' % g.size))
            continue
        Eh, info = edinterp.recover_native_e(EL, wL, g, return_info=True)
        cg = np.polyfit(x, g - wT, 1)
        rows_expect = dict(grid_err_shift=float(cg[1]), grid_err_stretch=float(cg[0]))
        rel = Eh / E - 1
        from scipy.interpolate import CubicSpline
        Er = CubicSpline(wT, E)(g)                     # the truth resampled at g's positions
        rel_r = Eh / Er - 1
        d2h = (Eh[2:] - 2 * Eh[1:-1] + Eh[:-2])
        d2t = (E[2:] - 2 * E[1:-1] + E[:-2])
        d2r = (Er[2:] - 2 * Er[1:-1] + Er[:-2])
        md = m[1:-1]
        sc = scan(EL, wL, g)
        nv = native_minus_via(wL, EL, g, Eh, vmodel)
        rows.append(dict(grid=label, dshift=nv['dshift'], dstretch=nv['dstretch'],
                         dshift_err=nv['dshift_err'],
                         rel_rms=float(np.sqrt(np.mean(rel[m] ** 2))),
                         rel_rms_vs_resampled=float(np.sqrt(np.mean(rel_r[m] ** 2))),
                         d2_err_frac_vs_resampled=float(
                             np.sqrt(np.mean((d2h - d2r)[m[1:-1]] ** 2)) /
                             np.sqrt(np.mean(d2r[m[1:-1]] ** 2))),
                         rel_max=float(np.abs(rel[m]).max()),
                         d2_err_frac=float(np.sqrt(np.mean((d2h - d2t)[md] ** 2)) /
                                           np.sqrt(np.mean(d2t[md] ** 2))),
                         sdp=float(edinterp.second_difference_power(Eh, g)),
                         sdp_truth=float(edinterp.second_difference_power(E, wT)),
                         ring=ringing(Eh, g), ring_truth=ringing(E, wT),
                         ring_min=float(sc.min()), delta_min=float(DELTAS[np.argmin(sc)]),
                         n_unconstrained=info['n_unconstrained'],
                         forward_rel_rms=info['forward_rel_rms'], **rows_expect))
    return pd.DataFrame(rows)


#: wrong-grid controls per instrument: other lab epochs (plain crop rule)
OTHER_EPOCH = {'120242': '2022_06', '121222': '2021_07_tst', '122302': '2024_09',
               '122304': '2024_09'}
LINES = (435.8328, 546.0735)


def _width_models(w, spec, ref, via, sh):
    """Gaussian and (if assigned) lab-p super-Gaussian template E models."""
    out = {}
    t = srf.fit_template_windows(w, spec, ref[0], ref[1], veil=False, via_grid=via)
    out['gauss'] = srf.SRFModel.from_lines('E', t, scale_cov=True)
    if sh['shape'] == 'supergauss':
        t = srf.fit_template_windows(w, spec, ref[0], ref[1], veil=False, via_grid=via,
                                     shape='supergauss', shape_par=sh['shape_p'],
                                     fix_shape_par=True)
        out['supergauss'] = srf.SRFModel.from_lines('E', t, scale_cov=True)
    return out


def native_minus_via(wL, EL, g, Eh, model, ref=None, err=None):
    """Wavelength-scale fits of the recovered native E (on g) and of the L2
    irradiance through g (via_grid); their difference measures the grid error."""
    ok = np.isfinite(Eh)
    a = srfcal.fit_wavelength_scale(g[ok], Eh[ok], srf=model, ref=ref,
                                    err=None if err is None else err[ok])
    b = srfcal.fit_wavelength_scale(wL, EL, srf=model, ref=ref, via_grid=g)
    return dict(shift_native=a['shift'], stretch_native=a['stretch'],
                shift_err_native=a['shift_err'], stretch_err_native=a['stretch_err'],
                shift_via=b['shift'], stretch_via=b['stretch'],
                dshift=a['shift'] - b['shift'], dstretch=a['stretch'] - b['stretch'],
                dshift_err=float(np.hypot(a['shift_err'], b['shift_err'])))


def one_sequence(args):
    r, ref = args
    t0 = time.time()
    try:
        d = wl.load_l2b(wl.release2_path(r['site'], r['sequence_time'], r['file']))
        m = d['meta']
        sysid = m['system_id']
        sn = sysid[-6:]
        irr = m['instrument_calibration_date_irr']
        cal = '%s/%s' % (m['instrument_calibration_date_rad'], irr)
        wL, EL = d['wave'], d['irradiance'][:, 0]
        grids = {'own': labcal.grid_for(sn, 'E'), 'own+0.05': labcal.grid_for(sn, 'E') + 0.05}
        if sn != '122304':
            grids['proxy_122304'] = labcal.grid_for('122304', 'E')
        if sn in OTHER_EPOCH:
            grids['epoch_' + OTHER_EPOCH[sn]] = labcal.grid_for(sn, 'E', epoch=OTHER_EPOCH[sn])
        keys = dict(instrument=sysid, cal_period=cal, site=r['site'], sequence_time=r['sequence_time'],
                    sza=float(np.nanmean(d['sza'])))
        srfs = twin.load_instrument_srfs()
        emod_srf = srfs.get('%s@%s' % (sysid, irr), srfs.get(sysid))['E']
        rows, curves = [], []
        for label, g in grids.items():
            Eh, info = edinterp.recover_native_e(EL, wL, g, return_info=True)
            band = (g >= BAND[0]) & (g <= BAND[1])
            sc = scan(EL, wL, g)
            d2 = (Eh[2:] - 2 * Eh[1:-1] + Eh[:-2]) / Eh[1:-1]
            alt = np.sign(d2[band[1:-1]])
            rows.append(dict(keys, grid=label, n_E=g.size,
                             n_unconstrained=info['n_unconstrained'],
                             forward_rel_rms=info['forward_rel_rms'],
                             min_E_band=float(np.nanmin(Eh[band])), positive=bool(np.nanmin(Eh[band]) > 0),
                             sign_flip_frac=float(np.mean(alt[1:] != alt[:-1])),
                             sdp=float(edinterp.second_difference_power(Eh, g)),
                             ring=ringing(Eh, g), ring_min=float(sc.min()),
                             delta_min=float(DELTAS[np.argmin(sc)])))
            sn_ = native_minus_via(wL, EL, g, Eh, emod_srf)
            rows[-1].update(sn_)
            curves.append(pd.DataFrame(dict(keys, grid=label, delta=DELTAS, ring=sc)))
        # widths: native-grid fit of the recovered E vs the via_grid route of task 3
        sh = labcal.srf_shape_for(sn, 'E', irr)
        g = grids['own']
        Eh = edinterp.recover_native_e(EL, wL, g)
        ok = np.isfinite(Eh)
        nat = _width_models(g[ok], Eh[ok], ref, None, sh)
        via = _width_models(wL, EL, ref, g, sh)
        ep, _ = labcal.physical_epoch(sn, irr)
        lab = labcal.load_lab_srf()
        lab = lab[(lab['sn'] == sn) & (lab['epoch'] == ep) & (lab['channel'] == 'E')]
        wrows = []
        for shape in nat:
            for lam in LINES:
                l = lab[np.isclose(lab['lam_ref'], lam)].iloc[0]
                wrows.append(dict(keys, shape=shape, lam=lam, fwhm_native=float(nat[shape].fwhm(lam)),
                                  fwhm_via=float(via[shape].fwhm(lam)),
                                  fwhm_lab=float(l['fwhm_sg'] if shape == 'supergauss' else l['lab_fwhm']),
                                  lab_epoch=ep, p=sh['shape_p'] if shape == 'supergauss' else np.nan))
        return rows, curves, wrows, round(time.time() - t0, 1)
    except Exception as ex:                       # noqa: BLE001 -- logged, not fatal
        return [dict(site=r['site'], sequence_time=r['sequence_time'],
                     error='%s: %s' % (type(ex).__name__, ex))], [], [], 0.0


def release2(ref, workers, per_group):
    from multiprocessing import Pool
    req = pd.read_csv(REQUEST, dtype=str)
    req = req[[os.path.exists(wl.release2_path(a, b, c))
               for a, b, c in zip(req['site'], req['sequence_time'], req['file'])]]
    req = req.groupby(['instrument', 'cal_dates_rad_irr'], group_keys=False).head(per_group)
    rows, curves, widths = [], [], []
    with Pool(workers) as pool:
        for r_, c_, w_, _ in pool.imap_unordered(one_sequence, [(r, ref) for _, r in req.iterrows()]):
            rows += r_
            curves += c_
            widths += w_
    return (pd.DataFrame(rows), pd.concat(curves, ignore_index=True) if curves else pd.DataFrame(),
            pd.DataFrame(widths))


#: VEIT's ringing index with the true grid (L1C, ``veit.csv``) and the factor
#: above it that a grossly wrong grid (the 122304 proxy) reaches on Release 2.
RING_VEIT = 0.0077
RING_GROSS = 1.3


def verdict(sn, d):
    """``verified`` where an L1A E grid confirms the epoch (122304).  Else
    ``consistent (gross); <=0.05 nm not testable`` if the own grid's ringing
    is at VEIT's level (< RING_GROSS x RING_VEIT), and ``suspect`` if it is
    above; ``check`` if the recovered E is not positive in 400-900 nm.  A
    0.05 nm offset of the grid is invisible from L2 (module docstring)."""
    if not d['positive_all']:
        return 'check'
    if labcal.PROCESSOR_EPOCH[sn]['E'] == 'verified':
        return 'verified'
    if d['ring_med'] > RING_GROSS * RING_VEIT:
        return 'suspect'
    return 'consistent (gross); <=0.05 nm not testable'


def summarise(t, wd):
    """Per instrument x period x grid: the ringing-minimum offset delta_min
    (median and range over the sequences) and the excess ringing of the grid
    over its minimum; the E-epoch verdict on the own grid.  Then the widths."""
    t = t.dropna(subset=['grid'])
    out = []
    for (inst, cal, grid), h in t.groupby(['instrument', 'cal_period', 'grid']):
        d = dict(instrument=inst, cal_period=cal, grid=grid, n_seq=h['sequence_time'].nunique(),
                 delta_min_med=float(h['delta_min'].median()),
                 delta_min_lo=float(h['delta_min'].min()), delta_min_hi=float(h['delta_min'].max()),
                 excess_med=float((h['ring'] / h['ring_min']).median() - 1),
                 ring_med=float(h['ring'].median()), sdp_med=float(h['sdp'].median()),
                 positive_all=bool(h['positive'].all()),
                 sign_flip_frac_med=float(h['sign_flip_frac'].median()),
                 n_unconstrained=int(h['n_unconstrained'].max()),
                 dshift_med=float(h['dshift'].median()), dstretch_med=float(h['dstretch'].median()),
                 dshift_err_med=float(h['dshift_err'].median()))
        if grid == 'own':
            d['verdict'] = verdict(inst[-6:], d)
            pr = t[(t['instrument'] == inst) & (t['cal_period'] == cal) & (t['grid'] == 'proxy_122304')]
            d['ring_proxy_over_own'] = (float(pr['ring'].median()) / d['ring_med']
                                        if len(pr) else np.nan)
        w = wd[(wd['instrument'] == inst) & (wd['cal_period'] == cal)] if grid == 'own' else wd.iloc[0:0]
        for shape, hh in w.groupby('shape'):
            for c in ('fwhm_native', 'fwhm_via', 'fwhm_lab'):
                d['%s_%s' % (c, shape)] = float(hh[c].median())
        out.append(d)
    return pd.DataFrame(out)


def figure(v, curves, summ, wd):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    fig, axs = plt.subplots(1, 3, figsize=(17, 5))
    a = axs[0]
    f = wl.sequence_files(**VEIT)
    irr, l1c = wl.load_l1a_irr(f['L1A_IRR']), wl.load_l1c(f['L1C_ALL'])
    EL = np.nanmean(l1c['irradiance'], axis=1)
    wT, E = irr['wave'], irr['mean']
    zoom = (wT > 425) & (wT < 445)
    for label, g, col in (('grid_for (true)', labcal.grid_for('122304', 'E'), '#1baf7a'),
                          ('2024_09 grid', labcal.grid_for('122304', 'E', epoch='2024_09'), '#e34948'),
                          ('true + 0.05 nm', labcal.grid_for('122304', 'E') + 0.05, '#eda100')):
        Eh = edinterp.recover_native_e(EL, l1c['wave'], g)
        a.plot(wT[zoom], 100 * (Eh / E - 1)[zoom], '-o', ms=2.5, lw=1, color=col,
               label='%s: rms %.1e (400-900 nm)' % (label, np.sqrt(np.mean(
                   (Eh / E - 1)[(wT > 400) & (wT < 900)] ** 2))))
    a.set_xlabel('wavelength (nm)')
    a.set_ylabel('recovered / L1A E − 1 (%)')
    a.set_title('(a) VEIT: a wrong E grid is visible only against the L1A truth', fontsize=10,
                loc='left')
    a.legend(fontsize=7.5, frameon=False)
    a = axs[1]
    insts = sorted(curves['instrument'].unique())
    cols = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7']
    for k, inst in enumerate(insts):
        for grid, ls in (('own', '-'), ('proxy_122304', ':')):
            c = curves[(curves['instrument'] == inst) & (curves['grid'] == grid)]
            if c.empty:
                continue
            q = c.groupby(['site', 'sequence_time', 'cal_period'])
            norm = c['ring'] / q['ring'].transform('min')
            m = c.assign(n=norm).groupby('delta')['n'].median()
            a.plot(m.index, m.values, ls, color=cols[k % 7], lw=1.5,
                   label=inst[-6:] + (' own' if grid == 'own' else ' 122304 proxy'))
    a.set_xlabel('shift δ (nm)')
    a.set_ylabel('ringing / minimum (median over sequences)')
    a.set_title('(b) Release 2 ringing scans (diagnostic, not a test): own (solid), proxy (dotted)',
                fontsize=9, loc='left')
    a.legend(fontsize=6.5, frameon=False, ncol=2)
    a = axs[2]
    w = wd[np.isclose(wd['lam'], 546.0735)]
    x = 0
    ticks = []
    for (inst, shape), h in w.groupby(['instrument', 'shape']):
        for c, col, mk in (('fwhm_native', '#2a78d6', 'o'), ('fwhm_via', '#eb6834', 's'),
                           ('fwhm_lab', '#333333', 'D')):
            a.plot(x + {'fwhm_native': -0.2, 'fwhm_via': 0, 'fwhm_lab': 0.2}[c], h[c].median(), mk,
                   color=col, ms=6, label={'fwhm_native': 'recovered native E', 'fwhm_via': 'via_grid',
                                           'fwhm_lab': 'lab'}[c] if x == 0 else None)
        ticks.append('%s\n%s' % (inst[-6:], 'sg' if shape == 'supergauss' else 'g'))
        x += 1
    a.set_xticks(range(len(ticks)))
    a.set_xticklabels(ticks, fontsize=7)
    a.set_ylabel('E FWHM at 546 nm (nm)')
    a.set_title('(c) E widths: native vs via_grid vs lab (g = Gaussian, sg = lab p)', fontsize=10,
                loc='left')
    a.legend(fontsize=8, frameon=False)
    for a in axs:
        a.grid(alpha=0.25, lw=0.5)
        a.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    p = os.path.join(FIG, 'inversion_check.png')
    fig.savefig(p, dpi=140)
    plt.close(fig)
    return p


def main(argv=None):
    ap = argparse.ArgumentParser(description='lab data task 5: the L2 inversion')
    ap.add_argument('--veit-only', action='store_true')
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--per-group', type=int, default=5)
    a = ap.parse_args(argv)
    os.makedirs(OUT, exist_ok=True)
    pd.set_option('display.width', 260)
    pd.set_option('display.max_columns', 60)
    v = veit()
    v.to_csv(os.path.join(OUT, 'veit.csv'), index=False, float_format='%.6g')
    print(v.to_string(index=False, float_format='%.4g'))
    if a.veit_only:
        return
    t0 = time.time()
    ref = whn_srf.load_reference('air')
    t, curves, wd = release2(ref, a.workers, a.per_group)
    t.to_csv(os.path.join(WIGGLES_DIR, 'lab_inversion_check.csv'), index=False, float_format='%.5g')
    curves.to_csv(os.path.join(OUT, 'ring_scans.csv'), index=False, float_format='%.6g')
    wd.to_csv(os.path.join(OUT, 'widths.csv'), index=False, float_format='%.4f')
    summ = summarise(t, wd)
    summ.to_csv(os.path.join(WIGGLES_DIR, 'lab_inversion_summary.csv'), index=False,
                float_format='%.4f')
    if 'error' in t:
        print('errors:', t['error'].dropna().tolist())
    print('\n%d sequences in %.0f s' % (t['sequence_time'].nunique(), time.time() - t0))
    print(summ.to_string(index=False, float_format='%.3f'))
    print('wrote', figure(v, curves, summ, wd))


if __name__ == '__main__':
    main()

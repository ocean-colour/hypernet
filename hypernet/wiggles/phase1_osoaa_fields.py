"""Phase 1, task 5: the smooth OSOAA fields for the twin experiment.

On a 5 nm grid from 380 to 1000 nm, for each case (SZA, aerosol, water), one
OSOAA run per wavelength (``hypernet.rt.osoaa``; ~10 s) gives:

- ``ed_dir``, ``ed_dif``: the direct and diffuse downwelling irradiance at
  0+, in units of the extraterrestrial irradiance E_sun (Flux.txt / pi);
- ``ld``: the sky radiance at 0+ seen by HYPSTAR (40 deg from zenith,
  relative azimuth 90), in units of E_sun (I / pi);
- ``lu0m``: the upwelling radiance at 0- along the refracted water view
  (asin(sin 40 / n)), and ``lw`` = lu0m (1 - rho_F(internal)) / n^2, the
  water-leaving radiance at 0+ (Phase 1 Q&A Q5);
- ``lu0p``: the total upward radiance at 0+ at VZA 40 (Lw plus surface
  reflection), kept as a check.  rho_eff = (lu0p - lw) / ld.

OSOAA has no gas absorption, so these are smooth.  The twin experiment
multiplies them by the high-resolution F0 x T of ``hypernet.emod``:
Ed = F0 (ed_dir T_direct + ed_dif T_diffuse), and similarly for Ld.

**Below 400 nm.**  OSOAA's phytoplankton absorption table
(``fic/OSOAA_SEA_PHYT_COEFFS.txt``, Bricaud) starts at 400 nm; the pure-water
table starts at 200 nm.  Below 400 nm OSOAA runs without pigment absorption,
so rho_w is far too high there, and through atmosphere-ocean coupling Ed and
Ld are a few per cent off (400 nm itself, on the table edge, is affected
too).  The assembled fields therefore replace every quantity at
lambda < 405 nm with a quadratic-in-log extrapolation fitted over 405-445 nm
(``extrapolated`` mask; the raw OSOAA values are kept as ``raw_*``).  The
fork is not edited.

Cases (Q&A Q4): the VEIT-like case first (SZA 40, AOT(550) 0.1, Chl 1), then
SZA {30, 50, 70} x AOT(550) {0.05, 0.25} x water {Chl 0.1, 1, 10 mg m-3;
turbid: Chl 1 + 5 mg L-1 sediment + YS 0.1 m-1}.  Each case is saved on its
own (resumable) under ``$OS_COLOR/hypernet/wiggles/phase1/osoaa_cases/``,
with its own Mie database (so that parallel cases do not race), then
assembled into ``osoaa_fields.npz``.  Run outputs are deleted after parsing.

Run from the repository root::

    python -m hypernet.wiggles.phase1_osoaa_fields --veit           # one case, timed
    python -m hypernet.wiggles.phase1_osoaa_fields --grid --workers 12
    python -m hypernet.wiggles.phase1_osoaa_fields --assemble        # npz + figure
"""
import argparse
import itertools
import os
import shutil
import time

import numpy as np

from hypernet.rt import osoaa
from hypernet.wiggles import FIGDIR

ROOT = os.path.join(os.getenv('OS_COLOR', '.'), 'hypernet', 'wiggles', 'phase1')
CASE_DIR = os.path.join(ROOT, 'osoaa_cases')
LAM = np.arange(380.0, 1000.0 + 1e-9, 5.0)
N_WATER = 1.34
VZA = 40.0                       # HYPSTAR water and sky views
PHI = 90.0                       # relative azimuth of every requested sequence
WATERS = {'chl0.1': dict(chl=0.1), 'chl1': dict(chl=1.0), 'chl10': dict(chl=10.0),
          'turbid': dict(chl=1.0, csed=5.0, ys440=0.1)}
VEIT_CASE = dict(name='veit_sza40_aot0.10_chl1', sza=40.0, aot550=0.10, water='chl1')


def grid_cases():
    out = []
    for sza, aot, w in itertools.product((30.0, 50.0, 70.0), (0.05, 0.25), WATERS):
        out.append(dict(name='sza%02.0f_aot%.2f_%s' % (sza, aot, w), sza=sza, aot550=aot,
                        water=w))
    return out


def fresnel(theta_i_deg, n1, n2):
    """Unpolarised Fresnel reflectance for incidence ``theta_i`` from medium n1."""
    ti = np.radians(theta_i_deg)
    st = n1 / n2 * np.sin(ti)
    if st >= 1:
        return 1.0
    tt = np.arcsin(st)
    rs = (n1 * np.cos(ti) - n2 * np.cos(tt)) / (n1 * np.cos(ti) + n2 * np.cos(tt))
    rp = (n1 * np.cos(tt) - n2 * np.cos(ti)) / (n1 * np.cos(tt) + n2 * np.cos(ti))
    return 0.5 * (rs ** 2 + rp ** 2)


def one_wavelength(lam, case, mie_dir, tmp_root):
    w = WATERS[case['water']]
    work = os.path.join(tmp_root, '%s_%06.1f' % (case['name'], lam))
    p = osoaa.default_params(lam, case['sza'], work, mie_dir=mie_dir, phi=PHI, level=3,
                             aot550=case['aot550'], chl=w.get('chl', 0.0),
                             csed=w.get('csed', 0.0), ys440=w.get('ys440', 0.0))
    try:
        r = osoaa.run(p)
    finally:
        pass
    f, up, dn = r['flux'], r['up'], r['down']
    l0p, l0m = dn['level_0plus'], dn['level_0minus']
    i0p = np.where(f['level'] == l0p)[0][0]
    th_w = np.degrees(np.arcsin(np.sin(np.radians(VZA)) / N_WATER))
    lu0m = osoaa.radiance_at(up, l0m, th_w) / np.pi
    out = dict(ed_dir=f['Ed_direct'][i0p] / np.pi, ed_dif=f['Ed_diffuse'][i0p] / np.pi,
               ld=osoaa.radiance_at(dn, l0p, VZA) / np.pi, lu0m=lu0m,
               lw=lu0m * (1 - fresnel(th_w, N_WATER, 1.0)) / N_WATER ** 2,
               lu0p=osoaa.radiance_at(up, l0p, VZA) / np.pi, seconds=r['seconds'])
    shutil.rmtree(work, ignore_errors=True)
    return out


def run_case(case, lam=LAM, verbose=False):
    """Run (or resume) one case; returns the path of its npz."""
    os.makedirs(CASE_DIR, exist_ok=True)
    path = os.path.join(CASE_DIR, case['name'] + '.npz')
    if os.path.exists(path):
        return path
    mie = os.path.join(ROOT, 'osoaa_db', case['name'])
    tmp = os.path.join(ROOT, 'osoaa_tmp')
    keys = ('ed_dir', 'ed_dif', 'ld', 'lu0m', 'lw', 'lu0p', 'seconds')
    res = {k: np.full(lam.size, np.nan) for k in keys}
    t0 = time.time()
    for i, l in enumerate(lam):
        try:
            o = one_wavelength(l, case, mie, tmp)
            for k in keys:
                res[k][i] = o[k]
        except Exception as ex:
            print('%s %.0f nm: %s' % (case['name'], l, str(ex)[:200]), flush=True)
        if verbose and (i < 3 or i % 25 == 0):
            print('%s %.0f nm: %.1f s (elapsed %.0f s)' % (case['name'], l, res['seconds'][i],
                                                          time.time() - t0), flush=True)
    np.savez(path, lam=lam, **res, **{'case_' + k: v for k, v in case.items()})
    print('%s done in %.0f s' % (case['name'], time.time() - t0), flush=True)
    return path


def assemble():
    cases = [VEIT_CASE] + grid_cases()
    found = [c for c in cases if os.path.exists(os.path.join(CASE_DIR, c['name'] + '.npz'))]
    keys = ('ed_dir', 'ed_dif', 'ld', 'lu0m', 'lw', 'lu0p', 'seconds')
    arr = {k: [] for k in keys}
    for c in found:
        d = np.load(os.path.join(CASE_DIR, c['name'] + '.npz'))
        for k in keys:
            arr[k].append(d[k])
    out = {k: np.array(v) for k, v in arr.items()}
    # below 400 nm: log-linear extrapolation from 400-440 nm (see the docstring)
    extrap = LAM < 405.0                     # 400 nm itself sits on the table edge
    fit = (LAM >= 405.0) & (LAM <= 445.0)
    for k in ('ed_dir', 'ed_dif', 'ld', 'lu0m', 'lw', 'lu0p'):
        out['raw_' + k] = out[k].copy()
        for i in range(out[k].shape[0]):
            y = out[k][i]
            # quadratic in ln(y): continuous curvature into the OSOAA values,
            # which matters because the twin experiment measures rho_w''
            c = np.polyfit(LAM[fit], np.log(np.clip(y[fit], 1e-30, None)), 2)
            out[k][i, extrap] = np.exp(np.polyval(c, LAM[extrap]))
    out['extrapolated'] = extrap
    meta = dict(name=np.array([c['name'] for c in found]),
                sza=np.array([c['sza'] for c in found]),
                aot550=np.array([c['aot550'] for c in found]),
                water=np.array([c['water'] for c in found]))
    path = os.path.join(ROOT, 'osoaa_fields.npz')
    np.savez(path, lam=LAM, vza=VZA, phi=PHI, n_water=N_WATER,
             units='E_sun (irradiance: Flux/pi; radiance: I/pi)', **out, **meta)
    print('assembled %d cases -> %s' % (len(found), path))
    return path


def clear_floor_check(path):
    """Modelled clear-sky Ld(750)/Ed(750) vs the observed clear floor 0.011-0.015."""
    d = np.load(path)
    j = np.argmin(np.abs(d['lam'] - 750))
    print('\nclear-sky check, Ld(750)/Ed(750) in sr-1 (observed clear floor 0.011-0.015):')
    for i, n in enumerate(d['name']):
        ed = d['ed_dir'][i, j] + d['ed_dif'][i, j]
        print('  %-28s %.4f   (diffuse fraction %.2f)' % (n, d['ld'][i, j] / ed,
                                                         d['ed_dif'][i, j] / ed))


def figure(path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    d = np.load(path)
    lam = d['lam']
    names = list(d['name'])
    sel = [n for n in names if n.startswith('sza50') or n.startswith('veit')]
    cols = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7',
            '#e34948', '#333333']
    fig, axs = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    for k, n in enumerate(sel):
        i = names.index(n)
        c = cols[k % len(cols)]
        ed = d['ed_dir'][i] + d['ed_dif'][i]
        axs[0, 0].plot(lam, ed, color=c, lw=1.5, label=n)
        axs[0, 1].plot(lam, d['ed_dif'][i] / ed, color=c, lw=1.5)
        axs[1, 0].plot(lam, d['ld'][i] / ed, color=c, lw=1.5)
        axs[1, 1].semilogy(lam, np.pi * d['lw'][i] / ed, color=c, lw=1.5)
    axs[1, 0].axhspan(0.011, 0.015, color='#888', alpha=0.25, lw=0)
    axs[1, 0].text(760, 0.016, 'observed clear floor of Ld/Ed at 750 nm', fontsize=8,
                   color='#555')
    labs = [('Ed(0+) / E_sun', 'total irradiance'), ('Ed_diffuse / Ed', 'diffuse fraction'),
            ('Ld / Ed (sr⁻¹)', 'sky radiance, VZA 40°, RAA 90°'),
            ('ρw = π Lw / Ed', 'water-leaving reflectance')]
    for ax, (yl, t) in zip(axs.ravel(), labs):
        ax.set_ylabel(yl)
        ax.set_title(t, fontsize=10, loc='left')
        ax.grid(alpha=0.25, lw=0.6)
        ax.spines[['top', 'right']].set_visible(False)
    for ax in axs[1]:
        ax.set_xlabel('wavelength (nm)')
    axs[0, 0].legend(fontsize=7, frameon=False, ncol=1)
    fig.suptitle('OSOAA smooth fields (no gas absorption): SZA 50° cases and the VEIT-like '
                 'case', fontsize=10.5, x=0.06, ha='left')
    fig.tight_layout()
    os.makedirs(FIGDIR.replace('phase0', 'phase1'), exist_ok=True)
    out = os.path.join(FIGDIR.replace('phase0', 'phase1'), 'osoaa_fields.png')
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print('wrote', out)


def main(argv=None):
    ap = argparse.ArgumentParser(description='Phase 1 task 5: OSOAA smooth fields')
    ap.add_argument('--veit', action='store_true', help='run the VEIT-like case (timed)')
    ap.add_argument('--grid', action='store_true', help='run the 24-case grid')
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--assemble', action='store_true', help='npz, check and figure')
    a = ap.parse_args(argv)
    if a.veit:
        run_case(VEIT_CASE, verbose=True)
    if a.grid:
        from multiprocessing import Pool
        cases = [c for c in grid_cases()
                 if not os.path.exists(os.path.join(CASE_DIR, c['name'] + '.npz'))]
        print('%d cases to run on %d workers' % (len(cases), a.workers), flush=True)
        with Pool(a.workers) as pool:
            for p in pool.imap_unordered(run_case, cases):
                print('saved', p, flush=True)
    if a.assemble:
        path = assemble()
        clear_floor_check(path)
        figure(path)


if __name__ == '__main__':
    main()

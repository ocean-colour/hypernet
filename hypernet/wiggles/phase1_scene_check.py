"""Phase 1, task 6 check: one twin scene, and the true rho_w with each control.

For the VEIT-like OSOAA case (SZA 40, AOT 0.10, Chl 1) and each water type,
compose the 0.01 nm scene (``hypernet.twin``) and the true rho_w on the VEIT L
grid with the shipped Ld SRF (``hypernet/data/veit_srf_model.json``).  Prints
the size of the fluorescence and Raman controls, and draws
``hypernet/wiggles/figs/phase1/scene_check.png``.

Run from the repository root: ``python -m hypernet.wiggles.phase1_scene_check``.
"""
import os

import numpy as np

from hypernet import emod, srf, twin
from hypernet import whn_l1a as wl
from hypernet.wiggles import DATA_DIR, FIGDIR

CASE = 'veit_sza40_aot0.10_chl1'


def main():
    fields = twin.load_fields()
    lib = twin.rhow_library(fields)
    e = emod.build_emod(40.0)
    f = wl.sequence_files(site='VEIT', seq_time='20260604T0845')
    grid_L = wl.load_l1a_rad(f['L1A_RAD'])['wave']
    grid_L = grid_L[(grid_L > 385) & (grid_L < 995)]
    srf_L = srf.load_srf_models(os.path.join(DATA_DIR, 'veit_srf_model.json'))['Ld']
    out = {}
    for w in lib:
        sc = twin.compose_scene(e, fields, CASE, water=w, lib=lib)
        rho = {}
        for tag, ctl in (('elastic', ()), ('fl', ('fl',)), ('raman', ('raman',)),
                         ('all', ('fl', 'raman'))):
            s2 = twin.compose_scene(e, fields, CASE, water=w, controls=ctl, lib=lib)
            rho[tag] = twin.rhow_true(s2, srf_L, grid_L)
        out[w] = (sc, rho)
        g = (grid_L > 540) & (grid_L < 560)
        print('%-7s rho_el(550) %.4f  fluorescence peak %.2e  Raman/total at 550 %.2f  '
              'Raman/total at 650 %.2f' % (
                  w, rho['elastic'][g].mean(), (rho['fl'] - rho['elastic']).max(),
                  ((rho['raman'] - rho['elastic']) / rho['raman'])[g].mean(),
                  ((rho['raman'] - rho['elastic']) / rho['raman'])[
                      (grid_L > 640) & (grid_L < 660)].mean()))
    figure(out, grid_L)


def figure(out, grid_L):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    sc, rho = out['chl1']
    lam = sc['lam']
    fig, axs = plt.subplots(2, 2, figsize=(11, 7.5))
    for ax, (lo, hi), t in ((axs[0, 0], (425, 440), 'G band'), (axs[0, 1], (755, 772), 'O₂-A')):
        m = (lam > lo) & (lam < hi)
        ax.plot(lam[m], sc['Ed'][m] / sc['Ed'][m].max(), color='#2a78d6', lw=0.6, label='Ed')
        ax.plot(lam[m], sc['Ld'][m] / sc['Ld'][m].max(), color='#eb6834', lw=0.6, label='Ld')
        ax.set_title('0.01 nm scene: %s (normalised)' % t, fontsize=10, loc='left')
        ax.set_xlabel('wavelength (nm)')
        ax.legend(fontsize=8, frameon=False)
    cols = {'chl0.1': '#2a78d6', 'chl1': '#eb6834', 'chl10': '#1baf7a', 'turbid': '#333333'}
    for w, (s_, r) in out.items():
        axs[1, 0].semilogy(grid_L, r['all'], color=cols[w], lw=1, label=w)
        axs[1, 1].plot(grid_L, r['fl'] - r['elastic'], color=cols[w], lw=1.2)
        axs[1, 1].plot(grid_L, r['raman'] - r['elastic'], color=cols[w], lw=1, ls=':')
    axs[1, 0].set_title('true ρw (L SRF, VEIT L grid), with controls', fontsize=10, loc='left')
    axs[1, 0].legend(fontsize=8, frameon=False)
    axs[1, 1].set_title('controls: fluorescence (solid), Raman (dotted)', fontsize=10, loc='left')
    for ax in axs[1]:
        ax.set_xlabel('wavelength (nm)')
    axs[1, 1].set_ylabel('Δρw')
    for ax in axs.ravel():
        ax.grid(alpha=0.25, lw=0.6)
        ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    path = os.path.join(FIGDIR.replace('phase0', 'phase1'), 'scene_check.png')
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print('wrote', path)


if __name__ == '__main__':
    main()

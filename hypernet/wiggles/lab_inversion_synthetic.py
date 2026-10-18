"""How the inversion's ringing index behaves on synthetic spectra (lab data
prompts, task 5): which line density makes a wrong E grid detectable, to set
up the Tier-1 test honestly.

For synthetic E (continuum x N absorption lines of 0.01-0.3 nm, seen at 2.3
nm FWHM; optionally white relative noise per E pixel, as the scan noise of a
real E) on 122304's grids: the shift scan delta_min and the change of the
plain and Nyquist second-difference power for a grid error of +-0.05 nm.
Noise-free, a wrong grid only re-interpolates a smooth spectrum and nothing
is detectable; the signal comes from the E pixel noise, whose Nyquist gain
through M_wrong^-1 M_true depends on the grid phase.

Run from the repository root:  ``python -m hypernet.wiggles.lab_inversion_synthetic``
"""
import numpy as np

from hypernet import edinterp, labcal
from hypernet.srf import convolve_to_grid

DELTAS = np.round(np.arange(-0.10, 0.1001, 0.01), 3)


def native(wE, n_lines, seed=3):
    lam = np.arange(340.0, 1110.0, 0.005)
    rng = np.random.default_rng(seed)
    y = 1.0 + 0.2 * np.sin(lam / 60.0)
    for mu in rng.uniform(345, 1105, n_lines):
        y *= 1 - rng.uniform(0.05, 0.9) * np.exp(-0.5 * ((lam - mu) / rng.uniform(0.01, 0.3)) ** 2)
    return convolve_to_grid(lam, y, wE, 2.3)


def main():
    wE, wL = labcal.grid_for('122304', 'E'), labcal.grid_for('122304', 'L')
    sdp = edinterp.second_difference_power
    for n, noise in ((2000, 0.0), (2000, 1e-4), (2000, 1e-3), (2000, 3e-3)):
        E = native(wE, n) * (1 + noise * np.random.default_rng(7).standard_normal(wE.size))
        EL = edinterp.interp_matrix(wE, wL) @ E
        def ring(g):
            return sdp(edinterp.recover_native_e(EL, wL, g), g, nyquist=True)
        out = []
        for err in (0.0, 0.05, -0.05):
            g = wE + err
            sc = np.array([ring(g + d) for d in DELTAS])
            Eg = edinterp.recover_native_e(EL, wL, g)
            out.append('err %+.2f: delta_min %+.2f, sdp x%.3f, nyquist x%.3f'
                       % (err, DELTAS[np.argmin(sc)], sdp(Eg, g) / sdp(E, wE),
                          sdp(Eg, g, nyquist=True) / sdp(E, wE, nyquist=True)))
        print('%5d lines, noise %.0e: ' % (n, noise) + ' | '.join(out))


if __name__ == '__main__':
    main()

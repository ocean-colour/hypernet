"""Tests for :mod:`hypernet.edinterp` (Phase 1 task 8).

All synthetic: a 0.01 nm 'solar' spectrum with many narrow lines, HYPSTAR-like
E and L grids with a drifting sub-pixel offset.
"""

import numpy as np
import pytest

from hypernet import edinterp
from hypernet.twin import convolve_to_grid

LAM = np.arange(390.0, 710.0, 0.01)
GE = np.arange(400.0, 700.0, 0.471) + 0.11                    # E grid
GL = np.arange(401.0, 698.0, 0.483)                           # L grid (drifting phase)


def _emod(seed=2):
    rng = np.random.default_rng(seed)
    f = 1.5 + 0.001 * (LAM - 550)
    for mu in rng.uniform(392, 708, 900):
        f = f * (1 - rng.uniform(0.05, 0.7) * np.exp(-0.5 * ((LAM - mu) / rng.uniform(0.01, 0.06)) ** 2))
    return LAM, f


def _measured(emod, fwhm_E, scale=None):
    """E observed by the instrument: Emod x a smooth 'atmosphere/calibration'."""
    lam, f = emod
    smooth = 1 + 0.2 * np.sin((lam - 400) / 90.0) if scale is None else scale
    return convolve_to_grid(lam, f * smooth, GE, fwhm_E)


def test_constant_emod_equals_linear():
    E = _measured(_emod(), 2.3)
    const = (LAM, np.full(LAM.size, 1.7))
    ref = np.interp(GL, GE, E)
    for m in ('linear', 'ruddick2023', 'srf'):
        out = edinterp.interpolate_ed_to_l(GE, E, GL, emod=const, srf_irr=2.3, srf_rad=2.8,
                                           method=m)
        np.testing.assert_allclose(out, ref, rtol=1e-12, atol=0)


def test_srf_with_equal_srfs_is_eq14():
    em = _emod()
    E = _measured(em, 2.3)
    out = edinterp.interpolate_ed_to_l(GE, E, GL, emod=em, srf_irr=2.3, srf_rad=2.3,
                                       method='srf')
    # direct implementation of Ruddick et al. (2023) eq. 14, pixel by pixel
    emE_at_E = convolve_to_grid(em[0], em[1], GE, 2.3)
    emE_at_L = convolve_to_grid(em[0], em[1], GL, 2.3)
    direct = np.empty(GL.size)
    for j, lj in enumerate(GL):
        i = np.searchsorted(GE, lj) - 1
        w = (lj - GE[i]) / (GE[i + 1] - GE[i])
        direct[j] = emE_at_L[j] * ((1 - w) * E[i] / emE_at_E[i] +
                                   w * E[i + 1] / emE_at_E[i + 1])
    np.testing.assert_allclose(out, direct, rtol=1e-12, atol=0)
    r23 = edinterp.interpolate_ed_to_l(GE, E, GL, emod=em, srf_irr=2.3, method='ruddick2023')
    np.testing.assert_allclose(out, r23, rtol=1e-14, atol=0)


def test_spectrum_equal_to_model_returns_model():
    em = _emod()
    emE = convolve_to_grid(em[0], em[1], GE, 2.3)
    out = edinterp.interpolate_ed_to_l(GE, emE, GL, emod=em, srf_irr=2.3, srf_rad=2.8,
                                       method='srf')
    np.testing.assert_allclose(out, convolve_to_grid(em[0], em[1], GL, 2.8), rtol=1e-12)
    # and alpha x Emod_E returns alpha x Emod_L (the paper's property)
    out2 = edinterp.interpolate_ed_to_l(GE, 3.0 * emE, GL, emod=em, srf_irr=2.3, srf_rad=2.8,
                                        method='srf')
    np.testing.assert_allclose(out2, 3.0 * convolve_to_grid(em[0], em[1], GL, 2.8), rtol=1e-12)


def test_srf_removes_the_mismatch_better_than_linear():
    """E at FWHM 2.3, L at 2.8: dividing an L-resolution spectrum by Ed_L
    leaves line residuals for linear interpolation, almost none for srf."""
    em = _emod()
    E = _measured(em, 2.3, scale=1.0)
    truth_L = convolve_to_grid(em[0], em[1], GL, 2.8)
    res = {}
    for m in ('linear', 'ruddick2023', 'srf'):
        out = edinterp.interpolate_ed_to_l(GE, E, GL, emod=em, srf_irr=2.3, srf_rad=2.8,
                                           method=m)
        res[m] = np.std(truth_L / out - 1)
    assert res['srf'] < 1e-6 < res['ruddick2023'] < res['linear'] * 1.5
    assert res['srf'] < 0.01 * res['linear']


def test_nulls_and_errors():
    x = np.sin(GE / 40.0) + 0.001 * GE
    smooth = np.sin(GL / 40.0) + 0.001 * GL
    for m in ('cubic', 'sinc'):
        out = edinterp.interpolate_ed_to_l(GE, x, GL, method=m)
        assert np.max(np.abs(out - smooth)[(GL > 410) & (GL < 690)]) < 2e-3
    with pytest.raises(ValueError):
        edinterp.interpolate_ed_to_l(GE, x, GL, method='srf')
    with pytest.raises(ValueError):
        edinterp.interpolate_ed_to_l(GE, x, GL, method='spline')


def test_per_pixel_srf_arrays():
    em = _emod()
    E = _measured(em, 2.3)
    fE = np.full(GE.size, 2.3)
    fL = np.full(GL.size, 2.8)
    a = edinterp.interpolate_ed_to_l(GE, E, GL, emod=em, srf_irr=fE, srf_rad=fL, method='srf')
    b = edinterp.interpolate_ed_to_l(GE, E, GL, emod=em, srf_irr=2.3, srf_rad=2.8, method='srf')
    np.testing.assert_allclose(a, b, rtol=1e-13)
    c = edinterp.interpolate_ed_to_l(GE, E, GL, emod=em, srf_irr=fE, method='ruddick2023')
    d = edinterp.interpolate_ed_to_l(GE, E, GL, emod=em, srf_irr=2.3, method='ruddick2023')
    np.testing.assert_allclose(c, d, rtol=1e-13)

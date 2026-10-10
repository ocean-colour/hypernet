"""Tests for :mod:`hypernet.edinterp` (Phase 1 task 8; hardened in Phase 2 task 3).

All synthetic: a 0.01 nm 'solar' spectrum with many narrow lines, HYPSTAR-like
E and L grids with a drifting sub-pixel offset.
"""

import numpy as np
import pytest

import time

from hypernet import edinterp, srf
from hypernet.srf import convolve_to_grid

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


# --- Phase 2 task 3 ----------------------------------------------------------------

def _model(ch, c0):
    return srf.SRFModel(ch, [c0, -0.05, 0.01], np.eye(3) * 1e-4)


def test_2d_irradiance_matches_columns():
    em = _emod()
    E1 = _measured(em, 2.3)
    E = np.stack([E1, 1.1 * E1, 0.9 * E1 + 0.01], axis=1)
    for m in ('linear', 'ruddick2023', 'srf'):
        out = edinterp.interpolate_ed_to_l(GE, E, GL, emod=em, srf_irr=2.3, srf_rad=2.8, method=m)
        assert out.shape == (GL.size, 3)
        for k in range(3):
            ref = edinterp.interpolate_ed_to_l(GE, E[:, k], GL, emod=em, srf_irr=2.3,
                                               srf_rad=2.8, method=m)
            np.testing.assert_allclose(out[:, k], ref, rtol=1e-14)


def test_validation():
    E = np.ones(GE.size)
    bad = [dict(wav_irr=GE[::-1]), dict(wav_irr=np.r_[GE[:5], GE[4], GE[6:]]),
           dict(irradiance=np.ones(GE.size - 1)), dict(irradiance=np.ones((GE.size, 2, 2))),
           dict(wav_rad=np.r_[GL[:-1], np.nan]), dict(method='spline')]
    for b in bad:
        kw = dict(wav_irr=GE, irradiance=E, wav_rad=GL, method='linear')
        kw.update(b)
        with pytest.raises(ValueError):
            edinterp.interpolate_ed_to_l(kw.pop('wav_irr'), kw.pop('irradiance'),
                                         kw.pop('wav_rad'), **kw)
    with pytest.raises(ValueError):          # uncertainties with a null method
        edinterp.interpolate_ed_to_l(GE, E, GL, method='cubic', u_rel_random=1.0)


def test_nan_stays_local():
    em = _emod()
    E = _measured(em, 2.3)
    k = 300
    E[k] = np.nan
    i, w = edinterp.linear_weights(GE, GL)
    uses = ((i == k) & (w < 1)) | ((i + 1 == k) & (w > 0))
    for m in ('linear', 'srf'):
        out = edinterp.interpolate_ed_to_l(GE, E, GL, emod=em, srf_irr=2.3, srf_rad=2.8, method=m)
        assert np.all(np.isnan(out[uses])) and np.all(np.isfinite(out[~uses]))
        assert 1 <= uses.sum() <= 2


def test_edges_fall_back_to_linear():
    """L pixels beyond the E grid, or whose SRF reaches past Emod, use linear."""
    em = _emod()
    E = _measured(em, 2.3)
    gl = np.r_[391.0, 395.0, GL, 699.9, 705.0]          # Emod spans 390-710
    out, info = edinterp.interpolate_ed_to_l(GE, E, gl, emod=em, srf_irr=2.3, srf_rad=2.8,
                                             method='srf', full_output=True)
    fb = info['fallback']
    assert fb[:2].all() and fb[-2:].all() and not fb[2:-2].any()
    np.testing.assert_allclose(out[fb], np.interp(gl[fb], GE, E), rtol=1e-12)
    inner = edinterp.interpolate_ed_to_l(GE, E, GL, emod=em, srf_irr=2.3, srf_rad=2.8)
    np.testing.assert_allclose(out[2:-2], inner, rtol=1e-14)


def test_choose_method_rule():
    assert edinterp.delta_fwhm(2.3, 2.8) == pytest.approx(0.5)
    assert edinterp.choose_method(2.3, 2.8) == 'srf'
    assert edinterp.choose_method(2.7, 2.8) == 'linear'
    assert edinterp.choose_method(2.65, 2.8) == 'srf'          # 0.15 is inclusive
    assert edinterp.choose_method(2.7, 2.8, threshold=0.05) == 'srf'
    assert edinterp.choose_method(_model('E', 2.3), _model('Ld', 2.8)) == 'srf'
    with pytest.raises(ValueError):
        edinterp.delta_fwhm(np.full(GE.size, 2.3), 2.8)


def test_choose_method_release2_instruments():
    from hypernet import twin
    s = twin.load_instrument_srfs()
    assert edinterp.choose_method(s['HYPSTAR_122304']['E'], s['HYPSTAR_122304']['Ld']) == 'srf'
    assert edinterp.choose_method(s['HYPSTAR_122305']['E'], s['HYPSTAR_122305']['Ld']) == 'linear'


def test_input_uncertainties_linear_analytic():
    E = 1.0 + 0.3 * np.sin(GE / 17.0)
    u = 0.5 + 0.001 * (GE - 400)                                # % per E pixel
    out, info = edinterp.interpolate_ed_to_l(GE, E, GL, method='linear', u_rel_random=u,
                                             u_rel_systematic_indep=u)
    i, w = edinterp.linear_weights(GE, GL)
    a, b = (1 - w) * E[i], w * E[i + 1]
    rand = np.sqrt((a * u[i]) ** 2 + (b * u[i + 1]) ** 2) / (a + b)
    sys_ = (a * u[i] + b * u[i + 1]) / (a + b)
    np.testing.assert_allclose(info['u_rel_random'], rand, rtol=1e-12)
    np.testing.assert_allclose(info['u_rel_systematic_indep'], sys_, rtol=1e-12)
    assert 'u_rel_systematic_corr_rad_irr' not in info
    assert np.all(info['u_rel_random'] <= info['u_rel_systematic_indep'] + 1e-12)


def test_constant_systematic_uncertainty_is_preserved():
    em = _emod()
    E = np.stack([_measured(em, 2.3)] * 2, axis=1)
    for m in ('linear', 'ruddick2023', 'srf'):
        _, info = edinterp.interpolate_ed_to_l(GE, E, GL, emod=em, srf_irr=2.3, srf_rad=2.8,
                                               method=m, u_rel_systematic_corr_rad_irr=2.0)
        assert info['u_rel_systematic_corr_rad_irr'].shape == (GL.size, 2)
        np.testing.assert_allclose(info['u_rel_systematic_corr_rad_irr'], 2.0, rtol=1e-12)


def test_convolution_matrix_cache():
    em = _emod()
    s = _model('E', 2.3)
    M1 = srf.convolution_matrix(em[0], GE, s)
    M2 = srf.convolution_matrix(em[0], GE, _model('E', 2.3))     # equal model, new object
    assert M1 is M2
    assert srf.convolution_matrix(em[0], GE, s.fwhm) is not M1   # callables are not cached
    np.testing.assert_allclose(M1 @ em[1], convolve_to_grid(em[0], em[1], GE, s.fwhm),
                               rtol=1e-14)


def test_timing_1536_to_1538():
    """A full-size srf call: < 1 s cold, < 50 ms with the matrices cached
    (measured 25 ms and 4.5 ms on an M-series Mac)."""
    lam = np.arange(380.0, 1000.0 + 1e-9, 0.01)
    f = 1.5 + 0.2 * np.sin(lam / 7.0)
    gE = np.linspace(389.0, 991.0, 1536)
    gL = np.linspace(390.2, 990.4, 1538)
    sE, sL = _model('E', 2.31), _model('Ld', 2.81)
    E = convolve_to_grid(lam, f, gE, sE.fwhm, cache=False)
    t0 = time.perf_counter()
    edinterp.interpolate_ed_to_l(gE, E, gL, emod=(lam, f), srf_irr=sE, srf_rad=sL)
    cold = time.perf_counter() - t0
    t0 = time.perf_counter()
    out = edinterp.interpolate_ed_to_l(gE, E, gL, emod=(lam, f), srf_irr=sE, srf_rad=sL)
    warm = time.perf_counter() - t0
    assert out.shape == (1538,) and np.all(np.isfinite(out))
    assert cold < 1.0 and warm < 0.05

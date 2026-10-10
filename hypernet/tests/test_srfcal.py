"""Tests for :mod:`hypernet.srfcal` (Phase 2 task 4).

Tier 1: a synthetic 'solar' spectrum with many narrow lines, observed with a
known FWHM(lambda) by :func:`hypernet.twin.observe`.  Tier 2 (``@needs_veit``,
``@needs_hsrs``): the VEIT E spectrum against the shipped template model.
"""
import os

import numpy as np
import pytest

from hypernet import refspec, srf, srfcal, twin
from hypernet import whn_l1a as wl

LAM = np.arange(375.0, 905.0, 0.01)
GRID = np.arange(385.0, 895.0, 0.47) + 0.07
TRUE = srf.SRFModel('E', [2.35, -0.12, 0.04], np.zeros((3, 3)))


def _ref(seed=3):
    rng = np.random.default_rng(seed)
    f = 1.5 + 0.001 * (LAM - 550)
    for mu in rng.uniform(376, 904, 1500):
        f = f * (1 - rng.uniform(0.05, 0.6) *
                 np.exp(-0.5 * ((LAM - mu) / rng.uniform(0.01, 0.05)) ** 2))
    return LAM, f


def _observed(ref, blue_scale=1.0, noise=1e-3, shift=0.0, seed=4):
    rng = np.random.default_rng(seed)
    obs = twin.observe(ref[0], ref[1], GRID, TRUE, true_grid=GRID + shift)
    obs = obs * np.where(GRID < 460, blue_scale, 1.0)
    e = np.full(obs.size, noise * np.median(obs[GRID > 500]))
    return obs + e * rng.standard_normal(obs.size), e


def test_recovers_known_fwhm_within_errors():
    ref = _ref()
    obs, e = _observed(ref)
    m, fits = srfcal.calibrate_srf(GRID, obs, err=e, ref=ref, return_fits=True)
    assert m.meta['ok'] and m.meta['flags'] == []
    assert m.meta['n_windows'] >= 20 and m.meta['n_blue'] >= 3
    lam = np.array([410.0, 500.0, 600.0, 700.0, 850.0])
    pull = (m.fwhm(lam) - TRUE.fwhm(lam)) / m.fwhm_err(lam)
    assert np.all(np.abs(pull) < 3), pull
    assert np.all(m.fwhm_err(lam) < 0.03)
    assert abs(m.offset) < 3 * m.offset_err
    assert len(fits) == len(srf.template_windows())


def test_recovers_centroid_offset():
    """Pixels that really sample at grid + 0.05 nm show the reference 0.05 nm
    early: dlam (measured - reference) = -0.05."""
    ref = _ref()
    obs, e = _observed(ref, shift=0.05)
    m = srfcal.calibrate_srf(GRID, obs, err=e, ref=ref)
    assert m.offset == pytest.approx(-0.05, abs=3 * m.offset_err)
    assert m.offset_err < 0.01


def test_low_blue_signal_is_flagged():
    ref = _ref()
    obs, e = _observed(ref, blue_scale=0.02)              # blue S/N ~ 20 (low sun)
    m = srfcal.calibrate_srf(GRID, obs, err=e, ref=ref)
    assert 'low_blue_signal' in m.meta['flags'] and not m.meta['ok']
    assert m.meta['blue_snr'] < srfcal.MIN_BLUE_SNR
    # without err the flag still follows from the blue windows / FWHM error
    obs2, _ = _observed(ref, blue_scale=0.002)
    m2 = srfcal.calibrate_srf(GRID, obs2, ref=ref)
    assert 'low_blue_signal' in m2.meta['flags']


def test_bad_input():
    with pytest.raises(ValueError):
        srfcal.calibrate_srf(GRID, np.ones(GRID.size - 1), ref=_ref())
    m = srfcal.calibrate_srf(GRID, np.full(GRID.size, np.nan), ref=_ref())
    assert 'fit_failed' in m.meta['flags'] and 'few_windows' in m.meta['flags']


def _veit():
    try:
        f = wl.sequence_files(site='VEIT', seq_time='20260604T0845')
        return f if 'L1A_IRR' in f else None
    except Exception:
        return None


needs_veit = pytest.mark.skipif(_veit() is None, reason='requires the VEIT sample')
needs_hsrs = pytest.mark.skipif(not refspec.hsrs_available(),
                                reason='requires the TSIS-1 HSRS file')


@needs_veit
@needs_hsrs
def test_veit_e_matches_shipped_model():
    """The same procedure on the VEIT E reproduces hypernet/data/veit_srf_model.json."""
    from hypernet import whn_srf
    irr = wl.load_l1a_irr(_veit()['L1A_IRR'])
    mean, e_mean, _ = srf.scan_errors(irr['scans'], flatten_px=whn_srf.FLATTEN_PX)
    m = srfcal.calibrate_srf(irr['wave'], mean, err=e_mean, instrument='HYPSTAR_122304')
    ship = srf.load_srf_models(os.path.join(os.path.dirname(srf.__file__), 'data',
                                            'veit_srf_model.json'))['E']
    lam = np.array([420.0, 500.0, 600.0, 660.0])
    np.testing.assert_allclose(m.fwhm(lam), ship.fwhm(lam), atol=0.01)
    assert m.meta['ok'], m.meta


# --- task 4b: wavelength scale -------------------------------------------------------

def _scaled(ref, shift, stretch, noise=2e-3, seed=5):
    """twin case iv: pixels really sample grid + shift + stretch (grid-630)/240."""
    rng = np.random.default_rng(seed)
    tg = GRID + shift + stretch * (GRID - srfcal.LAM0) / srfcal.HALF
    obs = twin.observe(ref[0], ref[1], GRID, TRUE, true_grid=tg)
    e = noise * obs
    return obs + e * rng.standard_normal(obs.size), e


@pytest.mark.parametrize('shift,stretch', [(0.1, 0.0), (-0.05, 0.08), (0.03, -0.05)])
def test_wavelength_scale_recovered(shift, stretch):
    """Known shift/stretch at VEIT-like noise (0.2 % per pixel): recovered
    within 3 sigma and to < 0.02 nm, with max_err well inside G1's 0.05 nm."""
    ref = _ref()
    obs, e = _scaled(ref, shift, stretch)
    r = srfcal.fit_wavelength_scale(GRID, obs, srf=TRUE, ref=ref, err=e)
    assert r['ok'] and r['n_windows'] >= 20
    for k, v in (('shift', shift), ('stretch', stretch)):
        assert abs(r[k] - v) < 3 * r[k + '_err'] and abs(r[k] - v) < 0.02, (k, r[k], v)
    assert r['max_err'] < srfcal.MAX_SCALE_ERR < srfcal.G1_TOLERANCE
    # correcting the grid removes the error
    corr = srfcal.correct_wavelengths(GRID, r)
    tg = GRID + shift + stretch * (GRID - srfcal.LAM0) / srfcal.HALF
    assert np.max(np.abs(corr - tg)) < 3 * r['max_err']
    # SRF width free gives the same answer
    r2 = srfcal.fit_wavelength_scale(GRID, obs, ref=ref, err=e)
    assert abs(r2['shift'] - r['shift']) < r['shift_err']


def test_wavelength_scale_via_grid():
    """An E resampled onto another grid (as in L2): the shift is that of the
    native grid."""
    ref = _ref()
    native = GRID
    obs, e = _scaled(ref, 0.06, 0.0, noise=1e-3)
    gL = np.arange(386.0, 893.0, 0.49) + 0.21
    r = srfcal.fit_wavelength_scale(gL, np.interp(gL, native, obs), srf=TRUE, ref=ref,
                                    via_grid=native)
    assert abs(r['shift'] - 0.06) < 0.02


def test_wavelength_scale_noisy_is_uncertain_and_falls_back():
    ref = _ref()
    obs, e = _scaled(ref, 0.0, 0.0, noise=0.05)
    r = srfcal.fit_wavelength_scale(GRID, obs, srf=TRUE, ref=ref, err=e)
    assert 'uncertain' in r['flags'] and not r['ok']
    fb = dict(shift=0.01, stretch=-0.02, ok=True)
    use, src = srfcal.scale_or_fallback(r, fb)
    assert use is fb and src == 'fallback'
    good = dict(shift=0.0, stretch=0.0, ok=True)
    assert srfcal.scale_or_fallback(good, fb) == (good, 'sequence')


def test_wavelength_scale_pattern_option():
    """A per-window pattern is subtracted from -dlam before the line fit: a
    constant pattern c moves the shift by exactly -c."""
    ref = _ref()
    obs, e = _scaled(ref, 0.04, 0.0)
    r0 = srfcal.fit_wavelength_scale(GRID, obs, srf=TRUE, ref=ref, err=e)
    const = {c: 0.02 for c in (0.5 * (a + b) for a, b in srf.template_windows())}
    r1 = srfcal.fit_wavelength_scale(GRID, obs, srf=TRUE, ref=ref, err=e, pattern=const)
    assert r1['shift'] == pytest.approx(r0['shift'] - 0.02, abs=1e-10)
    assert r1['stretch'] == pytest.approx(r0['stretch'], abs=1e-10)


def test_shipped_window_pattern():
    pat = srfcal.load_window_pattern()
    wins = {0.5 * (a + b) for a, b in srf.template_windows()}
    assert len(pat) >= 20 and set(pat) <= wins
    v = np.array(list(pat.values()))
    assert 0.05 < np.sqrt(np.mean(v ** 2)) < 0.3          # ~0.13 nm rms (task 4b)

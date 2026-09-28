"""Tests for the SRF line-fitting layer (:mod:`hypernet.srf`).

All data-independent: synthetic lines on a HYPSTAR-like grid (~0.48 nm).

Run from the repository root::

    pytest -q hypernet/tests/test_srf.py
"""

import numpy as np
import pytest

from hypernet import srf

# A HYPSTAR-like grid, deliberately not aligned with the line centres.
WAV = np.arange(350.0, 1100.0, 0.483) + 0.117


def _lines(x, comps, c0=100.0, c1=0.2, xref=600.0):
    """Linear continuum times Gaussian absorption lines (depth, mu, sigma)."""
    y = c0 + c1 * (x - xref)
    for a, mu, s in comps:
        y = y * (1.0 - a * np.exp(-0.5 * ((x - mu) / s) ** 2))
    return y


def test_air_to_vac():
    # H alpha and Na D2: air -> vacuum (NIST)
    assert srf.air_to_vac(656.281) == pytest.approx(656.461, abs=2e-3)
    assert srf.air_to_vac(588.995) == pytest.approx(589.158, abs=2e-3)
    assert np.all(srf.LINES['lam_vac'] > srf.LINES['lam_air'])


def test_lines_table():
    L = srf.LINES
    assert list(L.columns) == ['name', 'lam_air', 'lam_vac', 'half', 'group',
                               'blend', 'use_for_srf']
    # Ca H/K are one blended group; the bands are diagnostic only
    cahk = L[L['group'] == 'CaHK']
    assert set(cahk['name']) == {'Ca K', 'Ca H'} and cahk['blend'].all()
    assert not L.loc[L['name'].isin(['H2O', 'O2 A', 'O2 B']), 'use_for_srf'].any()
    assert L['name'].is_unique


def test_single_line_sloped_continuum():
    mu, sig = 500.13, 1.2
    y = _lines(WAV, [(0.3, mu, sig)], c1=0.5)
    r = srf.fit_line(WAV, y, 500.0)
    assert r['ok']
    assert r['mu'] == pytest.approx(mu, abs=0.01)
    assert r['fwhm'] == pytest.approx(srf.FWHM_PER_SIGMA * sig, abs=0.01)
    assert r['depth'] == pytest.approx(0.3, abs=1e-3)


def test_weighted_fit_errors_and_chi2():
    rng = np.random.default_rng(42)
    mu, sig = 656.35, 1.25
    y0 = _lines(WAV, [(0.25, mu, sig)])
    e = np.full_like(WAV, 0.1)
    y = y0 + rng.normal(0.0, 0.1, WAV.size)
    r = srf.fit_line(WAV, y, 656.281, err=e)
    assert r['ok']
    assert 0.3 < r['chi2_nu'] < 2.5
    assert abs(r['mu'] - mu) < 4 * r['mu_err']
    assert abs(r['sigma'] - sig) < 4 * r['sigma_err']
    assert r['fwhm_err'] == pytest.approx(srf.FWHM_PER_SIGMA * r['sigma_err'])


def test_blended_pair_recovers_both_centroids():
    comps = [(0.55, 393.42, 0.95), (0.45, 396.80, 1.05)]
    y = _lines(WAV, comps, c1=-0.3)
    res = srf.fit_blend(WAV, y, [393.366, 396.847], half=4.0)
    assert all(r['ok'] for r in res)
    for r, (a, mu, s) in zip(res, comps):
        assert r['mu'] == pytest.approx(mu, abs=0.01)
        assert r['sigma'] == pytest.approx(s, abs=0.01)
    # the continuum and chi2 are shared
    assert res[0]['c1'] == res[1]['c1']


def test_failures_return_nan():
    y = _lines(WAV, [(0.3, 500.0, 1.2)])
    # window off the grid
    r = srf.fit_line(WAV, y, 200.0)
    assert not r['ok'] and np.isnan(r['mu']) and np.isnan(r['fwhm'])
    # all-NaN spectrum
    r = srf.fit_line(WAV, np.full_like(WAV, np.nan), 500.0)
    assert not r['ok'] and np.isnan(r['sigma'])
    # an emission line is not an absorption line
    r = srf.fit_line(WAV, _lines(WAV, [(-0.3, 500.0, 1.2)]), 500.0)
    assert not r['ok'] and np.isnan(r['depth'])


def test_fit_lines_table():
    L = srf.LINES
    comps = [(0.3, lam + 0.05, 1.1 + 0.001 * (lam - 400.0))
             for lam in L['lam_air']]
    y = _lines(WAV, comps)
    # cut the grid at 800 nm: the red lines must come back NaN, not raise
    keep = WAV < 800.0
    df = srf.fit_lines(WAV[keep], y[keep])
    assert len(df) == len(L) and list(df['name']) == list(L['name'])
    red = df['lam_air'] > 800.0
    assert not df.loc[red, 'ok'].any() and df.loc[red, 'mu'].isna().all()
    np.testing.assert_allclose(df['fwhm'], srf.FWHM_PER_SIGMA * df['sigma'])
    good = df['ok']
    assert good.sum() >= 7
    np.testing.assert_allclose(df.loc[good, 'dmu'], 0.05, atol=0.03)


def test_scan_errors():
    rng = np.random.default_rng(1)
    scans = 10.0 + rng.normal(0.0, 0.5, (200, 6))  # (wavelength, scan)
    mean, e_mean, e_scan = srf.scan_errors(scans)
    assert mean.shape == (200,)
    np.testing.assert_allclose(e_scan, scans.std(axis=1, ddof=1))
    np.testing.assert_allclose(e_mean, e_scan / np.sqrt(6))

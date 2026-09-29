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


# --- FWHM(lambda) model -------------------------------------------------------

LAM_PTS = np.array([393.4, 430.8, 486.1, 517.3, 589.3, 656.3, 849.8, 866.2])
TRUE = np.array([3.0, 0.4, -0.25])  # FWHM at 600 nm, slope, curvature (x = (lam-600)/100)


def _quad(lam, c=TRUE):
    x = (lam - srf.LAM_REF) / srf.LAM_SCALE
    return c[0] + c[1] * x + c[2] * x ** 2


def test_fwhm_model_exact_quadratic():
    fw = _quad(LAM_PTS)
    coeffs, cov = srf.fit_fwhm_model(LAM_PTS, fw, np.full_like(fw, 0.05))
    np.testing.assert_allclose(coeffs, TRUE, atol=1e-10)
    np.testing.assert_allclose(srf.fwhm_at(LAM_PTS, coeffs), fw, atol=1e-10)
    assert isinstance(srf.fwhm_at(500.0, coeffs), float)
    # at lam_ref the error is sqrt(cov[0, 0])
    assert srf.fwhm_err_at(srf.LAM_REF, cov) == pytest.approx(np.sqrt(cov[0, 0]))


def test_fwhm_model_cov_shrinks_with_err():
    fw = _quad(LAM_PTS)
    _, cov1 = srf.fit_fwhm_model(LAM_PTS, fw, np.full_like(fw, 0.10))
    _, cov2 = srf.fit_fwhm_model(LAM_PTS, fw, np.full_like(fw, 0.05))
    np.testing.assert_allclose(cov2, cov1 / 4.0)
    assert np.all(np.diag(cov2) < np.diag(cov1))


def test_fwhm_model_weights_and_bad_points():
    rng = np.random.default_rng(3)
    err = rng.uniform(0.02, 0.2, LAM_PTS.size)
    fw = _quad(LAM_PTS) + rng.normal(0.0, err)
    fw_bad = np.append(fw, [np.nan, 9.0])
    err_bad = np.append(err, [0.1, 0.0])       # NaN value, zero error: dropped
    lam_bad = np.append(LAM_PTS, [700.0, 750.0])
    c1, _ = srf.fit_fwhm_model(LAM_PTS, fw, err)
    c2, _ = srf.fit_fwhm_model(lam_bad, fw_bad, err_bad)
    np.testing.assert_allclose(c1, c2)
    # too few points -> NaN, not an exception
    c, cov = srf.fit_fwhm_model(LAM_PTS[:2], fw[:2], err[:2])
    assert np.all(np.isnan(c)) and np.all(np.isnan(cov))


def test_srfmodel_json_roundtrip(tmp_path):
    fw = _quad(LAM_PTS)
    coeffs, cov = srf.fit_fwhm_model(LAM_PTS, fw, np.full_like(fw, 0.05))
    m = srf.SRFModel('E', coeffs, cov, offset=0.03, offset_err=0.01,
                     lam_min=393.4, lam_max=866.2, npts=8,
                     instrument='HYPSTAR_122304', meta={'seq': 'SEQ20260604T084543'})
    m2 = srf.SRFModel.from_json(m.to_json())
    np.testing.assert_array_equal(m2.coeffs, m.coeffs)
    np.testing.assert_array_equal(m2.cov, m.cov)
    assert (m2.channel, m2.offset, m2.instrument, m2.meta) == \
        (m.channel, m.offset, m.instrument, m.meta)
    assert np.isnan(m2.chi2_nu)          # NaN survives as JSON null
    assert 'NaN' not in m.to_json()
    # file round trip, single and several models
    p = tmp_path / 'e.json'
    m.to_json(p)
    assert srf.SRFModel.from_json(p).fwhm(500.0) == pytest.approx(m.fwhm(500.0))
    mL = srf.SRFModel('Ld', coeffs + 0.2, cov)
    srf.save_srf_models(tmp_path / 'all.json', [m, mL], meta={'site': 'VEIT'})
    got = srf.load_srf_models(tmp_path / 'all.json')
    assert set(got) == {'E', 'Ld'}
    assert got['Ld'].fwhm(600.0) == pytest.approx(TRUE[0] + 0.2)


def test_srfmodel_from_lines():
    # synthetic spectrum: every line at +0.04 nm, sigma following TRUE
    L = srf.LINES
    comps = [(0.3, lam + 0.04, _quad(lam) / srf.FWHM_PER_SIGMA) for lam in L['lam_air']]
    y = _lines(WAV, comps)
    df = srf.fit_lines(WAV, y, err=np.full_like(WAV, 0.01))
    m = srf.SRFModel.from_lines('E', df, instrument='test')
    use = df['ok'] & df['use_for_srf']
    assert m.npts == use.sum() >= 8
    np.testing.assert_allclose(m.coeffs, TRUE, atol=0.02)
    assert m.offset == pytest.approx(0.04, abs=0.01)
    assert m.lam_min == pytest.approx(393.366) and m.frame == 'air'
    assert m.sigma(600.0) == pytest.approx(TRUE[0] / srf.FWHM_PER_SIGMA, abs=0.01)


def test_scan_errors_flatten_removes_broadband_scatter():
    rng = np.random.default_rng(7)
    base = _lines(WAV, [(0.3, 500.0, 1.2)])
    scale = 1.0 + 0.05 * rng.normal(size=6)          # 5 % broadband changes
    tilt = 1.0 + 0.02 * rng.normal(size=6)[None, :] * (WAV[:, None] - 600) / 300
    noise = 1.0 + 0.004 * rng.normal(size=(WAV.size, 6))  # 0.4 % pixel noise
    scans = base[:, None] * scale[None, :] * tilt * noise
    mean, _, raw = srf.scan_errors(scans)
    mean_f, _, flat = srf.scan_errors(scans, flatten_px=20)
    np.testing.assert_array_equal(mean, mean_f)
    assert np.median(raw / mean) > 0.02
    assert np.median(flat / mean) == pytest.approx(0.004, rel=0.25)

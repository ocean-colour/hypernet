"""Tests for :mod:`hypernet.emod` (gas transmittance so far; Phase 1 task 4
adds ozone and ``build_emod``).

Tier 2: needs the cached HITRAN line lists in
``$OS_COLOR/hypernet/wiggles/ref/hitran`` (skips otherwise; no network use).
"""

import os

import numpy as np
import pytest

from hypernet import emod


def _have_lines():
    try:
        return all(os.path.exists(os.path.join(emod.hitran_dir(), m + '.data'))
                   for m in emod.MOLECULES)
    except Exception:
        return False


needs_hitran = pytest.mark.skipif(not _have_lines(), reason='requires cached HITRAN lines')


def test_columns():
    c = emod.columns(pwv_mm=15.0)
    assert c['O2'] == pytest.approx(0.2095 * 2.15e25, rel=0.02)    # molecules cm-2
    assert c['H2O'] == pytest.approx(5.01e22, rel=0.01)


@needs_hitran
def test_gas_transmittance_o2a():
    lam = np.arange(755.0, 775.0, 0.005)
    T = emod.gas_transmittance(lam, airmass=1.0)
    assert np.all((T >= 0) & (T <= 1))
    band = (lam > 759.5) & (lam < 770.0)
    assert T[band].min() < 0.01                       # saturated O2-A line cores
    assert 0.3 < T[band].mean() < 0.8                 # band-mean transmittance
    assert T[(lam > 755) & (lam < 758)].mean() > 0.99  # clear blueward of the band
    # more air mass, more absorption
    assert emod.gas_transmittance(lam, 2.0)[band].mean() < T[band].mean()


def _have_all():
    from hypernet import refspec
    return _have_lines() and refspec.hsrs_available() and os.path.exists(emod.ozone_path())


needs_emod = pytest.mark.skipif(not _have_all(),
                                reason='requires cached HITRAN lines, HSRS and the O3 file')


def test_kasten_young_airmass():
    assert emod.kasten_young_airmass(0.0) == pytest.approx(1.0, abs=1e-3)
    assert emod.kasten_young_airmass(60.0) == pytest.approx(1.995, abs=0.01)
    assert emod.kasten_young_airmass(85.0) > 10


def test_bin_average_conserves():
    x = np.linspace(0, 10, 10001)
    y = 1 + np.sin(5 * x)
    xo = np.arange(1, 9, 0.5)
    np.testing.assert_allclose(emod._bin_average(x, np.ones_like(x), xo), 1.0)
    # bin mean of sin over [a, b] analytically
    a, b = xo - 0.25, xo + 0.25
    exact = 1 + (np.cos(5 * a) - np.cos(5 * b)) / (5 * (b - a))
    np.testing.assert_allclose(emod._bin_average(x, y, xo), exact, atol=1e-5)  # trapezoid


@needs_emod
def test_ozone_transmittance():
    lam = np.array([450.0, 600.0, 900.0])
    T = emod.ozone_transmittance(lam, airmass=1.0, ozone_du=300.0)
    assert np.all((T > 0) & (T <= 1))
    assert T[1] == pytest.approx(0.96, abs=0.01)        # Chappuis peak, 300 DU, m = 1
    assert T[1] < T[0] and T[1] < T[2]


@needs_emod
def test_build_emod_and_cache(tmp_path):
    e = emod.build_emod(36.7, cache_dir=str(tmp_path))
    lam = e['lam']
    assert lam[0] == pytest.approx(380.0) and lam[-1] == pytest.approx(1000.0)
    np.testing.assert_allclose(np.diff(lam), 0.01, atol=1e-9)           # 0.01 nm grid
    for k in ('T_direct', 'T_diffuse'):
        assert np.all((e[k] >= 0) & (e[k] <= 1))
    j = np.argmin(np.abs(lam - 760.6))
    # O2-A at 760.6 nm, air mass 1.25, 0.01 nm bins: stated range 0.05-0.25
    assert 0.05 < e['T_direct'][j] < 0.25
    assert e['T_diffuse'][j] < e['T_direct'][j]                         # 1.66 > 1.25
    assert e['airmass_direct'] == pytest.approx(1.246, abs=0.005)
    assert 1.5 < e['F0'][np.argmin(np.abs(lam - 550))] < 2.2           # W m-2 nm-1
    assert not e['cached']
    e2 = emod.build_emod(36.7, cache_dir=str(tmp_path))
    assert e2['cached']
    for k in ('lam', 'F0', 'T_direct', 'T_diffuse'):
        np.testing.assert_array_equal(e2[k], e[k])

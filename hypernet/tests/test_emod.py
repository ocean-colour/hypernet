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


# --- Phase 2 task 5: the Emod cache --------------------------------------------------

def _fake_builder(calls):
    """A stand-in for build_emod that records its calls (tier 1)."""
    def build(sza, pwv_mm=15.0, ozone_du=300.0, pressure_hpa=1013.25, cache_dir=None):
        calls.append((sza, pwv_mm, ozone_du, pressure_hpa))
        lam = np.arange(380.0, 1000.0 + 1e-9, 0.01)
        line = 1 - 0.5 * np.exp(-0.5 * ((lam - 656.28) / 0.05) ** 2)
        return dict(lam=lam, F0=1.5 * line, T_direct=np.full(lam.size, 0.9 - 0.001 * sza),
                    T_diffuse=np.full(lam.size, 0.8), sza=sza, pwv_mm=pwv_mm, cached=False)
    return build


def test_snap_grid():
    assert emod.snap(36.7) == (35.0, 15.0)
    assert emod.snap(37.6, 22.0) == (40.0, 30.0)        # log-nearest: 22 is closer to 30
    assert emod.snap(37.6, 20.0) == (40.0, 15.0)
    assert emod.snap(-3.0, 1.0) == (0.0, 7.5)
    assert emod.snap(89.0, 90.0) == (85.0, 30.0)


def test_emod_for_same_bin_same_object():
    emod.clear_memory_cache()
    calls = []
    b = _fake_builder(calls)
    e1 = emod.emod_for(36.7, 14.0, builder=b, cache_dir='fake')
    e2 = emod.emod_for(35.9, 16.0, ozone_du=280.0, pressure_hpa=1000.0, builder=b,
                       cache_dir='fake')
    assert e1 is e2 and e1['E'] is e2['E'] and len(calls) == 1
    assert calls[0] == (35.0, 15.0, emod.OZONE_FIXED, emod.PRESSURE_FIXED)
    assert (e1['sza_bin'], e1['pwv_bin'], e1['ozone_du']) == (35.0, 15.0, 300.0)
    fd = emod.DIFFUSE_FRACTION
    np.testing.assert_allclose(e1['E'], e1['F0'] * ((1 - fd) * e1['T_direct'] + fd * 0.8))
    e3 = emod.emod_for(43.0, builder=b, cache_dir='fake')            # another bin: a miss
    assert e3 is not e1 and len(calls) == 2 and calls[1][0] == 45.0
    emod.clear_memory_cache()


def test_emod_on_grid_cached_and_used_by_edinterp():
    from hypernet import edinterp, srf
    emod.clear_memory_cache()
    e = emod.emod_for(36.7, builder=_fake_builder([]), cache_dir='fake')
    grid = np.arange(400.0, 900.0, 0.48)
    m = srf.SRFModel('E', [2.3, 0.0, 0.0], np.zeros((3, 3)))
    a = emod.emod_on_grid(e, grid, m)
    b = emod.emod_on_grid(e, grid, srf.SRFModel('E', [2.3, 0.0, 0.0], np.zeros((3, 3))))
    assert a is b and not a.flags.writeable
    np.testing.assert_allclose(a, srf.convolve_to_grid(e['lam'], e['E'], grid, 2.3), rtol=1e-13)
    assert edinterp.model_on_grid(e, grid, m) is a                  # edinterp uses the cache
    gl = grid[1:-1] + 0.2
    E = srf.convolve_to_grid(e['lam'], e['E'], grid, 2.3)
    out = edinterp.interpolate_ed_to_l(grid, E, gl, emod=e, srf_irr=m, srf_rad=2.8)
    np.testing.assert_allclose(out, srf.convolve_to_grid(e['lam'], e['E'], gl, 2.8), rtol=1e-10)
    emod.clear_memory_cache()


@needs_emod
def test_emod_for_builds_and_stores(tmp_path):
    """A cache miss builds with build_emod and stores the npz; a second process
    (memory cleared) reads it back."""
    emod.clear_memory_cache()
    e = emod.emod_for(36.7, cache_dir=str(tmp_path))
    files = os.listdir(tmp_path)
    assert len(files) == 1 and files[0].startswith('emod_sza35.00_pwv15.0')
    assert not e['cached']
    emod.clear_memory_cache()
    e2 = emod.emod_for(34.0, cache_dir=str(tmp_path))
    assert e2['cached'] and e2 is not e
    np.testing.assert_array_equal(e2['E'], e['E'])
    emod.clear_memory_cache()

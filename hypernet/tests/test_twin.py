"""Tests for the twin-experiment scene code (:mod:`hypernet.twin`).

Tier 1: synthetic flat fields and Emod.  Tier 2 (``needs_fields``): the real
OSOAA fields npz and the cached Emod.
"""

import os

import numpy as np
import pytest

from hypernet import twin

LAM5 = np.arange(380.0, 1000.0 + 1e-9, 5.0)
LAM = np.arange(380.0, 1000.0 + 1e-9, 0.01)
GRID = np.arange(390.0, 990.0, 0.48)


def _flat_fields(n_case=1):
    one = np.ones((n_case, LAM5.size))
    return dict(lam=LAM5, name=np.array(['flat_chl1']), water=np.array(['chl1']),
                sza=np.array([40.0]), ed_dir=0.6 * one, ed_dif=0.1 * one, ld=0.02 * one,
                lw=0.7 * 0.01 / np.pi * one, lu0p=(0.7 * 0.01 / np.pi + 0.02 * 0.028) * one)


def _flat_emod():
    return dict(lam=LAM, F0=np.full(LAM.size, 1.5), T_direct=np.ones(LAM.size),
                T_diffuse=np.ones(LAM.size))


def _flat_lib():
    return {'chl1': dict(lam5=LAM5, rho_el=np.full(LAM5.size, 0.01), chl=1.0,
                         fl_amp=float(twin.fluorescence_amplitude(1.0)),
                         raman_eff=twin.RAMAN_EFFICIENCY)}


def test_convolve_to_grid_constant_and_area():
    np.testing.assert_allclose(twin.convolve_to_grid(LAM, np.full(LAM.size, 3.0), GRID, 2.5),
                               3.0)
    # a narrow line keeps its equivalent width (area below the continuum)
    y = 1 - 0.5 * np.exp(-0.5 * ((LAM - 600.0) / 0.05) ** 2)
    g = np.arange(560.0, 640.0, 0.1)
    obs = twin.convolve_to_grid(LAM, y, g, 3.0)
    ew_true = 0.5 * 0.05 * np.sqrt(2 * np.pi)
    assert np.trapezoid(1 - obs, g) == pytest.approx(ew_true, rel=0.01)
    # variable FWHM through a callable
    obs2 = twin.convolve_to_grid(LAM, y, g, lambda l: np.full(np.shape(l), 3.0))
    np.testing.assert_allclose(obs, obs2)


def test_flat_scene_gives_flat_rhow():
    sc = twin.compose_scene(_flat_emod(), _flat_fields(), 'flat_chl1', controls=(),
                            lib=_flat_lib())
    np.testing.assert_allclose(sc['Ed'], 1.5 * 0.7)
    rho = twin.rhow_true(sc, 2.8, GRID)
    np.testing.assert_allclose(rho, 0.01, rtol=1e-10)


def test_fluorescence_control_peaks_at_683():
    base = twin.compose_scene(_flat_emod(), _flat_fields(), 'flat_chl1', controls=(),
                              lib=_flat_lib())
    fl = twin.compose_scene(_flat_emod(), _flat_fields(), 'flat_chl1', controls=('fl',),
                            lib=_flat_lib())
    d = twin.rhow_true(fl, 2.8, GRID) - twin.rhow_true(base, 2.8, GRID)
    assert GRID[np.argmax(d)] == pytest.approx(683.0, abs=0.5)
    assert d.max() == pytest.approx(twin.fluorescence_amplitude(1.0), rel=0.02)
    assert d[np.abs(GRID - 683) > 60].max() < 1e-3 * d.max()


def test_raman_excitation_and_shift():
    assert twin.raman_excitation(550.0) == pytest.approx(464.4, abs=0.2)
    # a single sharp 'line' in Ed reappears in the Raman radiance at the
    # emission wavelength of its excitation, broadened
    Ed = np.ones(LAM.size)
    Ed[np.argmin(np.abs(LAM - 464.4))] += 100.0
    r = twin.raman_radiance(LAM, Ed)
    assert LAM[np.argmax(r)] == pytest.approx(550.0, abs=0.5)


def _fields_ok():
    p = os.path.join(os.getenv('OS_COLOR', '.'), 'hypernet', 'wiggles', 'phase1',
                     'osoaa_fields.npz')
    from hypernet import refspec
    return os.path.exists(p) and refspec.hsrs_available()


needs_fields = pytest.mark.skipif(not _fields_ok(), reason='requires the OSOAA fields and HSRS')


@needs_fields
def test_real_scene_veit():
    from hypernet import emod as em
    fields = twin.load_fields()
    lib = twin.rhow_library(fields)
    assert set(lib) == {'chl0.1', 'chl1', 'chl10', 'turbid'}
    e = em.build_emod(40.0)
    sc = twin.compose_scene(e, fields, 'veit_sza40_aot0.10_chl1', lib=lib)
    for k in ('Ed', 'Ld', 'Lw', 'Lu'):
        assert np.all(np.isfinite(sc[k])) and np.all(sc[k] >= 0)
    # the elastic truth equals the smooth rho_w on the L grid (it carries Ed's
    # structure, which cancels in Lw / Ed observed with the same SRF)
    el = dict(sc, Lw=sc['Lw_elastic'])
    rho = twin.rhow_true(el, 3.0, GRID)
    ref = np.interp(GRID, sc['lam'], sc['rho_el'])
    m = (GRID > 420) & (GRID < 900)
    assert np.max(np.abs(rho[m] / ref[m] - 1)) < 0.02
    # Raman is a small fraction of Lw in the green for Chl 1
    j = (sc['lam'] > 540) & (sc['lam'] < 560)
    frac = sc['Lw_raman'][j].mean() / sc['Lw'][j].mean()
    assert 0.005 < frac < 0.5


def test_raman_efficiency_falls_in_the_red():
    from hypernet.rt import osoaa
    import os as _os
    if not _os.path.exists(_os.path.join(osoaa.DEFAULT_ROOT, 'fic')):
        pytest.skip('needs the OSOAA fic tables')
    e = twin.raman_efficiency(np.array([450.0, 550.0, 650.0, 750.0, 950.0]))
    assert e[1] == pytest.approx(twin.RAMAN_EFFICIENCY)
    assert e[2] < e[1] and e[3] < 0.3 * e[1] and e[4] < 0.02 * e[1]
    # pigments only lower it
    e10 = twin.raman_efficiency(np.array([450.0, 550.0]), chl=10.0)
    assert np.all(e10 < e[:2])


# --- instrument model (task 7) --------------------------------------------------------

def test_gaussian_srf_integrates_to_one():
    x = np.arange(-20, 20, 0.001)
    for f in (1.5, 2.3, 3.0):
        assert np.trapezoid(twin.gaussian_srf(x, f), x) == pytest.approx(1.0, abs=1e-9)


def test_delta_line_observed_fwhm3():
    from hypernet import srf
    lam = np.arange(560.0, 640.0, 0.01)
    y = 1 - 0.6 * np.exp(-0.5 * ((lam - 600.0) / 0.003) ** 2)      # ~delta line
    grid = np.arange(570.0, 630.0, 0.48) + 0.13
    obs = twin.observe(lam, y, grid, 3.0)
    r = srf.fit_line(grid, obs, 600.0)
    assert r['ok']
    assert r['fwhm'] == pytest.approx(3.0, abs=0.01)
    assert r['mu'] == pytest.approx(600.0, abs=0.01)


def test_observe_wavelength_error_and_noise():
    from hypernet import srf
    lam = np.arange(560.0, 640.0, 0.01)
    y = 1 - 0.5 * np.exp(-0.5 * ((lam - 600.0) / 0.5) ** 2)
    grid = np.arange(570.0, 630.0, 0.48)
    # pixels really sample 0.3 nm redward of what they report: the line
    # appears 0.3 nm blueward in reported wavelength
    # residual-free synthetic data leave curve_fit's covariance undefined, so
    # add a little seeded noise (1e-4) and give the fit that error
    e = np.full(grid.size, 1e-4)
    obs = twin.observe(lam, y, grid, 2.5, true_grid=grid + 0.3, noise=1e-4, seed=1)
    assert srf.fit_line(grid, obs, 600.0, err=e)['mu'] == pytest.approx(599.7, abs=0.01)
    a = twin.observe(lam, y, grid, 2.5, noise=0.01, seed=3)
    b = twin.observe(lam, y, grid, 2.5, noise=0.01, seed=3)
    np.testing.assert_array_equal(a, b)
    rel = a / twin.observe(lam, y, grid, 2.5) - 1
    assert 0.007 < rel.std() < 0.013
    # dfwhm widens the observed line
    w = srf.fit_line(grid, twin.observe(lam, y, grid, 2.5, dfwhm=0.3, noise=1e-4, seed=2),
                     600.0, err=e)['fwhm']
    assert w == pytest.approx(np.hypot(2.8, srf.FWHM_PER_SIGMA * 0.5), abs=0.02)


def test_grids_and_instrument_srfs():
    g = twin.load_grids()
    assert g['grid_E_122304'].size == 1536 and g['grid_L_122304'].size == 1538
    assert g['grid_L_122305'].size == 1536
    assert np.nanmedian(g['noise_E_scan']) < 0.02
    s = twin.load_instrument_srfs()
    for inst in ('HYPSTAR_122304', 'HYPSTAR_122305'):
        assert set(s[inst]) == {'E', 'Ld'}
    # the narrow-E / E ~ L split of Phase 0 task 8c
    d304 = s['HYPSTAR_122304']['Ld'].fwhm(450) - s['HYPSTAR_122304']['E'].fwhm(450)
    d305 = s['HYPSTAR_122305']['Ld'].fwhm(450) - s['HYPSTAR_122305']['E'].fwhm(450)
    assert d304 > 0.3 and abs(d305) < 0.15


def test_case_table():
    cases = twin.case_table('HYPSTAR_122304')
    names = [c.name for c in cases]
    for prefix in ('i_', 'ii_', 'iii_', 'iv_shift', 'iv_stretch', 'v_wrongsrf', 'vi_wrongemod',
                   'vii_'):
        assert any(n.startswith(prefix) for n in names), prefix
    c = {k.name: k for k in cases}
    assert c['i_offset'].srf_E is c['i_offset'].srf_L                 # equal SRFs
    assert np.array_equal(c['ii_mismatch'].grid_E, c['ii_mismatch'].grid_L)  # no offset
    st = c['iv_stretch+0.1']
    d = st.true_grid_E - st.grid_E
    assert d[np.argmin(np.abs(st.grid_E - 870))] == pytest.approx(0.1, abs=0.002)
    assert d[np.argmin(np.abs(st.grid_E - 390))] == pytest.approx(-0.1, abs=0.002)
    assert c['v_wrongsrf_E+0.3'].corr_dfwhm_E == 0.3
    assert c['vii_noise'].noise
    assert len(twin.case_table('HYPSTAR_122305')) == len(cases)

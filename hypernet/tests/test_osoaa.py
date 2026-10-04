"""Tests for the OSOAA wrapper (:mod:`hypernet.rt.osoaa`).

- Parsers: on fixtures from one real run (``hypernet/tests/data/osoaa/``:
  550 nm, SZA 40, RAA 90, AOT(550) 0.1, Chl 1).
- Live: ``@needs_osoaa`` runs OSOAA once at 550 nm; it skips when the
  executable is missing (``$OSOAA_ROOT`` or the fork's default path).

Run from the repository root::

    pytest -q hypernet/tests/test_osoaa.py
"""

import os

import numpy as np
import pytest

from hypernet.rt import osoaa

DATA = os.path.join(os.path.dirname(__file__), 'data', 'osoaa')
SZA = 40.0


def _have_osoaa():
    try:
        osoaa.osoaa_root()
        return True
    except FileNotFoundError:
        return False


needs_osoaa = pytest.mark.skipif(not _have_osoaa(), reason='requires the OSOAA executable')


@pytest.fixture(scope='module')
def fx():
    return dict(flux=osoaa.parse_flux(os.path.join(DATA, 'Flux.txt')),
                up=osoaa.parse_lum_advanced(os.path.join(DATA, 'LUM_Advanced_Up.txt')),
                down=osoaa.parse_lum_advanced_down(os.path.join(DATA, 'LUM_Advanced_Down.txt')),
                vsvza=osoaa.parse_lum_vsvza(os.path.join(DATA, 'LUM_vsVZA.txt')))


def test_parse_flux(fx):
    f = fx['flux']
    assert f['level'][0] == 0 and f['z'][0] > 0
    # TOA: the direct beam is pi cos(SZA) for an extraterrestrial irradiance of pi
    assert f['Ed_total'][0] == pytest.approx(np.pi * np.cos(np.radians(SZA)), rel=1e-4)
    np.testing.assert_allclose(f['Ed_direct'] + f['Ed_diffuse'], f['Ed_total'], rtol=2e-5)
    assert np.all(f['Ed_total'] > 0)
    assert np.all(np.diff(f['Ed_total'][27:]) <= 0)        # Ed falls with depth


def test_parse_lum_advanced(fx):
    up, dn = fx['up'], fx['down']
    assert up['direction'] == 'up' and dn['direction'] == 'down'
    assert (dn['level_0plus'], dn['level_0minus']) == (26, 27)
    assert dn['phi_pos'] == pytest.approx(90.0) and dn['phi_neg'] == pytest.approx(270.0)
    assert set(np.unique(dn['level'])) == {0, 25, 26, 27, 28}    # the trimmed fixture
    # no diffuse downward radiance enters at the top of the atmosphere
    assert np.all(dn['I'][dn['level'] == 0] == 0)
    with pytest.raises(ValueError):
        osoaa.parse_lum_advanced_down(os.path.join(DATA, 'LUM_Advanced_Up.txt'))


def test_sky_and_water_radiances(fx):
    f, up, dn = fx['flux'], fx['up'], fx['down']
    Ed0p = f['Ed_total'][f['level'] == dn['level_0plus']][0]
    # Ld at 0+, HYPSTAR sky view (40 deg from zenith, RAA 90): L/Ed = I/F
    ld_ed = osoaa.radiance_at(dn, '0+', 40.0) / Ed0p
    assert 0.01 < ld_ed < 0.1          # clear sky at 550 nm (the 750 nm index is ~0.02)
    # Lu just below the surface, at the refracted water-view angle
    th_w = np.degrees(np.arcsin(np.sin(np.radians(40.0)) / 1.34))
    lu = osoaa.radiance_at(up, '0-', th_w)
    assert 0 < lu / Ed0p < 0.05
    # interpolation returns the grid value at a grid node
    m = (dn['level'] == 26) & (dn['vza'] > 0)
    assert osoaa.radiance_at(dn, 26, dn['vza'][m][3]) == pytest.approx(dn['I'][m][3])


def test_parse_lum_vsvza(fx):
    v, f = fx['vsvza'], fx['flux']
    assert '0+' in v['level_label']
    assert v['vza'].min() < 0 < v['vza'].max()
    # REFL = pi L / Ed(0+) = pi I / F(0+)
    Ed0p = f['Ed_total'][f['level'] == 26][0]
    np.testing.assert_allclose(v['refl'], np.pi * v['I'] / Ed0p, rtol=2e-3)


def test_default_params_fix_aerosol_reference(tmp_path):
    p = osoaa.default_params(700.0, 30.0, str(tmp_path))
    assert p['AER.Waref'] == 0.55                       # fixed, not the run wavelength
    assert p['OSOAA.Wa'] == pytest.approx(0.7)
    assert p['OSOAA.View.Level'] == 3
    for k in ('MIE_AER', 'MIE_HYD', 'SURF'):
        assert os.path.isdir(tmp_path / k)


@needs_osoaa
def test_live_run_550(tmp_path):
    p = osoaa.default_params(550.0, SZA, str(tmp_path / 'run'), mie_dir=str(tmp_path / 'db'))
    r = osoaa.run(p)
    assert r['seconds'] > 0
    assert r['flux']['Ed_total'][0] == pytest.approx(np.pi * np.cos(np.radians(SZA)), rel=1e-4)
    for k in ('up', 'down', 'vsvza'):
        assert k in r
    # same physics as the fixture (same inputs)
    ref = osoaa.parse_flux(os.path.join(DATA, 'Flux.txt'))
    np.testing.assert_allclose(r['flux']['Ed_total'], ref['Ed_total'], rtol=1e-4)

"""Tests for the L1A/L1C/L2A readers (:mod:`hypernet.whn_l1a`).

Two tiers:

- **Tier 1** -- filename parsing and time interpolation, run everywhere.
- **Tier 2** -- ``@needs_wavelengths`` tests read Kevin's VEIT sample
  (SEQ20260604T084543) and skip when ``$OS_COLOR/WATERHYPERNET/Wavelengths``
  or the sample is not there.

Run from the repository root::

    pytest -q hypernet/tests/test_whn_l1a.py
"""

import numpy as np
import pytest

from hypernet import whn_l1a as wl

VEIT = dict(site='VEIT', seq_time='20260604T0845')


def _veit_files():
    try:
        files = wl.sequence_files(**VEIT)
    except Exception:
        return None
    return files if all(p in files for p in wl.PRODUCTS) else None


needs_wavelengths = pytest.mark.skipif(
    _veit_files() is None,
    reason='requires the VEIT sample in $OS_COLOR/WATERHYPERNET/Wavelengths')


# --- Tier 1 ------------------------------------------------------------------

def test_parse_name():
    a = wl.parse_name('/x/HYPERNETS_W_VEIT_L1A_IRR_20260604T0845_20260828T1547_v2.1.nc')
    assert a == dict(network='W', site='VEIT', product='L1A_IRR',
                     seq_time='20260604T0845', proc_time='20260828T1547',
                     azimuth=None, version='v2.1')
    c = wl.parse_name('HYPERNETS_W_VEIT_L1C_ALL_20260604T0845_20260828T1547_090_v2.1.nc')
    assert c['product'] == 'L1C_ALL' and c['azimuth'] == '090'
    assert wl.parse_name('plot_irradiance_HYPERNETS_W_VEIT_L1B_IRR.png') is None


def test_interp_in_time():
    wave = np.arange(5)
    s1, s2 = np.full(5, 10.0), np.full(5, 20.0)
    scans = np.stack([s1, s1, s2, s2, s2], axis=1)
    sid = np.array([4, 4, 10, 10, 10])
    t = np.array([100, 100, 200, 200, 200])
    np.testing.assert_allclose(wl.interp_in_time(scans, sid, t, 125.0), 12.5)
    np.testing.assert_allclose(wl.interp_in_time(scans, sid, t, 50.0), 10.0)   # clamped
    np.testing.assert_allclose(wl.interp_in_time(scans, sid, t, 300.0), 20.0)
    one = wl.interp_in_time(scans[:, :2], sid[:2], t[:2], 999.0)
    np.testing.assert_allclose(one, 10.0)
    means, times, ids = wl.series_means(scans, sid, t)
    assert means.shape == (wave.size, 2) and list(ids) == [4, 10]


def test_wavelengths_root_explicit(tmp_path):
    assert wl.wavelengths_root(str(tmp_path)) == str(tmp_path)


# --- Tier 2: the VEIT sample -------------------------------------------------

@pytest.fixture(scope='module')
def veit():
    f = _veit_files()
    return dict(irr=wl.load_l1a_irr(f['L1A_IRR']), rad=wl.load_l1a_rad(f['L1A_RAD']),
                l1c=wl.load_l1c(f['L1C_ALL']), l2a=wl.load_l2a(f['L2A_REF']))


@needs_wavelengths
def test_veit_grids_and_scans(veit):
    irr, rad, l1c, l2a = veit['irr'], veit['rad'], veit['l1c'], veit['l2a']
    assert irr['wave'].size == 1536 and irr['scans'].shape == (1536, 6)
    assert rad['wave'].size == 1538 and rad['n_scans'] == 12
    assert rad['Lu']['n_scans'] == 6 and rad['Ld']['n_scans'] == 6
    assert l1c['n'] == 6 and l1c['irradiance'].shape == (1538, 6)
    assert l2a['n'] == 1 and l2a['reflectance'].shape == (1538, 1)
    np.testing.assert_array_equal(l1c['wave'], rad['wave'])
    assert irr['meta']['system_id'] == 'HYPSTAR_122304'


@needs_wavelengths
def test_veit_vza_split(veit):
    rad, l1c = veit['rad'], veit['l1c']
    assert np.all(rad['Lu']['vza'] < 90) and np.all(rad['Ld']['vza'] >= 90)
    np.testing.assert_allclose(l1c['vza'], rad['Lu']['vza'])   # L1C is the water view
    assert np.all(veit['irr']['vza'] > 170)                      # E looks up
    # the sky is much brighter than the water view
    j = np.argmin(np.abs(rad['wave'] - 450))
    assert rad['Ld']['mean'][j] > 5 * rad['Lu']['mean'][j]


@needs_wavelengths
def test_veit_l1c_irradiance_is_linear_interp(veit):
    irr, l1c = veit['irr'], veit['l1c']
    w = l1c['wave']
    E_lin = np.interp(w, irr['wave'], irr['mean'])
    rel = np.abs(l1c['irradiance'] - E_lin[:, None]) / E_lin[:, None]
    # 1e-4 where the lines are; the departures grow at low-signal pixels
    # (grid ends, the 935 nm H2O band; wiggles/phase0a_l1c_consistency.py [8])
    assert np.nanmax(rel[(w > 400) & (w < 900)]) < 1e-4
    assert np.nanmax(rel[(w > 380) & (w < 1000)]) < 1e-3


@needs_wavelengths
def test_veit_check_against_l1c(veit):
    d = wl.check_against_l1c(veit['irr'], veit['rad'], veit['l1c'])
    assert d['E'] < 1e-4
    assert d['Lu'] == 0.0
    assert d['Ld'] < 1e-3


@needs_wavelengths
def test_veit_l2a_matches_l1c(veit):
    l1c, l2a = veit['l1c'], veit['l2a']
    np.testing.assert_allclose(l2a['sza'], np.mean(l1c['sza']), atol=0.1)
    band = (l2a['wave'] > 400) & (l2a['wave'] < 700)
    rho = l2a['reflectance'][band, 0]
    assert np.all(np.isfinite(rho)) and 0 < np.median(rho) < 0.1

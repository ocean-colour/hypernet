"""Tests for the WATERHYPERNET exploration layer (:mod:`hypernet.whn_explore`).

Two tiers:

- **Tier 1** -- data-independent tests of the pure logic (product selection,
  apportionment, shape normalisation) run everywhere.
- **Tier 2** -- ``@needs_whn`` tests touch the real archive and skip
  automatically when ``$OS_COLOR/WATERHYPERNET/RELEASE_2`` is not mounted, since
  the 9.8 GB release is not bundled with anything.

Run from the repository root::

    pytest -q hypernet/tests/test_whn_explore.py
"""

import os

import numpy as np
import pandas as pd
import pytest

from hypernet import whn_explore as wx


def _whn_available():
    """True if the WATERHYPERNET RELEASE_2 tree is reachable."""
    try:
        wx.whn_root()
        return True
    except Exception:
        return False


needs_whn = pytest.mark.skipif(
    not _whn_available(),
    reason='requires the WATERHYPERNET RELEASE_2 tree (unbundled, 9.8 GB)')


# --- Tier 1: data-independent ------------------------------------------------

def test_system_of():
    """Site-code suffix decides the measuring system."""
    assert wx.system_of('BEFR_H') == 'HYPSTAR'
    assert wx.system_of('O1BE_P') == 'PANTHYR'


def test_product_for_default_sites():
    """Non-turbid sites use the SimSpec-corrected product, per system spelling."""
    assert wx.product_for('BEFR_H') == ('reflectance', 'std_reflectance')
    assert wx.product_for('VEIT_P') == ('reflectance', 'reflectance_std')


def test_product_for_turbid_sites():
    """The three turbid sites fall back to ``reflectance_nosc``."""
    assert wx.NOSC_SITES == ('LPAR_H', 'MAFR_H', 'O1BE_P')
    assert wx.product_for('LPAR_H') == ('reflectance_nosc', 'std_reflectance_nosc')
    assert wx.product_for('MAFR_H') == ('reflectance_nosc', 'std_reflectance_nosc')
    assert wx.product_for('O1BE_P') == ('reflectance_nosc', 'reflectance_nosc_std')


def test_largest_remainder_totals_and_floor():
    """Apportionment hits the target exactly and never drops a cluster."""
    counts = pd.Series({0: 500, 1: 300, 2: 190, 3: 10})
    quota = wx._largest_remainder(counts, 100)
    assert sum(quota.values()) == 100
    assert set(quota) == set(counts.index)
    assert min(quota.values()) >= 1                 # rare cluster survives
    assert quota[0] > quota[3]                      # tracks occupancy


def test_largest_remainder_more_clusters_than_slots():
    """Asking for fewer spectra than clusters still returns exactly that many."""
    counts = pd.Series({0: 5, 1: 4, 2: 3})
    quota = wx._largest_remainder(counts, 2)
    assert sum(quota.values()) == 2


def test_shape_matrix_is_unit_norm_and_drops_bad_rows():
    """Shapes are L2-normalised; rows with NaN in the cluster range are dropped."""
    n = wx.ANALYSIS_WAVE.size
    grids = np.vstack([np.linspace(0.001, 0.01, n),        # fine
                       np.full(n, np.nan),                 # all NaN -> dropped
                       np.full(n, 0.005)])                 # fine
    shapes, ok = wx.shape_matrix(grids)
    assert ok.tolist() == [True, False, True]
    assert shapes.shape[0] == 2
    np.testing.assert_allclose(np.linalg.norm(shapes, axis=1), 1.0, rtol=1e-6)


def test_shape_matrix_handles_negative_reflectance():
    """Negative bands are kept -- L2 normalisation must not blow up on them."""
    n = wx.ANALYSIS_WAVE.size
    spec = np.full(n, 0.004)
    spec[wx.ANALYSIS_WAVE > 750] = -0.0005           # realistic NIR negatives
    shapes, ok = wx.shape_matrix(spec[None, :])
    assert ok.all()
    assert np.isfinite(shapes).all()


def test_analysis_grid_covers_both_systems():
    """The display grid stays inside the range both instruments actually measure."""
    assert wx.ANALYSIS_WAVE[0] >= 350.0
    assert wx.ANALYSIS_WAVE[-1] <= 945.0
    lo, hi = wx.CLUSTER_RANGE
    assert lo >= 400.0 and hi <= 900.0               # the notes' reliable range


# --- Tier 2: needs the archive -----------------------------------------------

@needs_whn
def test_whn_root_exists():
    """The archive resolves and contains the expected site directories."""
    root = wx.whn_root()
    sites = {d for d in os.listdir(root) if not d.startswith('0_')}
    assert {'BEFR_H', 'VEIT_H', 'O1BE_P', 'WRUK_H'} <= sites


@needs_whn
def test_load_spectrum_hypstar():
    """A HYPSTAR file reads as Rrs on its native grid with finite sigma."""
    root = wx.whn_root()
    import glob
    path = sorted(glob.glob(os.path.join(root, 'BEFR_H', '*', '*', '*', '*.nc')))[0]
    spec = wx.load_spectrum(path, 'BEFR_H')

    assert spec['wave'].size > 1500                  # ~1539 bands
    assert np.all(np.diff(spec['wave']) > 0)         # ascending
    assert spec['Rrs'].shape == spec['wave'].shape
    assert np.isfinite(spec['sigma']).all()          # std is always present
    assert spec['serial'].startswith('HYPSTAR_')
    assert 40 < spec['lat'] < 50                     # Berre lagoon, France


@needs_whn
def test_load_spectrum_panthyr_fixed_grid():
    """PANTHYR files share one fixed 237-band grid and read despite uint8 flags."""
    root = wx.whn_root()
    import glob
    path = sorted(glob.glob(os.path.join(root, 'O1BE_P', '*', '*', '*', '*.nc')))[0]
    spec = wx.load_spectrum(path, 'O1BE_P')

    assert spec['wave'].size == 237
    assert spec['wave'][0] == pytest.approx(355.0)
    assert spec['wave'][-1] == pytest.approx(945.0)
    assert np.isfinite(spec['Rrs']).any()


@needs_whn
def test_rrs_is_rhow_over_pi():
    """The stored quantity is Rrs, i.e. the file's rho_w divided by pi.

    Guards the single most consequential convention in this module: both
    archive products are water-leaving reflectance, not remote-sensing
    reflectance.
    """
    import glob
    import netCDF4

    root = wx.whn_root()
    path = sorted(glob.glob(os.path.join(root, 'BEFR_H', '*', '*', '*', '*.nc')))[0]
    spec = wx.load_spectrum(path, 'BEFR_H')

    ds = netCDF4.Dataset(path)
    raw = wx._filled(ds.variables['reflectance'])
    ds.close()

    finite = np.isfinite(raw) & np.isfinite(spec['Rrs'])
    np.testing.assert_allclose(spec['Rrs'][finite], raw[finite] / np.pi,
                               rtol=1e-6)


@needs_whn
def test_load_spectrum_analysis_grid():
    """The display grid is optional and matches :data:`ANALYSIS_WAVE`."""
    import glob
    root = wx.whn_root()
    path = sorted(glob.glob(os.path.join(root, 'VEIT_H', '*', '*', '*', '*.nc')))[0]

    spec = wx.load_spectrum(path, 'VEIT_H')
    assert spec['Rrs_grid'].shape == wx.ANALYSIS_WAVE.shape

    bare = wx.load_spectrum(path, 'VEIT_H', analysis_grid=False)
    assert bare['Rrs_grid'] is None


@needs_whn
@pytest.mark.parametrize('site', ['O1BE_P', 'THFR_H'])
def test_quality_flag_survives_the_panthyr_fill_value(site):
    """A passing ``quality_flag`` of 0 must not be read as missing.

    Regression test. PANTHYR declares ``_FillValue = 0`` on a bitmask whose zero
    value means *no flags set* -- a pass -- so a masked read turns every passing
    PANTHYR measurement into NaN, and anyone filtering ``quality_flag == 0``
    drops the whole PANTHYR half of the archive. HYPSTAR declares no
    ``_FillValue`` here and is checked alongside to show the fix changes nothing
    for it.
    """
    import glob
    import netCDF4

    root = wx.whn_root()
    path = sorted(glob.glob(os.path.join(root, site, '*', '*', '*', '*.nc')))[0]

    flag = wx.load_spectrum(path, site, analysis_grid=False)['quality_flag']

    ds = netCDF4.Dataset(path)
    ds.set_auto_mask(False)
    raw = float(np.asarray(ds.variables['quality_flag'][:]).ravel()[0])
    ds.close()

    assert np.isfinite(flag), f'{site}: quality_flag read as NaN, not {raw}'
    assert flag == raw


@needs_whn
def test_build_index_small():
    """The filename index parses both grammars and yields sane columns."""
    df = wx.build_index(save=False)
    assert len(df) > 50000
    assert set(df.columns) == {'site', 'system', 'datetime', 'azimuth', 'path'}
    assert df.datetime.notna().all()
    assert df.site.nunique() == 11
    assert set(df.system.unique()) == {'HYPSTAR', 'PANTHYR'}


# --- Similarity-Spectrum over-subtraction check ------------------------------

from hypernet import whn_simspec_check as sc        # noqa: E402


def test_flag_oversubtraction_logic():
    """Flagged only when the corrected product died and the raw one did not."""
    med_ref = np.array([-0.0002, 0.0050, -0.0003, 0.0000])
    med_nosc = np.array([0.0004, 0.0060, -0.0001, 0.0004])
    assert sc.flag_oversubtraction(med_ref, med_nosc).tolist() == \
        [True, False, False, True]      # a zero corrected median counts


def test_dark_sites_are_hypstar_only():
    """The check targets the three darkest sites, all HYPSTAR."""
    assert sc.DARK_SITES == ('THFR_H', 'BEFR_H', 'WRUK_H')
    assert all(s.endswith('_H') for s in sc.DARK_SITES)
    assert sc.TEST_RANGE == (440.0, 600.0)


@needs_whn
def test_check_file_on_a_known_oversubtracted_spectrum():
    """A file identified by the full scan still reads as over-subtracted.

    Regression guard on the whole chain: fill handling, the 440-600 nm window
    and the sign test.
    """
    path = os.path.join(
        wx.whn_root(), 'THFR_H', '2025', '11', '28',
        'HYPERNETS_W_THFR_L2B_REF_20251128T0730_20260528T2005_270_v2.1.nc')
    if not os.path.exists(path):
        pytest.skip('reference spectrum not present in this copy of the archive')

    stats = sc.check_file(path)
    assert stats is not None
    assert stats['med_ref'] <= 0 < stats['med_nosc']
    assert stats['quality_flag'] == 0          # it passed QC
    assert stats['n_valid'] == stats['n_total']
    assert sc.flag_oversubtraction(stats['med_ref'], stats['med_nosc'])


@needs_whn
def test_check_file_on_a_normal_spectrum():
    """A bright-site spectrum is not flagged."""
    import glob
    path = sorted(glob.glob(os.path.join(
        wx.whn_root(), 'GAIT_H', '*', '*', '*', '*.nc')))[0]
    stats = sc.check_file(path)
    assert stats is not None
    assert stats['med_ref'] > 0
    assert not sc.flag_oversubtraction(stats['med_ref'], stats['med_nosc'])

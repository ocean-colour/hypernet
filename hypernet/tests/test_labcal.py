"""Tests for :mod:`hypernet.labcal` and the committed lab tables
(lab data prompts, task 2).

- Tier 1: ``hypernet/data/lab_wl_coefs.csv`` and ``lab_srf.csv`` ship with
  the package (round trip, contents, crop rule, pixel counts).
- Tier 2: the lab tree (``needs_lab``), the VEIT sample (``needs_veit``) and
  the Release 2 archive (``needs_whn``) skip themselves when absent.
"""
import os

import numpy as np
import pandas as pd
import pytest

from hypernet import labcal
from hypernet import whn_l1a as wl

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REQUEST = os.path.join(REPO, 'docs', 'wiggles_data_request.csv')
EPOCHS = {'120242': ['2020_09', '2022_05', '2022_06', '2024_05'],
          '121222': ['2021_07', '2021_07_tst'], '121231': ['2021_07'],
          '122302': ['2023_01', '2024_09'], '122303': ['2023_01'],
          '122304': ['2023_01', '2024_09'], '122305': ['2023_01']}
# processor grids: (first px, n px); the Release 2 L lengths, and E by the rule
NPX = {('120242', 'L'): (411, 1540), ('120242', 'E'): (412, 1539),
       ('121222', 'L'): (412, 1538), ('121222', 'E'): (412, 1538),
       ('121231', 'L'): (410, 1538), ('121231', 'E'): (410, 1538),
       ('122302', 'L'): (412, 1539), ('122302', 'E'): (413, 1538),
       ('122303', 'L'): (413, 1538), ('122303', 'E'): (414, 1537),
       ('122304', 'L'): (411, 1538), ('122304', 'E'): (413, 1536),
       ('122305', 'L'): (412, 1536), ('122305', 'E'): (413, 1535)}


def _veit():
    try:
        f = wl.sequence_files(site='VEIT', seq_time='20260604T0845')
        return f if 'L1A_IRR' in f else None
    except Exception:
        return None


def _have_whn():
    try:
        from hypernet.whn_explore import whn_root
        return os.path.isdir(whn_root())
    except Exception:
        return False


needs_lab = pytest.mark.skipif(not labcal.lab_available(), reason='requires the lab tree')
needs_veit = pytest.mark.skipif(_veit() is None, reason='requires the VEIT sample')
needs_whn = pytest.mark.skipif(not _have_whn(), reason='requires the Release 2 archive')


# ---------------------------------------------------------------- Tier 1
def test_coef_table_contents():
    t = labcal.load_wl_coefs()
    assert list(t.columns) == labcal.COEF_COLUMNS
    assert len(t) == 26 and not t.duplicated(['sn', 'epoch', 'channel']).any()
    assert {sn: sorted(set(g['epoch'])) for sn, g in t.groupby('sn')} == EPOCHS
    assert (t['lines_source'] == 'assumed_all_deltawl').all()
    assert (t['n_lines'] == t['lines'].str.split().str.len()).all()
    # exactly one processor epoch per instrument and channel, as PROCESSOR_EPOCH says
    p = t[t['processor']]
    assert len(p) == 14
    for sn, d in labcal.PROCESSOR_EPOCH.items():
        for ch in labcal.CHANNELS:
            r = p[(p['sn'] == sn) & (p['channel'] == ch)]
            assert r['epoch'].tolist() == [d['epoch']]
            assert r['processor_status'].iloc[0] == d[ch]
    assert (t.loc[~t['processor'], 'processor_status'] == 'no').all()
    assert not t.loc[t['epoch'] == '2021_07_tst', 'processor'].any()
    assert set(labcal.PROCESSOR_EPOCH['122304'].values()) == {'2023_01', 'verified'}


def test_coef_table_round_trip(tmp_path):
    t = labcal.load_wl_coefs()
    p = labcal.write_wl_coefs(t, str(tmp_path / 'c.csv'))
    u = labcal.load_wl_coefs(p)
    pd.testing.assert_frame_equal(t, u)
    c = ['c0', 'c1', 'c2', 'c3', 'c4']
    assert np.array_equal(t[c].values, u[c].values)          # bit for bit


def test_crop_rule_and_pixel_counts():
    for (sn, ch), (start, n) in NPX.items():
        lam, px = labcal.grid_for('HYPSTAR_' + sn, ch, return_px=True)
        assert (px[0], len(px)) == (start, n), (sn, ch)
        assert np.all(np.diff(px) == 1) and np.all(np.diff(lam) > 0)
        assert lam[0] >= labcal.LAM_LO and lam[-1] < labcal.LAM_HI
        c = labcal.coefs(sn, ch)
        before = labcal.poly(c, px[0] - 1)
        assert before < labcal.LAM_LO                      # nothing left out at the blue end
        if sn not in labcal.CROP_EPOCH:
            assert labcal.poly(c, px[-1] + 1) >= labcal.LAM_HI
    # Ld and Lu share the L grid; the 120242 exception is the 2020_09 window
    assert np.array_equal(labcal.grid_for('122304', 'Ld'), labcal.grid_for('122304', 'Lu'))
    s, e = labcal.crop_window(labcal.coefs('120242', 'L', '2020_09'))
    assert (s, e) == (411, 1950)
    s, e = labcal.crop_window(labcal.coefs('120242', 'L', '2024_05'))     # the plain rule
    assert e - s + 1 == 1541
    assert len(labcal.grid_for('120242', 'L', epoch='2022_06')) == 1541   # non-processor epoch
    with pytest.raises(ValueError):
        labcal.grid_for('122304', 'X')
    with pytest.raises(KeyError):
        labcal.grid_for('122304', 'E', epoch='1999_01')


def test_lab_srf_table():
    s = labcal.load_lab_srf()
    assert list(s.columns) == labcal.SRF_COLUMNS
    assert not s.duplicated(['sn', 'epoch', 'channel', 'lam_ref']).any()
    t = labcal.load_wl_coefs()
    n = t[t['channel'] == 'L'].set_index(['sn', 'epoch'])['n_lines']
    assert (s.groupby(['sn', 'epoch']).size() == 2 * n).all()
    assert s['mu_px'].notna().all()
    assert set(s.loc[s['blend'], 'lam_ref']) == set(labcal.BLENDS)
    assert (s['low_snr'] == ~(s['snr'] >= labcal.SNR_MIN)).all()
    assert set(s.loc[s['outside_grid'], 'lam_ref']) == {296.728}
    # the two instrument classes: lab Delta FWHM at Hg 546 on the processor epochs
    p = s[s['lam_ref'] == 546.0735].pivot_table(index=['sn', 'epoch'], columns='channel',
                                                values='lab_fwhm')
    d = (p['L'] - p['E'])
    for sn, e in labcal.PROCESSOR_EPOCH.items():
        narrow = sn in ('120242', '122302', '122304')
        assert (d[(sn, e['epoch'])] > 0.3) if narrow else (abs(d[(sn, e['epoch'])]) < 0.1), sn
    # the polynomial residual at the bright lines is the stated sigma's size
    b = s[(s['snr'] > 100) & ~s['outside_grid']]
    assert b['resid'].abs().max() < 0.2


# ---------------------------------------------------------------- Tier 2
@needs_lab
def test_readers_and_committed_coefficients():
    t = labcal.load_wl_coefs()
    for sn, ep, d in labcal.epochs():
        e = labcal.read_epoch(d)
        assert e['lamp'][0].size == labcal.NPIX
        assert all(np.isclose(L.max(), 1) and np.isclose(E.max(), 1) for L, E in e['lamp'][1])
        assert len(e['lamps']) == len(e['lamp'][1]) and e['lamps'][0] == 'HgAr'
        for ch in labcal.CHANNELS:
            assert np.array_equal(labcal.coefs(sn, ch, ep, t), e['coefs'][ch])
    # the _tst coefficients differ from 2021_07's despite the shared file name
    a, b = (labcal.read_coefs(labcal.epoch_dir('121222', x)) for x in ('2021_07', '2021_07_tst'))
    assert a['coef_file'] == b['coef_file'] and not np.allclose(a['L'], b['L'])


@needs_veit
def test_grid_for_reproduces_veit_l1a():
    f = _veit()
    wE = wl.load_l1a_irr(f['L1A_IRR'])['wave']
    wL = wl.load_l1a_rad(f['L1A_RAD'])['wave']
    gE, pE = labcal.grid_for('HYPSTAR_122304', 'E', return_px=True)
    gL, pL = labcal.grid_for('HYPSTAR_122304', 'L', return_px=True)
    assert gE.shape == wE.shape and np.abs(gE - wE).max() <= 1e-4
    assert gL.shape == wL.shape and np.abs(gL - wL).max() <= 1e-4
    # guard: the 2024_09 lab polynomial is not what the processor uses
    for ch, px, w in (('E', pE, wE), ('L', pL, wL)):
        g = labcal.poly(labcal.coefs('122304', ch, '2024_09'), px)
        assert np.sqrt(np.mean((g - w) ** 2)) > 0.02, ch


@needs_whn
def test_grid_for_reproduces_release2_l_grids():
    req = pd.read_csv(REQUEST, dtype=str)
    seen = set()
    for _, r in req.iterrows():
        sn = r['instrument'].split('_')[1]
        if sn in seen:
            continue
        p = wl.release2_path(r['site'], r['sequence_time'], r['file'])
        if not os.path.exists(p):
            continue
        w = wl.load_l2b(p)['wave']
        g = labcal.grid_for(sn, 'Ld')
        assert g.shape == w.shape, sn
        assert np.abs(g - w).max() <= 1e-4, sn
        seen.add(sn)
    assert seen == set(labcal.INSTRUMENTS)


def test_npz_grids_equal_grid_for():
    """hypstar_grids.npz holds every instrument's processor grids, exactly
    labcal.grid_for (Setup Q11), so the two cannot drift."""
    from hypernet import twin
    g = twin.load_grids()
    for sn in labcal.INSTRUMENTS:
        for ch in labcal.CHANNELS:
            k = 'grid_%s_%s' % (ch, sn)
            assert k in g, k
            ref = labcal.grid_for(sn, ch)
            assert g[k].shape == ref.shape and np.abs(g[k] - ref).max() <= 1e-9, k
    assert {'noise_E_scan', 'noise_Ld_scan', 'noise_Lu_scan'} <= set(g)
    assert g['noise_E_scan'].size == g['grid_E_122304'].size


def test_srf_shape_for_classes():
    """Setup Q7: super-Gaussian E on 120242, 122302, 122304 only; L Gaussian;
    p from the physical epoch (Q8)."""
    for sn in labcal.INSTRUMENTS:
        e = labcal.srf_shape_for(sn, 'E', '2026-01-01')
        assert e['shape'] == ('supergauss' if sn in ('120242', '122302', '122304') else 'gauss')
        assert labcal.srf_shape_for(sn, 'Ld', '2026-01-01')['shape'] == 'gauss'
        if e['shape'] == 'supergauss':
            assert 2.4 < e['shape_p'] < 3.1
    a = labcal.srf_shape_for('122304', 'E', '2023-02-01')
    b = labcal.srf_shape_for('122304', 'E', '2024-11-12')
    assert (a['lab_epoch'], b['lab_epoch']) == ('2023_01', '2024_09')
    assert a['shape_p'] != b['shape_p']

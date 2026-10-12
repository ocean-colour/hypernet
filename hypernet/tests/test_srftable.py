"""Tests for :mod:`hypernet.srftable` and the committed ``hypernet/data/srf_table.csv``
(Phase 2 task 6).  Tier 1: the CSV ships with the package."""
import numpy as np
import pandas as pd
import pytest

from hypernet import edinterp, srfcal, srftable

INSTRUMENTS = ('HYPSTAR_120242', 'HYPSTAR_121222', 'HYPSTAR_121231', 'HYPSTAR_122302',
               'HYPSTAR_122303', 'HYPSTAR_122304', 'HYPSTAR_122305')


def test_table_contents():
    t = srftable.load_srf_table()
    need = {'instrument', 'channel', 'cal_period', 'valid_from', 'source', 'provisional',
            'c0', 'c1', 'c2', 'cov00', 'cov11', 'cov22', 'offset', 'offset_err',
            'n_sequences', 'delta_fwhm', 'method', 'scale_shift', 'scale_stretch'}
    assert need <= set(t.columns)
    assert set(t['instrument']) == set(INSTRUMENTS) and set(t['channel']) == {'E', 'Ld'}
    assert t['provisional'].all()                     # until the L1A delivery is fitted
    assert set(t['source']) <= set(srftable.SOURCE_RANK)
    assert t['fwhm_550'].between(1.8, 3.2).all()
    assert t['method'].isin(edinterp.METHODS).all()
    # the G1 rule, row by row, and one method per E/Ld pair
    rule = np.where(t['delta_fwhm'] >= edinterp.SRF_THRESHOLD - 1e-9, 'srf', 'linear')
    assert (t['method'] == rule).all()
    assert (t.groupby(['instrument', 'cal_period', 'source'])['method'].nunique() == 1).all()
    # valid_from is the channel's own calibration date
    e, ld = t[t['channel'] == 'E'], t[t['channel'] == 'Ld']
    assert (e['valid_from'] == pd.to_datetime(e['cal_date_irr'])).all()
    assert (ld['valid_from'] == pd.to_datetime(ld['cal_date_rad'])).all()


def test_covariances_are_valid():
    t = srftable.load_srf_table()
    for _, r in t.iterrows():
        m = srftable.model_from_row(r)
        assert np.all(np.linalg.eigvalsh(m.cov) > -1e-12), r['instrument']
        assert m.fwhm(550.0) == pytest.approx(r['fwhm_550'], rel=1e-5)


def test_srf_for_dates_and_preference():
    m = srftable.srf_for('HYPSTAR_122304', 'E', '20260604T0845')
    assert m.meta['source'] == 'veit_l1a'                       # L1A preferred over L2B
    assert m.meta['cal_period'] == '2024-11-08/2024-11-12'
    old = srftable.srf_for('HYPSTAR_122304', 'E', '2024-01-01')
    assert old.meta['cal_period'] == '2023-01-27/2023-01-25'
    # between the RAD (11-08) and IRR (11-12) calibrations: Ld new, E old
    assert srftable.srf_for('HYPSTAR_122304', 'Ld', '2024-11-10').meta['cal_period'] == \
        '2024-11-08/2024-11-12'
    assert srftable.srf_for('HYPSTAR_122304', 'E', '2024-11-10').meta['cal_period'] == \
        '2023-01-27/2023-01-25'
    with pytest.raises(KeyError):
        srftable.srf_for('HYPSTAR_122304', 'E', '2019-01-01')
    with pytest.raises(KeyError):
        srftable.srf_for('HYPSTAR_999999', 'E', '2025-01-01')


def test_method_for_matches_the_rule():
    assert srftable.method_for('HYPSTAR_122304', '2025-06-01') == 'srf'
    assert srftable.method_for('HYPSTAR_122305', '2025-06-01') == 'linear'
    for inst in INSTRUMENTS:
        e = srftable.srf_for(inst, 'E', '2026-01-01')
        ld = srftable.srf_for(inst, 'Ld', '2026-01-01')
        if e.meta['cal_period'] == ld.meta['cal_period']:
            assert edinterp.choose_method(e, ld) == srftable.method_for(inst, '2026-01-01')


def test_scale_for():
    s = srftable.scale_for('HYPSTAR_122302', 'E', '2024-06-01')
    assert s['ok'] and s['n'] >= 10 and abs(s['shift']) < 0.1 and abs(s['stretch']) < 0.15
    assert s['shift_err'] < 0.02
    wav = np.array([390.0, 630.0, 870.0])
    corr = srfcal.correct_wavelengths(wav, s)
    np.testing.assert_allclose(corr - wav, [s['shift'] - s['stretch'], s['shift'],
                                            s['shift'] + s['stretch']])
    use, src = srfcal.scale_or_fallback(dict(ok=False), s)
    assert use is s and src == 'fallback'


def test_g0_switches():
    assert srftable.TABLE_IS_DEFAULT and not srftable.APPLY_CENTROID_OFFSET

"""Tests for the per-sequence SRF glue (:mod:`hypernet.whn_srf`) and the
reference spectrum (:mod:`hypernet.refspec`).

Tier 2 only: ``@needs_veit`` reads Kevin's VEIT sample, ``@needs_hsrs`` the
TSIS-1 HSRS file in ``$OS_COLOR/hypernet/ref``; both skip when absent.  The
VEIT test stands in for task 10's "one delivered sequence fits without
error" until the delivery arrives.

Run from the repository root::

    pytest -q hypernet/tests/test_whn_srf.py
"""

import numpy as np
import pytest

from hypernet import refspec, whn_srf
from hypernet import whn_l1a as wl


def _veit():
    try:
        f = wl.sequence_files(site='VEIT', seq_time='20260604T0845')
        return f if 'L1A_IRR' in f and 'L1A_RAD' in f else None
    except Exception:
        return None


needs_veit = pytest.mark.skipif(_veit() is None, reason='requires the VEIT sample')
needs_hsrs = pytest.mark.skipif(not refspec.hsrs_available(),
                                reason='requires the TSIS-1 HSRS file')


@needs_hsrs
def test_hsrs_checksum_and_frames():
    assert refspec.fetch_hsrs(verify=True) == refspec.hsrs_path()
    wv, fv = refspec.load_hsrs(655, 658, frame='vac')
    wa, fa = refspec.load_hsrs(655, 658, frame='air')
    # H alpha: the minimum sits at 656.46 nm (vacuum) and 656.28 nm (air)
    assert wv[np.argmin(fv)] == pytest.approx(656.46, abs=0.02)
    assert wa[np.argmin(fa)] == pytest.approx(656.28, abs=0.02)


@needs_veit
@needs_hsrs
def test_fit_sequence_veit():
    res = whn_srf.fit_files(_veit(), veil=False)
    assert res['meta']['system_id'] == 'HYPSTAR_122304'
    assert set(res['lines']['channel']) == {'E', 'Ld', 'Lu'}
    methods = {(m.meta['method'], m.channel) for m in res['models']}
    assert methods == {(a, c) for a in ('empirical', 'template') for c in ('E', 'Ld', 'Lu')}
    tm = {m.channel: m for m in res['models'] if m.meta['method'] == 'template'}
    # the template SRF: E ~2.2-2.4 nm, Ld wider by ~0.5 nm in the blue (task 3b)
    assert 2.0 < tm['E'].fwhm(450) < 2.6
    assert 0.3 < tm['Ld'].fwhm(450) - tm['E'].fwhm(450) < 0.8
    assert res['template']['ok'].mean() > 0.8
    tab = whn_srf.models_table(res['models'])
    assert len(tab) == 6 and tab['fwhm_450'].notna().all()

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

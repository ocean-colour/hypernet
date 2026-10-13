"""Tests for :mod:`hypernet.snippets.interpolate_wav_srf` (Phase 2 task 8).

Tier 1: synthetic spectra and obsarray datasets shaped like the processor's.
Tier 2 (``@needs_fields``): the Phase 1 VEIT-like twin scene on the real
122304 grids with the SRF table.
"""
import warnings

import numpy as np
import punpy
import pytest
import scipy.interpolate
import xarray as xr

from hypernet import edinterp, srf
from hypernet.snippets.interpolate_wav_srf import (SRF_INPUT_NAMES, InterpolationWavSRF,
                                                   add_srf_variables)
from hypernet.tests.test_edinterp import GE, GL, _emod, _measured

warnings.filterwarnings('ignore', category=FutureWarning)


def _processor_linear(rad_wavs, irr_wavs, irr):
    """InterpolationWavLinear.meas_function at hypernets_processor 9a12819a
    (context/hypernets_processor/...), verbatim in effect."""
    return scipy.interpolate.interp1d(irr_wavs, irr, axis=0, fill_value='extrapolate')(rad_wavs)


def _models(sig=(0.02, 0.01, 0.005)):
    cov = np.diag(sig) ** 2
    return srf.SRFModel('E', [2.3, 0, 0], cov), srf.SRFModel('Ld', [2.8, 0, 0], cov)


def test_interface_matches_the_processor():
    mf = InterpolationWavSRF(prop=punpy.MCPropagation(10), corr_dims='wavelength',
                             yvariable='irradiance')
    assert InterpolationWavSRF.get_name() == 'InterpolationWavSRF'
    assert mf.get_argument_names() == ['radiance_wavelength', 'wavelength', 'irradiance']
    mf2 = InterpolationWavSRF(prop=punpy.MCPropagation(10), srf_inputs=True)
    assert mf2.get_argument_names()[3:] == SRF_INPUT_NAMES


def test_unwired_setup_falls_back_to_the_processor_linear():
    mf = InterpolationWavSRF(prop=punpy.MCPropagation(10))
    with pytest.warns(UserWarning):
        mf.setup(40.0, 'w', GE)                      # what the processor passes today
    E = np.stack([_measured(_emod(), 2.3)] * 2, axis=1)
    gl = np.r_[GL[0] - 2.0, GL, GL[-1] + 2.0]        # beyond the E grid: extrapolated
    np.testing.assert_array_equal(mf.meas_function(gl, GE, E), _processor_linear(gl, GE, E))
    # punpy's broadcast 2-D wavelength arrays are handled like the processor's
    np.testing.assert_array_equal(mf.meas_function(np.tile(gl[:, None], 2),
                                                   np.tile(GE[:, None], 2), E),
                                  _processor_linear(gl, GE, E))


def test_srf_method_calls_interpolate_ed_to_l():
    em = _emod()
    E = _measured(em, 2.3)
    sE, sL = _models()
    mf = InterpolationWavSRF(prop=punpy.MCPropagation(10))
    mf.setup(40.0, 'w', GE, srf_irr=sE, srf_rad=sL, emod=em)
    assert mf.method == 'srf'                                   # delta FWHM 0.5 nm
    gl = np.r_[GL[0] - 2.0, GL]
    out = mf.meas_function(gl, GE, E)
    ref = edinterp.interpolate_ed_to_l(GE, E, GL, emod=em, srf_irr=sE, srf_rad=sL)
    np.testing.assert_allclose(out[1:], ref, rtol=1e-14)
    assert out[0] == pytest.approx(_processor_linear(gl[:1], GE, E)[0])
    # explicit extra inputs equal to the nominal values give the same answer
    np.testing.assert_allclose(mf.meas_function(gl, GE, E, *[
        np.asarray(sE.coeffs), np.asarray(sL.coeffs), np.zeros(2), np.zeros(2)]), out,
        rtol=1e-12)


def _datasets(E, wE, wL, u_rand, u_sys):
    """Minimal processor-shaped datasets: the L1B irradiance (obsarray
    uncertainties, %) and the L1C radiance grid (``radiance_wavelength``)."""
    ns = E.shape[1]
    irr = xr.Dataset({'irradiance': (('wavelength', 'scan'), E, {'units': 'mW m^-2 nm^-1'})},
                     coords={'wavelength': wE, 'scan': np.arange(ns)})

    def ec(f1, f2):
        return [{'dim': 'wavelength', 'form': f1, 'params': [], 'units': []},
                {'dim': 'scan', 'form': f2, 'params': [], 'units': []}]
    irr.unc['irradiance']['u_rel_random_irradiance'] = (
        ('wavelength', 'scan'), np.full(E.shape, u_rand),
        {'err_corr': ec('random', 'random'), 'pdf_shape': 'gaussian', 'units': '%'})
    irr.unc['irradiance']['u_rel_systematic_indep_irradiance'] = (
        ('wavelength', 'scan'), np.full(E.shape, u_sys),
        {'err_corr': ec('systematic', 'systematic'), 'pdf_shape': 'gaussian', 'units': '%'})
    rad = xr.Dataset(coords={'radiance_wavelength': wL})
    return irr, rad


def test_dataset_route_agrees_with_propagate_ed_to_l():
    """Build #2 Q1: the snippet's dataset route (a), propagated by punpy as the
    processor does, agrees with edinterp.propagate_ed_to_l (b)."""
    em = _emod()
    E1 = _measured(em, 2.3)
    E = np.stack([E1, 1.02 * E1], axis=1)
    sE, sL = _models()
    sc = dict(shift=0.0, stretch=0.0, shift_err=0.01, stretch_err=0.01)
    n = 300
    mf = InterpolationWavSRF(prop=punpy.MCPropagation(n), srf_inputs=True,
                             corr_dims='radiance_wavelength', yvariable='irradiance',
                             ydims=['radiance_wavelength', 'scan'])
    mf.setup(40.0, 'w', GE, srf_irr=sE, srf_rad=sL, scale_irr=sc, scale_rad=sc, emod=em)
    irr, rad = _datasets(E, GE, GL, 0.5, 1.0)
    irr = add_srf_variables(irr, mf)
    np.random.seed(21)
    out = mf.propagate_ds_specific(['random', 'systematic_indep'], rad, irr,
                                   store_unc_percent=True)
    ed_b, info_b = edinterp.propagate_ed_to_l(GE, E, GL, emod=em, srf_irr=sE, srf_rad=sL,
                                              u_rel_random=0.5, u_rel_systematic_indep=1.0,
                                              scale_irr=sc, scale_rad=sc, mcsteps=n, seed=22)
    np.testing.assert_allclose(out['irradiance'].values, ed_b, rtol=1e-10)
    for comp in ('random', 'systematic_indep'):
        a = out['u_rel_%s_irradiance' % comp].values
        b = info_b['u_rel_' + comp]
        r = a / b - 1
        assert abs(np.median(r)) < 0.06 and np.percentile(np.abs(r), 95) < 0.25, \
            (comp, np.median(r))
    # the new (SRF + scale) terms alone (no input systematic): both routes agree,
    # and without the extra inputs there is nothing to propagate
    irr2, _ = _datasets(E, GE, GL, 0.5, 0.0)
    irr2 = add_srf_variables(irr2, mf)
    np.random.seed(23)
    out2 = mf.propagate_ds_specific(['systematic_indep'], rad, irr2, store_unc_percent=True)
    _, info2 = edinterp.propagate_ed_to_l(GE, E, GL, emod=em, srf_irr=sE, srf_rad=sL,
                                          scale_irr=sc, scale_rad=sc, mcsteps=n, seed=24)
    a2 = out2['u_rel_systematic_indep_irradiance'].values
    b2 = info2['u_rel_systematic_indep']
    assert np.median(b2) > 1e-3
    r = a2 / b2 - 1
    assert abs(np.median(r)) < 0.1, np.median(r)
    mf0 = InterpolationWavSRF(prop=punpy.MCPropagation(n), corr_dims='radiance_wavelength',
                              yvariable='irradiance', ydims=['radiance_wavelength', 'scan'])
    mf0.setup(40.0, 'w', GE, srf_irr=sE, srf_rad=sL, scale_irr=sc, scale_rad=sc, emod=em)
    irr0, _ = _datasets(E, GE, GL, 0.5, 0.0)
    out0 = mf0.propagate_ds_specific(['systematic_indep'], rad, irr0, store_unc_percent=True)
    name = 'u_rel_systematic_indep_irradiance'
    if name in out0:                 # punpy omits the variable when nothing propagates
        u0 = out0[name].values
        assert np.all(~np.isfinite(u0) | (u0 == 0))


def _have_fields():
    try:
        from hypernet import twin
        twin.load_fields()
        twin.load_grids()
        return True
    except Exception:
        return False


needs_fields = pytest.mark.skipif(not _have_fields(), reason='requires the Phase 1 fields')


@needs_fields
def test_phase1_scene_through_the_snippet():
    """The VEIT-like twin scene on the real 122304 grids, SRFs from the table:
    Ed_L from the snippet matches the truth (E seen through the L SRF) far
    better than linear interpolation near the lines."""
    from hypernet import emod as emodmod
    from hypernet import twin
    fields = twin.load_fields()
    grids = twin.load_grids()
    e = emodmod.build_emod(40.0)
    scene = twin.compose_scene(e, fields, 'veit_sza40_aot0.10_chl1')
    gE = grids['grid_E_122304']
    gE = gE[(gE > 389) & (gE < 991)]
    gL = grids['grid_L_122304']
    gL = gL[(gL > 392) & (gL < 988)]
    mf = InterpolationWavSRF(prop=punpy.MCPropagation(10))
    zero = dict(shift=0.0, stretch=0.0)                      # the twin has no scale error
    mf.setup(40.0, 'w', gE, system_id='HYPSTAR_122304', date='20260604T0845',
             scale_irr=zero, scale_rad=zero, emod=(scene['lam'], scene['Ed']))
    assert mf.method == 'srf' and mf.srf_irr.meta['source'] == 'veit_l1a'
    Eobs = twin.observe(scene['lam'], scene['Ed'], gE, mf.srf_irr)
    truth = twin.observe(scene['lam'], scene['Ed'], gL, mf.srf_rad)
    ed = mf.meas_function(gL, gE, Eobs)
    lin = _processor_linear(gL, gE, Eobs)
    err_srf = np.std(ed / truth - 1)
    err_lin = np.std(lin / truth - 1)
    assert err_srf < 1e-6 < 0.01 * err_lin

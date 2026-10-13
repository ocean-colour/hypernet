"""``InterpolationWavSRF``: a drop-in for ``hypernets_processor``'s
``InterpolationWavLinear`` that puts Ed on the radiance grid with the measured
E and L SRFs (the ``srf`` method of ``docs/wiggles_planning.md`` §3).

It is written against ``hypernets_processor`` commit ``9a12819a`` (copies in
``context/hypernets_processor/``): same base class (``punpy.MeasurementFunction``),
same ``meas_function(rad_wavs, irr_wavs, irr)`` and argument names, same
``setup(szas, network, irr_wavs)`` call, and a ``get_name()`` for the config.
It calls :func:`hypernet.edinterp.interpolate_ed_to_l` with:

- the SRFs and the gate G1 method (``srf`` where FWHM_L - FWHM_E >= 0.15 nm,
  else ``linear``) from the SRF table (:mod:`hypernet.srftable`);
- the wavelength-scale corrections of E and L (task 4b), from the table's
  calibration-period medians unless a per-sequence fit is passed;
- the Emod cache entry for the sequence's mean SZA (:func:`hypernet.emod.emod_for`).

Pixels outside the E grid are **extrapolated linearly**, as the processor's
linear class does (``interp1d(..., fill_value="extrapolate")``), not held.

**What RBINS must wire up** (nothing runs end to end here):

1. ``interpolation/measurement_functions/interpolation_factory.py``: import
   ``InterpolationWavSRF`` and add it to ``interpolators``.
2. ``etc/processor_water_defaults.config``:
   ``measurement_function_interpolate_wav = InterpolationWavSRF``.
3. ``interpolation/interpolate.py``, ``Interpolate.interpolate_irradiance``:
   pass the instrument and date to ``setup``, e.g.::

       interpolation_function_wav.setup(
           np.nanmean(dataset_l1b_irr["solar_zenith_angle"].values),
           self.context.get_config_value("network"),
           dataset_l1b_irr["wavelength"].values,
           system_id=dataset_l1b_irr.attrs["system_id"],
           date=str(dataset_l1b_irr["acquisition_time"].values[0]))

   Without ``system_id`` the class falls back to linear interpolation (with
   a warning), so an unwired install behaves like today's processor.
4. Optional, for the new uncertainty terms (Build #2 Q1, route (a)): build
   the class with ``srf_inputs=True`` and add the SRF and wavelength-scale
   variables to the irradiance dataset before propagation
   (:func:`add_srf_variables`).  punpy then propagates them with the
   existing ``propagate_ds_specific`` call once ``mcsteps > 0`` (the water
   default is 0).
5. Ship ``hypernet`` (for ``srftable``, ``emod``, ``edinterp`` and the data
   files) or vendor those modules.
"""
import warnings

import numpy as np
import scipy.interpolate
from punpy import MeasurementFunction

from hypernet import edinterp, srftable
from hypernet import srf as srfmod
from hypernet.srfcal import HALF, LAM0

#: Names of the optional extra input variables (``srf_inputs=True``).  None of
#: them may contain an existing input's name right after the component
#: (punpy finds ``u_*_<component>_<var>`` by substring: a name like
#: ``wavelength_scale_*`` would be taken as an uncertainty of ``wavelength``).
SRF_INPUT_NAMES = ['srf_coeffs_irradiance', 'srf_coeffs_radiance',
                   'wavescale_irradiance', 'wavescale_radiance']


def _first_col(w):
    w = np.asarray(w, dtype=float)
    return w[:, 0] if w.ndim > 1 else w


def _scale_vec(scale):
    if scale is None:
        return np.zeros(2)
    return np.array([scale.get('shift', 0.0), scale.get('stretch', 0.0)], dtype=float)


class InterpolationWavSRF(MeasurementFunction):
    """Ed on the L grid with the measured SRFs (see the module docstring)."""

    def __init__(self, *args, srf_inputs=False, **kwargs):
        self.srf_inputs = bool(srf_inputs)
        self.method = 'linear'
        self.srf_irr = self.srf_rad = None
        self.scale_irr = self.scale_rad = None
        self.emod = None
        super().__init__(*args, **kwargs)

    def setup(self, szas, network, irr_wavs, system_id=None, date=None, *, method=None,
              srf_irr=None, srf_rad=None, scale_irr=None, scale_rad=None, emod=None,
              pwv_mm=15.0, table=None):
        """Called by the processor before propagation.

        ``szas``, ``network`` and ``irr_wavs`` are what the processor passes
        today.  ``system_id`` and ``date`` select the SRF-table rows; the
        keyword overrides (``method``, ``srf_irr``, ``srf_rad``,
        ``scale_irr``, ``scale_rad``, ``emod``) replace table or cache values
        (e.g. a per-sequence wavelength-scale fit, Setup Q8).
        """
        if system_id is not None and date is not None:
            self.srf_irr = srf_irr or srftable.srf_for(system_id, 'E', date, table)
            self.srf_rad = srf_rad or srftable.srf_for(system_id, 'Ld', date, table)
            self.method = method or srftable.method_for(system_id, date, table)
            self.scale_irr = scale_irr or srftable.scale_for(system_id, 'E', date, table)
            self.scale_rad = scale_rad or srftable.scale_for(system_id, 'Ld', date, table)
        else:
            self.srf_irr, self.srf_rad = srf_irr, srf_rad
            self.scale_irr, self.scale_rad = scale_irr, scale_rad
            self.method = method or ('linear' if srf_irr is None else
                                     edinterp.choose_method(srf_irr, srf_rad))
            if srf_irr is None:
                warnings.warn('InterpolationWavSRF.setup: no system_id/date; falling back to '
                              'linear interpolation')
        if self.method != 'linear' and emod is None:
            from hypernet import emod as emodmod
            emod = emodmod.emod_for(float(np.nanmean(szas)), pwv_mm)
        self.emod = emod
        if self.srf_inputs and self.srf_irr is None:
            raise ValueError('srf_inputs=True needs SRFs (system_id/date or srf_irr/srf_rad)')

    def nominal_inputs(self):
        """Values of the extra inputs (``srf_inputs=True``), in
        :data:`SRF_INPUT_NAMES` order."""
        return [np.asarray(self.srf_irr.coeffs, float), np.asarray(self.srf_rad.coeffs, float),
                _scale_vec(self.scale_irr), _scale_vec(self.scale_rad)]

    def meas_function(self, rad_wavs, irr_wavs, irr, *srf_args):
        """Ed on ``rad_wavs`` from ``irr`` on ``irr_wavs`` ((n_E,) or (n_E, n_scan))."""
        wr, wi = _first_col(rad_wavs), _first_col(irr_wavs)
        irr = np.asarray(irr, dtype=float)
        if srf_args:
            cE, cL, sE, sL = (np.asarray(a, float).ravel() for a in srf_args)
        else:
            cE = cL = None
            sE, sL = _scale_vec(self.scale_irr), _scale_vec(self.scale_rad)
        wi_t = wi + sE[0] + sE[1] * (wi - LAM0) / HALF
        wr_t = wr + sL[0] + sL[1] * (wr - LAM0) / HALF
        if self.method == 'linear':
            return scipy.interpolate.interp1d(wi_t, irr, axis=0, fill_value='extrapolate')(wr_t)
        a = self._srf(self.srf_irr, cE)
        b = self._srf(self.srf_rad, cL)
        out = edinterp.interpolate_ed_to_l(wi_t, irr, wr_t, emod=self.emod, srf_irr=a,
                                           srf_rad=b, method=self.method)
        outside = (wr_t < wi_t[0]) | (wr_t > wi_t[-1])
        if outside.any():
            ext = scipy.interpolate.interp1d(wi_t, irr, axis=0, fill_value='extrapolate')
            out[outside] = ext(wr_t[outside])
        return out

    @staticmethod
    def _srf(model, coeffs):
        """The SRFModel itself (cached convolutions) unless drawn coefficients differ."""
        if coeffs is None or np.array_equal(coeffs, model.coeffs):
            return model
        return lambda lam: srfmod.fwhm_at(lam, coeffs, model.lam_ref, model.lam_scale)

    @staticmethod
    def get_name():
        return 'InterpolationWavSRF'

    def get_argument_names(self):
        base = ['radiance_wavelength', 'wavelength', 'irradiance']
        return base + SRF_INPUT_NAMES if getattr(self, 'srf_inputs', False) else base


def add_srf_variables(ds_irr, mf):
    """The irradiance dataset with the extra inputs of a set-up
    ``InterpolationWavSRF(srf_inputs=True)``: the E and L FWHM coefficients
    and the E and L wavelength-scale corrections, each with an absolute
    ``u_rel_systematic_indep_*`` uncertainty (units nm) and its
    error-correlation matrix (obsarray ``err_corr_matrix``)."""
    ds = ds_irr.copy()
    covs = [mf.srf_irr.cov, mf.srf_rad.cov,
            edinterp._scale_cov(mf.scale_irr)[1], edinterp._scale_cov(mf.scale_rad)[1]]
    for name, val, cov in zip(SRF_INPUT_NAMES, mf.nominal_inputs(), covs):
        dim = 'srf_coeff' if name.startswith('srf') else 'scale_param'
        u, corr = edinterp._u_corr_from_cov(cov)
        cname = 'err_corr_systematic_indep_' + name
        ds[name] = ((dim,), val, {'units': 'nm'})
        ds[cname] = ((dim, dim + '_2'), corr)
        ds.unc[name]['u_rel_systematic_indep_' + name] = (
            (dim,), u, {'err_corr': [{'dim': dim, 'form': 'err_corr_matrix',
                                      'params': [cname], 'units': []}],
                        'pdf_shape': 'gaussian', 'units': 'nm'})
    return ds

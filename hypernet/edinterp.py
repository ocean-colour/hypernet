"""Resample the downwelling irradiance Ed from the E grid onto the L grid.

One function, :func:`interpolate_ed_to_l`, implements the methods of the
wiggles plan (``docs/wiggles_planning.md`` §3).  With E_i measured at the E
pixel centres lambda_E^i, an L pixel lambda_L^j between lambda_E^i and
lambda_E^{i+1}, and the linear weight w_j:

- ``linear`` (paper eqs. 11-12; the processor's ``interpolate_wav_linear``)::

      E(lambda_L) = (1 - w) E_i + w E_{i+1}

- ``ruddick2023`` (Ruddick et al. 2023, eqs. 14-15): a high-resolution model
  irradiance Emod is convolved with the **E** SRF on both sides::

      E(lambda_L) = Emod_E(lambda_L) [(1 - w) E_i / Emod_E(lambda_E^i)
                                       + w E_{i+1} / Emod_E(lambda_E^{i+1})]

- ``srf`` (this work): the numerator uses the **L** SRF, so the output has
  the resolution of the channel it will be divided into::

      E(lambda_L) = Emod_L(lambda_L) [(1 - w) E_i / Emod_E(lambda_E^i)
                                       + w E_{i+1} / Emod_E(lambda_E^{i+1})]

  With srf_rad = srf_irr this *is* ``ruddick2023``; with a constant Emod it
  is ``linear``.
- ``cubic`` (scipy CubicSpline) and ``sinc`` (Lanczos-windowed sinc in
  fractional-pixel space, a = 8) are nulls: different interpolators that use
  no model.

The bracket is the measured-over-model ratio, which is smooth if the model is
good, so it is interpolated linearly.

This is the minimal tested version of Phase 1 task 8; Phase 2 adds
validation, caching and uncertainty.
"""

import numpy as np
from scipy.interpolate import CubicSpline

METHODS = ('linear', 'ruddick2023', 'srf', 'cubic', 'sinc')


def _emod_arrays(emod):
    """(lam, flux) from a tuple/list or a dict with 'lam' and 'E'/'Ed'/'Emod'."""
    if isinstance(emod, dict):
        for k in ('E', 'Ed', 'Emod'):
            if k in emod:
                return np.asarray(emod['lam'], float), np.asarray(emod[k], float)
        raise KeyError("emod dict needs 'lam' and one of 'E', 'Ed', 'Emod'")
    lam, flux = emod
    return np.asarray(lam, float), np.asarray(flux, float)


def model_on_grid(emod, grid, srf):
    """Emod convolved with the Gaussian ``srf`` (anything
    :func:`hypernet.twin.convolve_to_grid` accepts, or an SRFModel) at
    ``grid``."""
    from hypernet.twin import convolve_to_grid
    lam, flux = _emod_arrays(emod)
    f = srf.fwhm if hasattr(srf, 'fwhm') else srf
    return convolve_to_grid(lam, flux, grid, f)


def _srf_on(srf, own_grid, target_grid):
    """An SRF spec usable at ``target_grid``: a per-pixel FWHM array given on
    ``own_grid`` is interpolated; constants, callables and SRFModels pass."""
    if hasattr(srf, 'fwhm') or callable(srf) or np.ndim(srf) == 0:
        return srf
    a = np.asarray(srf, dtype=float)
    if a.size == np.size(target_grid):
        return a
    return np.interp(target_grid, own_grid, a)


def _lanczos(wav_irr, irradiance, wav_rad, a=8):
    """Windowed-sinc interpolation in fractional-pixel space of the E grid."""
    x = np.interp(wav_rad, wav_irr, np.arange(wav_irr.size))      # fractional index
    out = np.empty(x.size)
    n = wav_irr.size
    for j, xj in enumerate(x):
        k0 = int(np.floor(xj))
        k = np.arange(k0 - a + 1, k0 + a + 1)
        k = k[(k >= 0) & (k < n)]
        d = xj - k
        w = np.sinc(d) * np.sinc(d / a)
        out[j] = np.dot(w, irradiance[k]) / w.sum()
    return out


def interpolate_ed_to_l(wav_irr, irradiance, wav_rad, *, emod=None, srf_irr=None,
                        srf_rad=None, method='srf'):
    """Ed on the radiance (L) grid.

    Parameters
    ----------
    wav_irr, irradiance : array
        E-grid wavelengths (nm, ascending) and the measured irradiance.
    wav_rad : array
        L-grid wavelengths (nm).  Outside ``wav_irr`` the linear-weight
        methods hold the end values (as ``np.interp``).
    emod : tuple or dict, optional
        High-resolution model irradiance: ``(lam, flux)`` or a dict with
        ``'lam'`` and ``'E'``/``'Ed'``/``'Emod'`` (e.g. a ``hypernet.twin``
        scene).  Required for ``ruddick2023`` and ``srf``.
    srf_irr, srf_rad : float, array, callable or SRFModel
        Gaussian SRF FWHM (nm) of the E and L channels (an array must match
        its own grid).  ``ruddick2023`` uses ``srf_irr`` only.
    method : str
        One of :data:`METHODS`.

    Returns
    -------
    numpy.ndarray
        Ed on ``wav_rad``.
    """
    wav_irr = np.asarray(wav_irr, dtype=float)
    E = np.asarray(irradiance, dtype=float)
    wav_rad = np.asarray(wav_rad, dtype=float)
    if method == 'linear':
        return np.interp(wav_rad, wav_irr, E)
    if method == 'cubic':
        return CubicSpline(wav_irr, E)(wav_rad)
    if method == 'sinc':
        return _lanczos(wav_irr, E, wav_rad)
    if method not in ('ruddick2023', 'srf'):
        raise ValueError('method must be one of %s' % (METHODS,))
    if emod is None or srf_irr is None:
        raise ValueError('%s needs emod and srf_irr' % method)
    emod_E_at_E = model_on_grid(emod, wav_irr, srf_irr)
    ratio = np.interp(wav_rad, wav_irr, E / emod_E_at_E)
    if method == 'ruddick2023':
        return model_on_grid(emod, wav_rad, _srf_on(srf_irr, wav_irr, wav_rad)) * ratio
    if srf_rad is None:
        raise ValueError('srf needs srf_rad')
    return model_on_grid(emod, wav_rad, srf_rad) * ratio

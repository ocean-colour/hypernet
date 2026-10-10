"""Resample the downwelling irradiance Ed from the E grid onto the L grid.

One function, :func:`interpolate_ed_to_l`, implements the methods of the
wiggles plan (``docs/wiggles_planning.md`` §3).  Notation: E_i is the measured
irradiance at the E pixel centre lambda_E^i; an L pixel lambda_L^j lies
between lambda_E^i and lambda_E^{i+1}; w_j = (lambda_L^j - lambda_E^i) /
(lambda_E^{i+1} - lambda_E^i) is its linear weight; Emod is a high-resolution
model irradiance, and Emod_X(lambda_c) is Emod convolved with the X-channel
SRF centred on lambda_c:

    Emod_X(lambda_c) = int SRF_X(lambda; lambda_c) Emod(lambda) dlambda
                       / int SRF_X(lambda; lambda_c) dlambda

Methods (:data:`METHODS`):

- ``linear`` (paper eqs. 11-12; the processor's ``InterpolationWavLinear``)::

      E(lambda_L^j) = (1 - w_j) E_i + w_j E_{i+1}

- ``ruddick2023`` (Ruddick et al. 2023, eqs. 14-15), the model seen through
  the **E** SRF on both sides::

      E(lambda_L^j) = Emod_E(lambda_L^j) [(1 - w_j) E_i / Emod_E(lambda_E^i)
                                          + w_j E_{i+1} / Emod_E(lambda_E^{i+1})]

- ``srf`` (this work), the **L** SRF in front, so the output has the
  resolution of the channel it will be divided into::

      E(lambda_L^j) = Emod_L(lambda_L^j) [(1 - w_j) E_i / Emod_E(lambda_E^i)
                                          + w_j E_{i+1} / Emod_E(lambda_E^{i+1})]

  With SRF_L = SRF_E this *is* ``ruddick2023``; with a constant Emod it is
  ``linear``.  The bracket is the measured-over-model ratio, smooth if the
  model's line shapes are right, so it is interpolated linearly.

``cubic`` and ``sinc`` (:data:`NULL_METHODS`) are model-free nulls kept for
the Phase 1 twin experiment; they are not deliverables.

Gate G1 (``claude_prompts/wiggles/gate_G1.md``) chose ``srf`` where the
measured FWHM_L - FWHM_E is >= 0.15 nm and ``linear`` otherwise:
:func:`choose_method`.

Uncertainties follow the L1C convention: three relative components in %,
``u_rel_random``, ``u_rel_systematic_indep`` and
``u_rel_systematic_corr_rad_irr``, on (wavelength[, scan]).  This module
propagates the *input* components analytically (first order); Phase 2 task 7
adds the SRF, Emod and wavelength-scale terms with punpy.
"""

import numpy as np
from scipy.interpolate import CubicSpline

from hypernet import srf as srfmod

METHODS = ('linear', 'ruddick2023', 'srf')
#: Model-free nulls for the twin experiment; accepted, not deliverables.
NULL_METHODS = ('cubic', 'sinc')
#: The L1C relative-uncertainty components (``u_rel_<name>_irradiance``).
UNC_COMPONENTS = ('random', 'systematic_indep', 'systematic_corr_rad_irr')
#: Gate G1: ``srf`` where the median FWHM_L - FWHM_E over DFWHM_BAND >= this (nm).
SRF_THRESHOLD = 0.15
#: Band (nm) where both template SRF models are constrained (Setup Q7).
DFWHM_BAND = (400.0, 680.0)
#: An L pixel falls back to linear if the SRF (+- this many sigma) of it or of
#: its two E neighbours reaches beyond the model's wavelength range.
NSIG_EDGE = 3.0


# --- the G1 method rule ------------------------------------------------------------

def delta_fwhm(srf_irr, srf_rad, band=DFWHM_BAND, step=1.0):
    """Median of FWHM_L(lambda) - FWHM_E(lambda) over ``band`` (nm).

    ``srf_irr`` and ``srf_rad`` are SRFModels, callables or constants
    (per-pixel arrays have no wavelengths attached; pass a callable).
    """
    for s in (srf_irr, srf_rad):
        if np.ndim(s) > 0 and not hasattr(s, 'fwhm') and not callable(s):
            raise ValueError('delta_fwhm needs SRFModels, callables or constants')
    lam = np.arange(band[0], band[1] + 0.5 * step, step)
    return float(np.median(srfmod.fwhm_values(srf_rad, lam) - srfmod.fwhm_values(srf_irr, lam)))


def choose_method(srf_irr, srf_rad, threshold=SRF_THRESHOLD, band=DFWHM_BAND):
    """Gate G1's rule: ``'srf'`` if :func:`delta_fwhm` >= ``threshold`` (nm),
    else ``'linear'``.

    Below ~0.1 nm the wiggle excess is smaller than the true line structure
    and ``srf`` is fragile to SRF errors (Phase 1 tasks 10-11), so the
    processor's linear interpolation is kept there.
    """
    # 1e-9 nm tolerance, so that e.g. 2.8 - 2.65 counts as 0.15
    return 'srf' if delta_fwhm(srf_irr, srf_rad, band) >= threshold - 1e-9 else 'linear'


# --- the model on a grid -----------------------------------------------------------

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
    """Emod convolved with the Gaussian ``srf`` at ``grid``: Emod_X(grid).

    ``srf`` is anything :func:`hypernet.srf.convolve_to_grid` accepts.  The
    convolution matrix is cached per (Emod grid, ``grid``, SRF), so repeated
    calls with the same SRFModel cost one sparse product.
    """
    if isinstance(emod, dict) and 'key' in emod:      # an emod.emod_for entry
        from hypernet import emod as emodmod
        return emodmod.emod_on_grid(emod, grid, srf)
    lam, flux = _emod_arrays(emod)
    return srfmod.convolve_to_grid(lam, flux, grid, srf)


def _srf_on(srf, own_grid, target_grid):
    """An SRF spec usable at ``target_grid``: a per-pixel FWHM array given on
    ``own_grid`` is interpolated; constants, callables and SRFModels pass."""
    if hasattr(srf, 'fwhm') or callable(srf) or np.ndim(srf) == 0:
        return srf
    a = np.asarray(srf, dtype=float)
    if a.size == np.size(target_grid):
        return a
    return np.interp(target_grid, own_grid, a)


# --- interpolation kernels ----------------------------------------------------------

def linear_weights(wav_irr, wav_rad):
    """Left neighbour index and weight of each L pixel on the E grid.

    Returns ``(i, w)`` with E(lambda_L^j) = (1 - w_j) E_i + w_j E_{i+1}.
    Outside the E grid the end value is held (``np.interp`` semantics):
    below, i = 0 and w = 0; above, i = n - 2 and w = 1.
    """
    n = wav_irr.size
    i = np.clip(np.searchsorted(wav_irr, wav_rad, side='right') - 1, 0, n - 2)
    w = (wav_rad - wav_irr[i]) / (wav_irr[i + 1] - wav_irr[i])
    return i, np.clip(w, 0.0, 1.0)


def _combine(x, i, w):
    """(1 - w) x[i] + w x[i+1] along axis 0; a NaN reaches only the L pixels
    that use it with non-zero weight."""
    ww = w.reshape((-1,) + (1,) * (x.ndim - 1))
    a, b = x[i], x[i + 1]
    out = (1 - ww) * a + ww * b
    out = np.where(ww == 0, a, out)
    return np.where(ww == 1, b, out)


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
        wk = np.sinc(d) * np.sinc(d / a)
        out[j] = np.dot(wk, irradiance[k]) / wk.sum()
    return out


def _null(method, wav_irr, E, wav_rad):
    cols = E.reshape(E.shape[0], -1).T
    if method == 'cubic':
        out = [CubicSpline(wav_irr, c)(wav_rad) for c in cols]
    else:
        out = [_lanczos(wav_irr, c, wav_rad) for c in cols]
    return np.stack(out, axis=-1).reshape((wav_rad.size,) + E.shape[1:])


def _edge_mask(emod, wav_irr, wav_rad, i, srf_irr, srf_num):
    """L pixels whose model terms are unreliable: outside the E grid, or an
    SRF (+- NSIG_EDGE sigma) of the pixel or its E neighbours beyond Emod."""
    lam, _ = _emod_arrays(emod)
    lo, hi = lam[0], lam[-1]
    sL = NSIG_EDGE * srfmod.fwhm_values(srf_num, wav_rad) / srfmod.FWHM_PER_SIGMA
    sE = NSIG_EDGE * srfmod.fwhm_values(srf_irr, wav_irr) / srfmod.FWHM_PER_SIGMA
    bad = (wav_rad < wav_irr[0]) | (wav_rad > wav_irr[-1])
    bad |= (wav_rad - sL < lo) | (wav_rad + sL > hi)
    for k in (i, i + 1):
        bad |= (wav_irr[k] - sE[k] < lo) | (wav_irr[k] + sE[k] > hi)
    return bad


def _propagate_inputs(r, i, w, u, systematic):
    """First-order relative uncertainty (%) of (1-w) r_i + w r_{i+1}.

    r is E (linear) or E / Emod_E (model methods): the deterministic model
    factors cancel in relative terms.  Random components add in quadrature;
    systematic ones are taken as fully correlated between the two
    neighbouring E pixels (adjacent pixels; task 7 uses the err_corr
    matrices).
    """
    ww = w.reshape((-1,) + (1,) * (r.ndim - 1))
    t0, t1 = (1 - ww) * r[i], ww * r[i + 1]
    tot = t0 + t1
    a, b = t0 / tot, t1 / tot
    u = np.broadcast_to(u, r.shape)
    ui, uj = u[i], u[i + 1]
    if systematic:
        return np.abs(a * ui + b * uj)
    return np.sqrt((a * ui) ** 2 + (b * uj) ** 2)


# --- the function -------------------------------------------------------------------

def interpolate_ed_to_l(wav_irr, irradiance, wav_rad, *, emod=None, srf_irr=None,
                        srf_rad=None, method='srf', u_rel_random=None,
                        u_rel_systematic_indep=None, u_rel_systematic_corr_rad_irr=None,
                        full_output=None):
    """Ed on the radiance (L) grid.

    Parameters
    ----------
    wav_irr : array (n_E,)
        E-grid wavelengths, nm: finite and strictly increasing.
    irradiance : array (n_E,) or (n_E, n_scan)
        The measured irradiance.  NaNs reach only the L pixels that use them.
    wav_rad : array (n_L,)
        L-grid wavelengths, nm (finite).  Outside ``wav_irr`` the end value is
        held (``np.interp`` semantics; the processor extrapolates instead, see
        the task 8 snippet).
    emod : tuple or dict, optional
        High-resolution model irradiance: ``(lam, flux)`` or a dict with
        ``'lam'`` and ``'E'``/``'Ed'``/``'Emod'``.  Required for
        ``ruddick2023`` and ``srf``.
    srf_irr, srf_rad : float, array, callable or SRFModel
        Gaussian SRF FWHM (nm) of the E and L channels (an array is per pixel
        of its own grid).  ``ruddick2023`` uses ``srf_irr`` only.
    method : str
        One of :data:`METHODS` (or :data:`NULL_METHODS` for the twin).  See
        :func:`choose_method` for gate G1's rule.
    u_rel_random, u_rel_systematic_indep, u_rel_systematic_corr_rad_irr :
            array, optional
        Relative uncertainties of ``irradiance`` in % (L1C convention),
        broadcastable to its shape.  Propagated to first order (random in
        quadrature, systematic fully correlated between the two neighbours).
    full_output : bool, optional
        Return ``(Ed, info)``.  Defaults to True when any ``u_rel_*`` is given.

    Returns
    -------
    numpy.ndarray
        Ed on ``wav_rad``, shape (n_L,) or (n_L, n_scan).
    dict
        Only with ``full_output``: ``method``; ``fallback`` (bool (n_L,), L
        pixels set by ``linear`` because the model terms reach beyond Emod or
        the E grid); and ``u_rel_<component>`` (%) for each given input
        component.
    """
    wav_irr = np.asarray(wav_irr, dtype=float)
    E = np.asarray(irradiance, dtype=float)
    wav_rad = np.asarray(wav_rad, dtype=float)
    # validation
    if wav_irr.ndim != 1 or wav_irr.size < 2:
        raise ValueError('wav_irr must be 1-D with at least 2 pixels')
    if not np.all(np.isfinite(wav_irr)) or np.any(np.diff(wav_irr) <= 0):
        raise ValueError('wav_irr must be finite and strictly increasing')
    if wav_rad.ndim != 1 or not np.all(np.isfinite(wav_rad)):
        raise ValueError('wav_rad must be 1-D and finite')
    if E.ndim not in (1, 2) or E.shape[0] != wav_irr.size:
        raise ValueError('irradiance must be (n_E,) or (n_E, n_scan) with n_E = '
                         'wav_irr.size (got %s for %d)' % (E.shape, wav_irr.size))
    if method not in METHODS + NULL_METHODS:
        raise ValueError('method must be one of %s' % (METHODS,))
    if method in ('ruddick2023', 'srf') and (emod is None or srf_irr is None):
        raise ValueError('%s needs emod and srf_irr' % method)
    if method == 'srf' and srf_rad is None:
        raise ValueError('srf needs srf_rad')
    unc = {c: u for c, u in zip(UNC_COMPONENTS, (u_rel_random, u_rel_systematic_indep,
                                                 u_rel_systematic_corr_rad_irr))
           if u is not None}
    if full_output is None:
        full_output = bool(unc)
    if unc and method in NULL_METHODS:
        raise ValueError('uncertainties are not propagated for the null methods')

    i, w = linear_weights(wav_irr, wav_rad)
    fallback = np.zeros(wav_rad.size, dtype=bool)
    if method in NULL_METHODS:
        out = _null(method, wav_irr, E, wav_rad)
        r = None
    elif method == 'linear':
        r = E
        out = _combine(E, i, w)
    else:
        mE = model_on_grid(emod, wav_irr, srf_irr)
        shape = (-1,) + (1,) * (E.ndim - 1)
        r = E / mE.reshape(shape)
        num = (_srf_on(srf_irr, wav_irr, wav_rad) if method == 'ruddick2023' else srf_rad)
        mL = model_on_grid(emod, wav_rad, num)
        out = mL.reshape(shape) * _combine(r, i, w)
        fallback = _edge_mask(emod, wav_irr, wav_rad, i, srf_irr, num) | ~np.isfinite(mL)
        if fallback.any():
            out[fallback] = _combine(E, i, w)[fallback]
    if not full_output:
        return out
    info = dict(method=method, fallback=fallback)
    for c, u in unc.items():
        u = np.asarray(u, dtype=float)
        sys_ = c != 'random'
        v = _propagate_inputs(r, i, w, u, sys_)
        if fallback.any():
            v[fallback] = _propagate_inputs(E, i, w, u, sys_)[fallback]
        info['u_rel_' + c] = v
    return out, info

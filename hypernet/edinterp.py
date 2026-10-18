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


# --- uncertainty propagation with punpy (Phase 2 task 7) ------------------------------

def _u_corr_from_cov(cov):
    """1-sigma vector and correlation matrix of a covariance (zero-variance
    entries get unit diagonal)."""
    cov = np.atleast_2d(np.asarray(cov, dtype=float))
    u = np.sqrt(np.clip(np.diag(cov), 0, None))
    with np.errstate(invalid='ignore', divide='ignore'):
        corr = cov / np.outer(u, u)
    corr[~np.isfinite(corr)] = 0.0
    np.fill_diagonal(corr, 1.0)
    return u, corr


def _scale_cov(scale):
    """(values, cov) of a wavelength-scale correction dict (shift, stretch)."""
    if scale is None:
        return np.zeros(2), np.zeros((2, 2))
    v = np.array([scale.get('shift', 0.0), scale.get('stretch', 0.0)], dtype=float)
    if 'cov' in scale and np.all(np.isfinite(scale['cov'])):
        return v, np.asarray(scale['cov'], dtype=float)
    e = np.array([scale.get('shift_err', 0.0), scale.get('stretch_err', 0.0)], dtype=float)
    return v, np.diag(np.nan_to_num(e) ** 2)


def _corrected(wav, s):
    from hypernet.srfcal import HALF, LAM0
    return wav + s[0] + s[1] * (wav - LAM0) / HALF


def _fwhm_fn(model, coeffs):
    """A plain callable FWHM(lambda) with drawn coefficients (not cached by
    hypernet.srf, so MC draws do not fill the matrix cache), carrying the
    model's SRF shape."""
    shape, par = srfmod.srf_shape(model)
    return srfmod.ShapedFWHM(
        lambda lam: srfmod.fwhm_at(lam, coeffs, model.lam_ref, model.lam_scale), shape, par)


def propagate_ed_to_l(wav_irr, irradiance, wav_rad, *, emod=None, srf_irr=None, srf_rad=None,
                      method='srf', u_rel_random=None, u_rel_systematic_indep=None,
                      u_rel_systematic_corr_rad_irr=None, err_corr_systematic_indep=None,
                      err_corr_systematic_corr_rad_irr=None, scale_irr=None, scale_rad=None,
                      emod_perturbations=None, mcsteps=100, seed=None, return_corr=False):
    """Ed on the L grid with its three relative uncertainty components (%),
    propagated with ``punpy.MCPropagation``.

    One MC run per L1C component, each run perturbing:

    - ``random``: E with ``u_rel_random`` (uncorrelated in wavelength);
    - ``systematic_corr_rad_irr``: E with that component and its wavelength
      error-correlation matrix (default fully correlated);
    - ``systematic_indep``: E with that component and matrix, **and** the new
      terms, common to every scan of an instrument: the E and L FWHM(lambda)
      coefficients (``SRFModel.cov``) and the E and L wavelength-scale
      corrections (shift, stretch with covariance; task 4b), drawn jointly.

    The Emod inputs (SZA, water vapour, ozone) are added to
    ``systematic_indep`` in quadrature, linearised by finite differences
    (``emod_perturbations``, e.g. :func:`hypernet.emod.emod_perturbations`):
    rebuilding Emod per draw would take seconds each.  A pure F0 scale cancels
    in the model methods and is not a term.

    The nominal Ed is :func:`interpolate_ed_to_l` on the scale-corrected grids
    (``scale_irr`` / ``scale_rad`` applied to ``wav_irr`` / ``wav_rad`` for the
    model and the weights; the output is reported on ``wav_rad``).

    Parameters
    ----------
    wav_irr, irradiance, wav_rad, emod, method :
        As :func:`interpolate_ed_to_l`; ``irradiance`` 1-D or (n_E, n_scan)
        (2-D runs each scan separately).
    srf_irr, srf_rad : SRFModel
        Needed for the model methods; their ``cov`` gives the FWHM terms.
    u_rel_random, u_rel_systematic_indep, u_rel_systematic_corr_rad_irr :
            array or float, optional
        Input relative uncertainties of ``irradiance`` (%); missing = 0.
    err_corr_systematic_indep, err_corr_systematic_corr_rad_irr :
            array (n_E, n_E), optional
        Wavelength error-correlation of those components (default ``'syst'``).
    scale_irr, scale_rad : dict, optional
        Wavelength-scale corrections, ``shift`` and ``stretch`` (nm) with
        ``cov`` or ``shift_err`` / ``stretch_err`` (e.g.
        :func:`hypernet.srfcal.fit_wavelength_scale`,
        :func:`hypernet.srftable.scale_for`).
    emod_perturbations : list of (name, alt_emod, factor), optional
        Linearised Emod terms.
    mcsteps : int
        MC draws per run.
    seed : int, optional
        Seeds numpy's global RNG (which punpy / comet_maths use) for the call;
        the previous RNG state is restored afterwards.
    return_corr : bool
        Also return each component's wavelength error-correlation matrix
        on the L grid (1-D input only).

    Returns
    -------
    numpy.ndarray
        Ed on ``wav_rad``.
    dict
        ``u_rel_random``, ``u_rel_systematic_indep``,
        ``u_rel_systematic_corr_rad_irr`` (%), the breakdown ``u_rel_emod``
        (%) and ``u_rel_emod_<name>``, ``fallback``, ``method``, ``mcsteps``,
        and ``err_corr_<component>`` with ``return_corr``.
    """
    import punpy

    wav_irr = np.asarray(wav_irr, dtype=float)
    wav_rad = np.asarray(wav_rad, dtype=float)
    E = np.asarray(irradiance, dtype=float)
    if method not in METHODS:
        raise ValueError('propagate_ed_to_l: method must be one of %s' % (METHODS,))
    model = method != 'linear'
    if model and not (hasattr(srf_irr, 'cov') and (method == 'ruddick2023' or
                                                    hasattr(srf_rad, 'cov'))):
        raise ValueError('%s propagation needs SRFModel srf_irr/srf_rad' % method)
    if E.ndim == 2:
        cols = [propagate_ed_to_l(
            wav_irr, E[:, k], wav_rad, emod=emod, srf_irr=srf_irr, srf_rad=srf_rad,
            method=method, u_rel_random=_col(u_rel_random, k, E.shape),
            u_rel_systematic_indep=_col(u_rel_systematic_indep, k, E.shape),
            u_rel_systematic_corr_rad_irr=_col(u_rel_systematic_corr_rad_irr, k, E.shape),
            err_corr_systematic_indep=err_corr_systematic_indep,
            err_corr_systematic_corr_rad_irr=err_corr_systematic_corr_rad_irr,
            scale_irr=scale_irr, scale_rad=scale_rad, emod_perturbations=emod_perturbations,
            mcsteps=mcsteps, seed=None if seed is None else seed + k) for k in range(E.shape[1])]
        ed = np.stack([c[0] for c in cols], axis=1)
        info = {k: np.stack([c[1][k] for c in cols], axis=1) for k in cols[0][1]
                if k.startswith('u_rel')}
        info.update(method=method, mcsteps=mcsteps,
                    fallback=np.any([c[1]['fallback'] for c in cols], axis=0))
        return ed, info

    sv_irr, cov_irr = _scale_cov(scale_irr)
    sv_rad, cov_rad = _scale_cov(scale_rad)
    cE = np.asarray(srf_irr.coeffs, float) if model else np.zeros(3)
    cL = np.asarray(srf_rad.coeffs, float) if method == 'srf' else np.zeros(3)
    covE = np.asarray(srf_irr.cov, float) if model else np.zeros((3, 3))
    covL = np.asarray(srf_rad.cov, float) if method == 'srf' else np.zeros((3, 3))

    def nominal_srfs():
        a = srf_irr if model else None
        b = srf_rad if method == 'srf' else None
        return a, b

    wi0, wr0 = _corrected(wav_irr, sv_irr), _corrected(wav_rad, sv_rad)
    a0, b0 = nominal_srfs()
    ed, base_info = interpolate_ed_to_l(wi0, E, wr0, emod=emod, srf_irr=a0, srf_rad=b0,
                                        method=method, full_output=True)

    def f_e(e):                                     # E only: grids and SRFs nominal
        return interpolate_ed_to_l(wi0, e, wr0, emod=emod, srf_irr=a0, srf_rad=b0,
                                   method=method)

    def f_all(e, ce, cl, si, sr):                   # E, SRF coefficients, scales
        sa = _fwhm_fn(srf_irr, ce) if model else None
        sb = _fwhm_fn(srf_rad, cl) if method == 'srf' else None
        return interpolate_ed_to_l(_corrected(wav_irr, si), e, _corrected(wav_rad, sr),
                                   emod=emod, srf_irr=sa, srf_rad=sb, method=method)

    def absu(u):
        return np.zeros_like(E) if u is None else np.broadcast_to(
            np.asarray(u, float), E.shape) / 100.0 * np.abs(E)

    prop = punpy.MCPropagation(int(mcsteps))
    state = np.random.get_state()
    out = dict(method=method, mcsteps=int(mcsteps), fallback=base_info['fallback'])
    try:
        if seed is not None:
            np.random.seed(seed)
        runs = (('random', u_rel_random, 'rand'),
                ('systematic_corr_rad_irr', u_rel_systematic_corr_rad_irr,
                 'syst' if err_corr_systematic_corr_rad_irr is None
                 else np.asarray(err_corr_systematic_corr_rad_irr, float)))
        for name, u, corr in runs:
            uE = absu(u)
            if not np.any(uE > 0):
                out['u_rel_' + name] = np.zeros(wav_rad.size)
                if return_corr:
                    out['err_corr_' + name] = np.eye(wav_rad.size)
                continue
            res = prop.propagate_standard(f_e, [E], [uE], [corr], return_corr=return_corr)
            uy, cy = res if return_corr else (res, None)
            out['u_rel_' + name] = 100.0 * np.asarray(uy) / np.abs(ed)
            if return_corr:
                out['err_corr_' + name] = cy
        # systematic_indep: E's component plus SRF and wavelength-scale terms
        uE = absu(u_rel_systematic_indep)
        corrE = ('syst' if err_corr_systematic_indep is None
                 else np.asarray(err_corr_systematic_indep, float))
        x, ux, cx = [E], [uE], [corrE]
        for v, cov in ((cE, covE), (cL, covL), (sv_irr, cov_irr), (sv_rad, cov_rad)):
            u_, c_ = _u_corr_from_cov(cov)
            x.append(v)
            ux.append(u_)
            cx.append(c_)
        if any(np.any(np.asarray(u_) > 0) for u_ in ux):
            res = prop.propagate_standard(f_all, x, ux, cx, return_corr=return_corr)
            uy, cy = res if return_corr else (res, None)
            u_mc = 100.0 * np.asarray(uy) / np.abs(ed)
        else:
            u_mc, cy = np.zeros(wav_rad.size), np.eye(wav_rad.size)
    finally:
        np.random.set_state(state)
    # linearised Emod terms
    u_emod2 = np.zeros(wav_rad.size)
    for name, alt, factor in (emod_perturbations or []):
        if not model:
            out['u_rel_emod_' + name] = np.zeros(wav_rad.size)
            continue
        ed_alt = interpolate_ed_to_l(wi0, E, wr0, emod=alt, srf_irr=a0, srf_rad=b0,
                                     method=method)
        d = 100.0 * factor * (ed_alt - ed) / np.abs(ed)
        out['u_rel_emod_' + name] = np.abs(d)
        u_emod2 += d ** 2
    out['u_rel_emod'] = np.sqrt(u_emod2)
    out['u_rel_systematic_indep'] = np.sqrt(u_mc ** 2 + u_emod2)
    if return_corr:
        out['err_corr_systematic_indep'] = cy
    return ed, out


def _col(u, k, shape):
    if u is None:
        return None
    a = np.broadcast_to(np.asarray(u, float), shape)
    return a[:, k]

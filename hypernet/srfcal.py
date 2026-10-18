"""SRF self-calibration from a measured spectrum (Phase 2 task 4).

:func:`calibrate_srf` turns one spectrum (e.g. the scan-mean E of a sequence)
into an :class:`hypernet.srf.SRFModel`: Gaussian FWHM(lambda) as a quadratic
with covariance, and a mean centroid offset.  It uses the Phase 0 standard,
the TSIS-1 HSRS template fit (:func:`hypernet.srf.fit_template_windows`, 10 nm
windows avoiding telluric bands, then :meth:`SRFModel.from_lines`).  Gaussian
line fits (:func:`hypernet.srf.fit_lines`) measure SRF x intrinsic line
width and remain a diagnostic.

**Quality flags** (``model.meta['flags']``; ``model.meta['ok']`` is True only
with none):

- ``low_blue_signal``: the low-sun risk of plan §9.  Raised if the blue
  (``BLUE`` nm) has fewer than ``MIN_BLUE_WINDOWS`` good windows, or the
  median per-pixel S/N there is below ``MIN_BLUE_SNR`` (when ``err`` is
  given), or the model's FWHM error at ``LAM_BLUE_CHECK`` exceeds
  ``MAX_FWHM_ERR``;
- ``few_windows``: fewer than ``MIN_WINDOWS`` good windows overall;
- ``fit_failed``: the FWHM model could not be fitted (NaN coefficients).

``MAX_FWHM_ERR`` = 0.05 nm sits inside gate G1's FWHM_E tolerance
(≲ 0.07 nm for 80 % of the srf gain).

:func:`fit_wavelength_scale` (Phase 2 task 4b) measures a spectrum's
wavelength-scale error as a shift and a linear stretch, the binding
requirement of gate G1 (E to ≤ 0.05 nm relative to the model).
:func:`correct_wavelengths` applies it; :func:`scale_or_fallback` implements
Setup Q8's per-sequence-with-fallback rule.
"""

import numpy as np

from scipy.optimize import curve_fit

from hypernet import srf as srfmod

#: Blue band (nm) whose signal the low-sun flag checks.
BLUE = (390.0, 450.0)
MIN_BLUE_WINDOWS = 3
MIN_BLUE_SNR = 100.0
LAM_BLUE_CHECK = 410.0
MAX_FWHM_ERR = 0.05
MIN_WINDOWS = 8

_REF = {}


def default_reference(frame='air'):
    """TSIS-1 HSRS over the Phase 0 range (``hypernet.whn_srf.load_reference``),
    loaded once per frame."""
    if frame not in _REF:
        from hypernet import whn_srf
        _REF[frame] = whn_srf.load_reference(frame=frame)
    return _REF[frame]


def calibrate_srf(wav, spectrum, channel='E', err=None, ref=None, windows=None, veil=False,
                  via_grid=None, shape='gauss', scale_cov=True, instrument='',
                  return_fits=False, meta=None):
    """Fit an SRF model to one spectrum.

    Parameters
    ----------
    wav, spectrum : array
        Pixel wavelengths (nm, air scale for the default reference) and the
        spectrum (any units).
    channel : str
        ``'E'``, ``'Ld'`` or ``'Lu'``; stored in the model.
    err : array, optional
        1-sigma per-pixel error of ``spectrum`` (e.g. the flattened scan
        scatter / sqrt(n), :func:`hypernet.srf.scan_errors`).  Enables the S/N
        check and absolute weights.
    ref : (ref_wave, ref_flux), optional
        High-resolution reference; default TSIS-1 HSRS in air
        (:func:`default_reference`).
    windows : list of (lo, hi), optional
        Default :func:`hypernet.srf.template_windows` (390-880 nm, 10 nm,
        telluric bands excluded).
    veil : bool
        Fit an additive veil per window.  Default False, as in the shipped
        VEIT models (the free veil was consistent with 0 there).
    via_grid, shape : passed to :func:`hypernet.srf.fit_template_windows`.
    scale_cov : bool
        Inflate the FWHM covariance by chi2_nu when the scatter exceeds the
        window errors (:func:`hypernet.srf.fit_fwhm_model`).
    instrument : str
        Stored in the model.
    return_fits : bool
        Also return the per-window table.
    meta : dict, optional
        Extra provenance merged into ``model.meta``.

    Returns
    -------
    SRFModel
        With ``meta``: ``method='template'``, ``veil``, ``n_windows``,
        ``n_blue``, ``blue_snr``, ``fwhm_err_blue``, ``flags`` (list) and
        ``ok``.
    pandas.DataFrame
        Only with ``return_fits``: the window fits.
    """
    wav = np.asarray(wav, dtype=float)
    spectrum = np.asarray(spectrum, dtype=float)
    if wav.ndim != 1 or spectrum.shape != wav.shape:
        raise ValueError('wav and spectrum must be 1-D and the same length')
    if ref is None:
        ref = default_reference('air')
    fits = srfmod.fit_template_windows(wav, spectrum, ref[0], ref[1], windows=windows,
                                       err=err, veil=veil, via_grid=via_grid, shape=shape)
    model = srfmod.SRFModel.from_lines(channel, fits, scale_cov=scale_cov,
                                       instrument=instrument)
    good = fits['ok'].astype(bool)
    blue = good & (fits['lam_center'] >= BLUE[0]) & (fits['lam_center'] <= BLUE[1])
    blue_snr = np.nan
    if err is not None:
        e = np.asarray(err, dtype=float)
        mb = (wav >= BLUE[0]) & (wav <= BLUE[1]) & np.isfinite(e) & (e > 0)
        if mb.any():
            blue_snr = float(np.nanmedian(spectrum[mb] / e[mb]))
    fit_ok = np.all(np.isfinite(model.coeffs))
    fwhm_err_blue = float(model.fwhm_err(LAM_BLUE_CHECK)) if fit_ok else np.nan
    flags = []
    if not fit_ok:
        flags.append('fit_failed')
    if (blue.sum() < MIN_BLUE_WINDOWS or (np.isfinite(blue_snr) and blue_snr < MIN_BLUE_SNR)
            or not np.isfinite(fwhm_err_blue) or fwhm_err_blue > MAX_FWHM_ERR):
        flags.append('low_blue_signal')
    if good.sum() < MIN_WINDOWS:
        flags.append('few_windows')
    model.meta = dict(meta or {}, method='template', veil=veil, shape=shape,
                      n_windows=int(good.sum()), n_blue=int(blue.sum()), blue_snr=blue_snr,
                      fwhm_err_blue=fwhm_err_blue, flags=flags, ok=not flags)
    return (model, fits) if return_fits else model


# --- wavelength scale (Phase 2 task 4b) -----------------------------------------------

#: Pivot and half-range (nm) of the stretch: the scale error is
#: shift + stretch * (lambda - LAM0) / HALF, so ``stretch`` is the extra error
#: at 870 nm (and minus it at 390 nm), as in twin case (iv).
LAM0, HALF = 630.0, 240.0
#: Setup Q8: fall back to the calibration-period value above this (nm).
MAX_SCALE_ERR = 0.03
#: Gate G1: the srf method keeps 80 % of its gain within this (nm).
G1_TOLERANCE = 0.05


def _fit_window_shift(wav, spec, ref, lo, hi, sigma, err=None, via_grid=None):
    """Centroid shift (measured - reference, nm) in one window, with the SRF
    width fixed: spec = (c0 + c1 (x - xc)) R(x - dlam), R the reference
    through a Gaussian of ``sigma`` (evaluated on ``via_grid`` and
    resampled, if given)."""
    m = (wav >= lo) & (wav <= hi) & np.isfinite(spec)
    if err is not None:
        m &= np.isfinite(err) & (err > 0)
    x, y = wav[m], spec[m]
    e = err[m] if err is not None else None
    xc = 0.5 * (lo + hi)
    nan = dict(lo=lo, hi=hi, lam_center=xc, dlam=np.nan, dlam_err=np.nan, chi2_nu=np.nan,
               npix=int(x.size), ok=False)
    if x.size < 8 or ref[0][0] > lo - 6 or ref[0][-1] < hi + 6:
        return nan
    k0, k1 = np.searchsorted(ref[0], lo - 10.0), np.searchsorted(ref[0], hi + 10.0)
    rw, rf = ref[0][k0:k1], ref[1][k0:k1]
    g = None
    if via_grid is not None:
        g = via_grid[(via_grid > lo - 3.0) & (via_grid < hi + 3.0)]
        if g.size < 4:
            return nan

    def f(xx, c0, c1, dl):
        at = xx if g is None else g
        R = srfmod.convolve_gaussian(rw, rf, at - dl, sigma)
        R = R if g is None else np.interp(xx, g, R)
        return (c0 + c1 * (xx - xc)) * R

    R0 = f(x, 1.0, 0.0, 0.0)
    c0, c1 = np.linalg.lstsq(np.column_stack([R0, R0 * (x - xc)]), y, rcond=None)[0]
    try:
        with np.errstate(all='ignore'):
            p, cov = curve_fit(f, x, y, p0=[c0, c1, 0.0], sigma=e,
                               absolute_sigma=e is not None,
                               bounds=([-np.inf, -np.inf, -1.0], [np.inf, np.inf, 1.0]))
    except Exception:
        return nan
    if not np.all(np.isfinite(cov)) or abs(abs(p[2]) - 1.0) < 1e-6:
        return nan
    chi2 = float(np.sum(((y - f(x, *p)) / e) ** 2) / (x.size - 3)) if e is not None else np.nan
    return dict(lo=lo, hi=hi, lam_center=xc, dlam=float(p[2]),
                dlam_err=float(np.sqrt(cov[2, 2])), chi2_nu=chi2, npix=int(x.size), ok=True)


def fit_wavelength_scale(wav, spectrum, srf=None, ref=None, err=None, windows=None,
                         via_grid=None, pattern=None, return_fits=False):
    """Wavelength-scale error of a spectrum: shift and linear stretch.

    The true sampling wavelength of a pixel reported at lambda is::

        lambda + shift + stretch * (lambda - LAM0) / HALF

    (LAM0 = 630, HALF = 240 nm: ``stretch`` is the extra error at 870 nm, and
    minus it at 390 nm, as in twin case iv).  In each template window
    (:func:`hypernet.srf.template_windows`) the centroid offset dlam
    (measured - reference) is fitted, with the SRF width fixed by ``srf``
    (or free, via :func:`hypernet.srf.fit_srf_template`, if ``srf`` is None).
    Since a pixel that samples lambda + d shows the reference d early, the
    error is -dlam; a weighted straight line through -dlam(lambda) gives
    (shift, stretch) and their covariance, inflated by chi2_nu when the
    windows scatter more than their errors.

    Parameters
    ----------
    wav, spectrum : array
        Reported wavelengths (nm, air for the default reference) and spectrum.
    srf : float, callable or SRFModel, optional
        The channel's SRF FWHM (nm), held fixed in each window.
    ref : (ref_wave, ref_flux), optional
        High-resolution model; default TSIS-1 HSRS in air.  An Emod with gas
        absorption may be passed instead (the default windows avoid the
        telluric bands either way).
    err : array, optional
        Per-pixel 1-sigma error of ``spectrum`` (absolute weights).
    windows : list of (lo, hi), optional
        Default :func:`hypernet.srf.template_windows`.
    via_grid : array, optional
        The native grid the spectrum was observed on before resampling to
        ``wav`` (e.g. the E grid of an L2 irradiance); the shift then refers
        to ``via_grid``.
    pattern : dict or pandas.Series, optional
        Window centre (nm) -> a common offset (nm, in -dlam) to subtract
        before the straight-line fit, e.g. :func:`load_window_pattern`: the
        per-window structure shared by every HYPSTAR and channel, which is not
        an instrument scale error.  Windows without an entry are kept as is.
    return_fits : bool
        Also return the per-window table.

    Returns
    -------
    dict
        ``shift``, ``stretch`` (nm), ``cov`` (2x2, nm^2), ``shift_err``,
        ``stretch_err``, ``max_err`` (largest 1-sigma error of the correction
        over 390-870 nm), ``slope_per100nm`` (of dlam, as in Phase 0 8d),
        ``chi2_nu``, ``n_windows``, ``flags`` (``few_windows``, ``uncertain``
        when max_err > MAX_SCALE_ERR, ``fit_failed``) and ``ok``.
    pandas.DataFrame
        Only with ``return_fits``.
    """
    import pandas as pd

    wav = np.asarray(wav, dtype=float)
    spectrum = np.asarray(spectrum, dtype=float)
    if wav.ndim != 1 or spectrum.shape != wav.shape:
        raise ValueError('wav and spectrum must be 1-D and the same length')
    err = None if err is None else np.asarray(err, dtype=float)
    via_grid = None if via_grid is None else np.asarray(via_grid, dtype=float)
    ref = default_reference('air') if ref is None else ref
    windows = srfmod.template_windows() if windows is None else windows
    rows = []
    for lo, hi in windows:
        if srf is None:
            r = srfmod.fit_srf_template(wav, spectrum, ref[0], ref[1], lo, hi, err=err,
                                        veil=False, via_grid=via_grid)
            r = {k: r[k] for k in ('lo', 'hi', 'lam_center', 'dlam', 'dlam_err', 'chi2_nu',
                                   'npix', 'ok')}
        else:
            sig = float(srfmod.fwhm_values(srf, [0.5 * (lo + hi)])[0]) / srfmod.FWHM_PER_SIGMA
            r = _fit_window_shift(wav, spectrum, ref, lo, hi, sig, err=err, via_grid=via_grid)
        rows.append(r)
    fits = pd.DataFrame(rows)
    good = fits['ok'].astype(bool) & np.isfinite(fits['dlam_err']) & (fits['dlam_err'] > 0)
    out = dict(shift=np.nan, stretch=np.nan, cov=np.full((2, 2), np.nan), shift_err=np.nan,
               stretch_err=np.nan, max_err=np.nan, slope_per100nm=np.nan, chi2_nu=np.nan,
               n_windows=int(good.sum()))
    flags = []
    if good.sum() >= 3:
        g = fits[good]
        x = (g['lam_center'].values - LAM0) / HALF
        y = -g['dlam'].values
        if pattern is not None:
            pat = dict(pattern)
            y = y - np.array([pat.get(float(c), 0.0) for c in g['lam_center']])
        A = np.column_stack([np.ones_like(x), x])
        w = 1.0 / g['dlam_err'].values ** 2
        cov = np.linalg.inv(A.T @ (A * w[:, None]))
        coef = cov @ (A.T @ (w * y))
        chi2 = float(np.sum(w * (y - A @ coef) ** 2) / (len(g) - 2))
        cov = cov * max(chi2, 1.0)
        xx = (np.array([390.0, 870.0]) - LAM0) / HALF
        var = cov[0, 0] + xx ** 2 * cov[1, 1] + 2 * xx * cov[0, 1]
        out.update(shift=float(coef[0]), stretch=float(coef[1]), cov=cov,
                   shift_err=float(np.sqrt(cov[0, 0])), stretch_err=float(np.sqrt(cov[1, 1])),
                   max_err=float(np.sqrt(var.max())),
                   slope_per100nm=float(-coef[1] * 100.0 / HALF), chi2_nu=chi2)
    else:
        flags.append('fit_failed')
    if good.sum() < MIN_WINDOWS:
        flags.append('few_windows')
    if not np.isfinite(out['max_err']) or out['max_err'] > MAX_SCALE_ERR:
        flags.append('uncertain')
    out.update(flags=flags, ok=not flags)
    return (out, fits) if return_fits else out


def load_window_pattern(path=None):
    """The common per-window offset pattern (window centre -> nm, in -dlam),
    ``hypernet/data/wavescale_window_pattern.csv`` (written by
    ``hypernet/wiggles/phase2_wavescale.py``)."""
    import os

    import pandas as pd
    path = path or os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data',
                                'wavescale_window_pattern.csv')
    t = pd.read_csv(path)
    return dict(zip(t['lam_center'].astype(float), t['pattern'].astype(float)))


def correct_wavelengths(wav, scale):
    """True sampling wavelengths, ``wav + shift + stretch (wav - LAM0) / HALF``,
    from a :func:`fit_wavelength_scale` result (or any dict with ``shift`` and
    ``stretch``).  Apply to E's grid before :func:`hypernet.edinterp.interpolate_ed_to_l`."""
    wav = np.asarray(wav, dtype=float)
    return wav + scale['shift'] + scale['stretch'] * (wav - LAM0) / HALF


def scale_or_fallback(scale, fallback):
    """Setup Q8: the per-sequence ``scale`` if it is ``ok`` (max_err <=
    MAX_SCALE_ERR and enough windows), else ``fallback`` (e.g. the
    calibration-period median).  Returns (scale, source)."""
    if scale is not None and scale.get('ok'):
        return scale, 'sequence'
    return fallback, 'fallback'

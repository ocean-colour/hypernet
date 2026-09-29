"""Spectral response functions (SRFs) of hyperspectral radiometers from solar
and telluric lines.

The SRF of each channel is parameterised as a Gaussian whose FWHM varies
smoothly with wavelength.  This module measures it from the data: every line
in :data:`LINES` is fitted with a Gaussian absorption profile on a linear
continuum, which gives a centroid (wavelength calibration) and a width
(resolution) per line.  :func:`fit_fwhm_model` then fits a quadratic
FWHM(lambda) to the line widths of one channel, and :class:`SRFModel` holds
the result (with the channel's mean centroid offset) and reads and writes it
as JSON.

A fitted Gaussian measures the SRF convolved with the line's intrinsic
profile, not the SRF itself.  Blends and bands (Ca H/K, G band, Mg b, Na D,
O2) therefore give widths biased high.  The bias is common to the E and L
channels of one instrument, so the empirical fits here are good for relative
tests (E vs L, stability).  An absolute SRF needs a forward-model fit against
a high-resolution solar spectrum.

Wavelengths are in nm.  Laboratory wavelengths are given in air, with vacuum
values from :func:`air_to_vac`; which scale the HYPSTAR calibration uses is
not yet known.
"""

import dataclasses
import json
import os

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

#: FWHM / sigma for a Gaussian, 2 sqrt(2 ln 2).
FWHM_PER_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))


def air_to_vac(lam_air):
    """Convert air wavelengths (nm) to vacuum (nm).

    Uses the inverse of the Morton (2000) / Ciddor (1996) dispersion relation,
    as in VALD (N. Piskunov).  At 656.281 nm (air) this gives 656.461 nm.

    Args:
        lam_air (float or array): air wavelength(s) in nm.

    Returns:
        float or np.ndarray: vacuum wavelength(s) in nm.
    """
    lam_air = np.asarray(lam_air, dtype=float)
    s2 = (1e3 / lam_air) ** 2  # (1 / lambda[um])^2
    n = (1.0 + 8.336624212083e-5 + 2.408926869968e-2 / (130.1065924522 - s2)
         + 1.599740894897e-4 / (38.92568793293 - s2))
    return lam_air * n


# The line list.  Columns:
#   name        -- label
#   lam_air     -- laboratory (or nominal band) wavelength in air, nm
#   half        -- half-width of the fit window, nm
#   group       -- lines sharing a group are fitted jointly on one window
#   blend       -- the feature is itself a blend or band, so neither its
#                  absolute centroid nor its width is a clean SRF probe
#   use_for_srf -- enters the FWHM(lambda) fit
_LINE_ROWS = [
    # name          lam_air   half  group    blend  use_for_srf
    ('Ca K',        393.366,  4.0, 'CaHK',   True,  True),
    ('Ca H',        396.847,  4.0, 'CaHK',   True,  True),   # + H epsilon 397.0
    ('G band',      430.79,   4.0, 'G',      True,  True),   # CH band
    ('H beta',      486.133,  4.0, 'Hb',     False, True),
    ('Mg b',        517.27,   4.0, 'Mgb',    True,  True),   # 516.7/517.3/518.4
    ('Na D',        589.29,   4.0, 'NaD',    True,  True),   # 589.0/589.6
    ('H alpha',     656.281,  4.0, 'Ha',     False, True),
    ('O2 B',        687.0,    4.0, 'O2B',    True,  False),  # band
    ('O2 A',        760.6,    4.0, 'O2A',    True,  False),  # band
    ('Ca II 849.8', 849.802,  4.0, 'CaIR12', False, True),   # jointly with 854.2
    ('Ca II 854.2', 854.209,  4.0, 'CaIR12', False, True),
    ('Ca II 866.2', 866.214,  4.0, 'CaIR3',  False, True),
    ('H2O',         936.5,    4.0, 'H2O',    True,  False),  # band, diagnostic
]

#: The line table (see the comments above ``_LINE_ROWS`` for the columns).
LINES = pd.DataFrame(_LINE_ROWS, columns=['name', 'lam_air', 'half', 'group',
                                          'blend', 'use_for_srf'])
LINES.insert(2, 'lam_vac', air_to_vac(LINES['lam_air'].values))

#: Keys of one fitted component, in output order.
_FIT_KEYS = ['mu', 'mu_err', 'sigma', 'sigma_err', 'fwhm', 'fwhm_err',
             'depth', 'depth_err', 'c0', 'c0_err', 'c1', 'c1_err',
             'chi2_nu', 'npix', 'ok']


def flatten_scans(scans, axis=1, smooth_px=20.0):
    """Remove a smooth multiplicative change of each scan relative to the mean.

    Each scan is divided by a Gaussian-smoothed (``smooth_px`` pixels, along
    wavelength) version of its ratio to the mean spectrum.  This takes out
    broadband changes between scans (sky brightness, glint, pointing) while
    keeping pixel-to-pixel noise and line-scale changes (a few pixels).

    Args:
        scans (np.ndarray): spectra, with scans along ``axis``.
        axis (int): the scan axis.
        smooth_px (float): Gaussian sigma of the smoothing, pixels.

    Returns:
        np.ndarray: the flattened scans, same shape as ``scans``.
    """
    from scipy.ndimage import gaussian_filter1d
    scans = np.moveaxis(np.asarray(scans, dtype=float), axis, -1)
    with np.errstate(invalid='ignore', divide='ignore'):
        ratio = scans / np.nanmean(scans, axis=-1, keepdims=True)
    ratio = np.where(np.isfinite(ratio), ratio, 1.0)
    smooth = gaussian_filter1d(ratio, smooth_px, axis=0, mode='nearest')
    return np.moveaxis(scans / smooth, -1, axis)


def scan_errors(scans, axis=1, flatten_px=None):
    """Mean spectrum and its per-pixel uncertainty from repeated scans.

    L1A files carry no uncertainty variables (and the L1C ``u_rel_random_*``
    are zero), so the scan-to-scan scatter is the noise estimate.  With ``s``
    the sample standard deviation (``ddof=1``) over N scans:

    - the error of the mean spectrum is ``s / sqrt(N)``;
    - the error of a single scan is ``s``.

    The raw scatter is dominated by broadband changes between scans (in the
    VEIT sample 5-6 % in Lu and Ld against 0.4-0.8 % pixel noise; see
    ``wiggles/phase0a_l1c_consistency.py``).  With ``flatten_px`` the scatter
    is taken from :func:`flatten_scans` instead, which is the right noise for
    a line fit on a free continuum.  The mean is always that of the raw scans.

    Args:
        scans (np.ndarray): spectra, with scans along ``axis``.
        axis (int): the scan axis.  L1A arrays are (wavelength, scan), hence 1.
        flatten_px (float, optional): smoothing sigma (pixels) for
            :func:`flatten_scans`; None uses the raw scatter.

    Returns:
        tuple: (mean, err_mean, err_scan), each along wavelength.
    """
    scans = np.asarray(scans, dtype=float)
    n = np.sum(np.isfinite(scans), axis=axis)
    mean = np.nanmean(scans, axis=axis)
    work = scans if flatten_px is None else flatten_scans(scans, axis, flatten_px)
    with np.errstate(invalid='ignore', divide='ignore'):
        s = np.nanstd(work, axis=axis, ddof=1)
        return mean, s / np.sqrt(n), s


def _model(x, xref, ncomp, *p):
    """Linear continuum times ``ncomp`` Gaussian absorption profiles.

    ``p = (c0, c1, a_1, mu_1, s_1, ..., a_n, mu_n, s_n)``; the continuum is
    ``c0 + c1 (x - xref)``, so ``c0`` is its value at ``xref``.
    """
    y = p[0] + p[1] * (x - xref)
    for k in range(ncomp):
        a, mu, s = p[2 + 3 * k: 5 + 3 * k]
        y = y * (1.0 - a * np.exp(-0.5 * ((x - mu) / s) ** 2))
    return y


def _failed(npix=0):
    out = {k: np.nan for k in _FIT_KEYS}
    out['npix'] = npix
    out['ok'] = False
    return out


def fit_blend(wav, spec, lams, half=4.0, err=None, sigma0=1.3):
    """Fit one or more Gaussian absorption lines jointly on a linear continuum.

    The window is ``[min(lams) - half, max(lams) + half]``.  With ``err`` the
    fit is weighted and ``absolute_sigma=True``, so the parameter errors are
    the propagated measurement errors and ``chi2_nu`` tests the model.
    Without it, ``curve_fit`` rescales the covariance by the residuals and
    ``chi2_nu`` is NaN.

    A fit that fails (too few pixels, no convergence, non-finite covariance,
    or a result outside the window, with depth outside (0, 1) or width
    outside (0.05, half)) returns NaN for every component, with ``ok=False``.
    It never raises.

    Args:
        wav (np.ndarray): wavelength grid, nm.
        spec (np.ndarray): spectrum on ``wav``.
        lams (sequence of float): nominal line centres, nm.
        half (float): half-width of the window beyond the outer lines, nm.
        err (np.ndarray, optional): 1-sigma per-pixel uncertainty of ``spec``.
        sigma0 (float): initial Gaussian sigma, nm.

    Returns:
        list of dict: one per line, keys ``mu, mu_err, sigma, sigma_err, fwhm,
        fwhm_err, depth, depth_err, c0, c0_err, c1, c1_err, chi2_nu, npix, ok``
        (continuum and ``chi2_nu`` are shared by all lines in the window).
    """
    lams = np.atleast_1d(np.asarray(lams, dtype=float))
    ncomp = lams.size
    wav = np.asarray(wav, dtype=float)
    spec = np.asarray(spec, dtype=float)
    lo, hi = lams.min() - half, lams.max() + half
    m = (wav > lo) & (wav < hi) & np.isfinite(spec)
    if err is not None:
        err = np.asarray(err, dtype=float)
        m &= np.isfinite(err) & (err > 0)
    x, y = wav[m], spec[m]
    e = err[m] if err is not None else None
    npar = 2 + 3 * ncomp
    npix = int(x.size)
    if npix < npar + 3:
        return [_failed(npix) for _ in range(ncomp)]

    # Initial guesses: continuum through the window edges, depth from the
    # spectrum at each nominal centre.
    nedge = max(2, npix // 10)
    xl, yl = np.mean(x[:nedge]), np.median(y[:nedge])
    xr, yr = np.mean(x[-nedge:]), np.median(y[-nedge:])
    xref = 0.5 * (lo + hi)
    c1 = (yr - yl) / (xr - xl) if xr > xl else 0.0
    c0 = yl + c1 * (xref - xl)
    p0 = [c0, c1]
    for lam in lams:
        cont = c0 + c1 * (lam - xref)
        ymin = y[np.argmin(np.abs(x - lam))]
        a0 = np.clip(1.0 - ymin / cont, 0.05, 0.9) if cont > 0 else 0.2
        p0 += [a0, lam, sigma0]

    def f(xx, *p):
        return _model(xx, xref, ncomp, *p)

    try:
        with np.errstate(all='ignore'):
            p, cov = curve_fit(f, x, y, p0=p0, sigma=e,
                               absolute_sigma=e is not None, maxfev=20000)
    except Exception:
        return [_failed(npix) for _ in range(ncomp)]
    perr = np.sqrt(np.diag(cov)) if np.all(np.isfinite(cov)) else None
    if perr is None or not np.all(np.isfinite(p)) or not np.all(np.isfinite(perr)):
        return [_failed(npix) for _ in range(ncomp)]

    chi2_nu = np.nan
    if e is not None:
        chi2_nu = float(np.sum(((y - f(x, *p)) / e) ** 2) / (npix - npar))
    # c0 is reported at the nominal line centre, not at the window centre.
    out = []
    for k, lam in enumerate(lams):
        a, mu, s = p[2 + 3 * k: 5 + 3 * k]
        ea, emu, es = perr[2 + 3 * k: 5 + 3 * k]
        s = abs(s)
        if not (lo < mu < hi and 0.05 < s < half and 0.0 < a < 1.0):
            return [_failed(npix) for _ in range(ncomp)]
        c0_lam = p[0] + p[1] * (lam - xref)
        var_c0 = (cov[0, 0] + (lam - xref) ** 2 * cov[1, 1]
                  + 2 * (lam - xref) * cov[0, 1])
        out.append(dict(mu=mu, mu_err=emu, sigma=s, sigma_err=es,
                        fwhm=FWHM_PER_SIGMA * s, fwhm_err=FWHM_PER_SIGMA * es,
                        depth=a, depth_err=ea,
                        c0=c0_lam, c0_err=np.sqrt(max(var_c0, 0.0)),
                        c1=p[1], c1_err=perr[1],
                        chi2_nu=chi2_nu, npix=npix, ok=True))
    return out


def fit_line(wav, spec, lam, half=4.0, err=None, sigma0=1.3):
    """Fit a single Gaussian absorption line on a linear continuum.

    See :func:`fit_blend` for the model, weighting and failure behaviour.

    Args:
        wav (np.ndarray): wavelength grid, nm.
        spec (np.ndarray): spectrum on ``wav``.
        lam (float): nominal line centre, nm.
        half (float): half-width of the fit window, nm.
        err (np.ndarray, optional): 1-sigma per-pixel uncertainty of ``spec``.
        sigma0 (float): initial Gaussian sigma, nm.

    Returns:
        dict: ``mu, mu_err, sigma, sigma_err, fwhm, fwhm_err, depth,
        depth_err, c0, c0_err, c1, c1_err, chi2_nu, npix, ok``.
    """
    return fit_blend(wav, spec, [lam], half=half, err=err, sigma0=sigma0)[0]


def fit_lines(wav, spec, lines=LINES, err=None, lam_col='lam_air'):
    """Fit every line of a line table; lines sharing a ``group`` jointly.

    Args:
        wav (np.ndarray): wavelength grid, nm.
        spec (np.ndarray): spectrum on ``wav``.
        lines (pd.DataFrame): line table with the columns of :data:`LINES`.
        err (np.ndarray, optional): 1-sigma per-pixel uncertainty of ``spec``.
        lam_col (str): column of ``lines`` taken as the nominal centre
            (``'lam_air'`` or ``'lam_vac'``).

    Returns:
        pd.DataFrame: the rows of ``lines`` (same order) with the fit columns
        of :func:`fit_blend` appended, plus ``dmu = mu - lines[lam_col]``.
        Lines that fail, or fall off the grid, are NaN with ``ok=False``.
    """
    lines = lines.reset_index(drop=True)
    fits = [None] * len(lines)
    for _, idx in lines.groupby('group', sort=False).groups.items():
        sub = lines.loc[idx]
        res = fit_blend(wav, spec, sub[lam_col].values,
                        half=float(sub['half'].max()), err=err)
        for i, r in zip(idx, res):
            fits[i] = r
    out = pd.concat([lines, pd.DataFrame(fits, columns=_FIT_KEYS)], axis=1)
    out['dmu'] = out['mu'] - out[lam_col]
    out['ok'] = out['ok'].astype(bool)
    return out


# --- FWHM(lambda) model -------------------------------------------------------

#: Default reference wavelength and scale of the FWHM polynomial, nm.  The
#: polynomial is in x = (lam - LAM_REF) / LAM_SCALE, so coeffs[0] is the FWHM
#: at LAM_REF and the coefficients are well conditioned.
LAM_REF = 600.0
LAM_SCALE = 100.0


def _design(lam, deg, lam_ref, lam_scale):
    x = (np.asarray(lam, dtype=float) - lam_ref) / lam_scale
    return np.vander(np.atleast_1d(x), deg + 1, increasing=True)


def fit_fwhm_model(lam, fwhm, err, deg=2, lam_ref=LAM_REF, lam_scale=LAM_SCALE,
                   scale_cov=False):
    """Weighted least-squares polynomial (default quadratic) FWHM(lambda).

    ``FWHM = sum_k coeffs[k] x**k`` with ``x = (lam - lam_ref) / lam_scale``.
    The covariance is ``(A^T W A)^-1`` with ``W = 1/err**2``, i.e. it takes
    ``err`` as absolute 1-sigma errors.  Points with a non-finite value or a
    non-positive error are dropped.  With fewer than ``deg + 1`` points left,
    the coefficients and covariance are NaN.

    Args:
        lam (array): line wavelengths, nm.
        fwhm (array): fitted FWHM at each line, nm.
        err (array): 1-sigma error of each FWHM, nm.
        deg (int): polynomial degree.
        lam_ref (float): reference wavelength, nm.
        lam_scale (float): wavelength scale, nm.
        scale_cov (bool): multiply the covariance by ``chi2_nu`` when it
            exceeds 1 (when the scatter about the model exceeds ``err``).

    Returns:
        tuple: (coeffs, cov) -- arrays of shape (deg+1,) and (deg+1, deg+1).
    """
    lam, fwhm, err = (np.asarray(a, dtype=float).ravel() for a in (lam, fwhm, err))
    m = np.isfinite(lam) & np.isfinite(fwhm) & np.isfinite(err) & (err > 0)
    npar = deg + 1
    if m.sum() < npar:
        return np.full(npar, np.nan), np.full((npar, npar), np.nan)
    A = _design(lam[m], deg, lam_ref, lam_scale)
    w = 1.0 / err[m] ** 2
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    coeffs = cov @ (A.T @ (w * fwhm[m]))
    if scale_cov and m.sum() > npar:
        chi2_nu = np.sum(w * (fwhm[m] - A @ coeffs) ** 2) / (m.sum() - npar)
        cov = cov * max(chi2_nu, 1.0)
    return coeffs, cov


def fwhm_at(lam, coeffs, lam_ref=LAM_REF, lam_scale=LAM_SCALE):
    """Evaluate the FWHM(lambda) polynomial (nm) at ``lam`` (nm)."""
    coeffs = np.asarray(coeffs, dtype=float)
    y = _design(lam, coeffs.size - 1, lam_ref, lam_scale) @ coeffs
    return y if np.ndim(lam) else float(y[0])


def fwhm_err_at(lam, cov, lam_ref=LAM_REF, lam_scale=LAM_SCALE):
    """1-sigma error (nm) of the FWHM(lambda) polynomial at ``lam`` (nm)."""
    cov = np.asarray(cov, dtype=float)
    A = _design(lam, cov.shape[0] - 1, lam_ref, lam_scale)
    y = np.sqrt(np.einsum('ij,jk,ik->i', A, cov, A))
    return y if np.ndim(lam) else float(y[0])


def _nan_to_none(v):
    """Replace NaN by None, recursively, so the JSON is standard."""
    if isinstance(v, (list, tuple)):
        return [_nan_to_none(u) for u in v]
    if isinstance(v, float) and not np.isfinite(v):
        return None
    return v


def _none_to_nan(v):
    if isinstance(v, list):
        return [_none_to_nan(u) for u in v]
    return np.nan if v is None else v


@dataclasses.dataclass
class SRFModel:
    """Gaussian SRF of one channel: FWHM(lambda) and a centroid offset.

    Attributes:
        channel (str): e.g. ``'E'``, ``'Ld'``, ``'Lu'``.
        coeffs (np.ndarray): FWHM polynomial coefficients in
            ``x = (lam - lam_ref) / lam_scale``, nm.
        cov (np.ndarray): covariance of ``coeffs``, nm^2.
        offset (float): mean centroid offset, measured minus laboratory, nm.
        offset_err (float): its 1-sigma error, nm.
        lam_min, lam_max (float): wavelength range of the lines used, nm.
        lam_ref, lam_scale (float): polynomial reference and scale, nm.
        chi2_nu (float): reduced chi^2 of the FWHM fit.
        npts (int): number of lines in the FWHM fit.
        instrument (str): instrument id (e.g. ``'HYPSTAR_122304'``).
        frame (str): wavelength scale of ``offset`` (``'air'`` or ``'vac'``).
        meta (dict): free-form provenance (sequence, method, ...).
    """
    channel: str
    coeffs: np.ndarray
    cov: np.ndarray
    offset: float = np.nan
    offset_err: float = np.nan
    lam_min: float = np.nan
    lam_max: float = np.nan
    lam_ref: float = LAM_REF
    lam_scale: float = LAM_SCALE
    chi2_nu: float = np.nan
    npts: int = 0
    instrument: str = ''
    frame: str = 'air'
    meta: dict = dataclasses.field(default_factory=dict)

    def __post_init__(self):
        self.coeffs = np.asarray(self.coeffs, dtype=float)
        self.cov = np.asarray(self.cov, dtype=float)

    def fwhm(self, lam):
        """FWHM (nm) at ``lam`` (nm)."""
        return fwhm_at(lam, self.coeffs, self.lam_ref, self.lam_scale)

    def fwhm_err(self, lam):
        """1-sigma error of the FWHM (nm) at ``lam`` (nm)."""
        return fwhm_err_at(lam, self.cov, self.lam_ref, self.lam_scale)

    def sigma(self, lam):
        """Gaussian sigma (nm) at ``lam`` (nm)."""
        return self.fwhm(lam) / FWHM_PER_SIGMA

    @classmethod
    def from_lines(cls, channel, fits, deg=2, lam_col='lam_air',
                   scale_cov=False, **kwargs):
        """Build a model from a :func:`fit_lines` table.

        The FWHM(lambda) fit uses the rows with ``ok`` and ``use_for_srf``.
        The offset is the weighted mean of ``dmu`` over the rows that are
        also not ``blend`` (the clean lines).  Its error is the error of the
        weighted mean, inflated by ``sqrt(chi2_nu)`` when the scatter exceeds
        the errors.

        Args:
            channel (str): channel name.
            fits (pd.DataFrame): output of :func:`fit_lines`.
            deg (int): polynomial degree.
            lam_col (str): column used as the line wavelength.
            scale_cov (bool): passed to :func:`fit_fwhm_model`.
            **kwargs: other :class:`SRFModel` fields (instrument, meta, ...).

        Returns:
            SRFModel
        """
        use = fits['ok'].astype(bool) & fits['use_for_srf'].astype(bool)
        f = fits[use]
        lam_ref = kwargs.pop('lam_ref', LAM_REF)
        lam_scale = kwargs.pop('lam_scale', LAM_SCALE)
        coeffs, cov = fit_fwhm_model(f[lam_col], f['fwhm'], f['fwhm_err'],
                                     deg=deg, lam_ref=lam_ref,
                                     lam_scale=lam_scale, scale_cov=scale_cov)
        chi2_nu = np.nan
        if len(f) > deg + 1 and np.all(np.isfinite(coeffs)):
            r = (f['fwhm'] - fwhm_at(f[lam_col].values, coeffs, lam_ref,
                                     lam_scale)) / f['fwhm_err']
            chi2_nu = float(np.sum(r ** 2) / (len(f) - deg - 1))
        clean = f[~f['blend'].astype(bool)]
        clean = clean[np.isfinite(clean['mu_err']) & (clean['mu_err'] > 0)]
        offset = offset_err = np.nan
        if len(clean):
            w = 1.0 / clean['mu_err'].values ** 2
            offset = float(np.sum(w * clean['dmu']) / np.sum(w))
            offset_err = float(1.0 / np.sqrt(np.sum(w)))
            if len(clean) > 1:
                chi2_off = np.sum(w * (clean['dmu'] - offset) ** 2) / (len(clean) - 1)
                offset_err *= np.sqrt(max(chi2_off, 1.0))
        return cls(channel=channel, coeffs=coeffs, cov=cov, offset=offset,
                   offset_err=offset_err,
                   lam_min=float(f[lam_col].min()) if len(f) else np.nan,
                   lam_max=float(f[lam_col].max()) if len(f) else np.nan,
                   lam_ref=lam_ref, lam_scale=lam_scale, chi2_nu=chi2_nu,
                   npts=int(len(f)),
                   frame='vac' if lam_col == 'lam_vac' else 'air', **kwargs)

    def to_dict(self):
        """Plain-Python dict (lists, NaN as None)."""
        d = dataclasses.asdict(self)
        d['coeffs'] = self.coeffs.tolist()
        d['cov'] = self.cov.tolist()
        return {k: _nan_to_none(v) for k, v in d.items()}

    @classmethod
    def from_dict(cls, d):
        d = {k: _none_to_nan(v) for k, v in d.items()}
        return cls(**d)

    def to_json(self, path=None, indent=2):
        """JSON string; also written to ``path`` if given."""
        s = json.dumps(self.to_dict(), indent=indent)
        if path is not None:
            with open(path, 'w') as fh:
                fh.write(s + '\n')
        return s

    @classmethod
    def from_json(cls, s):
        """From a JSON string or the path of a JSON file."""
        if os.path.exists(str(s)):
            with open(s) as fh:
                s = fh.read()
        return cls.from_dict(json.loads(s))


def save_srf_models(path, models, meta=None, indent=2):
    """Write several :class:`SRFModel` to one JSON file, keyed by channel."""
    doc = {'meta': meta or {},
           'models': {m.channel: m.to_dict() for m in models}}
    with open(path, 'w') as fh:
        fh.write(json.dumps(doc, indent=indent) + '\n')


def load_srf_models(path):
    """Read a file written by :func:`save_srf_models`.

    Returns:
        dict: channel -> :class:`SRFModel`.
    """
    with open(path) as fh:
        doc = json.load(fh)
    return {k: SRFModel.from_dict(v) for k, v in doc['models'].items()}

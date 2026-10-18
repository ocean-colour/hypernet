"""HYPSTAR laboratory wavelength calibration (Joel's lamp data, 2026-10-17).

The lab tree is ``$OS_COLOR/hypernet/wiggles/lab_data/hypstar_<sn>/calibration/
wavelength/<epoch>/`` (:func:`lab_root`; the zip is never read).  Each of the
13 epochs (7 instruments, VNIR only) holds four files:

- ``hypstar_<sn>_wl_coefs_<yymmdd>.dat``: 4th-order pixel -> lambda
  polynomials (``c0 + c1 px + ... + c4 px^4``, px from 0 over 2048 px) for
  VNIR_L and VNIR_E, and the residual sigma (2 decimals).  The SWIR columns are
  zero.  The yymmdd in the name is *not* unique (121222 ``2021_07`` and
  ``2021_07_tst`` share a name but not coefficients): the epoch directory is
  the key.  :func:`read_coefs`.
- ``fwhm_vnir.dat``: lab FWHM_L and FWHM_E (nm) at the reference lines.
- ``deltawl_vnir.dat``: reference lambda minus the nominal-grid lambda at the
  fitted line centre.  Both are 4 columns (``wl L wl E``; the two ``wl``
  columns are identical) and the lines are not sorted.  :func:`read_line_file`.
- ``lamp_vnir.dat``: 2048 rows; a nominal per-instrument ``wl`` column (0.1 nm
  rounding, the same at every epoch) and 2-3 peak-normalised L/E column pairs,
  one per lamp (Hg-Ar and Xe; Kr in place of Xe on 120242 2020_09; Hg-Ar, Kr
  and Xe on 121222/121231 2021_07).  :func:`read_lamp`.

All reference lines are NIST **air** wavelengths.

Committed products (``hypernet/data/``), written by
``python -m hypernet.wiggles.lab_tables``:

- ``lab_wl_coefs.csv`` -- one row per instrument x epoch x channel (26 rows),
  coefficients at full precision (Python ``repr``, the shortest string that
  round-trips exactly; c4 ~ 3e-12 multiplies px^4 ~ 1.8e13, so a short format
  would break the 1e-4 nm grid tolerance).  :func:`load_wl_coefs`.
- ``lab_srf.csv`` -- one row per instrument x epoch x channel x line: the lab
  FWHM and deltawl plus our refits of the lamp lines.  :func:`load_lab_srf`.

Grids.  The processor's wavelength grid is the lab polynomial of
:data:`PROCESSOR_EPOCH` on the pixels with 350 <= lambda < 1100 nm
(:func:`grid_for`), with one exception (:data:`CROP_EPOCH`): on 120242 the
pixel window is the one selected by that instrument's 2020_09 polynomial,
while the wavelengths are 2024_05's.

Line fits (:func:`fit_line_px`, :func:`fit_epoch_lines`).  Each lamp line is
fitted in pixel space over +-7 px with a Gaussian, a super-Gaussian
(exponent p; same FWHM convention as :func:`hypernet.srf.srf_kernel`) and an
asymmetric Gaussian (blue/red widths sigma (1 -+ asym)), plus a constant.
FWHMs are converted to nm with the local dispersion of the epoch's own
polynomial.  The noise is a robust std of line-free parts of the lamp column
(450-530 and 600-680 nm on the nominal ``wl`` column), which misses the
photon noise of the bright lines: chi2_nu is 30-560 on Hg but 0.1-0.2 on Ar.
**Every ``*_err`` column is therefore the formal error multiplied by
sqrt(chi2_nu) of its own fit** (the same as ``curve_fit`` with
``absolute_sigma=False``).
"""
import glob
import os
import re

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

FWHM_PER_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))
NPIX = 2048
LAM_LO, LAM_HI = 350.0, 1100.0
CHANNELS = ('L', 'E')
INSTRUMENTS = ('120242', '121222', '121231', '122302', '122303', '122304', '122305')

#: The epoch whose polynomial the processor (Release 2) uses, one per
#: instrument, and how that is known: ``verified`` = a processor grid of that
#: channel is reproduced exactly (VEIT L1A E and L for 122304; the Release 2
#: L grid of every instrument); ``inferred_L`` = only the L grid was checked
#: and E is assumed to come from the same epoch.
PROCESSOR_EPOCH = {
    '120242': dict(epoch='2024_05', L='verified', E='inferred_L'),
    '121222': dict(epoch='2021_07', L='verified', E='inferred_L'),
    '121231': dict(epoch='2021_07', L='verified', E='inferred_L'),
    '122302': dict(epoch='2023_01', L='verified', E='inferred_L'),
    '122303': dict(epoch='2023_01', L='verified', E='inferred_L'),
    '122304': dict(epoch='2023_01', L='verified', E='verified'),
    '122305': dict(epoch='2023_01', L='verified', E='inferred_L'),
}

#: Per-instrument exception to the crop rule (Setup Q9).  The Release 2 L grid
#: of 120242 has 1540 px (px 411-1950); 350 <= lambda < 1100 on its 2024_05
#: polynomial gives 1541 (px 1951 = 1099.575 nm).  The window that the 2020_09
#: polynomial selects (px 1951 = 1100.30 nm) reproduces it exactly, so the
#: processor's pixel window was evidently fixed at the instrument's first
#: calibration (none of the 2022_05, 2022_06 or 2024_05 windows does).
#: Applied to E as well, which is unverified (``inferred_L``): the 2020_09 E
#: window is px 412-1950, 1539 px, against 1541 by the plain rule.  Checked by
#: ``hypernet.wiggles.lab_tables`` (``crop_hypothesis_120242.csv``).
CROP_EPOCH = {'120242': '2020_09'}

#: Species of the reference lines (identified by value; NIST air).
SPECIES = {296.728: 'Hg I', 404.6563: 'Hg I', 435.8328: 'Hg I', 546.0735: 'Hg I',
           1013.976: 'Hg I', 696.5431: 'Ar I', 763.5106: 'Ar I', 826.4522: 'Ar I',
           823.1635: 'Xe I', 881.9411: 'Xe I', 904.545: 'Xe I', 992.319: 'Xe I',
           1083.837: 'Xe I', 587.0916: 'Kr I', 785.48234: 'Kr I', 850.88728: 'Kr I',
           877.67505: 'Kr I', 892.86934: 'Kr I', 975.1761: 'Kr I'}
#: Lines whose lab width is raised by a neighbour.
BLENDS = (404.6563, 823.1635)
#: Below this line S/N a row is flagged ``low_snr``.
SNR_MIN = 10.0
#: The brightest line of each lamp, used to tell the lamp pairs apart.
LAMP_MARKERS = (('HgAr', 435.8328), ('Kr', 785.48234), ('Xe', 881.9411))

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
COEF_COLUMNS = ['instrument', 'sn', 'epoch', 'channel', 'c0', 'c1', 'c2', 'c3', 'c4',
                'sigma', 'coef_file', 'coef_date', 'file_created', 'lamps', 'n_lines',
                'lines', 'lines_source', 'processor', 'processor_status']
SRF_COLUMNS = ['instrument', 'sn', 'epoch', 'channel', 'lam_ref', 'species', 'lamp',
               'lab_fwhm', 'lab_deltawl', 'mu_px', 'mu_px_err', 'disp_nm_px', 'lam_fit',
               'resid', 'fwhm_gauss', 'fwhm_gauss_err', 'fwhm_halfmax', 'fwhm_sg',
               'fwhm_sg_err', 'sg_p', 'sg_p_err', 'asym', 'asym_err', 'chi2_gauss',
               'chi2_sg', 'chi2_asym', 'snr', 'blend', 'low_snr', 'outside_grid']


# ----------------------------------------------------------------------------
# the lab tree
# ----------------------------------------------------------------------------
def lab_root():
    """``$OS_COLOR/hypernet/wiggles/lab_data`` (the unzipped lab tree)."""
    return os.path.join(os.environ.get('OS_COLOR', ''), 'hypernet', 'wiggles', 'lab_data')


def lab_available(root=None):
    return bool(glob.glob(os.path.join(root or lab_root(), 'hypstar_*', 'calibration',
                                       'wavelength', '*')))


def epochs(root=None):
    """All ``(sn, epoch, directory)`` in the lab tree, sorted."""
    out = []
    for d in sorted(glob.glob(os.path.join(root or lab_root(), 'hypstar_*', 'calibration',
                                           'wavelength', '*'))):
        if os.path.isdir(d):
            out.append((re.search(r'hypstar_(\d+)', d).group(1), os.path.basename(d), d))
    return out


def epoch_dir(sn, epoch, root=None):
    return os.path.join(root or lab_root(), 'hypstar_%s' % sn, 'calibration', 'wavelength',
                        epoch)


def _sn(instrument):
    s = str(instrument)
    return s.split('_')[-1] if s.startswith('HYPSTAR') else s


# ----------------------------------------------------------------------------
# readers
# ----------------------------------------------------------------------------
def read_coefs(d):
    """The coefficient file of epoch directory ``d``.

    Returns:
        dict: ``L``, ``E`` (arrays c0..c4), ``sigma_L``, ``sigma_E``,
        ``coef_file``, ``coef_date`` (Timestamp from the yymmdd in the name),
        ``file_created`` (str from the header, '' if absent).
    """
    f = glob.glob(os.path.join(d, 'hypstar_*_wl_coefs_*.dat'))
    if len(f) != 1:
        raise FileNotFoundError('expected one wl_coefs file in %s, found %d' % (d, len(f)))
    f = f[0]
    ymd = re.search(r'coefs_(\d{6})', f).group(1)
    sig, created = {}, ''
    for line in open(f):
        m = re.match(r'#\s*sigma_(\w+)\s+([\d.]+)', line)
        if m:
            sig[m.group(1)] = float(m.group(2))
        m = re.match(r'#\s*File was created\s+(.+?)\s*$', line)
        if m:
            created = m.group(1)
    a = np.loadtxt(f, comments='#')
    if a.shape != (5, 4):
        raise ValueError('%s: expected 5 x 4 coefficients, got %s' % (f, a.shape))
    return dict(L=a[:, 0].copy(), E=a[:, 1].copy(), sigma_L=sig.get('VNIR_L', np.nan),
                sigma_E=sig.get('VNIR_E', np.nan), coef_file=os.path.basename(f),
                coef_date=pd.Timestamp('20%s-%s-%s' % (ymd[:2], ymd[2:4], ymd[4:])),
                file_created=created)


def read_line_file(path):
    """``fwhm_vnir.dat`` or ``deltawl_vnir.dat``: DataFrame ``wl, L, E`` (file
    order; the lines are not sorted)."""
    a = np.loadtxt(path, comments='#')
    if a.ndim != 2 or a.shape[1] != 4:
        raise ValueError('%s: expected 4 columns, got %s' % (path, a.shape))
    if not np.array_equal(a[:, 0], a[:, 2]):
        raise ValueError('%s: the L and E wl columns differ' % path)
    return pd.DataFrame(dict(wl=a[:, 0], L=a[:, 1], E=a[:, 3]))


def read_fwhm(d):
    return read_line_file(os.path.join(d, 'fwhm_vnir.dat'))


def read_deltawl(d):
    return read_line_file(os.path.join(d, 'deltawl_vnir.dat'))


def read_lamp(d):
    """``lamp_vnir.dat``: ``(wl, pairs)`` with ``pairs`` a list of (L, E)
    column arrays, one per lamp."""
    a = np.loadtxt(os.path.join(d, 'lamp_vnir.dat'), comments='#')
    if a.shape[0] != NPIX or (a.shape[1] - 1) % 2:
        raise ValueError('%s: unexpected lamp_vnir shape %s' % (d, a.shape))
    npair = (a.shape[1] - 1) // 2
    return a[:, 0].copy(), [(a[:, 1 + 2 * k].copy(), a[:, 2 + 2 * k].copy())
                            for k in range(npair)]


def read_epoch(d):
    """All four files of one epoch directory, plus the lamp names."""
    e = dict(dir=d, coefs=read_coefs(d), fwhm=read_fwhm(d), deltawl=read_deltawl(d),
             lamp=read_lamp(d))
    if not np.array_equal(e['fwhm']['wl'].values, e['deltawl']['wl'].values):
        raise ValueError('%s: fwhm and deltawl line lists differ' % d)
    e['lamps'] = lamp_names(e['lamp'][1], e['coefs']['L'])
    return e


# ----------------------------------------------------------------------------
# polynomials and grids
# ----------------------------------------------------------------------------
def poly(c, px):
    """lambda(px) for coefficients c0..c4."""
    return np.polyval(np.asarray(c, float)[::-1], px)


def dpoly(c, px):
    """d lambda / d px."""
    return np.polyval(np.polyder(np.poly1d(np.asarray(c, float)[::-1])), px)


def px_of(c, lam):
    """Fractional pixel of wavelength ``lam`` (Newton)."""
    c = np.asarray(c, float)
    px = (np.asarray(lam, float) - c[0]) / c[1]
    for _ in range(8):
        px = px - (poly(c, px) - lam) / dpoly(c, px)
    return px


def crop_window(c, lo=LAM_LO, hi=LAM_HI):
    """First and last pixel with ``lo <= poly(c, px) < hi`` (the processor's
    crop rule); the polynomial is monotonic over the detector."""
    lam = poly(c, np.arange(NPIX))
    m = np.where((lam >= lo) & (lam < hi))[0]
    return int(m[0]), int(m[-1])


_COEFS = {}


def load_wl_coefs(path=None):
    """``hypernet/data/lab_wl_coefs.csv`` (cached per path)."""
    path = path or os.path.join(DATA_DIR, 'lab_wl_coefs.csv')
    if path not in _COEFS:
        _COEFS[path] = pd.read_csv(path, dtype={'sn': str, 'epoch': str, 'coef_date': str,
                                                'file_created': str, 'lines': str},
                                   keep_default_na=False, na_values={'sigma': ['']},
                                   float_precision='round_trip')   # exact, not 'high'
    return _COEFS[path]


def write_wl_coefs(df, path=None):
    """Write the coefficient table with full float precision."""
    path = path or os.path.join(DATA_DIR, 'lab_wl_coefs.csv')
    out = df[COEF_COLUMNS].copy()
    for k in ('c0', 'c1', 'c2', 'c3', 'c4'):      # shortest exact round trip (repr)
        out[k] = [repr(float(v)) for v in out[k]]
    out.to_csv(path, index=False)
    _COEFS.pop(path, None)
    return path


def coefs(instrument, channel, epoch=None, table=None):
    """c0..c4 of ``instrument`` / ``channel`` at ``epoch`` (default
    :data:`PROCESSOR_EPOCH`) from the committed table."""
    sn = _sn(instrument)
    ch = _channel(channel)
    ep = epoch or PROCESSOR_EPOCH[sn]['epoch']
    t = load_wl_coefs() if table is None else table
    r = t[(t['sn'] == sn) & (t['epoch'] == ep) & (t['channel'] == ch)]
    if len(r) != 1:
        raise KeyError('no lab coefficients for %s %s %s' % (sn, ep, ch))
    return r[['c0', 'c1', 'c2', 'c3', 'c4']].values[0].astype(float)


def _channel(channel):
    ch = {'E': 'E', 'L': 'L', 'Ld': 'L', 'Lu': 'L'}.get(channel)
    if ch is None:
        raise ValueError('channel must be one of E, L, Ld, Lu (got %r)' % (channel,))
    return ch


def grid_for(instrument, channel, epoch=None, table=None, return_px=False):
    """The wavelength grid (nm, air) of ``instrument`` / ``channel``.

    Parameters
    ----------
    instrument : str
        ``'HYPSTAR_122304'`` or ``'122304'``.
    channel : str
        ``'E'``, or ``'L'`` / ``'Ld'`` / ``'Lu'`` (one L grid).
    epoch : str, optional
        Lab epoch; default :data:`PROCESSOR_EPOCH`, i.e. the processor's grid.
        The :data:`CROP_EPOCH` pixel window applies only to the processor
        epoch; any other epoch gets the plain 350 <= lambda < 1100 rule.
    return_px : bool
        Also return the pixel indices.
    """
    sn = _sn(instrument)
    ch = _channel(channel)
    pe = PROCESSOR_EPOCH[sn]['epoch']
    ep = epoch or pe
    c = coefs(sn, ch, ep, table)
    cw = coefs(sn, ch, CROP_EPOCH[sn], table) if (ep == pe and sn in CROP_EPOCH) else c
    s, e = crop_window(cw)
    px = np.arange(s, e + 1)
    lam = poly(c, px)
    return (lam, px) if return_px else lam


# ----------------------------------------------------------------------------
# line fitting in pixel space
# ----------------------------------------------------------------------------
def gauss(x, a, mu, sig, b):
    return b + a * np.exp(-0.5 * ((x - mu) / sig) ** 2)


def supergauss(x, a, mu, sig, b, p):
    """FWHM = FWHM_PER_SIGMA * sig whatever p (as :func:`hypernet.srf.srf_kernel`)."""
    f = FWHM_PER_SIGMA * sig
    s_ = f / (2.0 * (2.0 * np.log(2.0)) ** (1.0 / p))
    return b + a * np.exp(-0.5 * np.abs((x - mu) / s_) ** p)


def asym_gauss(x, a, mu, sig, b, asym):
    """Width sig (1 - asym) on the blue side and sig (1 + asym) on the red;
    asym < 0 is wider on the blue side."""
    s = np.where(x < mu, sig * (1 - asym), sig * (1 + asym))
    return b + a * np.exp(-0.5 * ((x - mu) / s) ** 2)


def halfmax_fwhm(x, y):
    """FWHM in x units from the interpolated half-maximum crossings above the
    window minimum."""
    yy = y - np.min(y)
    i = int(np.argmax(yy))
    h = 0.5 * yy[i]
    try:
        jl = np.where(yy[:i] < h)[0][-1]
        xl = x[jl] + (h - yy[jl]) / (yy[jl + 1] - yy[jl]) * (x[jl + 1] - x[jl])
        jr = i + np.where(yy[i:] < h)[0][0]
        xr = x[jr - 1] + (h - yy[jr - 1]) / (yy[jr] - yy[jr - 1]) * (x[jr] - x[jr - 1])
        return xr - xl
    except IndexError:
        return np.nan


def _fit(x, y, err, model, p0, bounds):
    """curve_fit with absolute errors; returns (p, perr scaled by
    sqrt(chi2_nu), chi2_nu) or (None, None, nan)."""
    try:
        p, cov = curve_fit(model, x, y, p0=p0, sigma=np.full(x.size, err), absolute_sigma=True,
                           bounds=bounds, maxfev=20000)
    except Exception:
        return None, None, np.nan
    chi2 = np.sum(((y - model(x, *p)) / err) ** 2) / (x.size - len(p))
    return p, np.sqrt(np.diag(cov)) * np.sqrt(chi2), chi2


def noise_of(y, wl):
    """Robust std of the lamp column in 450-530 and 600-680 nm (nominal wl)."""
    m = ((wl > 450) & (wl < 530)) | ((wl > 600) & (wl < 680))
    return 1.4826 * np.median(np.abs(y[m] - np.median(y[m])))


def fit_line_px(y, px0, noise, half=7):
    """Fit one lamp line around pixel ``px0`` with the three shapes.

    Returns a dict (pixel units; errors scaled by sqrt(chi2_nu)) or None if
    the Gaussian fit fails or the window leaves the detector.
    """
    lo, hi = int(round(px0)) - half, int(round(px0)) + half + 1
    if lo < 0 or hi > NPIX:
        return None
    x, yy = np.arange(lo, hi, dtype=float), np.asarray(y[lo:hi], float)
    out = dict(snr=(yy.max() - yy.min()) / noise, fwhm_halfmax_px=halfmax_fwhm(x, yy))
    b0, a0 = yy.min(), yy.max() - yy.min()
    p, e, chi = _fit(x, yy, noise, gauss, [a0, px0, 2.5, b0],
                     ([0, lo, 0.5, -np.inf], [np.inf, hi, 10, np.inf]))
    if p is None:
        return None
    out.update(mu_px=p[1], mu_px_err=e[1], fwhm_gauss_px=FWHM_PER_SIGMA * abs(p[2]),
               fwhm_gauss_px_err=FWHM_PER_SIGMA * e[2], chi2_gauss=chi)
    ps, es, chis = _fit(x, yy, noise, supergauss, [p[0], p[1], abs(p[2]), p[3], 2.0],
                        ([0, lo, 0.5, -np.inf, 1.0], [np.inf, hi, 10, np.inf, 8.0]))
    if ps is not None:
        out.update(fwhm_sg_px=FWHM_PER_SIGMA * abs(ps[2]), fwhm_sg_px_err=FWHM_PER_SIGMA * es[2],
                   sg_p=ps[4], sg_p_err=es[4], chi2_sg=chis)
    pa, ea, chia = _fit(x, yy, noise, asym_gauss, [p[0], p[1], abs(p[2]), p[3], 0.0],
                        ([0, lo, 0.5, -np.inf, -0.6], [np.inf, hi, 10, np.inf, 0.6]))
    if pa is not None:
        out.update(asym=pa[4], asym_err=ea[4], chi2_asym=chia)
    return out


def lamp_names(pairs, c_L):
    """Name each L/E pair by its brightest marker line (``HgAr``, ``Kr``, ``Xe``)."""
    names = []
    for L, _ in pairs:
        best = None
        for name, lam in LAMP_MARKERS:
            i = int(round(px_of(c_L, lam)))
            v = L[i - 3:i + 4].max() - np.median(L[i - 40:i + 40])
            if best is None or v > best[0]:
                best = (v, name)
        names.append(best[1])
    return names


def lamp_of_line(lam):
    """The lamp that carries reference line ``lam`` (Hg and Ar: the Hg-Ar pen
    lamp)."""
    sp = SPECIES.get(lam, '?').split()[0]
    return 'HgAr' if sp in ('Hg', 'Ar') else sp


def fit_epoch_lines(sn, epoch, e=None):
    """Rows of ``lab_srf.csv`` for one epoch (both channels, every line of
    ``fwhm_vnir.dat``; a line whose fit fails keeps its lab values and NaN
    fits)."""
    e = e or read_epoch(epoch_dir(sn, epoch))
    c = e['coefs']
    wl, pairs = e['lamp']
    rows = []
    for _, lf in e['fwhm'].iterrows():
        lam = float(lf['wl'])
        dw = e['deltawl'].loc[e['deltawl']['wl'] == lam].iloc[0]
        lamp = lamp_of_line(lam)
        k = e['lamps'].index(lamp) if lamp in e['lamps'] else None
        for col, ch in enumerate(CHANNELS):
            r = dict(instrument='HYPSTAR_' + sn, sn=sn, epoch=epoch, channel=ch, lam_ref=lam,
                     species=SPECIES.get(lam, '?'), lamp=lamp, lab_fwhm=float(lf[ch]),
                     lab_deltawl=float(dw[ch]), blend=lam in BLENDS,
                     outside_grid=not (LAM_LO <= lam < LAM_HI))
            f = None
            if k is not None:
                y = pairs[k][col]
                f = fit_line_px(y, float(px_of(c[ch], lam)), noise_of(y, wl))
            if f is not None:
                mu = f['mu_px']
                disp = float(dpoly(c[ch], mu))
                r.update(mu_px=mu, mu_px_err=f['mu_px_err'], disp_nm_px=disp,
                         lam_fit=float(poly(c[ch], mu)), resid=float(poly(c[ch], mu)) - lam,
                         fwhm_gauss=f['fwhm_gauss_px'] * disp,
                         fwhm_gauss_err=f['fwhm_gauss_px_err'] * disp,
                         fwhm_halfmax=f['fwhm_halfmax_px'] * disp,
                         fwhm_sg=f.get('fwhm_sg_px', np.nan) * disp,
                         fwhm_sg_err=f.get('fwhm_sg_px_err', np.nan) * disp,
                         sg_p=f.get('sg_p', np.nan), sg_p_err=f.get('sg_p_err', np.nan),
                         asym=f.get('asym', np.nan), asym_err=f.get('asym_err', np.nan),
                         chi2_gauss=f['chi2_gauss'], chi2_sg=f.get('chi2_sg', np.nan),
                         chi2_asym=f.get('chi2_asym', np.nan), snr=f['snr'])
            r['low_snr'] = not (r.get('snr', np.nan) >= SNR_MIN)
            rows.append(r)
    return pd.DataFrame(rows).reindex(columns=SRF_COLUMNS)


def coef_rows(sn, epoch, e=None):
    """Rows of ``lab_wl_coefs.csv`` for one epoch."""
    e = e or read_epoch(epoch_dir(sn, epoch))
    c = e['coefs']
    pe = PROCESSOR_EPOCH.get(sn, {})
    rows = []
    for ch in CHANNELS:
        proc = pe.get('epoch') == epoch
        rows.append(dict(instrument='HYPSTAR_' + sn, sn=sn, epoch=epoch, channel=ch,
                         **{'c%d' % i: float(c[ch][i]) for i in range(5)},
                         sigma=c['sigma_' + ch], coef_file=c['coef_file'],
                         coef_date=c['coef_date'].strftime('%Y-%m-%d'),
                         file_created=c['file_created'], lamps='+'.join(e['lamps']),
                         n_lines=len(e['deltawl']),
                         lines=' '.join(repr(float(w)) for w in e['deltawl']['wl']),
                         lines_source='assumed_all_deltawl', processor=proc,
                         processor_status=pe[ch] if proc else 'no'))
    return rows


_SRF = {}


def load_lab_srf(path=None):
    """``hypernet/data/lab_srf.csv`` (cached per path)."""
    path = path or os.path.join(DATA_DIR, 'lab_srf.csv')
    if path not in _SRF:
        _SRF[path] = pd.read_csv(path, dtype={'sn': str, 'epoch': str})
    return _SRF[path]


def write_lab_srf(df, path=None):
    path = path or os.path.join(DATA_DIR, 'lab_srf.csv')
    df[SRF_COLUMNS].to_csv(path, index=False, float_format='%.7g')
    _SRF.pop(path, None)
    return path


def physical_epoch(sn, date, table=None):
    """The latest non-``_tst`` lab epoch of ``sn`` with coefficient date <=
    ``date`` (the lab state of the instrument at that calibration), else the
    earliest; ``(epoch, before)``."""
    t = load_wl_coefs() if table is None else table
    g = t[(t['sn'] == sn) & ~t['epoch'].str.contains('tst')].drop_duplicates('epoch')
    d = pd.to_datetime(g['coef_date'])
    before = g[d <= pd.Timestamp(date)]
    if len(before):
        return before.loc[pd.to_datetime(before['coef_date']).idxmax(), 'epoch'], True
    return g.loc[d.idxmin(), 'epoch'], False


# ----------------------------------------------------------------------------
# the lab SRF shape (lab data prompts, Setup Q6-Q8)
# ----------------------------------------------------------------------------
#: Isolated lamp lines for the shape (as ``wavecal/lab_data_examine.py``).
CLEAN_LINES = (404.6563, 435.8328, 546.0735, 763.5106, 826.4522, 881.9411)
#: Setup Q7: a channel is super-Gaussian where the median chi2_nu of the
#: super-Gaussian fit is below this fraction of the Gaussian one.
SG_CHI2_RATIO = 0.5


def lab_shape(instrument, channel, epoch=None, srf=None):
    """Median lamp-line shape of ``instrument`` / ``channel`` at ``epoch``
    (default :data:`PROCESSOR_EPOCH`), over :data:`CLEAN_LINES` with
    S/N >= :data:`SNR_MIN`.

    Returns:
        dict: ``epoch``, ``n``, ``p`` (median super-Gaussian exponent),
        ``p_min``, ``p_max``, ``asym``, ``chi2_ratio`` (median chi2_sg /
        median chi2_gauss) and ``supergauss`` (chi2_ratio < SG_CHI2_RATIO).
    """
    sn = _sn(instrument)
    ch = _channel(channel)
    ep = epoch or PROCESSOR_EPOCH[sn]['epoch']
    s = load_lab_srf() if srf is None else srf
    g = s[(s['sn'] == sn) & (s['epoch'] == ep) & (s['channel'] == ch) &
          s['lam_ref'].isin(CLEAN_LINES) & ~s['low_snr']]
    if g.empty:
        raise KeyError('no clean lab lines for %s %s %s' % (sn, ep, ch))
    ratio = float(g['chi2_sg'].median() / g['chi2_gauss'].median())
    return dict(epoch=ep, n=int(len(g)), p=float(g['sg_p'].median()),
                p_min=float(g['sg_p'].min()), p_max=float(g['sg_p'].max()),
                asym=float(g['asym'].median()), chi2_ratio=ratio,
                supergauss=ratio < SG_CHI2_RATIO)


def srf_shape_for(instrument, channel, date, srf=None, table=None):
    """The SRF shape for the SRF table (Setup Q6-Q8).

    The class (super-Gaussian or Gaussian) is decided on the processor epoch
    (:func:`lab_shape`, rule :data:`SG_CHI2_RATIO`); L is always Gaussian.  A
    super-Gaussian's exponent p is the lab median at the epoch that describes
    the instrument at ``date`` (:func:`physical_epoch`).

    Returns:
        dict: ``shape`` ('supergauss' or 'gauss'), ``shape_p`` (NaN for
        Gaussian), ``shape_source`` ('lab'), ``lab_epoch``, ``chi2_ratio``.
    """
    sn = _sn(instrument)
    ch = _channel(channel)
    cls = lab_shape(sn, ch, srf=srf)
    ep, _ = physical_epoch(sn, date, table)
    if ch == 'L' or not cls['supergauss']:
        return dict(shape='gauss', shape_p=np.nan, shape_source='lab', lab_epoch=ep,
                    chi2_ratio=cls['chi2_ratio'])
    return dict(shape='supergauss', shape_p=lab_shape(sn, ch, ep, srf=srf)['p'],
                shape_source='lab', lab_epoch=ep, chi2_ratio=cls['chi2_ratio'])

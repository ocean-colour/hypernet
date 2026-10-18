"""High-resolution reference spectra: the TSIS-1 Hybrid Solar Reference
Spectrum (HSRS; Coddington et al. 2021, 2023).

The file is downloaded once, outside the repository, to
``$OS_COLOR/hypernet/ref/`` and verified by its SHA-256.  HSRS v2 spans
202-2730 nm at 0.005 nm resolution (sampled every 0.001 nm), in W m-2 nm-1,
on **vacuum** wavelengths.
"""

import hashlib
import os

import numpy as np

#: TSIS-1 HSRS v2, 0.005 nm resolution (LASP LISIRD).
HSRS_URL = ('https://lasp.colorado.edu/lisird/resources/lasp/hsrs/v2/'
            'hybrid_reference_spectrum_p005nm_resolution_c2022-11-30_with_unc.nc')
#: SHA-256 of that file as downloaded on 2026-09-29.
HSRS_SHA256 = 'dd9f62fb9b39433631013ebf052429f4daddb2bd7e0d970a6d292be6026f3e20'


def ref_root():
    """``$OS_COLOR/hypernet/ref`` (created if needed)."""
    root = os.path.join(os.getenv('OS_COLOR', '.'), 'hypernet', 'ref')
    os.makedirs(root, exist_ok=True)
    return root


def hsrs_path():
    """Local path of the HSRS file (it may not exist yet)."""
    return os.path.join(ref_root(), os.path.basename(HSRS_URL))


def _sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def fetch_hsrs(force=False, verify=True):
    """Download the HSRS file if it is not already present.

    Parameters
    ----------
    force : bool, optional
        Download again even if the file exists.
    verify : bool, optional
        Check the SHA-256 against :data:`HSRS_SHA256`.

    Returns
    -------
    str
        Path to the local file.

    Raises
    ------
    ValueError
        If the checksum does not match.
    """
    path = hsrs_path()
    if force or not os.path.exists(path):
        import urllib.request
        tmp = path + '.part'
        urllib.request.urlretrieve(HSRS_URL, tmp)
        os.replace(tmp, path)
    if verify and _sha256(path) != HSRS_SHA256:
        raise ValueError('HSRS checksum mismatch for %s' % path)
    return path


def hsrs_available():
    """True if the HSRS file is on disk."""
    return os.path.exists(hsrs_path())


def load_hsrs(wmin=350.0, wmax=1100.0, step=5, frame='vac', path=None):
    """Load the HSRS between ``wmin`` and ``wmax``.

    Parameters
    ----------
    wmin, wmax : float, optional
        Wavelength range, nm, in ``frame``.
    step : int, optional
        Block-average this many native samples (0.001 nm each).  The default
        of 5 gives 0.005 nm samples, the file's stated resolution.
    frame : {'vac', 'air'}, optional
        Wavelength scale of the returned grid.  ``'air'`` converts with
        :func:`hypernet.srf.vac_to_air`.
    path : str, optional
        File to read; default :func:`hsrs_path` (downloaded if absent).

    Returns
    -------
    tuple
        ``(wave, ssi)``: nm (in ``frame``) and W m-2 nm-1.
    """
    import xarray as xr
    from hypernet import srf

    path = path or fetch_hsrs(verify=False)
    ds = xr.open_dataset(path)
    try:
        w = ds['Vacuum Wavelength'].values
        # select on the vacuum scale with a margin; frame conversion after
        lo, hi = wmin - 1.0, wmax + 1.0
        i0, i1 = np.searchsorted(w, [lo, hi])
        w = w[i0:i1].astype(float)
        f = ds['SSI'].values[i0:i1].astype(float)
    finally:
        ds.close()
    if step > 1:
        n = (w.size // step) * step
        w = w[:n].reshape(-1, step).mean(axis=1)
        f = f[:n].reshape(-1, step).mean(axis=1)
    if frame == 'air':
        w = srf.vac_to_air(w)
    elif frame != 'vac':
        raise ValueError("frame must be 'vac' or 'air'")
    m = (w >= wmin) & (w <= wmax)
    return w[m], f[m]

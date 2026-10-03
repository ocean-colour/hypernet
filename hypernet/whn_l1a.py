"""Readers for WATERHYPERNET L1A, L1C and L2A sequence files.

Release 2 as distributed holds only L2A/L2B.  The L1A files, with the
irradiance (E) and radiance (L) on their own wavelength grids, come from RBINS
on request and are mirrored to ``$OS_COLOR/WATERHYPERNET/Wavelengths``.

Conventions established on the VEIT sample (SEQ20260604T084543;
``wavecal/check_vza_convention_veit.py``, ``hypernet/wiggles/phase0a_l1c_consistency.py``):

- ``viewing_zenith_angle`` is measured from **nadir**: the irradiance sensor
  looks up at 180 deg.  In L1A_RAD, scans with vza < 90 are the water view
  (Lu) and scans with vza >= 90 are the sky (Ld).
- A sequence has two sky series, one before and one after the water series.
  L1C ``downwelling_radiance`` is the two sky-series means interpolated
  linearly in time to the water-view time (:func:`interp_in_time`).  L1C
  ``upwelling_radiance`` is the L1A water scans unchanged, and L1C
  ``irradiance`` is ``np.interp`` of the time-interpolated E onto the L grid.
  E, too, has two series; in the sample they bracket the water view
  symmetrically, so the time interpolation equals the plain mean.
- L1C ``u_rel_random_*`` and ``std_downwelling_radiance`` are zero in the
  sample: placeholders, not uncertainties.
- ``instrument_calibration_file_rad`` points to the IRR file (Kevin's
  attribute bug); pair calibration files by instrument and date instead.

Every loader returns a plain dict (as :func:`hypernet.whn_explore.load_spectrum`
does).  Spectral arrays keep the file layout, (wavelength, scan) or
(wavelength, series).
"""

import glob
import os
import re

import numpy as np
import pandas as pd

#: The four products of a sequence.
PRODUCTS = ('L1A_IRR', 'L1A_RAD', 'L1C_ALL', 'L2A_REF')

#: vza below this is the water view (Lu), at or above it the sky (Ld).
VZA_SPLIT = 90.0

_NAME_RE = re.compile(
    r'HYPERNETS_(?P<network>[A-Z])_(?P<site>[A-Z0-9]+)_(?P<product>L\d[A-Z]_[A-Z]+)_'
    r'(?P<seq_time>\d{8}T\d{4})_(?P<proc_time>\d{8}T\d{4})'
    r'(?:_(?P<azimuth>\d{3}))?_(?P<version>v[\d.]+)\.nc$')

_META_ATTRS = ('sequence_id', 'system_id', 'instrument_id', 'site_id',
               'processor_version', 'product_level',
               'instrument_calibration_file_irr', 'instrument_calibration_date_irr',
               'instrument_calibration_file_rad', 'instrument_calibration_date_rad')


def wavelengths_root(path=None):
    """Return the directory holding the L1A/L1C/L2A sequence files.

    Resolution order: explicit ``path``, then ``$OS_COLOR/WATERHYPERNET/Wavelengths``.

    Parameters
    ----------
    path : str or None, optional
        Explicit directory.

    Returns
    -------
    str
        Absolute path to the existing directory.

    Raises
    ------
    FileNotFoundError
        If neither candidate exists on disk.
    """
    candidates = []
    if path is not None:
        candidates.append(path)
    os_color = os.getenv('OS_COLOR')
    if os_color is not None:
        candidates.append(os.path.join(os_color, 'WATERHYPERNET', 'Wavelengths'))
    for c in candidates:
        if os.path.isdir(c):
            return os.path.abspath(c)
    raise FileNotFoundError(
        'No Wavelengths directory found; tried: %s' % (candidates or
                                                        ['$OS_COLOR unset']))


def parse_name(path):
    """Parse a HYPERNETS product filename.

    Parameters
    ----------
    path : str
        File path or basename, e.g.
        ``HYPERNETS_W_VEIT_L1C_ALL_20260604T0845_20260828T1547_090_v2.1.nc``.

    Returns
    -------
    dict or None
        ``network, site, product, seq_time, proc_time, azimuth, version``
        (``azimuth`` is None for L1A), or None if the name does not match.
    """
    m = _NAME_RE.search(os.path.basename(path))
    return m.groupdict() if m else None


def find_products(root=None, recursive=True):
    """Index every HYPERNETS product file under ``root``.

    Parameters
    ----------
    root : str or None, optional
        Directory to search; default :func:`wavelengths_root`.
    recursive : bool, optional
        Search subdirectories too (deliveries may arrive in subfolders).

    Returns
    -------
    pandas.DataFrame
        One row per file: ``path`` plus the fields of :func:`parse_name`.
    """
    root = wavelengths_root(root)
    pattern = os.path.join(root, '**', 'HYPERNETS_*.nc') if recursive \
        else os.path.join(root, 'HYPERNETS_*.nc')
    rows = []
    for p in sorted(glob.glob(pattern, recursive=recursive)):
        info = parse_name(p)
        if info is not None:
            rows.append(dict(path=p, **info))
    cols = ['path', 'network', 'site', 'product', 'seq_time', 'proc_time',
            'azimuth', 'version']
    return pd.DataFrame(rows, columns=cols)


def sequence_files(site=None, seq_time=None, root=None, index=None):
    """The product files of one sequence, keyed by product.

    Parameters
    ----------
    site, seq_time : str or None, optional
        Site code (``'VEIT'``) and sequence time (``'20260604T0845'``).
        None matches anything, which is fine when the directory holds a
        single sequence.
    root : str or None, optional
        Passed to :func:`find_products` if ``index`` is not given.
    index : pandas.DataFrame, optional
        Output of :func:`find_products`.

    Returns
    -------
    dict
        product -> path for the products found (a subset of :data:`PRODUCTS`).

    Raises
    ------
    ValueError
        If a product matches more than one file.
    """
    df = find_products(root) if index is None else index
    if site is not None:
        df = df[df['site'] == site]
    if seq_time is not None:
        df = df[df['seq_time'] == seq_time]
    out = {}
    for prod, g in df.groupby('product'):
        if len(g) > 1:
            raise ValueError('%d %s files match site=%s seq_time=%s'
                             % (len(g), prod, site, seq_time))
        out[prod] = g['path'].iloc[0]
    return out


def _open(path):
    import xarray as xr
    return xr.open_dataset(path)


def _meta(ds, path):
    meta = {a: ds.attrs.get(a) for a in _META_ATTRS}
    meta['path'] = path
    return meta


def _geometry(ds, sel=None):
    """Per-scan geometry and bookkeeping, optionally for a subset of scans."""
    names = {'sza': 'solar_zenith_angle', 'saa': 'solar_azimuth_angle',
             'vza': 'viewing_zenith_angle', 'vaa': 'viewing_azimuth_angle',
             'paa': 'pointing_azimuth_angle', 'time': 'acquisition_time',
             'series_id': 'series_id', 'quality_flag': 'quality_flag'}
    out = {}
    for k, v in names.items():
        if v in ds:
            a = ds[v].values
            if k == 'time':
                a = a.astype(np.int64)
            out[k] = a if sel is None else a[sel]
    return out


def series_means(scans, series_id, time):
    """Mean spectrum and acquisition time of each series.

    Parameters
    ----------
    scans : numpy.ndarray
        (wavelength, scan).
    series_id, time : numpy.ndarray
        Per-scan series id and acquisition time (s).

    Returns
    -------
    tuple
        ``(means, times, ids)``: (wavelength, n_series), (n_series,),
        (n_series,), ordered by time.
    """
    ids = np.unique(series_id)
    times = np.array([np.mean(time[series_id == i]) for i in ids], dtype=float)
    order = np.argsort(times)
    ids, times = ids[order], times[order]
    means = np.stack([np.nanmean(scans[:, series_id == i], axis=1) for i in ids],
                     axis=1)
    return means, times, ids


def interp_in_time(scans, series_id, time, t):
    """Series means interpolated linearly in time to ``t``.

    This is how the processor builds L1C ``downwelling_radiance`` from the two
    sky series.  ``t`` outside the series times takes the nearest series.  A
    single series returns its mean.

    Parameters
    ----------
    scans : numpy.ndarray
        (wavelength, scan).
    series_id, time : numpy.ndarray
        Per-scan series id and acquisition time (s).
    t : float
        Target time (s), e.g. the water-view time.

    Returns
    -------
    numpy.ndarray
        (wavelength,).
    """
    means, times, _ = series_means(scans, series_id, time)
    if times.size == 1:
        return means[:, 0]
    j = int(np.clip(np.searchsorted(times, t), 1, times.size - 1))
    w = np.clip((t - times[j - 1]) / (times[j] - times[j - 1]), 0.0, 1.0)
    return (1.0 - w) * means[:, j - 1] + w * means[:, j]


def load_l1a_irr(path):
    """Read an L1A_IRR file: downwelling irradiance on the E grid.

    Parameters
    ----------
    path : str
        Path to the ``.nc`` file.

    Returns
    -------
    dict
        ``wave`` (n_E,), ``scans`` (n_E, n_scan), ``mean`` (n_E,),
        ``bandwidth``, ``n_scans``, the per-scan geometry (``sza, saa, vza,
        vaa, paa, time, series_id, quality_flag``) and ``meta`` (file
        attributes and path).
    """
    ds = _open(path)
    try:
        E = ds['irradiance'].values.astype(float)
        out = dict(wave=ds['wavelength'].values.astype(float), scans=E,
                   mean=np.nanmean(E, axis=1),
                   bandwidth=ds['bandwidth'].values.astype(float),
                   n_scans=E.shape[1], meta=_meta(ds, path))
        out.update(_geometry(ds))
    finally:
        ds.close()
    return out


def load_l1a_rad(path):
    """Read an L1A_RAD file, split into the water view (Lu) and the sky (Ld).

    Parameters
    ----------
    path : str
        Path to the ``.nc`` file.

    Returns
    -------
    dict
        ``wave`` (n_L,), ``bandwidth``, ``meta``, and for each of ``'Lu'``
        (vza < :data:`VZA_SPLIT`) and ``'Ld'`` (vza >= it) a sub-dict with
        ``scans`` (n_L, n), ``mean`` (n_L,), ``n_scans`` and the per-scan
        geometry.  ``Ld['at_lu_time']`` is the sky interpolated in time to
        the mean water-view time, as in L1C.
    """
    ds = _open(path)
    try:
        R = ds['radiance'].values.astype(float)
        vza = ds['viewing_zenith_angle'].values
        out = dict(wave=ds['wavelength'].values.astype(float),
                   bandwidth=ds['bandwidth'].values.astype(float),
                   n_scans=R.shape[1], meta=_meta(ds, path))
        for key, sel in (('Lu', vza < VZA_SPLIT), ('Ld', vza >= VZA_SPLIT)):
            sub = dict(scans=R[:, sel], mean=np.nanmean(R[:, sel], axis=1),
                       n_scans=int(sel.sum()))
            sub.update(_geometry(ds, sel))
            out[key] = sub
    finally:
        ds.close()
    t_lu = float(np.mean(out['Lu']['time']))
    out['Lu']['t_mean'] = t_lu
    Ld = out['Ld']
    Ld['at_lu_time'] = interp_in_time(Ld['scans'], Ld['series_id'], Ld['time'], t_lu)
    return out


_L1C_VARS = ('irradiance', 'downwelling_radiance', 'upwelling_radiance',
             'water_leaving_radiance', 'reflectance', 'reflectance_nosc')
_L2A_EXTRA = ('std_reflectance', 'std_reflectance_nosc', 'std_water_leaving_radiance')
_SCALARS = ('rhof', 'rhof_wind', 'epsilon', 'n_valid_scans', 'n_total_scans')


def _load_l1c_l2a(path, extra=()):
    ds = _open(path)
    try:
        dim = 'scan' if 'scan' in ds.dims else 'series'
        out = dict(wave=ds['wavelength'].values.astype(float), dim=dim,
                   n=ds.sizes[dim], meta=_meta(ds, path))
        for v in _L1C_VARS + tuple(extra):
            if v in ds:
                out[v] = ds[v].values.astype(float)
        for v in _SCALARS:
            if v in ds:
                out[v] = ds[v].values
        out.update(_geometry(ds))
    finally:
        ds.close()
    return out


def load_l1c(path):
    """Read an L1C_ALL file (all quantities on the L grid, per water scan).

    Parameters
    ----------
    path : str
        Path to the ``.nc`` file.

    Returns
    -------
    dict
        ``wave``, and (wavelength, scan) arrays ``irradiance``,
        ``downwelling_radiance``, ``upwelling_radiance``,
        ``water_leaving_radiance``, ``reflectance``, ``reflectance_nosc``;
        per-scan ``rhof``, ``rhof_wind``, ``epsilon`` and geometry;
        ``n`` (scans), ``dim`` and ``meta``.
    """
    return _load_l1c_l2a(path)


def load_l2a(path):
    """Read an L2A_REF file (the sequence average, per series).

    Parameters
    ----------
    path : str
        Path to the ``.nc`` file.

    Returns
    -------
    dict
        As :func:`load_l1c`, with (wavelength, series) arrays, plus
        ``std_reflectance``, ``std_reflectance_nosc``,
        ``std_water_leaving_radiance``, ``n_valid_scans``, ``n_total_scans``.
    """
    return _load_l1c_l2a(path, extra=_L2A_EXTRA)


def load_l2b(path):
    """Read a Release 2 L2B_REF file.

    In Release 2, L2B files carry the same variables as L2A (``product_level``
    is even ``W_L2A``): the sequence-mean ``downwelling_radiance`` (Ld),
    ``upwelling_radiance``, ``irradiance`` (E resampled onto the L grid),
    ``water_leaving_radiance`` and the reflectances, on the native L grid, with
    ``n_valid_scans`` but no per-scan spectra.  See :func:`load_l2a`.
    """
    return _load_l1c_l2a(path, extra=_L2A_EXTRA)


def release2_path(site, sequence_time, filename, root=None):
    """Path of a Release 2 file: ``RELEASE_2/<site>/<YYYY>/<MM>/<DD>/<filename>``.

    Parameters
    ----------
    site : str
        Site with suffix, e.g. ``'VEIT_H'``.
    sequence_time : str
        ``'YYYYMMDDTHHMM'``.
    filename : str
        Basename, e.g. from the ``file`` column of the data request.
    root : str, optional
        The ``RELEASE_2`` directory; default
        :func:`hypernet.whn_explore.whn_root`.
    """
    if root is None:
        from hypernet.whn_explore import whn_root
        root = whn_root()
    t = sequence_time
    return os.path.join(root, site, t[:4], t[4:6], t[6:8], filename)


def check_against_l1c(irr, rad, l1c, wmin=400.0, wmax=900.0):
    """Compare the L1A readers' products with the processor's L1C.

    Parameters
    ----------
    irr, rad, l1c : dict
        Outputs of :func:`load_l1a_irr`, :func:`load_l1a_rad`, :func:`load_l1c`.
    wmin, wmax : float, optional
        Wavelength range (nm) of the comparison.

    Returns
    -------
    dict
        Maximum relative differences: ``E`` (``np.interp`` onto the L grid
        of the L1A E interpolated in time to the water view, vs L1C
        ``irradiance``, all scans), ``Lu`` (L1A water
        scans vs L1C ``upwelling_radiance``) and ``Ld`` (time-interpolated sky
        vs L1C ``downwelling_radiance``).
    """
    w = rad['wave']
    band = (w >= wmin) & (w <= wmax)

    def maxrel(a, b):
        a, b = np.broadcast_arrays(a, b)
        r = np.abs(a[band] - b[band]) / np.abs(b[band])
        return float(np.nanmax(r))

    E_t = interp_in_time(irr['scans'], irr['series_id'], irr['time'],
                         rad['Lu']['t_mean'])
    E_L = np.interp(w, irr['wave'], E_t)
    Lu = rad['Lu']['scans'] if rad['Lu']['n_scans'] == l1c['n'] else rad['Lu']['mean'][:, None]
    return dict(E=maxrel(E_L[:, None], l1c['irradiance']),
                Lu=maxrel(Lu, l1c['upwelling_radiance']),
                Ld=maxrel(rad['Ld']['at_lu_time'][:, None], l1c['downwelling_radiance']))


def sequence_table(products):
    """One row per delivered sequence (site, seq_time), with a path column per
    product and the L1C/L2A azimuth.

    Parameters
    ----------
    products : pandas.DataFrame
        Output of :func:`find_products`.

    Returns
    -------
    pandas.DataFrame
        ``site, seq_time, azimuth, n_products, complete`` and one column per
        product in :data:`PRODUCTS` (path or None).  A product delivered more
        than once for a sequence keeps the latest ``proc_time``, and the count
        of duplicates goes in ``n_duplicates``.
    """
    rows = []
    for (site, t), g in products.groupby(['site', 'seq_time']):
        row = dict(site=site, seq_time=t, n_duplicates=0)
        az = g['azimuth'].dropna().unique()
        row['azimuth'] = az[0] if len(az) == 1 else (','.join(sorted(az)) if len(az) else None)
        for prod in PRODUCTS:
            h = g[g['product'] == prod].sort_values('proc_time')
            row[prod] = h['path'].iloc[-1] if len(h) else None
            row['n_duplicates'] += max(len(h) - 1, 0)
        row['n_products'] = sum(row[p] is not None for p in PRODUCTS)
        row['complete'] = row['n_products'] == len(PRODUCTS)
        rows.append(row)
    cols = ['site', 'seq_time', 'azimuth', 'n_products', 'complete', 'n_duplicates'] + \
        list(PRODUCTS)
    return pd.DataFrame(rows, columns=cols)


def match_request(sequences, request):
    """Match delivered sequences to the rows of the data request.

    The request (``docs/wiggles_data_request.csv``) names a sequence by
    ``site`` (``'VEIT_H'``), ``sequence_time`` (``'20260604T0845'``) and
    ``azimuth`` (``'090'``); product filenames drop the ``_H`` suffix, and
    only L1C/L2A carry the azimuth.

    Parameters
    ----------
    sequences : pandas.DataFrame
        Output of :func:`sequence_table`.
    request : pandas.DataFrame
        The request table, read with ``dtype=str``.

    Returns
    -------
    pandas.DataFrame
        One row per requested sequence plus one per unrequested delivered
        sequence.  ``status`` is ``'delivered'``, ``'missing'`` or ``'extra'``;
        ``azimuth_ok`` compares the requested and delivered azimuths (None
        when either is unknown); the request columns are carried over.
    """
    req = request.copy()
    req['site_code'] = req['site'].str.replace(r'_[HP]$', '', regex=True)
    seq = sequences.rename(columns={'site': 'site_code', 'seq_time': 'sequence_time',
                                    'azimuth': 'azimuth_delivered'})
    m = req.merge(seq, on=['site_code', 'sequence_time'], how='outer', indicator=True)
    m['status'] = m['_merge'].map({'both': 'delivered', 'left_only': 'missing',
                                   'right_only': 'extra'}).astype(str)
    m = m.drop(columns='_merge')

    def _az(r):
        a, b = r.get('azimuth'), r.get('azimuth_delivered')
        if not isinstance(a, str) or not isinstance(b, str):
            return None
        return a in b.split(',')
    m['azimuth_ok'] = m.apply(_az, axis=1)
    return m

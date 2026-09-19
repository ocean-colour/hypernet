"""Exploration of the WATERHYPERNET Release 2 archive (data layer).

This module is the data half of the WATERHYPERNET exploration: it indexes the
archive, reads single spectra with the agreed product/units conventions, and
draws an optically-diverse ~100-spectrum sample per site. Figures and the
summary table are produced by the companion script ``docs/whn_figures.py``.

Conventions agreed in ``claude_prompts/explore_prompts.md`` (Q&A rounds 1-2),
all implemented here:

- **Product.** ``reflectance`` (NIR Similarity-Spectrum corrected) at every site
  *except* the three the release notes call out as poor -- ``LPAR_H``,
  ``MAFR_H``, ``O1BE_P`` -- which use ``reflectance_nosc`` (see
  :data:`NOSC_SITES`). The correction is **NIR**-based for *both* systems; what
  is SWIR-based for HYPSTAR in Release 2 is the separate QC test, which is a
  different thing and easy to conflate.
- **Units.** Both products are water-leaving reflectance ``rho_w`` (verified
  numerically: ``reflectance_nosc == pi * Lw / Ed`` to machine precision), so
  everything here returns **Rrs = rho_w / pi** [1/sr] and **sigma = std / pi**.
- **Uncertainty.** ``sigma`` is the archive's own per-band standard deviation,
  used as-is with no floor. It is scan-to-scan variability only and therefore
  understates the true uncertainty (no calibration or glint-correction term).
- **No trimming, no resampling of the data itself.** Spectra are returned on
  their native grids. The common grid in :data:`ANALYSIS_WAVE` exists purely so
  that spectra from different instruments can be *displayed* and *clustered*
  together; nothing downstream consumes it as the record grid.
- **Negative reflectance is kept.** The release retains it deliberately, and so
  do we.

Stages, run from the repository root::

    python -m hypernet.whn_explore 1   # index the archive from filenames -> index.parquet
    python -m hypernet.whn_explore 2   # read a pool, cluster, sample ~100/site -> sample.*

Intermediates (parquet/npz) are written under ``$OS_COLOR/hypernet/whn_explore``
rather than into the repo; only the figures and the summary table land in the
repo, under ``docs/``.
"""

from __future__ import annotations

import os
import re
import glob

import numpy as np
import pandas as pd

# --- conventions -------------------------------------------------------------

#: Sites where the release notes say the SimSpec-corrected product is poor and
#: ``reflectance_nosc`` is "definitely recommended" (turbid waters).
NOSC_SITES = ('LPAR_H', 'MAFR_H', 'O1BE_P')

#: Common grid for *display and clustering only* -- never a record grid.
ANALYSIS_WAVE = np.arange(350.0, 900.0 + 0.1, 2.5)

#: Sub-range used for the spectral-shape clustering (both systems cover it well,
#: and it avoids the <400 / >900 nm ranges the release notes call unreliable).
CLUSTER_RANGE = (400.0, 800.0)

#: Bands for the time-series figure (blue / green / red).
TS_BANDS = (490.0, 560.0, 665.0)

N_CLUSTERS = 8          # optical water types for the pooled clustering
N_POOL = 400            # candidate spectra read per site before sampling
N_SAMPLE = 100          # spectra kept per site (the "~100 representative")
SEED = 1234             # matches the sweep config's seed


def whn_root(path=None):
    """Return the RELEASE_2 directory of the WATERHYPERNET archive.

    Resolution order: explicit ``path``, then ``$OS_COLOR/WATERHYPERNET/RELEASE_2``.

    Parameters
    ----------
    path : str or None, optional
        Explicit path to the ``RELEASE_2`` directory.

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
        candidates.append(os.path.join(os_color, 'WATERHYPERNET', 'RELEASE_2'))
    for cand in candidates:
        if cand is not None and os.path.isdir(cand):
            return os.path.abspath(cand)
    raise FileNotFoundError(
        'WATERHYPERNET RELEASE_2 not found. Tried: '
        + (', '.join(str(c) for c in candidates) or '<none>')
        + '. Pass path= or set $OS_COLOR.')


def out_root():
    """Directory for generated intermediates (index/sample tables, spectra).

    Kept outside the repo -- ``$OS_COLOR/hypernet/whn_explore`` -- so parquet and
    npz artifacts are not committed, matching the project's artifact split.
    """
    root = os.path.join(os.getenv('OS_COLOR', '.'), 'hypernet', 'whn_explore')
    os.makedirs(root, exist_ok=True)
    return root


def system_of(site):
    """``'HYPSTAR'`` for a ``*_H`` site code, ``'PANTHYR'`` for ``*_P``."""
    return 'HYPSTAR' if site.endswith('_H') else 'PANTHYR'


def product_for(site):
    """Return ``(reflectance_var, std_var)`` for ``site``.

    Implements the agreed per-site product choice: the SimSpec-corrected
    ``reflectance`` everywhere except the three turbid sites in
    :data:`NOSC_SITES`, which use ``reflectance_nosc``. The two systems spell
    the standard-deviation variable differently (``std_reflectance`` vs
    ``reflectance_std``), which is why this returns both names.
    """
    nosc = site in NOSC_SITES
    if system_of(site) == 'HYPSTAR':
        return ('reflectance_nosc', 'std_reflectance_nosc') if nosc else \
               ('reflectance', 'std_reflectance')
    return ('reflectance_nosc', 'reflectance_nosc_std') if nosc else \
           ('reflectance', 'reflectance_std')


# --- stage 1: index ----------------------------------------------------------

#: Acquisition timestamp in both filename grammars (HYPSTAR minute precision,
#: PANTHYR second precision).
_ACQ_RE = re.compile(r'_REF_(\d{8}T\d{4,6})_')


def _azimuth_token(basename, system):
    """Relative-azimuth token from a filename (field 7 HYPSTAR, 6 PANTHYR)."""
    parts = basename[:-3].split('_')
    idx = 7 if system == 'HYPSTAR' else 6
    return parts[idx] if len(parts) > idx else ''


def build_index(path=None, save=True):
    """Index every spectrum file in the archive from its **filename alone**.

    Opening 56k NetCDF files is unnecessary to enumerate the archive, so this
    parses the two filename grammars instead: site, system, acquisition time and
    relative-azimuth token. The result is the enumeration every later stage
    samples from.

    Parameters
    ----------
    path : str or None, optional
        Explicit ``RELEASE_2`` path (see :func:`whn_root`).
    save : bool, optional
        Write ``index.parquet`` under :func:`out_root` (default True).

    Returns
    -------
    pandas.DataFrame
        Columns ``site, system, datetime, azimuth, path``, sorted by site then
        acquisition time.
    """
    root = whn_root(path)
    sites = sorted(d for d in os.listdir(root)
                   if os.path.isdir(os.path.join(root, d)) and not d.startswith('0_'))

    rows = []
    for site in sites:
        system = system_of(site)
        for f in glob.glob(os.path.join(root, site, '*', '*', '*', '*.nc')):
            base = os.path.basename(f)
            m = _ACQ_RE.search(base)
            if m is None:                       # unparsable name -- record it
                continue
            stamp = m.group(1)
            fmt = '%Y%m%dT%H%M%S' if len(stamp) == 15 else '%Y%m%dT%H%M'
            rows.append((site, system, pd.to_datetime(stamp, format=fmt),
                         _azimuth_token(base, system), f))

    df = pd.DataFrame(rows, columns=['site', 'system', 'datetime',
                                     'azimuth', 'path'])
    df = df.sort_values(['site', 'datetime']).reset_index(drop=True)
    if save:
        df.to_parquet(os.path.join(out_root(), 'index.parquet'))
    return df


# --- reading one spectrum ----------------------------------------------------

def load_spectrum(path, site, analysis_grid=True):
    """Read one WATERHYPERNET file as Rrs on its native grid.

    Applies the agreed conventions: the site's product (:func:`product_for`),
    and the ``rho_w -> Rrs`` conversion (divide by pi) for both the reflectance
    and its standard deviation.

    Parameters
    ----------
    path : str
        Path to the ``.nc`` file.
    site : str
        Site code (selects the product and the std variable name).
    analysis_grid : bool, optional
        Also return ``Rrs_grid``, the spectrum linearly interpolated onto
        :data:`ANALYSIS_WAVE` for display/clustering (NaN outside the native
        range). Default True.

    Returns
    -------
    dict
        ``wave``, ``Rrs``, ``sigma`` (native grid), ``Rrs_grid`` (or None), plus
        metadata: ``sza``, ``vza``, ``saa``, ``vaa``, ``quality_flag``,
        ``serial``, ``lat``, ``lon``, ``n_bands``.
    """
    import netCDF4

    refl_var, std_var = product_for(site)
    ds = netCDF4.Dataset(path)
    try:
        wave = np.asarray(ds.variables['wavelength'][:], dtype=float).ravel()
        Rrs = _filled(ds.variables[refl_var]) / np.pi
        sigma = _filled(ds.variables[std_var]) / np.pi

        meta = {
            'sza': _scalar(ds, 'solar_zenith_angle'),
            'vza': _scalar(ds, 'viewing_zenith_angle'),
            'saa': _scalar(ds, 'solar_azimuth_angle'),
            'vaa': _scalar(ds, 'viewing_azimuth_angle'),
            'quality_flag': _quality_flag(ds),
            'serial': _serial(ds),
            'lat': _site_lat(ds),
            'lon': _site_lon(ds),
            'n_bands': wave.size,
        }
    finally:
        ds.close()

    order = np.argsort(wave)
    wave, Rrs, sigma = wave[order], Rrs[order], sigma[order]

    grid = None
    if analysis_grid:
        good = np.isfinite(Rrs)
        grid = (np.interp(ANALYSIS_WAVE, wave[good], Rrs[good],
                          left=np.nan, right=np.nan)
                if good.sum() > 2 else np.full(ANALYSIS_WAVE.size, np.nan))

    out = {'wave': wave, 'Rrs': Rrs, 'sigma': sigma, 'Rrs_grid': grid}
    out.update(meta)
    return out


def _filled(var):
    """Variable values as a float array with fill/masked entries as NaN.

    The cast to float happens *before* ``filled`` because several variables are
    integer-typed (PANTHYR's ``quality_flag`` is uint8 with ``_FillValue=0``,
    HYPSTAR's angles are scaled uint16) and NaN cannot be stored in those.
    """
    arr = var[:]
    if np.ma.isMaskedArray(arr):
        arr = arr.astype(float).filled(np.nan)
    return np.asarray(arr, dtype=float).ravel()


def _scalar(ds, name):
    """First element of a per-series variable as a float (NaN if absent/masked)."""
    if name not in ds.variables:
        return np.nan
    arr = _filled(ds.variables[name])
    return float(arr[0]) if arr.size else np.nan


def _quality_flag(ds):
    """``quality_flag`` read with masking disabled, so that 0 means *passed*.

    PANTHYR declares ``_FillValue = 0`` on this variable, but the variable is a
    bitmask (``flag_masks = 1, 2, 4 ...``) whose zero value means *no flags set*,
    i.e. the measurement **passed**. Under ordinary masked reads every passing
    PANTHYR measurement therefore comes back as missing, and filtering on
    ``quality_flag == 0`` silently discards all 12,080 PANTHYR spectra. Reading
    raw recovers the real value.

    The trade-off is deliberate: because the release declares the meaningful
    value as the fill value, a genuinely absent flag cannot be distinguished
    from a pass in these files, and a pass is what 0 means. HYPSTAR declares no
    ``_FillValue`` here at all, so raw and masked reads agree for it.

    Only ``quality_flag`` is read this way. PANTHYR's angles carry the same
    ``_FillValue = 0`` but are *genuinely* fill -- ``solar_azimuth_angle`` is
    empty in every file sampled -- so they keep the mask and stay NaN.
    """
    if 'quality_flag' not in ds.variables:
        return np.nan
    var = ds.variables['quality_flag']
    var.set_auto_mask(False)
    arr = np.asarray(var[:], dtype=float).ravel()
    return float(arr[0]) if arr.size else np.nan


def _serial(ds):
    """Instrument identity: HYPSTAR ``system_id``, PANTHYR radiance-sensor S/N.

    The HYPSTAR wavelength grid changes with the instrument, so this is what
    distinguishes the two grids found at GAIT, LPAR, MAFR and VEIT_H.
    """
    for attr in ('system_id', 'l_sensor_sn'):
        if hasattr(ds, attr):
            return str(getattr(ds, attr))
    return ''


def _site_lat(ds):
    """Site latitude (HYPSTAR fixed attribute, PANTHYR per-file GPS average)."""
    for attr in ('site_latitude', 'latitude_average'):
        if hasattr(ds, attr):
            return float(getattr(ds, attr))
    return np.nan


def _site_lon(ds):
    """Site longitude (see :func:`_site_lat`)."""
    for attr in ('site_longitude', 'longitude_average'):
        if hasattr(ds, attr):
            return float(getattr(ds, attr))
    return np.nan


# --- stage 2: pool, cluster, sample -----------------------------------------

def read_pool(index, site, n_pool=N_POOL, verbose=True):
    """Read an evenly-spaced candidate pool of spectra for one site.

    The pool is spread uniformly through the site's date-sorted record list so
    that the clustering downstream sees the site's whole history rather than one
    season.

    Returns
    -------
    tuple
        ``(rows, grids)`` -- a DataFrame of per-spectrum metadata and an
        ``(n, len(ANALYSIS_WAVE))`` array of interpolated spectra.
    """
    sub = index[index.site == site].reset_index(drop=True)
    step = max(1, len(sub) // n_pool)
    sel = sub.iloc[::step].head(n_pool)

    rows, grids = [], []
    for i, r in enumerate(sel.itertuples()):
        try:
            spec = load_spectrum(r.path, site)
        except Exception as exc:                # a corrupt file must not stop us
            print(f'  ! {os.path.basename(r.path)}: {type(exc).__name__}: {exc}')
            continue
        rows.append({
            'site': site, 'system': system_of(site), 'datetime': r.datetime,
            'azimuth': r.azimuth, 'path': r.path, 'serial': spec['serial'],
            'sza': spec['sza'], 'vza': spec['vza'], 'lat': spec['lat'],
            'lon': spec['lon'], 'n_bands': spec['n_bands'],
            'wave_min': spec['wave'][0], 'wave_max': spec['wave'][-1],
            'rel_sigma': _rel_sigma(spec),
            'any_negative': _any_negative(spec),
        })
        grids.append(spec['Rrs_grid'])
        if verbose and (i + 1) % 100 == 0:
            print(f'  {site}: {i + 1}/{len(sel)}')

    # Shape-stable even when every read failed, so the caller can stack blindly.
    arr = (np.asarray(grids) if grids
           else np.zeros((0, ANALYSIS_WAVE.size), dtype=float))
    return pd.DataFrame(rows), arr


def _rel_sigma(spec):
    """Median ``sigma/|Rrs|`` over 450-650 nm -- the measured relative error."""
    w, Rrs, sig = spec['wave'], spec['Rrs'], spec['sigma']
    k = (w > 450) & (w < 650) & np.isfinite(sig) & (np.abs(Rrs) > 1e-5)
    return float(np.median(sig[k] / np.abs(Rrs[k]))) if k.any() else np.nan


def _any_negative(spec):
    """True if any band in 400-900 nm is negative (kept, never filtered)."""
    w, Rrs = spec['wave'], spec['Rrs']
    k = (w >= 400) & (w <= 900) & np.isfinite(Rrs)
    return bool((Rrs[k] < 0).any())


def shape_matrix(grids):
    """L2-normalised spectra over :data:`CLUSTER_RANGE`, for shape clustering.

    L2 normalisation (rather than the more common area normalisation) is used
    because WATERHYPERNET deliberately retains negative reflectance, which can
    drive an integral towards zero and blow up an area-normalised shape.

    Returns
    -------
    tuple
        ``(shapes, ok)`` -- the normalised matrix for the usable rows and the
        boolean mask selecting them.
    """
    lo, hi = CLUSTER_RANGE
    band = (ANALYSIS_WAVE >= lo) & (ANALYSIS_WAVE <= hi)
    sub = grids[:, band]
    ok = np.isfinite(sub).all(axis=1)
    norm = np.linalg.norm(sub[ok], axis=1, keepdims=True)
    ok_idx = np.flatnonzero(ok)
    good = (norm.ravel() > 0)
    ok[ok_idx[~good]] = False
    shapes = sub[ok] / np.linalg.norm(sub[ok], axis=1, keepdims=True)
    return shapes, ok


def cluster_and_sample(index=None, n_pool=N_POOL, n_sample=N_SAMPLE,
                       n_clusters=N_CLUSTERS, seed=SEED, save=True):
    """Read a pool per site, cluster all sites together, sample ~100 per site.

    The clustering is *pooled across sites* so the cluster labels act like
    optical water types shared by the whole network -- that is what makes the
    per-site composition comparable. Each site's ``n_sample`` spectra are then
    allocated across the clusters that site actually occupies, in proportion to
    occupancy (largest remainder, at least one per occupied cluster), and within
    a cluster the picks are spread evenly in time.

    Returns
    -------
    tuple
        ``(pool, sample, grids)`` -- the full pool table (with ``cluster`` and
        ``selected`` columns), the selected subset, and the interpolated
        spectra for the pool rows.
    """
    from sklearn.cluster import KMeans

    if index is None:
        index = pd.read_parquet(os.path.join(out_root(), 'index.parquet'))

    pools, grid_list = [], []
    for site in sorted(index.site.unique()):
        print(f'reading pool: {site}')
        rows, grids = read_pool(index, site, n_pool=n_pool)
        pools.append(rows)
        grid_list.append(grids)

    pool = pd.concat(pools, ignore_index=True)
    grids = np.vstack(grid_list)

    # --- pooled spectral-shape clustering (optical water types)
    shapes, ok = shape_matrix(grids)
    km = KMeans(n_clusters=n_clusters, random_state=seed, n_init=10)
    labels = km.fit_predict(shapes)
    pool['cluster'] = -1
    pool.loc[np.flatnonzero(ok), 'cluster'] = labels

    # --- per-site allocation across the clusters that site occupies
    pool['selected'] = False
    rng = np.random.default_rng(seed)
    for site in pool.site.unique():
        idx = pool.index[(pool.site == site) & (pool.cluster >= 0)]
        if len(idx) == 0:
            continue
        take = min(n_sample, len(idx))
        counts = pool.loc[idx, 'cluster'].value_counts()
        quota = _largest_remainder(counts, take)
        chosen = []
        for cl, n in quota.items():
            members = pool.loc[idx][pool.loc[idx, 'cluster'] == cl]
            members = members.sort_values('datetime')
            if n >= len(members):
                chosen.extend(members.index.tolist())
            else:                                # spread evenly through time
                pos = np.linspace(0, len(members) - 1, n).round().astype(int)
                chosen.extend(members.index[np.unique(pos)].tolist())
        # top up if rounding/uniqueness left us short
        short = take - len(chosen)
        if short > 0:
            rest = [i for i in idx if i not in set(chosen)]
            if rest:
                chosen.extend(rng.choice(rest, size=min(short, len(rest)),
                                         replace=False).tolist())
        pool.loc[chosen, 'selected'] = True

    sample = pool[pool.selected].copy()
    if save:
        root = out_root()
        pool.to_parquet(os.path.join(root, 'pool.parquet'))
        sample.to_parquet(os.path.join(root, 'sample.parquet'))
        np.savez_compressed(os.path.join(root, 'pool_spectra.npz'),
                            wave=ANALYSIS_WAVE, grids=grids,
                            cluster_centers=km.cluster_centers_)
        print(f'wrote pool ({len(pool)}) + sample ({len(sample)}) to {root}')
    return pool, sample, grids


def _largest_remainder(counts, total):
    """Apportion ``total`` across ``counts`` proportionally, min 1 each.

    Largest-remainder (Hare) apportionment, so the per-site sample mirrors the
    site's cluster occupancy without dropping a rare cluster entirely.
    """
    counts = counts[counts > 0]
    if len(counts) == 0:
        return {}
    if total <= len(counts):                    # fewer slots than clusters
        return {cl: 1 for cl in counts.index[:total]}
    exact = counts / counts.sum() * (total - len(counts))
    base = np.floor(exact).astype(int) + 1      # the guaranteed 1 each
    short = total - int(base.sum())
    if short > 0:
        order = (exact - np.floor(exact)).sort_values(ascending=False).index
        for cl in list(order)[:short]:
            base[cl] += 1
    return dict(base)


# --- driver ------------------------------------------------------------------

def main(flg, n_pool=N_POOL, n_sample=N_SAMPLE):
    """Stage driver: ``1`` index, ``2`` pool+cluster+sample."""
    flg = int(flg)
    if flg == 1:
        df = build_index()
        print(df.groupby(['site', 'system']).size())
        print(f'total {len(df)} spectra')
    elif flg == 2:
        cluster_and_sample(n_pool=n_pool, n_sample=n_sample)


def _cli(argv=None):
    """CLI: ``whn_explore.py <flg> [--n-pool N] [--n-sample N]``."""
    import argparse

    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('flg', nargs='?', type=int, default=0,
                   help='stage: 1 index, 2 pool+cluster+sample (0 = no-op)')
    p.add_argument('--n-pool', type=int, default=N_POOL,
                   help='candidate spectra read per site (stage 2)')
    p.add_argument('--n-sample', type=int, default=N_SAMPLE,
                   help='spectra kept per site (stage 2)')
    a = p.parse_args(argv)
    main(a.flg, n_pool=a.n_pool, n_sample=a.n_sample)


if __name__ == '__main__':
    _cli()

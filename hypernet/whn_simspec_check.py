"""Full-archive check of Similarity-Spectrum over-subtraction at the dark sites.

The pooled spectral-shape clustering in :mod:`whn_figures` threw up a small
cluster (9 of 4,400) of spectra that are negative across most of the visible.
All of them were HYPSTAR records at the three darkest sites, all with
``quality_flag = 0``. This script tests the diagnosis against **every** file at
those sites rather than a sample, so the rate quoted to the data provider is a
real one.

The test, per file, over 440-600 nm (the visible peak, away from the <400 nm and
>900 nm ranges the release notes call unreliable):

- ``med_ref``  -- median of ``reflectance`` (Similarity-Spectrum corrected)
- ``med_nosc`` -- median of ``reflectance_nosc`` (uncorrected)
- **over-subtraction** is flagged when ``med_ref <= 0 < med_nosc``: the
  correction has removed more than the entire water signal, turning a plausible
  positive spectrum into a non-physical one.

Values are reported as rho_w straight from the file (divide by pi for Rrs); the
sign test is unaffected by that factor.

Run::

    python whn_simspec_check.py              # the three dark sites
    python whn_simspec_check.py --sites all  # every HYPSTAR site (slow, ~30 min)
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

from hypernet.whn_explore import out_root, whn_root  # noqa: F401  (whn_root: path check)

#: The three darkest sites by median Rrs(560) -- where cluster 5 came from.
DARK_SITES = ('THFR_H', 'BEFR_H', 'WRUK_H')

#: Visible band over which the sign test is made.
TEST_RANGE = (440.0, 600.0)


def flag_oversubtraction(med_ref, med_nosc):
    """True where the correction removed more than the entire water signal.

    The corrected product has gone non-positive across the visible while the
    uncorrected one is still positive -- which cannot be physical, and is the
    signature seen at the darkest sites.
    """
    return (np.asarray(med_ref) <= 0) & (np.asarray(med_nosc) > 0)


def check_file(path):
    """Return the per-file statistics used by the over-subtraction test.

    Returns
    -------
    dict or None
        ``med_ref``, ``med_nosc``, ``quality_flag``, ``sza``, ``neg_frac``
        (fraction of 400-900 nm bands where the corrected product is negative),
        and ``n_valid``/``n_total`` scans. ``None`` if the file cannot be read.
    """
    import netCDF4

    try:
        ds = netCDF4.Dataset(path)
    except Exception:
        return None
    try:
        ds.set_auto_mask(False)
        w = np.asarray(ds.variables['wavelength'][:], dtype=float).ravel()
        ref = np.asarray(ds.variables['reflectance'][:], dtype=float).ravel()
        nsc = np.asarray(ds.variables['reflectance_nosc'][:], dtype=float).ravel()
        qf = int(np.asarray(ds.variables['quality_flag'][:]).ravel()[0])
        sza = float(np.asarray(ds.variables['solar_zenith_angle'][:]).ravel()[0])
        nv = int(np.asarray(ds.variables['n_valid_scans'][:]).ravel()[0])
        nt = int(np.asarray(ds.variables['n_total_scans'][:]).ravel()[0])
    except Exception:
        return None
    finally:
        ds.close()

    fill = 9.96921e36                       # HYPSTAR _FillValue
    ref = np.where(np.abs(ref) > fill / 2, np.nan, ref)
    nsc = np.where(np.abs(nsc) > fill / 2, np.nan, nsc)

    lo, hi = TEST_RANGE
    k = (w >= lo) & (w <= hi)
    vis = (w >= 400) & (w <= 900)
    if not k.any():
        return None

    return {
        'med_ref': float(np.nanmedian(ref[k])),
        'med_nosc': float(np.nanmedian(nsc[k])),
        'quality_flag': qf,
        'sza': sza,
        'neg_frac': float(np.nanmean(ref[vis] < 0)),
        'n_valid': nv,
        'n_total': nt,
    }


def scan(sites=DARK_SITES, save=True, verbose=True):
    """Run :func:`check_file` over every file at ``sites``.

    Returns
    -------
    pandas.DataFrame
        One row per file, with an ``oversub`` boolean column.
    """
    index = pd.read_parquet(os.path.join(out_root(), 'index.parquet'))
    index = index[index.site.isin(sites)]

    rows = []
    for site, g in index.groupby('site'):
        if verbose:
            print(f'scanning {site}: {len(g)} files')
        for i, r in enumerate(g.itertuples()):
            stats = check_file(r.path)
            if stats is None:
                continue
            stats.update(site=site, datetime=r.datetime, path=r.path)
            rows.append(stats)
            if verbose and (i + 1) % 1000 == 0:
                print(f'  {site}: {i + 1}/{len(g)}')

    df = pd.DataFrame(rows)
    df['oversub'] = flag_oversubtraction(df.med_ref, df.med_nosc)
    if save:
        df.to_parquet(os.path.join(out_root(), 'simspec_check.parquet'))
    return df


def report(df):
    """Print the per-site over-subtraction rate and the QC cross-tab."""
    print(f'\n{"site":8s} {"n":>6s} {"oversub":>8s} {"rate %":>7s} '
          f'{"qf=0 of those":>14s} {"median sza":>11s}')
    for site, g in df.groupby('site'):
        bad = g[g.oversub]
        print(f'{site:8s} {len(g):6d} {len(bad):8d} {100 * len(bad) / len(g):7.3f} '
              f'{int((bad.quality_flag == 0).sum()):14d} '
              f'{(bad.sza.median() if len(bad) else np.nan):11.1f}')

    bad = df[df.oversub]
    print(f'\ntotal {len(bad)} of {len(df)} ({100 * len(bad) / len(df):.3f} %)')
    if len(bad):
        print(f'all quality_flag == 0 : {bool((bad.quality_flag == 0).all())}')
        print(f'all scans 6/6         : '
              f'{bool((bad.n_valid == bad.n_total).all())}')
        print(f'sza range of flagged  : {bad.sza.min():.1f} - {bad.sza.max():.1f} deg')
        print(f'median rho_w_nosc(440-600) of flagged : {bad.med_nosc.median():.5f}')


def main(sites=DARK_SITES):
    """Scan, save and report."""
    df = scan(sites)
    report(df)
    return df


if __name__ == '__main__':
    import argparse

    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--sites', default='dark',
                   help="'dark' (default: THFR_H, BEFR_H, WRUK_H) or 'all' HYPSTAR")
    a = p.parse_args()

    if a.sites == 'all':
        idx = pd.read_parquet(os.path.join(out_root(), 'index.parquet'))
        target = tuple(sorted(idx[idx.system == 'HYPSTAR'].site.unique()))
    else:
        target = DARK_SITES
    main(target)

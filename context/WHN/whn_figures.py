"""Figures and summary table for the WATERHYPERNET Release 2 exploration.

Consumes the tables written by :mod:`whn_explore` (stages 1-2) and produces, in
``context/WHN``:

- ``summary_table.csv`` + ``summary_table.md`` -- one row per site.
- ``figs/fig_sites_map.png``      -- where the 11 sites are.
- ``figs/fig_data_volume.png``    -- how much data each site has, and when.
- ``figs/fig_median_spectra.png`` -- per-site median Rrs with percentile
                                     envelopes, on one common y axis.
- ``figs/fig_median_spectra_autoscale.png`` -- the same, per-panel autoscaled.
- ``figs/fig_clusters.png``       -- pooled spectral-shape clusters (OWT-like)
                                     and each site's composition.
- ``figs/fig_band_timeseries.png``-- Rrs(490/560/665) against time, per site.

All spectra plotted here are on the display grid :data:`whn_explore.ANALYSIS_WAVE`
and expressed as **Rrs = rho_w/pi** [1/sr]; see :mod:`whn_explore` for the
product and units conventions. Run after stages 1-2::

    python whn_figures.py
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')                            # headless
import matplotlib.pyplot as plt                  # noqa: E402
from matplotlib import gridspec                  # noqa: E402

from whn_explore import (ANALYSIS_WAVE, CLUSTER_RANGE, N_CLUSTERS, NOSC_SITES,
                         TS_BANDS, fig_dir, out_root, product_for)  # noqa: E402

#: Consistent colour per measuring system across every figure.
SYS_COLOR = {'HYPSTAR': '#1f78b4', 'PANTHYR': '#e66101'}

#: Human-readable site names from the release notes.
SITE_NAME = {
    'BEFR_H': 'Berre lagoon (FR)',      'CBUS_P': 'Chesapeake Bay (US)',
    'GAIT_H': 'Lake Garda (IT)',        'LPAR_H': 'Rio de la Plata (AR)',
    'MAFR_H': 'Gironde/MAGEST (FR)',    'O1BE_P': 'Oostende RT1 (BE)',
    'TBBE_P': 'Thornton Bank (BE)',     'THFR_H': 'Thau lagoon (FR)',
    'VEIT_H': 'Acqua Alta AAOT (IT)',   'VEIT_P': 'Acqua Alta AAOT (IT)',
    'WRUK_H': 'Wraysbury reservoir (UK)',
}


def load_products():
    """Load the stage-1/2 artifacts.

    Returns
    -------
    tuple
        ``(index, pool, grids, centers)`` -- the full filename index, the pool
        table (with cluster labels and the ``selected`` flag), the interpolated
        spectra for the pool rows, and the KMeans cluster centres.
    """
    root = out_root()
    index = pd.read_parquet(os.path.join(root, 'index.parquet'))
    pool = pd.read_parquet(os.path.join(root, 'pool.parquet'))
    with np.load(os.path.join(root, 'pool_spectra.npz')) as npz:
        grids = npz['grids']
        centers = npz['cluster_centers']
    return index, pool, grids, centers


# --- summary table -----------------------------------------------------------

def summary_table(index, pool, grids, save=True):
    """One row per site: extent, geometry, product, grid, and Rrs statistics.

    Counts and date spans come from the full 56,669-file index; the optical
    statistics come from the sampled pool (the archive is far too large to read
    in full, and the pool is spread uniformly through each site's history).
    """
    rows = []
    for site, g in index.groupby('site'):
        p = pool[pool.site == site]
        sel = p[p.selected]
        gsel = grids[p.index[p.selected.values]] if len(sel) else np.zeros((0, len(ANALYSIS_WAVE)))
        refl_var, _ = product_for(site)
        i560 = int(np.argmin(np.abs(ANALYSIS_WAVE - 560.0)))
        i490 = int(np.argmin(np.abs(ANALYSIS_WAVE - 490.0)))
        with np.errstate(invalid='ignore', divide='ignore'):
            ratio = (np.nanmedian(gsel[:, i490] / gsel[:, i560])
                     if len(gsel) else np.nan)
        rows.append({
            'site': site,
            'name': SITE_NAME.get(site, ''),
            'system': p.system.iloc[0] if len(p) else '',
            'N_files': len(g),
            'first': g.datetime.min().date(),
            'last': g.datetime.max().date(),
            'months': g.datetime.dt.month.nunique(),
            'lat': round(float(p.lat.median()), 5) if len(p) else np.nan,
            'lon': round(float(p.lon.median()), 5) if len(p) else np.nan,
            'product': refl_var,
            'n_serial': int(p.serial.nunique()) if len(p) else 0,
            # A grid is a (length, start) pair: MAFR's two instruments share a
            # band count but start 0.18 nm apart, so length alone undercounts.
            'n_grids': (len(set(zip(p.n_bands, p.wave_min.round(2))))
                        if len(p) else 0),
            'n_bands': int(p.n_bands.median()) if len(p) else 0,
            'sza_min': round(float(p.sza.min()), 1) if len(p) else np.nan,
            'sza_max': round(float(p.sza.max()), 1) if len(p) else np.nan,
            'Rrs560_med': float(np.nanmedian(gsel[:, i560])) if len(gsel) else np.nan,
            'blue_green': float(ratio),
            'rel_sigma': round(float(p.rel_sigma.median()), 4) if len(p) else np.nan,
            'pct_negative': round(100 * float(p.any_negative.mean()), 1) if len(p) else np.nan,
            'n_pool': len(p),
            'n_sample': len(sel),
        })
    df = pd.DataFrame(rows).sort_values('site').reset_index(drop=True)

    if save:
        here = os.path.dirname(os.path.abspath(__file__))
        df.to_csv(os.path.join(here, 'summary_table.csv'), index=False)
        with open(os.path.join(here, 'summary_table.md'), 'w') as fh:
            fh.write(df.to_markdown(index=False, floatfmt='.4g'))
    return df


# --- figures -----------------------------------------------------------------

#: Extent of the Europe inset -- 9 of the 11 sites fall inside it.
EUROPE_EXTENT = (-7.0, 17.0, 41.0, 55.0)

#: Manual label placement (dlon, dlat, ha) per site group. The sites are far too
#: clustered -- Oostende and Thornton Bank are 33 km apart, Acqua Alta carries
#: two systems -- for automatic offsets to stay legible.
LABEL_OFFSET = {
    'WRUK_H':        (-0.5,  0.35, 'right'),
    'O1BE_P':        (-0.5, -0.75, 'right'),
    'TBBE_P':        ( 0.5,  0.35, 'left'),
    'MAFR_H':        (-0.5,  0.35, 'right'),
    'THFR_H':        (-0.5, -0.85, 'right'),
    'BEFR_H':        ( 0.5,  0.25, 'left'),
    'GAIT_H':        (-0.5,  0.45, 'right'),
    'VEIT_H+VEIT_P': (-0.6, -0.95, 'right'),
    'CBUS_P':        ( 4.0,  1.50, 'left'),
    'LPAR_H':        ( 4.0, -1.50, 'left'),
}


def _site_groups(summary):
    """Group co-located sites so overlapping markers get one shared label.

    VEIT_H and VEIT_P are the same platform (Acqua Alta) carrying both a HYPSTAR
    and a PANTHYR, so they plot on top of each other and must share a label.
    """
    groups = {}
    for r in summary.itertuples():
        if not np.isfinite(r.lat):
            continue
        key = (round(r.lat, 2), round(r.lon, 2))
        groups.setdefault(key, []).append(r)
    return groups


def fig_sites_map(summary):
    """World map of the 11 sites, with a Europe inset (9 sites are clustered there).

    Markers are only drawn where they fall inside the axis extent -- cartopy
    does not clip annotation text, so an out-of-extent label would otherwise be
    drawn outside the axes and distort the saved figure.
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    groups = _site_groups(summary)
    nmax = summary.N_files.max()

    fig = plt.figure(figsize=(13, 4.2))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.55, 1], wspace=0.06)

    # Scales are pinned explicitly: cartopy otherwise auto-selects a finer scale
    # for a zoomed extent and blocks on a Natural Earth download the first time.
    panels = ((0, (-170, 180, -58, 78), '110m',
               '(a) WATERHYPERNET Release 2 — 11 sites', True),
              (1, EUROPE_EXTENT, '50m', '(b) European sites', False))

    for ax_i, extent, scale, title, is_global in panels:
        ax = fig.add_subplot(gs[ax_i], projection=ccrs.PlateCarree())
        ax.set_extent(extent, crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND.with_scale(scale), facecolor='#eeeae4')
        ax.add_feature(cfeature.OCEAN.with_scale(scale), facecolor='#d6e6f2')
        ax.add_feature(cfeature.COASTLINE.with_scale(scale), linewidth=0.4)
        ax.add_feature(cfeature.BORDERS.with_scale(scale), linewidth=0.25,
                       edgecolor='0.6')
        ax.set_title(title, fontsize=11)

        lo_x, hi_x, lo_y, hi_y = extent
        ex0, ex1, ey0, ey1 = EUROPE_EXTENT
        for (lat, lon), members in groups.items():
            if not (lo_x < lon < hi_x and lo_y < lat < hi_y):
                continue                        # never label outside the axes
            for k, r in enumerate(members):     # co-located: nudge the marker
                ax.scatter(lon + (0.5 if is_global else 0.07) * k, lat,
                           s=22 + 300 * r.N_files / nmax,
                           c=SYS_COLOR[r.system], edgecolor='k', linewidth=0.5,
                           alpha=0.9, zorder=5, clip_on=True,
                           transform=ccrs.PlateCarree())
            label = '+'.join(m.site for m in members)
            # On the world map, label only what the inset does not cover.
            if is_global and (ex0 < lon < ex1 and ey0 < lat < ey1):
                continue
            dlon, dlat, ha = LABEL_OFFSET.get(label, (0.5, 0.3, 'left'))
            ax.text(lon + dlon, lat + dlat, label, fontsize=7.5, ha=ha,
                    transform=ccrs.PlateCarree(), zorder=6, clip_on=True,
                    bbox=dict(fc='white', ec='none', alpha=0.7, pad=0.7))

        if is_global:                           # show where the inset looks
            ax.add_patch(plt.Rectangle((ex0, ey0), ex1 - ex0, ey1 - ey0,
                                       fill=False, edgecolor='#c0392b',
                                       linewidth=1.2, zorder=7,
                                       transform=ccrs.PlateCarree()))
            ax.text(ex0, ey1 + 4, 'inset (b) — 9 sites', color='#c0392b',
                    fontsize=8, ha='left', transform=ccrs.PlateCarree(),
                    zorder=7)

    handles = [plt.Line2D([], [], marker='o', ls='', color=c, mec='k', label=s)
               for s, c in SYS_COLOR.items()]
    fig.legend(handles=handles, loc='lower center', ncol=2, frameon=False,
               fontsize=9, title='marker area ∝ number of spectra',
               title_fontsize=9, bbox_to_anchor=(0.5, 0.02))
    fig.subplots_adjust(bottom=0.12, top=0.92)
    _save(fig, 'fig_sites_map.png')


def fig_data_volume(index, summary):
    """How much data each site holds (bar) and when it was acquired (heatmap)."""
    fig = plt.figure(figsize=(12, 7))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1.35], hspace=0.32)

    # --- (a) total spectra per site
    ax = fig.add_subplot(gs[0])
    s = summary.sort_values('N_files', ascending=False)
    ax.bar(s.site, s.N_files, color=[SYS_COLOR[x] for x in s.system],
           edgecolor='k', linewidth=0.4)
    for x, (n, sys_) in enumerate(zip(s.N_files, s.system)):
        ax.text(x, n + 250, f'{n:,}', ha='center', fontsize=8)
    ax.set_ylabel('spectra in Release 2')
    ax.set_title(f'(a) Archive volume — {summary.N_files.sum():,} spectra, 11 sites')
    ax.set_ylim(0, summary.N_files.max() * 1.15)
    handles = [plt.Line2D([], [], marker='s', ls='', color=c, label=k)
               for k, c in SYS_COLOR.items()]
    ax.legend(handles=handles, frameon=False, loc='upper right')

    # --- (b) month x site coverage
    ax = fig.add_subplot(gs[1])
    idx = index.copy()
    idx['ym'] = idx.datetime.dt.to_period('M')
    piv = idx.pivot_table(index='site', columns='ym', values='path',
                          aggfunc='count').fillna(0)
    piv = piv.reindex(sorted(piv.index, reverse=True))
    im = ax.imshow(piv.values, aspect='auto', cmap='viridis',
                   norm=matplotlib.colors.LogNorm(vmin=1, vmax=piv.values.max()))
    ax.set_yticks(range(len(piv.index)), piv.index, fontsize=8)
    cols = [str(c) for c in piv.columns]
    step = max(1, len(cols) // 16)
    ax.set_xticks(range(0, len(cols), step), cols[::step], rotation=45,
                  ha='right', fontsize=7)
    ax.set_title('(b) Temporal coverage — spectra per month (log colour; white = no data)')
    fig.colorbar(im, ax=ax, pad=0.01, label='spectra / month')
    _save(fig, 'fig_data_volume.png')


def fig_median_spectra(pool, grids, shared=True, name=None):
    """Per-site median Rrs with 25-75 and 10-90 percentile envelopes.

    Parameters
    ----------
    shared : bool, optional
        True (default) puts every panel on one common y axis, so the relative
        brightness of the sites reads at a glance -- the network spans a factor
        of ~12 in median Rrs(560), which is the single most useful thing a
        reader can take from this figure. False lets each panel autoscale,
        which preserves the spectral shape of the dark sites; the two are
        produced as companion figures.
    name : str or None, optional
        Output filename; defaults to the appropriate one for ``shared``.
    """
    sites = sorted(pool.site.unique())
    fig, axes = plt.subplots(3, 4, figsize=(15, 9), sharex=True,
                             sharey=bool(shared))
    axes = axes.ravel()

    lims = [np.inf, -np.inf]
    for ax, site in zip(axes, sites):
        rows = pool.index[(pool.site == site).values & pool.selected.values]
        g = grids[rows]
        med = np.nanmedian(g, axis=0)
        c = SYS_COLOR[pool.loc[rows[0], 'system']]
        for lo, hi, a in ((10, 90, 0.18), (25, 75, 0.32)):
            band_lo = np.nanpercentile(g, lo, axis=0)
            band_hi = np.nanpercentile(g, hi, axis=0)
            ax.fill_between(ANALYSIS_WAVE, band_lo, band_hi, color=c, alpha=a,
                            lw=0)
            lims = [min(lims[0], np.nanmin(band_lo)),
                    max(lims[1], np.nanmax(band_hi))]
        ax.plot(ANALYSIS_WAVE, med, color=c, lw=1.6)
        ax.axhline(0, color='0.5', lw=0.6, ls=':')
        ax.set_title(f'{site} — {SITE_NAME.get(site, "")}\n'
                     f'n={len(rows)}, {product_for(site)[0]}', fontsize=8.5)
        ax.tick_params(labelsize=8)
        if site in NOSC_SITES:
            ax.patch.set_facecolor('#fff6ef')

    if shared:
        pad = 0.04 * (lims[1] - lims[0])
        for ax in axes[:len(sites)]:
            ax.set_ylim(lims[0] - pad, lims[1] + pad)
            ax.tick_params(labelleft=True)        # sharey hides inner labels

    ncol = axes.size // 3
    for j in range(ncol):                        # sharex hides all but the last
        col = [i for i in range(len(sites)) if i % ncol == j]
        if col:
            axes[col[-1]].tick_params(labelbottom=True)
    for ax in axes[len(sites):]:
        ax.axis('off')

    note = ('median with 25–75 and\n10–90 percentile envelopes\n\n'
            'shaded panels use\nreflectance_nosc\n(turbid sites)')
    note += ('\n\ncommon y axis:\nbrightness compares\ndirectly between sites'
             if shared else
             '\n\nper-panel y axis:\nshape detail at the\ndark sites')
    axes[len(sites)].text(0.05, 0.95, note, fontsize=9, va='top',
                          transform=axes[len(sites)].transAxes)

    fig.supxlabel('wavelength [nm]')
    fig.supylabel('R$_{rs}$ = ρ$_w$/π  [sr$^{-1}$]')
    kind = 'common y axis' if shared else 'per-panel y axis'
    fig.suptitle('Per-site Rrs distribution — WATERHYPERNET Release 2 '
                 f'(~100 optically-sampled spectra per site; {kind})',
                 fontsize=12)
    fig.tight_layout(rect=(0.01, 0.01, 1, 0.96))
    _save(fig, name or ('fig_median_spectra.png' if shared
                        else 'fig_median_spectra_autoscale.png'))


def fig_clusters(pool, grids, centers):
    """Pooled spectral-shape clusters and each site's composition."""
    lo, hi = CLUSTER_RANGE
    band = (ANALYSIS_WAVE >= lo) & (ANALYSIS_WAVE <= hi)
    wave = ANALYSIS_WAVE[band]
    cmap = plt.get_cmap('tab10')

    fig = plt.figure(figsize=(14, 9))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1.1],
                           hspace=0.32, wspace=0.22)

    # --- (a) cluster mean shapes
    ax = fig.add_subplot(gs[0, 0])
    for k in range(len(centers)):
        n = int((pool.cluster == k).sum())
        ax.plot(wave, centers[k], color=cmap(k), lw=1.8,
                label=f'C{k} (n={n})')
    ax.axhline(0, color='0.5', lw=0.6, ls=':')
    ax.set_xlabel('wavelength [nm]')
    ax.set_ylabel('L2-normalised Rrs shape')
    ax.set_title(f'(a) Pooled spectral-shape clusters (k={len(centers)}, {lo:.0f}–{hi:.0f} nm)')
    ax.legend(fontsize=7.5, ncol=2, frameon=False)

    # --- (b) site x cluster composition
    ax = fig.add_subplot(gs[0, 1])
    comp = (pool[pool.cluster >= 0]
            .pivot_table(index='site', columns='cluster', values='path',
                         aggfunc='count').fillna(0))
    frac = comp.div(comp.sum(axis=1), axis=0)
    bottom = np.zeros(len(frac))
    for k in frac.columns:
        ax.barh(frac.index, frac[k], left=bottom, color=cmap(k),
                edgecolor='w', linewidth=0.4, label=f'C{k}')
        bottom += frac[k].values
    ax.set_xlim(0, 1)
    ax.set_xlabel('fraction of the site’s pool')
    ax.set_title('(b) Optical composition by site — colours as in (a)')
    ax.tick_params(labelsize=8)

    # --- (c) all sampled spectra, coloured by cluster
    ax = fig.add_subplot(gs[1, 0])
    sel = pool.index[pool.selected.values]
    for k in range(len(centers)):
        rows = [i for i in sel if pool.loc[i, 'cluster'] == k]
        if not rows:
            continue
        ax.plot(ANALYSIS_WAVE, np.nanmedian(grids[rows], axis=0),
                color=cmap(k), lw=1.6, label=f'C{k}')
        ax.fill_between(ANALYSIS_WAVE,
                        np.nanpercentile(grids[rows], 25, axis=0),
                        np.nanpercentile(grids[rows], 75, axis=0),
                        color=cmap(k), alpha=0.15, lw=0)
    ax.axhline(0, color='0.5', lw=0.6, ls=':')
    ax.set_xlabel('wavelength [nm]')
    ax.set_ylabel('R$_{rs}$ [sr$^{-1}$]')
    ax.set_title('(c) Cluster median Rrs in absolute units (25–75 %)')
    ax.legend(fontsize=7.5, ncol=2, frameon=False)

    # --- (d) magnitude vs blue/green colour, coloured by cluster
    ax = fig.add_subplot(gs[1, 1])
    i560 = int(np.argmin(np.abs(ANALYSIS_WAVE - 560.0)))
    i490 = int(np.argmin(np.abs(ANALYSIS_WAVE - 490.0)))
    with np.errstate(invalid='ignore', divide='ignore'):
        x = grids[:, i490] / grids[:, i560]
    y = grids[:, i560]
    ok = np.isfinite(x) & np.isfinite(y) & (pool.cluster.values >= 0)
    ax.scatter(x[ok], y[ok], c=[cmap(k) for k in pool.cluster.values[ok]],
               s=7, alpha=0.5, lw=0)
    ax.set_xlim(-0.25, 2.6)                      # a handful of outliers clip
    ax.set_yscale('symlog', linthresh=1e-4)
    ax.set_xlabel('R$_{rs}$(490) / R$_{rs}$(560)   [blue-green ratio]')
    ax.set_ylabel('R$_{rs}$(560) [sr$^{-1}$]')
    n_clip = int((x[ok] > 2.6).sum() + (x[ok] < -0.25).sum())
    ax.set_title(f'(d) Optical spread of the pool ({n_clip} outliers clipped)')
    ax.grid(alpha=0.25, lw=0.5)

    fig.suptitle('Spectral-shape clustering — the axis used to draw the '
                 '~100 spectra per site', fontsize=12)
    _save(fig, 'fig_clusters.png')


def fig_band_timeseries(pool, grids):
    """Rrs at 490/560/665 nm against acquisition date, per site."""
    sites = sorted(pool.site.unique())
    fig, axes = plt.subplots(6, 2, figsize=(14, 13), sharex=True)
    axes = axes.ravel()
    colors = {490.0: '#2c7fb8', 560.0: '#41ab5d', 665.0: '#d7301f'}

    for ax, site in zip(axes, sites):
        rows = pool.index[(pool.site == site).values]
        sub = pool.loc[rows].sort_values('datetime')
        g = grids[sub.index]
        for b in TS_BANDS:
            i = int(np.argmin(np.abs(ANALYSIS_WAVE - b)))
            ax.plot(sub.datetime, g[:, i], '.', ms=2.6, color=colors[b],
                    label=f'{b:.0f} nm', alpha=0.8)
        ax.axhline(0, color='0.5', lw=0.6, ls=':')
        ax.set_title(f'{site} — {SITE_NAME.get(site, "")}', fontsize=9)
        ax.tick_params(labelsize=7.5)
        ax.set_yscale('symlog', linthresh=1e-4)
        # Negatives are real and kept, but they are rare -- do not give the
        # negative half of the symlog axis more room than the data needs.
        ax.set_ylim(bottom=-5e-4)

    for j in range(2):                           # sharex hides all but the last
        col = [i for i in range(len(sites)) if i % 2 == j]
        if col:
            axes[col[-1]].tick_params(labelbottom=True)
    for ax in axes[len(sites):]:
        ax.axis('off')
    axes[0].legend(fontsize=7.5, ncol=3, frameon=False)
    fig.supylabel('R$_{rs}$ [sr$^{-1}$]  (symlog)')
    fig.suptitle('Band time series over the uniformly-spread pool '
                 '(≤400 spectra per site)', fontsize=12)
    fig.tight_layout(rect=(0.01, 0.01, 1, 0.97))
    _save(fig, 'fig_band_timeseries.png')


def _save(fig, name):
    """Write ``fig`` into ``context/WHN/figs`` and report the path."""
    path = os.path.join(fig_dir(), name)
    fig.savefig(path, dpi=140, bbox_inches='tight')
    plt.close(fig)
    print(f'wrote {path}')


def main():
    """Build the summary table and every figure."""
    index, pool, grids, centers = load_products()
    summary = summary_table(index, pool, grids)
    print(summary.to_string(index=False))

    fig_sites_map(summary)
    fig_data_volume(index, summary)
    fig_median_spectra(pool, grids, shared=True)
    fig_median_spectra(pool, grids, shared=False)
    fig_clusters(pool, grids, centers)
    fig_band_timeseries(pool, grids)


if __name__ == '__main__':
    main()

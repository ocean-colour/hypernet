"""Pick candidate WATERHYPERNET sequences for the wiggles data request to Kevin.

Instrument periods (site x instrument x calibration) are derived from the
per-file attribute survey written by ``instrument_timeline.py``.  Within each
period, sequences are drawn from the RELEASE_2 index so that they span solar
zenith angle (tertiles, SZA <= 75 deg) and season (month).  SZA and the
instrument id are read from each L2 file, so every pick is verified to come
from the intended instrument.

Each request row gets N primary candidates plus N//2 spares, so Kevin can swap
or tag clear/overcast sky without us having to go round again.

Writes ``docs/wiggles_data_request.csv``.

    python wavecal/select_sequences.py
"""
import os

import netCDF4
import numpy as np
import pandas as pd

SEED = 20260928
SZA_MAX = 75.   # deg; beyond this Ed is low and the sky-glint term dominates
OUT_ROOT = os.path.join(os.environ['OS_COLOR'], 'hypernet')
OUT_CSV = os.path.join(os.path.dirname(__file__), '..', 'docs',
                       'wiggles_data_request.csv')

# (row label, site, instrument, which period of that instrument at that site
#  in time order (0 = first, -1 = last, None = all), N primary)
ROWS = [
    ('VEIT 122304 pre-recal', 'VEIT_H', 'HYPSTAR_122304', 0, 20),
    ('VEIT 122305', 'VEIT_H', 'HYPSTAR_122305', None, 20),
    ('VEIT 122304 post-recal', 'VEIT_H', 'HYPSTAR_122304', -1, 20),
    ('BEFR 122302', 'BEFR_H', 'HYPSTAR_122302', None, 20),
    ('THFR 122302', 'THFR_H', 'HYPSTAR_122302', None, 20),
    ('MAFR 121231', 'MAFR_H', 'HYPSTAR_121231', None, 22),
    ('MAFR 122303', 'MAFR_H', 'HYPSTAR_122303', None, 8),
    ('GAIT 121222 pre-recal', 'GAIT_H', 'HYPSTAR_121222', 0, 7),
    ('GAIT 120242', 'GAIT_H', 'HYPSTAR_120242', None, 6),
    ('GAIT 121222 post-recal', 'GAIT_H', 'HYPSTAR_121222', -1, 7),
]
OPTIONAL = {'GAIT_H'}


def periods(timeline):
    """Contiguous runs of (instrument, cal) per site -> start/end dates."""
    t = timeline.sort_values(['site', 'datetime']).copy()
    key = t['L'] + '|' + t['cal']
    t['run'] = (key != key.groupby(t['site']).shift()).cumsum()
    return (t.groupby('run')
            .agg(site=('site', 'first'), inst=('L', 'first'),
                 cal=('cal', 'first'), start=('datetime', 'min'),
                 end=('datetime', 'max'))
            .reset_index(drop=True))


def read_attrs(path):
    """Mean SZA and the file's own instrument id (to verify the period)."""
    try:
        with netCDF4.Dataset(path) as ds:
            return (float(np.nanmean(ds['solar_zenith_angle'][:])),
                    getattr(ds, 'system_id', ''))
    except Exception:
        return np.nan, ''


def pick(cands, n, rng):
    """Spread n picks over SZA tertiles x month, round-robin."""
    cands = cands.dropna(subset=['sza']).copy()
    cands['sza_bin'] = pd.qcut(cands['sza'], 3, labels=['low', 'mid', 'high'],
                               duplicates='drop')
    cands['month'] = cands['datetime'].dt.month
    cells = [g.sample(frac=1, random_state=rng.integers(1e9))
             for _, g in cands.groupby(['sza_bin', 'month'], observed=True)]
    rng.shuffle(cells)
    out, i = [], 0
    while len(out) < n and any(len(c) for c in cells):
        c = cells[i % len(cells)]
        if len(c):
            out.append(c.iloc[0])
            cells[i % len(cells)] = c.iloc[1:]
        i += 1
    return pd.DataFrame(out)


def main():
    rng = np.random.default_rng(SEED)
    idx = pd.read_parquet(os.path.join(OUT_ROOT, 'whn_explore', 'index.parquet'))
    idx = idx[idx.system == 'HYPSTAR']
    # One row per sequence: prefer the L2A file when both L2A and L2B exist.
    idx = (idx.assign(l2a=idx.path.str.contains('_L2A_'))
           .sort_values('l2a', ascending=False)
           .drop_duplicates(['site', 'datetime', 'azimuth']))
    tl = pd.read_csv(os.path.join(OUT_ROOT, 'wavecal', 'instrument_timeline.csv'),
                     parse_dates=['datetime'])
    per = periods(tl)

    rows = []
    for label, site, inst, which, n in ROWS:
        p = per[(per.site == site) & (per.inst == inst)].sort_values('start')
        if p.empty:
            print(f'!! no period for {label}')
            continue
        p = p if which is None else p.iloc[[which]]
        # The timeline samples every 20th file, so pad each period to the
        # neighbouring sample to avoid losing its edges.
        sel = np.zeros(len(idx), bool)
        for _, r in p.iterrows():
            sel |= ((idx.site == site) & (idx.datetime >= r.start)
                    & (idx.datetime <= r.end)).values
        cands = idx[sel]
        m = min(len(cands), 8 * n)
        cands = cands.sample(m, random_state=rng.integers(1e9)).copy()
        cands['sza'], cands['file_inst'] = zip(*[read_attrs(f) for f in cands.path])
        cands = cands[(cands.file_inst == inst) & (cands.sza <= SZA_MAX)]
        got = pick(cands, n + n // 2, rng)
        got['priority'] = ['primary'] * min(n, len(got)) + \
            ['spare'] * max(0, len(got) - n)
        got['row'] = label
        got['instrument'] = inst
        got['cal_dates_rad_irr'] = '; '.join(p.cal.unique())
        got['optional'] = site in OPTIONAL
        rows.append(got)
        print(f'{label:24s} period {p.start.min():%Y-%m-%d}..{p.end.max():%Y-%m-%d} '
              f'cands {len(cands):4d} picked {len(got)} '
              f'SZA {got.sza.min():.0f}-{got.sza.max():.0f} '
              f'months {got.datetime.dt.month.nunique()}')

    df = pd.concat(rows)
    df['sequence_time'] = df.datetime.dt.strftime('%Y%m%dT%H%M')
    df['sza'] = df.sza.round(1)
    df['file'] = df.path.map(os.path.basename)
    df['sky'] = ''   # for Kevin: clear / overcast / broken
    cols = ['row', 'priority', 'optional', 'site', 'instrument',
            'cal_dates_rad_irr', 'sequence_time', 'azimuth', 'sza', 'sky', 'file']
    df = df.sort_values(['row', 'priority', 'sequence_time'])[cols]
    df.to_csv(OUT_CSV, index=False)
    print(f'wrote {len(df)} rows ({(df.priority == "primary").sum()} primary) '
          f'to {os.path.normpath(OUT_CSV)}')


if __name__ == '__main__':
    main()

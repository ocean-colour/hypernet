"""Phase 0, task 8e: the sky index of every requested sequence, from Release 2.

Promoted from ``wavecal/ld_ed_750_request.py``.  For each row of
``docs/wiggles_data_request.csv`` (left as sent), reads the Release 2 L2B file
and computes ``ld_ed_750`` = Ld(750)/Ed(750) (``whn_l1a.sky_index``) and the
sky class (``whn_l1a.sky_class``: clear < 0.05, cloudy, overcast >= 0.25;
"broken" needs the L1A Ed scan CV, so it cannot be assigned here).

Writes ``hypernet/wiggles/phase0_sky_index.csv`` (request row, site,
sequence_time, azimuth, sza, raa, instrument, ld_ed_750, sky) and prints the
class counts per request row, the quantiles, the clear floor by SZA and site,
and whether every instrument has clear and cloudy cases across SZA.

Run from the repository root: ``python -m hypernet.wiggles.phase0_sky``.
"""
import os

import numpy as np
import pandas as pd

from hypernet import whn_l1a as wl
from hypernet.wiggles import REPO, WIGGLES_DIR

REQUEST = os.path.join(REPO, 'docs', 'wiggles_data_request.csv')
SZA_BINS = [0, 35, 50, 65, 80]


def build():
    req = pd.read_csv(REQUEST, dtype=str)
    rows = []
    for _, r in req.iterrows():
        path = wl.release2_path(r['site'], r['sequence_time'], r['file'])
        out = dict(row=r['row'], priority=r['priority'], site=r['site'],
                   sequence_time=r['sequence_time'], azimuth=r['azimuth'],
                   instrument=r['instrument'])
        if os.path.exists(path):
            d = wl.load_l2b(path)
            s = wl.sky_index(l2=d)
            saa, vaa = float(np.nanmean(d['saa'])), float(np.nanmean(d['vaa']))
            out.update(sza=float(np.nanmean(d['sza'])),
                       raa=abs(((vaa - saa) + 180) % 360 - 180),
                       ld_ed_750=s['ld_ed_750'], sky=wl.sky_class(s['ld_ed_750']))
        rows.append(out)
    return pd.DataFrame(rows)


def report(t):
    pd.set_option('display.width', 200)
    ok = t.dropna(subset=['ld_ed_750'])
    print('%d/%d sequences read; relative azimuths: %s'
          % (len(ok), len(t), np.unique(np.round(ok['raa'])).tolist()))
    print('ld_ed_750 quantiles (0, 10, 25, 50, 75, 90, 100 %):',
          np.round(ok['ld_ed_750'].quantile([0, .1, .25, .5, .75, .9, 1]).values, 4))
    order = pd.unique(t['row'])
    counts = (ok.groupby(['row', 'sky']).size().unstack(fill_value=0)
              .reindex(index=order, columns=['clear', 'cloudy', 'overcast'], fill_value=0))
    counts['n'] = counts.sum(axis=1)
    print('\nsky class counts per request row (broken needs L1A):')
    print(counts.to_string())
    ok = ok.assign(sza_bin=pd.cut(ok['sza'], SZA_BINS))
    print('\nclear floor (10th percentile) and median by SZA bin:')
    print(ok.groupby('sza_bin', observed=True)['ld_ed_750']
          .agg(n='size', p10=lambda x: x.quantile(.1), median='median').round(4).to_string())
    print('\nby site:')
    print(ok.groupby('site')['ld_ed_750']
          .agg(n='size', p10=lambda x: x.quantile(.1), median='median',
               frac_clear=lambda x: (x < wl.SKY_CLEAR).mean()).round(4).to_string())
    cov = (ok.assign(cloudy=ok['sky'] != 'clear')
           .groupby(['instrument', 'sza_bin', 'cloudy'], observed=True).size()
           .unstack(fill_value=0).rename(columns={False: 'clear', True: 'cloudy+'}))
    print('\nclear vs cloudy-or-overcast per instrument and SZA bin:')
    print(cov.to_string())
    gaps = []
    for inst, g in ok.groupby('instrument'):
        has = {k: set(g.loc[g['sky'] == k, 'sza_bin'].astype(str)) if k == 'clear'
               else set(g.loc[g['sky'] != 'clear', 'sza_bin'].astype(str))
               for k in ('clear', 'cloudy')}
        bins = set(g['sza_bin'].astype(str))
        missing = sorted(b for b in bins if b not in has['clear'] or b not in has['cloudy'])
        gaps.append(dict(instrument=inst, n=len(g), n_clear=int((g['sky'] == 'clear').sum()),
                         n_cloudy_or_overcast=int((g['sky'] != 'clear').sum()),
                         sza_bins=len(bins), bins_missing_a_class=', '.join(missing)))
    gaps = pd.DataFrame(gaps)
    print('\nper instrument: SZA bins without both clear and cloudy cases:')
    print(gaps.to_string(index=False))
    return counts, gaps


def main():
    t = build()
    t.to_csv(os.path.join(WIGGLES_DIR, 'phase0_sky_index.csv'), index=False,
             float_format='%.5g')
    report(t)
    print('\nwrote hypernet/wiggles/phase0_sky_index.csv')


if __name__ == '__main__':
    main()

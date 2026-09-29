"""Phase 0b, task 9: index the L1A/L1C/L2A delivery against the data request.

Before running, mirror the delivery (name per the delivery) with, e.g.::

    rclone copy AIOcean:data/Color/WATERHYPERNET/Wavelengths/<delivery> \\
        $OS_COLOR/WATERHYPERNET/Wavelengths/<delivery>

Then, from the repository root::

    python -m hypernet.wiggles.phase0b_index

For every product file under ``$OS_COLOR/WATERHYPERNET/Wavelengths`` (all
subfolders), groups the files into sequences (site, sequence time), matches
them to ``docs/wiggles_data_request.csv`` (site, ``sequence_time``, azimuth),
and for each delivered sequence reads the L1A attributes to check:

- the instrument (``system_id``) against the requested one;
- the calibration dates against the requested ``cal_dates_rad_irr``;
- Kevin's attribute bug: ``instrument_calibration_file_rad`` naming an IRR
  file;
- the SZA and the scan counts.

Per request row (e.g. "VEIT 122304 post-recal") it counts primaries
delivered, spares delivered and spares substituted for missing primaries.

Outputs: ``$OS_COLOR/hypernet/wiggles/phase0/request_index.parquet`` (one row
per requested or delivered sequence) and ``hypernet/wiggles/phase0_delivery_summary.csv``
(one row per request row, plus one for unrequested sequences).
"""
import argparse
import os

import numpy as np
import pandas as pd

from hypernet import whn_l1a as wl  # noqa: E402
from hypernet.wiggles import DATA_DIR, FIGDIR, OUT, REPO, WIGGLES_DIR  # noqa: F401

REQUEST = os.path.join(REPO, 'docs', 'wiggles_data_request.csv')
#: Water type per site (plan §5 / Setup #3 of wiggles_prompts.md).
WATER_TYPE = {'VEIT': 'clear', 'GAIT': 'clear', 'BEFR': 'dark', 'THFR': 'dark',
              'WRUK': 'dark', 'MAFR': 'turbid', 'LPAR': 'turbid', 'O1BE': 'turbid'}


def l1a_checks(row):
    """Attributes and counts from the L1A files of one delivered sequence."""
    out = {}
    if isinstance(row.get('L1A_IRR'), str):
        irr = wl.load_l1a_irr(row['L1A_IRR'])
        m = irr['meta']
        out.update(system_id=m['system_id'], sequence_id=m['sequence_id'],
                   cal_date_irr=m['instrument_calibration_date_irr'],
                   cal_file_irr=m['instrument_calibration_file_irr'],
                   sza_l1a=float(np.nanmean(irr['sza'])), n_scans_irr=irr['n_scans'],
                   n_irr_px=irr['wave'].size, time_l1a=int(np.nanmean(irr['time'])))
    if isinstance(row.get('L1A_RAD'), str):
        rad = wl.load_l1a_rad(row['L1A_RAD'])
        m = rad['meta']
        f_rad = m['instrument_calibration_file_rad'] or ''
        out.update(system_id_rad=m['system_id'], cal_date_rad=m['instrument_calibration_date_rad'],
                   cal_file_rad=f_rad, cal_file_rad_bug='_IRR_' in f_rad,
                   n_scans_lu=rad['Lu']['n_scans'], n_scans_ld=rad['Ld']['n_scans'],
                   n_rad_px=rad['wave'].size)
    return out


def build(root=None):
    request = pd.read_csv(REQUEST, dtype=str)
    products = wl.find_products(root)
    seqs = wl.sequence_table(products)
    idx = wl.match_request(seqs, request)
    checks = []
    for _, r in idx.iterrows():
        checks.append(l1a_checks(r) if r['status'] != 'missing' else {})
    idx = pd.concat([idx.reset_index(drop=True), pd.DataFrame(checks)], axis=1)
    # instrument and calibration consistency (requested rows only)
    idx['instrument_ok'] = np.where(idx['status'] == 'delivered',
                                    idx['instrument'] == idx.get('system_id'), None)
    cal_req = idx['cal_dates_rad_irr'].fillna('')
    cal_got = idx.get('cal_date_rad', pd.Series(index=idx.index, dtype=object)).fillna('') + '/' + \
        idx.get('cal_date_irr', pd.Series(index=idx.index, dtype=object)).fillna('')
    idx['cal_ok'] = np.where(idx['status'] == 'delivered', cal_req == cal_got, None)
    idx['cal_period'] = np.where(idx['status'] != 'missing', cal_got, cal_req)
    idx['water_type'] = idx['site_code'].map(WATER_TYPE)
    idx['month'] = idx['sequence_time'].str[4:6].astype(float)
    return idx, products


def summarise(idx):
    rows = []
    req = idx[idx['status'] != 'extra']
    order = pd.unique(pd.read_csv(REQUEST, dtype=str)['row'])
    for row in order:
        g = req[req['row'] == row]
        prim = g[g['priority'] == 'primary']
        spare = g[g['priority'] == 'spare']
        n_pd = int((prim['status'] == 'delivered').sum())
        n_sd = int((spare['status'] == 'delivered').sum())
        d = g[g['status'] == 'delivered']
        rows.append(dict(
            row=row, site=g['site'].iloc[0], instrument=g['instrument'].iloc[0],
            optional=g['optional'].iloc[0], n_primary=len(prim), n_primary_delivered=n_pd,
            n_spare=len(spare), n_spare_delivered=n_sd,
            n_spare_substituted=min(n_sd, len(prim) - n_pd),
            n_usable=n_pd + min(n_sd, len(prim) - n_pd),
            n_complete=int(d['complete'].eq(True).sum()),
            n_instrument_mismatch=int((d['instrument_ok'] == False).sum()),  # noqa: E712
            n_cal_mismatch=int((d['cal_ok'] == False).sum()),  # noqa: E712
            n_azimuth_mismatch=int((d['azimuth_ok'] == False).sum()),  # noqa: E712
            n_cal_attr_bug=int(d.get('cal_file_rad_bug', pd.Series(dtype=bool)).eq(True).sum()),
            n_sky_tagged=int(g['sky'].notna().sum())))
    ex = idx[idx['status'] == 'extra']
    rows.append(dict(row='(unrequested)', site=','.join(sorted(ex['site_code'].dropna().unique())),
                     instrument=','.join(sorted(ex.get('system_id', pd.Series(dtype=str))
                                                .dropna().unique())),
                     n_primary_delivered=len(ex), n_complete=int(ex['complete'].eq(True).sum()),
                     n_cal_attr_bug=int(ex.get('cal_file_rad_bug', pd.Series(dtype=bool))
                                        .eq(True).sum())))
    return pd.DataFrame(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--root', default=None, help='delivery root (default $OS_COLOR/...)')
    a = ap.parse_args(argv)
    os.makedirs(OUT, exist_ok=True)
    idx, products = build(a.root)
    idx.to_parquet(os.path.join(OUT, 'request_index.parquet'))
    summ = summarise(idx)
    summ.to_csv(os.path.join(WIGGLES_DIR, 'phase0_delivery_summary.csv'), index=False)
    pd.set_option('display.width', 220)
    print('%d product files; %d delivered sequences; status counts: %s'
          % (len(products), int((idx['status'] != 'missing').sum()),
             idx['status'].value_counts().to_dict()))
    print(summ.to_string(index=False))
    ex = idx[idx['status'] == 'extra']
    if len(ex):
        print('\nunrequested sequences:')
        print(ex[['site_code', 'sequence_time', 'azimuth_delivered', 'system_id', 'cal_period',
                  'complete', 'cal_file_rad_bug', 'sza_l1a']].to_string(index=False))
    print('\nwrote request_index.parquet and hypernet/wiggles/phase0_delivery_summary.csv')


if __name__ == '__main__':
    main()

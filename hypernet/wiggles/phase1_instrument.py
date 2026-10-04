"""Phase 1, task 7: the instrument data for the twin experiment.

Writes two small committed files read by ``hypernet.twin``:

- ``hypernet/data/hypstar_grids.npz``:
  - the real pixel grids: ``grid_E_122304`` (1536) and ``grid_L_122304``
    (1538) from the VEIT L1A sample (IRR cal 2024-11-12), and
    ``grid_L_122305`` from a Release 2 L2B of 122305 (it has 1536 pixels:
    pixel counts are instrument-specific).  No 122305 E
    grid exists without its L1A; ``twin.case_table`` uses 122304's;
  - the VEIT per-scan relative noise ``noise_{E,Ld,Lu}_scan``, the
    flattened scan scatter of Phase 0 Q9, divided by the mean.
- ``hypernet/data/release2_srf_models.json``: an E and an Ld SRFModel per
  instrument and calibration period, from the Release 2 template fits of
  Phase 0 task 8c (``release2_template.parquet``).
  - **Clear-sky sequences only** (sky index < 0.05; task 8e), because the
    fixed-veil Ld widths carry a Ring dependence on sky (task 8d).
  - Per 10 nm window: the median FWHM and dlam over the group's sequences,
    with error 1.25 MAD / sqrt(n).  Then ``SRFModel.from_lines``.
  - Keys: ``HYPSTAR_<id>@<IRR cal date>``, plus ``HYPSTAR_<id>`` for its
    latest calibration period.
  - E for every instrument but 122304 (2024-11) used a proxy E grid, and E
    from L2 reads ~0.1 nm wide (task 8c validation); see ``meta``.

Run from the repository root: ``python -m hypernet.wiggles.phase1_instrument``.
"""
import json
import os

import numpy as np
import pandas as pd

from hypernet import srf
from hypernet import whn_l1a as wl
from hypernet.wiggles import DATA_DIR, OUT, REPO, WIGGLES_DIR

REQUEST = os.path.join(REPO, 'docs', 'wiggles_data_request.csv')


def grids_and_noise():
    f = wl.sequence_files(site='VEIT', seq_time='20260604T0845')
    irr, rad = wl.load_l1a_irr(f['L1A_IRR']), wl.load_l1a_rad(f['L1A_RAD'])
    out = dict(grid_E_122304=irr['wave'], grid_L_122304=rad['wave'])
    for key, sc in (('E', irr['scans']), ('Ld', rad['Ld']['scans']), ('Lu', rad['Lu']['scans'])):
        mean, _, e_scan = srf.scan_errors(sc, flatten_px=20)
        out['noise_%s_scan' % key] = e_scan / mean
    req = pd.read_csv(REQUEST, dtype=str)
    r = req[req['row'] == 'VEIT 122305'].iloc[0]
    d = wl.load_l2b(wl.release2_path(r['site'], r['sequence_time'], r['file']))
    assert d['meta']['system_id'] == 'HYPSTAR_122305'
    out['grid_L_122305'] = d['wave']
    np.savez(os.path.join(DATA_DIR, 'hypstar_grids.npz'), **out)
    bE = (irr['wave'] > 400) & (irr['wave'] < 900)
    bL = (rad['wave'] > 400) & (rad['wave'] < 900)
    print('grids: E %d, L %d, L(122305) %d px; median per-scan noise 400-900 nm: E %.4f, '
          'Ld %.4f, Lu %.4f' % (out['grid_E_122304'].size, out['grid_L_122304'].size,
                                out['grid_L_122305'].size,
                                np.nanmedian(out['noise_E_scan'][bE]),
                                np.nanmedian(out['noise_Ld_scan'][bL]),
                                np.nanmedian(out['noise_Lu_scan'][bL])))
    for k in ('grid_E_122304', 'grid_L_122304', 'grid_L_122305'):
        g = out[k]
        print('  %-14s %4d px, %.2f-%.2f nm, spacing %.3f-%.3f nm'
              % (k, g.size, g[0], g[-1], np.diff(g).min(), np.diff(g).max()))
    return out


def release2_models():
    t = pd.read_parquet(os.path.join(OUT, 'release2_template.parquet'))
    sky = pd.read_csv(os.path.join(WIGGLES_DIR, 'phase0_sky_index.csv'),
                      dtype={'sequence_time': str})
    sky['site_code'] = sky['site'].str[:-2]
    t = t.merge(sky[['site_code', 'sequence_time', 'ld_ed_750']],
                on=['site_code', 'sequence_time'], how='left')
    t = t[t['ok'] & (t['ld_ed_750'] < 0.05)]
    models, rows = {}, []
    for (inst, cal, ch), g in t.groupby(['system_id', 'cal_period', 'channel']):
        n_seq = g['sequence_time'].nunique()
        agg = []
        for name, h in g.groupby('name'):
            fw, dl = h['fwhm'].values, h['dlam'].values
            mad = lambda v: 1.4826 * np.median(np.abs(v - np.median(v)))  # noqa: E731
            agg.append(dict(name=name, lam_air=h['lam_center'].iloc[0],
                            fwhm=np.median(fw), fwhm_err=max(1.25 * mad(fw), 0.005) /
                            np.sqrt(len(fw)), dmu=np.median(dl),
                            mu_err=max(1.25 * mad(dl), 0.005) / np.sqrt(len(dl)), n=len(fw)))
        a = pd.DataFrame(agg).assign(ok=True, use_for_srf=True, blend=False)
        a = a[a['n'] >= 3]
        irr_date = cal.split('/')[-1]
        m = srf.SRFModel.from_lines(ch, a.reset_index(drop=True), scale_cov=True, instrument=inst,
                                    meta=dict(cal_period=cal, n_sequences=int(n_seq),
                                              n_windows=int(len(a)), sky='clear (< 0.05)',
                                              method='median of Release 2 L2B template fits',
                                              e_grid='own' if (inst == 'HYPSTAR_122304' and
                                                               irr_date == '2024-11-12')
                                              else 'proxy (122304)'))
        key = '%s@%s' % (inst, irr_date)
        models.setdefault(key, {})['Ld' if ch == 'Ld' else ch] = m
        rows.append(dict(key=key, channel=ch, n_seq=n_seq, fwhm_450=m.fwhm(450),
                         fwhm_600=m.fwhm(600), fwhm_800=m.fwhm(800), chi2_nu=m.chi2_nu,
                         offset=m.offset))
    # plain keys -> the latest calibration period of each instrument
    latest = {}
    for key in models:
        inst, date = key.split('@')
        if inst not in latest or date > latest[inst]:
            latest[inst] = date
    doc = {'meta': dict(source='phase0c Release 2 template fits (task 8c), clear sky',
                        note='E from L2 irradiance reads ~0.1 nm wide (8c validation); '
                             'E grids other than 122304@2024-11-12 are a proxy',
                        script='hypernet/wiggles/phase1_instrument.py',
                        latest=latest),
           'models': {}}
    for key, chans in models.items():
        doc['models'][key] = {ch: m.to_dict() for ch, m in chans.items()}
    for inst, date in latest.items():
        doc['models'][inst] = doc['models']['%s@%s' % (inst, date)]
    with open(os.path.join(DATA_DIR, 'release2_srf_models.json'), 'w') as fh:
        json.dump(doc, fh, indent=1)
    s = pd.DataFrame(rows)
    pd.set_option('display.width', 200)
    print(s.to_string(index=False, float_format='%.3f'))
    return doc


def main():
    grids_and_noise()
    release2_models()
    print('wrote hypernet/data/hypstar_grids.npz and hypernet/data/release2_srf_models.json')


if __name__ == '__main__':
    main()

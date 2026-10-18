"""Phase 1, task 7: the instrument data for the twin experiment.

Writes two small committed files read by ``hypernet.twin``:

- ``hypernet/data/hypstar_grids.npz``:
  - the real pixel grids ``grid_E_<sn>`` and ``grid_L_<sn>`` of all seven
    instruments, the processor's grids from the lab wavelength polynomials
    (``hypernet.labcal.grid_for``; lab data prompts task 3).  Until then the
    file held only ``grid_E_122304`` / ``grid_L_122304`` (from the VEIT L1A)
    and ``grid_L_122305`` (from a Release 2 L2B), which these reproduce to
    <= 2.5e-5 nm (the L1A/L2B wavelengths are stored as float32); the old
    file is kept as ``hypstar_grids_pre_lab.npz``.  Pixel counts are
    instrument-specific (E 1535-1539, L 1536-1540);
  - the VEIT per-scan relative noise ``noise_{E,Ld,Lu}_scan``, the
    flattened scan scatter of Phase 0 Q9, divided by the mean (on the
    122304 grids).
- ``hypernet/data/release2_srf_models.json``: an E and an Ld SRFModel per
  instrument and calibration period, from the Release 2 template fits of
  Phase 0 task 8c (``release2_template.parquet``).
  - **Clear-sky sequences only** (sky index < 0.05; task 8e), because the
    fixed-veil Ld widths carry a Ring dependence on sky (task 8d).
  - Per 10 nm window: the median FWHM and dlam over the group's sequences,
    with error 1.25 MAD / sqrt(n).  Then ``SRFModel.from_lines``.
  - Keys: ``HYPSTAR_<id>@<IRR cal date>``, plus ``HYPSTAR_<id>`` for its
    latest calibration period.
  - E is fitted through each instrument's own E grid (lab data prompts task
    3; before, a 122304 proxy).  For the instruments whose lab E SRF is
    super-Gaussian (``labcal.srf_shape_for``), an extra channel ``E_sg``
    holds the same aggregation of the fixed-p super-Gaussian fits of
    ``hypernet.wiggles.lab_shape_fits`` (``shape``, ``shape_par`` set; the
    free-p diagnostic median in ``meta['shape_p_field']``).

Run from the repository root:
``python -m hypernet.wiggles.phase1_instrument [--grids-only | --models-only]``.
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

from hypernet import labcal, srf
from hypernet import whn_l1a as wl
from hypernet.wiggles import DATA_DIR, OUT, REPO, WIGGLES_DIR

REQUEST = os.path.join(REPO, 'docs', 'wiggles_data_request.csv')


def grids_and_noise():
    f = wl.sequence_files(site='VEIT', seq_time='20260604T0845')
    irr, rad = wl.load_l1a_irr(f['L1A_IRR']), wl.load_l1a_rad(f['L1A_RAD'])
    out = {}
    for sn in labcal.INSTRUMENTS:
        for ch in ('E', 'L'):
            out['grid_%s_%s' % (ch, sn)] = labcal.grid_for(sn, ch)
    # the 122304 grids are the VEIT L1A's (to the float32 storage)
    for k, w in (('grid_E_122304', irr['wave']), ('grid_L_122304', rad['wave'])):
        assert out[k].shape == w.shape and np.abs(out[k] - w).max() < 1e-4, k
    for key, sc in (('E', irr['scans']), ('Ld', rad['Ld']['scans']), ('Lu', rad['Lu']['scans'])):
        mean, _, e_scan = srf.scan_errors(sc, flatten_px=20)
        out['noise_%s_scan' % key] = e_scan / mean
    np.savez(os.path.join(DATA_DIR, 'hypstar_grids.npz'), **out)
    bE = (irr['wave'] > 400) & (irr['wave'] < 900)
    bL = (rad['wave'] > 400) & (rad['wave'] < 900)
    print('median per-scan noise 400-900 nm (VEIT): E %.4f, Ld %.4f, Lu %.4f'
          % (np.nanmedian(out['noise_E_scan'][bE]), np.nanmedian(out['noise_Ld_scan'][bL]),
             np.nanmedian(out['noise_Lu_scan'][bL])))
    for k in sorted(k for k in out if k.startswith('grid_')):
        g = out[k]
        print('  %-14s %4d px, %.2f-%.2f nm, spacing %.3f-%.3f nm'
              % (k, g.size, g[0], g[-1], np.diff(g).min(), np.diff(g).max()))
    return out


SG_TEMPLATE = os.path.join(labcal.lab_root(), 'lab_shape', 'release2_template_sg.parquet')


def _clear(t):
    sky = pd.read_csv(os.path.join(WIGGLES_DIR, 'phase0_sky_index.csv'),
                      dtype={'sequence_time': str})
    sky['site_code'] = sky['site'].str[:-2]
    t = t.merge(sky[['site_code', 'sequence_time', 'ld_ed_750']],
                on=['site_code', 'sequence_time'], how='left')
    return t[t['ok'] & (t['ld_ed_750'] < 0.05)]


def _aggregate(g):
    """Per 10 nm window: median FWHM and dlam over the sequences, error
    1.25 MAD / sqrt(n); windows with n >= 3."""
    mad = lambda v: 1.4826 * np.median(np.abs(v - np.median(v)))  # noqa: E731
    agg = []
    for name, h in g.groupby('name'):
        fw, dl = h['fwhm'].values, h['dlam'].values
        agg.append(dict(name=name, lam_air=h['lam_center'].iloc[0],
                        fwhm=np.median(fw), fwhm_err=max(1.25 * mad(fw), 0.005) /
                        np.sqrt(len(fw)), dmu=np.median(dl),
                        mu_err=max(1.25 * mad(dl), 0.005) / np.sqrt(len(dl)), n=len(fw)))
    a = pd.DataFrame(agg).assign(ok=True, use_for_srf=True, blend=False)
    return a[a['n'] >= 3].reset_index(drop=True)


def release2_models():
    t = _clear(pd.read_parquet(os.path.join(OUT, 'release2_template.parquet')))
    groups = [(k, g, {}) for k, g in t.groupby(['system_id', 'cal_period', 'channel'])]
    if os.path.exists(SG_TEMPLATE):
        sg = _clear(pd.read_parquet(SG_TEMPLATE))
        for (inst, cal), g in sg[sg['fit'] == 'fixed'].groupby(['system_id', 'cal_period']):
            free = sg[(sg['fit'] == 'free') & (sg['system_id'] == inst) &
                      (sg['cal_period'] == cal)]
            groups.append(((inst, cal, 'E_sg'), g, dict(
                shape='supergauss', shape_par=float(g['shape_par'].iloc[0]),
                shape_p_field=float(free['shape_par'].median()) if len(free) else None,
                shape_p_field_n=int(len(free)))))
    models, rows = {}, []
    for (inst, cal, ch), g, sh in groups:
        n_seq = g['sequence_time'].nunique()
        a = _aggregate(g)
        irr_date = cal.split('/')[-1]
        meta = dict(cal_period=cal, n_sequences=int(n_seq), n_windows=int(len(a)),
                    sky='clear (< 0.05)', method='median of Release 2 L2B template fits',
                    e_grid='own')
        kw = {}
        if sh:
            meta.update(shape_source='lab', shape_p_field=sh['shape_p_field'],
                        shape_p_field_n=sh['shape_p_field_n'],
                        method='median of Release 2 L2B super-Gaussian template fits, '
                               'p fixed at the lab value (hypernet.wiggles.lab_shape_fits)')
            kw = dict(shape='supergauss', shape_par=sh['shape_par'])
        m = srf.SRFModel.from_lines('E' if ch == 'E_sg' else ch, a, scale_cov=True,
                                    instrument=inst, meta=meta, **kw)
        key = '%s@%s' % (inst, irr_date)
        models.setdefault(key, {})[ch] = m
        rows.append(dict(key=key, channel=ch, n_seq=n_seq, fwhm_450=m.fwhm(450),
                         fwhm_600=m.fwhm(600), fwhm_800=m.fwhm(800), chi2_nu=m.chi2_nu,
                         offset=m.offset, shape_par=m.shape_par))
    # plain keys -> the latest calibration period of each instrument
    latest = {}
    for key in models:
        inst, date = key.split('@')
        if inst not in latest or date > latest[inst]:
            latest[inst] = date
    doc = {'meta': dict(source='phase0c Release 2 template fits (task 8c), clear sky',
                        note='E fitted through each instrument\'s own E grid '
                             '(hypernet.labcal.grid_for; lab data prompts task 3).  '
                             'E_sg: super-Gaussian E with p fixed at the lab value, '
                             'for the instruments labcal.srf_shape_for classes as '
                             'super-Gaussian; E stays the Gaussian fit',
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


def main(argv=None):
    ap = argparse.ArgumentParser(description='Phase 1 task 7: instrument data for the twin')
    g = ap.add_mutually_exclusive_group()
    g.add_argument('--grids-only', action='store_true')
    g.add_argument('--models-only', action='store_true')
    a = ap.parse_args(argv)
    if not a.models_only:
        grids_and_noise()
        print('wrote hypernet/data/hypstar_grids.npz')
    if not a.grids_only:
        release2_models()
        print('wrote hypernet/data/release2_srf_models.json')


if __name__ == '__main__':
    main()

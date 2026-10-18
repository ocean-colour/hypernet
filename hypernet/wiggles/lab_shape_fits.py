"""Field E SRFs with the lab super-Gaussian shape (lab data prompts, task 3;
Setup Q6-Q8).

For the instruments whose lab E SRF is super-Gaussian
(``labcal.srf_shape_for``: 120242, 122302, 122304), refit E with
``srf.fit_template_windows(shape='supergauss')``:

- ``fixed``: p held at the lab value of the epoch that describes the
  instrument at that calibration (``labcal.physical_epoch``); only the width,
  centroid and continuum are fitted.  This is the SRF-table E (Q6 (c)): its
  FWHM(lambda) belongs to the shape it is used with.  (The Gaussian template
  FWHM of a flat-topped SRF is not its FWHM.)
- ``free``: p fitted, started at the lab value.  A diagnostic only
  (``shape_p_field``): on L2 E, p trades off against the width and the
  E-grid resampling.

Inputs and conventions as Phase 0c (``phase0c_release2_srf``): veil 0, HSRS in
air, the 10 nm windows of ``srf.template_windows``.

- **VEIT L1A** (122304, IRR 2024-11-12): the scan-mean E with the flattened
  scan-scatter errors, as ``hypernet/data/veit_srf_model.json``.  Writes the
  committed ``hypernet/data/veit_srf_model_lab_shape.json`` (channel E: the
  fixed-p model; meta: the free-p median and the Gaussian model's FWHM for
  comparison).
- **Release 2 L2B** (``docs/wiggles_data_request.csv``): E through the own E
  grid (``via_grid`` = ``labcal.grid_for``).  Writes
  ``$OS_COLOR/hypernet/wiggles/lab_data/lab_shape/release2_template_sg.parquet``,
  which ``phase1_instrument.release2_models`` aggregates into the ``E_sg``
  models of ``release2_srf_models.json``.

Run from the repository root:
``python -m hypernet.wiggles.lab_shape_fits [--workers 12] [--limit N] [--veit-only]``
"""
import argparse
import os
import time

import numpy as np
import pandas as pd

from hypernet import labcal, srf, whn_srf
from hypernet import whn_l1a as wl
from hypernet.wiggles import DATA_DIR, REPO

REQUEST = os.path.join(REPO, 'docs', 'wiggles_data_request.csv')
OUT = os.path.join(labcal.lab_root(), 'lab_shape')
VEIT = dict(site='VEIT', seq_time='20260604T0845')
LAMS = (400, 450, 550, 700, 850)


def _fits(wave, spec, ref, p, err=None, via=None):
    """The fixed-p and free-p super-Gaussian template fits (one table)."""
    out = []
    for fit, fix in (('fixed', True), ('free', False)):
        t = srf.fit_template_windows(wave, spec, ref[0], ref[1], err=err, veil=False,
                                     via_grid=via, shape='supergauss', shape_par=p,
                                     fix_shape_par=fix)
        out.append(t.assign(fit=fit, p_lab=p))
    return pd.concat(out, ignore_index=True)


def veit(ref):
    f = wl.sequence_files(**VEIT)
    irr = wl.load_l1a_irr(f['L1A_IRR'])
    cal = irr['meta']['instrument_calibration_date_irr']
    sh = labcal.srf_shape_for('122304', 'E', cal)
    mean, e_mean, _ = srf.scan_errors(irr['scans'], flatten_px=whn_srf.FLATTEN_PX)
    t = _fits(irr['wave'], mean, ref, sh['shape_p'], err=e_mean)
    gauss = srf.load_srf_models(os.path.join(DATA_DIR, 'veit_srf_model.json'))['E']
    fx = t[t['fit'] == 'fixed'].reset_index(drop=True)
    fr = t[(t['fit'] == 'free') & t['ok']]
    m = srf.SRFModel.from_lines(
        'E', fx, scale_cov=True, instrument='HYPSTAR_122304', shape='supergauss',
        shape_par=sh['shape_p'],
        meta=dict(gauss.meta, method='template, super-Gaussian with p fixed at the lab value',
                  shape_source='lab', lab_epoch=sh['lab_epoch'],
                  shape_p_field=float(fr['shape_par'].median()),
                  shape_p_field_n=int(len(fr)),
                  gauss_fwhm={str(l): float(gauss.fwhm(l)) for l in LAMS},
                  script='hypernet/wiggles/lab_shape_fits.py'))
    srf.save_srf_models(os.path.join(DATA_DIR, 'veit_srf_model_lab_shape.json'), [m],
                        meta=dict(source='VEIT L1A (122304, IRR 2024-11-12)',
                                  note='E SRF with the lab super-Gaussian shape (Setup Q6 '
                                       '(c)); Ld/Lu are Gaussian: veit_srf_model.json',
                                  script='hypernet/wiggles/lab_shape_fits.py'))
    print('VEIT L1A E: lab p %.3f (%s); fixed-p FWHM %s vs Gaussian %s; free p median %.3f '
          '(%d windows), chi2 fixed/free/gauss %.2f/%.2f'
          % (sh['shape_p'], sh['lab_epoch'],
             '/'.join('%.3f' % m.fwhm(l) for l in LAMS),
             '/'.join('%.3f' % gauss.fwhm(l) for l in LAMS),
             m.meta['shape_p_field'], m.meta['shape_p_field_n'],
             fx['chi2_nu'].median(), fr['chi2_nu'].median()))
    return t.assign(source='VEIT L1A'), m


def _one(args):
    r, ref = args
    t0 = time.time()
    try:
        d = wl.load_l2b(wl.release2_path(r['site'], r['sequence_time'], r['file']))
        m = d['meta']
        sysid = m['system_id']
        irr = m['instrument_calibration_date_irr']
        sh = labcal.srf_shape_for(sysid, 'E', irr)
        if sh['shape'] != 'supergauss':
            return None
        t = _fits(d['wave'], d['irradiance'][:, 0], ref, sh['shape_p'],
                  via=labcal.grid_for(sysid, 'E'))
        return t.assign(source='Release 2 L2B', system_id=sysid, site_code=r['site'][:-2],
                        sequence_time=r['sequence_time'], lab_epoch=sh['lab_epoch'],
                        cal_period='%s/%s' % (m['instrument_calibration_date_rad'], irr),
                        seconds=round(time.time() - t0, 1))
    except Exception as ex:                      # noqa: BLE001 -- logged, not fatal
        return pd.DataFrame([dict(source='Release 2 L2B', site_code=r['site'][:-2],
                                  sequence_time=r['sequence_time'], ok=False,
                                  error='%s: %s' % (type(ex).__name__, ex))])


def release2(ref, workers, limit=None):
    from multiprocessing import Pool
    req = pd.read_csv(REQUEST, dtype=str)
    req = req[req['instrument'].str[-6:].isin(
        [sn for sn in labcal.INSTRUMENTS
         if labcal.srf_shape_for(sn, 'E', '2026-01-01')['shape'] == 'supergauss'])]
    if limit:
        req = req.groupby('row', group_keys=False).head(limit)
    out = []
    with Pool(workers) as pool:
        for t in pool.imap_unordered(_one, [(r, ref) for _, r in req.iterrows()]):
            if t is not None:
                out.append(t)
    t = pd.concat(out, ignore_index=True)
    os.makedirs(OUT, exist_ok=True)
    t.to_parquet(os.path.join(OUT, 'release2_template_sg.parquet'))
    return t


def summarise(t):
    ok = t[t['ok'].astype(bool)]
    s = ok.groupby(['system_id', 'cal_period', 'fit']).agg(
        n_seq=('sequence_time', 'nunique'), p_lab=('p_lab', 'first'),
        p_med=('shape_par', 'median'), fwhm_med=('fwhm', 'median'),
        chi2_med=('chi2_nu', 'median')).reset_index()
    return s


def main(argv=None):
    ap = argparse.ArgumentParser(description='lab-shape (super-Gaussian) E template fits')
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--veit-only', action='store_true')
    a = ap.parse_args(argv)
    ref = whn_srf.load_reference('air')
    t0 = time.time()
    veit(ref)
    if a.veit_only:
        return
    t = release2(ref, a.workers, a.limit)
    pd.set_option('display.width', 200)
    print(summarise(t).to_string(index=False, float_format='%.3f'))
    err = t['error'].notna().sum() if 'error' in t else 0
    print('%d sequences, %d errors, %.0f s' % (t['sequence_time'].nunique(), err,
                                              time.time() - t0))


if __name__ == '__main__':
    main()

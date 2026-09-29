"""SRF fits of one WATERHYPERNET sequence (E, Ld, Lu), by both methods.

Glue between :mod:`hypernet.whn_l1a` (readers) and :mod:`hypernet.srf`
(fitting), shared by the Phase 0 scripts (``hypernet/wiggles/phase0a_*.py``,
``hypernet/wiggles/phase0b_fit_all.py``).

- **empirical** -- :func:`hypernet.srf.fit_lines` over :data:`hypernet.srf.LINES`:
  a Gaussian per line.  The widths include each line's intrinsic profile.
- **template** -- :func:`hypernet.srf.fit_template_windows` against TSIS-1
  HSRS in air: sigma is the SRF width itself.

Conventions (Q&A Setup #1 and Build #4-#8 of
``claude_prompts/wiggles/wiggles_phase0_prompts.md``):

- E is the mean of the L1A_IRR scans; Ld the mean of the sky scans (vza >= 90)
  and Lu the mean of the water scans (vza < 90).
- Errors are the flattened scan scatter of the mean
  (``scan_errors(..., flatten_px=20)``).
- Ld is the reference L channel; Lu is carried as a check.
"""

import numpy as np
import pandas as pd

from hypernet import srf
from hypernet import whn_l1a as wl

#: Smoothing (pixels) for the flattened scan-scatter errors (Q&A Q9).
FLATTEN_PX = 20
#: Channels, in order; Ld is the reference L channel (Q&A Q11).
CHANNELS = ('E', 'Ld', 'Lu')
#: HSRS range loaded for the template fits, nm (air).
REF_RANGE = (375.0, 900.0)


def load_reference(frame='air'):
    """TSIS-1 HSRS over :data:`REF_RANGE`, in ``frame``."""
    from hypernet import refspec
    return refspec.load_hsrs(*REF_RANGE, frame=frame)


def channel_spectra(irr, rad):
    """``{channel: (wave, scans)}`` for E, Ld and Lu."""
    return {'E': (irr['wave'], irr['scans']),
            'Ld': (rad['wave'], rad['Ld']['scans']),
            'Lu': (rad['wave'], rad['Lu']['scans'])}


def _meta(irr, rad):
    m = irr['meta']
    return dict(site=m.get('site_id'), sequence_id=m.get('sequence_id'),
                system_id=m.get('system_id'),
                cal_date_irr=m.get('instrument_calibration_date_irr'),
                cal_date_rad=rad['meta'].get('instrument_calibration_date_rad'),
                sza=float(np.nanmean(irr['sza'])))


def fit_sequence(irr, rad, ref=None, template=True, empirical=True, veil=True,
                 per_scan=False, scale_cov=True):
    """Fit E, Ld and Lu of one sequence.

    Parameters
    ----------
    irr, rad : dict
        :func:`hypernet.whn_l1a.load_l1a_irr` / ``load_l1a_rad`` outputs.
    ref : tuple, optional
        ``(wave, flux)`` reference in air; default :func:`load_reference`.
        Only needed with ``template``.
    template, empirical : bool, optional
        Which methods to run.
    veil : bool, optional
        Fit the additive veil in the template fits.
    per_scan : bool, optional
        Also fit every single scan (empirical method only; slow-ish).
    scale_cov : bool, optional
        Inflate the FWHM(lambda) covariance by chi2_nu (see
        :func:`hypernet.srf.fit_fwhm_model`).

    Returns
    -------
    dict
        ``meta``; ``lines`` (empirical fits, one row per channel x line);
        ``template`` (one row per channel x window); ``scans`` (per-scan
        empirical fits, if ``per_scan``); ``models`` -- a list of
        :class:`hypernet.srf.SRFModel`, with ``meta['method']`` set to
        ``'empirical'`` or ``'template'``.
    """
    meta = _meta(irr, rad)
    chans = channel_spectra(irr, rad)
    if template and ref is None:
        ref = load_reference()
    out = dict(meta=meta, lines=None, template=None, scans=None, models=[])
    lines, temps, scans = [], [], []
    for ch, (w, sc) in chans.items():
        mean, e_mean, e_scan = srf.scan_errors(sc, flatten_px=FLATTEN_PX)
        mm = dict(meta, channel=ch)
        if empirical:
            df = srf.fit_lines(w, mean, err=e_mean)
            df.insert(0, 'channel', ch)
            lines.append(df)
            out['models'].append(srf.SRFModel.from_lines(
                ch, df, scale_cov=scale_cov, instrument=meta['system_id'],
                meta=dict(mm, method='empirical')))
            if per_scan:
                for k in range(sc.shape[1]):
                    d = srf.fit_lines(w, sc[:, k], err=e_scan)
                    d.insert(0, 'scan', k)
                    d.insert(0, 'channel', ch)
                    scans.append(d)
        if template:
            dt = srf.fit_template_windows(w, mean, ref[0], ref[1], err=e_mean, veil=veil)
            dt.insert(0, 'channel', ch)
            temps.append(dt)
            out['models'].append(srf.SRFModel.from_lines(
                ch, dt, scale_cov=scale_cov, instrument=meta['system_id'],
                meta=dict(mm, method='template', veil=veil)))
    for key, rows in (('lines', lines), ('template', temps), ('scans', scans)):
        if rows:
            df = pd.concat(rows, ignore_index=True)
            for k in ('site', 'sequence_id', 'system_id'):
                df.insert(0, k, meta[k])
            out[key] = df
    return out


def models_table(models):
    """Flatten a list of SRFModel into a DataFrame (one row per model)."""
    rows = []
    for m in models:
        d = dict(m.meta)
        d.update(channel=m.channel, instrument=m.instrument,
                 offset=m.offset, offset_err=m.offset_err, chi2_nu=m.chi2_nu,
                 npts=m.npts, lam_min=m.lam_min, lam_max=m.lam_max,
                 lam_ref=m.lam_ref, lam_scale=m.lam_scale, frame=m.frame)
        for i, c in enumerate(m.coeffs):
            d['c%d' % i] = c
            d['c%d_err' % i] = np.sqrt(m.cov[i, i])
        d['cov'] = m.cov.ravel().tolist()
        for lam in (400, 450, 500, 600, 700, 850):
            d['fwhm_%d' % lam] = m.fwhm(lam)
            d['fwhm_%d_err' % lam] = m.fwhm_err(lam)
        rows.append(d)
    return pd.DataFrame(rows)


def fit_files(files, **kw):
    """:func:`fit_sequence` from a ``{product: path}`` dict (L1A_IRR, L1A_RAD)."""
    irr = wl.load_l1a_irr(files['L1A_IRR'])
    rad = wl.load_l1a_rad(files['L1A_RAD'])
    return fit_sequence(irr, rad, **kw)

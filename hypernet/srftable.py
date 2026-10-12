"""The per-instrument SRF table (Phase 2 task 6).

``hypernet/data/srf_table.csv``, written by
``hypernet/wiggles/phase2_srf_table.py``: one row per instrument x channel
(``E``, ``Ld``) x calibration period x source, with the Gaussian FWHM(lambda)
quadratic and its covariance, the centroid offset, the wavelength-scale
correction (shift + stretch, calibration-period median, task 4b), the
E-vs-L width difference and the method gate G1 selects.

Lookups pick the row in force at a date: the latest ``valid_from`` (the
channel's own calibration date: IRR for E, RAD for Ld) on or before it, and
among sources the ``preferred`` one (L1A over Release 2 L2B).

**G0 is not closed** (Phase 0 task 14 waits on the L1A delivery).  Until it
is, the table follows the preliminary G0 findings, recorded here so that the
snippet can switch when G0 closes:

- :data:`TABLE_IS_DEFAULT` (G0(b), SRFs stable: preliminary pass): the table
  is the default SRF source; per-sequence :func:`hypernet.srfcal.calibrate_srf`
  is the alternative.  If G0(b) fails, flip it: self-calibration first, the
  table as fallback.
- :data:`APPLY_CENTROID_OFFSET` (G0(c), relative offset: preliminary pass):
  the SRF centroid offset is *not* applied as a recalibration.  The
  wavelength-scale shift + stretch of task 4b is applied instead
  (:func:`scale_for`, :func:`hypernet.srfcal.correct_wavelengths`).

Every row is ``provisional`` until the table is rebuilt from the delivery.
"""
import os

import numpy as np
import pandas as pd

from hypernet import srf as srfmod

TABLE_IS_DEFAULT = True
APPLY_CENTROID_OFFSET = False
#: Source preference: lower is preferred.
SOURCE_RANK = {'phase0b_l1a': 0, 'veit_l1a': 1, 'release2_l2b': 2}

_TABLE = {}


def table_path():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'srf_table.csv')


def load_srf_table(path=None):
    """The table as a DataFrame (cached per path); ``valid_from`` as Timestamps."""
    path = path or table_path()
    if path not in _TABLE:
        t = pd.read_csv(path, dtype={'cal_period': str, 'note': str})
        t['valid_from'] = pd.to_datetime(t['valid_from'])
        _TABLE[path] = t
    return _TABLE[path]


def _ts(date):
    """A Timestamp from a date, datetime, 'YYYY-MM-DD' or a sequence time
    'YYYYMMDDTHHMM'."""
    if isinstance(date, str) and len(date) >= 8 and date[:8].isdigit():
        return pd.Timestamp('%s-%s-%s' % (date[:4], date[4:6], date[6:8]))
    return pd.Timestamp(date)


def row_for(instrument, channel, date, table=None):
    """The table row in force for ``instrument`` / ``channel`` at ``date``.

    Raises KeyError if the instrument or channel is absent, or if ``date`` is
    before every calibration period in the table.
    """
    t = load_srf_table() if table is None else table
    g = t[(t['instrument'] == instrument) & (t['channel'] == channel)]
    if g.empty:
        raise KeyError('no SRF table rows for %s %s' % (instrument, channel))
    when = _ts(date)
    g = g[g['valid_from'] <= when]
    if g.empty:
        raise KeyError('%s %s: %s is before the first calibration period in the table'
                       % (instrument, channel, when.date()))
    g = g[g['valid_from'] == g['valid_from'].max()]
    g = g.assign(_rank=g['source'].map(SOURCE_RANK).fillna(99)).sort_values('_rank')
    return g.iloc[0].drop(labels='_rank')


def model_from_row(r):
    """The :class:`hypernet.srf.SRFModel` of a table row."""
    c = np.array([r['c0'], r['c1'], r['c2']], dtype=float)
    cov = np.array([[r['cov00'], r['cov01'], r['cov02']],
                    [r['cov01'], r['cov11'], r['cov12']],
                    [r['cov02'], r['cov12'], r['cov22']]], dtype=float)
    return srfmod.SRFModel(channel=r['channel'], coeffs=c, cov=cov, offset=r['offset'],
                           offset_err=r['offset_err'], lam_min=r['lam_min'],
                           lam_max=r['lam_max'], lam_ref=r['lam_ref'],
                           lam_scale=r['lam_scale'], chi2_nu=r['chi2_nu'],
                           instrument=r['instrument'], frame='air',
                           meta=dict(source=r['source'], cal_period=r['cal_period'],
                                     provisional=bool(r['provisional']),
                                     n_sequences=int(r['n_sequences'])))


def srf_for(instrument, channel, date, table=None):
    """The SRFModel of ``instrument`` / ``channel`` (``'E'`` or ``'Ld'``) in force
    at ``date``."""
    return model_from_row(row_for(instrument, channel, date, table))


def scale_for(instrument, channel, date, table=None):
    """The calibration-period wavelength-scale correction (task 4b), a dict
    for :func:`hypernet.srfcal.correct_wavelengths` and
    :func:`hypernet.srfcal.scale_or_fallback`: ``shift``, ``stretch``, their
    errors, ``n`` and ``ok``."""
    r = row_for(instrument, channel, date, table)
    return dict(shift=float(r['scale_shift']), stretch=float(r['scale_stretch']),
                shift_err=float(r['scale_shift_err']),
                stretch_err=float(r['scale_stretch_err']), n=int(r['scale_n']),
                ok=bool(np.isfinite(r['scale_shift']) and np.isfinite(r['scale_stretch'])),
                source='period median (%s)' % r['source'])


def method_for(instrument, date, table=None):
    """Gate G1's method for the instrument at ``date`` (from its E row)."""
    return str(row_for(instrument, 'E', date, table)['method'])

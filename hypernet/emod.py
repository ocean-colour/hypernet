"""High-resolution model irradiance (Emod) for the wiggles twin experiment.

Started for Phase 0 task 8d(i), where the telluric bands are fitted in the
template SRF fit; Phase 1 task 4 adds ozone, ``build_emod`` and the full
380-1000 nm cache.  Here: O2 and H2O line-by-line absorption from HITRAN via
HAPI (Kochanov et al. 2016; Gordon et al. 2022).

Approximations (good at the 3 nm HYPSTAR resolution, stated for Phase 1):

- One homogeneous layer at ``pressure_hpa`` and ``temperature_k`` (Voigt
  profiles, air-broadened).  The real column is pressure-weighted; lines
  from high, low-pressure air are narrower.
- Columns for one air mass: O2 = 0.2095 x the dry-air column for the
  surface pressure; H2O from the precipitable water vapour.
- The HITRAN line lists are fetched once with ``hapi.fetch`` (main
  isotopologue, 10000-26400 cm-1, i.e. 379-1000 nm) into
  ``$OS_COLOR/hypernet/wiggles/ref/hitran``.  Cross-sections are cached as
  npz next to them.
"""

import os

import numpy as np

from hypernet import srf

#: HITRAN molecule numbers (main isotopologue only).
MOLECULES = {'O2': 7, 'H2O': 1}
#: Wavenumber range of the cached line lists, cm-1.
NU_RANGE = (10000.0, 26400.0)
AVOGADRO = 6.02214076e23
#: Dry-air column for 1013.25 hPa, molecules cm-2 (p / (m_air g)).
AIR_COLUMN_1ATM = 101325.0 / (28.964e-3 / AVOGADRO * 9.80665) * 1e-4
O2_VMR = 0.2095


def hitran_dir():
    """``$OS_COLOR/hypernet/wiggles/ref/hitran`` (created if needed)."""
    d = os.path.join(os.getenv('OS_COLOR', '.'), 'hypernet', 'wiggles', 'ref', 'hitran')
    os.makedirs(d, exist_ok=True)
    return d


def _hapi():
    import hapi
    cwd = os.getcwd()
    try:
        os.chdir(hitran_dir())              # HAPI writes its cache in the cwd
        hapi.db_begin(hitran_dir())
    finally:
        os.chdir(cwd)
    return hapi


def fetch_lines(molecule):
    """Fetch the HITRAN lines of ``molecule`` once (needs network)."""
    if os.path.exists(os.path.join(hitran_dir(), molecule + '.data')):
        return
    hapi = _hapi()
    hapi.fetch(molecule, MOLECULES[molecule], 1, *NU_RANGE)


def cross_section(molecule, lam_min, lam_max, step_cm=0.02, pressure_hpa=1013.25,
                  temperature_k=296.0, frame='air'):
    """Absorption cross-section (cm2 / molecule) on a wavelength grid.

    Computed with ``hapi.absorptionCoefficient_Voigt`` on a uniform wavenumber
    grid (``step_cm``) and cached as npz.

    Returns
    -------
    tuple
        ``(lam, sigma)``: wavelength (nm, ascending; air or vacuum) and
        cross-section.
    """
    tag = '%s_%.1f_%.1f_%.3f_%.1f_%.1f.npz' % (molecule, lam_min, lam_max, step_cm,
                                                pressure_hpa, temperature_k)
    path = os.path.join(hitran_dir(), 'xsec_' + tag)
    if os.path.exists(path):
        d = np.load(path)
        nu, xs = d['nu'], d['xs']
    else:
        fetch_lines(molecule)
        hapi = _hapi()
        numin, numax = 1e7 / lam_max, 1e7 / lam_min
        nu, xs = hapi.absorptionCoefficient_Voigt(
            SourceTables=molecule, WavenumberRange=[numin, numax], WavenumberStep=step_cm,
            Environment={'p': pressure_hpa / 1013.25, 'T': temperature_k},
            Diluent={'air': 1.0}, HITRAN_units=True)
        nu, xs = np.asarray(nu), np.asarray(xs)
        np.savez(path, nu=nu, xs=xs)
    lam = 1e7 / nu[::-1]                     # vacuum nm, ascending
    if frame == 'air':
        lam = srf.vac_to_air(lam)
    return lam, xs[::-1]


def columns(pwv_mm=15.0, pressure_hpa=1013.25):
    """Vertical columns (molecules cm-2) for one air mass: O2 and H2O."""
    air = AIR_COLUMN_1ATM * pressure_hpa / 1013.25
    h2o = pwv_mm * 0.1 / 18.01528 * AVOGADRO        # mm -> g cm-2 -> molecules cm-2
    return {'O2': O2_VMR * air, 'H2O': h2o}


def optical_depth(lam, molecule, pwv_mm=15.0, pressure_hpa=1013.25, temperature_k=296.0,
                  step_cm=0.02, frame='air'):
    """Vertical (air mass 1) optical depth of ``molecule`` interpolated onto
    ``lam`` (nm)."""
    lam = np.asarray(lam, dtype=float)
    lo, hi = float(np.floor(lam.min())) - 1.0, float(np.ceil(lam.max())) + 1.0
    lg, xs = cross_section(molecule, lo, hi, step_cm, pressure_hpa, temperature_k, frame)
    col = columns(pwv_mm, pressure_hpa)[molecule]
    return np.interp(lam, lg, xs * col, left=0.0, right=0.0)


def gas_transmittance(lam, airmass, pwv_mm=15.0, pressure_hpa=1013.25, temperature_k=296.0,
                      molecules=('O2', 'H2O'), step_cm=0.02, frame='air'):
    """Slant-path transmittance exp(-airmass * sum tau) of O2 and H2O on ``lam``."""
    tau = sum(optical_depth(lam, m, pwv_mm, pressure_hpa, temperature_k, step_cm, frame)
              for m in molecules)
    return np.exp(-airmass * tau)

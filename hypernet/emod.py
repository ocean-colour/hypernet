"""High-resolution model irradiance (Emod) for the wiggles twin experiment.

``build_emod(sza, ...)`` returns, on a regular 0.01 nm air-wavelength grid
over 380-1000 nm:

- F0, the TSIS-1 HSRS v2 (``hypernet.refspec``), bin-averaged;
- the direct-beam transmittance, for the Kasten & Young (1989) air mass of
  the SZA;
- the diffuse transmittance, for an effective diffuse air mass
  (:data:`DIFFUSE_AIRMASS`, see below).

The transmittances include O2 and H2O line-by-line absorption from HITRAN via
HAPI (Kochanov et al. 2016; Gordon et al. 2022) and O3 from Serdyuchenko et
al. (2014).  Rayleigh, aerosol and the angular structure of the sky are
*not* included: they are smooth, and OSOAA supplies them on a 5 nm grid
(Phase 1 task 5); the twin experiment multiplies the two.

Started for Phase 0 task 8d(i), where the telluric bands are fitted in the
template SRF fit; completed in Phase 1 task 4.

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
- O2-O2 collision-induced absorption (477, 577, 630 nm) is not in the
  HITRAN line lists, so it is not modelled.  It is broad (several nm) and
  weak (a few per cent).
- O3 at a single temperature (default 293 K, Phase 1 Q&A Q7).  The Chappuis
  band that matters here (450-750 nm) changes by only a few per cent between
  223 and 293 K.
- **Diffuse air mass.**  Skylight is scattered mostly at altitude, above
  most of the water vapour, and reaches the surface from all directions.  Its
  gas path is approximated by the two-stream diffusivity factor, 1.66 (the
  mean secant of an isotropic radiance field), independent of the SZA.  The
  direct beam's path above the scattering height is ignored.  Case (vi) of
  the twin experiment (a wrong Emod) bounds the effect.
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
#: Effective air mass of the diffuse (sky) irradiance: the two-stream
#: diffusivity factor (see the module docstring).
DIFFUSE_AIRMASS = 1.66
#: One Dobson unit, molecules cm-2.
DU = 2.6867e16
#: Serdyuchenko et al. (2014) O3 cross-sections, IUP Bremen, 213-1100 nm at
#: 0.01 nm, vacuum wavelengths, 11 temperatures 193-293 K.
OZONE_URL = ('https://www.iup.uni-bremen.de/gruppen/molspec/downloads/'
             'serdyuchenkogorshelev5digits.dat')
#: SHA-256 of that file as downloaded on 2026-10-04.
OZONE_SHA256 = '4dfbf021b746512c192df5f0d43c54cee6ea3b4365bb490bcf6ed347f0ce7092'
OZONE_TEMPS = (293, 283, 273, 263, 253, 243, 233, 223, 213, 203, 193)


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


def _bin_average(x, y, x_out):
    """Mean of y(x) in bins centred on ``x_out`` (edges at the midpoints).

    ``x`` must be ascending and finer than ``x_out``.  Flux-conserving: it
    integrates the trapezoid rule on ``x``.
    """
    x_out = np.asarray(x_out, dtype=float)
    edges = np.empty(x_out.size + 1)
    edges[1:-1] = 0.5 * (x_out[1:] + x_out[:-1])
    edges[0] = x_out[0] - 0.5 * (x_out[1] - x_out[0])
    edges[-1] = x_out[-1] + 0.5 * (x_out[-1] - x_out[-2])
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))])
    c = np.interp(edges, x, cum)
    return np.diff(c) / np.diff(edges)


def gas_transmittance(lam, airmass, pwv_mm=15.0, pressure_hpa=1013.25, temperature_k=296.0,
                      molecules=('O2', 'H2O'), step_cm=0.02, frame='air', bin_average=True):
    """Slant-path O2 + H2O transmittance on ``lam`` (nm).

    exp(-airmass * sum tau) is computed on the fine HITRAN grid (``step_cm``)
    and, with ``bin_average``, averaged into bins centred on ``lam``, which
    is right for an output grid coarser than the lines (e.g. 0.01 nm).
    Otherwise tau is interpolated onto ``lam``.
    """
    lam = np.asarray(lam, dtype=float)
    if not bin_average:
        tau = sum(optical_depth(lam, m, pwv_mm, pressure_hpa, temperature_k, step_cm, frame)
                  for m in molecules)
        return np.exp(-airmass * tau)
    lo, hi = float(np.floor(lam.min())) - 1.0, float(np.ceil(lam.max())) + 1.0
    col = columns(pwv_mm, pressure_hpa)
    lg, tau = None, 0.0
    for m in molecules:
        lgm, xs = cross_section(m, lo, hi, step_cm, pressure_hpa, temperature_k, frame)
        if lg is None:
            lg = lgm
        tau = tau + np.interp(lg, lgm, xs * col[m])
    return _bin_average(lg, np.exp(-airmass * tau), lam)


# --- ozone ----------------------------------------------------------------------

def ozone_path():
    d = os.path.join(os.getenv('OS_COLOR', '.'), 'hypernet', 'wiggles', 'ref', 'ozone')
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, os.path.basename(OZONE_URL))


def fetch_ozone(verify=True):
    """Download the Serdyuchenko O3 file once; check its SHA-256."""
    import hashlib
    path = ozone_path()
    if not os.path.exists(path):
        import urllib.request
        urllib.request.urlretrieve(OZONE_URL, path + '.part')
        os.replace(path + '.part', path)
    if verify:
        h = hashlib.sha256(open(path, 'rb').read()).hexdigest()
        if h != OZONE_SHA256:
            raise ValueError('O3 cross-section checksum mismatch: %s' % path)
    return path


def ozone_cross_section(temperature_k=293.0, frame='air'):
    """O3 cross-section (cm2 / molecule) at the tabulated temperature nearest
    ``temperature_k``.

    Returns
    -------
    tuple
        ``(lam, sigma)``: nm (air or vacuum, ascending) and cm2 / molecule.
    """
    npz = ozone_path() + '.npz'
    if os.path.exists(npz):
        a = np.load(npz)['a']
    else:
        a = np.loadtxt(fetch_ozone(verify=False), skiprows=45)
        np.savez(npz, a=a)
    j = 1 + int(np.argmin(np.abs(np.array(OZONE_TEMPS) - temperature_k)))
    lam = a[:, 0]
    if frame == 'air':
        lam = srf.vac_to_air(lam)
    return lam, np.clip(a[:, j], 0.0, None)


def ozone_transmittance(lam, airmass, ozone_du=300.0, temperature_k=293.0, frame='air'):
    """Slant-path O3 transmittance on ``lam`` (nm) for ``ozone_du`` Dobson units."""
    lo, xs = ozone_cross_section(temperature_k, frame)
    return np.exp(-airmass * ozone_du * DU * np.interp(lam, lo, xs))


# --- Emod ---------------------------------------------------------------------------

def kasten_young_airmass(sza):
    """Relative optical air mass (Kasten & Young 1989) for SZA in degrees."""
    z = np.asarray(sza, dtype=float)
    return 1.0 / (np.cos(np.radians(z)) + 0.50572 * (96.07995 - z) ** -1.6364)


def f0_on_grid(lam, frame='air'):
    """TSIS-1 HSRS (W m-2 nm-1) bin-averaged onto ``lam`` (nm)."""
    from hypernet import refspec
    lam = np.asarray(lam, dtype=float)
    rw, rf = refspec.load_hsrs(lam.min() - 1.0, lam.max() + 1.0, step=1, frame=frame)
    return _bin_average(rw, rf, lam)


def emod_dir():
    d = os.path.join(os.getenv('OS_COLOR', '.'), 'hypernet', 'wiggles', 'ref', 'emod')
    os.makedirs(d, exist_ok=True)
    return d


def build_emod(sza, pwv_mm=15.0, ozone_du=300.0, pressure_hpa=1013.25, temperature_k=296.0,
               ozone_temperature_k=293.0, diffuse_airmass=DIFFUSE_AIRMASS, lam_min=380.0,
               lam_max=1000.0, step=0.01, cache=True, cache_dir=None):
    """High-resolution model irradiance pieces for one geometry and atmosphere.

    Parameters
    ----------
    sza : float
        Solar zenith angle, degrees; the direct air mass is
        :func:`kasten_young_airmass`.
    pwv_mm, ozone_du, pressure_hpa : float
        Precipitable water (mm), ozone column (DU), surface pressure (hPa).
    temperature_k, ozone_temperature_k : float
        Temperature of the O2/H2O layer, and of the O3 cross-section.
    diffuse_airmass : float
        Effective air mass of the diffuse irradiance.
    lam_min, lam_max, step : float
        The air-wavelength grid, nm (inclusive, regular).
    cache : bool
        Read and write ``<cache_dir>/*.npz``.
    cache_dir : str, optional
        Default ``$OS_COLOR/hypernet/wiggles/ref/emod`` (:func:`emod_dir`).

    Returns
    -------
    dict
        ``lam``, ``F0`` (W m-2 nm-1, 1 AU), ``T_direct``, ``T_diffuse`` and
        the separate ``T_gas_direct``, ``T_o3_direct``; scalars
        ``airmass_direct``, ``airmass_diffuse``, the inputs, and ``cached``.
        Ed = F0 cos(SZA) (T_direct t_dir + T_diffuse t_diff) once OSOAA's
        smooth direct and diffuse transmittances are multiplied in.
    """
    key = 'emod_sza%.2f_pwv%.1f_o3%.0f_p%.1f_t%.0f_to3%.0f_md%.2f_%.0f-%.0f_%.3f.npz' % (
        sza, pwv_mm, ozone_du, pressure_hpa, temperature_k, ozone_temperature_k,
        diffuse_airmass, lam_min, lam_max, step)
    path = os.path.join(cache_dir or emod_dir(), key)
    if cache and os.path.exists(path):
        d = dict(np.load(path))
        out = {k: (v.item() if v.ndim == 0 else v) for k, v in d.items()}
        out['cached'] = True
        return out
    n = int(round((lam_max - lam_min) / step)) + 1
    lam = lam_min + step * np.arange(n)
    m_dir = float(kasten_young_airmass(sza))
    out = dict(lam=lam, F0=f0_on_grid(lam), airmass_direct=m_dir,
               airmass_diffuse=float(diffuse_airmass), sza=float(sza), pwv_mm=float(pwv_mm),
               ozone_du=float(ozone_du), pressure_hpa=float(pressure_hpa))
    for tag, m in (('direct', m_dir), ('diffuse', float(diffuse_airmass))):
        tg = gas_transmittance(lam, m, pwv_mm, pressure_hpa, temperature_k)
        to = ozone_transmittance(lam, m, ozone_du, ozone_temperature_k)
        out['T_' + tag] = tg * to
        if tag == 'direct':
            out['T_gas_direct'], out['T_o3_direct'] = tg, to
    if cache:
        np.savez(path, **out)
    out['cached'] = False
    return out


# --- the Emod cache (Phase 2 task 5) --------------------------------------------------

#: The cache grid (Setup Q9): SZA in 5 deg bins, water vapour snapped to these
#: columns (mm), ozone and pressure fixed.  Gate G1 found the srf correction
#: insensitive to SZA +- 10 deg and water vapour x 0.5-3, so this is generous.
SZA_STEP = 5.0
SZA_RANGE = (0.0, 85.0)
PWV_GRID = (7.5, 15.0, 30.0)
OZONE_FIXED = 300.0
PRESSURE_FIXED = 1013.25
#: Diffuse fraction of Ed used to weight T_direct and T_diffuse in the model
#: spectrum (clear sky).  Only the line shapes matter to the interpolation
#: (the measured/model ratio absorbs the broad level), and they depend on it
#: weakly.
DIFFUSE_FRACTION = 0.15

_EMOD_MEM = {}
_CONV_MEM = {}
CONV_CACHE_SIZE = 64


def snap(sza, pwv_mm=15.0):
    """The cache-grid point for (sza, pwv): SZA to the nearest 5 deg (clipped to
    0-85), water vapour to the nearest of :data:`PWV_GRID` in log."""
    s = float(np.clip(SZA_STEP * np.round(float(sza) / SZA_STEP), *SZA_RANGE))
    g = np.asarray(PWV_GRID)
    w = float(g[np.argmin(np.abs(np.log(g) - np.log(max(float(pwv_mm), 1e-3))))])
    return s, w


def emod_for(sza, pwv_mm=15.0, ozone_du=None, pressure_hpa=None,
             diffuse_fraction=DIFFUSE_FRACTION, cache_dir=None, builder=None):
    """The model irradiance for a sequence, from the nearest cache-grid entry.

    Parameters
    ----------
    sza : float
        Solar zenith angle, deg (snapped to 5 deg bins).
    pwv_mm : float
        Water vapour column, mm (snapped to :data:`PWV_GRID`).
    ozone_du, pressure_hpa : float, optional
        Accepted for the interface; the grid fixes them at 300 DU and
        1013.25 hPa (Setup Q9), and the values used are returned.
    diffuse_fraction : float
        Weight of ``T_diffuse`` in the model spectrum.
    cache_dir : str, optional
        Where :func:`build_emod` stores its npz (default :func:`emod_dir`).
    builder : callable, optional
        Replaces :func:`build_emod` (tests).

    Returns
    -------
    dict
        ``lam``, ``E`` = F0 [(1 - f_d) T_direct + f_d T_diffuse] (W m-2 nm-1,
        a model *shape*: no cos(SZA) or aerosol), the ``build_emod`` arrays,
        and ``sza_bin``, ``pwv_bin``, ``ozone_du``, ``pressure_hpa``,
        ``diffuse_fraction``, ``key``.  The same object is returned for every
        call in the same bin (in-memory cache); a miss builds it (or reads
        ``build_emod``'s disk cache under ``ref/emod/``).  Usable directly as
        ``emod=`` in :func:`hypernet.edinterp.interpolate_ed_to_l`.
    """
    s, w = snap(sza, pwv_mm)
    key = ('emod', s, w, OZONE_FIXED, PRESSURE_FIXED, float(diffuse_fraction),
           cache_dir or '')
    if key in _EMOD_MEM:
        return _EMOD_MEM[key]
    build = builder or build_emod
    e = dict(build(s, pwv_mm=w, ozone_du=OZONE_FIXED, pressure_hpa=PRESSURE_FIXED,
                   cache_dir=cache_dir))
    fd = float(diffuse_fraction)
    e['E'] = e['F0'] * ((1.0 - fd) * e['T_direct'] + fd * e['T_diffuse'])
    e.update(sza_bin=s, pwv_bin=w, ozone_du=OZONE_FIXED, pressure_hpa=PRESSURE_FIXED,
             diffuse_fraction=fd, key=key)
    _EMOD_MEM[key] = e
    return e


def emod_on_grid(entry, grid, srf):
    """``entry['E']`` through ``srf`` at ``grid`` (Emod_X of plan §3), cached per
    (Emod entry, grid, SRF).  ``srf``: an SRFModel, constant or per-pixel
    array (anything :func:`hypernet.srf.srf_key` can key; a plain callable is
    computed but not cached)."""
    from hypernet import srf as srfmod
    grid = np.asarray(grid, dtype=float)
    k = srfmod.srf_key(srf)
    key = None if k is None else (entry.get('key'), hash(grid.tobytes()), k)
    if key is not None and key in _CONV_MEM:
        return _CONV_MEM[key]
    out = srfmod.convolve_to_grid(entry['lam'], entry['E'], grid, srf)
    out.setflags(write=False)
    if key is not None:
        _CONV_MEM[key] = out
        while len(_CONV_MEM) > CONV_CACHE_SIZE:
            _CONV_MEM.pop(next(iter(_CONV_MEM)))
    return out


def clear_memory_cache():
    """Drop the in-memory Emod entries and convolutions (the disk cache stays)."""
    _EMOD_MEM.clear()
    _CONV_MEM.clear()

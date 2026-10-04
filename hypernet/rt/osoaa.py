"""Thin wrapper around OSOAA (Chami et al. 2015; JXP's fork).

Merged from the ``OSOAASimulation`` class in the fork's
``tests/test_water_Ed.py`` and its tutorial notebooks.  Nothing is imported
from the fork.  The executable is ``$OSOAA_ROOT/exe/OSOAA_MAIN.exe``.  One call
is one wavelength.

OSOAA facts used here (Phase 1 task 2; the fork's ``src/OSOAA_MAIN.F``):

- ``-OSOAA.View.Level``: 1 = TOA, 2 = sea bed, 3 = sea surface 0+,
  4 = sea surface 0-, 5 = user ``-OSOAA.View.Z``.  (The fork's Sphinx docs,
  which give +-1, are wrong.)
- ``Advanced_outputs/Flux.txt``: direct/diffuse/total down- and upwelling
  flux at every level, for a TOA solar irradiance of pi.
- ``-OSOAA.ResFile.Adv.Up`` / ``.Adv.Down``: upward / downward radiance at
  every level and viewing zenith angle, for the azimuth ``-OSOAA.View.Phi``.
  So one run gives Ld at 0+ and Lu at 0- together.
- ``Standard_outputs/<vsVZA>``: radiance vs VZA at the one level set by
  ``View.Level``.
- Radiances are normalised, I = pi L / E_sun (TOA).

**Aerosol reference wavelength.**  The fork's wrappers set ``AER.Waref`` to the
simulated wavelength, which makes AOT the same at every wavelength and loses
the aerosol's spectral shape.  Here ``AER.Waref`` is fixed (default 0.55 um),
and OSOAA scales AOT(lambda) with its Mie model.
"""

import os
import shutil
import subprocess
import tempfile
import time

import numpy as np

#: Fallback location of JXP's fork when $OSOAA_ROOT is unset (Phase 1 Q&A Q3).
DEFAULT_ROOT = '/Users/xavier/Oceanography/python/RadiativeTransferCode-OSOAA'
#: View.Level codes.
LEVELS = dict(toa=1, seabed=2, surface_plus=3, surface_minus=4, user=5)
ADV_UP, ADV_DOWN, VS_VZA = 'LUM_Advanced_Up.txt', 'LUM_Advanced_Down.txt', 'LUM_vsVZA.txt'


def osoaa_root(path=None):
    """The OSOAA root: ``path``, else ``$OSOAA_ROOT``, else :data:`DEFAULT_ROOT`.

    Raises
    ------
    FileNotFoundError
        If ``exe/OSOAA_MAIN.exe`` is not there.
    """
    root = path or os.getenv('OSOAA_ROOT') or DEFAULT_ROOT
    exe = os.path.join(root, 'exe', 'OSOAA_MAIN.exe')
    if not os.path.exists(exe):
        raise FileNotFoundError('OSOAA executable not found: %s' % exe)
    return os.path.abspath(root)


def exe_path(root=None):
    return os.path.join(osoaa_root(root), 'exe', 'OSOAA_MAIN.exe')


def default_params(wavelength_nm, sza, work_dir, mie_dir=None, phi=90.0, level=3,
                   aot550=0.1, aer_waref_um=0.55, chl=1.0, csed=0.0, ys440=0.0,
                   det440=0.0, wind=5.0, sea_depth=100.0, pressure=1013.0):
    """OSOAA keyword dict for one wavelength.

    Parameters
    ----------
    wavelength_nm, sza : float
        Wavelength (nm) and solar zenith angle (deg).
    work_dir : str
        ``-OSOAA.ResRoot`` (outputs go to its Standard_ and Advanced_outputs).
    mie_dir : str, optional
        Directory for the Mie and surface-matrix databases (``MIE_AER``,
        ``MIE_HYD``, ``SURF``).  Reuse it across runs: OSOAA caches per
        wavelength.  Default: inside ``work_dir``.
    phi : float
        Relative azimuth (deg) of the radiance outputs.  HYPSTAR uses 90
        (the VEIT sample and every requested sequence).
    level : int
        ``View.Level`` for the vsVZA file (default 3 = 0+).
    aot550, aer_waref_um : float
        AOT at the fixed reference wavelength (default 0.55 um; see the module
        docstring).
    chl, csed, ys440, det440 : float
        Chlorophyll (mg m-3), sediment (mg L-1), and yellow-substance and
        detritus absorption at 440 nm (m-1).
    wind, sea_depth, pressure : float
        Wind (m s-1), sea depth (m), surface pressure (hPa).
    """
    mie_dir = mie_dir or work_dir
    d = {k: os.path.join(mie_dir, k) for k in ('MIE_AER', 'MIE_HYD', 'SURF')}
    for p in d.values():
        os.makedirs(p, exist_ok=True)
    return {
        'OSOAA.ResRoot': work_dir,
        'OSOAA.Wa': wavelength_nm / 1000.0,
        'ANG.Thetas': sza,
        'OSOAA.View.Phi': phi,
        'OSOAA.View.Level': level,
        'OSOAA.View.Z': 0.0,
        'AP.Pressure': pressure, 'AP.HR': 8.0, 'AP.HA': 2.0,
        'AER.DirMie': d['MIE_AER'], 'AER.Waref': aer_waref_um, 'AER.AOTref': aot550,
        'AER.Model': 0, 'AER.MMD.MRwa': 1.45, 'AER.MMD.MIwa': -0.001,
        'AER.MMD.SDtype': 1, 'AER.MMD.LNDradius': 0.10, 'AER.MMD.LNDvar': 0.46,
        'SEA.Depth': sea_depth,
        'HYD.DirMie': d['MIE_HYD'], 'HYD.Model': 1,
        'PHYTO.Chl': chl, 'PHYTO.ProfilType': 1,
        'PHYTO.JD.slope': 4.0, 'PHYTO.JD.rmin': 0.01, 'PHYTO.JD.rmax': 200.0,
        'PHYTO.JD.MRwa': 1.05, 'PHYTO.JD.MIwa': 0.0, 'PHYTO.JD.rate': 1.0,
        'SED.Csed': csed, 'YS.Abs440': ys440, 'DET.Abs440': det440,
        'SEA.Dir': d['SURF'], 'SEA.Ind': 1.34, 'SEA.Wind': wind, 'SEA.SurfAlb': 0.0,
        'SEA.BotType': 1, 'SEA.BotAlb': 0.30,
        'OSOAA.ResFile.vsVZA': VS_VZA,
        'OSOAA.ResFile.Adv.Up': ADV_UP,
        'OSOAA.ResFile.Adv.Down': ADV_DOWN,
    }


def build_command(params, root=None):
    """``[exe, -KEY, value, ...]`` for :func:`subprocess.run`."""
    cmd = [exe_path(root)]
    for k, v in params.items():
        cmd += ['-' + k, str(v)]
    return cmd


def run(params, work_dir=None, root=None, timeout=3600):
    """Run OSOAA once (one wavelength) and parse its outputs.

    Parameters
    ----------
    params : dict
        From :func:`default_params`; ``OSOAA.ResRoot`` is set to ``work_dir``
        if given.
    work_dir : str, optional
        Output directory; default ``params['OSOAA.ResRoot']``.
    root : str, optional
        OSOAA root (see :func:`osoaa_root`).
    timeout : float
        Seconds.

    Returns
    -------
    dict
        ``flux`` (:func:`parse_flux`), ``up`` and ``down``
        (:func:`parse_lum_advanced`), ``vsvza`` (:func:`parse_lum_vsvza`),
        ``seconds`` and ``work_dir``.

    Raises
    ------
    RuntimeError
        If OSOAA exits non-zero (stdout/stderr tail in the message).
    """
    root = osoaa_root(root)
    if work_dir is not None:
        params = dict(params, **{'OSOAA.ResRoot': work_dir})
    work_dir = params['OSOAA.ResRoot']
    os.makedirs(work_dir, exist_ok=True)
    t0 = time.time()
    res = subprocess.run(build_command(params, root), cwd=work_dir, capture_output=True,
                         text=True, timeout=timeout, env=dict(os.environ, OSOAA_ROOT=root))
    dt = time.time() - t0
    if res.returncode != 0:
        raise RuntimeError('OSOAA failed (%d): %s %s' % (res.returncode, res.stdout[-800:],
                                                         res.stderr[-800:]))
    adv = os.path.join(work_dir, 'Advanced_outputs')
    std = os.path.join(work_dir, 'Standard_outputs')
    out = dict(seconds=dt, work_dir=work_dir,
               flux=parse_flux(os.path.join(adv, 'Flux.txt')))
    for key, fn, folder in (('up', params.get('OSOAA.ResFile.Adv.Up'), adv),
                            ('down', params.get('OSOAA.ResFile.Adv.Down'), adv),
                            ('vsvza', params.get('OSOAA.ResFile.vsVZA'), std)):
        p = os.path.join(folder, fn) if fn else None
        if p and os.path.exists(p):
            out[key] = parse_lum_vsvza(p) if key == 'vsvza' else parse_lum_advanced(p)
    return out


# --- parsers -------------------------------------------------------------------

def _numeric_rows(lines, ncol):
    """Rows of ``lines`` that parse as at least ``ncol`` floats."""
    rows = []
    for ln in lines:
        parts = ln.split()
        if len(parts) < ncol:
            continue
        try:
            rows.append([float(x) for x in parts[:ncol]])
        except ValueError:
            continue
    return np.array(rows, dtype=float).reshape(-1, ncol)


def parse_flux(path):
    """``Advanced_outputs/Flux.txt``: fluxes at every level (TOA irradiance = pi).

    Returns
    -------
    dict
        Arrays ``level`` (int), ``z`` (m; > 0 atmosphere, <= 0 sea),
        ``Ed_direct``, ``Ed_diffuse``, ``Ed_total``, ``Eu_direct``,
        ``Eu_diffuse``, ``Eu_total``, ``ratio`` (Eu/Ed).
    """
    with open(path, encoding='latin-1') as fh:
        a = _numeric_rows(fh.readlines(), 9)
    keys = ['level', 'z', 'Ed_direct', 'Ed_diffuse', 'Ed_total', 'Eu_direct', 'Eu_diffuse',
            'Eu_total', 'ratio']
    out = {k: a[:, i] for i, k in enumerate(keys)}
    out['level'] = out['level'].astype(int)
    return out


def _header_value(lines, key):
    for ln in lines:
        if key in ln:
            try:
                return float(ln.split()[-1])
            except ValueError:
                return None
    return None


def parse_lum_advanced(path):
    """``Advanced_outputs/LUM_Advanced_{Up,Down}.txt``: radiance at every level
    and viewing zenith angle, for one relative azimuth.

    VZA < 0 is the azimuth ``phi_neg`` (phi + 180), VZA > 0 the azimuth
    ``phi_pos`` (the requested ``OSOAA.View.Phi``).  I, Q, U and LPOL are
    normalised to pi L / E_sun.

    Returns
    -------
    dict
        ``direction`` ('up' / 'down'), ``level_0plus``, ``level_0minus``,
        ``phi_neg``, ``phi_pos``, and arrays ``level`` (int), ``z``, ``vza``,
        ``sca``, ``I``, ``Q``, ``U``, ``pol_ang``, ``pol_rate``, ``lpol``.
    """
    with open(path, encoding='latin-1') as fh:
        lines = fh.readlines()
    head = ''.join(lines[:40]).upper()
    direction = 'down' if 'DOWNWARD' in head else ('up' if 'UPWARD' in head else None)
    out = dict(direction=direction)
    for ln in lines[:60]:
        if '0+ (over surface) is level' in ln:
            out['level_0plus'] = int(ln.split()[-1])
        if '0- (under surface) is level' in ln:
            out['level_0minus'] = int(ln.split()[-1])
        if 'for VZA < 0' in ln:
            out['phi_neg'] = float(ln.split()[-1])
        if 'for VZA > 0' in ln:
            out['phi_pos'] = float(ln.split()[-1])
    a = _numeric_rows(lines, 10)
    keys = ['level', 'z', 'vza', 'sca', 'I', 'Q', 'U', 'pol_ang', 'pol_rate', 'lpol']
    out.update({k: a[:, i] for i, k in enumerate(keys)})
    out['level'] = out['level'].astype(int)
    return out


def parse_lum_advanced_down(path):
    """:func:`parse_lum_advanced` for the downward file (sky radiance, Ld)."""
    d = parse_lum_advanced(path)
    if d['direction'] not in (None, 'down'):
        raise ValueError('%s is not a downward-radiance file' % path)
    return d


def radiance_at(adv, level, vza):
    """I at ``level``, linearly interpolated in VZA (signed; see
    :func:`parse_lum_advanced`).  Accepts ``'0+'`` / ``'0-'`` as levels."""
    if level == '0+':
        level = adv['level_0plus']
    elif level == '0-':
        level = adv['level_0minus']
    m = adv['level'] == level
    v, i = adv['vza'][m], adv['I'][m]
    o = np.argsort(v)
    return np.interp(vza, v[o], i[o])


def parse_lum_vsvza(path):
    """``Standard_outputs/LUM_vsVZA.txt``: radiance vs VZA at the ``View.Level``.

    Returns
    -------
    dict
        ``level_label`` (the header line naming the level), and arrays ``vza``,
        ``sca``, ``I`` (pi L / E_sun), ``refl`` (pi L / Ed(z)), ``pol_rate``,
        ``lpol``, ``refl_pol``.
    """
    with open(path, encoding='latin-1') as fh:
        lines = fh.readlines()
    label = next((ln.strip() for ln in lines[:40] if 'level' in ln.lower() and
                  ('surface' in ln.lower() or 'TOA' in ln or 'bottom' in ln.lower()
                   or 'Altitude' in ln or 'Depth' in ln)), None)
    a = _numeric_rows(lines, 7)
    keys = ['vza', 'sca', 'I', 'refl', 'pol_rate', 'lpol', 'refl_pol']
    out = {k: a[:, i] for i, k in enumerate(keys)}
    out['level_label'] = label
    return out

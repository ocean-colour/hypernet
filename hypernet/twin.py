"""Twin experiment (Phase 1): high-resolution scenes, the rho_w library, and
the true rho_w.

A scene combines the 0.01 nm model irradiance of :mod:`hypernet.emod` (F0 x
the gas transmittances) with the smooth OSOAA fields of
``hypernet/wiggles/phase1_osoaa_fields.py`` (5 nm, no gas absorption):

- Ed = F0 (ed_dir T_direct + ed_dif T_diffuse): separate direct and diffuse
  air masses;
- Ld = F0 ld T_diffuse: skylight takes the diffuse gas path;
- Lw = rho_el Ed / pi + Lw_fl + Lw_raman.  The elastic part carries Ed's line
  structure (it is Ed reflected by the water body).  The controls are real
  spectral features that a correction must **not** remove:
  - **chlorophyll fluorescence**: a Gaussian at 683 nm, 25 nm FWHM, excited
    by broadband light, so it has no Fraunhofer structure;
  - **water Raman**: Ed redistributed by the 3357 cm-1 O-H stretch, broadened
    by the Raman band (~200 cm-1 FWHM).  Its efficiency follows a simple
    clear-water model, rho_R(lambda) proportional to
    b_R(lambda_exc) / (2 a_w(lambda_exc) + a_w(lambda)), with
    b_R ~ lambda_exc^-5.5 (Bartlett et al. 1998) and a_w from OSOAA's
    total absorption a_t = a_w + a_ph(Chl) + a_ys (OSOAA's pure-water and
    Bricaud tables).  It is normalised so that pure water gives 5e-4 at
    550 nm; absorption at the emission wavelength suppresses it in the red
    and NIR, and pigments and CDOM reduce it in greener water.  It carries
    the solar lines, shifted.
  - Both controls are emitted in the water, below the sensor's atmospheric
    path, so rho_w = pi Lw / Ed genuinely *peaks* in the gas bands (O2-A,
    O2-B, H2O), where Ed is absorbed but they are not.  That is the real
    "structure at low-transmittance wavelengths" that a correction must not
    erase.
- Lu = Lw + rho_eff Ld: the water-view radiance HYPSTAR measures, with the
  OSOAA effective surface reflectance rho_eff = (Lu(0+) - Lw) / Ld.

The smooth fields are interpolated from 5 nm onto 0.01 nm with a cubic spline
in ln(value).

The **true** rho_w (:func:`rhow_true`) is what a perfect instrument with the
L SRF would report: Lw and Ed both observed with the L SRF on the L grid,
then pi Lw / Ed.

Units: F0 in W m-2 nm-1, so Ed in W m-2 nm-1 and radiances in W m-2 nm-1 sr-1.
"""

import numpy as np
from scipy.interpolate import CubicSpline

from hypernet import srf as srfmod

#: Chlorophyll fluorescence control: centre and FWHM, nm.
FL_CENTRE, FL_FWHM = 683.0, 25.0
#: Water Raman shift (cm-1) and band FWHM (cm-1).
RAMAN_SHIFT, RAMAN_FWHM = 3357.0, 200.0


# --- observation -----------------------------------------------------------------

def convolve_to_grid(lam_hr, spec_hr, grid, fwhm, nsig=5.0):
    """Observe ``spec_hr`` (on ``lam_hr``, nm) with a Gaussian SRF on ``grid``.

    Parameters
    ----------
    lam_hr, spec_hr : array
        High-resolution spectrum (ascending ``lam_hr``; uniform or not).
    grid : array
        Pixel centres, nm.
    fwhm : float, array or callable
        SRF FWHM (nm): a constant, one value per pixel, or a function of
        wavelength (e.g. ``SRFModel.fwhm``).
    nsig : float
        Kernel half-width in sigma.

    Returns
    -------
    numpy.ndarray
        The SRF-weighted mean of ``spec_hr`` at each pixel (weights include
        the sample spacing, so a constant maps to the same constant).
    """
    lam_hr = np.asarray(lam_hr, dtype=float)
    spec_hr = np.asarray(spec_hr, dtype=float)
    grid = np.asarray(grid, dtype=float)
    f = fwhm(grid) if callable(fwhm) else np.broadcast_to(np.asarray(fwhm, float), grid.shape)
    sig = np.asarray(f, dtype=float) / srfmod.FWHM_PER_SIGMA
    dw = np.gradient(lam_hr)
    out = np.empty(grid.size)
    lo = np.searchsorted(lam_hr, grid - nsig * sig)
    hi = np.searchsorted(lam_hr, grid + nsig * sig)
    for k in range(grid.size):
        sl = slice(lo[k], hi[k])
        w = np.exp(-0.5 * ((lam_hr[sl] - grid[k]) / sig[k]) ** 2) * dw[sl]
        out[k] = np.dot(w, spec_hr[sl]) / w.sum() if w.size else np.nan
    return out


# --- fields and the rho_w library ---------------------------------------------------

def load_fields(path=None):
    """The OSOAA fields npz as a dict (default ``$OS_COLOR/.../osoaa_fields.npz``)."""
    import os
    path = path or os.path.join(os.getenv('OS_COLOR', '.'), 'hypernet', 'wiggles', 'phase1',
                                'osoaa_fields.npz')
    d = np.load(path, allow_pickle=False)
    return {k: d[k] for k in d.files}


def _case_index(fields, case):
    if isinstance(case, (int, np.integer)):
        return int(case)
    names = list(np.asarray(fields['name']).astype(str))
    return names.index(case)


def _spline_log(lam5, y5, lam):
    """Cubic spline in ln(y) from the 5 nm grid onto ``lam`` (y > 0)."""
    return np.exp(CubicSpline(lam5, np.log(np.clip(y5, 1e-30, None)))(lam))


def fluorescence_amplitude(chl):
    """Peak fluorescence reflectance (pi Lw_fl / Ed) for a chlorophyll
    concentration (mg m-3): 3e-4 chl^0.7, i.e. ~6e-5 (Chl 0.1), 3e-4 (Chl 1)
    and 1.5e-3 (Chl 10).  These are typical in situ magnitudes, scaled
    sub-linearly (fluorescence quantum yield falls with Chl)."""
    return 3e-4 * np.asarray(chl, dtype=float) ** 0.7


#: Raman reflectance at 550 nm, pi Lw_R / Ed(lambda_exc) ~ 5e-4: the clear-
#: water magnitude in the green.
RAMAN_EFFICIENCY = 5e-4


def _osoaa_table(name):
    import os
    from hypernet.rt import osoaa
    try:
        root = osoaa.osoaa_root()
    except FileNotFoundError:
        root = osoaa.DEFAULT_ROOT
    return np.loadtxt(os.path.join(root, 'fic', name))


def water_absorption(lam_nm):
    """Pure-water absorption (m-1) from OSOAA's table
    (``fic/OSOAA_SEA_MOL_COEFFS_JUNE_2013.txt``, column 2)."""
    t = _osoaa_table('OSOAA_SEA_MOL_COEFFS_JUNE_2013.txt')
    return np.interp(lam_nm, t[:, 0], t[:, 1])


def total_absorption(lam_nm, chl=0.0, ys440=0.0):
    """a_w + a_ph + a_ys (m-1).  a_ph = A(lam) Chl^E(lam) from OSOAA's Bricaud
    table (``fic/OSOAA_SEA_PHYT_COEFFS.txt``; held at its 400 nm values below
    400); a_ys = ys440 exp(-0.014 (lam - 440))."""
    lam = np.asarray(lam_nm, dtype=float)
    a = water_absorption(lam)
    if chl > 0:
        t = _osoaa_table('OSOAA_SEA_PHYT_COEFFS.txt')
        A, E = np.interp(lam, t[:, 0], t[:, 1]), np.interp(lam, t[:, 0], t[:, 2])
        a = a + A * chl ** E
    if ys440 > 0:
        a = a + ys440 * np.exp(-0.014 * (lam - 440.0))
    return a


def raman_efficiency(lam_nm, at550=RAMAN_EFFICIENCY, chl=0.0, ys440=0.0,
                     shift_cm=RAMAN_SHIFT):
    """Raman reflectance factor at emission ``lam_nm``:
    C b_R(exc) / (2 a_t(exc) + a_t(lam)), b_R ~ exc^-5.5, with C set so that
    *pure water* gives ``at550`` at 550 nm.  Absorption by phytoplankton and
    yellow substance (``chl``, ``ys440``) therefore only lowers it."""
    lam = np.asarray(lam_nm, dtype=float)

    def f(l, c, y):
        x = raman_excitation(l, shift_cm)
        return x ** -5.5 / (2 * total_absorption(x, c, y) + total_absorption(l, c, y))
    return at550 * f(lam, chl, ys440) / f(550.0, 0.0, 0.0)


def rhow_library(fields, ref_case=None):
    """Elastic rho_w per water type, from the OSOAA water grid, plus controls.

    Parameters
    ----------
    fields : dict
        :func:`load_fields`.
    ref_case : str, optional
        Geometry/aerosol prefix for the reference rho_w of each water type
        (default ``'sza50_aot0.05'``; rho_w depends only weakly on them).

    Returns
    -------
    dict
        water -> dict(``lam5``, ``rho_el`` (5 nm), ``chl``, ``fl_amp``,
        ``raman_eff``).
    """
    from hypernet.wiggles.phase1_osoaa_fields import WATERS
    ref_case = ref_case or 'sza50_aot0.05'
    lam5 = fields['lam']
    out = {}
    for w, props in WATERS.items():
        i = _case_index(fields, '%s_%s' % (ref_case, w))
        ed = fields['ed_dir'][i] + fields['ed_dif'][i]
        out[w] = dict(lam5=lam5, rho_el=np.pi * fields['lw'][i] / ed,
                      chl=props.get('chl', 0.0), ys440=props.get('ys440', 0.0),
                      fl_amp=float(fluorescence_amplitude(props.get('chl', 0.0))),
                      raman_eff=RAMAN_EFFICIENCY)
    return out


def raman_excitation(lam_nm, shift_cm=RAMAN_SHIFT):
    """Excitation wavelength (nm) whose Raman emission lands at ``lam_nm``."""
    return 1.0 / (1.0 / np.asarray(lam_nm, dtype=float) + shift_cm * 1e-7)


def raman_radiance(lam, Ed, efficiency=RAMAN_EFFICIENCY, shift_cm=RAMAN_SHIFT,
                   fwhm_cm=RAMAN_FWHM):
    """Water-Raman Lw: Ed at the excitation wavelengths, broadened by the
    Raman band (Gaussian in wavenumber), times ``efficiency / pi``.

    ``lam`` must extend ~(shift) blueward of the band of interest; outside the
    excitation range the result is NaN-free but uses the edge of ``Ed``.
    """
    lam = np.asarray(lam, dtype=float)
    nu = 1e7 / lam                                  # cm-1, descending
    # Ed per unit wavenumber would add a lam^2 factor; a constant efficiency
    # absorbs it at this level of approximation
    nu_exc = nu + shift_cm
    order = np.argsort(nu)
    nu_s, ed_s = nu[order], np.asarray(Ed, dtype=float)[order]
    step = np.median(np.diff(nu_s))
    sig = fwhm_cm / srfmod.FWHM_PER_SIGMA
    half = int(np.ceil(4 * sig / step))
    kern = np.exp(-0.5 * (np.arange(-half, half + 1) * step / sig) ** 2)
    kern /= kern.sum()
    # resample Ed on a uniform wavenumber grid, broaden, then read at nu_exc
    nu_u = np.arange(nu_s[0], nu_s[-1], step)
    ed_u = np.interp(nu_u, nu_s, ed_s)
    ed_b = np.convolve(ed_u, kern, mode='same')
    return efficiency / np.pi * np.interp(nu_exc, nu_u, ed_b)


def compose_scene(emod, fields, case, water=None, controls=('fl', 'raman'), lib=None):
    """High-resolution Ed, Ld, Lw (and Lu) for one OSOAA case.

    Parameters
    ----------
    emod : dict
        :func:`hypernet.emod.build_emod` output (``lam``, ``F0``,
        ``T_direct``, ``T_diffuse``); its SZA should match the case's.
    fields : dict
        :func:`load_fields`.
    case : str or int
        OSOAA case (name or index): geometry, aerosol, and default water.
    water : str, optional
        Water type for Lw (default: the case's own).
    controls : tuple
        Any of ``'fl'`` (fluorescence) and ``'raman'``.
    lib : dict, optional
        :func:`rhow_library` (computed if absent).

    Returns
    -------
    dict
        ``lam``, ``Ed``, ``Ed_dir``, ``Ed_dif``, ``Ld``, ``Lw``,
        ``Lw_elastic``, ``Lw_fl``, ``Lw_raman``, ``Lu``, ``rho_eff`` (smooth,
        on ``lam``), ``rho_el`` (smooth elastic rho_w on ``lam``), and the
        case metadata.
    """
    i = _case_index(fields, case)
    lam5, lam = fields['lam'], emod['lam']
    F0 = emod['F0']
    s = {k: _spline_log(lam5, fields[k][i], lam) for k in ('ed_dir', 'ed_dif', 'ld')}
    Ed_dir = F0 * s['ed_dir'] * emod['T_direct']
    Ed_dif = F0 * s['ed_dif'] * emod['T_diffuse']
    Ed = Ed_dir + Ed_dif
    Ld = F0 * s['ld'] * emod['T_diffuse']
    lib = lib or rhow_library(fields)
    water = water or str(np.asarray(fields['water'])[i])
    w = lib[water]
    rho_el = _spline_log(lam5, w['rho_el'], lam)
    Lw_el = rho_el * Ed / np.pi
    Lw_fl = np.zeros_like(lam)
    Lw_r = np.zeros_like(lam)
    if 'fl' in controls:
        sig = FL_FWHM / srfmod.FWHM_PER_SIGMA
        # fluorescence is emission: its radiance follows the smooth (line-free)
        # irradiance level -- Ed smoothed with a 5 nm Gaussian -- not Ed's lines
        from scipy.ndimage import gaussian_filter1d
        step = float(np.median(np.diff(lam)))
        ed_smooth = gaussian_filter1d(Ed, 5.0 / step, mode='nearest')
        Lw_fl = w['fl_amp'] * np.exp(-0.5 * ((lam - FL_CENTRE) / sig) ** 2) * ed_smooth / np.pi
    if 'raman' in controls:
        Lw_r = raman_radiance(lam, Ed, efficiency=raman_efficiency(
            lam, w['raman_eff'], chl=w.get('chl', 0.0), ys440=w.get('ys440', 0.0)))
    Lw = Lw_el + Lw_fl + Lw_r
    rho_eff = _spline_log(lam5, (fields['lu0p'][i] - fields['lw'][i]) / fields['ld'][i], lam)
    return dict(lam=lam, Ed=Ed, Ed_dir=Ed_dir, Ed_dif=Ed_dif, Ld=Ld, Lw=Lw, Lw_elastic=Lw_el,
                Lw_fl=Lw_fl, Lw_raman=Lw_r, Lu=Lw + rho_eff * Ld, rho_eff=rho_eff,
                rho_el=rho_el, case=str(np.asarray(fields['name'])[i]), water=water,
                sza=float(np.asarray(fields['sza'])[i]), controls=tuple(controls))


def rhow_true(scene, srf_L, grid_L):
    """The true rho_w: Lw and Ed both observed with the L SRF on the L grid.

    Parameters
    ----------
    scene : dict
        :func:`compose_scene`.
    srf_L : float, array, callable or SRFModel
        L-channel FWHM (nm): anything :func:`convolve_to_grid` accepts, or an
        :class:`hypernet.srf.SRFModel` (its ``fwhm`` is used).
    grid_L : array
        L pixel centres, nm.

    Returns
    -------
    numpy.ndarray
        pi Lw_L / Ed_L on ``grid_L``.
    """
    f = srf_L.fwhm if hasattr(srf_L, 'fwhm') else srf_L
    lw = convolve_to_grid(scene['lam'], scene['Lw'], grid_L, f)
    ed = convolve_to_grid(scene['lam'], scene['Ed'], grid_L, f)
    return np.pi * lw / ed


# --- instrument model (Phase 1 task 7) ---------------------------------------------

import dataclasses  # noqa: E402
import os  # noqa: E402


def data_path(name):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', name)


def gaussian_srf(dlam, fwhm):
    """The Gaussian SRF as a density in nm-1 (it integrates to 1 over dlam)."""
    sig = np.asarray(fwhm, dtype=float) / srfmod.FWHM_PER_SIGMA
    return np.exp(-0.5 * (np.asarray(dlam) / sig) ** 2) / (sig * np.sqrt(2 * np.pi))


def load_grids(path=None):
    """``hypernet/data/hypstar_grids.npz``: the real pixel grids and the VEIT
    per-scan relative noise (written by ``hypernet/wiggles/phase1_instrument.py``).

    Keys: ``grid_E_122304`` (1536), ``grid_L_122304`` (1538),
    ``grid_L_122305`` (1536), and the relative per-scan noise
    ``noise_E_scan``, ``noise_Ld_scan``, ``noise_Lu_scan`` on the 122304 grids
    (flattened scan scatter, Phase 0 Q9).
    """
    d = np.load(path or data_path('hypstar_grids.npz'))
    return {k: d[k] for k in d.files}


def load_instrument_srfs(path=None):
    """``hypernet/data/release2_srf_models.json``: instrument -> {E, Ld}
    :class:`hypernet.srf.SRFModel` (clear-sky Release 2 template fits)."""
    import json
    path = path or data_path('release2_srf_models.json')
    doc = json.load(open(path))
    return {inst: {ch: srfmod.SRFModel.from_dict(m) for ch, m in chans.items()}
            for inst, chans in doc['models'].items()}


def observe(lam_hr, spec_hr, grid, srf_model, seed=None, noise=None, true_grid=None,
            dfwhm=0.0):
    """Observe a high-resolution spectrum with a Gaussian SRF on a pixel grid.

    Parameters
    ----------
    lam_hr, spec_hr : array
        The 0.01 nm scene spectrum.
    grid : array
        The *reported* pixel wavelengths (what the instrument says).
    srf_model : SRFModel, callable, array or float
        The SRF FWHM (nm), as in :func:`convolve_to_grid`; for an
        :class:`hypernet.srf.SRFModel` its ``fwhm(lambda)`` is used.
    seed : int, optional
        RNG seed for the noise.
    noise : float or array, optional
        Relative 1-sigma noise per pixel (e.g. ``noise_E_scan`` /
        sqrt(n_scans)).  None = noise-free.
    true_grid : array, optional
        Where the pixels really sample, if not at ``grid`` (a wavelength-
        calibration error: case iv).
    dfwhm : float
        Added to the FWHM (nm).

    Returns
    -------
    numpy.ndarray
        The observed spectrum on ``grid``.
    """
    f = srf_model.fwhm if hasattr(srf_model, 'fwhm') else srf_model
    at = grid if true_grid is None else true_grid
    if dfwhm:
        f0 = f
        f = (lambda l: f0(l) + dfwhm) if callable(f0) else np.asarray(f0) + dfwhm
    obs = convolve_to_grid(lam_hr, spec_hr, at, f)
    if noise is not None:
        rng = np.random.default_rng(seed)
        obs = obs * (1.0 + np.asarray(noise) * rng.standard_normal(obs.size))
    return obs


@dataclasses.dataclass
class Case:
    """One twin-experiment instrument case.

    The truth is always observed with ``srf_L`` on ``grid_L``.  E is observed
    with ``srf_E`` at ``true_grid_E`` (``grid_E`` plus any wavelength error)
    but reported on ``grid_E``.  The fields ``corr_*`` describe what the
    correction is *told*.
    """
    name: str
    label: str
    instrument: str
    srf_E: object
    srf_L: object
    grid_E: np.ndarray
    grid_L: np.ndarray
    true_grid_E: np.ndarray = None
    corr_dfwhm_E: float = 0.0
    corr_dfwhm_L: float = 0.0
    emod_variant: dict = dataclasses.field(default_factory=dict)
    noise: bool = False
    n_scans: int = 6


def case_table(instrument='HYPSTAR_122304', grids=None, srfs=None):
    """The cases (i)-(vii) of the plan for one instrument.

    - (i) ``offset``: equal SRFs (both = the L SRF), real E and L grids.
    - (ii) ``mismatch``: the E and L SRFs, E observed on the L grid (no
      offset).
    - (iii) ``both``: the E and L SRFs, real grids.
    - (iv) ``wavecal_*``: as (iii), with E's true grid shifted rigidly by
      +-0.1 and +-0.3 nm, or stretched linearly by +-0.1 nm at 390/870 nm
      (Phase 0 8d found a linear term).
    - (v) ``wrongsrf_*``: as (iii), with the correction told FWHM_E or
      FWHM_L +-0.3 nm.
    - (vi) ``wrongemod_*``: as (iii), with the correction's Emod built for SZA
      + 10, water vapour x 2, or the low-aerosol fields.
    - (vii) ``noise``: as (iii), with photon noise at the VEIT per-scan level
      for a mean of 6 scans.
    """
    grids = grids or load_grids()
    srfs = srfs or load_instrument_srfs()
    sE, sL = srfs[instrument]['E'], srfs[instrument]['Ld']
    tag = instrument.replace('HYPSTAR_', '')
    gE = grids['grid_E_122304']
    gL = grids.get('grid_L_%s' % tag, grids['grid_L_122304'])
    C = []
    C.append(Case('i_offset', '(i) grid offset only, equal SRFs', instrument, sL, sL, gE, gL))
    C.append(Case('ii_mismatch', '(ii) SRF mismatch, no offset', instrument, sE, sL, gL, gL))
    C.append(Case('iii_both', '(iii) offset + SRF mismatch', instrument, sE, sL, gE, gL))
    for d in (-0.3, -0.1, 0.1, 0.3):
        C.append(Case('iv_shift%+.1f' % d, '(iv) E wavelength error %+.1f nm (rigid)' % d,
                      instrument, sE, sL, gE, gL, true_grid_E=gE + d))
    for d in (-0.1, 0.1):
        stretch = d * (gE - 630.0) / 240.0           # +-d at 390 / 870 nm
        C.append(Case('iv_stretch%+.1f' % d, '(iv) E wavelength stretch %+.1f nm at the '
                      'ends' % d, instrument, sE, sL, gE, gL, true_grid_E=gE + stretch))
    for ch in ('E', 'L'):
        for d in (-0.3, 0.3):
            kw = {'corr_dfwhm_%s' % ch: d}
            C.append(Case('v_wrongsrf_%s%+.1f' % (ch, d), '(v) correction told FWHM_%s %+.1f '
                          'nm' % (ch, d), instrument, sE, sL, gE, gL, **kw))
    for key, v, lab in (('sza_offset', 10.0, 'SZA + 10 deg'), ('pwv_factor', 2.0,
                                                               'water vapour x 2'),
                        ('aerosol', 'low', 'no aerosol (AOT 0.05 fields)')):
        C.append(Case('vi_wrongemod_%s' % key, '(vi) wrong Emod: %s' % lab, instrument, sE, sL,
                      gE, gL, emod_variant={key: v}))
    C.append(Case('vii_noise', '(vii) VEIT per-scan noise, mean of 6 scans', instrument, sE,
                  sL, gE, gL, noise=True))
    return C

"""Phase 0a, task 7: preliminary H1 / H2 wiggle budget at each line (VEIT).

At each line (or joint-fit group: Ca H/K, Ca II 849.8/854.2) the relative
residual profile r(x) = y / continuum - 1 is computed within +-5 nm, the
continuum being a straight line through the pixels 5-7 nm either side.  For:

- **measured** -- Ld / Ed_L and Lu / Ed_L, with Ld interpolated in time to the
  water view and Ed_L = np.interp of E onto the L grid, i.e. exactly the
  quantities the processor divides (L1C), and rho_w from L2A;
- **H1 (interpolation)** -- Ed_cubic / Ed_linear - 1 on the L grid: what the
  ratio gains from linear instead of cubic interpolation of the same Ed (as in
  ``wavecal/sanity_checks_veit.py``);
- **H2 (SRF mismatch)** -- the E line re-observed at the L width with its
  equivalent width conserved, over the E line:
  ``prod_k (1 - a_k s_E,k/s_L,k g(x; mu_k, s_L,k)) / (1 - a_k g(x; mu_k, s_E,k)) - 1``,
  from the task 5 fits of E (depth a, sigma s_E, centroid mu) and of Ld
  (sigma s_L).  Measured L depths are not used, so the Ring effect does not
  enter the prediction.

Reported per line:

- the rms over +-5 nm of each profile, and of the noise (flattened errors);
- the task's closed form, depth x (FWHM_L^2 - FWHM_E^2) / FWHM_L^2, next to
  the H2 profile's central value.  The closed form is about twice the central
  value, because to first order the central value is a (1 - s_E/s_L)
  ~ a (s_L^2 - s_E^2) / (2 s_L^2);
- weighted least-squares projections of each measured profile: onto H2
  alone (alpha2) and H1 alone (alpha1) -- alpha = 1 means that hypothesis
  alone accounts for the feature -- and jointly onto (H1, H2, shift), where
  shift = d ln(E line model)/d lambda, so that alpha3 is an effective L - E
  wavelength offset in nm.  Without the shift term H1 absorbs any residual
  offset and the joint fit is degenerate;
- rho_w in absolute units (rho_w is ~0 in the NIR, so relative residuals mean
  nothing there): rms of the measured residual (rho_w - continuum) against
  the predictions rho_w x rms(H2) and rho_w x rms(H1), and the noise from
  L2A std_reflectance / sqrt(n_valid_scans);
- rms(rho_w'') within +-5 nm.

Outputs: ``hypernet/wiggles/phase0_veit_budget.csv`` and
``hypernet/wiggles/figs/phase0/veit_budget.png``.

This is the template for task 12: :func:`budget` works on any sequence given
the loaded L1A files, the L2A dict and the task-5-style line fits.

Run from the repository root: ``python -m hypernet.wiggles.phase0a_budget_veit``.
"""
import os

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

from hypernet import srf  # noqa: E402
from hypernet import whn_l1a as wl  # noqa: E402
from hypernet.wiggles import phase0a_veit as p5  # noqa: E402
from hypernet.wiggles import DATA_DIR, FIGDIR, OUT, REPO, WIGGLES_DIR  # noqa: F401

HALF = 5.0        # nm, the budget window
EDGE = 2.0        # nm, continuum pixels from HALF to HALF + EDGE

# Colours: palette slots in fixed order (measured Ld = the reference, dark).
C = {'meas': '#333333', 'H1': '#2a78d6', 'H2': '#eb6834', 'rho': '#1baf7a',
     'Lu': '#888888', 'noise': '#bbbbbb'}


def _residual(x, y, lo, hi):
    """y / (linear continuum through the edge pixels) - 1, and the in-window mask."""
    edge = ((x > lo - EDGE) & (x < lo)) | ((x > hi) & (x < hi + EDGE))
    edge &= np.isfinite(y)
    inwin = (x >= lo) & (x <= hi)
    if edge.sum() < 2:
        return np.full_like(y, np.nan), inwin
    c = np.polyfit(x[edge], y[edge], 1)
    return y / np.polyval(c, x) - 1.0, inwin


def _h2_profile(x, comps):
    """EW-conserving re-observation of the E line(s) at the L width, over E."""
    num = np.ones_like(x)
    den = np.ones_like(x)
    for a, mu, sE, sL in comps:
        num *= 1.0 - a * sE / sL * np.exp(-0.5 * ((x - mu) / sL) ** 2)
        den *= 1.0 - a * np.exp(-0.5 * ((x - mu) / sE) ** 2)
    return num / den - 1.0


def _shift_profile(x, comps):
    """d ln(E line model) / d lambda (per nm): the ratio residual of a unit
    L - E wavelength offset, up to sign."""
    lnE = np.zeros_like(x)
    for a, mu, sE, _ in comps:
        lnE += np.log(1.0 - a * np.exp(-0.5 * ((x - mu) / sE) ** 2))
    return np.gradient(lnE, x)


def _project(r, e, basis):
    """Weighted LS of r on the basis columns; returns coefficients and errors."""
    m = np.isfinite(r) & np.isfinite(e) & (e > 0) & np.all(np.isfinite(basis), axis=1)
    if m.sum() <= basis.shape[1]:
        return np.full(basis.shape[1], np.nan), np.full(basis.shape[1], np.nan)
    A = basis[m] / e[m, None]
    b = r[m] / e[m]
    cov = np.linalg.pinv(A.T @ A)
    coef = cov @ (A.T @ b)
    chi2 = np.sum((b - A @ coef) ** 2) / max(m.sum() - basis.shape[1], 1)
    return coef, np.sqrt(np.diag(cov) * max(chi2, 1.0))


def _rms(r, m):
    v = r[m]
    v = v[np.isfinite(v)]
    return float(np.sqrt(np.mean(v ** 2))) if v.size else np.nan


def budget(irr, rad, l2a, fits_E, fits_L, lines=srf.LINES, return_profiles=False):
    """The H1/H2 budget of one sequence.

    Parameters
    ----------
    irr, rad : dict
        :func:`hypernet.whn_l1a.load_l1a_irr` / ``load_l1a_rad`` outputs.
    l2a : dict
        :func:`hypernet.whn_l1a.load_l2a` output.
    fits_E, fits_L : pandas.DataFrame
        :func:`hypernet.srf.fit_lines` of mean E and mean Ld.
    lines : pandas.DataFrame
        Line table (groups define the windows).
    return_profiles : bool
        Also return the per-group residual profiles (for plotting).

    Returns
    -------
    pandas.DataFrame (and dict of profiles if ``return_profiles``)
    """
    wL, wE = rad['wave'], irr['wave']
    t_lu = rad['Lu']['t_mean']
    E_t = wl.interp_in_time(irr['scans'], irr['series_id'], irr['time'], t_lu)
    E_lin = np.interp(wL, wE, E_t)
    E_cub = CubicSpline(wE, E_t)(wL)
    Ld = rad['Ld']['at_lu_time']
    Lu = rad['Lu']['mean']
    rho = l2a['reflectance'][:, 0]
    n_rho = float(np.atleast_1d(l2a.get('n_valid_scans', [1]))[0])
    e_rho = l2a['std_reflectance'][:, 0] / np.sqrt(max(n_rho, 1.0))
    # errors of the ratios, from the flattened scan scatter (relative, quadrature)
    _, eE, _ = srf.scan_errors(irr['scans'], flatten_px=p5.FLATTEN_PX)
    eE_L = np.interp(wL, wE, eE / E_t)
    _, eLd, _ = srf.scan_errors(rad['Ld']['scans'], flatten_px=p5.FLATTEN_PX)
    _, eLu, _ = srf.scan_errors(rad['Lu']['scans'], flatten_px=p5.FLATTEN_PX)
    e_ld = np.hypot(eLd / Ld, eE_L)
    e_lu = np.hypot(eLu / Lu, eE_L)
    rho2 = np.gradient(np.gradient(rho, wL), wL)
    h1_full = E_cub / E_lin - 1.0

    fE = fits_E.set_index('name')
    fL = fits_L.set_index('name')
    rows, prof = [], {}
    for grp, g in lines.groupby('group', sort=False):
        names = list(g['name'])
        lo, hi = g['lam_air'].min() - HALF, g['lam_air'].max() + HALF
        ok = all(fE.loc[n, 'ok'] and fL.loc[n, 'ok'] for n in names)
        comps = [(fE.loc[n, 'depth'], fE.loc[n, 'mu'], fE.loc[n, 'sigma'],
                  fL.loc[n, 'sigma']) for n in names] if ok else []
        r_ld, m = _residual(wL, Ld / E_lin, lo, hi)
        r_lu, _ = _residual(wL, Lu / E_lin, lo, hi)
        r_rho, _ = _residual(wL, rho, lo, hi)
        r_h1, _ = _residual(wL, 1.0 + h1_full, lo, hi)
        r_h2 = _residual(wL, 1.0 + _h2_profile(wL, comps), lo, hi)[0] if ok \
            else np.full_like(wL, np.nan)
        r_sh = _shift_profile(wL, comps) if ok else np.full_like(wL, np.nan)
        # rho_w absolute residual and its continuum
        rho_cont = rho / (1.0 + r_rho)
        d_rho = rho - rho_cont
        e_rho_rel = e_rho / rho_cont
        one = lambda b: np.column_stack([b])[m]  # noqa: E731
        three = np.column_stack([r_h1, r_h2, r_sh])[m]
        res = {}
        for key, r, e in (('LdEd', r_ld, e_ld), ('LuEd', r_lu, e_lu),
                          ('rho', r_rho, e_rho_rel)):
            a2, ea2 = _project(r[m], e[m], one(r_h2))
            a1, ea1 = _project(r[m], e[m], one(r_h1))
            a3, ea3 = _project(r[m], e[m], three)
            res[key] = dict(alpha2=a2[0], alpha2_err=ea2[0],
                            alpha1=a1[0], alpha1_err=ea1[0],
                            j_alpha1=a3[0], j_alpha1_err=ea3[0],
                            j_alpha2=a3[1], j_alpha2_err=ea3[1],
                            j_shift=a3[2], j_shift_err=ea3[2])
        fin = m & np.isfinite(r_h1) & np.isfinite(r_h2)
        corr12 = float(np.corrcoef(r_h1[fin], r_h2[fin])[0, 1]) if fin.sum() > 2 \
            else np.nan
        label = 'Ca H/K' if grp == 'CaHK' else ('Ca II 849.8/854.2' if grp == 'CaIR12'
                                                else names[0])
        # the closed form of the task, and the H2 profile's central value, at
        # the strongest component
        k = int(np.argmax([c[0] for c in comps])) if ok else 0
        if ok:
            a, mu, sE, sL = comps[k]
            closed = a * (sL ** 2 - sE ** 2) / sL ** 2
            central = float(_h2_profile(np.array([mu]), [comps[k]])[0])
        else:
            closed = central = np.nan
        rows.append(dict(
            group=grp, line=label, lam_air=float(g['lam_air'].mean()),
            blend=bool(g['blend'].any()), use_for_srf=bool(g['use_for_srf'].all()),
            fit_ok=ok, depth_E=comps[k][0] if ok else np.nan,
            fwhm_E=FW * comps[k][2] if ok else np.nan,
            fwhm_L=FW * comps[k][3] if ok else np.nan,
            h2_closed_form=closed, h2_central=central,
            rms_meas_LdEd=_rms(r_ld, m), rms_meas_LuEd=_rms(r_lu, m),
            rms_noise_LdEd=_rms(e_ld, m), rms_H1=_rms(r_h1, m), rms_H2=_rms(r_h2, m),
            rms_rho_rel=_rms(r_rho, m), rms_rho2=_rms(rho2, m),
            rho_median=float(np.nanmedian(rho[m])),
            rms_rho_abs=_rms(d_rho, m),
            rms_rho_pred_H2=_rms(rho_cont * r_h2, m),
            rms_rho_pred_H1=_rms(rho_cont * r_h1, m),
            rms_rho_noise=_rms(e_rho, m), corr_H1_H2=corr12,
            **{'%s_%s' % (k2, key): v for key, d in res.items() for k2, v in d.items()}))
        prof[label] = dict(x=wL[m], ld=r_ld[m], lu=r_lu[m], rho=r_rho[m],
                           h1=r_h1[m], h2=r_h2[m], e=e_ld[m])
    df = pd.DataFrame(rows)
    return (df, prof) if return_profiles else df


FW = srf.FWHM_PER_SIGMA


def figure(df, prof, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(10, 13))
    gs = fig.add_gridspec(4, 4, height_ratios=[1.1, 1.2, 1.0, 1.0], hspace=0.35,
                          wspace=0.3)
    # top: example profiles
    for i, lab in enumerate(['G band', 'H beta', 'Na D', 'H alpha']):
        ax = fig.add_subplot(gs[0, i])
        p = prof[lab]
        ax.fill_between(p['x'], -p['e'], p['e'], color=C['noise'], alpha=0.5, lw=0)
        ax.plot(p['x'], p['ld'], color=C['meas'], lw=1.8, label='Ld/Ed measured')
        ax.plot(p['x'], p['h2'], color=C['H2'], lw=2, label='H2 predicted')
        ax.plot(p['x'], p['h1'], color=C['H1'], lw=2, label='H1 predicted')
        ax.axhline(0, color='#888', lw=0.6)
        ax.set_title(lab, fontsize=9.5, loc='left')
        ax.tick_params(labelsize=8)
        ax.set_xlabel('λ (nm)', fontsize=8.5)
        if i == 0:
            ax.set_ylabel('relative residual', fontsize=8.5)
        ax.spines[['top', 'right']].set_visible(False)
    fig.axes[0].legend(fontsize=7.5, frameon=False, loc='lower left')

    x = np.arange(len(df))
    # row 2: rms per line, relative ratios
    ax = fig.add_subplot(gs[1, :])
    series = [('rms_meas_LdEd', 'Ld/Ed measured', C['meas'], 'o'),
              ('rms_meas_LuEd', 'Lu/Ed measured', C['Lu'], 'v'),
              ('rms_H2', 'H2 predicted (SRF mismatch)', C['H2'], 's'),
              ('rms_H1', 'H1 predicted (linear interpolation)', C['H1'], '^'),
              ('rms_noise_LdEd', 'noise (Ld/Ed)', C['noise'], '_')]
    for j, (col, lab, c, mk) in enumerate(series):
        ax.plot(x + (j - 2) * 0.1, df[col], ls='', marker=mk, ms=8 if mk != '_' else 14,
                color=c, mew=2 if mk == '_' else 1, label=lab)
    ax.set_yscale('log')
    ax.set_ylabel('rms relative residual\nwithin ±5 nm')
    ax.legend(fontsize=8, frameon=False, ncol=3, loc='lower left')
    ax.set_title('VEIT 2026-06-04, HYPSTAR 122304: residual line structure vs the '
                 'H1 and H2 predictions', fontsize=10.5, loc='left')
    # row 3: rho_w, absolute
    cx = fig.add_subplot(gs[2, :], sharex=ax)
    for j, (col, lab, c, mk) in enumerate([
            ('rms_rho_abs', 'ρw measured', C['rho'], 'D'),
            ('rms_rho_pred_H2', 'ρw × H2', C['H2'], 's'),
            ('rms_rho_pred_H1', 'ρw × H1', C['H1'], '^'),
            ('rms_rho_noise', 'noise (L2A std / √n)', C['noise'], '_')]):
        cx.plot(x + (j - 1.5) * 0.1, df[col], ls='', marker=mk, ms=8 if mk != '_' else 14,
                color=c, mew=2 if mk == '_' else 1, label=lab)
    cx.set_yscale('log')
    cx.set_ylabel('rms ρw residual\nwithin ±5 nm')
    cx.legend(fontsize=8, frameon=False, ncol=4, loc='lower left')
    # row 4: single-hypothesis projections
    bx = fig.add_subplot(gs[3, :], sharex=ax)
    for j, (pre, lab, c, mk) in enumerate([('LdEd', 'Ld/Ed', C['meas'], 'o'),
                                           ('LuEd', 'Lu/Ed', C['Lu'], 'v'),
                                           ('rho', 'ρw', C['rho'], 'D')]):
        bx.errorbar(x + (j - 1) * 0.15, df['alpha2_' + pre], yerr=df['alpha2_err_' + pre],
                    fmt=mk, ms=7, color=c, lw=1.2, label='%s' % lab)
    bx.axhline(1, color=C['H2'], lw=0.8, ls='--')
    bx.axhline(0, color='#888', lw=0.8)
    bx.set_ylim(-1, 2.5)
    bx.set_ylabel('α₂: fraction of the\nH2 profile present')
    bx.legend(fontsize=8, frameon=False, ncol=3, loc='upper left')
    bx.set_xticks(x)
    bx.set_xticklabels(df['line'], rotation=30, ha='right', fontsize=8.5)
    for a in (ax, cx):
        plt.setp(a.get_xticklabels(), visible=False)
    for a in (ax, cx, bx):
        a.grid(alpha=0.25, lw=0.6)
        a.spines[['top', 'right']].set_visible(False)
        a.tick_params(labelsize=9)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)


def main():
    pd.set_option('display.width', 220)
    irr, rad, chans = p5.load()
    files = wl.sequence_files(site='VEIT', seq_time='20260604T0845')
    l2a = wl.load_l2a(files['L2A_REF'])
    lines = p5.fit_all(chans)
    f = lines[lines['err_kind'] == 'flat']
    fits_E = f[f['channel'] == 'E'].reset_index(drop=True)
    fits_L = f[f['channel'] == 'Ld'].reset_index(drop=True)
    df, prof = budget(irr, rad, l2a, fits_E, fits_L, return_profiles=True)
    df.to_csv(os.path.join(WIGGLES_DIR, 'phase0_veit_budget.csv'), index=False,
              float_format='%.6g')
    print('\nrms relative residual within +-5 nm (and rho_w\'\' rms, 1/nm^2):')
    print(df[['line', 'rms_meas_LdEd', 'rms_meas_LuEd', 'rms_noise_LdEd', 'rms_H1',
              'rms_H2', 'rms_rho_rel', 'rms_rho2']].to_string(index=False,
                                                                float_format='%.2e'))
    print('\nH2 closed form (task) vs the H2 profile central value:')
    print(df[['line', 'depth_E', 'fwhm_E', 'fwhm_L', 'h2_closed_form',
              'h2_central']].to_string(index=False, float_format='%.3f'))
    print("\nrho_w absolute: measured residual rms vs rho_w x H2 / H1, noise:")
    print(df[['line', 'rho_median', 'rms_rho_abs', 'rms_rho_pred_H2', 'rms_rho_pred_H1',
              'rms_rho_noise']].to_string(index=False, float_format='%.2e'))
    print('\nsingle-hypothesis projections alpha2 (H2 alone), alpha1 (H1 alone):')
    cols = ['line']
    for pre in ('LdEd', 'LuEd', 'rho'):
        cols += ['alpha2_' + pre, 'alpha2_err_' + pre, 'alpha1_' + pre, 'alpha1_err_' + pre]
    print(df[cols].to_string(index=False, float_format='%.2f'))
    print('\njoint projection onto (H1, H2, shift); shift in nm; corr(H1, H2) of the bases:')
    cols = ['line', 'corr_H1_H2']
    for pre in ('LdEd', 'rho'):
        cols += ['j_alpha1_' + pre, 'j_alpha1_err_' + pre, 'j_alpha2_' + pre,
                 'j_alpha2_err_' + pre, 'j_shift_' + pre, 'j_shift_err_' + pre]
    print(df[cols].to_string(index=False, float_format='%.2f'))
    clean = df[df['use_for_srf'] & df['fit_ok']]
    for pre in ('LdEd', 'LuEd', 'rho'):
        for a in ('alpha2', 'alpha1', 'j_alpha1', 'j_alpha2', 'j_shift'):
            w = 1 / clean['%s_err_%s' % (a, pre)] ** 2
            print('  %-4s %-8s weighted mean over the SRF lines: %+.3f +- %.3f' % (
                pre, a, np.sum(w * clean['%s_%s' % (a, pre)]) / w.sum(),
                1 / np.sqrt(w.sum())))
    r = clean['rms_H2'] / clean['rms_H1']
    print('  rms_H2 / rms_H1 over the SRF lines: median %.1f, range %.1f-%.1f'
          % (r.median(), r.min(), r.max()))
    figure(df, prof, os.path.join(FIGDIR, 'veit_budget.png'))
    print('\nwrote hypernet/wiggles/phase0_veit_budget.csv, hypernet/wiggles/figs/phase0/veit_budget.png')


if __name__ == '__main__':
    main()

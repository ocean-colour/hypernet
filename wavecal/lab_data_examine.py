"""Examine Joel's HYPSTAR laboratory wavelength-calibration files (lamp lines)
against our field-derived SRF table and wavelength grids.

Per instrument x lab epoch, `$OS_COLOR/hypernet/wiggles/lab_data/hypstar_<sn>/
calibration/wavelength/<YYYY_MM>/` holds `deltawl_vnir.dat`, `fwhm_vnir.dat`,
`hypstar_<sn>_wl_coefs_<yymmdd>.dat` and `lamp_vnir.dat`.  This script

1. inventories the epochs against our calibration periods (srf_table, timeline);
2. works out the file semantics (what `deltawl` is, the lamp pairs in
   `lamp_vnir.dat`, the wavelength grid of its `wl` column, air vs vacuum);
3. evaluates the polynomials on px 0..2047 and finds the crops that reproduce
   the VEIT L1A E/L grids and every Release 2 L grid, then builds each
   instrument's own E grid (replacing the 122304 proxy);
4. compares the lab FWHM_E / FWHM_L at the lines with our srf_table model;
5. fits the lamp lines with Gaussian / super-Gaussian / pseudo-Voigt shapes;
6. looks at the wings far from the lines (stray light) and at whether the
   polynomial between the lamp lines could carry the 0.13 nm window pattern;
7. checks 121222 `2021_07_tst` vs `2021_07` and 120242's four epochs.

Outputs (CSV/PNG, not committed) go to
`$OS_COLOR/hypernet/wiggles/lab_data/derived/`; a text summary is printed and
saved as `summary.txt`.

Run from the repository root:
    conda run -n ocean14 python wavecal/lab_data_examine.py
"""
import glob
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from hypernet import srf as srfmod
from hypernet.srftable import load_srf_table
from hypernet.whn_l1a import release2_path, load_l2b, load_l1a_irr, load_l1a_rad

OS_COLOR = os.environ['OS_COLOR']
LAB = os.path.join(OS_COLOR, 'hypernet', 'wiggles', 'lab_data')
OUT = os.path.join(LAB, 'derived')
os.makedirs(OUT, exist_ok=True)
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VEIT_DIR = os.path.join(OS_COLOR, 'WATERHYPERNET', 'Wavelengths', '')
NPIX = 2048
FWHM_PER_SIGMA = srfmod.FWHM_PER_SIGMA
# fixed categorical order (dataviz reference palette)
COL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']

# Species of the reference lines in deltawl/fwhm files (NIST air wavelengths,
# identified by value; "inferred" in the report).
SPECIES = {296.728: 'Hg I', 404.6563: 'Hg I', 435.8328: 'Hg I', 546.0735: 'Hg I',
           1013.976: 'Hg I', 696.5431: 'Ar I', 763.5106: 'Ar I', 826.4522: 'Ar I',
           823.1635: 'Xe I', 881.9411: 'Xe I', 904.545: 'Xe I', 992.319: 'Xe I',
           1083.837: 'Xe I', 587.0916: 'Kr I', 785.48234: 'Kr I', 850.88728: 'Kr I',
           877.67505: 'Kr I', 892.86934: 'Kr I', 975.1761: 'Kr I'}
# independent NIST ASD air / vacuum values for a few lines (typed from NIST;
# these are the check against the lab reference values)
NIST = {'Hg 546': (546.0735, 546.2268), 'Hg 404': (404.6563, 404.7710),
        'Hg 435': (435.8328, 435.9560), 'Ar 696': (696.5431, 696.7352),
        'Ar 763': (763.5106, 763.7208), 'Ar 826': (826.4522, 826.6794)}

LOG = []


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s)
    LOG.append(s)


# ----------------------------------------------------------------------------
# loaders
# ----------------------------------------------------------------------------
def epochs():
    """All (sn, epoch, dir) in the lab tree."""
    out = []
    for d in sorted(glob.glob(os.path.join(LAB, 'hypstar_*', 'calibration', 'wavelength', '*'))):
        sn = re.search(r'hypstar_(\d+)', d).group(1)
        out.append((sn, os.path.basename(d), d))
    return out


def read_coefs(d):
    f = glob.glob(os.path.join(d, 'hypstar_*_wl_coefs_*.dat'))[0]
    date = re.search(r'coefs_(\d{6})', f).group(1)
    date = pd.Timestamp('20%s-%s-%s' % (date[:2], date[2:4], date[4:6]))
    sig = {}
    for line in open(f):
        m = re.match(r'#\s*sigma_(\w+)\s+([\d.]+)', line)
        if m:
            sig[m.group(1)] = float(m.group(2))
    a = np.loadtxt(f, comments='#')
    return dict(L=a[:, 0], E=a[:, 1], date=date, sigma_L=sig.get('VNIR_L', np.nan),
                sigma_E=sig.get('VNIR_E', np.nan), file=os.path.basename(f))


def read_two_col(d, name):
    a = np.loadtxt(os.path.join(d, name))
    return pd.DataFrame(dict(wl=a[:, 0], L=a[:, 1], E=a[:, 3]))


def read_lamp(d):
    a = np.loadtxt(os.path.join(d, 'lamp_vnir.dat'))
    npair = (a.shape[1] - 1) // 2
    pairs = [(a[:, 1 + 2 * k], a[:, 2 + 2 * k]) for k in range(npair)]
    return a[:, 0], pairs


def poly(c, px):
    return np.polyval(c[::-1], px)


def dpoly(c, px):
    return np.polyval(np.polyder(np.poly1d(c[::-1])), px)


def px_of(c, lam):
    """Pixel of wavelength lam under polynomial c (Newton)."""
    px = (lam - c[0]) / c[1]
    for _ in range(6):
        px = px - (poly(c, px) - lam) / dpoly(c, px)
    return px


EP = {}
for sn, ep, d in epochs():
    EP[(sn, ep)] = dict(dir=d, coefs=read_coefs(d), deltawl=read_two_col(d, 'deltawl_vnir.dat'),
                        fwhm=read_two_col(d, 'fwhm_vnir.dat'), lamp=read_lamp(d))

# ----------------------------------------------------------------------------
# 1. inventory and matching
# ----------------------------------------------------------------------------
say('=' * 78)
say('1. INVENTORY')
tab = load_srf_table()
rows = []
for (sn, ep), e in EP.items():
    c = e['coefs']
    rows.append(dict(sn=sn, epoch=ep, coef_date=c['date'].date(), sigma_L=c['sigma_L'],
                     sigma_E=c['sigma_E'], n_lines=len(e['deltawl']),
                     n_lamp_pairs=len(e['lamp'][1]),
                     lines=' '.join('%.1f' % w for w in e['deltawl']['wl'])))
inv = pd.DataFrame(rows)
inv.to_csv(os.path.join(OUT, 'inventory.csv'), index=False)
say(inv.drop(columns='lines').to_string(index=False))

# match srf_table rows (instrument x cal period) to lab epochs
match = []
for (inst, cp), g in tab.groupby(['instrument', 'cal_period']):
    sn = inst.split('_')[1]
    r = g.iloc[0]
    d_rad, d_irr = pd.Timestamp(r['cal_date_rad']), pd.Timestamp(r['cal_date_irr'])
    labs = [(ep, EP[(s, ep)]['coefs']['date']) for (s, ep) in EP if s == sn and 'tst' not in ep]
    before = [(ep, dt) for ep, dt in labs if dt <= d_rad]
    if before:
        ep, dt = max(before, key=lambda t: t[1])
        lag = (d_rad - dt).days
        status = 'lab precedes cal by %d d' % lag
    elif labs:
        ep, dt = min(labs, key=lambda t: t[1])
        status = 'NO lab epoch before cal; nearest is %s (%d d after)' % (ep, (dt - d_rad).days)
    else:
        ep, dt, status = None, None, 'instrument absent from lab set'
    match.append(dict(instrument=inst, cal_period=cp, cal_date_rad=d_rad.date(),
                      cal_date_irr=d_irr.date(), lab_epoch=ep,
                      lab_date=None if dt is None else dt.date(), status=status,
                      sources=','.join(sorted(set(g['source'])))))
match = pd.DataFrame(match)
match.to_csv(os.path.join(OUT, 'match_cal_periods.csv'), index=False)
say('\nsrf_table cal periods vs lab epochs:')
say(match.to_string(index=False))
used = set(match['lab_epoch'].dropna())
say('\nLab epochs with no srf_table cal period: ' + ', '.join(
    '%s %s (%s)' % (sn, ep, EP[(sn, ep)]['coefs']['date'].date()) for (sn, ep) in EP
    if ep not in used or not ((match['instrument'] == 'HYPSTAR_' + sn) & (match['lab_epoch'] == ep)).any()))
BEST_EPOCH = {(r.instrument.split('_')[1], r.cal_period): r.lab_epoch for r in match.itertuples()}

# ----------------------------------------------------------------------------
# line fitting in pixel space
# ----------------------------------------------------------------------------
def gauss(x, a, mu, sig, b):
    return b + a * np.exp(-0.5 * ((x - mu) / sig) ** 2)


def sgauss(x, a, mu, sig, b, p):
    # super-Gaussian with the same FWHM convention as hypernet.srf.srf_kernel
    f = FWHM_PER_SIGMA * sig
    s_ = f / (2.0 * (2.0 * np.log(2.0)) ** (1.0 / p))
    return b + a * np.exp(-0.5 * np.abs((x - mu) / s_) ** p)


def pvoigt(x, a, mu, sig, b, eta):
    f = FWHM_PER_SIGMA * sig
    gam = 0.5 * f
    g = np.exp(-0.5 * ((x - mu) / sig) ** 2)
    lor = gam ** 2 / ((x - mu) ** 2 + gam ** 2)
    return b + a * ((1 - eta) * g + eta * lor)


def skew_gauss(x, a, mu, sig, b, asym):
    # width differs left/right of the peak: sig*(1 -/+ asym)
    s = np.where(x < mu, sig * (1 - asym), sig * (1 + asym))
    return b + a * np.exp(-0.5 * ((x - mu) / s) ** 2)


def halfmax_fwhm(x, y):
    """FWHM in x units from linear interpolation of the half-maximum crossings
    (background = min of the window)."""
    b = np.min(y)
    yy = y - b
    i = int(np.argmax(yy))
    h = 0.5 * yy[i]
    try:
        jl = np.where(yy[:i] < h)[0][-1]
        xl = x[jl] + (h - yy[jl]) / (yy[jl + 1] - yy[jl]) * (x[jl + 1] - x[jl])
        jr = i + np.where(yy[i:] < h)[0][0]
        xr = x[jr - 1] + (h - yy[jr - 1]) / (yy[jr] - yy[jr - 1]) * (x[jr] - x[jr - 1])
        return xr - xl
    except IndexError:
        return np.nan


def fit_shape(x, y, err, model, p0, bounds):
    try:
        p, cov = curve_fit(model, x, y, p0=p0, sigma=np.full(x.size, err), absolute_sigma=True,
                           bounds=bounds, maxfev=20000)
        chi2 = np.sum(((y - model(x, *p)) / err) ** 2) / (x.size - len(p))
        return p, np.sqrt(np.diag(cov)), chi2
    except Exception:
        return None, None, np.nan


def noise_of(y, wl):
    """Noise from line-free windows (450-530 and 600-680 nm), robust std."""
    m = ((wl > 450) & (wl < 530)) | ((wl > 600) & (wl < 680))
    return 1.4826 * np.median(np.abs(y[m] - np.median(y[m])))


def fit_line_px(y, px0, half=7, noise=1e-3):
    """Fit one emission line around pixel px0; returns dict of all shapes."""
    px = np.arange(NPIX)
    lo, hi = int(round(px0)) - half, int(round(px0)) + half + 1
    x, yy = px[lo:hi].astype(float), y[lo:hi]
    out = dict(px_lo=lo, peak=yy.max(), noise=noise, snr=(yy.max() - yy.min()) / noise)
    # saturation / flat-top diagnostics
    top = np.sort(yy)[::-1]
    out['n_top_98'] = int((yy > 0.98 * yy.max()).sum())
    out['top2_ratio'] = top[1] / top[0]
    out['fwhm_halfmax_px'] = halfmax_fwhm(x, yy)
    b0, a0 = yy.min(), yy.max() - yy.min()
    p, e, chi = fit_shape(x, yy, noise, gauss, [a0, px0, 2.5, b0],
                          ([0, lo, 0.5, -np.inf], [np.inf, hi, 10, np.inf]))
    if p is None:
        return None
    out.update(mu_px=p[1], mu_px_err=e[1], sig_px=abs(p[2]), fwhm_gauss_px=FWHM_PER_SIGMA * abs(p[2]),
               fwhm_gauss_px_err=FWHM_PER_SIGMA * e[2], chi2_gauss=chi, amp=p[0], bkg=p[3])
    # residual sign at the peak (flat top => negative residual at centre)
    out['resid_peak'] = (yy - gauss(x, *p))[np.argmax(yy)] / p[0]
    ps, es, chis = fit_shape(x, yy, noise, sgauss, [p[0], p[1], abs(p[2]), p[3], 2.0],
                             ([0, lo, 0.5, -np.inf, 1.0], [np.inf, hi, 10, np.inf, 8.0]))
    if ps is not None:
        out.update(fwhm_sg_px=FWHM_PER_SIGMA * abs(ps[2]), sg_p=ps[4], sg_p_err=es[4], chi2_sg=chis)
    pv, ev, chiv = fit_shape(x, yy, noise, pvoigt, [p[0], p[1], abs(p[2]), p[3], 0.1],
                             ([0, lo, 0.5, -np.inf, 0.0], [np.inf, hi, 10, np.inf, 1.0]))
    if pv is not None:
        out.update(fwhm_pv_px=FWHM_PER_SIGMA * abs(pv[2]), pv_eta=pv[4], pv_eta_err=ev[4], chi2_pv=chiv)
    pa, ea, chia = fit_shape(x, yy, noise, skew_gauss, [p[0], p[1], abs(p[2]), p[3], 0.0],
                             ([0, lo, 0.5, -np.inf, -0.6], [np.inf, hi, 10, np.inf, 0.6]))
    if pa is not None:
        out.update(asym=pa[4], asym_err=ea[4], chi2_asym=chia, mu_asym_px=pa[1])
    out['x'], out['y'], out['gfit'] = x, yy, gauss(x, *p)
    return out


def lamp_pairs(e, c):
    """Identify the lamp of each L/E pair by its brightest line: the pair
    peaking at Hg 435.83 is the Hg-Ar pen lamp, at Kr 785.48 the Kr lamp, at
    Xe 881.94 the Xe lamp."""
    wl, pairs = e['lamp']
    out = {}
    for k, (L, _) in enumerate(pairs):
        best = None
        for name, lam in (('Hg', 435.8328), ('Kr', 785.48234), ('Xe', 881.9411)):
            i = int(round(px_of(c, lam)))
            v = L[i - 3:i + 4].max() - np.median(L[i - 40:i + 40])
            if best is None or v > best[0]:
                best = (v, name)
        out[k] = best[1]
    return out


def line_pair(e, lam, c):
    """The lamp pair holding the line at lam (by species: Hg and Ar lines come
    from the Hg-Ar pen lamp, Kr and Xe lines from their own lamps) and the
    expected pixel."""
    lamps = lamp_pairs(e, c)
    sp = SPECIES.get(lam, '?').split()[0]
    want = 'Hg' if sp in ('Hg', 'Ar') else sp
    k = [k for k, v in lamps.items() if v == want]
    return (k[0] if k else 0), px_of(c, lam)


# fit every line of every epoch, both channels
say('\n' + '=' * 78)
say('2. FILE SEMANTICS')
fits = []
for (sn, ep), e in EP.items():
    c = e['coefs']
    wl, pairs = e['lamp']
    for _, r in e['deltawl'].iterrows():
        lam = r['wl']
        k, px0 = line_pair(e, lam, c['L'])
        for ch, col in (('L', 0), ('E', 1)):
            y = pairs[k][col]
            if ch == 'E':
                px0 = px_of(c['E'], lam)
            if px0 < 10 or px0 > NPIX - 10:
                continue
            f = fit_line_px(y, px0, noise=noise_of(y, wl))
            if f is None:
                continue
            mu = f['mu_px']
            disp = dpoly(c[ch], mu)
            fits.append(dict(sn=sn, epoch=ep, channel=ch, lam_ref=lam, species=SPECIES.get(lam, '?'),
                             pair=k, px_expected=px0, mu_px=mu, mu_px_err=f['mu_px_err'],
                             lam_fit=poly(c[ch], mu), lam_fit_on_L=poly(c['L'], mu),
                             resid_this=poly(c[ch], mu) - lam, disp_nm_px=disp,
                             fwhm_gauss_nm=f['fwhm_gauss_px'] * disp,
                             fwhm_gauss_nm_err=f['fwhm_gauss_px_err'] * disp,
                             fwhm_halfmax_nm=f['fwhm_halfmax_px'] * disp,
                             fwhm_sg_nm=f.get('fwhm_sg_px', np.nan) * disp, sg_p=f.get('sg_p', np.nan),
                             sg_p_err=f.get('sg_p_err', np.nan),
                             fwhm_pv_nm=f.get('fwhm_pv_px', np.nan) * disp, pv_eta=f.get('pv_eta', np.nan),
                             pv_eta_err=f.get('pv_eta_err', np.nan), asym=f.get('asym', np.nan),
                             asym_err=f.get('asym_err', np.nan),
                             chi2_gauss=f['chi2_gauss'], chi2_sg=f.get('chi2_sg', np.nan),
                             chi2_pv=f.get('chi2_pv', np.nan), chi2_asym=f.get('chi2_asym', np.nan),
                             peak=f['peak'], snr=f['snr'], n_top_98=f['n_top_98'],
                             top2_ratio=f['top2_ratio'], resid_peak=f['resid_peak'],
                             lab_fwhm=float(e['fwhm'].loc[e['fwhm']['wl'] == lam, ch].iloc[0]),
                             lab_deltawl=float(r[ch]), _x=f['x'], _y=f['y'], _g=f['gfit']))
fits = pd.DataFrame(fits)
fits.drop(columns=['_x', '_y', '_g']).to_csv(os.path.join(OUT, 'line_fits.csv'), index=False)

# 2a. the wl column of lamp_vnir.dat
say('\n2a. lamp_vnir.dat "wl" column vs polynomials (rms of wl - poly(px), nm; wl is rounded to 0.1):')
wlrows = []
for (sn, ep), e in EP.items():
    wl = e['lamp'][0]
    c = e['coefs']
    px = np.arange(NPIX)
    d = dict(sn=sn, epoch=ep, rms_L_this=np.sqrt(np.mean((wl - poly(c['L'], px)) ** 2)),
             rms_E_this=np.sqrt(np.mean((wl - poly(c['E'], px)) ** 2)))
    # previous epoch of the same instrument
    prev = [(k, EP[k]) for k in EP if k[0] == sn and EP[k]['coefs']['date'] < c['date']]
    if prev:
        kp, ep_ = max(prev, key=lambda t: t[1]['coefs']['date'])
        d['prev'] = kp[1]
        d['rms_L_prev'] = np.sqrt(np.mean((wl - poly(ep_['coefs']['L'], px)) ** 2))
        d['rms_E_prev'] = np.sqrt(np.mean((wl - poly(ep_['coefs']['E'], px)) ** 2))
    # the wl column as its own polynomial (nominal grid)
    cn = np.polyfit(px, wl, 4)[::-1]
    d['rms_wl_own_poly4'] = np.sqrt(np.mean((wl - poly(cn, px)) ** 2))
    d['nom_c0'], d['nom_c1'] = cn[0], cn[1]
    e['nominal'] = cn
    for lam0 in (400, 550, 700, 850):
        pxx = px_of(c['L'], lam0)
        d['L-nom_%d' % lam0] = poly(c['L'], pxx) - poly(cn, pxx)
        d['E-nom_%d' % lam0] = poly(c['E'], pxx) - poly(cn, pxx)
    wlrows.append(d)
wlrows = pd.DataFrame(wlrows)
wlrows.to_csv(os.path.join(OUT, 'lamp_wl_column.csv'), index=False)
say(wlrows.round(3).to_string(index=False))
say('   wl column identical across the epochs of one instrument: ' + ', '.join(
    '%s %s' % (sn, all(np.array_equal(EP[(sn, ep)]['lamp'][0], EP[(s, e)]['lamp'][0]) for (s, e) in EP if s == sn))
    for sn in sorted(set(s for s, _ in EP)) for ep in [next(e for (s, e) in EP if s == sn)]))
say('   (rounding to 0.1 nm alone gives rms 0.029; the wl column is a 4th-order polynomial, the same at every epoch')
say('    of an instrument, and is neither the L nor the E polynomial of the epoch: a nominal per-instrument grid)')

# 2b. lamp pairs: which lines peak in which pair
say('\n2b. lamp pairs: lines (species) found in each L/E pair:')
pr = fits[fits.channel == 'L'].groupby(['sn', 'epoch', 'pair']).apply(
    lambda g: ' '.join('%s%.0f' % (s.split()[0], l) for s, l in zip(g.species, g.lam_ref))).reset_index(name='lines')
say(pr.to_string(index=False))
# are L and E of a pair the same lamp (same lines) and are they normalised?
say('column maxima (all 1 => each column normalised to its peak):',
    sorted(set(np.round([p[j].max() for e in EP.values() for p in e['lamp'][1] for j in (0, 1)], 3))))

# 2c. deltawl hypotheses
say('\n2c. deltawl hypotheses (rms over lines of [candidate - deltawl], nm):')
hyp = []
for (sn, ep), e in EP.items():
    c = e['coefs']
    prev = [(k, EP[k]) for k in EP if k[0] == sn and EP[k]['coefs']['date'] < c['date']]
    cp = max(prev, key=lambda t: t[1]['coefs']['date'])[1]['coefs'] if prev else None
    for ch in ('L', 'E'):
        g = fits[(fits.sn == sn) & (fits.epoch == ep) & (fits.channel == ch) & (fits.snr > 10)]
        dw = g['lab_deltawl'].values
        mu = g['mu_px'].values
        lam = g['lam_ref'].values
        cn = e['nominal']
        cand = {'(a) poly_this(mu)-ref': poly(c[ch], mu) - lam,
                '(a-) ref-poly_this(mu)': lam - poly(c[ch], mu),
                '(c) poly_L(mu)-ref': poly(c['L'], mu) - lam,
                '(c) poly_E(mu)-ref': poly(c['E'], mu) - lam,
                '(c) poly_E-poly_L at mu': poly(c['E'], mu) - poly(c['L'], mu),
                '(d) poly_this-nominal at mu': poly(c[ch], mu) - poly(cn, mu),
                '(d2) ref-nominal(mu)': lam - poly(cn, mu)}
        if cp is not None:
            cand['(b) poly_prev(mu)-ref'] = poly(cp[ch], mu) - lam
            cand['(b) poly_prev-poly_this at mu'] = poly(cp[ch], mu) - poly(c[ch], mu)
            cand['(b-) poly_this-poly_prev at mu'] = poly(c[ch], mu) - poly(cp[ch], mu)
        d = dict(sn=sn, epoch=ep, channel=ch, n=len(g), deltawl_mean=dw.mean(), deltawl_rms=np.sqrt(np.mean(dw ** 2)))
        for k, v in cand.items():
            d[k] = np.sqrt(np.mean((v - dw) ** 2))
        hyp.append(d)
hyp = pd.DataFrame(hyp)
hyp.to_csv(os.path.join(OUT, 'deltawl_hypotheses.csv'), index=False)
pd.set_option('display.width', 250)
say(hyp.round(3).to_string(index=False))
# stated sigma vs rms of our residuals and of deltawl
say('\n   stated sigma vs rms(deltawl) vs rms(our Gaussian-centre residual of this poly):')
for (sn, ep), e in EP.items():
    for ch in ('L', 'E'):
        g = fits[(fits.sn == sn) & (fits.epoch == ep) & (fits.channel == ch)]
        g = g[(g.lam_ref > 380) & (g.lam_ref < 1100)]
        say('   %s %-12s %s sigma=%.2f rms(deltawl)=%.3f mean(deltawl)=%+.3f rms(resid_this)=%.3f mean(resid_this)=%+.3f  n=%d'
            % (sn, ep, ch, e['coefs']['sigma_' + ch], np.sqrt(np.mean(g.lab_deltawl ** 2)), g.lab_deltawl.mean(),
               np.sqrt(np.mean(g.resid_this ** 2)), g.resid_this.mean(), len(g)))

# 2d. air vs vacuum
say('\n2d. air vs vacuum: lab reference values vs NIST air / vacuum (nm):')
all_refs = sorted(set(np.concatenate([e['deltawl']['wl'].values for e in EP.values()])))
for name, (air, vac) in NIST.items():
    near = min(all_refs, key=lambda w: abs(w - air))
    say('   %-7s lab ref %.4f  NIST air %.4f (d=%+.4f)  NIST vac %.4f (d=%+.4f)  air_to_vac(lab)=%.4f'
        % (name, near, air, near - air, vac, near - vac, srfmod.air_to_vac(near)))
say('   all lab reference lines (species inferred by value): ' +
    ', '.join('%.4f %s' % (w, SPECIES.get(w, '?')) for w in all_refs))

# ----------------------------------------------------------------------------
# 3. wavelength grids
# ----------------------------------------------------------------------------
say('\n' + '=' * 78)
say('3. WAVELENGTH GRIDS')


def best_crop(c, target):
    lam = poly(c, np.arange(NPIX))
    n = len(target)
    best = (np.inf, None)
    for s in range(NPIX - n + 1):
        r = np.sqrt(np.mean((lam[s:s + n] - target) ** 2))
        if r < best[0]:
            best = (r, s)
    r, s = best
    seg = lam[s:s + n]
    return dict(start=s, end=s + n - 1, rms=r, maxabs=np.max(np.abs(seg - target)),
                lam_first=seg[0], lam_before=lam[s - 1] if s > 0 else np.nan,
                lam_last=seg[-1], lam_after=lam[s + n] if s + n < NPIX else np.nan)


targets = []
irr = load_l1a_irr(glob.glob(VEIT_DIR + '*L1A_IRR*.nc')[0])
rad = load_l1a_rad(glob.glob(VEIT_DIR + '*L1A_RAD*.nc')[0])
targets.append(('122304', 'VEIT L1A E (cal irr %s)' % irr['meta']['instrument_calibration_date_irr'], 'E', irr['wave']))
targets.append(('122304', 'VEIT L1A L', 'L', rad['wave']))
req = pd.read_csv(os.path.join(REPO, 'docs', 'wiggles_data_request.csv'))
R2GRID = {}
for (inst, cal), g in req.groupby(['instrument', 'cal_dates_rad_irr']):
    for _, row in g.iterrows():
        p = release2_path(row['site'], row['sequence_time'], row['file'])
        if os.path.exists(p):
            w = load_l2b(p)['wave']
            R2GRID[(inst.split('_')[1], cal)] = w
            targets.append((inst.split('_')[1], 'Release 2 L2B L, cal %s (%s)' % (cal, row['site']), 'L', w))
            break
crop = []
for sn, label, ch, w in targets:
    for (s2, ep), e in EP.items():
        if s2 != sn:
            continue
        for cch in ('L', 'E'):
            b = best_crop(e['coefs'][cch], w)
            crop.append(dict(sn=sn, target=label, n=len(w), epoch=ep, poly=cch, **b))
crop = pd.DataFrame(crop)
crop.to_csv(os.path.join(OUT, 'grid_crops.csv'), index=False)
say('best crop of each epoch polynomial reproducing each target grid (rms, max |d| in nm):')
say(crop[crop.poly == crop.target.str.extract(r'L1A (E|L)|L2B (L)').bfill(axis=1).iloc[:, 0]]
    .round(4).to_string(index=False))
say('\n(cross-channel check: E poly against an L grid, best rms per target)')
say(crop[crop.poly == 'E'][crop.target.str.contains('L2B|L1A L')].groupby('target').rms.min().round(3).to_string())

# the crop rule: first/last wavelengths of the matched crops
good = crop[(crop.rms < 0.01)]
say('\nmatched crops (rms < 0.01 nm): start px, lam[start-1] < 350 <= lam[start]? lam[end] <= 1100 < lam[end+1]?')
say(good[['sn', 'target', 'epoch', 'poly', 'start', 'end', 'n', 'lam_before', 'lam_first', 'lam_last', 'lam_after', 'rms']]
    .round(4).to_string(index=False))


def crop_rule(c, lo=350.0, hi=1100.0, select=None):
    """Pixels with lo <= lam < hi, where lam is the polynomial ``c`` (or the
    polynomial ``select`` when given, e.g. the nominal grid), and the
    calibrated wavelengths of ``c`` on them."""
    px = np.arange(NPIX)
    lam = poly(c, px)
    sel = lam if select is None else poly(select, px)
    m = np.where((sel >= lo) & (sel < hi))[0]
    return m[0], m[-1], lam[m]


# which selection rule reproduces the matched crops: the calibrated polynomial
# itself, or the instrument's nominal grid (the wl column of lamp_vnir.dat)?
say('\ncrop rule test: observed (start, n) vs 350 <= lam < 1100 on the calibrated poly, and on the nominal grid:')
for r in good.itertuples():
    c = EP[(r.sn, r.epoch)]['coefs']
    s1, e1, _ = crop_rule(c[r.poly])
    s2, e2, _ = crop_rule(c[r.poly], select=EP[(r.sn, r.epoch)]['nominal'])
    say('   %s %-52s obs %d,%d  calpoly %d,%d %s  nominal %d,%d %s' % (
        r.sn, r.target[:52], r.start, r.n, s1, e1 - s1 + 1, 'OK' if (s1, e1 - s1 + 1) == (r.start, r.n) else 'X',
        s2, e2 - s2 + 1, 'OK' if (s2, e2 - s2 + 1) == (r.start, r.n) else 'X'))


# own E grids for every instrument, from the epoch matched to each cal period
say('\nown E grids (crop rule 350 <= lam < 1100 applied to the E polynomial of the matched epoch):')
egrids = {}
erows = []
proxy_c = EP[('122304', '2024_09')]['coefs']['E']
for r in match.itertuples():
    if r.lab_epoch is None:
        continue
    sn = r.instrument.split('_')[1]
    c = EP[(sn, r.lab_epoch)]['coefs']
    s_e, e_e, lamE = crop_rule(c['E'])
    s_l, e_l, lamL = crop_rule(c['L'])
    egrids[(sn, r.cal_period)] = lamE
    px = np.arange(NPIX)
    pxE = px[s_e:e_e + 1]
    # proxy (122304 2024_09 E poly, same crop rule) vs own E, by pixel index and by wavelength
    s_p, e_p, lamP = crop_rule(proxy_c)
    d = dict(instrument=r.instrument, cal_period=r.cal_period, lab_epoch=r.lab_epoch, nE=len(lamE),
             E_start=s_e, nL=len(lamL), L_start=s_l, E_first=lamE[0], E_last=lamE[-1])
    for lam0 in (400, 550, 700, 850):
        pxe = px_of(c['E'], lam0)
        d['E-L_samepx_%d' % lam0] = poly(c['E'], pxe) - poly(c['L'], pxe)   # E minus L at the E pixel
        d['proxy-own_samepx_%d' % lam0] = poly(proxy_c, pxe) - poly(c['E'], pxe)
        # proxy minus own E, compared at the same *grid index* (what the resampling model sees)
        j = int(round(pxe)) - s_e
        d['proxy-own_sameidx_%d' % lam0] = (lamP[j] - lamE[j]) if 0 <= j < min(len(lamP), len(lamE)) else np.nan
        # phase of the E pixels relative to the L pixels, as a fraction of the E pitch:
        # the distance from the L pixel at lam0 to the nearest own / proxy E pixel
        pxl = int(round(px_of(c['L'], lam0)))
        laml = poly(c['L'], pxl)
        pitch = dpoly(c['E'], pxe)
        d['E_phase_own_%d' % lam0] = (laml - lamE[np.argmin(np.abs(lamE - laml))]) / pitch
        d['E_phase_proxy_%d' % lam0] = (laml - lamP[np.argmin(np.abs(lamP - laml))]) / pitch
    # L2 R2 grid vs own L poly
    if (sn, r.cal_period) in R2GRID:
        w = R2GRID[(sn, r.cal_period)]
        b = best_crop(c['L'], w)
        d['R2_L_rms'] = b['rms']
        d['R2_L_start'] = b['start']
        d['R2_nL'] = len(w)
    erows.append(d)
    pd.DataFrame(dict(px=pxE, wavelength_E=lamE)).to_csv(
        os.path.join(OUT, 'egrid_%s_%s.csv' % (sn, r.cal_period.replace('/', '_'))), index=False)
erows = pd.DataFrame(erows)
erows.to_csv(os.path.join(OUT, 'own_e_grids.csv'), index=False)
say(erows.round(3).to_string(index=False))

# figure: E - L (same pixel) and proxy - own vs wavelength
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for i, r in enumerate(match.itertuples()):
    if r.lab_epoch is None:
        continue
    sn = r.instrument.split('_')[1]
    c = EP[(sn, r.lab_epoch)]['coefs']
    px = np.arange(300, 2000)
    lamE = poly(c['E'], px)
    ax[0].plot(lamE, lamE - poly(c['L'], px), color=COL[i % 8], lw=1.5, label='%s %s' % (sn, r.lab_epoch))
    ax[1].plot(lamE, poly(proxy_c, px) - lamE, color=COL[i % 8], lw=1.5)
ax[0].set_xlabel('wavelength (nm)'); ax[0].set_ylabel('lam_E - lam_L at the same pixel (nm)')
ax[1].set_xlabel('wavelength (nm)'); ax[1].set_ylabel('proxy (122304 2024_09 E) - own E, same pixel (nm)')
for a in ax:
    a.set_xlim(350, 1100); a.grid(alpha=0.3)
ax[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'grid_offsets.png'), dpi=130); plt.close(fig)

# ----------------------------------------------------------------------------
# 4. FWHM comparison with srf_table
# ----------------------------------------------------------------------------
say('\n' + '=' * 78)
say('4. FWHM: lab vs srf_table')
cmp_rows = []
for r in tab.itertuples():
    sn = r.instrument.split('_')[1]
    ep = BEST_EPOCH.get((sn, r.cal_period))
    if ep is None:
        continue
    e = EP[(sn, ep)]
    ch = 'L' if r.channel == 'Ld' else 'E'
    coeffs = np.array([r.c0, r.c1, r.c2])
    for _, l in e['fwhm'].iterrows():
        lam = l['wl']
        if not (395 <= lam <= 865):
            continue
        model = srfmod.fwhm_at(lam, coeffs, r.lam_ref, r.lam_scale)
        f = fits[(fits.sn == sn) & (fits.epoch == ep) & (fits.channel == ch) & (fits.lam_ref == lam)]
        cmp_rows.append(dict(instrument=r.instrument, channel=r.channel, cal_period=r.cal_period, source=r.source,
                             lab_epoch=ep, lam=lam, species=SPECIES.get(lam, '?'), lab_fwhm=l[ch],
                             model_fwhm=model, lab_minus_model=l[ch] - model,
                             our_gauss_fwhm=float(f.fwhm_gauss_nm.iloc[0]) if len(f) else np.nan,
                             our_halfmax_fwhm=float(f.fwhm_halfmax_nm.iloc[0]) if len(f) else np.nan,
                             snr=float(f.snr.iloc[0]) if len(f) else np.nan))
cmp = pd.DataFrame(cmp_rows)
cmp['weak'] = cmp.snr < 10
cmp.to_csv(os.path.join(OUT, 'fwhm_compare.csv'), index=False)
say('per line (lab FWHM from fwhm_vnir.dat; model = srf_table FWHM(lambda); our = Gaussian refit of lamp line; weak = S/N < 10):')
say(cmp.round(3).to_string(index=False))

say('\nsummary per srf_table row (lines 395-865 nm with S/N >= 10, i.e. without Hg 297 / Ar 697 / Ar 826):')
summ = cmp[~cmp.weak].groupby(['instrument', 'channel', 'cal_period', 'source', 'lab_epoch']).agg(
    n=('lam', 'size'), lab_med=('lab_fwhm', 'median'), model_med=('model_fwhm', 'median'),
    lab_minus_model_med=('lab_minus_model', 'median'), lab_minus_model_min=('lab_minus_model', 'min'),
    lab_minus_model_max=('lab_minus_model', 'max')).reset_index()
cmp['ours_minus_model'] = cmp.our_gauss_fwhm - cmp.model_fwhm
summ['ourgauss_minus_model_med'] = cmp[~cmp.weak].groupby(
    ['instrument', 'channel', 'cal_period', 'source', 'lab_epoch']).ours_minus_model.median().values
# without the two lines whose lab width is raised by a neighbour (Hg 404.66, Xe 823.16)
cl2 = cmp[~cmp.weak & ~cmp.lam.isin([404.6563, 823.1635])]
summ['lab_minus_model_med_noblend'] = cl2.groupby(
    ['instrument', 'channel', 'cal_period', 'source', 'lab_epoch']).lab_minus_model.median().values
summ['n_noblend'] = cl2.groupby(['instrument', 'channel', 'cal_period', 'source', 'lab_epoch']).lam.size().values
say(summ.round(3).to_string(index=False))

# lab delta FWHM (L - E), median over lines 400-680, vs our delta_fwhm
say('\nlab dFWHM (L - E) at each epoch: median over lines in 400-680 nm (and 400-865), vs srf_table delta_fwhm:')
drows = []
for (sn, ep), e in EP.items():
    f = e['fwhm']
    m = (f.wl >= 400) & (f.wl <= 680)
    m2 = (f.wl >= 400) & (f.wl <= 865)
    inst = 'HYPSTAR_' + sn
    ours = tab[(tab.instrument == inst) & (tab.channel == 'E')]
    ours = ours[[BEST_EPOCH.get((sn, cp)) == ep for cp in ours.cal_period]]
    d = dict(sn=sn, epoch=ep, n_400_680=int(m.sum()), lab_dfwhm_400_680=float((f.L - f.E)[m].median()),
             lab_dfwhm_400_865=float((f.L - f.E)[m2].median()),
             lab_lines_400_680=' '.join('%.0f:%+.2f' % (w, d_) for w, d_ in zip(f.wl[m], (f.L - f.E)[m])),
             our_delta_fwhm=' / '.join('%.3f (%s %s)' % (x, cp, s) for x, cp, s in zip(ours.delta_fwhm, ours.cal_period, ours.source)) or '-',
             our_method=' / '.join(ours.method) or '-')
    d['lab_method_G1'] = 'srf' if d['lab_dfwhm_400_680'] >= 0.15 else 'linear'
    drows.append(d)
drows = pd.DataFrame(drows)
drows.to_csv(os.path.join(OUT, 'delta_fwhm_compare.csv'), index=False)
say(drows.to_string(index=False))

# FWHM definition check: lab fwhm vs our Gaussian / half-max refits of the same lines
say('\nFWHM definition: lab value minus our refit of the same lamp line (median over lines 395-865 with S/N >= 10, nm):')
g = fits[(fits.lam_ref >= 395) & (fits.lam_ref <= 865) & (fits.snr >= 10)]
say('   per strong line (S/N >= 50), median over epochs of lab - Gaussian / lab - half-max / lab - super-Gaussian:')
for (lam, ch), gg in g[g.snr >= 50].groupby(['lam_ref', 'channel']):
    say('     %8.3f %s n=%2d  lab %.3f  gauss %+.3f  halfmax %+.3f  sgauss %+.3f' % (
        lam, ch, len(gg), gg.lab_fwhm.median(), (gg.lab_fwhm - gg.fwhm_gauss_nm).median(),
        (gg.lab_fwhm - gg.fwhm_halfmax_nm).median(), (gg.lab_fwhm - gg.fwhm_sg_nm).median()))
say('   lab - Gaussian fit: median %+.3f, lab - half-max crossing: median %+.3f, lab - super-Gaussian: %+.3f, lab - pseudo-Voigt: %+.3f'
    % ((g.lab_fwhm - g.fwhm_gauss_nm).median(), (g.lab_fwhm - g.fwhm_halfmax_nm).median(),
       (g.lab_fwhm - g.fwhm_sg_nm).median(), (g.lab_fwhm - g.fwhm_pv_nm).median()))
say('   rms: Gaussian %.3f, half-max %.3f, super-Gaussian %.3f, pseudo-Voigt %.3f'
    % tuple(np.sqrt(np.mean((g.lab_fwhm - g[k]) ** 2)) for k in ('fwhm_gauss_nm', 'fwhm_halfmax_nm', 'fwhm_sg_nm', 'fwhm_pv_nm')))

# figure: lab FWHM vs model per instrument
insts = sorted(set(cmp.instrument))
fig, axs = plt.subplots(2, 4, figsize=(15, 7), sharex=True, sharey=True)
for a, inst in zip(axs.ravel(), insts):
    g = cmp[cmp.instrument == inst]
    for ch, col, mk in (('E', COL[0], 'o'), ('Ld', COL[1], 's')):
        gg = g[g.channel == ch]
        a.scatter(gg.lam, gg.lab_fwhm, s=28, color=col, marker=mk, label='lab %s' % ch, zorder=3)
        for (cp, src), g3 in gg.groupby(['cal_period', 'source']):
            r = tab[(tab.instrument == inst) & (tab.channel == ch) & (tab.cal_period == cp) & (tab.source == src)].iloc[0]
            lam = np.linspace(395, 865, 100)
            a.plot(lam, srfmod.fwhm_at(lam, [r.c0, r.c1, r.c2], r.lam_ref, r.lam_scale), color=col, lw=1.5,
                   ls='-' if src == 'release2_l2b' else '--', alpha=0.8,
                   label='table %s %s%s' % (ch, cp[:7], '' if src == 'release2_l2b' else ' L1A'))
    a.set_title(inst, fontsize=10); a.grid(alpha=0.3); a.legend(fontsize=6)
for a in axs[1]:
    a.set_xlabel('wavelength (nm)')
for a in axs[:, 0]:
    a.set_ylabel('FWHM (nm)')
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'fwhm_lab_vs_table.png'), dpi=130); plt.close(fig)

# ----------------------------------------------------------------------------
# 5. line shapes
# ----------------------------------------------------------------------------
say('\n' + '=' * 78)
say('5. LINE SHAPES (lamp lines fitted in pixel space; FWHM in nm via the local dispersion)')
say('saturation / S/N per line (L channel, all epochs): n_top_98 = pixels within 2% of the peak; top2 = 2nd/1st pixel')
sat = fits.groupby(['lam_ref', 'species', 'channel']).agg(peak_med=('peak', 'median'), snr_med=('snr', 'median'),
                                                           n_top98_max=('n_top_98', 'max'), top2_med=('top2_ratio', 'median'),
                                                           resid_peak_med=('resid_peak', 'median'),
                                                           sg_p_med=('sg_p', 'median'), chi2_g=('chi2_gauss', 'median'),
                                                           chi2_sg=('chi2_sg', 'median'), chi2_pv=('chi2_pv', 'median'),
                                                           chi2_asym=('chi2_asym', 'median'),
                                                           pv_eta_med=('pv_eta', 'median'), asym_med=('asym', 'median'),
                                                           fwhm_med=('fwhm_gauss_nm', 'median'), lab_fwhm_med=('lab_fwhm', 'median')).reset_index()
sat.to_csv(os.path.join(OUT, 'line_shape_summary.csv'), index=False)
say(sat.round(3).to_string(index=False))

CLEAN = [404.6563, 435.8328, 546.0735, 763.5106, 826.4522, 881.9411]
say('\nclean lines, per instrument/epoch/channel: p (super-Gaussian), eta (pseudo-Voigt), asym, chi2_nu Gaussian -> best:')
cl = fits[fits.lam_ref.isin(CLEAN)].copy()
cl['best'] = cl[['chi2_gauss', 'chi2_sg', 'chi2_pv', 'chi2_asym']].idxmin(axis=1).str.replace('chi2_', '')
agg = cl.groupby(['sn', 'epoch', 'channel']).agg(n=('lam_ref', 'size'), p=('sg_p', 'median'), p_min=('sg_p', 'min'),
                                                p_max=('sg_p', 'max'), eta=('pv_eta', 'median'), asym=('asym', 'median'),
                                                chi2_g=('chi2_gauss', 'median'), chi2_sg=('chi2_sg', 'median'),
                                                chi2_pv=('chi2_pv', 'median'), chi2_asym=('chi2_asym', 'median'),
                                                fwhm_gauss_minus_lab=('fwhm_gauss_nm', lambda s: np.nan)).reset_index()
agg['fwhm_gauss_minus_lab'] = cl.groupby(['sn', 'epoch', 'channel']).apply(lambda g: (g.fwhm_gauss_nm - g.lab_fwhm).median()).values
agg['best_counts'] = cl.groupby(['sn', 'epoch', 'channel']).best.agg(
    lambda s: ' '.join('%s:%d' % kv for kv in s.value_counts().items())).values
agg.to_csv(os.path.join(OUT, 'line_shapes_clean.csv'), index=False)
say(agg.round(3).to_string(index=False))
say('\nclean-line p by species x channel (median, [min,max]):')
for (sp, ch), g in cl.groupby(['species', 'channel']):
    say('   %-5s %s: p=%.2f [%.2f,%.2f]  eta=%.2f  asym=%+.3f  chi2_g/sg/pv/asym=%.1f/%.1f/%.1f/%.1f  n=%d'
        % (sp, ch, g.sg_p.median(), g.sg_p.min(), g.sg_p.max(), g.pv_eta.median(), g.asym.median(),
           g.chi2_gauss.median(), g.chi2_sg.median(), g.chi2_pv.median(), g.chi2_asym.median(), len(g)))

# figure: profiles of clean lines for 122304 2024_09 and 121222 2021_07, L and E, with residuals
fig, axs = plt.subplots(2, 6, figsize=(17, 6))
for row, key in enumerate([('122304', '2024_09'), ('121222', '2021_07')]):
    for j, lam in enumerate(CLEAN):
        a = axs[row, j]
        for ch, col in (('L', COL[0]), ('E', COL[1])):
            f = fits[(fits.sn == key[0]) & (fits.epoch == key[1]) & (fits.channel == ch) & (fits.lam_ref == lam)]
            if not len(f):
                continue
            f = f.iloc[0]
            x = (f['_x'] - f['mu_px']) * f['disp_nm_px']
            yn = (f['_y'] - f['_y'].min()) / (f['_y'].max() - f['_y'].min())
            a.semilogy(x, np.clip(yn, 1e-4, None), 'o-', ms=3, lw=1, color=col,
                       label='%s FWHM %.2f (lab %.2f) p=%.2f' % (ch, f['fwhm_gauss_nm'], f['lab_fwhm'], f['sg_p']))
            gn = (f['_g'] - f['_y'].min()) / (f['_y'].max() - f['_y'].min())
            a.semilogy(x, np.clip(gn, 1e-4, None), '--', lw=1, color=col, alpha=0.6)
        a.set_title('%s %s  %s %.1f' % (key[0], key[1], SPECIES.get(lam, ''), lam), fontsize=8)
        a.set_ylim(1e-4, 2); a.legend(fontsize=6); a.grid(alpha=0.3)
        a.set_xlabel('nm from centre')
fig.suptitle('lamp line profiles (points) and Gaussian fits (dashed), normalised, log scale', fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'line_profiles.png'), dpi=130); plt.close(fig)

# ----------------------------------------------------------------------------
# 6. wings / stray light and the common window pattern
# ----------------------------------------------------------------------------
say('\n' + '=' * 78)
say('6. WINGS FAR FROM LINES, AND THE WINDOW PATTERN')
say('Hg-Ar pair: level relative to the 546.07 peak at offsets from the line (median of L and E columns over epochs):')
offs = [2, 3, 5, 10, 20, 50]
wrows = []
for (sn, ep), e in EP.items():
    wl, pairs = e['lamp']
    c = e['coefs']
    for ch, col in (('L', 0), ('E', 1)):
        y = pairs[0][col]
        p546 = px_of(c[ch], 546.0735)
        i = int(round(p546))
        pk = y[i - 3:i + 4].max()
        d = dict(sn=sn, epoch=ep, channel=ch, peak=pk)
        for o in offs:
            dpx = o / dpoly(c[ch], p546)
            for sgn, lab in ((-1, 'blue'), (1, 'red')):
                j = int(round(p546 + sgn * dpx))
                d['%s_%dnm' % (lab, o)] = np.median(y[j - 1:j + 2]) / pk
        # broad background between lines (560-690 nm, no lines) relative to peak, and its noise
        m = (wl > 560) & (wl < 690)
        d['bkg_560_690'] = np.median(y[m]) / pk
        d['noise_560_690'] = noise_of(y, wl) / pk
        wrows.append(d)
wrows = pd.DataFrame(wrows)
wrows.to_csv(os.path.join(OUT, 'wings_546.csv'), index=False)
say(wrows.groupby('channel').median(numeric_only=True).T.round(5).to_string())
say('   (a Gaussian of FWHM 2.7 nm is at 1e-4 of the peak by 3.5 nm; 1e-9 at 5 nm)')
say('Xe pair (lower noise floor): level relative to the Xe 881.94 peak, median over epochs (L and E):')
xrows = []
for (sn, ep), e in EP.items():
    wl, pairs = e['lamp']
    c = e['coefs']
    kx = [k for k, v in lamp_pairs(e, c['L']).items() if v == 'Xe']
    if not kx:
        continue
    for ch, col in (('L', 0), ('E', 1)):
        y = pairs[kx[0]][col]
        p0 = px_of(c[ch], 881.9411)
        i = int(round(p0))
        pk = y[i - 3:i + 4].max()
        d = dict(sn=sn, epoch=ep, channel=ch, peak=pk)
        for o in (3, 5, 10, 20, 50):
            dpx = o / dpoly(c[ch], p0)
            for sgn, lab in ((-1, 'blue'), (1, 'red')):
                j = int(round(p0 + sgn * dpx))
                d['%s_%dnm' % (lab, o)] = np.median(y[j - 1:j + 2]) / pk
        m = (wl > 600) & (wl < 800)   # no Xe lines
        d['bkg_600_800'] = np.median(y[m]) / pk
        d['noise_600_800'] = noise_of(y, wl) / pk
        xrows.append(d)
xrows = pd.DataFrame(xrows)
xrows.to_csv(os.path.join(OUT, 'wings_xe882.csv'), index=False)
say(xrows.groupby('channel').median(numeric_only=True).T.round(5).to_string())

# a common component in the polynomial residuals at the lines (line x epoch matrix)
say('\nPolynomial residual at each line (our Gaussian centre through this epoch\'s poly, minus ref; nm):')
say('  mean over instruments and its sem, vs the scatter: a common component would show as mean >> sem.')
for ch in ('L', 'E'):
    g = fits[(fits.channel == ch) & (fits.snr >= 10) & (~fits.epoch.str.contains('tst'))]
    # one epoch per instrument (the one matched to our cal periods, or the latest)
    keep = []
    for sn in sorted(set(g.sn)):
        eps_ = sorted(set(g[g.sn == sn].epoch))
        keep.append((sn, eps_[-1]))
    g = g[[(a, b) in keep for a, b in zip(g.sn, g.epoch)]]
    piv = g.pivot_table(index='lam_ref', columns='sn', values='resid_this')
    piv['mean'] = piv.mean(axis=1)
    piv['sem'] = piv.std(axis=1) / np.sqrt(piv.notna().sum(axis=1))
    piv['scatter'] = piv.drop(columns=['mean', 'sem']).std(axis=1)
    piv.round(3).to_csv(os.path.join(OUT, 'poly_resid_matrix_%s.csv' % ch))
    say('  channel %s (latest epoch per instrument):' % ch)
    say(piv.round(3).to_string())

# figure: full Hg-Ar lamp spectra, log scale, for all instruments (L)
fig, ax = plt.subplots(figsize=(12, 4.5))
for i, ((sn, ep), e) in enumerate(EP.items()):
    if 'tst' in ep:
        continue
    wl, pairs = e['lamp']
    ax.semilogy(wl, np.clip(pairs[0][0], 1e-5, None), lw=0.8, color=COL[i % 8], alpha=0.8, label='%s %s' % (sn, ep))
for lam0 in (415, 505, 565, 605, 665):
    ax.axvline(lam0, color='k', lw=0.6, ls=':')
ax.set_xlim(380, 900); ax.set_ylim(1e-5, 2); ax.set_xlabel('wavelength (nm)'); ax.set_ylabel('Hg-Ar lamp, L column (normalised)')
ax.legend(fontsize=6, ncol=4); ax.grid(alpha=0.3)
ax.set_title('Hg-Ar lamp spectra, L channel; dotted: windows of the 0.13 nm common pattern', fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'lamp_spectra_log.png'), dpi=130); plt.close(fig)

# polynomial-gap test: refit the wavelength polynomial with and without the Kr lines
# (121222, 121231) and compare at the window-pattern wavelengths.
say('\nPolynomial gap test (121222, 121231 have Kr lines at 587/785/851/878/893/975 nm):')
say('  4th-order poly fitted to our Gaussian centres with all lines vs the standard Hg/Ar/Xe set only;')
say('  difference at the window-pattern wavelengths and the residual of the all-line fit at the Kr lines (nm):')
pat = pd.read_csv(os.path.join(REPO, 'hypernet', 'data', 'wavescale_window_pattern.csv'))
PATLAM = [415, 505, 565, 605, 665]
gap = []
for (sn, ep), e in EP.items():
    for ch in ('L', 'E'):
        g = fits[(fits.sn == sn) & (fits.epoch == ep) & (fits.channel == ch)]
        g = g[(g.species != '?')]
        if not (g.species == 'Kr I').any():
            continue
        std = g[g.species != 'Kr I']
        c_all = np.polyfit(g.mu_px, g.lam_ref, 4)
        c_std = np.polyfit(std.mu_px, std.lam_ref, 4)
        kr = g[g.species == 'Kr I']
        d = dict(sn=sn, epoch=ep, channel=ch, n_all=len(g), n_std=len(std))
        for lam0 in PATLAM:
            pxx = px_of(e['coefs'][ch], lam0)
            d['std-all_%d' % lam0] = np.polyval(c_std, pxx) - np.polyval(c_all, pxx)
            d['file-all_%d' % lam0] = poly(e['coefs'][ch], pxx) - np.polyval(c_all, pxx)
        d['std_resid_at_Kr'] = ' '.join('%.0f:%+.3f' % (l, np.polyval(c_std, m) - l) for l, m in zip(kr.lam_ref, kr.mu_px))
        d['all_resid_at_Kr'] = ' '.join('%.0f:%+.3f' % (l, np.polyval(c_all, m) - l) for l, m in zip(kr.lam_ref, kr.mu_px))
        d['file_resid_at_Kr'] = ' '.join('%.0f:%+.3f' % (l, poly(e['coefs'][ch], m) - l) for l, m in zip(kr.lam_ref, kr.mu_px))
        gap.append(d)
gap = pd.DataFrame(gap)
gap.to_csv(os.path.join(OUT, 'poly_gap_test.csv'), index=False)
say(gap.drop(columns=['std_resid_at_Kr', 'all_resid_at_Kr', 'file_resid_at_Kr']).round(3).to_string(index=False))
say('  residual at the Kr lines of: the standard-set fit (what a 13-line instrument\'s poly would do) / the all-line fit / the file poly:')
for r in gap.itertuples():
    say('   %s %s %s  std: %s | all: %s | file: %s' % (r.sn, r.epoch, r.channel, r.std_resid_at_Kr, r.all_resid_at_Kr, r.file_resid_at_Kr))
say('\nCommon lab residual at the lines vs the window pattern (pattern is in -dlam, dlam = measured - reference;')
say('  the lab residual poly(mu) - ref is a dlam, so compare the pattern with MINUS the mean residual):')
common = fits[(fits.snr >= 10) & (~fits.epoch.str.contains('tst'))].groupby('lam_ref').resid_this.agg(['mean', 'std', 'size'])
for lam, r in common.iterrows():
    near = pat.iloc[(pat.lam_center - lam).abs().argmin()]
    if abs(near.lam_center - lam) <= 5 and 395 <= lam <= 865:
        say('   line %.2f: -mean resid %+.3f (std %.3f, n=%d)   window %.0f pattern %+.3f (scatter %.3f)' % (
            lam, -r['mean'], r['std'], r['size'], near.lam_center, near.pattern, near.scatter))
common.round(3).to_csv(os.path.join(OUT, 'poly_resid_common.csv'))
say('  window pattern (hypernet/data/wavescale_window_pattern.csv) at those windows: ' +
    ', '.join('%d:%+.3f' % (l, pat.loc[pat.lam_center == l, 'pattern'].iloc[0]) for l in PATLAM))
say('  pattern rms over all windows: %.3f nm' % np.sqrt(np.mean(pat.pattern ** 2)))
# where are the lamp lines relative to the pattern windows?  the only lines in 440-690 nm
say('  lab lines in 440-690 nm: ' + ', '.join('%.1f' % w for w in all_refs if 440 < w < 690) +
    '  (the standard set has only 546.07 between 435.83 and 696.54)')

# ----------------------------------------------------------------------------
# 7. stability: 121222 tst vs 2021_07, and 120242's four epochs
# ----------------------------------------------------------------------------
say('\n' + '=' * 78)
say('7. STABILITY')
a, b = EP[('121222', '2021_07')], EP[('121222', '2021_07_tst')]
say('121222 2021_07 vs 2021_07_tst:')
say('   coef files identical: %s (L max |dc|=%.2e, E %.2e)' % (
    np.allclose(a['coefs']['L'], b['coefs']['L']) and np.allclose(a['coefs']['E'], b['coefs']['E']),
    np.max(np.abs(a['coefs']['L'] - b['coefs']['L'])), np.max(np.abs(a['coefs']['E'] - b['coefs']['E']))))
say('   lines: %d vs %d; lamp pairs: %d vs %d' % (len(a['deltawl']), len(b['deltawl']), len(a['lamp'][1]), len(b['lamp'][1])))
mm = a['fwhm'].merge(b['fwhm'], on='wl', suffixes=('', '_tst'))
say('   common lines %d: max |dFWHM_L|=%.4f, max |dFWHM_E|=%.4f' % (len(mm), (mm.L - mm.L_tst).abs().max(), (mm.E - mm.E_tst).abs().max()))
md = a['deltawl'].merge(b['deltawl'], on='wl', suffixes=('', '_tst'))
say('   max |d deltawl_L|=%.4f, |d deltawl_E|=%.4f' % ((md.L - md.L_tst).abs().max(), (md.E - md.E_tst).abs().max()))
say('   lines only in 2021_07: ' + ' '.join('%.1f' % w for w in set(a['fwhm'].wl) - set(b['fwhm'].wl)))
say('   poly tst - 2021_07 at the pixels of 400/550/700/850 nm: L ' + '/'.join(
    '%+.3f' % (poly(b['coefs']['L'], px_of(a['coefs']['L'], l)) - l) for l in (400, 550, 700, 850)) + '  E ' + '/'.join(
    '%+.3f' % (poly(b['coefs']['E'], px_of(a['coefs']['E'], l)) - l) for l in (400, 550, 700, 850)))
wa, pa = a['lamp']; wb, pb = b['lamp']
say('   lamp: wl columns identical %s; pair0 (Hg) identical %s; last pair (Xe) of 2021_07 == pair1 of tst: %s' % (
    np.allclose(wa, wb), np.allclose(pa[0][0], pb[0][0]) and np.allclose(pa[0][1], pb[0][1]),
    np.allclose(pa[-1][0], pb[-1][0]) and np.allclose(pa[-1][1], pb[-1][1])))

say('\n120242, four epochs: FWHM at common lines and the polynomial differences')
eps = [ep for (sn, ep) in EP if sn == '120242']
f0 = None
for ep in eps:
    f = EP[('120242', ep)]['fwhm']
    f0 = f if f0 is None else f0.merge(f, on='wl', suffixes=('', '_' + ep))
f0 = f0.rename(columns={'L': 'L_' + eps[0], 'E': 'E_' + eps[0]})
f0 = f0[(f0.wl > 380) & (f0.wl < 1100)]
say(f0.round(3).to_string(index=False))
f0.to_csv(os.path.join(OUT, 'stability_120242_fwhm.csv'), index=False)
say('   median over lines of FWHM(epoch) - FWHM(2020_09):')
for ep in eps[1:]:
    say('     %s: L %+.3f  E %+.3f' % (ep, (f0['L_' + ep] - f0['L_' + eps[0]]).median(), (f0['E_' + ep] - f0['E_' + eps[0]]).median()))
say('   polynomial differences lam_epoch(px) - lam_2020_09(px) at the pixels of 400/550/700/850 nm (nm):')
c0 = EP[('120242', eps[0])]['coefs']
srows = []
for ep in eps:
    c = EP[('120242', ep)]['coefs']
    d = dict(epoch=ep, date=c['date'].date())
    for ch in ('L', 'E'):
        for lam0 in (400, 550, 700, 850):
            pxx = px_of(c0[ch], lam0)
            d['%s_%d' % (ch, lam0)] = poly(c[ch], pxx) - poly(c0[ch], pxx)
    srows.append(d)
srows = pd.DataFrame(srows)
srows.to_csv(os.path.join(OUT, 'stability_120242_poly.csv'), index=False)
say(srows.round(3).to_string(index=False))
# same for 122302 and 122304 (two epochs each)
say('   two-epoch instruments, lam_2nd(px) - lam_1st(px) at the same pixels, and dFWHM median (2nd - 1st):')
for sn in ('122302', '122304'):
    e1, e2 = [EP[(s, ep)] for (s, ep) in EP if s == sn]
    f = e1['fwhm'].merge(e2['fwhm'], on='wl', suffixes=('_1', '_2'))
    f = f[(f.wl > 380) & (f.wl < 1100)]
    s = '     %s %s -> %s:' % (sn, [ep for (s_, ep) in EP if s_ == sn][0], [ep for (s_, ep) in EP if s_ == sn][1])
    for ch in ('L', 'E'):
        s += ' %s poly d(400/550/700/850)=' % ch + '/'.join(
            '%+.3f' % (poly(e2['coefs'][ch], px_of(e1['coefs'][ch], l)) - l) for l in (400, 550, 700, 850))
        s += ' dFWHM=%+.3f;' % (f[ch + '_2'] - f[ch + '_1']).median()
    say(s)

# ----------------------------------------------------------------------------
say('\n' + '=' * 78)
say('OUTPUT FILES in %s:' % OUT)
for f in sorted(os.listdir(OUT)):
    say('   ' + f)
with open(os.path.join(OUT, 'summary.txt'), 'w') as fh:
    fh.write('\n'.join(LOG) + '\n')

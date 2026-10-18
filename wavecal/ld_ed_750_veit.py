"""Ld(750)/Ed(750) sky index (Ruddick et al. 2006, eqs 23-24) for the
VEIT sample sequence, plus a variable inventory of the RELEASE_2 L2B and
PANTHYR L2A files, to see which products carry Ld and Ed.

Run from the repo root:  conda run -n ocean14 python wavecal/ld_ed_750_veit.py
"""
import glob
import os

import numpy as np
import xarray as xr

OS_COLOR = os.environ["OS_COLOR"]
WAV = os.path.join(OS_COLOR, "WATERHYPERNET", "Wavelengths")
REL2 = os.path.join(OS_COLOR, "WATERHYPERNET", "RELEASE_2")

LAM0 = 750.0


def band(ds, name, lam0=LAM0, half=5.0):
    """Mean of variable `name` over |wavelength - lam0| <= half nm."""
    wl = ds["wavelength"].values
    sel = np.abs(wl - lam0) <= half
    v = ds[name]
    wdim = "wavelength"
    return v.isel({wdim: np.where(sel)[0]}).mean(wdim)


def inventory(path, label):
    ds = xr.open_dataset(path)
    print(f"\n=== {label}: {os.path.basename(path)}")
    print("dims:", dict(ds.sizes))
    for k in ds.data_vars:
        print(f"  {k:40s} {ds[k].dims} {ds[k].attrs.get('units','')}")
    return ds


# ---- VEIT sample: L1A_IRR, L1A_RAD, L1C_ALL, L2A_REF ----------------------
f_irr = glob.glob(os.path.join(WAV, "*L1A_IRR*.nc"))[0]
f_rad = glob.glob(os.path.join(WAV, "*L1A_RAD*.nc"))[0]
f_l1c = glob.glob(os.path.join(WAV, "*L1C_ALL*.nc"))[0]
f_l2a = glob.glob(os.path.join(WAV, "*L2A_REF*.nc"))[0]

irr = inventory(f_irr, "L1A_IRR")
rad = inventory(f_rad, "L1A_RAD")
l1c = inventory(f_l1c, "L1C_ALL")
l2a = inventory(f_l2a, "L2A_REF")

# L1A: per-scan; Ld = sky view, vza >= 90 (vza from nadir)
vza = rad["viewing_zenith_angle"].values
sky = vza >= 90
print("\nL1A_RAD vza:", np.round(vza, 2), " sky scans:", sky.sum())
Ld_scans = band(rad, "radiance").values  # (scan,)
Ed_scans = band(irr, "irradiance").values
Ld_1a = Ld_scans[sky].mean()
Ed_1a = Ed_scans.mean()
print(f"L1A  Ld750 (sky mean) = {Ld_1a:.4f}  Ed750 mean = {Ed_1a:.3f}  "
      f"Ld/Ed = {Ld_1a/Ed_1a:.4f} sr^-1")
print("  per-scan Ld/Ed (sky):", np.round(Ld_scans[sky] / Ed_1a, 4))
print("  Ed scan CV:", f"{Ed_scans.std()/Ed_scans.mean():.4f}",
      " Ld(sky) scan CV:", f"{Ld_scans[sky].std()/Ld_scans[sky].mean():.4f}")

# L1C: already split
Ld_1c = band(l1c, "downwelling_radiance").values
Ed_1c = band(l1c, "irradiance").values
print(f"L1C  Ld750 = {np.atleast_1d(Ld_1c).mean():.4f}  Ed750 = "
      f"{np.atleast_1d(Ed_1c).mean():.3f}  Ld/Ed = "
      f"{np.atleast_1d(Ld_1c).mean()/np.atleast_1d(Ed_1c).mean():.4f}")
for k in ("solar_zenith_angle", "solar_azimuth_angle", "viewing_azimuth_angle",
          "pointing_azimuth_angle", "relative_azimuth_angle"):
    if k in l1c:
        print(f"  {k}: {np.unique(np.round(l1c[k].values, 1))}")

# L2A: does it carry Ld and Ed?
have = [k for k in ("downwelling_radiance", "irradiance", "upwelling_radiance",
                    "water_leaving_radiance", "rhof", "reflectance",
                    "reflectance_nosc") if k in l2a]
print("L2A carries:", have)
if "downwelling_radiance" in l2a and "irradiance" in l2a:
    Ld_2a = float(band(l2a, "downwelling_radiance").mean())
    Ed_2a = float(band(l2a, "irradiance").mean())
    print(f"L2A  Ld750 = {Ld_2a:.4f}  Ed750 = {Ed_2a:.3f}  Ld/Ed = {Ld_2a/Ed_2a:.4f}")
if "rhof" in l2a:
    print("  rhof:", float(l2a["rhof"].mean()))
print("  L2A attrs of interest:",
      {k: l2a.attrs[k] for k in l2a.attrs if "rho" in k.lower() or "sky" in k.lower()})

# ---- RELEASE_2: one HYPSTAR L2B and one PANTHYR L2A --------------------------
f_l2b = sorted(glob.glob(os.path.join(REL2, "VEIT_H", "2025", "11", "03", "*L2B*.nc")))[0]
l2b = inventory(f_l2b, "RELEASE_2 HYPSTAR L2B")
have = [k for k in ("downwelling_radiance", "irradiance", "upwelling_radiance",
                    "water_leaving_radiance", "rhof") if k in l2b]
print("L2B carries:", have)
if "downwelling_radiance" in l2b and "irradiance" in l2b:
    Ld = float(band(l2b, "downwelling_radiance").mean())
    Ed = float(band(l2b, "irradiance").mean())
    print(f"L2B  Ld750/Ed750 = {Ld/Ed:.4f}")
print("  L2B global attrs mentioning rho/sky/cloud:",
      {k: v for k, v in l2b.attrs.items()
       if any(s in k.lower() for s in ("rho", "sky", "cloud"))})

f_pan = sorted(glob.glob(os.path.join(REL2, "O1BE_P", "2026", "08", "03", "*L2A*.nc")))[0]
pan = inventory(f_pan, "RELEASE_2 PANTHYR L2A")

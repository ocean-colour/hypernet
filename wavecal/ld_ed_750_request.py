"""Ld(750)/Ed(750) for every sequence in docs/wiggles_data_request.csv, read
from the local RELEASE_2 L2B files, binned by SZA and relative azimuth.

Run from the repo root:  conda run -n ocean14 python wavecal/ld_ed_750_request.py
Writes $OS_COLOR/hypernet/wavecal/ld_ed_750_request.csv.
"""
import os

import numpy as np
import pandas as pd
import xarray as xr

from hypernet.whn_l1a import release2_path

LAM0, HALF = 750.0, 5.0

req = pd.read_csv("docs/wiggles_data_request.csv", dtype={"azimuth": str})
rows = []
for _, r in req.iterrows():
    path = release2_path(r.site, r.sequence_time, r.file)
    if not os.path.exists(path):
        rows.append(dict(ratio=np.nan))
        continue
    with xr.open_dataset(path) as ds:
        wl = ds["wavelength"].values
        sel = np.abs(wl - LAM0) <= HALF
        Ld = float(ds["downwelling_radiance"].values[sel].mean())
        Ed = float(ds["irradiance"].values[sel].mean())
        Ld_sza = float(ds["solar_zenith_angle"].values.mean())
        saa = float(ds["solar_azimuth_angle"].values.mean())
        vaa = float(ds["viewing_azimuth_angle"].values.mean())
        raa = abs(((vaa - saa) + 180) % 360 - 180)
        rows.append(dict(ratio=Ld / Ed, Ld750=Ld, Ed750=Ed, sza_file=Ld_sza,
                         raa=raa, n_valid=int(ds["n_valid_scans"].values.mean()),
                         rhof=float(ds["rhof"].values.mean())))
out = pd.concat([req, pd.DataFrame(rows)], axis=1)
odir = os.path.join(os.environ["OS_COLOR"], "hypernet", "wavecal")
os.makedirs(odir, exist_ok=True)
out.to_csv(os.path.join(odir, "ld_ed_750_request.csv"), index=False)

ok = out.dropna(subset=["ratio"])
print(f"{len(ok)}/{len(out)} sequences read")
print("relative azimuth values:", np.unique(np.round(ok.raa)))
print("\nratio quantiles:", np.round(ok.ratio.quantile([0, .1, .25, .5, .75, .9, 1]).values, 4))
print(f"fraction < 0.05: {(ok.ratio < 0.05).mean():.2f}   "
      f"< 0.03: {(ok.ratio < 0.03).mean():.2f}   > 0.10: {(ok.ratio > 0.10).mean():.2f}")

ok = ok.assign(sza_bin=pd.cut(ok.sza, [0, 35, 50, 65, 80]))
print("\nmedian / 10th-pct ratio by SZA bin (10th pct ~ clear-sky floor):")
print(ok.groupby("sza_bin", observed=True).ratio
        .agg(n="size", p10=lambda x: x.quantile(.1), median="median",
             frac_lt_005=lambda x: (x < 0.05).mean()).round(4))
print("\nby site:")
print(ok.groupby("site").ratio
        .agg(n="size", p10=lambda x: x.quantile(.1), median="median",
             frac_lt_005=lambda x: (x < 0.05).mean()).round(4))
print("\nby relative azimuth (rounded):")
print(ok.assign(raa_r=np.round(ok.raa / 5) * 5).groupby("raa_r").ratio
        .agg(n="size", p10=lambda x: x.quantile(.1), median="median").round(4))
# a rough clear-sky envelope: p10 vs 1/cos(sza)
clear = ok[ok.ratio < 0.05]
print("\nclear subset: corr(ratio, 1/cos sza) =",
      f"{np.corrcoef(clear.ratio, 1/np.cos(np.radians(clear.sza)))[0,1]:.2f}")

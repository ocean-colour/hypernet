"""Per-site instrument timeline from every Nth file in the whn_explore index."""
import os, sys, pandas as pd, netCDF4
STEP = int(sys.argv[1]) if len(sys.argv) > 1 else 20
root = os.path.join(os.environ['OS_COLOR'], 'hypernet', 'whn_explore')
idx = pd.read_parquet(os.path.join(root, 'index.parquet')).sort_values(['site', 'datetime'])
rows = []
for site, g in idx.groupby('site'):
    sub = g.iloc[::STEP]
    for _, r in sub.iterrows():
        try:
            ds = netCDF4.Dataset(r.path)
        except Exception as e:
            continue
        a = ds.__dict__
        if r.system == 'HYPSTAR':
            inst = a.get('system_id', ''); e_inst = inst
            cal = f"{a.get('instrument_calibration_date_rad','')}/{a.get('instrument_calibration_date_irr','')}"
            pv = a.get('processor_version', '')
        else:
            inst = a.get('l_sensor_sn', ''); e_inst = a.get('e_sensor_sn', '')
            cal = f"{a.get('l_sensor_cal','')}/{a.get('e_sensor_cal','')}"
            pv = a.get('processor_version', '')
        ds.close()
        rows.append(dict(site=site, datetime=r.datetime, L=inst, E=e_inst, cal=cal, proc=pv))
df = pd.DataFrame(rows)
out = os.path.join(os.environ['OS_COLOR'], 'hypernet', 'wavecal')
os.makedirs(out, exist_ok=True)
df.to_csv(os.path.join(out, 'instrument_timeline.csv'), index=False)
# contiguous runs
for site, g in df.groupby('site'):
    g = g.sort_values('datetime')
    key = g.L + '|' + g.E + '|' + g.cal
    run_id = (key != key.shift()).cumsum()
    print(f"\n== {site}  (sampled {len(g)} files, every {STEP}th)")
    for _, rr in g.groupby(run_id):
        print(f"   {rr.datetime.min():%Y-%m-%d} -> {rr.datetime.max():%Y-%m-%d}  L={rr.L.iloc[0]:16s} E={rr.E.iloc[0]:16s} cal={rr.cal.iloc[0]}  n={len(rr)}  proc={rr.proc.iloc[0]}")
# cross-site instrument reuse
print("\n== instruments seen at more than one site")
for inst, g in df.groupby('L'):
    sites = sorted(g.site.unique())
    if len(sites) > 1:
        for s in sites:
            gg = g[g.site == s]
            print(f"   {inst}: {s} {gg.datetime.min():%Y-%m-%d} -> {gg.datetime.max():%Y-%m-%d}")

#!/usr/bin/env python3
"""Does the wavefront actually move calibration skill?

Byte-identity is a mechanical property. The question that matters for a CALIBRATED model is
whether the parallel routing mode shifts sim-vs-obs skill (NSE / KGE / PBIAS) enough to change
what a modeller would conclude. This runs one calibrated model in three modes and scores each
against the same observed record.

  serial   OMP=1                              -> the reference the calibration was fitted to
  rtser    OMP=N, SWATPLUS_ROUTING_SERIAL=1   -> land-parallel, routing serial (reproducible)
  wave     OMP=N, full DAG wavefront          -> routing parallel (the fast, non-identical mode)

Scoring follows scripts/wb_inventory_explorer/score_wb_apply_experiment.py: monthly-mean
flo_out (m3/s) from channel_sd, observed cfs -> m3/s, months kept at >=90% daily coverage.

usage: calibration_impact.py <engine> <TxtInOut> <obs_csv> <gis_id> <workdir> [threads]
"""
import os, shutil, subprocess, sys
from pathlib import Path

import numpy as np, pandas as pd, xarray as xr

CFS = 0.0283168
NETCDF_LIB = "/data/SWATGenXApp/codes/lib/netcdf-ifx"
IFX_LIB = "/opt/intel/oneapi/compiler/2026.0/lib"
MODES = {"serial": (1, "1"), "rtser": (None, "1"), "wave": (None, "0")}


def run(engine, src, dst, threads, routing_serial):
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = f"{NETCDF_LIB}:{IFX_LIB}"
    env["OMP_NUM_THREADS"] = str(threads)
    env["SWATPLUS_ROUTING_SERIAL"] = routing_serial
    with open(dst / "engine.log", "w") as fh:
        subprocess.run([str(engine)], cwd=dst, env=env, stdout=fh,
                       stderr=subprocess.STDOUT, timeout=14400, check=False)


def sim_monthly(run_dir, gis_id):
    """Monthly-mean flo_out (m3/s) for one channel, from daily channel_sd NetCDF."""
    for cand in ("channel_sd_day.nc", "channel_sd_mon.nc"):
        p = run_dir / cand
        if p.exists():
            break
    else:
        return None
    with xr.open_dataset(p) as ds:
        gis = ds["gis_id"].values.astype(int)
        if gis_id not in gis:
            return None
        col = int(np.where(gis == gis_id)[0][0])
        flo = ds["flo_out"].values if "flo_out" in ds else ds["v41"].values
        flo = np.asarray(flo, dtype=np.float64)
        if flo.shape[0] != ds.sizes["time"]:
            flo = flo.T
        idx = pd.to_datetime(dict(year=ds["yrc"].values.astype(int),
                                  month=ds["mo"].values.astype(int),
                                  day=ds["day"].values.astype(int)
                                  if "day" in ds else 1))
        s = pd.Series(flo[:, col], index=idx)
    return s.resample("MS").mean()


def obs_monthly(path):
    df = pd.read_csv(path, parse_dates=["date"])
    df.loc[df["streamflow"] <= -90, "streamflow"] = np.nan
    q = df.set_index("date")["streamflow"] * CFS
    g = q.resample("MS")
    return g.mean()[(g.count() / g.size()) >= 0.90].dropna()


def skill(o, s):
    j = pd.concat([o, s], axis=1, join="inner").dropna()
    if len(j) < 12:
        return dict(n=len(j), nse=np.nan, kge=np.nan, pbias=np.nan)
    ov, sv = j.iloc[:, 0].values, j.iloc[:, 1].values
    den = ((ov - ov.mean()) ** 2).sum()
    nse = 1 - ((ov - sv) ** 2).sum() / den if den > 0 else np.nan
    pb = 100 * (ov - sv).sum() / ov.sum() if abs(ov.sum()) > 1e-12 else np.nan
    r = np.corrcoef(ov, sv)[0, 1]
    a, b = sv.std() / ov.std(), sv.mean() / ov.mean()
    kge = 1 - np.sqrt((r - 1) ** 2 + (a - 1) ** 2 + (b - 1) ** 2)
    return dict(n=len(j), nse=nse, kge=kge, pbias=pb)


def main():
    engine, src, obs_csv = Path(sys.argv[1]).resolve(), Path(sys.argv[2]), Path(sys.argv[3])
    gis_id, work = int(sys.argv[4]), Path(sys.argv[5])
    threads = int(sys.argv[6]) if len(sys.argv) > 6 else 8
    work.mkdir(parents=True, exist_ok=True)

    print(f"CALIBRATION IMPACT  model={src.parents[2].name}  gis_id={gis_id}  threads={threads}")
    sims = {}
    for mode, (t, rs) in MODES.items():
        n = t or threads
        print(f"  running {mode:6s} (OMP={n}, ROUTING_SERIAL={rs}) ...", flush=True)
        run(engine, src, work / mode, n, rs)
        sims[mode] = sim_monthly(work / mode, gis_id)

    o = obs_monthly(obs_csv)
    print(f"\n  observed months (>=90% coverage): {len(o)}\n")
    print(f"  {'mode':8s} {'n':>4s} {'NSE':>8s} {'KGE':>8s} {'PBIAS%':>8s}")
    base = None
    rows = {}
    for mode in MODES:
        if sims[mode] is None:
            print(f"  {mode:8s}   -- channel not found in output --")
            continue
        m = skill(o, sims[mode])
        rows[mode] = m
        print(f"  {mode:8s} {m['n']:4d} {m['nse']:8.4f} {m['kge']:8.4f} {m['pbias']:8.2f}")
        if base is None:
            base = m

    print(f"\n  DELTA vs serial (the reference the calibration was fitted to):")
    for mode in ("rtser", "wave"):
        if mode in rows and "serial" in rows:
            d = rows[mode]
            b = rows["serial"]
            print(f"  {mode:8s} dNSE={d['nse']-b['nse']:+.5f}  dKGE={d['kge']-b['kge']:+.5f}  "
                  f"dPBIAS={d['pbias']-b['pbias']:+.3f} pp")

    # pure sim-vs-sim: how far apart are the hydrographs themselves?
    print(f"\n  sim-vs-sim vs serial (hydrograph displacement, not skill):")
    for mode in ("rtser", "wave"):
        if sims.get(mode) is None or sims.get("serial") is None:
            continue
        j = pd.concat([sims["serial"], sims[mode]], axis=1, join="inner").dropna()
        b, t = j.iloc[:, 0].values, j.iloc[:, 1].values
        dmean = 100 * (t.mean() - b.mean()) / b.mean()
        dmax = 100 * np.max(np.abs(t - b) / np.maximum(np.abs(b), 1e-12))
        print(f"  {mode:8s} mean flow {dmean:+.4f}%   worst month {dmax:.4f}%   "
              f"corr={np.corrcoef(b, t)[0,1]:.8f}")


if __name__ == "__main__":
    main()

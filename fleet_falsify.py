#!/usr/bin/env python3
"""Try to BREAK two claims across the calibrated fleet, on adversarially chosen basins.

  CLAIM A  HRU-parallel + serial routing is identical to serial   <- the shipped mode's guarantee
  CLAIM B  the full wavefront materially perturbs routing/nutrients

Basins are ordered by reservoir+wetland density, because the race that was fixed lived in
reservoir/wetland release logic reached from the parallel land phase — so that is where
CLAIM A is most likely to fail. Every basin is run under stress: in-stream water quality
forced ON, channel routing + sediment ON, full daily NetCDF output for every object.

Results are appended to a JSON-lines file after each basin so partial runs are usable.

usage: fleet_falsify.py <engine> <out.jsonl> [n_basins] [years] [threads]
"""
import json, os, shutil, subprocess, sys, time
from pathlib import Path

import numpy as np, pandas as pd, xarray as xr

HERE = Path(__file__).resolve().parent
SCRATCH = HERE / "scratch" / "fleet"
NETCDF_LIB = "/data/SWATGenXApp/codes/lib/netcdf-ifx"
IFX_LIB = "/opt/intel/oneapi/compiler/2026.0/lib"
INV = HERE / "scratch" / "fleet_inventory.json"

_BASE = ["flo", "sed", "orgn", "sedp", "no3", "solp", "chla", "nh3", "no2",
         "cbod", "dox", "san", "sil", "cla", "sag", "lag", "grv", "temp"]
SD = (["area", "precip", "evap", "seep"] + [f"{c}_stor" for c in _BASE]
      + [f"{c}_in" for c in _BASE] + [f"{c}_out" for c in _BASE])


def run(engine, src, dst, threads, routing_serial):
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = f"{NETCDF_LIB}:{IFX_LIB}"
    env["OMP_NUM_THREADS"] = str(threads)
    env["SWATPLUS_ROUTING_SERIAL"] = routing_serial
    t0 = time.time()
    with open(dst / "engine.log", "w") as fh:
        rc = subprocess.run([str(engine)], cwd=dst, env=env, stdout=fh,
                            stderr=subprocess.STDOUT, timeout=36000).returncode
    return rc, time.time() - t0


def compare_streaming(a_dir, b_dir):
    """Max abs difference + differing-cell count across every output, ONE FILE AT A TIME.

    The first version of this loaded every variable of every NetCDF for BOTH runs into
    memory at once and was OOM-killed on the first basin (1.2 GB of output per run, x2,
    widened to float64). Streaming file-by-file keeps the footprint to a single variable
    pair and makes big-basin runs survivable.
    """
    worst, ncells, nvars, total = 0.0, 0, 0, 0
    worst_key = ""
    for pa in sorted(Path(a_dir).glob("*.nc")):
        pb = Path(b_dir) / pa.name
        if not pb.exists():
            continue
        try:
            da, db = xr.open_dataset(pa), xr.open_dataset(pb)
        except Exception:
            continue
        try:
            idn = next((n for n in ("obj_id", "id", "gis_id", "unit") if n in da.variables), None)
            order = None
            if idn is not None:
                ids = np.asarray(da.variables[idn][:]).ravel()
                if np.isfinite(ids).all():
                    order = np.argsort(ids)
            for v in da.variables:
                if da[v].dtype.kind != "f" or v == idn or v not in db.variables:
                    continue
                x = np.asarray(da[v].values, dtype=np.float64)
                y = np.asarray(db[v].values, dtype=np.float64)
                if x.shape != y.shape:
                    continue
                total += x.size
                d = np.abs(x - y)
                nz = int((d > 0).sum())
                if nz:
                    nvars += 1
                    ncells += nz
                    m = float(d.max())
                    if m > worst:
                        worst, worst_key = m, f"{pa.name}:{v}"
                del x, y, d
        finally:
            da.close(); db.close()
    return dict(max_abs=worst, worst_var=worst_key, n_cells=ncells,
                n_vars=nvars, n_total=total)


def chan(d):
    with xr.open_dataset(Path(d) / "channel_sd_day.nc") as ds:
        out = {}
        for v in ds.variables:
            if not v.startswith("v"):
                continue
            i = int(v[1:]) - 1
            if i < len(SD):
                a = np.asarray(ds[v].values, dtype=np.float64)
                out[SD[i]] = a[:, 0] if (a.ndim == 2 and a.shape[1] == 1) else a
    return out


def dev(S, W, key):
    if key not in S or key not in W:
        return None
    s, w = np.ravel(S[key]), np.ravel(W[key])
    m = s > (np.nanmax(s) * 1e-4) if np.nanmax(s) > 0 else np.zeros_like(s, bool)
    if not m.any():
        return None
    e = 100 * np.abs(w[m] - s[m]) / s[m]
    net = 100 * (np.nansum(w) - np.nansum(s)) / np.nansum(s) if np.nansum(s) else np.nan
    return dict(median=float(np.median(e)), p90=float(np.percentile(e, 90)),
                worst=float(e.max()), net=float(net))


def main():
    engine = Path(sys.argv[1]).resolve()
    outp = Path(sys.argv[2])
    nbas = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    years = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    threads = int(sys.argv[5]) if len(sys.argv) > 5 else 8

    inv = json.load(open(INV))
    inv.sort(key=lambda r: -(r["res"] + r["wet"]))
    seen, picks = set(), []
    for r in inv:
        if r["s"] in seen:
            continue
        seen.add(r["s"])
        picks.append(r)
        if len(picks) >= nbas:
            break

    SCRATCH.mkdir(parents=True, exist_ok=True)
    for i, r in enumerate(picks, 1):
        tag = f"{r['s']}_{r['mn']}"
        work = SCRATCH / tag
        print(f"\n[{i}/{len(picks)}] {r['s']}  HRU {r['hru']} cha {r['cha']} "
              f"res {r['res']} wet {r['wet']}", flush=True)
        model = work / "model"
        subprocess.run([sys.executable, str(HERE / "prep_test_model.py"),
                        r["t"], str(model), str(years)], check=True,
                       stdout=subprocess.DEVNULL)
        rec = dict(site=r["s"], model=r["mn"], hru=r["hru"], cha=r["cha"],
                   res=r["res"], wet=r["wet"], aqu=r["aqu"], years=years,
                   threads=threads, engine=engine.name)
        try:
            for lab, n, rs in (("serial", 1, "1"), ("rtser", threads, "1"),
                               ("wave", threads, "0")):
                rc, el = run(engine, model, work / lab, n, rs)
                rec[f"{lab}_rc"], rec[f"{lab}_s"] = rc, round(el, 1)
                print(f"    {lab:7s} rc={rc} {el:7.1f}s", flush=True)
            rec["A_hru_parallel_vs_serial"] = compare_streaming(work / "serial", work / "rtser")
            rec["B_wave_vs_serial"] = compare_streaming(work / "serial", work / "wave")
            S, W = chan(work / "serial"), chan(work / "wave")
            R = chan(work / "rtser")
            rec["wave_dev"] = {k: dev(S, W, k) for k in
                               ("flo_out", "orgn_out", "sedp_out", "nh3_out", "no3_out")}
            rec["hru_dev"] = {k: dev(S, R, k) for k in ("flo_out", "orgn_out", "sedp_out")}
            a = rec["A_hru_parallel_vs_serial"]
            print(f"    CLAIM A  max|Δ|={a['max_abs']:.3e}  cells={a['n_cells']}/"
                  f"{a['n_total']}  -> {'HOLDS' if a['n_cells']==0 else 'BROKEN'}", flush=True)
            wd = rec["wave_dev"].get("orgn_out")
            if wd:
                print(f"    CLAIM B  orgn median {wd['median']:.1f}% net {wd['net']:+.1f}%",
                      flush=True)
        except Exception as e:
            rec["error"] = repr(e)
            print(f"    ERROR {e!r}", flush=True)
        with open(outp, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        # keep disk in check: drop the raw runs, keep the record
        for lab in ("serial", "rtser", "wave"):
            shutil.rmtree(work / lab, ignore_errors=True)
        shutil.rmtree(model, ignore_errors=True)


if __name__ == "__main__":
    main()

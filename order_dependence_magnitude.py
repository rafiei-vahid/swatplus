#!/usr/bin/env python3
"""How much does the SEQUENTIAL answer change if the order-dependence defects stay in?

The engine paper reports that two of the mechanisms it removed are defects in the RELEASED
SEQUENTIAL engine, not artifacts of parallelism. That claim is about existence. This measures
magnitude: run ONE model, single-threaded, under two binaries built from the SAME tree with the
SAME flags, differing ONLY in whether the two defects are present, and report how far apart the
answers are.

  arm A (fixed)  -- enratio imported from hru_module; gra zeroed before the algcon guard
  arm B (defect) -- enratio shadowed by a saved local; gra left to carry over

Both run at OMP_NUM_THREADS=1, so no parallelism is involved anywhere in the comparison. Any
difference is the defect changing the sequential answer.

usage: order_dependence_magnitude.py <engineA> <engineB> <TxtInOut> <workdir> [out.json]
"""
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import xarray as xr

NETCDF_LIB = "/data/SWATGenXApp/codes/lib/netcdf-ifx"
IFX_LIB = "/opt/intel/oneapi/compiler/2026.0/lib"

# The two defects have named consumers in the upstream report; flow is the control that should
# not move, because neither mechanism reaches the water balance.
EXPECTED = {
    "flo_out": "control - neither defect reaches the water balance",
    "chla_out": "gra (ch_watqual4) -> factk -> alg_m -> algcon_out",
    "solp_out": "gra (ch_watqual4) -> zz",
    "cbod_out": "enratio (varinit) -> swr_subwq org_c -> cbodu",
    "dox_out": "enratio (varinit) -> cbodu -> doxq",
}


def run(engine, src, dst):
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = f"{NETCDF_LIB}:{IFX_LIB}"
    env["OMP_NUM_THREADS"] = "1"          # single-threaded on BOTH arms
    env["SWATPLUS_ROUTING_SERIAL"] = "1"  # and no routing wavefront either
    with open(dst / "engine.log", "w") as fh:
        p = subprocess.run([str(engine)], cwd=dst, env=env, stdout=fh,
                           stderr=subprocess.STDOUT, timeout=14400, check=False)
    return p.returncode


def load(run_dir):
    p = run_dir / "channel_sd_day.nc"
    if not p.exists():
        raise SystemExit(f"no channel_sd_day.nc in {run_dir}")
    ds = xr.open_dataset(p)
    return ds


def relstats(a, b, floor_frac=1e-4):
    """Relative difference on cells that carry signal.

    Cells at or near zero make percentages meaningless (a 1e-12 vs 2e-12 is 100%), so cells are
    kept only where the reference exceeds floor_frac of that variable's peak. The retained count
    IS the denominator and is reported with every statistic.
    """
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    finite = np.isfinite(a) & np.isfinite(b)
    a, b = a[finite], b[finite]
    peak = np.nanmax(np.abs(a)) if a.size else 0.0
    if peak <= 0:
        return None
    keep = np.abs(a) > floor_frac * peak
    n_kept = int(keep.sum())
    if n_kept == 0:
        return None
    rel = np.abs(b[keep] - a[keep]) / np.abs(a[keep]) * 100.0
    ndiff = int((b != a).sum())
    return {
        "n_total": int(a.size),
        "n_differing_bitwise": ndiff,
        "pct_cells_differing": round(100.0 * ndiff / a.size, 3),
        "n_scored": n_kept,
        "median_abs_pct": float(np.median(rel)),
        "p90_abs_pct": float(np.percentile(rel, 90)),
        "max_abs_pct": float(np.max(rel)),
        "sum_A": float(np.nansum(a)),
        "sum_B": float(np.nansum(b)),
        "net_pct": float((np.nansum(b) - np.nansum(a)) / np.nansum(a) * 100.0)
        if np.nansum(a) else None,
    }


def main():
    engA, engB, txt, work = (Path(x) for x in sys.argv[1:5])
    out = Path(sys.argv[5]) if len(sys.argv) > 5 else work / "order_dependence_magnitude.json"
    work.mkdir(parents=True, exist_ok=True)

    result = {
        "fixture": str(txt),
        "engine_fixed": str(engA),
        "engine_defect": str(engB),
        "threads": 1,
        "routing_serial": True,
    }

    for name, eng in (("A_fixed", engA), ("B_defect", engB)):
        rc = run(eng, txt, work / name)
        result[f"rc_{name}"] = rc
        print(f"{name}: rc={rc}", flush=True)

    dsA, dsB = load(work / "A_fixed"), load(work / "B_defect")
    shared = [v for v in dsA.data_vars if v in dsB.data_vars
              and np.issubdtype(dsA[v].dtype, np.floating)]
    result["n_variables_compared"] = len(shared)

    per_var, n_var_differ = {}, 0
    for v in sorted(shared):
        st = relstats(dsA[v].values, dsB[v].values)
        if st is None:
            continue
        if st["n_differing_bitwise"]:
            n_var_differ += 1
        per_var[v] = st
    result["n_variables_differing"] = n_var_differ
    result["variables"] = per_var
    result["named_consumers"] = {k: {"mechanism": m, "stats": per_var.get(k)}
                                 for k, m in EXPECTED.items()}

    out.write_text(json.dumps(result, indent=2))
    print(f"\nwrote {out}")

    print(f"\n{len(shared)} float variables compared; {n_var_differ} differ bitwise\n")
    hdr = f"{'variable':16} {'cells differ %':>14} {'median %':>10} {'p90 %':>9} {'worst %':>10} {'net %':>9}"
    print(hdr)
    print("-" * len(hdr))
    for v, st in sorted(per_var.items(), key=lambda kv: -kv[1]["median_abs_pct"])[:20]:
        net = f"{st['net_pct']:9.3f}" if st["net_pct"] is not None else "        -"
        print(f"{v:16} {st['pct_cells_differing']:14.3f} {st['median_abs_pct']:10.3f} "
              f"{st['p90_abs_pct']:9.2f} {st['max_abs_pct']:10.2f} {net}")


if __name__ == "__main__":
    main()

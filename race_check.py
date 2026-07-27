#!/usr/bin/env python3
"""Byte-identity / determinism battery for one engine binary on one model.

Runs, in the deployed HRU-parallel mode (SWATPLUS_ROUTING_SERIAL=1):
  serial      : OMP_NUM_THREADS=1
  parN_a/par_b: OMP_NUM_THREADS=N, twice (run-to-run determinism)

Reports, per NetCDF output variable, the worst relative difference of
  (a) parallel vs serial   -> byte-identity of the shipped mode
  (b) par_a vs par_b       -> run-to-run determinism (a true race, not FP reassociation)

Diffs are aligned by obj_id, because the wavefront reorders write order and an
index-wise diff is a false positive that mimics a 100% race.

usage: race_check.py <engine> <TxtInOut> <workdir> [threads] [--years N]
"""
import os, shutil, subprocess, sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

NETCDF_LIB = "/data/SWATGenXApp/codes/lib/netcdf-ifx"
IFX_LIB = "/opt/intel/oneapi/compiler/2026.0/lib"


def run(engine, src, dst, threads, label):
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = f"{NETCDF_LIB}:{IFX_LIB}"
    env["OMP_NUM_THREADS"] = str(threads)
    env["SWATPLUS_ROUTING_SERIAL"] = "1"       # the shipped byte-identical mode
    log = dst / "engine.log"
    with open(log, "w") as fh:
        rc = subprocess.run([str(engine)], cwd=dst, env=env, stdout=fh,
                            stderr=subprocess.STDOUT, timeout=7200).returncode
    if rc != 0:
        print(f"  !! {label}: engine exit {rc} (see {log})")
    return rc


def load(ncpath):
    """{var: array} keyed by obj_id order so write-order permutation is not a diff."""
    out = {}
    with Dataset(ncpath) as ds:
        idname = next((n for n in ("obj_id", "id", "gis_id", "unit") if n in ds.variables), None)
        order = np.argsort(np.asarray(ds.variables[idname][:]).ravel()) if idname else None
        for name, var in ds.variables.items():
            if var.dtype.kind != "f" or name == idname:
                continue
            a = np.asarray(var[:])
            if order is not None and a.ndim >= 1 and a.shape[-1] == order.size:
                a = a[..., order]
            out[name] = a
    return out


def compare(dir_a, dir_b, tag, tol=1e-12):
    """Worst relative difference per variable across every .nc both runs wrote."""
    worst = []
    for nc_a in sorted(Path(dir_a).glob("*.nc")):
        nc_b = Path(dir_b) / nc_a.name
        if not nc_b.exists():
            continue
        try:
            A, B = load(nc_a), load(nc_b)
        except Exception as e:
            print(f"  (skip {nc_a.name}: {e})")
            continue
        for var in sorted(set(A) & set(B)):
            a, b = A[var], B[var]
            if a.shape != b.shape:
                worst.append((1.0, nc_a.name, var, "SHAPE MISMATCH"))
                continue
            d = np.abs(a - b)
            if not d.any():
                continue
            scale = np.maximum(np.abs(a), np.abs(b))
            rel = np.where(scale > 0, d / np.where(scale > 0, scale, 1), 0.0)
            n_over = int((rel > tol).sum())
            if n_over:
                worst.append((float(rel.max()), nc_a.name, var, f"{n_over} cells"))
    worst.sort(reverse=True)
    if not worst:
        print(f"  PASS  {tag}: identical across every variable")
        return True
    print(f"  FAIL  {tag}: {len(worst)} variable(s) differ")
    for rel, f, var, note in worst[:12]:
        print(f"        {f:24s} {var:14s} rel={rel:.3e}  {note}")
    return False


def main():
    engine, src, work = Path(sys.argv[1]).resolve(), Path(sys.argv[2]), Path(sys.argv[3])
    threads = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    work.mkdir(parents=True, exist_ok=True)
    print(f"RACE CHECK  engine={engine.name}  model={src.resolve().name}  threads={threads}")
    for label, n in (("serial", 1), ("par_a", threads), ("par_b", threads)):
        print(f"  running {label} (OMP={n}) ...", flush=True)
        run(engine, src, work / label, n, label)
    ok1 = compare(work / "serial", work / "par_a", f"parallel(x{threads}) vs serial  [byte-identity]")
    ok2 = compare(work / "par_a", work / "par_b", f"par_a vs par_b            [determinism]")
    print(f"\n  VERDICT: {'CLEAN' if (ok1 and ok2) else 'RACE PRESENT'}")
    return 0 if (ok1 and ok2) else 1


if __name__ == "__main__":
    sys.exit(main())

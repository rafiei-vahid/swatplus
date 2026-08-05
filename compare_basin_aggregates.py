#!/usr/bin/env python3
"""Serial vs parallel comparison of basin-level NetCDF aggregates.

Basin aggregate outputs have never been compared by any gate, because every fixture
examined ships with `basin_*` disabled in print.prt. This compares them for the first
time, at mon/yr/aa (the daily interval is a separate open question).

Reports in the form the paper demands: files, variables compared, variables carrying
NON-ZERO data, cells, differing. A comparison over all-zero variables proves nothing
regardless of the differing count, so the non-zero column is not optional.

PRECONDITION (s-302d7da0's free control): channel_sd_day.nc is what the whole
certification rests on, so if it has zero records the run printed nothing and no other
file in the directory carries information. Abort rather than compare.

usage: compare_basin_aggregates.py <ref_dir> <cmp_dir>
"""
import sys
from pathlib import Path

import numpy as np
import xarray as xr


def records(path):
    if not path.exists():
        return None
    with xr.open_dataset(path) as ds:
        return int(ds.sizes.get("time", 0))


def main():
    ref, cmp_ = Path(sys.argv[1]), Path(sys.argv[2])

    # --- precondition: did the run print anything at all? ---
    print("PRECONDITION  channel_sd_day.nc records")
    for lbl, d in (("reference", ref), ("compare  ", cmp_)):
        n = records(d / "channel_sd_day.nc")
        print(f"  {lbl}: {n if n is not None else 'FILE MISSING'}")
        if not n:
            print("  ABORT: the run printed no records (check nyskip vs simulated years).")
            print("  No basin result from this directory carries information.")
            return 1
    print("  precondition PASSED\n")

    hdr = f"{'file':22} {'vars':>5} {'non-zero':>9} {'cells':>10} {'differing':>10}"
    print(hdr); print("-" * len(hdr))
    tot_v = tot_nz = tot_c = tot_d = 0
    for name in ("basin_wb_mon.nc", "basin_wb_yr.nc", "basin_wb_aa.nc"):
        a_p, b_p = ref / name, cmp_ / name
        if not (a_p.exists() and b_p.exists()):
            print(f"{name:22} {'MISSING':>5}")
            continue
        A, B = xr.open_dataset(a_p), xr.open_dataset(b_p)
        shared = [v for v in A.data_vars if v in B.data_vars
                  and np.issubdtype(A[v].dtype, np.floating)]
        nz = cells = diff = 0
        for v in shared:
            a = np.asarray(A[v].values, dtype=np.float64)
            b = np.asarray(B[v].values, dtype=np.float64)
            if a.shape != b.shape:
                print(f"  !! {name}:{v} shape mismatch {a.shape} vs {b.shape}")
                continue
            if a.size and np.nanmax(np.abs(a)) > 0:
                nz += 1
            cells += a.size
            diff += int((a != b).sum())
        print(f"{name:22} {len(shared):5} {nz:9} {cells:10,} {diff:10,}")
        tot_v += len(shared); tot_nz += nz; tot_c += cells; tot_d += diff
    print("-" * len(hdr))
    print(f"{'TOTAL':22} {tot_v:5} {tot_nz:9} {tot_c:10,} {tot_d:10,}")

    print()
    if tot_nz == 0:
        print("!! VERDICT VOID: every compared variable is all-zero. This proves nothing,")
        print("   whatever the differing count says. Diagnose the fixture, do not explain the result.")
        return 1
    if tot_d == 0:
        print(f"VERDICT: basin aggregates are BIT-IDENTICAL serial vs parallel")
        print(f"         over {tot_c:,} cells, {tot_nz} of {tot_v} variables carrying non-zero data.")
    else:
        print(f"VERDICT: {tot_d:,} of {tot_c:,} cells DIFFER ({100*tot_d/tot_c:.4f}%),")
        print(f"         {tot_nz} of {tot_v} variables carrying non-zero data.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Make a small, fast, race-exposing test model from a real TxtInOut.

The race only surfaces with in-stream water quality active over a multi-year run
(publication/engine-acceleration/REVISION_NOTES_parallel_nutrient.md), and every
process must be observable, so this turns on wq_cha + full daily NetCDF print and
trims the run to a few years.

usage: prep_test_model.py <src TxtInOut> <dst> [nyears]
"""
import shutil, sys
from pathlib import Path

SRC, DST = Path(sys.argv[1]), Path(sys.argv[2])
NYEARS = int(sys.argv[3]) if len(sys.argv) > 3 else 4

if DST.exists():
    shutil.rmtree(DST)
shutil.copytree(SRC, DST)

# --- time.sim: last NYEARS of the available weather record -------------------
t = (DST / "time.sim").read_text().splitlines()
f = t[2].split()
yr_end = int(f[3])
f[0], f[1], f[2], f[3] = "1", str(yr_end - NYEARS + 1), "365", str(yr_end)
t[2] = "  " + "  ".join(f)
(DST / "time.sim").write_text("\n".join(t) + "\n")

# --- codes.bsn: turn ON in-stream WQ + channel routing/sediment --------------
c = (DST / "codes.bsn").read_text().splitlines()
hdr, val = c[1].split(), c[2].split()
for key, want in (("wq_cha", "1"), ("rte_cha", "1"), ("sed_cha", "1")):
    if key in hdr:
        val[hdr.index(key)] = want
c[2] = " " + " ".join(val)
(DST / "codes.bsn").write_text("\n".join(c) + "\n")

# --- print.prt: nyskip=0 and DAILY on for every object ----------------------
p = (DST / "print.prt").read_text().splitlines()
n = p[2].split()
n[0] = "0"
p[2] = "  " + "  ".join(n)
out, seen_objects = [], False
for line in p:
    if line.strip().startswith("objects"):
        seen_objects = True
        out.append(line)
        continue
    if seen_objects and line.strip():
        parts = line.split()
        if len(parts) == 5 and parts[1] in ("y", "n"):
            parts[1] = "y"                     # daily on, every object
            out.append(f"{parts[0]:<12}{parts[1]:<5}{parts[2]:<5}{parts[3]:<5}{parts[4]}")
            continue
    out.append(line)
(DST / "print.prt").write_text("\n".join(out) + "\n")

print(f"prepared {DST}")
print(f"  years   : {f[1]}-{f[3]}")
print(f"  wq_cha  : on")
print(f"  print   : daily, all objects, cdfout")

#!/usr/bin/env python3
"""Serial vs parallel comparison for the management-race experiment.

Evaluates two gates on one basin:
  G2  arm C (threadprivate(mgt) removed) must FAIL  -- the fixture has power
  G3  arm A (ship engine) must PASS                 -- the shipped engine is correct

mgt_out.txt is compared SORTED. Under threads the same operations are emitted in a
different sequence because HRUs complete out of order; that is a write-order artifact,
established on two independent models, and comparing it unsorted reports a race that
does not exist. Every other output is compared byte for byte.

usage: compare_mgt_experiment.py <run_dir_serial> <run_dir_parallel> <label>
"""
import sys
from pathlib import Path

SORTED_COMPARE = {"mgt_out.txt"}


def compare(a: Path, b: Path):
    rows, n_diff_files, n_diff_lines = [], 0, 0
    names = sorted({p.name for p in a.glob("*.txt")} & {p.name for p in b.glob("*.txt")})
    for name in names:
        la = (a / name).read_text(errors="replace").split("\n")
        lb = (b / name).read_text(errors="replace").split("\n")
        note = ""
        if name in SORTED_COMPARE:
            la_c, lb_c = sorted(la), sorted(lb)
            raw_same = la == lb
            note = "sorted" if raw_same else "sorted (raw order differs)"
        else:
            la_c, lb_c = la, lb
        if len(la_c) != len(lb_c):
            d = abs(len(la_c) - len(lb_c))
        else:
            d = sum(1 for x, y in zip(la_c, lb_c) if x != y)
        if d:
            n_diff_files += 1
            n_diff_lines += d
        rows.append((name, len(la_c), d, note))
    return rows, n_diff_files, n_diff_lines


def main():
    a, b, label = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
    rows, nf, nl = compare(a, b)
    print(f"=== {label}: {a.name} vs {b.name} ===")
    print(f"{'file':26} {'lines':>10} {'differing':>10}  note")
    for name, n, d, note in rows:
        flag = "  <<<" if d else ""
        print(f"{name:26} {n:10,} {d:10,}  {note}{flag}")
    print(f"{'':26} {'':>10} {'':>10}")
    print(f"files compared {len(rows)} | files differing {nf} | lines differing {nl:,}")
    return nf


if __name__ == "__main__":
    sys.exit(0 if main() == 0 else 1)

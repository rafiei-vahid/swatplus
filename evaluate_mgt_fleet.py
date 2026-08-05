#!/usr/bin/env python3
"""Fleet gate evaluation for the management-race experiment.

Consumes the `out/` directory fetched from the AWS run (5 basins x 2 arms x 2 thread counts)
and reports, per basin:

  G2  arm C (threadprivate(mgt) REMOVED) must DIFFER serial-vs-parallel.
      This is the POSITIVE CONTROL. A fixture that passes G2 has demonstrated it can express
      the race. Without it, arm A passing means nothing -- an empty schedule also "passes".
  G3  arm A (the shipped engine) must be IDENTICAL serial-vs-parallel.

The gates are asymmetric on purpose. G3 is the claim; G2 is what makes G3 admissible.
A basin where G2 fails is reported as INADMISSIBLE, not as a pass, because its arm-A result
carries no information either way.

mgt_out.txt is compared SORTED: under threads the same operations are emitted in a different
sequence because HRUs complete out of order. That is a write-order artifact established on two
independent models; comparing it unsorted reports a race that does not exist. Every other output
is compared byte for byte. The RAW order difference is still reported, so the artifact stays
visible rather than being quietly normalised away.

usage: evaluate_mgt_fleet.py <out_dir>
"""
import sys
from collections import OrderedDict
from pathlib import Path

SORTED_COMPARE = {"mgt_out.txt"}
EXPECT_RUNS = 20


def compare(a: Path, b: Path):
    """Return (n_files, n_diff_files, n_diff_lines, raw_order_only)."""
    names = sorted({p.name for p in a.glob("*.txt")} & {p.name for p in b.glob("*.txt")})
    n_diff_files = n_diff_lines = 0
    raw_order_only = []
    for name in names:
        la = (a / name).read_text(errors="replace").split("\n")
        lb = (b / name).read_text(errors="replace").split("\n")
        if name in SORTED_COMPARE:
            if la != lb:
                raw_order_only.append(name)
            la, lb = sorted(la), sorted(lb)
        if len(la) != len(lb):
            d = abs(len(la) - len(lb))
        else:
            d = sum(1 for x, y in zip(la, lb) if x != y)
        if d:
            n_diff_files += 1
            n_diff_lines += d
    return len(names), n_diff_files, n_diff_lines, raw_order_only


def rc_of(d: Path) -> str:
    f = d / "RC"
    return f.read_text().strip() if f.exists() else "MISSING"


def main() -> int:
    out = Path(sys.argv[1])
    dirs = sorted(p for p in out.iterdir() if p.is_dir())

    # Denominator first. A gate table over an unknown number of runs is not a result.
    bad_rc = [(p.name, rc_of(p)) for p in dirs if rc_of(p) != "0"]
    print(f"runs present: {len(dirs)}/{EXPECT_RUNS} | nonzero-or-missing rc: {len(bad_rc)}")
    for n, rc in bad_rc:
        print(f"   !! {n} rc={rc}")
    if len(dirs) != EXPECT_RUNS or bad_rc:
        print("\nINADMISSIBLE: the run set is incomplete or contains failures.")
        print("Fix the run before reading any gate below -- a missing run is not a passing run.")
        if not dirs:
            return 1

    basins = list(OrderedDict.fromkeys(p.name.rsplit("_", 2)[0] for p in dirs))
    print(f"\n{'basin':12} {'arm':4} {'files':>6} {'difffiles':>10} {'difflines':>11}  gate")
    print("-" * 68)

    verdicts = {}
    for basin in basins:
        for arm, gate, want_diff in (("C", "G2", True), ("A", "G3", False)):
            s, p = out / f"{basin}_{arm}_t1", out / f"{basin}_{arm}_t8"
            if not (s.is_dir() and p.is_dir()):
                print(f"{basin:12} {arm:4} {'--':>6} {'--':>10} {'--':>11}  {gate} NO DATA")
                verdicts[(basin, gate)] = "NO DATA"
                continue
            nf, ndf, ndl, raw = compare(s, p)
            differs = ndf > 0
            ok = differs if want_diff else not differs
            verdicts[(basin, gate)] = "PASS" if ok else "FAIL"
            note = f"{gate} {'PASS' if ok else 'FAIL'}"
            if raw:
                note += f"  [{','.join(raw)} raw order differs, compared sorted]"
            print(f"{basin:12} {arm:4} {nf:6} {ndf:10} {ndl:11,}  {note}")

    print("\nsummary")
    for basin in basins:
        g2, g3 = verdicts.get((basin, "G2")), verdicts.get((basin, "G3"))
        # NO DATA is not FAIL. Collapsing them is the same defect these gates exist to catch:
        # a basin whose arm-A runs never landed was being announced as "the shipped engine
        # DIFFERS under threads", which is a fabricated finding, not a cautious one.
        if g2 == "NO DATA" or g3 == "NO DATA":
            missing = [g for g, v in (("G2", g2), ("G3", g3)) if v == "NO DATA"]
            print(f"  {basin:12} NOT MEASURED — {'/'.join(missing)} run(s) absent; "
                  f"no verdict either way")
        elif g2 != "PASS":
            print(f"  {basin:12} INADMISSIBLE — positive control failed (G2={g2}); "
                  f"arm A result carries no information")
        elif g3 == "PASS":
            print(f"  {basin:12} shipped engine reproduces the serial result, "
                  f"on a fixture proven able to break it")
        else:
            print(f"  {basin:12} *** G3 FAIL — the shipped engine DIFFERS under threads ***")

    measured = [b for b in basins
                if "NO DATA" not in (verdicts.get((b, "G2")), verdicts.get((b, "G3")))]
    admissible = [b for b in measured if verdicts.get((b, "G2")) == "PASS"]
    passed = [b for b in admissible if verdicts.get((b, "G3")) == "PASS"]
    print(f"\n{len(passed)}/{len(admissible)} admissible basins pass G3 | "
          f"{len(basins) - len(measured)} not measured, "
          f"{len(measured) - len(admissible)} measured but inadmissible "
          f"(of {len(basins)} basins, {EXPECT_RUNS} expected runs)")
    # A clean exit requires the FULL fleet measured and admissible. Returning 0 because the one
    # basin that happened to land looked fine is how a partial run gets read as a result.
    complete = len(measured) == len(basins) == len(admissible) and len(dirs) == EXPECT_RUNS
    return 0 if complete and len(passed) == len(admissible) else 1


if __name__ == "__main__":
    sys.exit(main())

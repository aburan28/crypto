"""Run ONE cryptanalysis experiments/ic-bench cell and print its receipt as JSON.

    python3 cryptanalysis_cell.py --root ../cryptanalysis --cell '{"n": 13, ...}' --limits '{...}'

bench.py itself runs only whole suites, through an unpinned process pool, and
derives run numbers from its history file; ``--record`` rewrites its
baseline.  This runner calls the same ``run_cell(cell, calibration)`` in this
one process (which ICMS has already pinned), with the frozen calibration,
and writes nothing into the cryptanalysis checkout except the C kernel that
``kernel.py`` builds on first use.  It refuses to run if bench.py's solver
limits are not the ones the spec pinned.
"""
from __future__ import annotations

import argparse
import json
import os
import sys


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--cell", required=True)
    ap.add_argument("--limits", required=True)
    args = ap.parse_args()
    bench_dir = os.path.join(os.path.abspath(args.root), "experiments", "ic-bench")
    sys.path.insert(0, bench_dir)
    import bench  # noqa: E402  (bench.py puts pdp-degree-heuristics on sys.path itself)

    want = json.loads(args.limits)
    if dict(bench.LIMITS) != want:
        print(json.dumps({"error": f"bench.LIMITS is {bench.LIMITS}, the spec pinned {want}"}))
        return 3
    cell = json.loads(args.cell)
    with open(bench.CALIBRATION, encoding="utf-8") as fh:
        calibration = json.load(fh)
    receipt = bench.run_cell(cell, calibration)
    receipt["icms_runner"] = {"calibration_id": calibration.get("calibration_id"),
                              "limits": dict(bench.LIMITS), "phases": list(bench.PHASES)}
    print(json.dumps(receipt, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())

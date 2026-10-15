#!/usr/bin/env python3
"""Secondary amortized S and descriptive wall panel from ecbench resources JSON.

Does not invent GAE prices. Amortized figures share factor_base+oracle_setup
across N targets and are labelled secondary. Wall figures are whatever the
session recorded (often L0 / descriptive).
"""
from __future__ import annotations

import argparse
import json
import statistics
from collections import defaultdict
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--records", type=Path, required=True)
    ap.add_argument("--ns", type=int, nargs="+", default=[1, 8, 64])
    args = ap.parse_args()
    by = defaultdict(list)
    for line in args.records.open():
        r = json.loads(line)
        if r.get("outcome", {}).get("status") != "verified":
            continue
        phases = {p["name"]: p for p in r["phases"]}
        setup = sum(
            phases.get(n, {}).get("gae", 0.0)
            for n in ("factor_base", "oracle_setup")
        )
        online = sum(
            phases.get(n, {}).get("gae", 0.0)
            for n in ("relations", "linear_algebra", "verify")
        )
        total = r["cost"]["total_gae"]
        s = r["cost"]["s"]
        sqrt_r = total / s if s else 1.0
        wall = r.get("time", {}).get("solve_wall_ns") or 0
        by[r["arm"]].append((setup, online, s, sqrt_r, wall))

    print(f"{'arm':32} {'coldS':>8}" + "".join(f"{'amS@'+str(n):>10}" for n in args.ns) + f"{'wall_ms':>10}")
    for arm in sorted(by):
        rows = by[arm]
        cold = statistics.mean(r[2] for r in rows)
        wall_ms = statistics.mean(r[4] for r in rows) / 1e6
        ams = []
        for n in args.ns:
            vals = [(su + n * on) / n / sr for su, on, _, sr, _ in rows]
            ams.append(statistics.mean(vals))
        print(
            f"{arm:32} {cold:8.3f}"
            + "".join(f"{v:10.3f}" for v in ams)
            + f"{wall_ms:10.2f}"
        )


if __name__ == "__main__":
    main()

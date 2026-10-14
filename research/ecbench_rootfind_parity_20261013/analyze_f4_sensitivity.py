#!/usr/bin/env python3
"""Secondary sensitivity: add pinned word-XOR GAE to F4 arms.

Does NOT change admitted S. Uses docs/ic/calibration.json ns_per_word_xor
with ns_per_add=1 (ecbench pinned_only_calibration). Labelled secondary.
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
    ap.add_argument(
        "--calibration",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "docs/ic/calibration.json",
    )
    ap.add_argument("--curve-key", default="koblitz/K_0 / GF(2^13)")
    args = ap.parse_args()
    cal = json.loads(args.calibration.read_text())["instances"][args.curve_key]
    ns_xor = float(cal["ns_per_word_xor"])
    by: dict[str, list[tuple[float, float, float]]] = defaultdict(list)
    for line in args.records.open():
        r = json.loads(line)
        if r.get("outcome", {}).get("status") != "verified":
            continue
        s = r["cost"]["s"]
        tot = r["cost"]["total_gae"]
        sr = tot / s if s else 1.0
        ops = (r.get("solver") or {}).get("ops") or 0
        xor_gae = ops * ns_xor  # ns_per_add = 1
        by[r["arm"]].append((s, (tot + xor_gae) / sr if ops else s, xor_gae))

    print(
        f"# F4 XOR sensitivity (secondary): ns_per_word_xor={ns_xor} "
        f"from {args.curve_key}; ns_per_add=1"
    )
    print(f"{'arm':32} {'coldS':>8} {'S+xor':>8} {'xorGAE':>10}")
    for arm in sorted(by):
        rows = by[arm]
        cold = statistics.mean(r[0] for r in rows)
        priced = statistics.mean(r[1] for r in rows)
        xg = statistics.mean(r[2] for r in rows)
        print(f"{arm:32} {cold:8.3f} {priced:8.3f} {xg:10.1f}")


if __name__ == "__main__":
    main()

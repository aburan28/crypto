#!/usr/bin/env python3
"""Reduce calibration_run{1,2,3}.json to the registered prices (calibration_summary.json).

Bounds, per RESEARCH_CACHEGRIND_N41_20260929.md:
  L1            median over runs of latency at 16 KiB
  P_mid_low     min over runs of latency at 512 KiB  - L1
  P_mid_high    max over runs of latency at 2 MiB    - L1
  P_far_low     min over runs of latency at 8 MiB    - L1
  P_far_high    max over runs of latency at 256 MiB  - L1
  MLP_max       max over runs of the best MLP over k
  branch_ns     median over runs of ns per mispredicted branch
"""
import glob
import json
import statistics

runs = [json.load(open(p)) for p in sorted(glob.glob("calibration_run[0-9]*.json"))]
KIB, MIB = 1024, 1024 * 1024
lat = lambda r, w: r["latency_ns_by_bytes"][str(w)]["ns_per_load_median_of_runs"]
l1 = statistics.median(lat(r, 16 * KIB) for r in runs)
out = {
    "n_runs": len(runs),
    "f_hz_median": statistics.median(r["f_hz"] for r in runs),
    "L1_ns": l1,
    "P_mid_low_ns": min(lat(r, 512 * KIB) for r in runs) - l1,
    "P_mid_high_ns": max(lat(r, 2 * MIB) for r in runs) - l1,
    "P_far_low_ns": min(lat(r, 8 * MIB) for r in runs) - l1,
    "P_far_high_ns": max(lat(r, 256 * MIB) for r in runs) - l1,
    "MLP_max": max(r["mlp_max"] for r in runs),
    "branch_ns": statistics.median(r["branch_ns_per_mispredict"] for r in runs),
    "latency_ns_by_bytes_min_max": {
        w: [min(lat(r, int(w)) for r in runs), max(lat(r, int(w)) for r in runs)]
        for w in runs[0]["latency_ns_by_bytes"]
    },
}
json.dump(out, open("calibration_summary.json", "w"), indent=2)
print(json.dumps({k: v for k, v in out.items() if k != "latency_ns_by_bytes_min_max"}, indent=1))
print(out["latency_ns_by_bytes_min_max"])

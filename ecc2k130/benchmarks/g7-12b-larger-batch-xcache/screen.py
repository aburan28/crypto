#!/usr/bin/env python3
import json
import math
from pathlib import Path
import statistics

import runbench

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
VARIANTS = {
    "selected": (ROOT / "build/ecc2k130-local-packed", 524286, 16),
    "batch24-cache4-min3": (ROOT / "build/g7-split-batch24", 349524, 24),
    "batch24-cache4-min2": (ROOT / "build/g7-batch24-cache4-min2", 349524, 24),
    "batch24-cache6-min2": (ROOT / "build/g7-batch24-cache6-min2", 349524, 24),
    "batch32-cache4-min3": (ROOT / "build/g7-split-batch32", 262143, 32),
    "batch32-cache4-min2": (ROOT / "build/g7-batch32-cache4-min2", 262143, 32),
    "batch32-cache8-min2": (ROOT / "build/g7-batch32-cache8-min2", 262143, 32),
}
ORDERS = [
    list(VARIANTS),
    list(reversed(VARIANTS)),
    ["batch24-cache6-min2", "batch32-cache4-min2", "selected",
     "batch24-cache4-min3", "batch32-cache8-min2",
     "batch24-cache4-min2", "batch32-cache4-min3"],
]

rows = []
for repetition, order in enumerate(ORDERS):
    for name in order:
        binary, workers, batch = VARIANTS[name]
        rows.append(runbench.bench(binary, f"timing-r{repetition}-{name}",
                                   workers=workers, batch=batch, run_id=62177))

rates = {name: [] for name in VARIANTS}
for row in rows:
    name = row["label"].split("-", 2)[2]
    rates[name].append(row["million_updates_per_second"] / 1000)

def paired(candidate, reference):
    ratios = [candidate[i] / reference[i] for i in range(3)]
    logs = [math.log(value) for value in ratios]
    mean = statistics.mean(logs)
    standard_error = statistics.stdev(logs) / math.sqrt(3)
    tcritical = 4.3026527297
    return {"paired_geomean_speedup": math.exp(mean),
            "paired_speedup_95_ci": [math.exp(mean - tcritical * standard_error),
                                      math.exp(mean + tcritical * standard_error)]}

summary = []
for name, values in rates.items():
    comparison = ({"paired_geomean_speedup": 1.0,
                   "paired_speedup_95_ci": [1.0, 1.0]}
                  if name == "selected" else paired(values, rates["selected"]))
    summary.append({"variant": name,
                    "median_billion_updates_per_second": statistics.median(values),
                    "rate_over_12b": statistics.median(values) / 12,
                    **comparison, "rates": values,
                    "decision": "reference" if name == "selected" else
                    ("screen qualifier" if comparison["paired_speedup_95_ci"][0] > 1
                     else "unconfirmed or regression")})

contrasts = {
    "batch24_cache6_over_cache4_min2": paired(
        rates["batch24-cache6-min2"], rates["batch24-cache4-min2"]),
    "batch32_cache8_over_cache4_min2": paired(
        rates["batch32-cache8-min2"], rates["batch32-cache4-min2"]),
    "batch24_min2_over_min3": paired(
        rates["batch24-cache4-min2"], rates["batch24-cache4-min3"]),
    "batch32_min2_over_min3": paired(
        rates["batch32-cache4-min2"], rates["batch32-cache4-min3"]),
}
result = {"passed": all(row["passed"] for row in rows),
          "unit": "billion_complete_scalar_updates_per_second",
          "scalar_walks_per_arm": 8388576,
          "updates_per_sample": 34359607296,
          "timing_samples": len(rows), "rows": summary,
          "contrasts": contrasts,
          "raw_labels": [row["label"] for row in rows],
          "goal_billion_updates_per_second": 12,
          "goal_met": max(statistics.median(values) for values in rates.values()) >= 12,
          "generic_work_boundary": "sqrt(n/262)",
          "generic_work_ratio": 1, "full_dlp_S": None}
(OUT / "comparison.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))

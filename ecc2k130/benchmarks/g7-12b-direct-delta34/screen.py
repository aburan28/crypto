#!/usr/bin/env python3
import json, math, pathlib, statistics

import runbench

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = pathlib.Path(__file__).resolve().parent
COLLISION_WORK_RATIO = 1.0130803718
VARIANTS = {
    "selected": ROOT / "build/ecc2k130-local-packed",
    "two-jump-indexed": ROOT / "build/g7-two-jump-indexed",
    "direct-delta34": ROOT / "build/g7-direct-delta34",
}
ORDERS = [
    ["selected", "two-jump-indexed", "direct-delta34"],
    ["direct-delta34", "two-jump-indexed", "selected"],
    ["two-jump-indexed", "selected", "direct-delta34"],
]

rows = []
for repetition, order in enumerate(ORDERS):
    for name in order:
        rows.append(
            runbench.bench(
                VARIANTS[name],
                f"timing-r{repetition}-{name}",
                run_id=57341,
            )
        )

rates = {name: [] for name in VARIANTS}
for row in rows:
    name = row["label"].split("-", 2)[2]
    rates[name].append(row["million_updates_per_second"] / 1000)

reference = rates["selected"]
t_critical_95_df2 = 4.3026527297
summary = []
for name, values in rates.items():
    paired = [values[i] / reference[i] for i in range(3)]
    logs = [math.log(value) for value in paired]
    mean = statistics.mean(logs)
    standard_error = (
        statistics.stdev(logs) / math.sqrt(3) if name != "selected" else 0
    )
    lower = math.exp(mean - t_critical_95_df2 * standard_error)
    upper = math.exp(mean + t_critical_95_df2 * standard_error)
    work_ratio = 1 if name == "selected" else COLLISION_WORK_RATIO
    adjusted = [value / work_ratio for value in values]
    adjusted_lower = lower / work_ratio
    adjusted_upper = upper / work_ratio
    summary.append(
        {
            "variant": name,
            "median_billion_updates_per_second": statistics.median(values),
            "rate_over_12b": statistics.median(values) / 12,
            "collision_work_ratio": work_ratio,
            "collision_adjusted_median_billion_iterations_per_second": statistics.median(adjusted),
            "paired_geomean_raw_speedup": math.exp(mean),
            "paired_raw_speedup_95_ci": [lower, upper],
            "paired_geomean_collision_adjusted_speedup": math.exp(mean) / work_ratio,
            "paired_collision_adjusted_speedup_95_ci": [adjusted_lower, adjusted_upper],
            "rates": values,
            "decision": (
                "reference"
                if name == "selected"
                else (
                    "screen qualifier"
                    if adjusted_lower > 1
                    else "unconfirmed or regression"
                )
            ),
        }
    )

out = {
    "passed": all(row["passed"] for row in rows),
    "unit": "billion_complete_scalar_updates_per_second",
    "updates_per_sample": 34_359_738_368,
    "timing_samples": len(rows),
    "rows": summary,
    "raw_labels": [row["label"] for row in rows],
    "goal_billion_updates_per_second": 12,
    "goal_met": max(statistics.median(values) for values in rates.values()) >= 12,
    "generic_work_boundary": "sqrt(n/262)",
    "generic_work_ratio": {
        "selected": 1,
        "two_jump": COLLISION_WORK_RATIO,
    },
    "full_dlp_S": None,
}
(OUT / "comparison.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))

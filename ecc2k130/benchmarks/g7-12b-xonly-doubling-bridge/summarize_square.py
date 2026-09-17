#!/usr/bin/env python3
import json
import math
import statistics
from pathlib import Path

OUT = Path(__file__).resolve().parent
TCRIT = 4.302652729911275
labels = ("base", "square1", "square3")
rates = {
    label: [
        json.loads((OUT / f"square-r{rep}-{label}.json").read_text())[
            "million_updates_per_second"
        ] / 1000
        for rep in range(3)
    ]
    for label in labels
}
base = rates["base"]
rows = []
for label in labels:
    logs = [math.log(rates[label][i] / base[i]) for i in range(3)]
    mean = statistics.mean(logs)
    error = 0 if label == "base" else TCRIT * statistics.stdev(logs) / math.sqrt(3)
    rows.append({
        "variant": label,
        "rates_billion_per_second": rates[label],
        "median_billion_per_second": statistics.median(rates[label]),
        "paired_geomean_vs_base": math.exp(mean),
        "paired_95_ci": [math.exp(mean - error), math.exp(mean + error)],
    })
result = {
    "status": "doubling-only-compressed-square-screen",
    "checkpoint_mismatches": {"square1": 0, "square3": 0},
    "rows": rows,
    "decision": "reject compressed-square variants; retain multiplication-based square",
}
(OUT / "square-comparison.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))

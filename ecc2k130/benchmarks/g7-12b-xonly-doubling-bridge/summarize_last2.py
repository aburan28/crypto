#!/usr/bin/env python3
import json
import math
import statistics
from pathlib import Path

OUT = Path(__file__).resolve().parent
TCRIT = 4.302652729911275
labels = ("base", "last2")
rates = {
    label: [
        json.loads((OUT / f"last2-r{rep}-{label}.json").read_text())[
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
    "status": "nonpromotable-common-kernel-last-prefix-diagnostic",
    "checkpoint_mismatches": {"last2": 0},
    "rows": rows,
    "decision": "reject last2 spill regression; retain base B16",
}
(OUT / "last2-comparison.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))

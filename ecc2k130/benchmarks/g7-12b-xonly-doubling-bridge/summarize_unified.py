#!/usr/bin/env python3
import json
import math
import statistics
from pathlib import Path

OUT = Path(__file__).resolve().parent
TCRIT = 4.302652729911275
labels = ("base", "complete")
rates = {
    label: [
        json.loads((OUT / f"unified-r{rep}-{label}.json").read_text())[
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
    "status": "complete-map-outlined-unified-diagnostic",
    "checkpoint_mismatches": {"complete": 0},
    "rows": rows,
    "decision": "reject outlined unified bridge: device-call stack traffic dominates",
}
(OUT / "unified-comparison.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))

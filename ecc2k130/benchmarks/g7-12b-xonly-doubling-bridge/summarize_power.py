#!/usr/bin/env python3
import json
import math
import statistics
from pathlib import Path

OUT = Path(__file__).resolve().parent
TCRIT = 4.302652729911275
labels = ("block", "power")
rates = {
    label: [
        json.loads((OUT / f"power-r{rep}-{label}.json").read_text())[
            "million_updates_per_second"
        ] / 1000
        for rep in range(3)
    ]
    for label in labels
}
base = rates["block"]
rows = []
for label in labels:
    logs = [math.log(rates[label][i] / base[i]) for i in range(3)]
    mean = statistics.mean(logs)
    error = 0 if label == "block" else TCRIT * statistics.stdev(logs) / math.sqrt(3)
    rows.append({
        "variant": label,
        "rates_billion_per_second": rates[label],
        "median_billion_per_second": statistics.median(rates[label]),
        "paired_geomean_vs_block": math.exp(mean),
        "paired_95_ci": [math.exp(mean - error), math.exp(mean + error)],
    })
result = {
    "status": "complete-map-bridge-power-chain-followup",
    "checkpoint_mismatches": {"power": 0},
    "rows": rows,
    "decision": "reject unconfirmed power-chain change; retain block-local reference",
}
(OUT / "power-comparison.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))

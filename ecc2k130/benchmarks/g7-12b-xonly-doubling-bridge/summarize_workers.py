#!/usr/bin/env python3
import json, math, statistics
from pathlib import Path
OUT = Path(__file__).resolve().parent
TCRIT = 4.302652729911275
labels = ("w131", "w262", "w524", "w1048")
rates = {label: [json.loads((OUT / f"workers-r{rep}-{label}.json").read_text())["million_updates_per_second"] / 1000 for rep in range(3)] for label in labels}
base = rates["w524"]
rows = []
for label in labels:
    logs = [math.log(rates[label][i] / base[i]) for i in range(3)]
    mean = statistics.mean(logs)
    error = 0 if label == "w524" else TCRIT * statistics.stdev(logs) / math.sqrt(3)
    rows.append({"variant": label, "rates_billion_per_second": rates[label], "median_billion_per_second": statistics.median(rates[label]), "paired_geomean_vs_w524": math.exp(mean), "paired_95_ci": [math.exp(mean-error), math.exp(mean+error)]})
result = {"status": "complete-map-worker-count-screen", "equal_scalar_updates_per_sample": 34359738368, "rows": rows, "decision": "retain 524288 workers; no alternative has a paired 95% interval above one"}
(OUT / "workers-comparison.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))

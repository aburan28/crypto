#!/usr/bin/env python3
import json, math, statistics
from pathlib import Path

OUT = Path(__file__).resolve().parent
TCRIT = 4.302652729911275


def load_rates(prefix, labels):
    rates = {}
    for label in labels:
        rates[label] = [
            json.loads((OUT / f"{prefix}-r{rep}-{label}.json").read_text())["million_updates_per_second"] / 1000
            for rep in range(3)
        ]
    return rates


def summarize(prefix, labels, base_label, status):
    rates = load_rates(prefix, labels)
    base = rates[base_label]
    rows = []
    for label in labels:
        logs = [math.log(rates[label][i] / base[i]) for i in range(3)]
        mean = statistics.mean(logs)
        error = 0 if label == base_label else TCRIT * statistics.stdev(logs) / math.sqrt(3)
        rows.append({
            "variant": label,
            "rates_billion_per_second": rates[label],
            "median_billion_per_second": statistics.median(rates[label]),
            f"paired_geomean_vs_{base_label}": math.exp(mean),
            "paired_95_ci": [math.exp(mean - error), math.exp(mean + error)],
            "rate_over_10b": statistics.median(rates[label]) / 10.0,
        })
    result = {
        "status": status,
        "collision_work_ratio_vs_selected": 1.0066958832826605,
        "generic_work_boundary": "sqrt(n/262) times collision work",
        "full_dlp_S": None,
        "rows": rows,
    }
    candidate = next(row for row in rows if row["variant"] != base_label)
    if candidate["paired_95_ci"][0] > 1:
        result["decision"] = "promote empty-queue skip; paired interval exceeds one"
        result["class"] = "engineering"
    else:
        result["decision"] = "retain modulus-72 pipeline; empty-queue skip not established"
        result["class"] = "engineering"
    return result


if __name__ == "__main__":
    bench = summarize(
        "skipempty",
        ("control", "skipempty"),
        "control",
        "complete-map-empty-queue-skip-benchmark-screen",
    )
    (OUT / "skipempty-comparison.json").write_text(json.dumps(bench, indent=2) + "\n")
    print(json.dumps(bench, indent=2), flush=True)
    dp34_files = list(OUT.glob("skipempty-dp34-r0-control.json"))
    if dp34_files:
        dp34 = summarize(
            "skipempty-dp34",
            ("control", "skipempty"),
            "control",
            "complete-map-empty-queue-skip-dp34-screen",
        )
        (OUT / "skipempty-dp34-comparison.json").write_text(json.dumps(dp34, indent=2) + "\n")
        print(json.dumps(dp34, indent=2), flush=True)

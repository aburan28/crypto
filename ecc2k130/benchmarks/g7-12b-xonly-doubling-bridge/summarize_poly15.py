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
        median = statistics.median(rates[label])
        rows.append({
            "variant": label,
            "rates_billion_per_second": rates[label],
            "median_billion_per_second": median,
            f"paired_geomean_vs_{base_label}": math.exp(mean),
            "paired_95_ci": [math.exp(mean - error), math.exp(mean + error)],
            "rate_over_10b": median / 10.0,
            "rate_over_12b": median / 12.0,
            "rate_over_15b": median / 15.0,
        })
    result = {
        "status": status,
        "generic_work_boundary": "sqrt(n/262) times collision work",
        "full_dlp_S": None,
        "rows": rows,
    }
    poly = next(row for row in rows if row["variant"] == "poly12")
    arith = next(row for row in rows if row["variant"] == "arith")
    if poly["median_billion_per_second"] >= 15.0 and poly["paired_95_ci"][0] > 1:
        result["decision"] = "promote poly12 complete map; median reaches 15 B/s and paired interval exceeds one"
        result["class"] = "engineering"
    elif arith["median_billion_per_second"] < 15.5:
        result["decision"] = "product+inverse ceiling below 15.5 B/s; 15 B/s is not available on this common path"
        result["class"] = "engineering"
    else:
        result["decision"] = "retain modulus-72 pipeline; poly12 did not establish 15 B/s"
        result["class"] = "engineering"
    return result


if __name__ == "__main__":
    bench = summarize(
        "poly15",
        ("control", "arith", "poly12"),
        "control",
        "complete-map-poly12-15b-benchmark-screen",
    )
    (OUT / "poly15-comparison.json").write_text(json.dumps(bench, indent=2) + "\n")
    print(json.dumps(bench, indent=2), flush=True)

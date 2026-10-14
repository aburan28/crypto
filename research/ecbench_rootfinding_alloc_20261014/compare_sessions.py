#!/usr/bin/env python3
"""Check deterministic ecbench records and describe paired L0 timings."""

import argparse
import json
import math
import random
import statistics
from pathlib import Path


def records(path):
    with (path / "records.jsonl").open() as handle:
        items = [json.loads(line) for line in handle]
    return {(x["arm"], x["workload"]["target_index"], x["round"]): x for x in items}


def deterministic(record):
    keys = (
        "arm", "round", "warmup", "role", "method", "workload", "algorithm_seed",
        "outcome", "cost", "boundaries", "counters", "factor_base",
        "detail", "solver", "field_ops",
    )
    result = {key: record.get(key) for key in keys}
    result["phases"] = [
        {key: value for key, value in phase.items() if key != "wall_ns"}
        for phase in record["phases"]
    ]
    return result


def cluster_ratio(a, b, arm, phase=None, draws=5000):
    def elapsed(record):
        if phase is None:
            return record["time"]["solve_wall_ns"]
        return next(row["wall_ns"] for row in record["phases"] if row["name"] == phase)

    target_logs = []
    for target in range(8):
        log_ratios = []
        for round_id in (1, 2, 3):
            key = arm, target, round_id
            old = elapsed(a[key])
            new = elapsed(b[key])
            log_ratios.append(math.log(new / old))
        target_logs.append(statistics.mean(log_ratios))
    rng = random.Random(20261014)
    samples = sorted(
        math.exp(statistics.mean(rng.choices(target_logs, k=8)))
        for _ in range(draws)
    )
    return {
        "candidate_over_baseline_geomean": math.exp(statistics.mean(target_logs)),
        "target_cluster_bootstrap_95_percent": [samples[125], samples[4874]],
        "target_log_ratios": target_logs,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    a, b = records(args.baseline), records(args.candidate)
    missing_a = sorted(set(b) - set(a))
    missing_b = sorted(set(a) - set(b))
    mismatches = []
    for key in sorted(set(a) & set(b)):
        av, bv = deterministic(a[key]), deterministic(b[key])
        for field in av:
            if av[field] != bv[field]:
                mismatches.append({"key": key, "field": field})
    result = {
        "baseline": str(args.baseline),
        "candidate": str(args.candidate),
        "records_baseline": len(a),
        "records_candidate": len(b),
        "missing_baseline": missing_a,
        "missing_candidate": missing_b,
        "deterministic_mismatches": mismatches,
        "all_verified": all(x["outcome"]["status"] == "verified" for x in list(a.values()) + list(b.values())),
    }
    if not (missing_a or missing_b or mismatches):
        result["l0_timing"] = {}
        for arm in ("rho-neg", "b48-serial-8", "b48-batch-8", "rho-neg-aa"):
            result["l0_timing"][arm] = {"solve_wall": cluster_ratio(a, b, arm)}
            if arm.startswith("b48-"):
                result["l0_timing"][arm]["relations_wall"] = cluster_ratio(
                    a, b, arm, phase="relations"
                )
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "l0_timing"}, indent=2))
    if missing_a or missing_b or mismatches or not result["all_verified"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

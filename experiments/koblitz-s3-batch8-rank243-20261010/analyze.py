"""Analyze the frozen paired single-target batch-inversion screen."""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import random
import statistics


HERE = Path(__file__).resolve().parent


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def geometric(values: list[float]) -> float:
    assert values and all(value > 0 for value in values)
    return math.exp(statistics.mean(math.log(value) for value in values))


def analyze(ledger: Path, output: Path) -> None:
    freeze_path = HERE / "FREEZE.json"
    freeze = json.loads(freeze_path.read_bytes())
    raw = ledger.read_bytes() if ledger.exists() else b""
    rows = [json.loads(line) for line in raw.splitlines()]
    target_map = {(target["n"], target["label"]): target
                  for target in freeze["selected_targets"]}
    expected = []
    for index in freeze["selected_indices"]:
        for n in ((41, 53) if index % 2 else (53, 41)):
            label = f"T{index:03}"
            arms = (("baseline", "batch8", "batch8", "baseline")
                    if (index + (n == 53)) % 2 else
                    ("batch8", "baseline", "baseline", "batch8"))
            repeats = defaultdict(int)
            for arm in arms:
                repeats[arm] += 1
                expected.append((n, label, arm, repeats[arm]))
    actual = [(row["n"], row["label"], row["arm"], row["repetition"])
              for row in rows]
    order_prefix_ok = actual == expected[:len(actual)]
    verified = []
    for row in rows:
        target = target_map[(row["n"], row["label"])]
        arm = row["arm"]
        key = f"n{row['n']}_{arm}"
        report = row["report"]
        valid = bool(row["verified_native_and_fixture"] and report and row["phase_ok"]
                     and row["candidate_id"] ==
                     freeze["candidate_manifests"][key]["candidate_id"]
                     and row["binary_sha256"] == freeze["binary_sha256"][arm]
                     and row["source_sha256"] == freeze["source_sha256"][arm]
                     and row["workload_id"] == target["workload_id"]
                     and row["point_sha256"] == target["point_sha256"]
                     and report["target"] == target["public_point"]
                     and report["recovered_scalar"] == target["fixture_scalar"]
                     and report["group_verified"] is True)
        verified.append(valid)
    complete = bool(order_prefix_ok and len(rows) == len(expected)
                    and all(verified))
    result = {
        "kind": "s3_rank243_one_worker_batch8_screen_analysis_v1",
        "freeze_sha256": sha(freeze_path.read_bytes()),
        "ledger_sha256": sha(raw), "planned_runs": len(expected),
        "observed_runs": len(rows), "verified_runs": sum(verified),
        "order_prefix_ok": order_prefix_ok, "complete": complete,
        "development_gate_per_curve": 2.0,
        "development_gate_passed": None,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
        "curves": {},
    }
    if complete:
        for n in (41, 53):
            curve_rows = [row for row in rows if row["n"] == n]
            target_results = []
            for index in freeze["selected_indices"]:
                label = f"T{index:03}"
                target_rows = [row for row in curve_rows if row["label"] == label]
                arms = {arm: [row for row in target_rows if row["arm"] == arm]
                        for arm in ("baseline", "batch8")}
                assert all(len(group) == 2 for group in arms.values())
                times = {arm: [row["report"]["timing_ms"]["target_online_after_reusable_setup"]
                               for row in group] for arm, group in arms.items()}
                pdp = {arm: [row["report"]["timing_ms"]["target_pdp_charged"]
                             for row in group] for arm, group in arms.items()}
                calls = {arm: [row["report"]["target_s3_calls"] for row in group]
                         for arm, group in arms.items()}
                target_results.append({
                    "label": label, "workload_id": target_map[(n, label)]["workload_id"],
                    "online_ms": times, "target_pdp_charged_ms": pdp,
                    "target_s3_calls": calls,
                    "online_ratio": geometric(times["baseline"]) /
                                    geometric(times["batch8"]),
                    "pdp_ratio": geometric(pdp["baseline"]) /
                                 geometric(pdp["batch8"]),
                })
            ratios = [item["online_ratio"] for item in target_results]
            rng = random.Random(20261010 + n)
            boot = sorted(geometric([ratios[rng.randrange(len(ratios))]
                                     for _ in ratios]) for _ in range(20000))
            arm_rows = {arm: [row for row in curve_rows if row["arm"] == arm]
                        for arm in ("baseline", "batch8")}
            aa_abs_log_p95 = {}
            for arm, group in arm_rows.items():
                by_target = defaultdict(list)
                for row in group:
                    by_target[row["label"]].append(
                        row["report"]["timing_ms"]["target_online_after_reusable_setup"])
                logs = sorted(abs(math.log(pair[0] / pair[1])) for pair in by_target.values())
                aa_abs_log_p95[arm] = logs[min(len(logs) - 1, math.ceil(.95 * len(logs)) - 1)]
            result["curves"][str(n)] = {
                "target_count": len(target_results),
                "online_baseline_over_batch8_geometric": geometric(ratios),
                "online_target_bootstrap_95pct": [boot[499], boot[19499]],
                "median_online_ms": {arm: statistics.median(
                    row["report"]["timing_ms"]["target_online_after_reusable_setup"]
                    for row in group) for arm, group in arm_rows.items()},
                "median_target_pdp_charged_ms": {arm: statistics.median(
                    row["report"]["timing_ms"]["target_pdp_charged"]
                    for row in group) for arm, group in arm_rows.items()},
                "aa_p95_absolute_log_online_ratio": aa_abs_log_p95,
                "targets": target_results,
            }
        result["development_gate_passed"] = all(
            curve["online_baseline_over_batch8_geometric"] >= 2.0
            for curve in result["curves"].values())
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"complete": complete,
                      "verified_runs": result["verified_runs"],
                      "ratios": {n: curve["online_baseline_over_batch8_geometric"]
                                 for n, curve in result["curves"].items()},
                      "development_gate_passed": result["development_gate_passed"]},
                     sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.ledger.resolve(), args.output.resolve())

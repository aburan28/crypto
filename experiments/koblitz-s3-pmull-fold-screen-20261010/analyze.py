"""Audit source-matched ARM PMULL folding rows and estimate target-paired costs."""

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


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def geometric(values: list[float]) -> float:
    assert values and all(value > 0 for value in values)
    return math.exp(statistics.mean(math.log(value) for value in values))


def analyze(ledger: Path, output: Path) -> dict:
    freeze_path = HERE / "FREEZE.json"
    freeze = json.loads(freeze_path.read_bytes())
    raw = ledger.read_bytes() if ledger.exists() else b""
    rows = [json.loads(line) for line in raw.splitlines()]
    actual = [(row["n"], row["label"], row["arm"], row["repetition"]) for row in rows]
    expected = [(row["n"], row["label"], row["arm"], row["repetition"])
                for row in freeze["run_order"]]
    order_prefix_ok = actual == expected[:len(actual)]
    targets = {(t["n"], t["label"]): t for t in freeze["selected_targets"]}
    audited = []
    for row in rows:
        n, label, arm = row["n"], row["label"], row["arm"]
        config = freeze["curves"][str(n)]
        frozen = config["arms"][arm]
        target = targets[(n, label)]
        stem = f"N{n}_{label}_{arm}_R1"
        report_path = ledger.parent / f"{stem}.json"
        stdout_path = ledger.parent / f"{stem}.stdout"
        stderr_path = ledger.parent / f"{stem}.stderr"
        report_raw = report_path.read_bytes() if report_path.exists() else b""
        stdout_raw = stdout_path.read_bytes() if stdout_path.exists() else b""
        stderr_raw = stderr_path.read_bytes() if stderr_path.exists() else b""
        try:
            report = json.loads(report_raw) if report_raw else None
        except ValueError:
            report = None
        timing = report.get("timing_ms", {}) if report else {}
        valid = bool(row["verified_native_and_fixture"] and row["phase_ok"] and report
                     and row["exit_code"] == 0 and row["error"] is None
                     and row["candidate_id"] == frozen["candidate_id"]
                     and row["run_id"] == f"{frozen['candidate_id']}W{target['workload_id']}R1"
                     and row["workload_id"] == target["workload_id"]
                     and row["point_sha256"] == target["point_sha256"]
                     and row["source_sha256"] == freeze["source_sha256"][frozen["binary_name"]]
                     and row["binary_sha256"] == freeze["binary_sha256"][frozen["binary_name"]]
                     and row["report_sha256"] == sha(report_raw)
                     and row["stdout_sha256"] == sha(stdout_raw)
                     and row["stderr_sha256"] == sha(stderr_raw)
                     and report_raw == stdout_raw
                     and report["group_verified"] is True
                     and report["target"] == target["public_point"]
                     and report["recovered_scalar"] == target["fixture_scalar"]
                     and report["factor_base_points"] == config["factor_base_points"]
                     and report["factor_base_digest"] == config["factor_base_digest"]
                     and report["prepared_rank"] == config["prepare_rank"]
                     and report["target_query_batch_width"] == 8
                     and row["gf2_arm_pmull_fold"] == freeze["arms_env"][arm]
                     and report["target_s3_calls"] == row["target_s3_calls"]
                     and timing["target_online_after_reusable_setup"] == row["online_ms"]
                     and timing["target_pdp_charged"] == row["target_pdp_charged_ms"])
        audited.append({"n": n, "label": label, "arm": arm, "valid": valid})
    complete = bool(order_prefix_ok and len(rows) == len(expected)
                    and all(item["valid"] for item in audited))
    result = {
        "kind": "s3_pmull_fold_screen_analysis_v1",
        "freeze_sha256": sha(freeze_path.read_bytes()),
        "ledger_sha256": sha(raw),
        "planned_runs": len(expected), "observed_runs": len(rows),
        "verified_runs": sum(item["valid"] for item in audited),
        "order_prefix_ok": order_prefix_ok, "complete": complete,
        "audited_rows": audited, "curves": {},
        "host_exclusive_partition_receipt": None, "controlled_speedup": None,
    }
    if complete:
        by_target = defaultdict(dict)
        for row in rows:
            by_target[(row["n"], row["label"])][row["arm"]] = row
        for n in (41, 53):
            pairs = []
            for target in (item for item in freeze["selected_targets"] if item["n"] == n):
                baseline = by_target[(n, target["label"])]["table"]
                candidate = by_target[(n, target["label"])]["fold"]
                pairs.append({
                    "label": target["label"], "workload_id": target["workload_id"],
                    "table_online_ms": baseline["online_ms"],
                    "fold_online_ms": candidate["online_ms"],
                    "table_target_pdp_ms": baseline["target_pdp_charged_ms"],
                    "fold_target_pdp_ms": candidate["target_pdp_charged_ms"],
                    "table_s3_calls": baseline["target_s3_calls"],
                    "fold_s3_calls": candidate["target_s3_calls"],
                    "online_ratio": baseline["online_ms"] / candidate["online_ms"],
                })
            ratios = [pair["online_ratio"] for pair in pairs]
            calls_equal = all(pair["table_s3_calls"] == pair["fold_s3_calls"] for pair in pairs)
            rng = random.Random(20261010 + n)
            boot = sorted(geometric([ratios[rng.randrange(len(ratios))] for _ in ratios])
                          for _ in range(20000))
            curve_rows = [row for row in rows if row["n"] == n]
            result["curves"][str(n)] = {
                "targets": pairs,
                "geometric_online_table_over_fold": geometric(ratios),
                "target_bootstrap_95pct": [boot[499], boot[19499]],
                "all_s3_calls_equal": calls_equal,
                "big_speedup_gate_at_least_twofold": calls_equal and geometric(ratios) >= 2.0,
                "median_online_ms": {
                    arm: statistics.median(row["online_ms"] for row in curve_rows if row["arm"] == arm)
                    for arm in ("table", "fold")},
                "median_target_pdp_ms": {
                    arm: statistics.median(row["target_pdp_charged_ms"] for row in curve_rows if row["arm"] == arm)
                    for arm in ("table", "fold")},
                "median_setup_ms": {
                    arm: statistics.median(row["reusable_setup_ms"] for row in curve_rows if row["arm"] == arm)
                    for arm in ("table", "fold")},
                "max_peak_rss_bytes": max(row["peak_rss_bytes"] for row in curve_rows),
                "major_page_faults": sum(row["resource_counters"]["major_page_faults"]
                                         for row in curve_rows),
            }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"complete": complete, "verified_runs": result["verified_runs"],
                      "ratios": {n: item["geometric_online_table_over_fold"]
                                 for n, item in result["curves"].items()}}, sort_keys=True))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.ledger.resolve(), args.output.resolve())

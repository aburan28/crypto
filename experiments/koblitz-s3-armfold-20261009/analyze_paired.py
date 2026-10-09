"""Analyze frozen N41/N53 ARM table versus PMULL folding target blocks."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import random
import statistics


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
BATCHED = ROOT / "experiments/koblitz-s3-batched-target-20261009/FREEZE.json"
ARMS = ("table", "fold")
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def geomean(values: list[float]) -> float:
    return math.exp(statistics.mean(map(math.log, values)))


def bootstrap_95(values: list[float], seed: int) -> list[float]:
    rng = random.Random(seed)
    distribution = sorted(geomean([values[rng.randrange(len(values))]
                                   for _ in values]) for _ in range(10_000))
    return [distribution[249], distribution[9749]]


def p95(values: list[float]) -> float:
    selected = sorted(values)
    return selected[math.ceil(.95 * len(selected)) - 1]


def main(ledger: Path, output: Path) -> None:
    raw = ledger.read_bytes()
    rows = [json.loads(line) for line in raw.splitlines()]
    plan_bytes = (HERE / "FREEZE.json").read_bytes()
    plan = json.loads(plan_bytes)
    parent = json.loads(PARENT.read_bytes())
    batched = json.loads(BATCHED.read_bytes())
    assert plan["parent_plan_sha256"] == sha(PARENT.read_bytes())
    assert plan["batched_plan_sha256"] == sha(BATCHED.read_bytes())
    assert batched["source_sha256"] == plan["example_source_sha256"]
    targets = {(item["n"], item["label"]): item for item in parent["targets"]}
    frozen = {(item["n"], item["label"]): item for item in plan["targets"]}
    by_target = defaultdict(list)
    seen_run_ids = set()
    for row in rows:
        n, label, arm = row["n"], row["label"], row["arm"]
        assert row["kind"] == "s3_arm_pmull_fold_paired_run_v1"
        assert (n, label) in frozen and arm in ARMS
        assert row["target_sha256"] == targets[n, label]["point_sha256"]
        assert row["target_sha256"] == frozen[n, label]["point_sha256"]
        assert row["example_source_sha256"] == plan["example_source_sha256"]
        assert row["field_arithmetic_source_sha256"] == plan["field_arithmetic_source_sha256"]
        assert row["binary_sha256"] == plan["local_binary_sha256"]
        assert row["GF2_ARM_PMULL_FOLD"] == ("1" if arm == "fold" else "unset")
        candidate = plan["curves"][str(n)][f"{arm}_candidate_id"]
        assert row["candidate_id"] == candidate
        assert row["run_id"] == f"{candidate}W{row['workload_id']}R{row['repetition']}"
        assert row["run_id"] not in seen_run_ids
        seen_run_ids.add(row["run_id"])
        by_target[n, label].append(row)
    assert len(rows) == 4 * len(frozen)
    assert set(by_target) == set(frozen)
    blocks = []
    for (n, label), selected in sorted(by_target.items()):
        assert len(selected) == 4
        assert Counter(row["arm"] for row in selected) == {arm: 2 for arm in ARMS}
        by_arm = {arm: sorted((row for row in selected if row["arm"] == arm),
                              key=lambda row: row["repetition"])
                  for arm in ARMS}
        successful = all(row["verified_native_and_fixture"] and row["phase_ok"]
                         for row in selected)
        block = {"n": n, "label": label,
                 "all_four_verified": successful,
                 "statuses": [{"arm": row["arm"], "repetition": row["repetition"],
                               "exit_code": row["exit_code"], "error": row["error"],
                               "verified": row["verified_native_and_fixture"]}
                              for row in selected]}
        if successful:
            time = {arm: [row["report"]["timing_ms"]["target_online_after_reusable_setup"]
                          for row in by_arm[arm]] for arm in ARMS}
            calls = {arm: [row["report"]["target_s3_calls"] for row in by_arm[arm]]
                     for arm in ARMS}
            relation = {arm: [row["report"]["target_relation_indices"]
                              for row in by_arm[arm]] for arm in ARMS}
            if n == 41:
                assert all(item == relation["table"][0]
                           for arm in ARMS for item in relation[arm])
                assert all(item == calls["table"][0]
                           for arm in ARMS for item in calls[arm])
            block.update({
                "online_ms_by_arm": time,
                "target_s3_calls_by_arm": calls,
                "target_relation_indices_by_arm": relation,
                "table_over_fold_online_ratio": (
                    geomean(time["table"]) / geomean(time["fold"])),
                "table_over_fold_s3_call_ratio": (
                    geomean(calls["table"]) / geomean(calls["fold"])),
                "aa_abs_log_ratio_by_arm": {
                    arm: abs(math.log(time[arm][0] / time[arm][1]))
                    for arm in ARMS}})
        blocks.append(block)
    curves = {}
    for n in (41, 53):
        selected = [block for block in blocks if block["n"] == n]
        complete = [block for block in selected if block["all_four_verified"]]
        ratios = [block["table_over_fold_online_ratio"] for block in complete]
        metrics = {"target_count": len(selected),
                   "verified_target_blocks": len(complete),
                   "run_count": sum(len(by_target[n, block["label"]]) for block in selected),
                   "controlled_speedup": None}
        if ratios:
            metrics.update({
                "exploratory_table_over_fold_online_geomean": geomean(ratios),
                "exploratory_target_bootstrap_95": bootstrap_95(ratios, 20261009 + n),
                "points_at_least_2x": sum(value >= 2 for value in ratios),
                "points_above_1x": sum(value > 1 for value in ratios),
                "table_over_fold_s3_call_geomean": geomean([
                    block["table_over_fold_s3_call_ratio"] for block in complete]),
                "aa_p95_abs_log_ratio_by_arm": {
                    arm: p95([block["aa_abs_log_ratio_by_arm"][arm]
                              for block in complete]) for arm in ARMS},
                "median_online_ms_by_arm": {
                    arm: statistics.median(value for block in complete
                                           for value in block["online_ms_by_arm"][arm])
                    for arm in ARMS},
                "median_target_s3_calls_by_arm": {
                    arm: statistics.median(value for block in complete
                                           for value in block["target_s3_calls_by_arm"][arm])
                    for arm in ARMS},
                "median_charged_phase_ms_by_arm": {
                    arm: {phase: statistics.median(
                        row["report"]["timing_ms"][phase]
                        for block in selected
                        for row in by_target[n, block["label"]]
                        if row["arm"] == arm and row["verified_native_and_fixture"])
                          for phase in PHASES} for arm in ARMS}})
        curves[str(n)] = metrics
    result = {"kind": "s3_arm_pmull_fold_paired_analysis_v1",
              "plan_sha256": sha(plan_bytes), "ledger_sha256": sha(raw),
              "run_count": len(rows), "all_native_verified": all(
                  row["verified_native_and_fixture"] for row in rows),
              "all_phase_sums_valid": all(row["phase_ok"] for row in rows),
              "controlled_speedup": None, "curves": curves, "blocks": blocks}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(output), "curves": curves}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    main(args.ledger, args.output)

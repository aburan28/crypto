"""Summarize the N53 delayed-worker parameter sweep by frozen target."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import statistics


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
PARALLEL = ROOT / "experiments/koblitz-s3-parallel-target-20261009/FREEZE.json"
ARMS = ("parallel", "hedged_0", "hedged_1000", "hedged_3000")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def geomean(values: list[float]) -> float:
    return math.exp(statistics.mean(map(math.log, values)))


def main(ledger: Path, output: Path) -> None:
    raw = ledger.read_bytes()
    rows = [json.loads(line) for line in raw.splitlines()]
    plan_bytes = (HERE / "FREEZE.json").read_bytes()
    plan = json.loads(plan_bytes)
    parent = json.loads(PARENT.read_bytes())
    parallel = json.loads(PARALLEL.read_bytes())
    assert plan["parent_plan_sha256"] == sha(PARENT.read_bytes())
    assert plan["parallel_plan_sha256"] == sha(PARALLEL.read_bytes())
    targets = {(item["n"], item["label"]): item for item in parent["targets"]}
    by_label = defaultdict(list)
    seen_run_ids = set()
    for row in rows:
        assert row["kind"] == "s3_hedged_delay_exploratory_probe_v1"
        assert row["n"] == 53 and row["arm"] in ARMS
        assert row["target_sha256"] == targets[53, row["label"]]["point_sha256"]
        assert row["source_sha256"] == plan["source_sha256"]
        assert row["hedged_binary_sha256"] == plan["local_binary_sha256"]
        expected_binary = (parallel["local_binary_sha256"] if row["delay_us"] is None
                           else plan["local_binary_sha256"])
        assert row["binary_sha256"] == expected_binary
        candidate = (plan["reference_candidate_id"] if row["delay_us"] is None
                     else plan["hedged_candidate_ids_by_delay_us"][str(row["delay_us"])])
        assert row["candidate_id"] == candidate
        assert row["run_id"] == f"{candidate}W{row['workload_id']}R{row['repetition']}"
        assert row["run_id"] not in seen_run_ids
        seen_run_ids.add(row["run_id"])
        assert row["verified_native_and_fixture"] and row["phase_ok"]
        by_label[row["label"]].append(row)
    blocks = []
    for label, selected in sorted(by_label.items()):
        assert len(selected) == 8
        assert Counter(row["arm"] for row in selected) == {arm: 2 for arm in ARMS}
        by_arm = {arm: sorted((row for row in selected if row["arm"] == arm),
                              key=lambda row: row["repetition"])
                  for arm in ARMS}
        time = {arm: [row["report"]["timing_ms"]["target_online_after_reusable_setup"]
                      for row in values] for arm, values in by_arm.items()}
        calls = {arm: [row["report"]["target_s3_calls"] for row in values]
                 for arm, values in by_arm.items()}
        block = {"label": label, "online_ms_by_arm": time,
                 "target_s3_calls_by_arm": calls,
                 "aa_abs_log_ratio_by_arm": {
                     arm: abs(math.log(values[0] / values[1]))
                     for arm, values in time.items()},
                 "parallel_over_hedged_online_ratio": {
                     arm: geomean(time["parallel"]) / geomean(time[arm])
                     for arm in ARMS if arm != "parallel"},
                 "parallel_over_hedged_root_call_ratio": {
                     arm: geomean(calls["parallel"]) / geomean(calls[arm])
                     for arm in ARMS if arm != "parallel"}}
        blocks.append(block)
    curves = {}
    baseline_aas = sorted(block["aa_abs_log_ratio_by_arm"]["parallel"]
                          for block in blocks)
    for arm in ARMS:
        if arm == "parallel":
            continue
        ratios = [block["parallel_over_hedged_online_ratio"][arm] for block in blocks]
        calls = [block["parallel_over_hedged_root_call_ratio"][arm] for block in blocks]
        aas = sorted(block["aa_abs_log_ratio_by_arm"][arm] for block in blocks)
        root_aas = sorted(abs(math.log(block["target_s3_calls_by_arm"][arm][0] /
                                        block["target_s3_calls_by_arm"][arm][1]))
                          for block in blocks)
        curves[arm] = {"target_count": len(blocks),
                       "exploratory_online_geomean_ratio": geomean(ratios),
                       "exploratory_online_min_ratio": min(ratios),
                       "online_ratio_above_one_count": sum(value > 1 for value in ratios),
                       "root_call_geomean_ratio": geomean(calls),
                       "aa_p95_abs_log_ratio": aas[math.ceil(.95 * len(aas)) - 1],
                       "root_call_aa_p95_abs_log_ratio": root_aas[
                           math.ceil(.95 * len(root_aas)) - 1]}
    result = {"kind": "s3_hedged_delay_tuning_analysis_v1",
              "plan_sha256": sha(plan_bytes), "ledger_sha256": sha(raw),
              "run_count": len(rows),
              "target_count": len(by_label), "controlled_speedup": None,
              "parallel_baseline_aa_p95_abs_log_ratio": baseline_aas[
                  math.ceil(.95 * len(baseline_aas)) - 1],
              "arms": curves, "blocks": blocks}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(output), "arms": curves}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    main(args.ledger, args.output)

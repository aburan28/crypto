"""Compare persistent tiles with per-tile spawning and first-winner search."""

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
TILED = ROOT / "experiments/koblitz-s3-tiled-target-20261009/FREEZE.json"
ARMS = ("parallel", "tiled", "persistent")


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
    tiled = json.loads(TILED.read_bytes())
    assert plan["parent_plan_sha256"] == sha(PARENT.read_bytes())
    assert plan["parallel_plan_sha256"] == sha(PARALLEL.read_bytes())
    assert plan["tiled_plan_sha256"] == sha(TILED.read_bytes())
    targets = {(item["n"], item["label"]): item for item in parent["targets"]}
    by_label = defaultdict(list)
    seen_run_ids = set()
    for row in rows:
        assert row["kind"] == "s3_persistent_tiled_paired_run_v1"
        assert row["n"] == 53 and row["arm"] in ARMS
        assert row["target_sha256"] == targets[53, row["label"]]["point_sha256"]
        assert row["source_sha256"] == plan["source_sha256"]
        assert row["persistent_binary_sha256"] == plan["local_binary_sha256"]
        expected_binary = {"parallel": parallel["local_binary_sha256"],
                           "tiled": tiled["local_binary_sha256"],
                           "persistent": plan["local_binary_sha256"]}[row["arm"]]
        assert row["binary_sha256"] == expected_binary
        candidate = plan["candidate_ids_by_arm"][row["arm"]]
        assert row["candidate_id"] == candidate
        assert row["run_id"] == f"{candidate}W{row['workload_id']}R{row['repetition']}"
        assert row["run_id"] not in seen_run_ids
        seen_run_ids.add(row["run_id"])
        assert row["verified_native_and_fixture"] and row["phase_ok"]
        by_label[row["label"]].append(row)
    blocks = []
    for label, selected in sorted(by_label.items()):
        assert len(selected) == 6
        assert Counter(row["arm"] for row in selected) == {arm: 2 for arm in ARMS}
        by_arm = {arm: sorted((row for row in selected if row["arm"] == arm),
                              key=lambda row: row["repetition"])
                  for arm in ARMS}
        time = {arm: [row["report"]["timing_ms"]["target_online_after_reusable_setup"]
                      for row in values] for arm, values in by_arm.items()}
        calls = {arm: [row["report"]["target_s3_calls"] for row in values]
                 for arm, values in by_arm.items()}
        tiled_relations = [row["report"]["target_relation_indices"]
                           for arm in ("tiled", "persistent") for row in by_arm[arm]]
        assert all(item == tiled_relations[0] for item in tiled_relations)
        assert calls["tiled"] == calls["persistent"]
        assert calls["tiled"][0] == calls["tiled"][1]
        block = {"label": label, "online_ms_by_arm": time,
                 "target_s3_calls_by_arm": calls,
                 "same_first_relation_and_root_calls_for_tile_implementations": True,
                 "aa_abs_log_ratio_by_arm": {
                     arm: abs(math.log(values[0] / values[1]))
                     for arm, values in time.items()},
                 "tiled_over_persistent_online_ratio": (
                     geomean(time["tiled"]) / geomean(time["persistent"])),
                 "parallel_over_persistent_online_ratio": (
                     geomean(time["parallel"]) / geomean(time["persistent"])),
                 "parallel_over_persistent_root_call_ratio": (
                     geomean(calls["parallel"]) / geomean(calls["persistent"]))}
        blocks.append(block)
    comparisons = {}
    baseline_aas = sorted(block["aa_abs_log_ratio_by_arm"]["parallel"]
                          for block in blocks)
    for name in ("tiled_over_persistent_online_ratio",
                 "parallel_over_persistent_online_ratio"):
        ratios = [block[name] for block in blocks]
        comparisons[name] = {"exploratory_geomean": geomean(ratios),
                             "minimum": min(ratios),
                             "count_above_one": sum(value > 1 for value in ratios)}
    aas = {arm: sorted(block["aa_abs_log_ratio_by_arm"][arm] for block in blocks)
           for arm in ARMS}
    result = {"kind": "s3_persistent_tiled_paired_analysis_v1",
              "plan_sha256": sha(plan_bytes), "ledger_sha256": sha(raw),
              "run_count": len(rows),
              "target_count": len(by_label), "controlled_speedup": None,
              "parallel_baseline_aa_p95_abs_log_ratio": baseline_aas[
                  math.ceil(.95 * len(baseline_aas)) - 1],
              "aa_p95_abs_log_ratio_by_arm": {
                  arm: values[math.ceil(.95 * len(values)) - 1]
                  for arm, values in aas.items()},
              "all_tile_relations_and_root_counts_match": True,
              "comparisons": comparisons, "blocks": blocks}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(output), "comparisons": comparisons}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    main(args.ledger, args.output)

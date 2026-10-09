"""Analyze verified one-target batch4/paired S3 comparison blocks."""

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


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def geomean(values: list[float]) -> float | None:
    return math.exp(statistics.mean(map(math.log, values))) if values else None


def analyze(ledger: Path, output: Path) -> dict:
    plan_bytes = (HERE / "FREEZE.json").read_bytes()
    plan = json.loads(plan_bytes)
    parent_bytes = PARENT.read_bytes()
    parent = json.loads(parent_bytes)
    assert plan["parent_plan_sha256"] == sha(parent_bytes)
    assert plan["source_sha256"] == sha(
        (ROOT / "examples/s3_partial_prepared_batch4_targetspan.rs").read_bytes())
    target_by_key = {(item["n"], item["label"]): item for item in parent["targets"]}
    ledger_bytes = ledger.read_bytes()
    runs = [json.loads(line) for line in ledger_bytes.splitlines()]
    assert runs and len(runs) % 8 == 0
    by_target = defaultdict(list)
    seen_run_ids = set()
    for row in runs:
        n, label, arm = row["n"], row["label"], row["arm"]
        assert row["kind"] == "s3_batch4_single_target_run_v1"
        assert n in (41, 53) and arm in ("prepared", "batch4")
        target = target_by_key[(n, label)]
        assert 1 <= int(label[1:]) <= 24
        assert row["target_sha256"] == target["point_sha256"]
        assert row["workload_id"] == target["workload_id"]
        cid = plan["curves"][str(n)][f"{arm}_candidate_id"]
        assert row["candidate_id"] == cid
        assert row["run_id"].startswith(f"{cid}W{target['workload_id']}R")
        assert row["run_id"] not in seen_run_ids
        seen_run_ids.add(row["run_id"])
        expected_binary = (parent["partial_binary_sha256_local"] if arm == "prepared"
                           else plan["local_binary_sha256"])
        assert row["binary_sha256"] == expected_binary
        if row["verified_native_and_fixture"]:
            report = row["report"]
            assert row["exit_code"] == 0 and row["phase_ok"]
            assert report["group_verified"] and report["target_count"] == 1
            assert report["target"] == target["public_point"]
            assert report["recovered_scalar"] == target["fixture_scalar"]
            assert report["prepared_rank_reached"] and report["prepared_rank"] == 243
        by_target[n, label].append(row)
    labels = {int(label[1:]) for _, label in by_target}
    assert len(by_target) == 2 * len(labels)
    assert set(by_target) == {(n, f"T{index:03}") for index in labels for n in (41, 53)}
    blocks = []
    for (n, label), rows in sorted(by_target.items()):
        assert len(rows) == 4 and Counter(row["arm"] for row in rows) == {
            "prepared": 2, "batch4": 2}
        assert {row["run_id"] for row in rows} == {
            f"{plan['curves'][str(n)][f'{arm}_candidate_id']}W{target_by_key[n, label]['workload_id']}R{rep}"
            for arm in ("prepared", "batch4") for rep in (1, 2)}
        good = all(row["verified_native_and_fixture"] and row["phase_ok"] for row in rows)
        block = {"n": n, "label": label, "all_verified": good,
                 "run_ids": [row["run_id"] for row in rows],
                 "failed_run_ids": [row["run_id"] for row in rows
                                    if not row["verified_native_and_fixture"]],
                 "online_ratio_prepared_over_batch4": None,
                 "target_pdp_ratio_prepared_over_batch4": None,
                 "cold_ratio_prepared_over_batch4": None,
                 "root_call_ratio_prepared_over_batch4": None,
                 "prepared_aa_abs_log_ratio": None,
                 "batch4_aa_abs_log_ratio": None,
                 "relation_witness_equations_match": None}
        if good:
            arm_rows = {arm: [row for row in rows if row["arm"] == arm]
                        for arm in ("prepared", "batch4")}
            reports = {arm: [row["report"] for row in selected]
                       for arm, selected in arm_rows.items()}
            first_witnesses = {
                arm: [[(w["relation_scalar"], w["point_indices"])
                       for w in report["rank_relation_witnesses"]]
                      for report in values]
                for arm, values in reports.items()}
            matching = (all(w == first_witnesses["prepared"][0]
                            for arm in first_witnesses for w in first_witnesses[arm])
                        and len({tuple(report["target_relation_indices"])
                                 for values in reports.values() for report in values}) == 1)
            assert matching, f"first-witness policy changed for n{n} {label}"
            metric = {}
            for name, key in (("online", "target_online_after_reusable_setup"),
                              ("target_pdp", "target_pdp"),
                              ("cold", "cold_after_launch")):
                values = {arm: [report["timing_ms"][key] for report in selected]
                          for arm, selected in reports.items()}
                assert all(value > 0 for pair in values.values() for value in pair)
                metric[name] = math.sqrt(math.prod(values["prepared"]) /
                                         math.prod(values["batch4"]))
            calls = {arm: [report["target_s3_calls"] + report["fallback_s3_calls"]
                           for report in selected] for arm, selected in reports.items()}
            assert all(value > 0 for pair in calls.values() for value in pair)
            online_times = {arm: [report["timing_ms"]["target_online_after_reusable_setup"]
                                  for report in selected]
                            for arm, selected in reports.items()}
            block.update({
                "online_ratio_prepared_over_batch4": metric["online"],
                "target_pdp_ratio_prepared_over_batch4": metric["target_pdp"],
                "cold_ratio_prepared_over_batch4": metric["cold"],
                "root_call_ratio_prepared_over_batch4": (
                    math.sqrt(math.prod(calls["prepared"]) / math.prod(calls["batch4"]))),
                "prepared_aa_abs_log_ratio": abs(math.log(
                    online_times["prepared"][0] / online_times["prepared"][1])),
                "batch4_aa_abs_log_ratio": abs(math.log(
                    online_times["batch4"][0] / online_times["batch4"][1])),
                "relation_witness_equations_match": True,
                "online_ms_by_arm": online_times,
                "target_s3_calls_by_arm": calls,
            })
        blocks.append(block)
    curves = {}
    for n in (41, 53):
        selected = [block for block in blocks if block["n"] == n]
        valid = [block for block in selected if block["all_verified"]]
        ratios = [block["online_ratio_prepared_over_batch4"] for block in valid]
        curves[str(n)] = {"target_count": len(selected),
                          "verified_target_count": len(valid),
                          "failed_run_count": sum(len(block["failed_run_ids"])
                                                  for block in selected),
                          "online_ratio_geomean_exploratory": geomean(ratios),
                          "online_ratio_min_exploratory": min(ratios) if ratios else None,
                          "at_least_2_online_ratio_count": sum(value >= 2 for value in ratios),
                          "target_pdp_ratio_geomean_exploratory": geomean([
                              block["target_pdp_ratio_prepared_over_batch4"] for block in valid]),
                          "cold_ratio_geomean_exploratory": geomean([
                              block["cold_ratio_prepared_over_batch4"] for block in valid]),
                          "root_call_ratio_geomean": geomean([
                              block["root_call_ratio_prepared_over_batch4"] for block in valid]),
                          "prepared_aa_p95_abs_log_ratio": (
                              sorted(block["prepared_aa_abs_log_ratio"] for block in valid)[
                                  math.ceil(.95 * len(valid)) - 1] if valid else None),
                          "batch4_aa_p95_abs_log_ratio": (
                              sorted(block["batch4_aa_abs_log_ratio"] for block in valid)[
                                  math.ceil(.95 * len(valid)) - 1] if valid else None)}
    result = {"kind": "s3_batch4_paired_pilot_analysis_v1",
              "plan_sha256": sha(plan_bytes), "ledger_sha256": sha(ledger_bytes),
              "run_count": len(runs), "matched_index_count": len(labels),
              "strict_host_isolation_receipt": None, "controlled_speedup": None,
              "curves": curves, "blocks": blocks}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(output), "curves": curves}, sort_keys=True))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.ledger, args.output)

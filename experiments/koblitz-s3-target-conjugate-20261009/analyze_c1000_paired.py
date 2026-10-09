"""Validate the 396-column baseline and 1000-column conjugate panel."""

from __future__ import annotations

import argparse
from collections import defaultdict
import gzip
import hashlib
import json
import math
from pathlib import Path
import random
import statistics


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical(value: dict) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def geomean(values: list[float]) -> float:
    assert values and all(value > 0 for value in values)
    return math.exp(statistics.mean(map(math.log, values)))


def bootstrap_95(values: list[float]) -> list[float]:
    rng = random.Random(20261009)
    samples = sorted(geomean(rng.choices(values, k=len(values)))
                     for _ in range(20_000))
    return [samples[499], samples[19_499]]


def p95(values: list[float]) -> float:
    ordered = sorted(values)
    return ordered[math.ceil(.95 * len(ordered)) - 1]


def analyze(ledger_path: Path, output_path: Path) -> dict:
    if output_path.exists():
        raise FileExistsError(output_path)
    plan_bytes = (HERE / "C1000_FREEZE.json").read_bytes()
    plan = json.loads(plan_bytes)
    parent_bytes = PARENT.read_bytes()
    assert plan["kind"] == "s3_base1000_oneworker_paired_v1"
    assert plan["parent_plan_sha256"] == sha(parent_bytes)
    base1000_compressed = (HERE / "local/n53_base_c1000.json.gz").read_bytes()
    assert plan["base1000_raw_sha256"] == sha(gzip.decompress(base1000_compressed))
    assert {arm: plan["arms"][arm]["columns"] for arm in ("baseline", "expanded")} == {
        "baseline": 396, "expanded": 1000}
    assert plan["target_workers"] == 1 and plan["relation_workers"] == 14
    assert plan["target_labels"] == [f"T{i:03}" for i in range(25, 33)]
    assert plan["repetitions_per_arm_per_target"] == 2
    assert plan["controlled_speedup"] is None
    targets = {(row["n"], row["label"]): row
               for row in json.loads(parent_bytes)["targets"]}
    candidates = {}
    for arm in ("baseline", "expanded"):
        frozen = plan["arms"][arm]
        candidate = frozen["candidate_id"]
        path = HERE / "identities/candidates" / f"{candidate}.json"
        record = json.loads(path.read_bytes())
        assert record.pop("candidate_id") == candidate
        assert sha(canonical(record))[:12] == candidate[-12:]
        source = ROOT / "examples" / (
            "s3_partial_prepared_parallel_targetspan.rs" if arm == "baseline"
            else "s3_partial_prepared_target_conjugate_targetspan.rs")
        assert frozen["source_sha256"] == sha(source.read_bytes())
        assert record["point_decomposition"]["root_table_value_bits"] == 64
        assert record["factor_base"]["actual_usable_points_B"] == frozen["factor_base_points"]
        assert record["factor_base"]["effective_relation_columns"] == frozen["columns"]
        assert record["factor_base"]["enumerated_set_digest"] == frozen["factor_base_digest"]
        assert record["point_decomposition"]["root_table_slots"] == frozen["root_table_slots"]
        assert record["implementation"]["rust_library_source_tree_git_oid"] == plan["source_tree_git_oid"]
        candidates[arm] = record
    payload = ledger_path.read_bytes()
    raw = gzip.decompress(payload) if payload[:2] == b"\x1f\x8b" else payload
    rows = [json.loads(line) for line in raw.splitlines()]
    assert len(rows) == 32
    blocks = defaultdict(list)
    failures = []
    trace_hashes = set()
    for row in rows:
        assert row["kind"] == "s3_base1000_oneworker_paired_run_v1"
        assert row["n"] == 53 and row["arm"] in candidates
        arm = row["arm"]
        assert row["label"] in plan["target_labels"]
        target = targets[(53, row["label"])]
        frozen = plan["arms"][arm]
        candidate = frozen["candidate_id"]
        assert row["candidate_id"] == candidate
        assert row["run_id"] == f"{candidate}W{target['workload_id']}R{row['repetition']}"
        assert row["workload_id"] == target["workload_id"]
        assert row["point_sha256"] == target["point_sha256"]
        assert row["source_sha256"] == frozen["source_sha256"]
        assert row["binary_sha256"] == frozen["binary_sha256"]
        assert row["argv"][3] == str(frozen["columns"])
        assert row["argv"][-3:] == ["14", str(frozen["requested_pre_target_rank"]), "1"]
        report = row["report"]
        if row["error"] or row["exit_code"] != 0 or report is None or not row["verified_native_and_fixture"]:
            failures.append({"label": row["label"], "arm": arm,
                             "repetition": row["repetition"],
                             "error": row["error"], "exit_code": row["exit_code"]})
        else:
            times = report["timing_ms"]
            assert all(isinstance(times[k], (int, float)) and times[k] >= 0 for k in PHASES)
            assert abs(sum(times[k] for k in PHASES) -
                       times["target_online_after_reusable_setup"]) < 0.02
            assert row["phase_ok"] and row["verified_native_and_fixture"]
            assert report["group_verified"] and report["target_count"] == 1
            assert report["target"] == target["public_point"]
            assert report["recovered_scalar"] == target["fixture_scalar"]
            assert report["factor_base_points"] == frozen["factor_base_points"]
            assert report["factor_base_digest"] == frozen["factor_base_digest"]
            assert report["orbit_columns"] == frozen["columns"]
            assert report["root_table_slots"] == frozen["root_table_slots"]
            assert frozen["requested_pre_target_rank"] <= report["prepared_rank"] <= frozen["columns"]
            assert report["target_query_workers"] == 1
            if arm == "expanded":
                assert report["target_query_batch_width"] == 8
                assert report["target_conjugate_count"] == 53
                assert report["index_precomputed_poly_point_bytes"] == (
                    frozen["indexed_state_extra_poly_point_bytes"])
                assert report["target_query_orientation_policy"] == (
                    "precomputed_inverse_frobenius_target_conjugates_with_unrotated_index_points")
            trace_hashes.add((arm, sha(canonical(report["rank_relation_witnesses"]))))
        blocks[(row["label"], arm)].append(row)
    assert len(blocks) == 16
    assert all(len(block) == 2 and {row["repetition"] for row in block} == {1, 2}
               for block in blocks.values())
    summary = {"kind": "s3_base1000_oneworker_paired_analysis_v1",
               "freeze_sha256": sha(plan_bytes), "raw_ledger_sha256": sha(raw),
               "run_count": len(rows), "failure_rows": failures,
               "all_native_verified": not failures and all(row["verified_native_and_fixture"] for row in rows),
               "distinct_prepared_relation_traces": len(trace_hashes),
               "controlled_speedup": None,
               "factor_base_is_controlled_variable": True,
               "expanded_index_extra_poly_point_bytes":
                   plan["arms"]["expanded"]["indexed_state_extra_poly_point_bytes"],
               "root_table_payload_bytes": {
                   arm: 16 * plan["arms"][arm]["root_table_slots"]
                   for arm in ("baseline", "expanded")}}
    if not failures:
        target_blocks = {}
        for label in plan["target_labels"]:
            target_blocks[label] = {}
            for arm in ("baseline", "expanded"):
                block = sorted(blocks[(label, arm)], key=lambda row: row["repetition"])
                target_blocks[label][arm] = {
                    "online_ms": [row["report"]["timing_ms"]["target_online_after_reusable_setup"] for row in block],
                    "target_s3_calls": [row["report"]["target_s3_calls"] for row in block],
                }
        summary["target_blocks"] = target_blocks
        online_ratios, calls_ratios = [], []
        for block in target_blocks.values():
            online_ratios.append(geomean(block["baseline"]["online_ms"]) /
                                 geomean(block["expanded"]["online_ms"]))
            calls_ratios.append(geomean(block["baseline"]["target_s3_calls"]) /
                                geomean(block["expanded"]["target_s3_calls"]))
        summary["exploratory_baseline_over_expanded_online_geomean"] = geomean(online_ratios)
        summary["target_bootstrap_95"] = bootstrap_95(online_ratios)
        summary["points_at_least_2x"] = sum(value >= 2 for value in online_ratios)
        summary["points_above_1x"] = sum(value > 1 for value in online_ratios)
        summary["baseline_over_expanded_s3_call_geomean"] = geomean(calls_ratios)
        summary["arms"] = {}
        for arm in ("baseline", "expanded"):
            selected = [row for row in rows if row["arm"] == arm]
            times = [row["report"]["timing_ms"] for row in selected]
            aa = [abs(math.log(block[arm]["online_ms"][0] / block[arm]["online_ms"][1]))
                  for block in target_blocks.values()]
            summary["arms"][arm] = {
                "candidate_id": plan["arms"][arm]["candidate_id"],
                "factor_base_points": plan["arms"][arm]["factor_base_points"],
                "folded_columns": plan["arms"][arm]["columns"],
                "root_table_entries": sorted({row["report"]["root_table_entries"] for row in selected}),
                "root_table_slots": plan["arms"][arm]["root_table_slots"],
                "actual_prepared_ranks": sorted({row["report"]["prepared_rank"] for row in selected}),
                "fallback_attempts": sum(row["report"]["fallback_attempts"] for row in selected),
                "median_rank_attempts": statistics.median(row["report"]["rank_attempts"] for row in selected),
                "median_rank_verified_relations": statistics.median(row["report"]["rank_verified_relations"] for row in selected),
                "median_rank_novel_rows": statistics.median(row["report"]["rank_new_rows"] for row in selected),
                "median_online_ms": statistics.median(t["target_online_after_reusable_setup"] for t in times),
                "median_charged_phase_ms": {key: statistics.median(t[key] for t in times) for key in PHASES},
                "median_target_s3_calls": statistics.median(row["report"]["target_s3_calls"] for row in selected),
                "median_reusable_setup_ms_supplementary": statistics.median(t["reusable_setup_total"] for t in times),
                "median_setup_phase_ms_supplementary": {
                    key: statistics.median(t[key] for t in times)
                    for key in ("factor_base_construct", "index_build", "rank_collection_total",
                                "rank_matrix_build", "rank_linear_algebra_incremental",
                                "rank_linear_algebra_final")},
                "peak_rss_bytes": [row["report"]["peak_rss_bytes"] for row in selected],
                "aa_p95_abs_log_online_ratio": p95(aa),
            }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, default=HERE / "local/raw_paired32.jsonl.gz")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.ledger, args.output)
    print(json.dumps({"runs": result["run_count"], "failures": len(result["failure_rows"]),
                      "online_ratio": result.get("exploratory_baseline_over_expanded_online_geomean"),
                      "bootstrap_95": result.get("target_bootstrap_95")}, sort_keys=True))

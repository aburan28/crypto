"""Validate the target-conjugate S3 panel and report paired target costs."""

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
BASE = ROOT / "experiments/koblitz-s3-base-size-20261009/FREEZE.json"
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
    plan_bytes = (HERE / "FREEZE.json").read_bytes()
    plan = json.loads(plan_bytes)
    base_bytes = BASE.read_bytes()
    parent_bytes = PARENT.read_bytes()
    assert plan["kind"] == "s3_target_conjugate_paired_v1"
    assert plan["base_plan_sha256"] == sha(base_bytes)
    assert plan["parent_plan_sha256"] == sha(parent_bytes)
    assert plan["factor_base_columns"] == 396
    assert plan["factor_base_points"] == 41976
    assert plan["root_table_slots"] == 33554432
    assert plan["target_workers"] == plan["relation_workers"] == 14
    assert plan["requested_pre_target_rank"] == 395
    assert plan["target_labels"] == [f"T{i:03}" for i in range(25, 33)]
    assert plan["repetitions_per_arm_per_target"] == 2
    assert plan["controlled_speedup"] is None
    targets = {(row["n"], row["label"]): row
               for row in json.loads(parent_bytes)["targets"]}
    candidates = {}
    for arm in ("baseline", "conjugate"):
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
        assert row["kind"] == "s3_conjugate_target_paired_run_v1"
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
        assert row["argv"][3] == "396" and row["argv"][-3:] == ["14", "395", "14"]
        report = row["report"]
        if row["error"] or row["exit_code"] != 0 or report is None:
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
            assert report["factor_base_points"] == plan["factor_base_points"]
            assert report["factor_base_digest"] == plan["factor_base_digest"]
            assert report["root_table_slots"] == plan["root_table_slots"]
            assert report["prepared_rank"] == 396 and report["target_query_workers"] == 14
            assert report["fallback_attempts"] == 0
            if arm == "conjugate":
                assert report["target_query_batch_width"] == 8
                assert report["target_query_inversion_policy"] == (
                    "one normal-basis field inversion per eight regular state pair products")
                assert report["target_conjugate_count"] == 53
                assert report["index_precomputed_poly_point_bytes"] == (
                    plan["indexed_state_extra_poly_point_bytes"])
                assert report["target_query_orientation_policy"] == (
                    "precomputed_inverse_frobenius_target_conjugates_with_unrotated_index_points")
            trace_hashes.add(sha(canonical(report["rank_relation_witnesses"])))
        blocks[(row["label"], arm)].append(row)
    assert len(blocks) == 16
    assert all(len(block) == 2 and {row["repetition"] for row in block} == {1, 2}
               for block in blocks.values())
    summary = {"kind": "s3_conjugate_target_paired_analysis_v1",
               "freeze_sha256": sha(plan_bytes), "raw_ledger_sha256": sha(raw),
               "run_count": len(rows), "failure_rows": failures,
               "all_native_verified": not failures and all(row["verified_native_and_fixture"] for row in rows),
               "distinct_prepared_relation_traces": len(trace_hashes),
               "controlled_speedup": None,
               "indexed_state_extra_poly_point_bytes": plan["indexed_state_extra_poly_point_bytes"],
               "root_table_payload_bytes": {"baseline": 16 * plan["root_table_slots"],
                                            "conjugate": 16 * plan["root_table_slots"]}}
    if not failures:
        target_blocks = {}
        for label in plan["target_labels"]:
            target_blocks[label] = {}
            for arm in ("baseline", "conjugate"):
                block = sorted(blocks[(label, arm)], key=lambda row: row["repetition"])
                target_blocks[label][arm] = {
                    "online_ms": [row["report"]["timing_ms"]["target_online_after_reusable_setup"] for row in block],
                    "target_s3_calls": [row["report"]["target_s3_calls"] for row in block],
                }
        summary["target_blocks"] = target_blocks
        online_ratios, calls_ratios = [], []
        for block in target_blocks.values():
            online_ratios.append(geomean(block["baseline"]["online_ms"]) /
                                 geomean(block["conjugate"]["online_ms"]))
            calls_ratios.append(geomean(block["baseline"]["target_s3_calls"]) /
                                geomean(block["conjugate"]["target_s3_calls"]))
        summary["exploratory_baseline_over_conjugate_online_geomean"] = geomean(online_ratios)
        summary["target_bootstrap_95"] = bootstrap_95(online_ratios)
        summary["points_at_least_2x"] = sum(value >= 2 for value in online_ratios)
        summary["points_above_1x"] = sum(value > 1 for value in online_ratios)
        summary["baseline_over_conjugate_s3_call_geomean"] = geomean(calls_ratios)
        summary["arms"] = {}
        for arm in ("baseline", "conjugate"):
            selected = [row for row in rows if row["arm"] == arm]
            times = [row["report"]["timing_ms"] for row in selected]
            aa = [abs(math.log(block[arm]["online_ms"][0] / block[arm]["online_ms"][1]))
                  for block in target_blocks.values()]
            summary["arms"][arm] = {
                "candidate_id": plan["arms"][arm]["candidate_id"],
                "median_online_ms": statistics.median(t["target_online_after_reusable_setup"] for t in times),
                "median_charged_phase_ms": {key: statistics.median(t[key] for t in times) for key in PHASES},
                "median_target_s3_calls": statistics.median(row["report"]["target_s3_calls"] for row in selected),
                "median_reusable_setup_ms_supplementary": statistics.median(t["reusable_setup_total"] for t in times),
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
                      "online_ratio": result.get("exploratory_baseline_over_conjugate_online_geomean"),
                      "bootstrap_95": result.get("target_bootstrap_95")}, sort_keys=True))

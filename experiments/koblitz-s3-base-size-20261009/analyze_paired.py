"""Validate and summarize the frozen N53 S3 factor-base-size panel."""

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
BASELINE = ROOT / "experiments/koblitz-s3-parallel-target-20261009/FREEZE.json"
SOURCE = ROOT / "examples/s3_partial_prepared_parallel_targetspan.rs"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical(value: dict) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    return ordered[max(0, math.ceil(fraction * len(ordered)) - 1)]


def geomean(values: list[float]) -> float:
    assert values and all(value > 0 for value in values)
    return math.exp(statistics.mean(map(math.log, values)))


def bootstrap_95(values: list[float], seed: int) -> list[float]:
    rng = random.Random(seed)
    samples = sorted(geomean(rng.choices(values, k=len(values)))
                     for _ in range(20_000))
    return [samples[499], samples[19_499]]


def load_ledger(path: Path) -> tuple[bytes, list[dict]]:
    payload = path.read_bytes()
    raw = gzip.decompress(payload) if payload[:2] == b"\x1f\x8b" else payload
    rows = [json.loads(line) for line in raw.splitlines()]
    assert len(rows) == 48
    return raw, rows


def analyze(ledger: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    freeze_bytes = (HERE / "FREEZE.json").read_bytes()
    freeze = json.loads(freeze_bytes)
    baseline_bytes = BASELINE.read_bytes()
    parent_bytes = PARENT.read_bytes()
    assert freeze["kind"] == "s3_base_size_paired_v1"
    assert freeze["baseline_plan_sha256"] == sha(baseline_bytes)
    assert freeze["parent_plan_sha256"] == sha(parent_bytes)
    baseline = json.loads(baseline_bytes)
    parent = json.loads(parent_bytes)
    assert freeze["source_sha256"] == baseline["source_sha256"] == sha(SOURCE.read_bytes())
    assert freeze["binary_sha256"] == baseline["local_binary_sha256"]
    assert freeze["curve_id"] == parent["curves"]["53"]["curve_id"]
    assert freeze["columns"] == [244, 320, 396]
    assert freeze["target_labels"] == [f"T{i:03}" for i in range(25, 33)]
    assert freeze["repetitions_per_candidate_per_target"] == 2
    assert freeze["target_query_workers"] == freeze["relation_workers"] == 14
    assert freeze["controlled_speedup"] is None
    targets = {(item["n"], item["label"]): item for item in parent["targets"]}
    candidates = {}
    for columns in freeze["columns"]:
        candidate = freeze["candidate_ids"][str(columns)]
        if columns == 244:
            assert candidate == baseline["curves"]["53"]["parallel_candidate_id"]
            path = BASELINE.parent / "identities/candidates" / f"{candidate}.json"
        else:
            path = HERE / "identities/candidates" / f"{candidate}.json"
        record = json.loads(path.read_bytes())
        assert record.pop("candidate_id") == candidate
        assert sha(canonical(record))[:12] == candidate[-12:]
        base = record["factor_base"]
        assert base["effective_relation_columns"] == columns
        assert base["actual_usable_points_B"] == 106 * columns
        assert record["implementation"]["algorithm_flags"]["prepare_rank"] == columns - 1
        assert record["point_decomposition"]["target_query_workers"] == 14
        if columns != 244:
            dump = json.loads(gzip.decompress(
                (HERE / "local" / f"n53_base_c{columns}.json.gz").read_bytes()))
            assert dump["base_hash"] == base["enumerated_set_digest"]
            assert dump["factor_base_scanned_x"] == base["last_scanned_x"]
            assert len(dump["factor_base_point_coordinates"]) == base["actual_usable_points_B"]
        candidates[columns] = record
    raw, rows = load_ledger(ledger)
    blocks: dict[tuple[str, int], list[dict]] = defaultdict(list)
    failures = []
    for row in rows:
        columns = row["columns"]
        label = row["label"]
        assert row["kind"] == "s3_base_size_paired_single_target_run_v1"
        assert row["n"] == 53 and columns in candidates
        assert label in freeze["target_labels"]
        target = targets[(53, label)]
        candidate = freeze["candidate_ids"][str(columns)]
        assert row["candidate_id"] == candidate
        assert row["run_id"] == f"{candidate}W{target['workload_id']}R{row['repetition']}"
        assert row["workload_id"] == target["workload_id"]
        assert row["point_sha256"] == target["point_sha256"]
        assert row["source_sha256"] == freeze["source_sha256"]
        assert row["binary_sha256"] == freeze["binary_sha256"]
        assert row["argv"][-3:] == ["14", str(columns - 1), "14"]
        assert row["argv"][3] == str(columns)
        report = row["report"]
        if row["error"] or row["exit_code"] != 0 or not report:
            failures.append({"label": label, "columns": columns,
                             "repetition": row["repetition"], "error": row["error"],
                             "exit_code": row["exit_code"]})
        else:
            times = report["timing_ms"]
            assert all(isinstance(times[k], (int, float)) and times[k] >= 0 for k in PHASES)
            assert abs(sum(times[k] for k in PHASES) -
                       times["target_online_after_reusable_setup"]) < 0.02
            assert row["phase_ok"] and row["verified_native_and_fixture"]
            assert report["group_verified"] and report["target_count"] == 1
            assert report["target"] == target["public_point"]
            assert report["recovered_scalar"] == target["fixture_scalar"]
            assert report["factor_base_points"] == candidates[columns]["factor_base"]["actual_usable_points_B"]
            assert report["factor_base_digest"] == candidates[columns]["factor_base"]["enumerated_set_digest"]
            assert report["root_table_slots"] == candidates[columns]["point_decomposition"]["root_table_slots"]
            assert report["prepared_rank_reached"] and report["prepared_rank"] >= columns - 1
            assert report["target_query_workers"] == 14 and report["fallback_attempts"] == 0
        blocks[(label, columns)].append(row)
    assert len(blocks) == 24
    assert all(len(block) == 2 and {row["repetition"] for row in block} == {1, 2}
               for block in blocks.values())
    result = {"kind": "s3_base_size_paired_analysis_v1",
              "freeze_sha256": sha(freeze_bytes), "raw_ledger_sha256": sha(raw),
              "run_count": len(rows), "failure_rows": failures,
              "all_native_verified": not failures and all(row["verified_native_and_fixture"] for row in rows),
              "controlled_speedup": None, "factor_base_is_controlled_variable": True,
              "curves": {"53": {}}}
    if failures:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        return result
    by_target = {}
    for label in freeze["target_labels"]:
        arms = {}
        for columns in freeze["columns"]:
            runs = sorted(blocks[(label, columns)], key=lambda row: row["repetition"])
            arms[str(columns)] = {
                "online_ms": [row["report"]["timing_ms"]["target_online_after_reusable_setup"] for row in runs],
                "target_s3_calls": [row["report"]["target_s3_calls"] for row in runs],
                "reusable_setup_ms": [row["report"]["timing_ms"]["reusable_setup_total"] for row in runs],
                "prepared_rank": [row["report"]["prepared_rank"] for row in runs],
            }
        by_target[label] = arms
    result["target_blocks"] = by_target
    curve = result["curves"]["53"]
    for columns in (320, 396):
        online_ratios, call_ratios, setup_ratios = [], [], []
        for label, arms in by_target.items():
            old, new = arms["244"], arms[str(columns)]
            online_ratios.append(geomean(old["online_ms"]) / geomean(new["online_ms"]))
            call_ratios.append(geomean(old["target_s3_calls"]) / geomean(new["target_s3_calls"]))
            setup_ratios.append(geomean(old["reusable_setup_ms"]) / geomean(new["reusable_setup_ms"]))
        curve[f"244_over_{columns}"] = {
            "exploratory_online_geomean": geomean(online_ratios),
            "online_target_bootstrap_95": bootstrap_95(online_ratios, 20261009 + columns),
            "online_points_at_least_2x": sum(value >= 2 for value in online_ratios),
            "online_points_above_1x": sum(value > 1 for value in online_ratios),
            "target_s3_call_geomean_ratio": geomean(call_ratios),
            "s3_call_target_bootstrap_95": bootstrap_95(call_ratios, 20261109 + columns),
            "reusable_setup_geomean_ratio_supplementary": geomean(setup_ratios),
        }
    for columns in freeze["columns"]:
        selected = [row for row in rows if row["columns"] == columns]
        times = [row["report"]["timing_ms"] for row in selected]
        aa = [abs(math.log(arms[str(columns)]["online_ms"][0] /
                           arms[str(columns)]["online_ms"][1]))
              for arms in by_target.values()]
        curve[str(columns)] = {
            "candidate_id": freeze["candidate_ids"][str(columns)],
            "factor_base_points": candidates[columns]["factor_base"]["actual_usable_points_B"],
            "folded_columns": columns,
            "root_table_slots": candidates[columns]["point_decomposition"]["root_table_slots"],
            "run_count": len(selected),
            "median_online_ms": statistics.median(t["target_online_after_reusable_setup"] for t in times),
            "median_charged_phase_ms": {key: statistics.median(t[key] for t in times) for key in PHASES},
            "median_target_s3_calls": statistics.median(row["report"]["target_s3_calls"] for row in selected),
            "median_reusable_setup_ms_supplementary": statistics.median(t["reusable_setup_total"] for t in times),
            "aa_p95_abs_log_online_ratio": percentile(aa, 0.95),
            "actual_prepared_ranks": sorted({row["report"]["prepared_rank"] for row in selected}),
        }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, default=HERE / "local/raw_paired48.jsonl.gz")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.ledger, args.output)
    print(json.dumps({"runs": result["run_count"], "failures": len(result["failure_rows"]),
                      "curves": result["curves"]}, sort_keys=True))

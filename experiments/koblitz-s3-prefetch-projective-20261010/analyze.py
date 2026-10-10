"""Audit raw paired screen records and summarize complete one-target pairs."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import random
import statistics


HERE = Path(__file__).resolve().parent
PLAN = HERE / "SCREEN_FREEZE.json"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def geometric(values: list[float]) -> float:
    return math.exp(sum(math.log(value) for value in values) / len(values))


def bootstrap95(ratios: list[float]) -> list[float]:
    rng = random.Random(20261010)
    draws = sorted(geometric([ratios[rng.randrange(len(ratios))]
                              for _ in ratios]) for _ in range(20000))
    return [draws[499], draws[19499]]


def schedule(plan: dict) -> list[tuple[str, str]]:
    result = []
    for index, label in enumerate(plan["target_labels"]):
        order = ("baseline", "combined") if index % 2 == 0 else ("combined", "baseline")
        result.extend((label, arm) for arm in order)
    return result


def analyze(run_dir: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    frozen = PLAN.read_bytes()
    plan = json.loads(frozen)
    ledger_raw = (run_dir / "ledger.jsonl").read_bytes()
    ledger = [json.loads(line) for line in ledger_raw.splitlines()]
    completion_raw = (run_dir / "completion.json").read_bytes()
    completion = json.loads(completion_raw)
    assert completion["freeze_sha256"] == sha(frozen)
    assert completion["completed_run_count"] == len(ledger)
    assert completion["planned_run_count"] == len(schedule(plan)) == 16
    host = json.loads((run_dir / "host.json").read_bytes())
    assert host["host_exclusive_partition_receipt"] is None
    target_by_label = {row["label"]: row for row in plan["targets"]}
    verified = []
    failures = []
    for index, row in enumerate(ledger):
        label, arm = schedule(plan)[index]
        assert (row["index"], row["label"], row["arm"]) == (index, label, arm)
        receipt_path = run_dir / row["relative_receipt"]
        receipt_raw = receipt_path.read_bytes()
        assert sha(receipt_raw) == row["receipt_sha256"]
        receipt = json.loads(receipt_raw)
        assert receipt["freeze_sha256"] == sha(frozen)
        assert receipt["candidate_id"] == plan["arms"][arm]["candidate_id"]
        assert receipt["workload_id"] == target_by_label[label]["workload_id"]
        assert receipt["run_id"] == (receipt["candidate_id"] + "W" +
                                     receipt["workload_id"] + "R1")
        assert receipt["source_sha256"] == plan["arms"][arm]["source_sha256"]
        assert receipt["binary_sha256"] == plan["arms"][arm]["binary_sha256"]
        archive = receipt["archive"]
        unpacked = {}
        for name, meta in archive.items():
            packed = (receipt_path.parent / (name + ".gz")).read_bytes()
            assert sha(packed) == meta["gzip_sha256"]
            raw = gzip.decompress(packed)
            assert sha(raw) == meta["raw_sha256"] and len(raw) == meta["raw_bytes"]
            unpacked[name] = raw
        assert (unpacked["report.json"] == unpacked["stdout.bin"]) == (
            receipt["stdout_matches_report"])
        try:
            report = json.loads(unpacked["report.json"]) if unpacked["report.json"] else None
        except ValueError:
            report = None
        assert row["verified"] == receipt["verified_native_and_fixture"]
        assert row["memory_within_envelope"] == receipt["memory_within_envelope"]
        if not row["verified"] or not row["memory_within_envelope"]:
            failures.append({"label": label, "arm": arm,
                             "exit_code": receipt["exit_code"], "error": receipt["error"],
                             "receipt_sha256": row["receipt_sha256"]})
            continue
        assert receipt["phase_accounting_ok"] and receipt["exit_code"] == 0
        assert report["group_verified"] is True
        assert report["target_query_workers"] == 1
        assert report["target"] == target_by_label[label]["public_point"]
        assert report["recovered_scalar"] == target_by_label[label]["fixture_scalar"]
        timing = report["timing_ms"]
        online = timing["target_online_after_reusable_setup"]
        assert online > 0 and abs(sum(timing[key] for key in PHASES) - online) < .02
        assert receipt["online_ms"] == row["online_ms"] == online
        assert receipt["target_s3_calls"] == row["target_s3_calls"] == report["target_s3_calls"]
        assert receipt["reusable_setup_ms"] == timing["reusable_setup_total"]
        assert receipt["peak_rss_bytes"] == row["peak_rss_bytes"]
        assert receipt["peak_rss_bytes"] <= plan["memory_envelope_bytes"]
        verified.append({
            "label": label, "arm": arm, "candidate_id": receipt["candidate_id"],
            "workload_id": receipt["workload_id"], "run_id": receipt["run_id"],
            "online_ms": online, "target_phase_ms": {key: timing[key] for key in PHASES},
            "target_s3_calls": report["target_s3_calls"],
            "target_query_workers": report["target_query_workers"],
            "target_relation_probes": report["target_relation_probes"],
            "target_relation_indices": report["target_relation_indices"],
            "prepared_rank": report["prepared_rank"],
            "rank_attempts": report["rank_attempts"],
            "rank_verified_relations": report["rank_verified_relations"],
            "rank_new_rows": report["rank_new_rows"],
            "root_table_entries": report["root_table_entries"],
            "factor_base_points": report["factor_base_points"],
            "setup_ms": timing["reusable_setup_total"],
            "index_build_ms": timing["index_build"],
            "prepared_relations_ms": timing["prepared_relations"],
            "final_LA_ms": timing["rank_linear_algebra_final"],
            "total_s3_calls": report["total_s3_calls"],
            "peak_rss_bytes": receipt["peak_rss_bytes"],
            "involuntary_context_switches":
                receipt["resource_counters"]["involuntary_context_switches"],
            "receipt_sha256": row["receipt_sha256"],
        })
    by_pair = {}
    for row in verified:
        by_pair.setdefault(row["label"], {})[row["arm"]] = row
    pairs = []
    for label in plan["target_labels"]:
        item = by_pair.get(label, {})
        if set(item) != {"baseline", "combined"}:
            continue
        a, b = item["baseline"], item["combined"]
        assert a["factor_base_points"] == b["factor_base_points"] == 127200
        assert a["prepared_rank"] == b["prepared_rank"] == 1200
        assert a["target_relation_indices"] == b["target_relation_indices"]
        pairs.append({"label": label, "workload_id": a["workload_id"],
                      "baseline_online_ms": a["online_ms"],
                      "combined_online_ms": b["online_ms"],
                      "online_ratio": a["online_ms"] / b["online_ms"],
                      "baseline_target_s3_calls": a["target_s3_calls"],
                      "combined_target_s3_calls": b["target_s3_calls"],
                      "root_call_ratio": a["target_s3_calls"] / b["target_s3_calls"],
                      "baseline_setup_ms": a["setup_ms"],
                      "combined_setup_ms": b["setup_ms"],
                      "same_first_witness":
                          a["target_relation_indices"] == b["target_relation_indices"],
                      "baseline_peak_rss_bytes": a["peak_rss_bytes"],
                      "combined_peak_rss_bytes": b["peak_rss_bytes"]})
    all_complete = len(ledger) == 16 and len(verified) == 16 and len(pairs) == 8
    ratios = [pair["online_ratio"] for pair in pairs]
    call_ratios = [pair["root_call_ratio"] for pair in pairs]
    setup_control = {}
    for arm in ("baseline", "combined"):
        setup = [row["setup_ms"] for row in verified if row["arm"] == arm]
        median = statistics.median(setup) if setup else None
        setup_control[arm] = {
            "median_ms": median,
            "ratios_to_median": [value / median for value in setup] if median else [],
            "all_within_10_percent": (all(abs(value / median - 1) <= .10 for value in setup)
                                      if median else None),
        }
    helper_failures = [json.loads(path.read_bytes()) for path in
                       sorted(run_dir.glob("*_helper_failure.json"))]
    result = {
        "kind": "s3_n53_c1200_prefetch_projective_paired_screen_analysis_v1",
        "freeze_sha256": sha(frozen), "ledger_sha256": sha(ledger_raw),
        "completion_sha256": sha(completion_raw),
        "run_count": len(ledger), "verified_run_count": len(verified),
        "failed_runs": failures, "helper_failures": helper_failures,
        "unstarted_run_count": 16 - len(ledger) - len(helper_failures),
        "complete_pair_count": len(pairs), "all_complete": all_complete,
        "runs": verified, "pairs": pairs,
        "online_geomean_baseline_over_combined": geometric(ratios) if ratios else None,
        "online_bootstrap95_target_resample": bootstrap95(ratios) if ratios else None,
        "root_call_geomean_baseline_over_combined":
            geometric(call_ratios) if call_ratios else None,
        "root_call_bootstrap95_target_resample": bootstrap95(call_ratios) if call_ratios else None,
        "setup_negative_control": setup_control,
        "exploratory_twofold_stage_gate_pass":
            bool(all_complete and geometric(ratios) >= 2),
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: result[key] for key in
                      ("run_count", "verified_run_count", "complete_pair_count",
                       "online_geomean_baseline_over_combined",
                       "root_call_geomean_baseline_over_combined",
                       "exploratory_twofold_stage_gate_pass")},
                     sort_keys=True))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    analyze(args.run_dir, args.output)

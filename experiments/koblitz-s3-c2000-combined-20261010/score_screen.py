"""Audit C2000 screen archives and score deterministic S3-call savings."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
PLAN = HERE / "SCREEN_FREEZE.json"
REFERENCE = HERE.parent / "koblitz-s3-c1800-combined-20261010/analysis.json"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def load_archive(run_path: Path, receipt: dict, name: str) -> bytes:
    packed = (run_path / (name + ".gz")).read_bytes()
    expected = receipt["archive"][name]
    assert sha(packed) == expected["gzip_sha256"]
    raw = gzip.decompress(packed)
    assert len(raw) == expected["raw_bytes"]
    assert len(packed) == expected["gzip_bytes"]
    assert sha(raw) == expected["raw_sha256"]
    return raw


def score(run_dir: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    frozen = PLAN.read_bytes()
    plan = json.loads(frozen)
    assert plan["kind"] == "s3_n53_c2000_combined_screen_freeze_v1"
    reference_raw = REFERENCE.read_bytes()
    assert sha(reference_raw) == plan["reference_analysis_sha256"]
    reference = json.loads(reference_raw)
    reference_runs = {row["label"]: row for row in reference["runs"]
                      if row["arm"] == "combined" and row["repetition"] == 1}
    assert set(reference_runs) == set(plan["target_labels"])
    completion_raw = (run_dir / "completion.json").read_bytes()
    completion = json.loads(completion_raw)
    assert completion["freeze_sha256"] == sha(frozen)
    assert completion["stopped_reason"] is None
    assert completion["completed_run_count"] == len(plan["target_labels"]) == 8
    assert completion["verified_run_count"] == 8
    assert completion["memory_pass_count"] == 8
    ledger_raw = (run_dir / "ledger.jsonl").read_bytes()
    ledger = [json.loads(line) for line in ledger_raw.splitlines()]
    assert [row["label"] for row in ledger] == plan["run_order"]
    targets = {row["label"]: row for row in plan["targets"]}
    points = []
    for row in ledger:
        label = row["label"]
        assert row["verified"] and row["memory_within_envelope"]
        receipt_path = run_dir / row["relative_receipt"]
        receipt_raw = receipt_path.read_bytes()
        assert sha(receipt_raw) == row["receipt_sha256"]
        receipt = json.loads(receipt_raw)
        assert receipt["verified_native_and_fixture"]
        assert receipt["memory_within_envelope"]
        assert receipt["freeze_sha256"] == sha(frozen)
        assert receipt["source_sha256"] == plan["arms"]["combined"]["source_sha256"]
        assert receipt["binary_sha256"] == plan["arms"]["combined"]["binary_sha256"]
        assert receipt["candidate_id"] == row["candidate_id"] == plan["arms"]["combined"]["candidate_id"]
        assert receipt["workload_id"] == row["workload_id"] == targets[label]["workload_id"]
        assert receipt["run_id"] == row["run_id"]
        run_path = receipt_path.parent
        report_raw = load_archive(run_path, receipt, "report.json")
        assert load_archive(run_path, receipt, "stdout.bin") == report_raw
        assert load_archive(run_path, receipt, "stderr.txt") == b""
        report = json.loads(report_raw)
        assert report["target"] == targets[label]["public_point"]
        assert report["recovered_scalar"] == targets[label]["fixture_scalar"]
        assert report["target_s3_calls"] == receipt["target_s3_calls"] == row["target_s3_calls"]
        assert receipt["online_ms"] == row["online_ms"]
        assert receipt["peak_rss_bytes"] == row["peak_rss_bytes"]
        reference_calls = plan["reference_target_s3_calls"][label]
        assert isinstance(reference_calls, int) and reference_calls > 0
        assert isinstance(row["target_s3_calls"], int) and row["target_s3_calls"] > 0
        old = reference_runs[label]
        assert old["target_s3_calls"] == reference_calls
        assert old["workload_id"] == row["workload_id"]
        points.append({
            "label": label,
            "workload_id": row["workload_id"],
            "reference_c1800_target_s3_calls": reference_calls,
            "c2000_target_s3_calls": row["target_s3_calls"],
            "c1800_over_c2000_s3_call_ratio": reference_calls / row["target_s3_calls"],
            "same_first_relation_witness_as_c1800": (
                old["target_relation_indices"] == report["target_relation_indices"]),
            "same_state_probes_as_c1800": (
                old["target_state_probes"] == report["target_state_probes"]),
            "c2000_online_ms_exploratory": row["online_ms"],
            "c2000_setup_ms": receipt["reusable_setup_ms"],
            "c2000_peak_rss_bytes": row["peak_rss_bytes"],
            "receipt_sha256": sha(receipt_raw),
            "report_raw_sha256": sha(report_raw),
        })
    ratio = math.exp(sum(math.log(point["c1800_over_c2000_s3_call_ratio"])
                         for point in points) / len(points))
    threshold = float(plan["advance_root_call_geomean_ratio_min"])
    result = {
        "kind": "s3_n53_c2000_combined_screen_analysis_v1",
        "freeze_sha256": sha(frozen),
        "completion_sha256": sha(completion_raw),
        "ledger_sha256": sha(ledger_raw),
        "candidate_id": plan["arms"]["combined"]["candidate_id"],
        "reference_candidate_id": plan["reference_candidate_id"],
        "verified_run_count": len(points),
        "root_call_geomean_c1800_over_c2000": ratio,
        "advance_root_call_geomean_ratio_min": threshold,
        "advance_to_timed_confirmation": ratio >= threshold,
        "points": points,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: result[key] for key in
                      ("verified_run_count", "root_call_geomean_c1800_over_c2000",
                       "advance_to_timed_confirmation")}, sort_keys=True))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    score(args.run_dir.resolve(), args.output.resolve())

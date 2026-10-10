"""Audit the matched batch-8/batch-32 N53 screen and retain all runs."""

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
PLAN = HERE / "PAIR_FREEZE.json"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def archive(run_path: Path, receipt: dict, name: str) -> bytes:
    packed = (run_path / (name + ".gz")).read_bytes()
    expected = receipt["archive"][name]
    assert sha(packed) == expected["gzip_sha256"]
    raw = gzip.decompress(packed)
    assert sha(raw) == expected["raw_sha256"]
    assert len(raw) == expected["raw_bytes"]
    return raw


def quantile(sorted_values: list[float], proportion: float) -> float:
    position = (len(sorted_values) - 1) * proportion
    lower = int(position)
    upper = min(lower + 1, len(sorted_values) - 1)
    return sorted_values[lower] + (sorted_values[upper] - sorted_values[lower]) * (
        position - lower)


def analyze(run_dir: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    frozen = PLAN.read_bytes()
    plan = json.loads(frozen)
    assert plan["kind"] == "s3_n53_c1800_batch32_pair_freeze_v1"
    completion_raw = (run_dir / "completion.json").read_bytes()
    completion = json.loads(completion_raw)
    assert completion["freeze_sha256"] == sha(frozen)
    assert completion["stopped_reason"] is None
    assert completion["planned_run_count"] == completion["completed_run_count"] == 16
    assert completion["verified_run_count"] == completion["memory_pass_count"] == 16
    ledger_raw = (run_dir / "ledger.jsonl").read_bytes()
    ledger = [json.loads(line) for line in ledger_raw.splitlines()]
    assert len(ledger) == 16
    by_label: dict[str, dict[str, dict]] = {}
    run_rows = []
    target_by_label = {target["label"]: target for target in plan["targets"]}
    for index, entry in enumerate(ledger):
        label, arm = entry["label"], entry["arm"]
        assert entry["index"] == index and entry["repetition"] == 1
        assert arm in ("baseline", "candidate")
        assert entry["verified"] and entry["memory_within_envelope"]
        receipt_path = run_dir / entry["relative_receipt"]
        receipt_raw = receipt_path.read_bytes()
        assert sha(receipt_raw) == entry["receipt_sha256"]
        receipt = json.loads(receipt_raw)
        frozen_arm = plan["arms"][arm]
        assert receipt["freeze_sha256"] == sha(frozen)
        assert receipt["verified_native_and_fixture"]
        assert receipt["phase_accounting_ok"]
        assert receipt["memory_within_envelope"]
        assert receipt["candidate_id"] == entry["candidate_id"] == frozen_arm["candidate_id"]
        assert receipt["source_sha256"] == frozen_arm["source_sha256"]
        assert receipt["binary_sha256"] == frozen_arm["binary_sha256"]
        assert receipt["workload_id"] == entry["workload_id"] == target_by_label[label]["workload_id"]
        assert receipt["online_ms"] == entry["online_ms"]
        assert receipt["target_s3_calls"] == entry["target_s3_calls"]
        assert receipt["peak_rss_bytes"] == entry["peak_rss_bytes"]
        report_raw = archive(receipt_path.parent, receipt, "report.json")
        assert archive(receipt_path.parent, receipt, "stdout.bin") == report_raw
        assert archive(receipt_path.parent, receipt, "stderr.txt") == b""
        report = json.loads(report_raw)
        target = target_by_label[label]
        assert report["target"] == target["public_point"]
        assert report["recovered_scalar"] == target["fixture_scalar"]
        assert report["target_query_batch_width"] == frozen_arm["target_query_batch_width"]
        assert report["target_s3_calls"] == entry["target_s3_calls"]
        assert report["target_relation_indices"] == receipt["target_relation_indices"]
        assert arm not in by_label.setdefault(label, {})
        by_label[label][arm] = {"entry": entry, "receipt": receipt, "report": report,
                                "receipt_sha256": sha(receipt_raw),
                                "report_raw_sha256": sha(report_raw)}
        run_rows.append({
            "label": label, "arm": arm,
            "candidate_id": receipt["candidate_id"],
            "workload_id": receipt["workload_id"],
            "run_id": receipt["run_id"],
            "online_ms": entry["online_ms"],
            "target_phase_ms": receipt["target_phase_ms"],
            "target_s3_calls": entry["target_s3_calls"],
            "peak_rss_bytes": entry["peak_rss_bytes"],
            "reusable_setup_ms": receipt["reusable_setup_ms"],
            "receipt_sha256": sha(receipt_raw),
            "report_raw_sha256": sha(report_raw),
        })
    assert list(by_label) == plan["target_labels"]
    points = []
    for label in plan["target_labels"]:
        pair = by_label[label]
        assert set(pair) == {"baseline", "candidate"}
        baseline = pair["baseline"]
        candidate = pair["candidate"]
        bms = baseline["entry"]["online_ms"]
        cms = candidate["entry"]["online_ms"]
        assert bms > 0 and cms > 0
        points.append({
            "label": label,
            "workload_id": baseline["entry"]["workload_id"],
            "baseline_online_ms": bms,
            "candidate_online_ms": cms,
            "baseline_over_candidate_online_ratio": bms / cms,
            "baseline_target_s3_calls": baseline["entry"]["target_s3_calls"],
            "candidate_target_s3_calls": candidate["entry"]["target_s3_calls"],
            "same_first_relation_witness": (
                baseline["report"]["target_relation_indices"] ==
                candidate["report"]["target_relation_indices"]),
            "same_target_state_probes": (
                baseline["report"]["target_state_probes"] ==
                candidate["report"]["target_state_probes"]),
        })
    logs = [math.log(point["baseline_over_candidate_online_ratio"]) for point in points]
    ratio = math.exp(statistics.mean(logs))
    rng = random.Random(20261010)
    bootstrap = sorted(math.exp(sum(rng.choice(logs) for _ in logs) / len(logs))
                       for _ in range(20000))
    result = {
        "kind": "s3_n53_c1800_batch32_pair_analysis_v1",
        "freeze_sha256": sha(frozen),
        "completion_sha256": sha(completion_raw),
        "ledger_sha256": sha(ledger_raw),
        "verified_run_count": len(run_rows),
        "distinct_target_count": len(points),
        "online_geomean_baseline_over_candidate": ratio,
        "online_bootstrap95_target_resample": [quantile(bootstrap, .025),
                                                quantile(bootstrap, .975)],
        "advance_online_geomean_ratio_min": float(plan["advance_online_geomean_ratio_min"]),
        "advance_to_repeat_confirmation": ratio >= float(plan["advance_online_geomean_ratio_min"]),
        "same_first_relation_witness_count": sum(point["same_first_relation_witness"]
                                                  for point in points),
        "same_target_state_probes_count": sum(point["same_target_state_probes"]
                                              for point in points),
        "points": points,
        "runs": run_rows,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: result[key] for key in
                      ("verified_run_count", "online_geomean_baseline_over_candidate",
                       "online_bootstrap95_target_resample",
                       "advance_to_repeat_confirmation")}, sort_keys=True))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    analyze(args.run_dir.resolve(), args.output.resolve())

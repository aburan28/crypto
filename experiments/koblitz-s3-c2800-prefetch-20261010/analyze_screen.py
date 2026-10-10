"""Audit the complete frozen C2800 mixed/prefetch same-point screen."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import random


HERE = Path(__file__).resolve().parent
FREEZE = HERE / "SCREEN_FREEZE.json"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def geometric(values: list[float]) -> float:
    return math.exp(sum(math.log(x) for x in values) / len(values))


def bootstrap(values: list[float]) -> list[float]:
    rng = random.Random(20261010)
    draws = sorted(geometric([values[rng.randrange(len(values))]
                              for _ in values]) for _ in range(20000))
    return [draws[499], draws[19499]]


def analyze(run_dir: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    frozen = FREEZE.read_bytes()
    plan = json.loads(frozen)
    assert plan["kind"] == "s3_n53_c2800_prefetch_paired_screen_freeze_v1"
    host = json.loads((run_dir / "host.json").read_bytes())
    assert host["host_exclusive_partition_receipt"] is None
    rows = [json.loads(line) for line in (run_dir / "raw.jsonl").read_text().splitlines()]
    complete = len(rows) == 2 * len(plan["target_labels"])
    by_target: dict[str, dict] = {}
    failed = []
    for row in rows:
        assert row["kind"] == "s3_n53_c2800_prefetch_paired_screen_run_v1"
        assert row["freeze_sha256"] == sha(frozen)
        arm = row["arm"]
        assert row["candidate_id"] == plan["arms"][arm]["candidate_id"]
        assert row["binary_sha256"] == plan["arms"][arm]["binary_sha256"]
        assert row["source_sha256"] == plan["arms"][arm]["source_sha256"]
        assert row["run_id"] == (row["candidate_id"] + "W" +
                                  row["workload_id"] + "R1")
        report_path = run_dir / f"{row['label']}_{arm}_R1.json"
        raw = report_path.read_bytes() if report_path.exists() else b""
        assert row["raw_sha256"] == (sha(raw) if raw else None)
        report = json.loads(raw) if raw else None
        if not (row["verified"] and row["memory_within_envelope"] and report):
            failed.append({"label": row["label"], "arm": arm,
                           "error": row["error"], "exit_code": row["exit_code"]})
        if report:
            assert report["timing_ms"]["target_online_after_reusable_setup"] == row["target_online_ms"]
            assert report["target_s3_calls"] == row["target_s3_calls"]
            assert report["target_relation_indices"] == row["target_relation_indices"]
            row["report"] = report
        assert arm not in by_target.setdefault(row["label"], {})
        by_target[row["label"]][arm] = row
    assert all(label in plan["target_labels"] for label in by_target)
    assert all(set(pair) <= {"baseline", "prefetch"} for pair in by_target.values())
    targets = []
    ratios = []
    pdp_ratios = []
    setup_ratios = []
    same_witness = True
    semantic_differences = {}
    allowed = {"timing_ms", "peak_rss_bytes", "target_state_probes",
               "target_s3_calls", "total_s3_calls"}
    for label in plan["target_labels"]:
        pair = by_target.get(label, {})
        if set(pair) != {"baseline", "prefetch"} or any(
            arm.get("report") is None for arm in pair.values()
        ):
            targets.append({"label": label, "status": "incomplete_or_failed"})
            continue
        baseline, prefetch = pair["baseline"], pair["prefetch"]
        a, b = baseline["report"], prefetch["report"]
        differences = sorted(key for key in set(a) | set(b)
                             if a.get(key) != b.get(key))
        unexpected = [key for key in differences if key not in allowed]
        semantic_differences[label] = unexpected
        same_witness &= (a.get("target_relation_indices") == b.get("target_relation_indices")
                         and a.get("target_relation_probes") == b.get("target_relation_probes"))
        ratio = baseline["target_online_ms"] / prefetch["target_online_ms"]
        pdp_ratio = baseline["target_pdp_charged_ms"] / prefetch["target_pdp_charged_ms"]
        setup_ratio = (a["timing_ms"]["reusable_setup_total"] /
                       b["timing_ms"]["reusable_setup_total"])
        ratios.append(ratio)
        pdp_ratios.append(pdp_ratio)
        setup_ratios.append(setup_ratio)
        targets.append({
            "label": label, "workload_id": baseline["workload_id"],
            "baseline_online_ms": baseline["target_online_ms"],
            "prefetch_online_ms": prefetch["target_online_ms"],
            "baseline_pdp_charged_ms": baseline["target_pdp_charged_ms"],
            "prefetch_pdp_charged_ms": prefetch["target_pdp_charged_ms"],
            "baseline_online_phases_ms": {key: a["timing_ms"][key] for key in PHASES},
            "prefetch_online_phases_ms": {key: b["timing_ms"][key] for key in PHASES},
            "baseline_s3_calls": baseline["target_s3_calls"],
            "prefetch_s3_calls": prefetch["target_s3_calls"],
            "baseline_pdp_us_per_s3_call":
                1000 * baseline["target_pdp_charged_ms"] / baseline["target_s3_calls"],
            "prefetch_pdp_us_per_s3_call":
                1000 * prefetch["target_pdp_charged_ms"] / prefetch["target_s3_calls"],
            "baseline_reusable_setup_ms": a["timing_ms"]["reusable_setup_total"],
            "prefetch_reusable_setup_ms": b["timing_ms"]["reusable_setup_total"],
            "baseline_peak_rss_bytes": baseline["peak_rss_bytes"],
            "prefetch_peak_rss_bytes": prefetch["peak_rss_bytes"],
            "online_ratio": ratio, "pdp_ratio": pdp_ratio,
            "reusable_setup_negative_control_ratio": setup_ratio,
            "reusable_setup_negative_control_within_10_percent":
                .90 <= setup_ratio <= 1.10,
            "first_witness_equal": (a.get("target_relation_indices") ==
                                    b.get("target_relation_indices") and
                                    a.get("target_relation_probes") ==
                                    b.get("target_relation_probes")),
            "unexpected_report_differences": unexpected,
        })
    valid = (complete and not failed and len(ratios) == len(plan["target_labels"])
             and same_witness and all(not x for x in semantic_differences.values()))
    result = {
        "kind": "s3_n53_c2800_prefetch_paired_screen_analysis_v1",
        "freeze_sha256": sha(frozen),
        "ledger_sha256": sha((run_dir / "raw.jsonl").read_bytes()),
        "run_count": len(rows), "target_count": len(targets),
        "complete": complete,
        "verified_and_memory_ok_count": len(rows) - len(failed),
        "failed_runs": failed,
        "same_first_witness_all": same_witness if valid else None,
        "unexpected_report_differences": semantic_differences,
        "online_geomean_baseline_over_prefetch": geometric(ratios) if valid else None,
        "online_bootstrap95_target_resample": bootstrap(ratios) if valid else None,
        "pdp_geomean_baseline_over_prefetch": geometric(pdp_ratios) if valid else None,
        "pdp_bootstrap95_target_resample": bootstrap(pdp_ratios) if valid else None,
        "setup_negative_control_geomean_ratio":
            geometric(setup_ratios) if valid else None,
        "setup_negative_control_all_within_10_percent":
            all(.90 <= ratio <= 1.10 for ratio in setup_ratios) if valid else None,
        "exploratory_threshold_met": (geometric(ratios) >=
                                      float(plan["exploratory_geomean_speedup_threshold"])
                                      if valid else None),
        "peak_rss_bytes_max": (max(row["peak_rss_bytes"] for row in rows)
                               if rows else None),
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
        "targets": targets,
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: result[key] for key in
                      ("run_count", "verified_and_memory_ok_count",
                       "same_first_witness_all",
                       "online_geomean_baseline_over_prefetch",
                       "online_bootstrap95_target_resample",
                       "exploratory_threshold_met")}, sort_keys=True))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.run_dir, args.output)

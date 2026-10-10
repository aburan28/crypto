"""Audit same-index lookup modes with balanced triplet and CPU-time summaries."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import random
import statistics


HERE = Path(__file__).resolve().parent
PLAN = HERE / "DIAGNOSTIC_FREEZE.json"
SCHEDULE = HERE / "SCHEDULE.json"
NAMES = ("none", "late", "early")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def geometric(values: list[float]) -> float:
    return math.exp(sum(math.log(value) for value in values) / len(values))


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, max(0, round(fraction * (len(ordered) - 1))))]


def bootstrap_ratios(ratios: list[float]) -> list[float]:
    rng = random.Random(20261010)
    draws = sorted(geometric([ratios[rng.randrange(len(ratios))]
                              for _ in ratios]) for _ in range(20000))
    return [draws[499], draws[19499]]


def analyze(run_dir: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    frozen = PLAN.read_bytes()
    plan = json.loads(frozen)
    receipt = json.loads((run_dir / "receipt.json").read_bytes())
    assert receipt["kind"] == "s3_c2800_lookup_ab_t025_diagnostic_receipt_v1"
    assert receipt["freeze_sha256"] == sha(frozen)
    assert receipt["source_sha256"] == plan["source_sha256"]
    assert receipt["binary_sha256"] == plan["binary_sha256"]
    assert receipt["candidate_ids"] == plan["candidate_ids"]
    primary_raw = (run_dir / "report.json").read_bytes()
    diagnostic_raw = (run_dir / "diagnostic.jsonl").read_bytes()
    assert receipt["primary_report_sha256"] == sha(primary_raw)
    assert receipt["diagnostic_jsonl_sha256"] == sha(diagnostic_raw)
    primary = json.loads(primary_raw)
    rows = [json.loads(line) for line in diagnostic_raw.splitlines()]
    schedule = json.loads(SCHEDULE.read_bytes())
    assert sha(SCHEDULE.read_bytes()) == plan["schedule_sha256"]
    assert len(rows) == len(schedule) == 54
    assert receipt["primary_verified"] and receipt["diagnostic_verified"]
    assert receipt["memory_within_envelope"]
    assert receipt["unexpected_reference_differences"] == []
    assert primary["target"] == plan["public_point"]
    assert primary["recovered_scalar"] == plan["fixture_scalar"]
    witness = primary["target_relation_indices"]
    probes = primary["target_relation_probes"]
    calls = primary["target_s3_calls"]
    blocks = []
    for base in range(0, len(rows), 3):
        triplet = rows[base:base + 3]
        assert sorted(row["mode"] for row in triplet) == [0, 1, 2]
        block = {}
        for index, row in enumerate(triplet, base):
            mode = row["mode"]
            assert row["sequence_index"] == index
            assert mode == schedule[index]
            assert row["mode_name"] == NAMES[mode]
            assert row["cache_flush_bytes"] == plan["cache_flush_bytes_per_query"]
            assert row["target"] == plan["public_point"]
            assert row["point_indices"] == witness
            assert row["partner_probes"] == probes
            assert row["s3_calls"] == calls
            assert row["recovered_scalar"] == plan["fixture_scalar"]
            assert row["group_verified"] is True
            assert row["wall_ns"] > 0
            pdp_ns = max(1.0, row["wall_ns"] - row["relation_check_ms"] * 1_000_000)
            block[mode] = {
                "wall_ns": row["wall_ns"],
                "thread_cpu_ns": row["thread_cpu_ns"],
                "wall_pdp_ns": pdp_ns,
            }
        blocks.append(block)
    cpu_complete = all(block[mode]["thread_cpu_ns"] is not None
                       and block[mode]["thread_cpu_ns"] > 0
                       for block in blocks for mode in range(3))
    modes = {}
    for mode, name in enumerate(NAMES):
        wall = [block[mode]["wall_ns"] / 1_000_000 for block in blocks]
        pdp = [block[mode]["wall_pdp_ns"] / 1_000_000 for block in blocks]
        cpu = ([block[mode]["thread_cpu_ns"] / 1_000_000 for block in blocks]
               if cpu_complete else None)
        gaps = ([block[mode]["wall_ns"] / block[mode]["thread_cpu_ns"]
                 for block in blocks] if cpu_complete else None)
        modes[name] = {
            "candidate_id": plan["candidate_ids"][name],
            "queries": len(blocks),
            "wall_ms_median": statistics.median(wall),
            "wall_ms_geomean": geometric(wall),
            "wall_pdp_ms_median": statistics.median(pdp),
            "thread_cpu_ms_median": statistics.median(cpu) if cpu else None,
            "thread_cpu_ms_geomean": geometric(cpu) if cpu else None,
            "wall_over_thread_cpu_ratio_p95": percentile(gaps, .95) if gaps else None,
            "wall_over_thread_cpu_ratio_gt_1_10":
                sum(ratio > 1.10 for ratio in gaps) if gaps else None,
        }
    comparisons = {}
    for numerator, denominator in ((0, 1), (0, 2), (1, 2)):
        label = NAMES[numerator] + "_over_" + NAMES[denominator]
        metrics = {}
        for key in ("wall_ns", "wall_pdp_ns", "thread_cpu_ns"):
            if key == "thread_cpu_ns" and not cpu_complete:
                metrics[key] = None
                continue
            ratios = [block[numerator][key] / block[denominator][key]
                      for block in blocks]
            metrics[key] = {
                "paired_triplet_geomean_ratio": geometric(ratios),
                "bootstrap95_triplet_resample": bootstrap_ratios(ratios),
                "median_triplet_ratio": statistics.median(ratios),
            }
        comparisons[label] = metrics
    result = {
        "kind": "s3_c2800_lookup_ab_t025_diagnostic_analysis_v1",
        "freeze_sha256": sha(frozen),
        "receipt_sha256": sha((run_dir / "receipt.json").read_bytes()),
        "primary_report_sha256": sha(primary_raw),
        "diagnostic_jsonl_sha256": sha(diagnostic_raw),
        "primary_run_id": receipt["primary_run_id"],
        "primary_verified": True,
        "primary_online_ms": primary["timing_ms"]["target_online_after_reusable_setup"],
        "primary_candidate_id": plan["candidate_ids"]["late"],
        "workload_id": plan["workload_id"],
        "query_count": len(rows), "balanced_triplet_count": len(blocks),
        "all_query_scalars_and_witnesses_verified": True,
        "thread_cpu_complete": cpu_complete,
        "cache_flush_bytes_per_query": plan["cache_flush_bytes_per_query"],
        "peak_rss_bytes": receipt["peak_rss_bytes"],
        "modes": modes,
        "comparisons": comparisons,
        "shared_preparation": True,
        "repeated_known_target_stage_diagnostic": True,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"primary_verified": True,
                      "query_count": len(rows),
                      "thread_cpu_complete": cpu_complete,
                      "comparisons": {name: metrics["wall_ns"]["paired_triplet_geomean_ratio"]
                                      for name, metrics in comparisons.items()}},
                     sort_keys=True))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.run_dir, args.output)

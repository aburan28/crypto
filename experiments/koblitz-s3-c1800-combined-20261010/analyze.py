"""Audit every confirmation run and summarize target-level ABBA ratios."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import random
import statistics

from run_paired import planned_rows


HERE = Path(__file__).resolve().parent
PLAN = HERE / "CONFIRM_FREEZE.json"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def geometric(values: list[float]) -> float:
    assert values and all(value > 0 for value in values)
    return math.exp(statistics.mean(map(math.log, values)))


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    low = int(position)
    return ordered[low] + (ordered[min(low + 1, len(ordered) - 1)] -
                           ordered[low]) * (position - low)


def bootstrap95(ratios: list[float]) -> list[float]:
    rng = random.Random(20261010)
    draws = [geometric([ratios[rng.randrange(len(ratios))]
                        for _ in ratios]) for _ in range(20000)]
    return [percentile(draws, .025), percentile(draws, .975)]


def analyze(run_dir: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    frozen = PLAN.read_bytes()
    plan = json.loads(frozen)
    schedule = planned_rows(plan)
    assert len(schedule) == 32
    ledger_raw = (run_dir / "ledger.jsonl").read_bytes()
    ledger = [json.loads(line) for line in ledger_raw.splitlines()]
    completion_raw = (run_dir / "completion.json").read_bytes()
    completion = json.loads(completion_raw)
    assert completion["freeze_sha256"] == sha(frozen)
    assert completion["completed_run_count"] == len(ledger)
    assert completion["planned_run_count"] == 32
    host = json.loads((run_dir / "host.json").read_bytes())
    assert host["host_exclusive_partition_receipt"] is None
    target_by_label = {row["label"]: row for row in plan["targets"]}
    verified = []
    failures = []
    for index, row in enumerate(ledger):
        label, arm, repeat = schedule[index]
        assert (row["index"], row["label"], row["arm"], row["repetition"]) == (
            index, label, arm, repeat)
        receipt_path = run_dir / row["relative_receipt"]
        receipt_raw = receipt_path.read_bytes()
        assert sha(receipt_raw) == row["receipt_sha256"]
        receipt = json.loads(receipt_raw)
        assert receipt["freeze_sha256"] == sha(frozen)
        assert (receipt["label"], receipt["arm"], receipt["repetition"]) == (
            label, arm, repeat)
        assert receipt["candidate_id"] == plan["arms"][arm]["candidate_id"]
        assert receipt["workload_id"] == target_by_label[label]["workload_id"]
        assert receipt["run_id"] == (receipt["candidate_id"] + "W" +
                                     receipt["workload_id"] + "R" + str(repeat))
        assert receipt["source_sha256"] == plan["arms"][arm]["source_sha256"]
        assert receipt["binary_sha256"] == plan["arms"][arm]["binary_sha256"]
        unpacked = {}
        for name, meta in receipt["archive"].items():
            packed = (receipt_path.parent / (name + ".gz")).read_bytes()
            assert sha(packed) == meta["gzip_sha256"]
            raw = gzip.decompress(packed)
            assert sha(raw) == meta["raw_sha256"]
            assert len(raw) == meta["raw_bytes"] and len(packed) == meta["gzip_bytes"]
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
            failures.append({"label": label, "arm": arm, "repetition": repeat,
                             "exit_code": receipt["exit_code"], "error": receipt["error"],
                             "receipt_sha256": row["receipt_sha256"]})
            continue
        target = target_by_label[label]
        assert receipt["phase_accounting_ok"] and receipt["exit_code"] == 0
        assert report["group_verified"] is True
        assert report["target"] == target["public_point"]
        assert report["recovered_scalar"] == target["fixture_scalar"]
        assert report["target_count"] == report["target_query_workers"] == 1
        frozen_arm = plan["arms"][arm]
        assert report["factor_base_digest"] == frozen_arm["factor_base_digest"]
        assert report["factor_base_points"] == frozen_arm["factor_base_points"]
        assert report["prepared_rank"] == frozen_arm["prepare_rank"]
        assert report["fallback_attempts"] == 0
        assert report["target_full_rank_direct_recovery"] is True
        timing = report["timing_ms"]
        online = timing["target_online_after_reusable_setup"]
        assert online > 0 and abs(sum(timing[key] for key in PHASES) - online) < .02
        assert receipt["online_ms"] == row["online_ms"] == online
        assert receipt["target_s3_calls"] == row["target_s3_calls"] == report["target_s3_calls"]
        assert receipt["reusable_setup_ms"] == timing["reusable_setup_total"]
        assert receipt["peak_rss_bytes"] == row["peak_rss_bytes"]
        assert receipt["peak_rss_bytes"] <= plan["memory_envelope_bytes"]
        verified.append({
            "label": label, "arm": arm, "repetition": repeat,
            "candidate_id": receipt["candidate_id"],
            "workload_id": receipt["workload_id"], "run_id": receipt["run_id"],
            "online_ms": online, "target_phase_ms": {key: timing[key] for key in PHASES},
            "target_s3_calls": report["target_s3_calls"],
            "target_state_probes": report["target_state_probes"],
            "target_relation_indices": report["target_relation_indices"],
            "setup_ms": timing["reusable_setup_total"],
            "rank_attempts": report["rank_attempts"],
            "rank_verified_relations": report["rank_verified_relations"],
            "rank_new_rows": report["rank_new_rows"],
            "final_LA_ms": timing["rank_linear_algebra_final"],
            "peak_rss_bytes": receipt["peak_rss_bytes"],
            "involuntary_context_switches":
                receipt["resource_counters"]["involuntary_context_switches"],
            "report_raw_sha256": sha(unpacked["report.json"]),
            "receipt_sha256": row["receipt_sha256"],
        })
    by_key = {(row["label"], row["arm"], row["repetition"]): row for row in verified}
    points = []
    aa = {"baseline": [], "combined": []}
    for label in plan["target_labels"]:
        keys = [(label, arm, repeat) for arm in ("baseline", "combined")
                for repeat in (1, 2)]
        if any(key not in by_key for key in keys):
            continue
        runs = [by_key[key] for key in keys]
        witnesses = {
            arm: {tuple(by_key[label, arm, repeat]["target_relation_indices"])
                  for repeat in (1, 2)}
            for arm in ("baseline", "combined")
        }
        assert all(len(values) == 1 for values in witnesses.values())
        samples = {arm: [by_key[label, arm, repeat]["online_ms"] for repeat in (1, 2)]
                   for arm in ("baseline", "combined")}
        for arm in samples:
            aa[arm].append(abs(math.log(samples[arm][0] / samples[arm][1])))
        base = geometric(samples["baseline"])
        combined = geometric(samples["combined"])
        points.append({
            "label": label, "workload_id": runs[0]["workload_id"],
            "same_first_relation_witness":
                next(iter(witnesses["baseline"])) == next(iter(witnesses["combined"])),
            "baseline_online_ms": samples["baseline"],
            "combined_online_ms": samples["combined"],
            "baseline_online_geomean_ms": base,
            "combined_online_geomean_ms": combined,
            "online_ratio": base / combined,
            "baseline_target_s3_calls": by_key[label, "baseline", 1]["target_s3_calls"],
            "combined_target_s3_calls": by_key[label, "combined", 1]["target_s3_calls"],
            "root_call_ratio": (by_key[label, "baseline", 1]["target_s3_calls"] /
                                by_key[label, "combined", 1]["target_s3_calls"]),
        })
    all_complete = len(ledger) == len(verified) == 32 and len(points) == 8
    ratios = [point["online_ratio"] for point in points]
    calls = [point["root_call_ratio"] for point in points]
    online_interval = bootstrap95(ratios) if all_complete else None
    call_interval = bootstrap95(calls) if all_complete else None
    online_geo = geometric(ratios) if all_complete else None
    setup_control = {}
    for arm in ("baseline", "combined"):
        setups = [row["setup_ms"] for row in verified if row["arm"] == arm]
        median = statistics.median(setups) if setups else None
        setup_control[arm] = {
            "median_ms": median,
            "ratios_to_median": [value / median for value in setups] if median else [],
            "all_within_10_percent": (all(abs(value / median - 1) <= .10
                                          for value in setups) if median else None),
        }
    aa_p95 = {arm: percentile(values, .95) if values else None
              for arm, values in aa.items()}
    replay_path = run_dir / "independent_sage_replay.json"
    if replay_path.exists():
        replay = json.loads(replay_path.read_bytes())
        assert replay["freeze_sha256"] == sha(frozen)
        assert replay["verified_run_count"] == 32
        assert replay["generator_order_verified"]
        assert {(row["label"], row["arm"], row["report_raw_sha256"])
                for row in replay["runs"]} == {
                    (row["label"], row["arm"], row["report_raw_sha256"])
                    for row in verified}
        sage_replay_verified = True
    else:
        sage_replay_verified = False
    helper_failures = [json.loads(path.read_bytes()) for path in
                       sorted(run_dir.glob("*_helper_failure.json"))]
    result = {
        "kind": "s3_n53_c1800_combined_confirmation_analysis_v1",
        "freeze_sha256": sha(frozen), "ledger_sha256": sha(ledger_raw),
        "completion_sha256": sha(completion_raw),
        "run_count": len(ledger), "verified_run_count": len(verified),
        "failed_runs": failures, "helper_failures": helper_failures,
        "unstarted_run_count": 32 - len(ledger) - len(helper_failures),
        "complete_point_count": len(points), "all_complete": all_complete,
        "runs": verified, "points": points,
        "online_geomean_baseline_over_combined": online_geo,
        "online_bootstrap95_target_resample": online_interval,
        "root_call_geomean_baseline_over_combined": geometric(calls) if all_complete else None,
        "root_call_bootstrap95_target_resample": call_interval,
        "aa_p95_absolute_log_ratio": aa_p95,
        "setup_negative_control": setup_control,
        "independent_sage_replay_verified": sage_replay_verified,
        "local_confirmation_gate_pass": bool(
            all_complete and sage_replay_verified and online_geo >= 2
            and online_interval[0] > 1),
        "twofold_threshold_95_lower_pass": bool(
            all_complete and sage_replay_verified and online_interval[0] > 2),
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: result[key] for key in
                      ("run_count", "verified_run_count", "complete_point_count",
                       "online_geomean_baseline_over_combined",
                       "online_bootstrap95_target_resample",
                       "local_confirmation_gate_pass",
                       "twofold_threshold_95_lower_pass")}, sort_keys=True))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    analyze(args.run_dir, args.output)

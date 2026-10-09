"""Validate the dense-root paired screen and preserve its A/A noise."""

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
ROOT = HERE.parents[2]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
ARMS = ("original", "dense")
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")
CAPACITY_FIELDS = {"root_table_slots", "root_table_slots_per_shard",
                   "root_table_slot_bytes", "root_table_value_bits",
                   "root_table_capacity_policy"}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def geomean(values: list[float]) -> float:
    assert values and all(value > 0 for value in values)
    return math.exp(statistics.mean(map(math.log, values)))


def p95(values: list[float]) -> float:
    ordered = sorted(values)
    return ordered[math.ceil(.95 * len(ordered)) - 1]


def bootstrap_95(values: list[float]) -> list[float]:
    rng = random.Random(20261009)
    samples = sorted(geomean(rng.choices(values, k=len(values)))
                     for _ in range(20_000))
    return [samples[499], samples[19_499]]


def analyze(ledger: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    plan_bytes = (HERE / "FREEZE.json").read_bytes()
    plan = json.loads(plan_bytes)
    parent_bytes = PARENT.read_bytes()
    assert plan["kind"] == "s3_dense_root_oneworker_paired_screen_v1"
    assert plan["parent_freeze_sha256"] == sha(parent_bytes)
    assert plan["runner_sha256"] == sha((HERE / "run_paired.py").read_bytes())
    targets = {(row["n"], row["label"]): row
               for row in json.loads(parent_bytes)["targets"]}
    payload = ledger.read_bytes()
    raw = gzip.decompress(payload) if payload[:2] == b"\x1f\x8b" else payload
    rows = [json.loads(line) for line in raw.splitlines()]
    assert len(rows) == 32
    blocks = defaultdict(list)
    failures = []
    expected_order = []
    for label_index, label in enumerate(plan["target_labels"]):
        for n in (41, 53):
            order = (["original", "dense", "dense", "original"] if
                     (label_index + (n == 53)) % 2 == 0 else
                     ["dense", "original", "original", "dense"])
            expected_order.extend((n, label, arm) for arm in order)
    assert [(row["n"], row["label"], row["arm"]) for row in rows] == expected_order
    for row in rows:
        assert row["kind"] == "s3_dense_root_oneworker_paired_screen_run_v1"
        n, arm, label = row["n"], row["arm"], row["label"]
        assert n in (41, 53) and arm in ARMS and label in plan["target_labels"]
        frozen = plan["curves"][str(n)][arm]
        target = targets[(n, label)]
        assert row["candidate_id"] == frozen["candidate_id"]
        assert row["workload_id"] == target["workload_id"]
        assert row["run_id"] == f"{frozen['candidate_id']}W{target['workload_id']}R{row['repetition']}"
        assert row["point_sha256"] == target["point_sha256"]
        assert row["source_sha256"] == frozen["source_sha256"]
        assert row["binary_sha256"] == frozen["binary_sha256"]
        assert row["argv"][3] == str(frozen["columns"])
        assert row["argv"][-3:] == ["14", str(frozen["prepare_rank"]), "1"]
        report = row["report"]
        if row["error"] or row["exit_code"] != 0 or report is None or not row["verified_native_and_fixture"]:
            failures.append({"n": n, "label": label, "arm": arm,
                             "repetition": row["repetition"], "error": row["error"],
                             "exit_code": row["exit_code"]})
        else:
            times = report["timing_ms"]
            assert all(isinstance(times[k], (int, float)) and times[k] >= 0 for k in PHASES)
            assert abs(sum(times[k] for k in PHASES) -
                       times["target_online_after_reusable_setup"]) < .02
            assert row["phase_ok"] and report["group_verified"]
            assert report["target"] == target["public_point"]
            assert report["recovered_scalar"] == target["fixture_scalar"]
            assert report["prepared_rank"] == frozen["prepare_rank"]
            assert report["factor_base_digest"] == frozen["factor_base_digest"]
            assert report["root_table_slots"] == frozen["root_table_slots"]
            assert report["fallback_attempts"] == 0
        blocks[(n, label, arm)].append(row)
    assert len(blocks) == 16
    assert all(len(group) == 2 and {row["repetition"] for row in group} == {1, 2}
               for group in blocks.values())
    host = json.loads((ledger.parent / "host.json").read_bytes())
    assert host["host_exclusive_partition_receipt"] is None
    result = {"kind": "s3_dense_root_oneworker_paired_screen_analysis_v1",
              "freeze_sha256": sha(plan_bytes), "raw_ledger_sha256": sha(raw),
              "run_count": len(rows), "failure_rows": failures,
              "all_native_verified": not failures,
              "host": host, "controlled_speedup": None,
              "curves": {}}
    if not failures:
        for n in (41, 53):
            ratios = []
            arm_times = {arm: [] for arm in ARMS}
            aa = {arm: [] for arm in ARMS}
            semantic_equal = True
            for label in plan["target_labels"]:
                groups = {arm: sorted(blocks[(n, label, arm)],
                                      key=lambda row: row["repetition"]) for arm in ARMS}
                for repetition in (0, 1):
                    original = groups["original"][repetition]["report"]
                    dense = groups["dense"][repetition]["report"]
                    semantic_equal &= ({k: v for k, v in original.items()
                                        if k != "timing_ms" and k not in CAPACITY_FIELDS} ==
                                       {k: v for k, v in dense.items()
                                        if k != "timing_ms" and k not in CAPACITY_FIELDS})
                point_times = {}
                for arm in ARMS:
                    pair = [row["report"]["timing_ms"]["target_online_after_reusable_setup"]
                            for row in groups[arm]]
                    point_times[arm] = pair
                    arm_times[arm].extend(row["report"]["timing_ms"] for row in groups[arm])
                    aa[arm].append(abs(math.log(pair[0] / pair[1])))
                ratios.append(geomean(point_times["original"]) /
                              geomean(point_times["dense"]))
            result["curves"][str(n)] = {
                "target_count": len(ratios),
                "all_non_timing_semantics_equal": semantic_equal,
                "exploratory_original_over_dense_online_geomean": geomean(ratios),
                "target_bootstrap_95": bootstrap_95(ratios),
                "point_ratios": dict(zip(plan["target_labels"], ratios)),
                "points_dense_faster": sum(value > 1 for value in ratios),
                "aa_p95_abs_log_online_ratio": {arm: p95(values) for arm, values in aa.items()},
                "aa_noise_gate_pass": all(p95(values) <= math.log(1.10)
                                          for values in aa.values()),
                "arms": {arm: {
                    "candidate_id": plan["curves"][str(n)][arm]["candidate_id"],
                    "root_table_slots": plan["curves"][str(n)][arm]["root_table_slots"],
                    "root_payload_bytes": (16 if arm == "original" else 13) *
                                          plan["curves"][str(n)][arm]["root_table_slots"],
                    "median_online_ms": statistics.median(
                        t["target_online_after_reusable_setup"] for t in arm_times[arm]),
                    "median_charged_phase_ms": {key: statistics.median(t[key]
                                                for t in arm_times[arm]) for key in PHASES},
                    "median_reusable_setup_ms_supplementary": statistics.median(
                        t["reusable_setup_total"] for t in arm_times[arm]),
                    "median_index_build_ms_supplementary": statistics.median(
                        t["index_build"] for t in arm_times[arm]),
                } for arm in ARMS},
            }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({str(n): result["curves"].get(str(n), {}).get(
        "exploratory_original_over_dense_online_geomean") for n in (41, 53)},
        sort_keys=True))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.ledger, args.output)

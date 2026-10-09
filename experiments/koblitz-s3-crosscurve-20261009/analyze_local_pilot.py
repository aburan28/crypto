"""Summarize frozen local S3 ABBA runs and their exact online root-call counts.

Local wall ratios are exploratory. A separate host isolation receipt is
required before interpreting them as controlled CPU speedups.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import statistics


HERE = Path(__file__).resolve().parent
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")
Z95_ONE_SIDED = 1.6448536269514722
Z95_TWO_SIDED = 1.959963984540054
Z80 = 0.8416212335729143


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def root_calls(row: dict) -> int:
    report = row["report"]
    if row["arm"] == "reference":
        online_rank = report["rank_target_s3_calls"]
    else:
        online_rank = report["fallback_s3_calls"]
    calls = online_rank + report["target_s3_calls"]
    assert isinstance(calls, int) and calls >= 0
    return calls


def analyze(ledger: Path, output: Path) -> dict:
    plan_bytes = (HERE / "FREEZE.json").read_bytes()
    plan = json.loads(plan_bytes)
    raw_bytes = ledger.read_bytes()
    rows = [json.loads(line) for line in raw_bytes.splitlines()]
    assert rows and len(rows) % 8 == 0
    targets = {(target["n"], target["label"]): target for target in plan["targets"]}
    grouped = defaultdict(list)
    binary_hashes = defaultdict(set)
    run_ids = set()
    for row in rows:
        n, label, arm = row["n"], row["label"], row["arm"]
        assert n in (41, 53) and arm in ("reference", "prepared")
        target = targets[(n, label)]
        assert row["kind"] == "s3_crosscurve_single_target_run_v1"
        assert row["target_sha256"] == target["point_sha256"]
        assert row["workload_id"] == target["workload_id"]
        candidate = plan["curves"][str(n)][
            "reference_candidate_id" if arm == "reference" else "rank243_candidate_id"]
        assert row["candidate_id"] == candidate
        assert row["run_id"].startswith(f"{candidate}W{target['workload_id']}R")
        assert row["run_id"] not in run_ids
        run_ids.add(row["run_id"])
        binary_hashes[arm].add(row["binary_sha256"])
        report = row["report"]
        if report is not None:
            assert row["raw_sha256"] is not None
            timing = report.get("timing_ms", {})
            actual_phase_ok = (all(isinstance(timing.get(phase), (int, float))
                                   and timing[phase] >= 0 for phase in PHASES)
                               and isinstance(timing.get("target_online_after_reusable_setup"),
                                              (int, float))
                               and abs(sum(timing[phase] for phase in PHASES)
                                       - timing["target_online_after_reusable_setup"]) < 0.02)
            assert row["phase_ok"] == actual_phase_ok
            if row["verified_native_and_fixture"]:
                assert row["exit_code"] == 0 and actual_phase_ok
                assert report["group_verified"] and report["target_count"] == 1
                assert report["target"] == target["public_point"]
                assert report["factor_base_points"] == plan["curves"][str(n)]["factor_base_points"]
                assert report["recovered_scalar"] == target["fixture_scalar"]
        else:
            assert not row["verified_native_and_fixture"]
        grouped[(n, label)].append(row)
    assert all(len(hashes) == 1 for hashes in binary_hashes.values())
    indices = {int(label[1:]) for _, label in grouped}
    assert len(grouped) == 2 * len(indices)
    assert set(grouped) == {(n, f"T{i:03}") for i in indices for n in (41, 53)}
    blocks = []
    for (n, label), block in sorted(grouped.items()):
        assert len(block) == 4 and Counter(row["arm"] for row in block) == {
            "reference": 2, "prepared": 2}
        expected_ids = {f"{plan['curves'][str(n)]['reference_candidate_id']}W{targets[(n, label)]['workload_id']}R{i}"
                        for i in (1, 2)} | {
            f"{plan['curves'][str(n)]['rank243_candidate_id']}W{targets[(n, label)]['workload_id']}R{i}"
            for i in (1, 2)}
        assert {row["run_id"] for row in block} == expected_ids
        good = all(row["verified_native_and_fixture"] and row["phase_ok"] for row in block)
        result = {"n": n, "label": label, "workload_id": targets[(n, label)]["workload_id"],
                  "all_verified": good, "run_ids": [row["run_id"] for row in block],
                  "failed_run_ids": [row["run_id"] for row in block
                                     if not row["verified_native_and_fixture"]],
                  "exploratory_online_ratio": None, "online_root_call_ratio": None,
                  "prepared_aa_abs_log_ratio": None, "fallback_attempts": None,
                  "online_ms_by_arm": None, "online_root_calls_by_arm": None}
        if good:
            by_arm = {arm: [row for row in block if row["arm"] == arm]
                      for arm in ("reference", "prepared")}
            times = {arm: [row["report"]["timing_ms"]["target_online_after_reusable_setup"]
                           for row in selected] for arm, selected in by_arm.items()}
            calls = {arm: [root_calls(row) for row in selected]
                     for arm, selected in by_arm.items()}
            assert all(value > 0 for values in times.values() for value in values)
            assert all(value > 0 for values in calls.values() for value in values)
            result.update({
                "exploratory_online_ratio": math.sqrt(math.prod(times["reference"]) /
                                                      math.prod(times["prepared"])),
                "online_root_call_ratio": math.sqrt(math.prod(calls["reference"]) /
                                                    math.prod(calls["prepared"])),
                "prepared_aa_abs_log_ratio": abs(math.log(times["prepared"][0] /
                                                       times["prepared"][1])),
                "fallback_attempts": [row["report"]["fallback_attempts"]
                                      for row in by_arm["prepared"]],
                "online_ms_by_arm": times, "online_root_calls_by_arm": calls,
            })
        blocks.append(result)
    curves = {}
    for n in (41, 53):
        selected = [block for block in blocks if block["n"] == n]
        valid = [block for block in selected if block["all_verified"]]
        ratios = [block["exploratory_online_ratio"] for block in valid]
        ops = [block["online_root_call_ratio"] for block in valid]
        aas = sorted(block["prepared_aa_abs_log_ratio"] for block in valid)
        margins = [math.log(ratio / 2) for ratio in ratios]
        mean_margin = statistics.mean(margins) if margins else None
        needed_for_two = (max(24, math.ceil(((Z95_ONE_SIDED + Z80) *
                                            statistics.stdev(margins) / mean_margin) ** 2))
                          if len(margins) >= 4 and mean_margin > 0 else None)
        curves[str(n)] = {
            "target_count": len(selected), "verified_target_count": len(valid),
            "failed_run_count": sum(len(block["failed_run_ids"]) for block in selected),
            "at_least_2_online_ratio_count": sum(ratio >= 2 for ratio in ratios),
            "exploratory_online_ratio_geomean": (
                math.exp(statistics.mean(map(math.log, ratios))) if ratios else None),
            "exploratory_online_ratio_min": min(ratios) if ratios else None,
            "online_root_call_ratio_geomean": (
                math.exp(statistics.mean(map(math.log, ops))) if ops else None),
            "online_root_call_ratio_min": min(ops) if ops else None,
            "prepared_aa_p95_abs_log_ratio": (
                aas[max(0, math.ceil(.95 * len(aas)) - 1)] if aas else None),
            "fallback_attempts_total": sum(sum(block["fallback_attempts"]) for block in valid),
            "provisional_pairs_for_geomean_above_2_one_sided_95_80_percent_power":
                needed_for_two,
        }
    by_index = {(block["n"], int(block["label"][1:])): block for block in blocks}
    matched_differences = [
        math.log(by_index[41, index]["exploratory_online_ratio"] /
                 by_index[53, index]["exploratory_online_ratio"])
        for index in sorted(indices)
        if by_index[41, index]["all_verified"] and by_index[53, index]["all_verified"]]
    mean_difference = statistics.mean(matched_differences) if matched_differences else None
    needed_for_difference = (
        max(24, math.ceil(((Z95_TWO_SIDED + Z80) *
                           statistics.stdev(matched_differences) / abs(mean_difference)) ** 2))
        if len(matched_differences) >= 4 and mean_difference != 0 else None)
    result = {"kind": "s3_crosscurve_local_pilot_analysis_v1",
              "plan_sha256": sha(plan_bytes), "ledger_sha256": sha(raw_bytes),
              "run_count": len(rows), "matched_index_count": len(indices),
              "binary_sha256": {arm: next(iter(hashes)) for arm, hashes in binary_hashes.items()},
              "strict_host_isolation_receipt": None, "controlled_speedup": None,
              "sample_size_estimates_are_exploratory": True,
              "matched_curve_difference": {
                  "verified_pairs": len(matched_differences),
                  "ratio_of_curve_geomean_online_ratios": (
                      math.exp(mean_difference) if mean_difference is not None else None),
                  "provisional_pairs_for_difference_two_sided_95_80_percent_power":
                      needed_for_difference,
              },
              "curves": curves, "blocks": blocks}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(output), "runs": len(rows), "curves": curves}, sort_keys=True))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.ledger, args.output)

"""Analyze frozen N53 rank-242 fresh targets and full-rank stress controls."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import statistics

from analyze_holdout import bootstrap_interval, sha, wilson


HERE = Path(__file__).resolve().parent


def summarize(pairs: list[dict]) -> dict:
    count = len(pairs)
    complete = [pair for pair in pairs if pair["all_verified"]]
    result = {"target_count": count, "verified_count": len(complete),
              "all_verified": len(complete) == count}
    for arm in ("shallow", "deep"):
        ratios = [pair[f"{arm}_online_ratio"] for pair in complete]
        at_least_2 = sum(ratio >= 2 for ratio in ratios)
        result[arm] = {
            "fallback_count": sum(pair[f"{arm}_fallback_attempts"] > 0 for pair in complete),
            "fallback_free_wilson_95": wilson(
                sum(pair[f"{arm}_fallback_attempts"] == 0 for pair in complete),
                len(complete)) if complete else None,
            "geometric_mean_online_ratio":
                math.exp(statistics.mean(map(math.log, ratios))) if result["all_verified"] else None,
            "target_bootstrap_95": bootstrap_interval(list(map(math.log, ratios)))
                if result["all_verified"] else None,
            "minimum_online_ratio": min(ratios) if result["all_verified"] else None,
            "at_least_2_count": at_least_2,
            "at_least_2_wilson_95": wilson(at_least_2, len(complete)) if complete else None,
            "below_2_labels": [pair["label"] for pair in complete
                               if pair[f"{arm}_online_ratio"] < 2],
            "minimum_online_s3_root_call_ratio": min(
                pair[f"{arm}_online_s3_root_call_ratio"] for pair in complete)
                if result["all_verified"] else None,
            "reference_cold_s3_root_calls_equal_count": sum(
                pair[f"{arm}_cold_s3_root_calls_equal_reference"] for pair in complete),
            "median_online_ms": statistics.median(
                pair[f"{arm}_online_ms"] for pair in complete) if result["all_verified"] else None,
            "median_cold_ms": statistics.median(
                pair[f"{arm}_cold_ms"] for pair in complete) if result["all_verified"] else None,
        }
    result["median_reference_online_ms"] = statistics.median(
        pair["reference_online_ms"] for pair in complete) if result["all_verified"] else None
    result["median_reference_cold_ms"] = statistics.median(
        pair["reference_cold_ms"] for pair in complete) if result["all_verified"] else None
    return result


def main(ledger: Path, archive: Path, output: Path, plan_path: Path) -> None:
    if output.exists():
        raise FileExistsError(output)
    plan = json.loads(plan_path.read_text())
    rows = [json.loads(line) for line in ledger.read_text().splitlines()]
    assert len(rows) == 81
    by_target = {}
    for row in rows:
        by_target.setdefault(row["label"], {})[row["arm"]] = row
    assert set(by_target) == {target["label"] for target in plan["targets"]}
    pairs = []
    for target in plan["targets"]:
        label = target["label"]
        arms = by_target[label]
        assert set(arms) == {"reference", "shallow", "deep"}
        verified = all(row["verified_native_and_fixture"] for row in arms.values())
        if verified:
            reference = arms["reference"]["report"]
            shallow = arms["shallow"]["report"]
            deep = arms["deep"]["report"]
            assert len({row["target_sha256"] for row in arms.values()}) == 1
            assert all(row["target_sha256"] == target["point_sha256"] for row in arms.values())
            assert reference["target"] == shallow["target"] == deep["target"] == target["public_point"]
            assert reference["recovered_scalar"] == shallow["recovered_scalar"] == deep["recovered_scalar"] == target["fixture_scalar"]
            assert len({report["factor_base_digest"] for report in (reference, shallow, deep)}) == 1
            assert shallow["prepare_rank_requested"] == plan["rank_policies"]["shallow"]
            assert deep["prepare_rank_requested"] == plan["rank_policies"]["deep"]
            ref_trace = reference["rank_relation_witnesses"]
            ref_online_calls = reference["rank_target_s3_calls"] + reference["target_s3_calls"]
            pair = {"label": label, "family": target["family"], "all_verified": True,
                    "reference_online_ms": reference["timing_ms"]["target_online_after_reusable_setup"],
                    "reference_cold_ms": reference["timing_ms"]["cold_after_launch"],
                    "reference_rank": reference["rank"],
                    "reference_relation_attempts": reference["rank_attempts"],
                    "reference_online_s3_root_calls": ref_online_calls,
                    "reference_raw_sha256": arms["reference"]["raw_sha256"]}
            for name, report in (("shallow", shallow), ("deep", deep)):
                trace = report["rank_relation_witnesses"] + report["fallback_relation_witnesses"]
                assert ref_trace[:min(len(ref_trace), len(trace))] == trace[:min(len(ref_trace), len(trace))]
                online_calls = report["fallback_s3_calls"] + report["target_s3_calls"]
                cold_calls = report["rank_target_s3_calls"] + online_calls
                online = report["timing_ms"]["target_online_after_reusable_setup"]
                pair.update({f"{name}_online_ms": online,
                             f"{name}_cold_ms": report["timing_ms"]["cold_after_launch"],
                             f"{name}_online_ratio": pair["reference_online_ms"] / online,
                             f"{name}_online_s3_root_call_ratio": ref_online_calls / online_calls,
                             f"{name}_cold_s3_root_calls_equal_reference": cold_calls == ref_online_calls,
                             f"{name}_fallback_attempts": report["fallback_attempts"],
                             f"{name}_prepared_attempts": report["rank_attempts"],
                             f"{name}_prepared_rank": report["rank"],
                             f"{name}_final_rank": report["final_rank"],
                             f"{name}_raw_sha256": arms[name]["raw_sha256"]})
        else:
            pair = {"label": label, "family": target["family"], "all_verified": False,
                    "raw_sha256": {name: row["raw_sha256"] for name, row in arms.items()},
                    "errors": {name: row["error"] for name, row in arms.items()}}
        pairs.append(pair)
    fresh = [pair for pair in pairs if pair["family"] == "fresh"]
    stress = [pair for pair in pairs if pair["family"] == "full_rank_stress_control"]
    assert len(fresh) == 24 and len(stress) == 3
    result = {"kind": "s3_partial_rank_n53_fresh_and_stress_analysis_v1",
              "frozen_plan_sha256": sha(plan_path),
              "ledger_sha256": sha(ledger), "raw_archive_sha256": sha(archive),
              "run_count": len(rows),
              "failure_count": sum(not row["verified_native_and_fixture"] for row in rows),
              "strict_host_isolation_receipt": None, "controlled_speedup": None,
              "fresh": summarize(fresh), "stress": summarize(stress), "pairs": pairs}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"fresh": result["fresh"], "stress": result["stress"],
                      "output": str(output)}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", required=True, type=Path)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--plan", type=Path,
                        default=HERE / "rank242/FREEZE.json")
    args = parser.parse_args()
    main(args.ledger, args.archive, args.output, args.plan)

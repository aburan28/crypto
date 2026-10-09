"""Analyze verified one-target group-lift/parallel comparison blocks."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import random
import statistics


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
GROUPLIFT = ROOT / "experiments/koblitz-s3-grouplift-20261009/FREEZE.json"
Z95 = 1.959963984540054
Z80 = 0.8416212335729143


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def geomean(values: list[float]) -> float | None:
    return math.exp(statistics.mean(map(math.log, values))) if values else None


def bootstrap_geomean_95(values: list[float], seed: int) -> list[float] | None:
    if len(values) < 4:
        return None
    logs = list(map(math.log, values))
    rng = random.Random(seed)
    samples = sorted(math.exp(statistics.mean(rng.choices(logs, k=len(logs))))
                     for _ in range(10_000))
    return [samples[249], samples[9749]]


def wilson_95(successes: int, total: int) -> list[float] | None:
    if not total:
        return None
    p = successes / total
    z2 = Z95 * Z95
    center = (p + z2 / (2 * total)) / (1 + z2 / total)
    radius = (Z95 * math.sqrt(p * (1 - p) / total + z2 / (4 * total * total)) /
              (1 + z2 / total))
    return [center - radius, center + radius]


def exploratory_pairs_for_2x(log_ratios: list[float]) -> int | None:
    if len(log_ratios) < 4:
        return None
    effect = statistics.mean(log_ratios) - math.log(2)
    if effect <= 0:
        return None
    return max(24, math.ceil(((Z95 + Z80) * statistics.stdev(log_ratios) / effect) ** 2))


def analyze(ledger: Path, output: Path) -> dict:
    plan_bytes = (HERE / "FREEZE.json").read_bytes()
    plan = json.loads(plan_bytes)
    parent_bytes = PARENT.read_bytes()
    parent = json.loads(parent_bytes)
    assert plan["parent_plan_sha256"] == sha(parent_bytes)
    assert plan["source_sha256"] == sha(
        (ROOT / "examples/s3_partial_prepared_parallel_targetspan.rs").read_bytes())
    old_bytes = GROUPLIFT.read_bytes()
    old = json.loads(old_bytes)
    assert plan["grouplift_plan_sha256"] == sha(old_bytes)
    target_by_key = {(item["n"], item["label"]): item for item in parent["targets"]}
    ledger_bytes = ledger.read_bytes()
    runs = [json.loads(line) for line in ledger_bytes.splitlines()]
    assert runs and len(runs) % 8 == 0
    by_target = defaultdict(list)
    seen_run_ids = set()
    for row in runs:
        n, label, arm = row["n"], row["label"], row["arm"]
        assert row["kind"] == "s3_grouplift_parallel_single_target_run_v1"
        assert n in (41, 53) and arm in ("grouplift", "parallel")
        target = target_by_key[(n, label)]
        assert 1 <= int(label[1:]) <= 24
        assert row["target_sha256"] == target["point_sha256"]
        assert row["workload_id"] == target["workload_id"]
        cid = plan["curves"][str(n)][f"{arm}_candidate_id"]
        assert row["candidate_id"] == cid
        assert row["run_id"].startswith(f"{cid}W{target['workload_id']}R")
        assert row["run_id"] not in seen_run_ids
        seen_run_ids.add(row["run_id"])
        expected_binary = (old["local_binary_sha256"] if arm == "grouplift"
                           else plan["local_binary_sha256"])
        assert row["binary_sha256"] == expected_binary
        if row["verified_native_and_fixture"]:
            report = row["report"]
            assert row["exit_code"] == 0 and row["phase_ok"]
            assert report["group_verified"] and report["target_count"] == 1
            assert report["target"] == target["public_point"]
            assert report["recovered_scalar"] == target["fixture_scalar"]
            assert report["prepared_rank_reached"] and report["prepared_rank"] == 243
            if arm in ("grouplift", "parallel"):
                assert report["group_lift_valid_states"] == report["regular_states"]
                assert report["group_lift_fallback_states"] == 0
            if arm == "parallel":
                workers = plan["target_workers_by_curve"][str(n)]
                assert report["target_query_workers"] == workers
                assert 0 <= report["target_query_winner_worker"] < workers
        by_target[n, label].append(row)
    labels = {int(label[1:]) for _, label in by_target}
    assert len(by_target) == 2 * len(labels)
    assert set(by_target) == {(n, f"T{index:03}") for index in labels for n in (41, 53)}
    blocks = []
    for (n, label), rows in sorted(by_target.items()):
        assert len(rows) == 4 and Counter(row["arm"] for row in rows) == {
            "grouplift": 2, "parallel": 2}
        assert {row["run_id"] for row in rows} == {
            f"{plan['curves'][str(n)][f'{arm}_candidate_id']}W{target_by_key[n, label]['workload_id']}R{rep}"
            for arm in ("grouplift", "parallel") for rep in (1, 2)}
        good = all(row["verified_native_and_fixture"] and row["phase_ok"] for row in rows)
        block = {"n": n, "label": label, "all_verified": good,
                 "run_ids": [row["run_id"] for row in rows],
                 "failed_run_ids": [row["run_id"] for row in rows
                                    if not row["verified_native_and_fixture"]],
                 "online_ratio_grouplift_over_parallel": None,
                 "target_pdp_ratio_grouplift_over_parallel": None,
                 "cold_ratio_grouplift_over_parallel": None,
                 "root_call_ratio_grouplift_over_parallel": None,
                 "grouplift_aa_abs_log_ratio": None,
                 "parallel_aa_abs_log_ratio": None,
                 "rank_relation_witnesses_match": None,
                 "target_relation_indices_equal": None}
        if good:
            arm_rows = {arm: [row for row in rows if row["arm"] == arm]
                        for arm in ("grouplift", "parallel")}
            reports = {arm: [row["report"] for row in selected]
                       for arm, selected in arm_rows.items()}
            first_witnesses = {
                arm: [[(w["relation_scalar"], w["point_indices"])
                       for w in report["rank_relation_witnesses"]]
                      for report in values]
                for arm, values in reports.items()}
            matching = all(w == first_witnesses["grouplift"][0]
                           for arm in first_witnesses for w in first_witnesses[arm])
            assert matching, f"rank witnesses changed for n{n} {label}"
            target_relations_equal = len({tuple(report["target_relation_indices"])
                                          for values in reports.values()
                                          for report in values}) == 1
            metric = {}
            for name, key in (("online", "target_online_after_reusable_setup"),
                              ("target_pdp", "target_pdp"),
                              ("cold", "cold_after_launch")):
                values = {arm: [report["timing_ms"][key] for report in selected]
                          for arm, selected in reports.items()}
                assert all(value > 0 for pair in values.values() for value in pair)
                metric[name] = math.sqrt(math.prod(values["grouplift"]) /
                                         math.prod(values["parallel"]))
            calls = {arm: [report["target_s3_calls"] + report["fallback_s3_calls"]
                           for report in selected] for arm, selected in reports.items()}
            assert all(value > 0 for pair in calls.values() for value in pair)
            online_times = {arm: [report["timing_ms"]["target_online_after_reusable_setup"]
                                  for report in selected]
                            for arm, selected in reports.items()}
            block.update({
                "online_ratio_grouplift_over_parallel": metric["online"],
                "target_pdp_ratio_grouplift_over_parallel": metric["target_pdp"],
                "cold_ratio_grouplift_over_parallel": metric["cold"],
                "root_call_ratio_grouplift_over_parallel": (
                    math.sqrt(math.prod(calls["grouplift"]) / math.prod(calls["parallel"]))),
                "grouplift_aa_abs_log_ratio": abs(math.log(
                    online_times["grouplift"][0] / online_times["grouplift"][1])),
                "parallel_aa_abs_log_ratio": abs(math.log(
                    online_times["parallel"][0] / online_times["parallel"][1])),
                "rank_relation_witnesses_match": True,
                "target_relation_indices_equal": target_relations_equal,
                "online_ms_by_arm": online_times,
                "target_s3_calls_by_arm": calls,
            })
        blocks.append(block)
    curves = {}
    for n in (41, 53):
        selected = [block for block in blocks if block["n"] == n]
        valid = [block for block in selected if block["all_verified"]]
        ratios = [block["online_ratio_grouplift_over_parallel"] for block in valid]
        two_x_count = sum(value >= 2 for value in ratios)
        curves[str(n)] = {"target_count": len(selected),
                          "verified_target_count": len(valid),
                          "failed_run_count": sum(len(block["failed_run_ids"])
                                                  for block in selected),
                          "online_ratio_geomean_exploratory": geomean(ratios),
                          "online_ratio_min_exploratory": min(ratios) if ratios else None,
                          "online_ratio_target_bootstrap_95_descriptive": (
                              bootstrap_geomean_95(ratios, 20261009 + n)),
                          "at_least_2_online_ratio_count": two_x_count,
                          "at_least_2_wilson_95_descriptive": wilson_95(two_x_count, len(selected)),
                          "exploratory_pairs_for_geomean_above_2x": (
                              exploratory_pairs_for_2x(list(map(math.log, ratios)))
                              if len(valid) == len(selected) else None),
                          "target_pdp_ratio_geomean_exploratory": geomean([
                              block["target_pdp_ratio_grouplift_over_parallel"] for block in valid]),
                          "cold_ratio_geomean_exploratory": geomean([
                              block["cold_ratio_grouplift_over_parallel"] for block in valid]),
                          "root_call_ratio_geomean": geomean([
                              block["root_call_ratio_grouplift_over_parallel"] for block in valid]),
                          "grouplift_aa_p95_abs_log_ratio": (
                              sorted(block["grouplift_aa_abs_log_ratio"] for block in valid)[
                                  math.ceil(.95 * len(valid)) - 1] if valid else None),
                          "parallel_aa_p95_abs_log_ratio": (
                              sorted(block["parallel_aa_abs_log_ratio"] for block in valid)[
                                  math.ceil(.95 * len(valid)) - 1] if valid else None)}
    blocks_by_key = {(block["n"], block["label"]): block for block in blocks}
    paired_logs = [
        math.log(blocks_by_key[53, f"T{index:03}"]["online_ratio_grouplift_over_parallel"])
        - math.log(blocks_by_key[41, f"T{index:03}"]["online_ratio_grouplift_over_parallel"])
        for index in sorted(labels)
        if blocks_by_key[41, f"T{index:03}"]["all_verified"]
        and blocks_by_key[53, f"T{index:03}"]["all_verified"]
    ]
    crosscurve_effect = statistics.mean(paired_logs) if paired_logs else None
    crosscurve = {
        "paired_label_count": len(paired_logs),
        "n53_over_n41_geometric_mean_ratio_of_arm_effects_exploratory": (
            math.exp(crosscurve_effect) if crosscurve_effect is not None else None),
        "paired_label_bootstrap_95_descriptive": (
            bootstrap_geomean_95(list(map(math.exp, paired_logs)), 20261009)),
        "exploratory_pairs_for_positive_difference": (
            max(24, math.ceil(((Z95 + Z80) * statistics.stdev(paired_logs) /
                                crosscurve_effect) ** 2))
            if len(paired_logs) >= 4 and crosscurve_effect is not None
            and crosscurve_effect > 0 else None),
        "worker_policies_differ_by_curve": True,
    }
    result = {"kind": "s3_grouplift_parallel_paired_pilot_analysis_v1",
              "plan_sha256": sha(plan_bytes), "ledger_sha256": sha(ledger_bytes),
              "run_count": len(runs), "matched_index_count": len(labels),
              "strict_host_isolation_receipt": None, "controlled_speedup": None,
              "curves": curves, "crosscurve": crosscurve, "blocks": blocks}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(output), "curves": curves,
                      "crosscurve": crosscurve}, sort_keys=True))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.ledger, args.output)

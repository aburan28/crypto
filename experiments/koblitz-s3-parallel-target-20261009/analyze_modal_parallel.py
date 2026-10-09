"""Analyze paired cooperative S3 Modal panels with guest fidelity and power gates.

This keeps every raw failure in the input files. Numerical inference is
descriptive until an independent host-wide CPU isolation receipt is supplied.
"""

from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import random
import statistics


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FIDELITY = ROOT / "experiments/koblitz-s3-partial-prepared-20261009/analyze_modal_fidelity.py"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
GROUPLIFT = ROOT / "experiments/koblitz-s3-grouplift-20261009/FREEZE.json"
spec = importlib.util.spec_from_file_location("s3_guest_fidelity", FIDELITY)
fidelity = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fidelity)
Z95_TWO_SIDED = 1.959963984540054
Z95_ONE_SIDED = 1.6448536269514722
Z80 = 0.8416212335729143
BOOTSTRAPS = 20_000


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_panel(path: Path, plan: dict, plan_sha256: str) -> tuple[dict, dict]:
    payload = path.read_bytes()
    document = json.loads(gzip.decompress(payload) if payload[:2] == b"\x1f\x8b" else payload)
    assert document["kind"] == "modal_vm_s3_parallel_target_panel_v1"
    assert document["frozen_plan_sha256"] == plan_sha256
    assert document["source_sha256"]["parallel"] == plan["source_sha256"]
    assert sha(GROUPLIFT.read_bytes()) == plan["grouplift_plan_sha256"]
    assert document["source_sha256"]["grouplift"] == json.loads(
        GROUPLIFT.read_text())["source_sha256"]
    assert document["strict_host_isolation_receipt"] is None
    assert document["controlled_speedup"] is None
    count = document["target_count_per_curve"]
    assert 1 <= count <= 24 and 1 <= document["first_target_index"] <= 24
    assert document["first_target_index"] + count <= 25
    assert len(document["blocks"]) == 2 * count
    assert len(document["runs"]) == 8 * count
    return document, {"path": str(path), "sha256": sha(payload),
                      "first_target_index": document["first_target_index"],
                      "target_count_per_curve": count}


def hardware_profile(document: dict) -> dict:
    host = document["host"]
    lscpu = host.get("lscpu", {}).get("stdout")
    fields = {}
    if lscpu:
        try:
            fields = {row["field"].strip().rstrip(":"): row["data"]
                      for row in json.loads(lscpu).get("lscpu", [])}
        except (ValueError, KeyError, TypeError):
            pass
    selected = host["selection"]["selected"]
    frequencies = host.get("frequency_by_cpu", {})
    policies = sorted({
        ((frequencies.get(str(cpu), {}).get("scaling_governor") or "").strip(),
         (frequencies.get(str(cpu), {}).get("scaling_min_freq") or "").strip(),
         (frequencies.get(str(cpu), {}).get("scaling_max_freq") or "").strip())
        for cpu in selected})
    return {"model": fields.get("Model name"),
            "architecture": fields.get("Architecture"),
            "numa_nodes": fields.get("NUMA node(s)"),
            "platform": host.get("platform"),
            "frequency_policies": policies,
            "rustc": host.get("rustc", {}).get("stdout"),
            "binary_sha256": document["binary_sha256"]}


def percentile(sorted_values: list[float], fraction: float) -> float:
    return sorted_values[max(0, math.ceil(fraction * len(sorted_values)) - 1)]


def bootstrap_geomean(values: list[float], seed: int) -> list[float] | None:
    if len(values) < 4:
        return None
    logs = list(map(math.log, values))
    rng = random.Random(seed)
    samples = sorted(math.exp(statistics.mean(rng.choices(logs, k=len(logs))))
                     for _ in range(BOOTSTRAPS))
    return [samples[499], samples[19_499]]


def hierarchical_bootstrap_geomean(blocks: list[dict], seed: int) -> list[float] | None:
    """Resample public points and the two executions within each selected arm."""
    if len(blocks) < 4:
        return None
    rng = random.Random(seed)
    samples = []
    for _ in range(BOOTSTRAPS):
        effects = []
        for block in rng.choices(blocks, k=len(blocks)):
            logs = {arm: list(map(math.log, block["online_ms_by_arm"][arm]))
                    for arm in ("grouplift", "parallel")}
            effects.append(statistics.mean(rng.choices(logs["grouplift"], k=2))
                           - statistics.mean(rng.choices(logs["parallel"], k=2)))
        samples.append(math.exp(statistics.mean(effects)))
    samples.sort()
    return [samples[499], samples[19_499]]


def repeat_count_plan(blocks: list[dict]) -> dict | None:
    """Estimate target and within-target log-time variance for repeat choices."""
    if len(blocks) < 4:
        return None
    effects = []
    within = {arm: [] for arm in ("grouplift", "parallel")}
    for block in blocks:
        logs = {arm: list(map(math.log, block["online_ms_by_arm"][arm]))
                for arm in within}
        effects.append(statistics.mean(logs["grouplift"])
                       - statistics.mean(logs["parallel"]))
        for arm in within:
            within[arm].append(statistics.variance(logs[arm]))
    observed = statistics.variance(effects)
    within_by_arm = {arm: statistics.mean(values) for arm, values in within.items()}
    within_sum = sum(within_by_arm.values())
    between = max(0.0, observed - within_sum / 2)
    margin = statistics.mean(effects) - math.log(2)
    choices = {}
    for repeats in (2, 4, 8, 16):
        required = (max(24, math.ceil(((Z95_TWO_SIDED + Z80) ** 2)
                                      * (between + within_sum / repeats) / margin ** 2))
                    if margin > 0 else None)
        choices[str(repeats)] = {
            "target_pairs_required_exploratory": required,
            "one_target_solves_for_this_curve": (
                2 * repeats * required if required is not None else None),
        }
    return {"effect_mean_log_ratio": statistics.mean(effects),
            "margin_above_log_2": margin,
            "observed_between_target_effect_variance_at_two_repeats": observed,
            "estimated_true_between_target_variance": between,
            "within_target_log_time_variance_by_arm": within_by_arm,
            "assumption": "independent targets and repetitions, stable target law and host",
            "95_percent_two_sided_80_percent_power_by_repeats": choices}


def required_pairs_for_log_effect(log_effects: list[float], one_sided: bool) -> int | None:
    if len(log_effects) < 4:
        return None
    effect = statistics.mean(log_effects)
    if effect <= 0:
        return None
    z = Z95_ONE_SIDED if one_sided else Z95_TWO_SIDED
    return max(24, math.ceil(((z + Z80) * statistics.stdev(log_effects) / effect) ** 2))


def analyze(paths: list[Path], output: Path) -> dict:
    assert paths
    plan_path = HERE / "FREEZE.json"
    plan_bytes = plan_path.read_bytes()
    plan = json.loads(plan_bytes)
    assert plan["kind"] == "s3_parallel_matched_curve_pilot_v1"
    parent_bytes = PARENT.read_bytes()
    assert plan["parent_plan_sha256"] == sha(parent_bytes)
    parent = json.loads(parent_bytes)
    assert parent["target_count_per_curve"] == 200
    parent_targets = {(item["n"], item["label"]): item for item in parent["targets"]}
    targets = {}
    for item in plan["targets"]:
        target = parent_targets[item["n"], item["label"]]
        assert item["point_sha256"] == target["point_sha256"]
        assert item["workload_id"] == target["workload_id"]
        targets[item["n"], item["label"]] = target
    assert len(targets) == 48
    panel_records = [load_panel(path, plan, sha(plan_bytes)) for path in paths]
    documents = [document for document, _ in panel_records]
    provenance = [record for _, record in panel_records]
    hardware = [hardware_profile(document) for document in documents]
    same_hardware = (all(item == hardware[0] for item in hardware)
                     and all(item["model"] and item["architecture"] for item in hardware))
    seen_indices = set()
    analyzed_runs = []
    analyzed_blocks = []
    all_perf = []
    for document_index, document in enumerate(documents):
        first = document["first_target_index"]
        count = document["target_count_per_curve"]
        indices = set(range(first, first + count))
        assert not seen_indices.intersection(indices), "overlapping panel indices"
        seen_indices.update(indices)
        by_tag = {row["tag"]: row for row in document["runs"]}
        assert len(by_tag) == len(document["runs"])
        used_tags = []
        analyzed_runs.extend({"panel_index": document_index, "n": row["n"],
                              **fidelity.run_quality(row, document)}
                             for row in document["runs"])
        for perf in document.get("perf_whole_process", []):
            one = fidelity.perf_ipc({"perf_whole_process": [perf]})[0]
            all_perf.append({"panel_index": document_index, "n": perf["n"], **one})
        for block in document["blocks"]:
            n, label = block["n"], block["label"]
            assert n in (41, 53)
            target = targets[(n, label)]
            assert first <= int(label[1:]) < first + count
            assert block["workload_id"] == target["workload_id"]
            assert block["point"] == target["public_point"]
            rows = [by_tag[tag] for tag in block["run_tags"]]
            used_tags.extend(block["run_tags"])
            assert len(rows) == 4 and Counter(row["arm"] for row in rows) == {
                "grouplift": 2, "parallel": 2}
            assert tuple(row["arm"] for row in rows) == tuple(block["order"])
            assert tuple(block["order"]) in (("grouplift", "parallel", "parallel", "grouplift"),
                                             ("parallel", "grouplift", "grouplift", "parallel"))
            assert len({row["run_id"] for row in rows}) == 4
            repetitions = Counter()
            for row in rows:
                assert row["n"] == n and row["label"] == label
                assert row["point"] == target["public_point"]
                assert row["target_sha256"] == target["point_sha256"]
                assert row["workload_id"] == target["workload_id"]
                expected_binary = document["binary_sha256"][row["arm"]]
                assert row["binary_sha256"] == expected_binary
                if row["raw_text"] is not None:
                    assert row["raw_sha256"] == sha(row["raw_text"].encode())
                repetitions[row["arm"]] += 1
                candidate = plan["curves"][str(n)][
                    "grouplift_candidate_id" if row["arm"] == "grouplift" else "parallel_candidate_id"]
                assert row["candidate_id"] == candidate
                assert row["run_id"] == f"{candidate}W{target['workload_id']}R{repetitions[row['arm']]}"
                if row["verified_native_and_fixture"]:
                    report = json.loads(row["raw_text"])
                    assert report["group_verified"] and report["target_count"] == 1
                    assert report["recovered_scalar"] == target["fixture_scalar"]
                    assert report["target"] == target["public_point"]
                    assert report["prepared_rank_reached"] and report["prepared_rank"] == 243
                    if row["arm"] == "parallel":
                        assert report["target_query_workers"] == (
                            plan["target_workers_by_curve"][str(n)])
            good = all(row["verified_native_and_fixture"] and row["phase_ok"] for row in rows)
            ratio = aa = by_arm = calls_by_arm = root_aa = None
            if good:
                reports = [(row["arm"], json.loads(row["raw_text"])) for row in rows]
                by_arm = {arm: [report["timing_ms"]["target_online_after_reusable_setup"]
                                for name, report in reports if name == arm]
                          for arm in ("grouplift", "parallel")}
                calls_by_arm = {arm: [report["target_s3_calls"] for name, report in reports
                                      if name == arm]
                                for arm in by_arm}
                assert all(value > 0 for values in by_arm.values() for value in values)
                assert all(value > 0 for values in calls_by_arm.values() for value in values)
                ratio = math.sqrt(math.prod(by_arm["grouplift"]) / math.prod(by_arm["parallel"]))
                aa = {arm: abs(math.log(values[0] / values[1]))
                      for arm, values in by_arm.items()}
                root_aa = {arm: abs(math.log(values[0] / values[1]))
                           for arm, values in calls_by_arm.items()}
                assert math.isclose(ratio, block["exploratory_online_ratio"], rel_tol=1e-10)
                assert all(math.isclose(aa[arm], block["aa_abs_log_ratio_by_arm"][arm],
                                        rel_tol=1e-10) for arm in aa)
            analyzed_blocks.append({"panel_index": document_index, "n": n, "label": label,
                                    "workload_id": target["workload_id"],
                                    "all_verified": good, "online_ratio": ratio,
                                    "aa_abs_log_ratio_by_arm": aa,
                                    "online_ms_by_arm": by_arm,
                                    "target_s3_calls_by_arm": calls_by_arm,
                                    "aa_abs_log_root_calls_by_arm": root_aa,
                                    "run_ids": [row["run_id"] for row in rows]})
        assert len(used_tags) == len(set(used_tags)) == len(by_tag)
        assert set(used_tags) == set(by_tag)
        assert {(block["n"], block["label"]) for block in document["blocks"]} == {
            (n, f"T{index:03}") for index in indices for n in (41, 53)}
    assert len(analyzed_blocks) == 2 * len(seen_indices)
    curves = {}
    for n in (41, 53):
        blocks = [block for block in analyzed_blocks if block["n"] == n]
        verified_blocks = [block for block in blocks if block["all_verified"]]
        scored_runs = [row for row in analyzed_runs if row["n"] == n]
        ratios = [block["online_ratio"] for block in blocks if block["online_ratio"] is not None]
        aas = {arm: sorted(block["aa_abs_log_ratio_by_arm"][arm] for block in blocks
                           if block["aa_abs_log_ratio_by_arm"] is not None)
               for arm in ("grouplift", "parallel")}
        root_aas = {arm: sorted(block["aa_abs_log_root_calls_by_arm"][arm]
                                for block in verified_blocks)
                    for arm in ("grouplift", "parallel")}
        log_margin = [math.log(value) - math.log(2) for value in ratios]
        curves[str(n)] = {
            "candidate_id": plan["curves"][str(n)]["parallel_candidate_id"],
            "reference_candidate_id": plan["curves"][str(n)]["grouplift_candidate_id"],
            "selected_plan_targets_per_curve": plan["target_count_per_curve"],
            "parent_frozen_pool_targets_per_curve": parent["target_count_per_curve"],
            "aggregate_poolable_on_recorded_hardware": same_hardware,
            "guest_fidelity_score_0_to_100": (
                min(row["guest_score_0_to_100"] for row in scored_runs) if scored_runs else None),
            "target_count": len(blocks), "verified_count": len(ratios),
            "at_least_2_count": sum(value >= 2 for value in ratios),
            "at_least_2_wilson_95": fidelity.wilson(sum(value >= 2 for value in ratios), len(blocks)),
            "minimum_online_ratio": min(ratios) if ratios else None,
            "geometric_mean_online_ratio": (
                math.exp(statistics.mean(map(math.log, ratios)))
                if ratios and same_hardware else None),
            "target_bootstrap_95_descriptive": (
                bootstrap_geomean(ratios, 20261009 + n) if same_hardware else None),
            "target_and_repeat_hierarchical_bootstrap_95_descriptive": (
                hierarchical_bootstrap_geomean(verified_blocks, 20261109 + n)
                if same_hardware else None),
            "aa_p95_abs_log_ratio_by_arm": {
                arm: percentile(values, .95) if values else None
                for arm, values in aas.items()},
            "root_call_aa_p95_abs_log_ratio_by_arm": {
                arm: percentile(values, .95) if values else None
                for arm, values in root_aas.items()},
            "aa_10_percent_environment_gate": (
                percentile(aas["grouplift"], .95) <= fidelity.AA_ABS_LOG_CAP
                if aas["grouplift"] else None),
            "aa_environment_gate_arm": "grouplift_deterministic_target_query",
            "pilot_pairs_for_geomean_above_2_one_sided_95_80_percent_power":
                required_pairs_for_log_effect(log_margin, one_sided=True)
                if same_hardware else None,
            "target_and_repeat_count_plan_exploratory": (
                repeat_count_plan(verified_blocks) if same_hardware else None),
        }
    by_index = {(block["n"], int(block["label"][1:])): block for block in analyzed_blocks}
    paired = []
    for index in sorted(seen_indices):
        left, right = by_index[(41, index)], by_index[(53, index)]
        if left["online_ratio"] is not None and right["online_ratio"] is not None:
            paired.append({"index": index,
                           "log_ratio_53_minus_41":
                               math.log(right["online_ratio"]) - math.log(left["online_ratio"])})
    differences = [item["log_ratio_53_minus_41"] for item in paired]
    signed_differences = (differences if differences and statistics.mean(differences) > 0
                          else [-value for value in differences])
    required_difference = (required_pairs_for_log_effect(signed_differences, one_sided=False)
                           if same_hardware else None)
    curve_difference = {
        "verified_matched_index_count": len(paired),
        "n53_over_n41_ratio_of_geometric_mean_speedups":
            math.exp(statistics.mean(differences)) if differences and same_hardware else None,
        "matched_index_bootstrap_95_descriptive":
            bootstrap_geomean([math.exp(value) for value in differences], 20261009)
            if same_hardware else None,
        "pilot_pairs_per_curve_for_difference_two_sided_95_80_percent_power": required_difference,
        "selected_plan_targets_per_curve": plan["target_count_per_curve"],
        "parent_frozen_pool_targets_per_curve": parent["target_count_per_curve"],
    }
    required = curve_difference["pilot_pairs_per_curve_for_difference_two_sided_95_80_percent_power"]
    curve_difference["selected_plan_covers_estimate"] = (
        required <= plan["target_count_per_curve"] if required is not None else None)
    curve_difference["parent_pool_covers_estimate"] = (
        required <= parent["target_count_per_curve"] if required is not None else None)
    guest_score = min(row["guest_score_0_to_100"] for row in analyzed_runs)
    result = {
        "kind": "s3_parallel_modal_fidelity_analysis_v1",
        "plan_sha256": sha(plan_bytes),
        "panel_inputs": provenance,
        "same_recorded_hardware_and_binaries": same_hardware,
        "hardware_profiles": hardware,
        "guest_fidelity_score_0_to_100": guest_score,
        "guest_score_definition": "minimum run score; each run scores passed checks among the 14 guest checks, with unavailable checks scored zero",
        "strict_host_isolation_receipt": None,
        "controlled_speedup": None,
        "curves": curves,
        "curve_difference": curve_difference,
        "runs": analyzed_runs,
        "blocks": analyzed_blocks,
        "whole_process_perf": all_perf,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(output), "score": guest_score,
                      "targets_per_curve": len(seen_indices),
                      "verified_per_curve": {n: curves[n]["verified_count"] for n in curves},
                      "difference_pilot_pairs": required}, sort_keys=True))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.inputs, args.output)

"""Analyze paired N41/N53 Modal panels with guest fidelity and power gates.

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
    assert document["kind"] == "modal_vm_s3_crosscurve_rank243_panel_v1"
    assert document["frozen_plan_sha256"] == plan_sha256
    assert document["source_sha256"]["reference"] == plan["reference_source_sha256"]
    assert document["source_sha256"]["partial"] == plan["partial_source_sha256"]
    assert document["strict_host_isolation_receipt"] is None
    assert document["controlled_speedup"] is None
    count = document["target_count_per_curve"]
    assert 1 <= count <= 24 and 1 <= document["first_target_index"] <= 200
    assert document["first_target_index"] + count <= 201
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
    assert plan["kind"] == "s3_rank243_matched_curve_fresh_panel_v1"
    targets = {(item["n"], item["label"]): item for item in plan["targets"]}
    assert len(targets) == 400
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
                "reference": 2, "prepared": 2}
            assert tuple(row["arm"] for row in rows) == tuple(block["order"])
            assert tuple(block["order"]) in (("reference", "prepared", "prepared", "reference"),
                                             ("prepared", "reference", "reference", "prepared"))
            assert len({row["run_id"] for row in rows}) == 4
            repetitions = Counter()
            for row in rows:
                assert row["n"] == n and row["label"] == label
                assert row["point"] == target["public_point"]
                assert row["target_sha256"] == target["point_sha256"]
                assert row["workload_id"] == target["workload_id"]
                expected_binary = document["binary_sha256"][
                    "reference" if row["arm"] == "reference" else "partial"]
                assert row["binary_sha256"] == expected_binary
                if row["raw_text"] is not None:
                    assert row["raw_sha256"] == sha(row["raw_text"].encode())
                repetitions[row["arm"]] += 1
                candidate = plan["curves"][str(n)][
                    "reference_candidate_id" if row["arm"] == "reference" else "rank243_candidate_id"]
                assert row["candidate_id"] == candidate
                assert row["run_id"] == f"{candidate}W{target['workload_id']}R{repetitions[row['arm']]}"
                if row["verified_native_and_fixture"]:
                    report = json.loads(row["raw_text"])
                    assert report["group_verified"] and report["target_count"] == 1
                    assert report["recovered_scalar"] == target["fixture_scalar"]
                    assert report["target"] == target["public_point"]
            good = all(row["verified_native_and_fixture"] and row["phase_ok"] for row in rows)
            ratio = aa = None
            if good:
                times = [(row["arm"], json.loads(row["raw_text"])["timing_ms"][
                    "target_online_after_reusable_setup"]) for row in rows]
                by_arm = {arm: [value for name, value in times if name == arm]
                          for arm in ("reference", "prepared")}
                assert all(value > 0 for values in by_arm.values() for value in values)
                ratio = math.sqrt(math.prod(by_arm["reference"]) / math.prod(by_arm["prepared"]))
                aa = abs(math.log(by_arm["prepared"][0] / by_arm["prepared"][1]))
                assert math.isclose(ratio, block["exploratory_online_ratio"], rel_tol=1e-10)
                assert math.isclose(aa, block["prepared_aa_abs_log_ratio"], rel_tol=1e-10)
            analyzed_blocks.append({"panel_index": document_index, "n": n, "label": label,
                                    "workload_id": target["workload_id"],
                                    "all_verified": good, "online_ratio": ratio,
                                    "prepared_aa_abs_log_ratio": aa,
                                    "run_ids": [row["run_id"] for row in rows]})
        assert len(used_tags) == len(set(used_tags)) == len(by_tag)
        assert set(used_tags) == set(by_tag)
        assert {(block["n"], block["label"]) for block in document["blocks"]} == {
            (n, f"T{index:03}") for index in indices for n in (41, 53)}
    assert len(analyzed_blocks) == 2 * len(seen_indices)
    curves = {}
    for n in (41, 53):
        blocks = [block for block in analyzed_blocks if block["n"] == n]
        scored_runs = [row for row in analyzed_runs if row["n"] == n]
        ratios = [block["online_ratio"] for block in blocks if block["online_ratio"] is not None]
        aas = sorted(block["prepared_aa_abs_log_ratio"] for block in blocks
                     if block["prepared_aa_abs_log_ratio"] is not None)
        log_margin = [math.log(value) - math.log(2) for value in ratios]
        curves[str(n)] = {
            "candidate_id": plan["curves"][str(n)]["rank243_candidate_id"],
            "reference_candidate_id": plan["curves"][str(n)]["reference_candidate_id"],
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
            "prepared_aa_p95_abs_log_ratio": percentile(aas, .95) if aas else None,
            "aa_10_percent_gate": percentile(aas, .95) <= fidelity.AA_ABS_LOG_CAP if aas else None,
            "pilot_pairs_for_geomean_above_2_one_sided_95_80_percent_power":
                required_pairs_for_log_effect(log_margin, one_sided=True)
                if same_hardware else None,
        }
    by_index = {(block["n"], int(block["label"][1:])): block for block in analyzed_blocks}
    paired = []
    for index in sorted(seen_indices):
        left, right = by_index[(41, index)], by_index[(53, index)]
        if left["online_ratio"] is not None and right["online_ratio"] is not None:
            paired.append({"index": index,
                           "log_ratio_41_minus_53":
                               math.log(left["online_ratio"]) - math.log(right["online_ratio"])})
    differences = [item["log_ratio_41_minus_53"] for item in paired]
    signed_differences = (differences if differences and statistics.mean(differences) > 0
                          else [-value for value in differences])
    required_difference = (required_pairs_for_log_effect(signed_differences, one_sided=False)
                           if same_hardware else None)
    curve_difference = {
        "verified_matched_index_count": len(paired),
        "ratio_of_curve_geometric_mean_speedups":
            math.exp(statistics.mean(differences)) if differences and same_hardware else None,
        "matched_index_bootstrap_95_descriptive":
            bootstrap_geomean([math.exp(value) for value in differences], 20261009)
            if same_hardware else None,
        "pilot_pairs_per_curve_for_difference_two_sided_95_80_percent_power": required_difference,
        "frozen_pool_size_per_curve": plan["target_count_per_curve"],
    }
    required = curve_difference["pilot_pairs_per_curve_for_difference_two_sided_95_80_percent_power"]
    curve_difference["frozen_pool_covers_estimate"] = required <= 200 if required is not None else None
    guest_score = min(row["guest_score_0_to_100"] for row in analyzed_runs)
    result = {
        "kind": "s3_crosscurve_modal_fidelity_analysis_v1",
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

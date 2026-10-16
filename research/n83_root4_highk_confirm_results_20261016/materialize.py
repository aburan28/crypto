#!/usr/bin/env python3
"""Seal the replayed N83 root4 high-K confirmation."""

import argparse
import hashlib
import importlib.util
import json
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
COMMON = REPO / "research/n83_best_e2e_results_20261015/materialize.py"
PANEL = REPO / "research/n83_factor_base_100_results_20261014/materialize.py"
COMMON_SHA256 = "dc7c42bbc7f5c1d49e11150c16ff89c3bfb322ebd7c7fa0177e8306184b9efeb"
PANEL_SHA256 = "843745c6c74721f35292bea3767b08afbdbd5915ed65e46d3f8ed06809bd7a5e"
SOURCE_COMMIT = "adf9ca9b254ac7de8bb777446e74788d25e369a8"
RUN_ID = 38835419672
SOURCE_SHA256 = "f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc"
GENERATOR_SHA256 = "fdd2e6ad3a9260d7e66f235fa01e2c41a90f35d6cb752e7b1ab2bfcaf3bb8699"
POINT_SHA256 = "1dcb2665d71cf0979ab7df205cd43281369f24db3b4b696fe6d97fb03f603579"
ANSWER_SHA256 = "da9cd25bfe87eedd2a2c74d39f5c7c7ed5eaff296026a580f07aba76d2e444a8"
ADMISSION_SHA256 = "5c14104122027a0877f2bc02d2aa8a56bcd3b9e613ea414d9f26336f1570879d"
RSS_CAP_BYTES = 61 * 1024**3 // 4
ARMS = (("root4k3500std_a", 3500, "root4", "standard"),
        ("root4k3596std_a", 3596, "root4", "standard"),
        ("root4k3596std_b", 3596, "root4", "standard"),
        ("root4k3500std_b", 3500, "root4", "standard"))

if (hashlib.sha256(COMMON.read_bytes()).hexdigest() != COMMON_SHA256
        or hashlib.sha256(PANEL.read_bytes()).hexdigest() != PANEL_SHA256):
    raise ValueError("reusable N83 result-sealing implementation changed")
spec = importlib.util.spec_from_file_location("n83_direct_materializer", COMMON)
common = importlib.util.module_from_spec(spec)
spec.loader.exec_module(common)
common.SOURCE_COMMIT = SOURCE_COMMIT
common.panel.SOURCE_COMMIT = SOURCE_COMMIT
common.SOURCE_SHA256 = SOURCE_SHA256
common.panel.SOURCE_SHA256 = SOURCE_SHA256
common.BASE_HASHES = {}
common.ROOT_KEYS = {}


def frozen_points():
    generator = common.git_file("research/n83_root4_highk_confirm_20261016/generate.py")
    if hashlib.sha256(generator).hexdigest() != GENERATOR_SHA256:
        raise ValueError("executed fixture generator differs")
    files = (REPO / "research/n83_root4_highk_confirm_20261016/targets.jsonl",
             REPO / "research/n83_root4_highk_confirm_20261016/answers.jsonl")
    blobs = [path.read_bytes() for path in files]
    if [hashlib.sha256(blob).hexdigest() for blob in blobs] != [POINT_SHA256, ANSWER_SHA256]:
        raise ValueError("regenerated fixture differs from source run")
    points, answers = ([json.loads(line) for line in blob.splitlines() if line]
                       for blob in blobs)
    if (len(points) != 8000 or len(answers) != 8000
            or len({tuple(point) for point in points}) != 8000
            or any(row["fixture_index"] != index or row["point"] != point
                   for index, (point, row) in enumerate(zip(points, answers)))):
        raise ValueError("public-point count, uniqueness or sidecar pairing differs")
    return points


def candidate(k, variant, profile, summary, dump_path, field, curve):
    if profile == "wide3":
        return common.candidate(k, variant, profile, summary, dump_path,
                                field, curve)
    if variant != "root4" or profile != "standard":
        raise ValueError("unexpected root4 Bloom pipeline")
    panel = common.panel
    record = panel.candidate(k, variant, summary, panel.read(dump_path),
                             dump_path, field, curve)
    if (summary["sparse_canonical_enabled"] is not True
            or summary["root_bloom_enabled"] is not True
            or (summary["root_bloom_profile"], summary["root_bloom_hashes"],
                summary["root_bloom_bits_per_request"], summary["root_bloom_bytes"])
            != ("standard", 3, 8, 1073741824)
            or summary["base_hash"] != common.BASE_HASHES[k]
            or summary["root_table_entries"] != common.ROOT_KEYS[k]):
        raise ValueError(f"K{k}: standard Bloom, base or index differs")
    record["implementation"]["algorithm_flags"].update(
        {"sparse_canonical": True, "root_bloom": True,
         "root_bloom_profile": "standard"})
    record["point_decomposition"]["canonical_root_selector"] = "sparse_linear_map_exact"
    record["point_decomposition"]["root_bloom_filter"] = {
        "profile": "standard", "hashes": 3, "bits_per_request": 8,
        "allocated_bytes": 1073741824,
        "policy": "membership_prefilter_then_exact_root_lookup"}
    body = {key: value for key, value in record.items() if key != "candidate_id"}
    record["candidate_id"] = (
        f"IC1N83Ckb1a1fb{166*k}PDP4rootRCguidedLAgaussTDdirectISO0h"
        + panel.digest(body)[:12])
    return record


def seal(root, artifact_zip, artifact_sha256, analysis_sha256, artifact_id):
    panel = common.panel
    if (artifact_id <= 0 or panel.file_sha(artifact_zip) != artifact_sha256
            or panel.file_sha(root / "analysis.json") != analysis_sha256):
        raise ValueError("archived ZIP or local replay digest differs")
    analysis = panel.read(root / "analysis.json")
    source = panel.read(root / "source.json")
    host = panel.read(root / "host.json")
    receipt = panel.read(root.parent / "fixture-check.json")
    names = [arm for arm, _, _, _ in ARMS]
    hashes = source["hashes"]
    if (host["git_sha"] != SOURCE_COMMIT
            or source["schema"] != "n83.root4-highk-confirm-source/v1"
            or source["arm_order"] != names
            or hashes["examples/koblitz_orbit_dlp_fast.rs"] != SOURCE_SHA256
            or hashes["Cargo.lock"] != common.LOCK_SHA256
            or hashes["research/n83_root4_highk_confirm_20261016/targets.jsonl"] != POINT_SHA256
            or hashes["research/n83_root4_highk_confirm_20261016/answers.jsonl"] != ANSWER_SHA256
            or hashes["research/n83_root4_memory_frontier_results_20261015/selection.json"] != ADMISSION_SHA256
            or receipt != {"count": 8000, "targets_sha256": POINT_SHA256,
                           "answers_sha256": ANSWER_SHA256}
            or source["panel_runner_sha256"] != panel.file_sha(
                REPO / "research/n83_root4_highk_confirm_20261016/run.py")
            or analysis["schema"] != "n83.root4-highk-confirm-analysis/v1"
            or analysis["hardware_tier"] != "cloud_instance_exploratory"
            or analysis["target_count"] != 8000
            or analysis["arm_order"] != names
            or len(analysis["arms"]) != len(ARMS)):
        raise ValueError("source, input, host or arm identity differs")
    for k in (3500, 3596):
        accepted = [row for row in analysis["arms"]
                    if row.get("k") == k and row["status"] == "PASS_INDEPENDENT_REPLAY"]
        if accepted:
            bases = {row["base_hash"] for row in accepted}
            counts = {row["root_table_entries"] for row in accepted}
            if len(bases) != 1 or len(counts) != 1:
                raise ValueError(f"K{k}: repeated base or root-index identity differs")
            common.BASE_HASHES[k] = bases.pop()
            common.ROOT_KEYS[k] = counts.pop()
    points = frozen_points()
    workload = {
        "schema_version": 1, "curve_id": panel.CURVE_ID,
        "target_count": 8000, "public_targets": points,
        "input_law": {
            "domain": "n83-root4-highk-confirm-fresh-8000-20261016-v1",
            "scalar_rule": "1+(little_endian(SHA256(domain||u32_le(index))) mod (r-1))",
            "answer_sidecar_sha256": ANSWER_SHA256},
        "target_sha256": POINT_SHA256,
        "cold_or_warm": "four_mirrored_fresh_processes_with_empty_index_and_logs; reusable_setup_outside_per_target_online_intervals",
    }
    workload_id = panel.digest(workload)[:12]
    workload_record = common.write_gzip_once(HERE / "workload.json.gz",
                                             {"workload_id": workload_id, **workload})
    field, curve = panel.curve_record()
    verified = {row["arm"]: row for row in analysis["arms"]}
    if set(verified) != set(names):
        raise ValueError("verifier arm set differs")
    rows = []
    for arm, k, variant, profile in ARMS:
        item = verified[arm]
        arm_root = root / "runs" / arm
        usage_path = arm_root / "usage.json"
        usage = panel.read(usage_path) if usage_path.exists() else {}
        row = {
            "arm": arm, "candidate_id": None, "workload_id": workload_id,
            "run_id": None, "status": item["status"],
            "hardware_tier": "cloud_instance_exploratory",
            "workload_role": "secondary_fresh_8000_point_balanced_highk_root4_confirmation",
            "controlled_speedup": None, "source_commit": SOURCE_COMMIT,
            "executed_binary_sha256": source["binary_sha256"],
            "host_receipt_sha256": panel.file_sha(root / "host.json"),
            "source_receipt_sha256": panel.file_sha(root / "source.json"),
            "correctness_receipt_sha256": analysis_sha256,
            "full_process_wall_seconds": usage.get("wall_seconds"),
            "peak_rss_bytes": usage.get("peak_rss_bytes"),
            "sampled_peak_rss_bytes": usage.get("sampled_peak_rss_bytes"),
            "rss_guard_bytes": RSS_CAP_BYTES,
            "exit_status": usage.get("exit_status"),
            "stop_reason": usage.get("stop_reason"),
            "major_faults": usage.get("major_faults"),
            "involuntary_context_switches": usage.get("involuntary_context_switches"),
        }
        if item["status"] == "PASS_INDEPENDENT_REPLAY":
            summary = panel.read(arm_root / "stdout.jsonl")
            dump_path = arm_root / "base_dump.jsonl"
            manifest = candidate(k, variant, profile, summary, dump_path,
                                 field, curve)
            manifest_path = HERE / f"candidate_K{k}_{variant}_{profile}.json"
            panel.write_once(manifest_path, manifest)
            targets = [json.loads(line) for line in
                       (arm_root / "targets_result.jsonl").read_text().splitlines() if line]
            if len(targets) != 8000:
                raise ValueError(f"{arm}: expected 8000 target records")
            candidate_id = manifest["candidate_id"]
            boundary = [float(target["target_ms"]) - sum(float(target[key]) for key in
                        ("target_query_ms", "target_pdp_ms",
                         "target_relation_check_ms", "target_descent_ms",
                         "target_recovery_check_ms")) for target in targets]
            if any(not 0 <= value <= 1 for value in boundary):
                raise ValueError(f"{arm}: target phase boundary overhead differs")
            row.update({
                "candidate_id": candidate_id,
                "run_id": candidate_id + "W" + workload_id +
                          ("R1" if arm.endswith("_a") else "R2"),
                "candidate_manifest_sha256": panel.file_sha(manifest_path),
                "actual_factor_base_points": item["actual_factor_base_points"],
                "effective_columns": k, "final_rank": summary["rank"],
                "source_base_hash_blake3": item["base_hash"],
                "root_table_entries": item["root_table_entries"],
                "root_table_allocated_bytes": item["root_table_allocated_bytes"],
                "root_state_allocated_bytes": item["root_state_allocated_bytes"],
                "root_bloom_bytes": summary["root_bloom_bytes"],
                "target_probes": [target["probes"] for target in targets],
                "ic_online_ms_per_target_exploratory":
                    [float(target["target_ms"]) for target in targets],
                "ic_online_mean_ms_exploratory": item["mean_online_ms_exploratory"],
                "ic_online_phases_ms_per_target": [
                    {key: target[key] for key in
                     ("target_query_ms", "target_pdp_ms", "target_relation_check_ms",
                      "target_descent_ms", "target_recovery_check_ms")}
                    for target in targets],
                "ic_online_boundary_overhead_ms_per_target": [
                    float(target["target_ms"]) - sum(float(target[key]) for key in
                        ("target_query_ms", "target_pdp_ms",
                         "target_relation_check_ms", "target_descent_ms",
                         "target_recovery_check_ms"))
                    for target in targets],
                "setup_ms": item["setup_ms"],
                "rank_ms": summary["timing_ms"]["rank_stage"],
                "final_la_ms": summary["timing_ms"]["linear_algebra"],
                "cpu_throttled_periods_delta": item["cpu_throttled_periods_delta"],
                "oom_kills_delta": item["oom_kills_delta"],
            })
        rows.append(row)
    measurements_record = common.write_gzip_once(HERE / "measurements.json.gz", rows)
    eligible = [row for row in rows if row["status"] == "PASS_INDEPENDENT_REPLAY"
                and row["peak_rss_bytes"] <= RSS_CAP_BYTES
                and row.get("cpu_throttled_periods_delta") == 0
                and row.get("oom_kills_delta") == 0]
    complete = len(eligible) == len(ARMS) and analysis["status"] == "COMPLETE"
    comparison = analysis["balanced_comparison"]
    if complete != (comparison is not None):
        raise ValueError("balanced selection eligibility differs")
    selected_pipeline = comparison["selected_pipeline"] if complete else None
    if complete:
        by_arm = {row["arm"]: row for row in eligible}
        means = {}
        for pipeline in ("root4k3500std", "root4k3596std"):
            a = by_arm[pipeline + "_a"]["ic_online_ms_per_target_exploratory"]
            b = by_arm[pipeline + "_b"]["ic_online_ms_per_target_exploratory"]
            mean = statistics.mean((x + y) / 2 for x, y in zip(a, b))
            if abs(mean - comparison["pipeline_mean_ms"][pipeline]) > 1e-9:
                raise ValueError(f"{pipeline}: paired mean differs")
            means[pipeline] = mean
        if selected_pipeline is not None and selected_pipeline != min(means, key=means.get):
            raise ValueError("selected pipeline is not lowest mean")
    selected_row = next((row for row in eligible
                         if row["arm"] == f"{selected_pipeline}_a"), None)
    panel.write_once(HERE / "selection.json", {
        "schema": "n83.root4-highk-confirm-selection/v1",
        "status": "PANEL_SELECTED" if selected_row else
                  "BOUNDED_UNKNOWN" if complete else "INCOMPLETE",
        "selected_pipeline": selected_pipeline,
        "descriptive_lowest_pipeline": (
            min(comparison["pipeline_mean_ms"],
                key=comparison["pipeline_mean_ms"].get) if complete else None),
        "candidate_id": selected_row["candidate_id"] if selected_row else None,
        "candidate_manifest_sha256": selected_row["candidate_manifest_sha256"] if selected_row else None,
        "source_base_hash_blake3": selected_row["source_base_hash_blake3"] if selected_row else None,
        "balanced_comparison": comparison,
        "workload_id": workload_id, "artifact_run_id": RUN_ID,
        "artifact_id": artifact_id, "artifact_sha256": artifact_sha256,
        "analysis_sha256": analysis_sha256, "source_commit": SOURCE_COMMIT,
        "source_sha256": SOURCE_SHA256,
    })
    panel.write_once(HERE / "records_index.json", {
        "schema": "n83.root4-highk-confirm-compressed-records/v1",
        "workload_id": workload_id, "measurement_rows": len(rows),
        "workload": workload_record, "measurements": measurements_record,
        "artifact_run_id": RUN_ID, "artifact_id": artifact_id,
        "artifact_sha256": artifact_sha256, "analysis_sha256": analysis_sha256,
    })
    print(json.dumps({"workload_id": workload_id,
                      "candidate_ids": [row["candidate_id"] for row in rows],
                      "selected_pipeline": selected_pipeline}, sort_keys=True))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-zip", required=True, type=Path)
    parser.add_argument("--artifact-sha256", required=True)
    parser.add_argument("--analysis-sha256", required=True)
    parser.add_argument("--artifact-id", required=True, type=int)
    parser.add_argument("artifact_results", type=Path)
    args = parser.parse_args()
    seal(args.artifact_results.resolve(), args.artifact_zip.resolve(),
         args.artifact_sha256, args.analysis_sha256, args.artifact_id)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Seal the replayed N83 root4 Bloom factor-base screen."""

import argparse
import hashlib
import importlib.util
import json
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
COMMON = REPO / "research/n83_best_e2e_results_20261015/materialize.py"
SOURCE_COMMIT = "afea7179c40c9e00f8701a32c596e47cd8129f62"
SOURCE_SHA256 = "f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc"
GENERATOR_SHA256 = "327a50e945927c870ee8539228d423e88287605259dbd38ac51ca8145f695269"
POINT_SHA256 = "71b93fd3a397ee2c1204dc14fc451894eee7537b99f50d75c98d65dee6b1e3d1"
ANSWER_SHA256 = "5417c1a10a4c2240eaec304627eddfbd287e1d005f674b1a2b17b3b957aa53e5"
ARTIFACT_SHA256 = "8687a0ad20ec1a0525a437f172801278a1a931771e88bef404c3559c77ed7626"
ANALYSIS_SHA256 = "15b47de167f99e004bccada699f739ff7d0d7e7f7e7960239d63b4b00183e547"
RSS_CAP_BYTES = 29 * 1024**3 // 2
ARMS = (("root4k3350std_a", 3350, "root4", "standard"),
        ("root4k3400std_a", 3400, "root4", "standard"),
        ("root4k3425std_a", 3425, "root4", "standard"),
        ("root11k3300wide3_a", 3300, "root11", "wide3"),
        ("root11k3300wide3_b", 3300, "root11", "wide3"),
        ("root4k3425std_b", 3425, "root4", "standard"),
        ("root4k3400std_b", 3400, "root4", "standard"),
        ("root4k3350std_b", 3350, "root4", "standard"))

spec = importlib.util.spec_from_file_location("n83_direct_materializer", COMMON)
common = importlib.util.module_from_spec(spec)
spec.loader.exec_module(common)
common.SOURCE_COMMIT = SOURCE_COMMIT
common.panel.SOURCE_COMMIT = SOURCE_COMMIT
common.SOURCE_SHA256 = SOURCE_SHA256
common.panel.SOURCE_SHA256 = SOURCE_SHA256
common.BASE_HASHES = {
    3300: "a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b",
    3350: "dc2e0c5185b88b2b0c74ede0bbc7b889d756f2f21558e9b932850d97a1bdd07b",
    3400: "baa4ee439a1571543c9294cde96e4acddb930c38cf81d0b916302dfc5e44aeb4",
    3425: "47e4eb7a1ef7e596f3505e9a1997b60dab56f9f306ad6b969f5c452923f772c6",
}
common.ROOT_KEYS = {3300: 903852075, 3350: 931448963,
                    3400: 959460734, 3425: 973622246}


def frozen_points():
    generator = common.git_file("research/n83_root4_bloom_20261015/generate.py")
    if hashlib.sha256(generator).hexdigest() != GENERATOR_SHA256:
        raise ValueError("executed fixture generator differs")
    files = (REPO / "research/n83_root4_bloom_20261015/targets.jsonl",
             REPO / "research/n83_root4_bloom_20261015/answers.jsonl")
    blobs = [path.read_bytes() for path in files]
    if [hashlib.sha256(blob).hexdigest() for blob in blobs] != [POINT_SHA256, ANSWER_SHA256]:
        raise ValueError("regenerated fixture differs from source run")
    points, answers = ([json.loads(line) for line in blob.splitlines() if line]
                       for blob in blobs)
    if (len(points) != 1000 or len(answers) != 1000
            or len({tuple(point) for point in points}) != 1000
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


def seal(root, artifact_zip):
    panel = common.panel
    if (panel.file_sha(artifact_zip) != ARTIFACT_SHA256
            or panel.file_sha(root / "analysis.json") != ANALYSIS_SHA256):
        raise ValueError("archived ZIP or local replay digest differs")
    analysis = panel.read(root / "analysis.json")
    source = panel.read(root / "source.json")
    host = panel.read(root / "host.json")
    receipt = panel.read(root.parent / "fixture-check.json")
    names = [arm for arm, _, _, _ in ARMS]
    hashes = source["hashes"]
    if (host["git_sha"] != SOURCE_COMMIT
            or source["schema"] != "n83.root4-bloom-source/v1"
            or source["arm_order"] != names
            or hashes["examples/koblitz_orbit_dlp_fast.rs"] != SOURCE_SHA256
            or hashes["Cargo.lock"] != common.LOCK_SHA256
            or hashes["research/n83_root4_bloom_20261015/targets.jsonl"] != POINT_SHA256
            or hashes["research/n83_root4_bloom_20261015/answers.jsonl"] != ANSWER_SHA256
            or receipt != {"count": 1000, "targets_sha256": POINT_SHA256,
                           "answers_sha256": ANSWER_SHA256}
            or source["panel_runner_sha256"] != panel.file_sha(
                REPO / "research/n83_root4_bloom_20261015/run.py")
            or analysis["schema"] != "n83.root4-bloom-analysis/v1"
            or analysis["hardware_tier"] != "cloud_instance_exploratory"
            or analysis["target_count"] != 1000
            or analysis["arm_order"] != names
            or len(analysis["arms"]) != len(ARMS)):
        raise ValueError("source, input, host or arm identity differs")
    points = frozen_points()
    workload = {
        "schema_version": 1, "curve_id": panel.CURVE_ID,
        "target_count": 1000, "public_targets": points,
        "input_law": {
            "domain": "n83-root4-bloom-fresh-1000-20261015-v1",
            "scalar_rule": "1+(little_endian(SHA256(domain||u32_le(index))) mod (r-1))",
            "answer_sidecar_sha256": ANSWER_SHA256},
        "target_sha256": POINT_SHA256,
        "cold_or_warm": "eight_mirrored_fresh_processes_with_empty_index_and_logs; reusable_setup_outside_per_target_online_intervals",
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
        usage = panel.read(arm_root / "usage.json")
        row = {
            "arm": arm, "candidate_id": None, "workload_id": workload_id,
            "run_id": None, "status": item["status"],
            "hardware_tier": "cloud_instance_exploratory",
            "workload_role": "secondary_fresh_1000_point_balanced_root4_bloom_factor_base_screen",
            "controlled_speedup": None, "source_commit": SOURCE_COMMIT,
            "executed_binary_sha256": source["binary_sha256"],
            "host_receipt_sha256": panel.file_sha(root / "host.json"),
            "source_receipt_sha256": panel.file_sha(root / "source.json"),
            "correctness_receipt_sha256": ANALYSIS_SHA256,
            "full_process_wall_seconds": usage["wall_seconds"],
            "peak_rss_bytes": usage["peak_rss_bytes"],
            "sampled_peak_rss_bytes": usage["sampled_peak_rss_bytes"],
            "rss_guard_bytes": RSS_CAP_BYTES,
            "exit_status": usage["exit_status"],
            "stop_reason": usage["stop_reason"],
            "major_faults": usage["major_faults"],
            "involuntary_context_switches": usage["involuntary_context_switches"],
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
            if len(targets) != 1000:
                raise ValueError(f"{arm}: expected 1000 target records")
            candidate_id = manifest["candidate_id"]
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
        for pipeline in ("root4k3350std", "root4k3400std",
                         "root4k3425std", "root11k3300wide3"):
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
        "schema": "n83.root4-bloom-selection/v1",
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
        "workload_id": workload_id, "artifact_sha256": ARTIFACT_SHA256,
        "analysis_sha256": ANALYSIS_SHA256, "source_commit": SOURCE_COMMIT,
        "source_sha256": SOURCE_SHA256,
    })
    panel.write_once(HERE / "records_index.json", {
        "schema": "n83.root4-bloom-compressed-records/v1",
        "workload_id": workload_id, "measurement_rows": len(rows),
        "workload": workload_record, "measurements": measurements_record,
        "artifact_sha256": ARTIFACT_SHA256, "analysis_sha256": ANALYSIS_SHA256,
    })
    print(json.dumps({"workload_id": workload_id,
                      "candidate_ids": [row["candidate_id"] for row in rows],
                      "selected_pipeline": selected_pipeline}, sort_keys=True))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-zip", required=True, type=Path)
    parser.add_argument("artifact_results", type=Path)
    args = parser.parse_args()
    seal(args.artifact_results.resolve(), args.artifact_zip.resolve())


if __name__ == "__main__":
    main()

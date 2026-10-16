#!/usr/bin/env python3
"""Seal the independently replayed N83 root11 K3500/K3550 edge screen."""

import argparse
import hashlib
import importlib.util
import json
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
COMMON = REPO / "research/n83_best_e2e_results_20261015/materialize.py"
SOURCE_COMMIT = "10c7c80d9489ef4fd1f302c9063c7452444aeaee"
SOURCE_SHA256 = "f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc"
PRODUCER_SHA256 = "decd62935e518ff769cb1198b9da4d368f7b98ec8169d4f9e73c15d0e902ece9"
GENERATOR_SHA256 = "0fe13bf072f42a0ab9a7811ed76b534bf06030c2dd6dc2de2a01b548d7fa159b"
POINT_SHA256 = "c75cbb92ef1a23057272cf3481cc649b12495135c75affd1c2610205a42b7aef"
ANSWER_SHA256 = "ca95d57ff091b8fc7de8aa2fd22120cadac94a7c7a697e3c757aebc53115a159"
ARTIFACT_SHA256 = "a7bbcefaff56322e234db96ea799771941dab020252e0adc89a6c66190ba62e3"
ANALYSIS_SHA256 = "f50192442355a64fd1c4cfd07147f5b35f8aef4ae77d783f1c7614035992c3eb"
ARTIFACT_ID = 12005586706
RUN_ID = 38794757719
RSS_CAP_BYTES = 14 * 1024**3
ARMS = (("root11k3450wide3_a", 3450), ("root11k3500wide3_a", 3500),
        ("root11k3550wide3_a", 3550), ("root11k3550wide3_b", 3550),
        ("root11k3500wide3_b", 3500), ("root11k3450wide3_b", 3450))

spec = importlib.util.spec_from_file_location("n83_direct_materializer", COMMON)
common = importlib.util.module_from_spec(spec)
spec.loader.exec_module(common)
common.SOURCE_COMMIT = SOURCE_COMMIT
common.panel.SOURCE_COMMIT = SOURCE_COMMIT
common.SOURCE_SHA256 = SOURCE_SHA256
common.panel.SOURCE_SHA256 = SOURCE_SHA256
common.BASE_HASHES = {
    3450: "ae75c86565f5fa12fb6f30a5539f6f4638bcf1ca8b347a5cefa20f5073a9f572",
    3500: "3387d0c40335c49256ad2048542455572980bef501c3735fd8a26f42ac035e67",
    3550: "1dfb31b3e371e6b846947fa28aa1ceb1ed440fd2c89135d9f8b5eb1505c32548",
}
common.ROOT_KEYS = {3450: 987887478, 3500: 1016729255, 3550: 1045986026}


def points():
    generator = common.git_file("research/n83_root11_edgek_20261015/generate.py")
    if hashlib.sha256(generator).hexdigest() != GENERATOR_SHA256:
        raise ValueError("executed fixture generator differs")
    files = (REPO / "research/n83_root11_edgek_20261015/targets.jsonl",
             REPO / "research/n83_root11_edgek_20261015/answers.jsonl")
    blobs = [path.read_bytes() for path in files]
    if [hashlib.sha256(blob).hexdigest() for blob in blobs] != [POINT_SHA256, ANSWER_SHA256]:
        raise ValueError("frozen public points or sidecar differs")
    public, answers = ([json.loads(line) for line in blob.splitlines() if line]
                       for blob in blobs)
    if (len(public) != 1000 or len(answers) != 1000
            or len({tuple(point) for point in public}) != 1000
            or any(row["fixture_index"] != i or row["point"] != point
                   for i, (point, row) in enumerate(zip(public, answers)))):
        raise ValueError("public-point count, uniqueness or sidecar pairing differs")
    return public


def seal(root, artifact_zip):
    panel = common.panel
    if (panel.file_sha(artifact_zip) != ARTIFACT_SHA256
            or panel.file_sha(root / "analysis.json") != ANALYSIS_SHA256):
        raise ValueError("archived ZIP or independent replay digest differs")
    analysis = panel.read(root / "analysis.json")
    source = panel.read(root / "source.json")
    host = panel.read(root / "host.json")
    fixture = panel.read(root.parent / "fixture-check.json")
    names = [arm for arm, _ in ARMS]
    hashes = source["hashes"]
    if (host["git_sha"] != SOURCE_COMMIT
            or source["schema"] != "n83.root11-edgek-source/v1"
            or source["binary_sha256"] != PRODUCER_SHA256
            or source["arm_order"] != names
            or hashes["examples/koblitz_orbit_dlp_fast.rs"] != SOURCE_SHA256
            or hashes["Cargo.lock"] != common.LOCK_SHA256
            or hashes["research/n83_root11_edgek_20261015/targets.jsonl"] != POINT_SHA256
            or hashes["research/n83_root11_edgek_20261015/answers.jsonl"] != ANSWER_SHA256
            or fixture != {"count": 1000, "targets_sha256": POINT_SHA256,
                           "answers_sha256": ANSWER_SHA256}
            or source["panel_runner_sha256"] != panel.file_sha(
                REPO / "research/n83_root11_edgek_20261015/run.py")
            or analysis["schema"] != "n83.root11-edgek-analysis/v1"
            or analysis["hardware_tier"] != "cloud_instance_exploratory"
            or analysis["target_count"] != 1000
            or analysis["arm_order"] != names
            or len(analysis["arms"]) != len(ARMS)):
        raise ValueError("source, public fixture, host or arm identity differs")
    public = points()
    workload = {
        "schema_version": 1, "curve_id": panel.CURVE_ID,
        "target_count": 1000, "public_targets": public,
        "input_law": {
            "domain": "n83-root11-edgek-fresh-1000-20261015-v1",
            "scalar_rule": "1+(little_endian(SHA256(domain||u32_le(index))) mod (r-1))",
            "answer_sidecar_sha256": ANSWER_SHA256},
        "target_sha256": POINT_SHA256,
        "cold_or_warm": "six_mirrored_fresh_processes_with_empty_index_and_logs; reusable_setup_outside_per_target_online_intervals",
    }
    workload_id = panel.digest(workload)[:12]
    workload_record = common.write_gzip_once(HERE / "workload.json.gz",
                                             {"workload_id": workload_id, **workload})
    field, curve = panel.curve_record()
    verified = {row["arm"]: row for row in analysis["arms"]}
    if set(verified) != set(names):
        raise ValueError("verifier arm set differs")
    rows = []
    for arm, k in ARMS:
        item = verified[arm]
        arm_root = root / "runs" / arm
        usage = panel.read(arm_root / "usage.json")
        row = {
            "arm": arm, "candidate_id": None, "workload_id": workload_id,
            "run_id": None, "status": item["status"],
            "hardware_tier": "cloud_instance_exploratory",
            "workload_role": "secondary_fresh_1000_point_balanced_order_root11_edgek_factor_base_screen",
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
            dump = arm_root / "base_dump.jsonl"
            manifest = common.candidate(k, "root11", "wide3", summary,
                                        dump, field, curve)
            manifest_path = HERE / f"candidate_K{k}_root11_wide3.json"
            panel.write_once(manifest_path, manifest)
            targets = [json.loads(line) for line in
                       (arm_root / "targets_result.jsonl").read_text().splitlines() if line]
            if len(targets) != 1000:
                raise ValueError(f"{arm}: expected 1000 target rows")
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
        raise ValueError("balanced comparison eligibility differs")
    selected = comparison["selected_pipeline"] if complete else None
    if complete:
        by_arm = {row["arm"]: row for row in eligible}
        means = {}
        for pipeline in ("root11k3450wide3", "root11k3500wide3",
                         "root11k3550wide3"):
            a = by_arm[pipeline + "_a"]["ic_online_ms_per_target_exploratory"]
            b = by_arm[pipeline + "_b"]["ic_online_ms_per_target_exploratory"]
            mean = statistics.mean((x + y) / 2 for x, y in zip(a, b))
            if abs(mean - comparison["pipeline_mean_ms"][pipeline]) > 1e-9:
                raise ValueError(f"{pipeline}: paired mean differs")
            means[pipeline] = mean
        if selected is not None and selected != min(means, key=means.get):
            raise ValueError("selected pipeline is not the lowest mean")
    no_new_challenger = (complete and
                         min(comparison["pipeline_mean_ms"][name]
                             for name in ("root11k3500wide3", "root11k3550wide3")) >=
                         comparison["pipeline_mean_ms"]["root11k3450wide3"])
    panel.write_once(HERE / "selection.json", {
        "schema": "n83.root11-edgek-selection/v1",
        "status": "PANEL_SELECTED" if selected else
                  "BOUNDED_UNKNOWN" if complete else "INCOMPLETE",
        "selected_pipeline": selected,
        "no_new_challenger_under_frozen_confirmation_rule": no_new_challenger,
        "balanced_comparison": comparison,
        "workload_id": workload_id,
        "artifact_run_id": RUN_ID, "artifact_id": ARTIFACT_ID,
        "artifact_sha256": ARTIFACT_SHA256, "analysis_sha256": ANALYSIS_SHA256,
        "source_commit": SOURCE_COMMIT, "source_sha256": SOURCE_SHA256,
    })
    panel.write_once(HERE / "records_index.json", {
        "schema": "n83.root11-edgek-compressed-records/v1",
        "workload_id": workload_id, "measurement_rows": len(rows),
        "workload": workload_record, "measurements": measurements_record,
        "artifact_run_id": RUN_ID, "artifact_id": ARTIFACT_ID,
        "artifact_sha256": ARTIFACT_SHA256, "analysis_sha256": ANALYSIS_SHA256,
    })
    print(json.dumps({"workload_id": workload_id,
                      "candidate_ids": [row["candidate_id"] for row in rows],
                      "selected_pipeline": selected,
                      "no_new_challenger": no_new_challenger}, sort_keys=True))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-zip", required=True, type=Path)
    parser.add_argument("artifact_results", type=Path)
    args = parser.parse_args()
    seal(args.artifact_results.resolve(), args.artifact_zip.resolve())


if __name__ == "__main__":
    main()

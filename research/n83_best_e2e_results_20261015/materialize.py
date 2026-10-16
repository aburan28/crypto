#!/usr/bin/env python3
"""Seal the replayed six-arm N83 factor-base comparison."""

import argparse
import gzip
import hashlib
import importlib.util
import io
import json
import statistics
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
PANEL = REPO / "research/n83_factor_base_100_results_20261014/materialize.py"
SOURCE_COMMIT = "ae6a10b641e17629bd282d8205787ede129dbf98"
SOURCE_SHA256 = "f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc"
GENERATOR_SHA256 = "bc4d8dca6d56fa4cfcaac33390b1321218c158f1c05424cd479a48f1c3185551"
LOCK_SHA256 = "ac5cd7a3459be776180afe1eb2975c76b1a88f84bd1e99aa20475927ec260a90"
POINT_SHA256 = "bbb3a2fffdcb51db6105d00f12f9dc58450520a5bae646bd3dfb08d8dc72eab4"
ANSWER_SHA256 = "3bc11ae82d81c86b6ba15491e96cdbefb874b7da6eb63de04f5b745e7bfd432c"
ARTIFACT_SHA256 = "b77c196d6c2fcc3227cf3d0a1e0cacb893210fb53e95d23a2dbd01b398a54851"
ANALYSIS_SHA256 = "59e2a1cd6bf404ed46a4e82694d5ecde61f191f5e97142ad27282948cae568a7"
RSS_CAP_BYTES = 14 * 1024**3
ARMS = (("wide3_a", 3300, "root11", "wide3"),
        ("root4k3550_a", 3550, "root4", "off"),
        ("root4k3596_a", 3596, "root4", "off"),
        ("root4k3596_b", 3596, "root4", "off"),
        ("root4k3550_b", 3550, "root4", "off"),
        ("wide3_b", 3300, "root11", "wide3"))
BASE_HASHES = {
    3300: "a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b",
    3550: "1dfb31b3e371e6b846947fa28aa1ceb1ed440fd2c89135d9f8b5eb1505c32548",
    3596: "9bdadb868c18ccaad3446a62316356e3fe21ecc76512c0c27ccc1fa24d661bc1",
}
ROOT_KEYS = {3300: 903852075, 3550: 1045986026, 3596: 1073268803}

spec = importlib.util.spec_from_file_location("n83_panel_materializer", PANEL)
panel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(panel)
panel.SOURCE_COMMIT = SOURCE_COMMIT
panel.SOURCE_SHA256 = SOURCE_SHA256


def write_gzip_once(path, value):
    raw = (json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n").encode()
    buffer = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buffer,
                       compresslevel=9, mtime=0) as stream:
        stream.write(raw)
    compressed = buffer.getvalue()
    if path.exists():
        if path.read_bytes() != compressed:
            raise ValueError(f"existing sealed record differs: {path}")
    else:
        path.write_bytes(compressed)
    return {"sha256": hashlib.sha256(compressed).hexdigest(),
            "json_sha256": hashlib.sha256(raw).hexdigest(),
            "compressed_bytes": len(compressed), "json_bytes": len(raw)}


def git_file(relative):
    return subprocess.check_output(["git", "show", SOURCE_COMMIT + ":" + relative], cwd=REPO)


def frozen_points():
    names = ("research/n83_best_e2e_20261015/targets.jsonl",
             "research/n83_best_e2e_20261015/answers.jsonl")
    if hashlib.sha256(git_file("research/n83_best_e2e_20261015/generate.py")).hexdigest() != GENERATOR_SHA256:
        raise ValueError("executed fixture generator differs")
    blobs = [(REPO / name).read_bytes() for name in names]
    if ([hashlib.sha256(blob).hexdigest() for blob in blobs]
            != [POINT_SHA256, ANSWER_SHA256]):
        raise ValueError("executed public points or answer sidecar changed")
    points, answers = ([json.loads(line) for line in blob.splitlines() if line]
                       for blob in blobs)
    if (len(points) != 4000 or len(answers) != 4000
            or len({tuple(point) for point in points}) != 4000
            or any(row["fixture_index"] != index or row["point"] != point
                   for index, (point, row) in enumerate(zip(points, answers)))):
        raise ValueError("frozen point count, uniqueness or sidecar pairing differs")
    return points


def candidate(k, variant, profile, summary, dump_path, field, curve):
    record = panel.candidate(k, variant, summary, panel.read(dump_path),
                             dump_path, field, curve)
    if (summary["sparse_canonical_enabled"] is not True
            or summary["base_hash"] != BASE_HASHES[k]
            or summary["root_table_entries"] != ROOT_KEYS[k]
            or summary["root_bloom_enabled"] is not (profile == "wide3")):
        raise ValueError(f"K{k}: source base, root index or selector differs")
    record["implementation"]["algorithm_flags"]["sparse_canonical"] = True
    record["point_decomposition"]["canonical_root_selector"] = "sparse_linear_map_exact"
    if profile == "wide3":
        if (summary["root_bloom_profile"], summary["root_bloom_hashes"],
                summary["root_bloom_bits_per_request"], summary["root_bloom_bytes"]) != (
                    "wide3", 3, 16, 1073741824):
            raise ValueError("wide3 Bloom allocation differs")
        record["implementation"]["algorithm_flags"].update(
            {"root_bloom": True, "root_bloom_profile": "wide3"})
        record["point_decomposition"]["root_bloom_filter"] = {
            "profile": "wide3", "hashes": 3, "bits_per_request": 16,
            "allocated_bytes": 1073741824,
            "policy": "membership_prefilter_then_exact_root_lookup"}
    elif (summary["root_bloom_profile"], summary["root_bloom_hashes"],
          summary["root_bloom_bits_per_request"], summary["root_bloom_bytes"]) != (
              "none", 0, 0, 0):
        raise ValueError("Bloom-off arm has an allocated filter")
    body = {key: value for key, value in record.items() if key != "candidate_id"}
    record["candidate_id"] = (
        f"IC1N83Ckb1a1fb{166*k}PDP4rootRCguidedLAgaussTDdirectISO0h"
        + panel.digest(body)[:12])
    return record


def seal(root, artifact_zip):
    if (panel.file_sha(artifact_zip) != ARTIFACT_SHA256
            or panel.file_sha(root / "analysis.json") != ANALYSIS_SHA256):
        raise ValueError("archived ZIP or independently replayed analysis digest differs")
    analysis = panel.read(root / "analysis.json")
    host = panel.read(root / "host.json")
    source = panel.read(root / "source.json")
    fixture_receipt = panel.read(root.parent / "fixture-check.json")
    hashes = source["hashes"]
    arm_names = [arm for arm, _, _, _ in ARMS]
    if (host["git_sha"] != SOURCE_COMMIT
            or source["schema"] != "n83.best-e2e-source/v1"
            or source["arm_order"] != arm_names
            or hashes["examples/koblitz_orbit_dlp_fast.rs"] != SOURCE_SHA256
            or hashes["Cargo.lock"] != LOCK_SHA256
            or hashes["research/n83_best_e2e_20261015/targets.jsonl"] != POINT_SHA256
            or hashes["research/n83_best_e2e_20261015/answers.jsonl"] != ANSWER_SHA256
            or fixture_receipt != {"count": 4000, "targets_sha256": POINT_SHA256,
                                   "answers_sha256": ANSWER_SHA256}
            or source["panel_runner_sha256"] != panel.file_sha(
                REPO / "research/n83_best_e2e_20261015/run.py")
            or analysis["schema"] != "n83.best-e2e-analysis/v1"
            or analysis["hardware_tier"] != "cloud_instance_exploratory"
            or analysis["target_count"] != 4000
            or analysis["arm_order"] != arm_names
            or len(analysis["arms"]) != len(ARMS)):
        raise ValueError("source, fixture, allocation or arm identity differs")
    points = frozen_points()
    workload = {
        "schema_version": 1, "curve_id": panel.CURVE_ID,
        "target_count": 4000, "public_targets": points,
        "input_law": {
            "domain": "n83-best-e2e-fresh-4000-20261015-v1",
            "scalar_rule": "1+(little_endian(SHA256(domain||u32_le(index))) mod (r-1))",
            "answer_sidecar_sha256": ANSWER_SHA256},
        "target_sha256": POINT_SHA256,
        "cold_or_warm": "six_mirrored_fresh_processes_with_empty_index_and_logs; reusable_setup_outside_per_target_online_intervals",
    }
    workload_id = panel.digest(workload)[:12]
    workload_record = write_gzip_once(HERE / "workload.json.gz",
                                      {"workload_id": workload_id, **workload})
    field, curve = panel.curve_record()
    verified = {row["arm"]: row for row in analysis["arms"]}
    if set(verified) != set(arm_names):
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
            "workload_role": "secondary_fresh_4000_point_balanced_order_three_pipeline_factor_base_comparison",
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
            target_rows = [json.loads(line) for line in
                           (arm_root / "targets_result.jsonl").read_text().splitlines() if line]
            if len(target_rows) != 4000:
                raise ValueError(f"{arm}: expected 4000 target rows")
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
                "target_probes": [target["probes"] for target in target_rows],
                "ic_online_ms_per_target_exploratory":
                    [float(target["target_ms"]) for target in target_rows],
                "ic_online_mean_ms_exploratory": item["mean_online_ms_exploratory"],
                "ic_online_phases_ms_per_target": [
                    {key: target[key] for key in
                     ("target_query_ms", "target_pdp_ms", "target_relation_check_ms",
                      "target_descent_ms", "target_recovery_check_ms")}
                    for target in target_rows],
                "setup_ms": item["setup_ms"],
                "rank_ms": summary["timing_ms"]["rank_stage"],
                "final_la_ms": summary["timing_ms"]["linear_algebra"],
                "cpu_throttled_periods_delta": item["cpu_throttled_periods_delta"],
                "oom_kills_delta": item["oom_kills_delta"],
            })
        rows.append(row)
    measurements_record = write_gzip_once(HERE / "measurements.json.gz", rows)
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
        means = {}
        by_arm = {row["arm"]: row for row in eligible}
        for pipeline in ("wide3", "root4k3550", "root4k3596"):
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
    selection = {
        "schema": "n83.best-e2e-selection/v1",
        "status": "PANEL_SELECTED" if selected_row else
                  "BOUNDED_UNKNOWN" if complete else "INCOMPLETE",
        "selected_pipeline": selected_pipeline,
        "candidate_id": selected_row["candidate_id"] if selected_row else None,
        "candidate_manifest_sha256": selected_row["candidate_manifest_sha256"] if selected_row else None,
        "source_base_hash_blake3": selected_row["source_base_hash_blake3"] if selected_row else None,
        "balanced_comparison": comparison,
        "workload_id": workload_id,
        "artifact_sha256": ARTIFACT_SHA256,
        "analysis_sha256": ANALYSIS_SHA256,
        "source_commit": SOURCE_COMMIT,
        "source_sha256": SOURCE_SHA256,
    }
    panel.write_once(HERE / "selection.json", selection)
    panel.write_once(HERE / "records_index.json", {
        "schema": "n83.best-e2e-compressed-records/v1",
        "workload_id": workload_id, "measurement_rows": len(rows),
        "workload": workload_record, "measurements": measurements_record,
        "artifact_sha256": ARTIFACT_SHA256, "analysis_sha256": ANALYSIS_SHA256,
    })
    print(json.dumps({"workload_id": workload_id,
                      "candidate_ids": [row["candidate_id"] for row in rows],
                      "selected_pipeline": selected_pipeline,
                      "selection_status": selection["status"]}, sort_keys=True))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-zip", required=True, type=Path)
    parser.add_argument("artifact_results", type=Path)
    args = parser.parse_args()
    seal(args.artifact_results.resolve(), args.artifact_zip.resolve())


if __name__ == "__main__":
    main()

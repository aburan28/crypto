#!/usr/bin/env python3
"""Seal the balanced-order N83 K3300 root Bloom profile comparison."""

import argparse
import gzip
import hashlib
import importlib.util
import io
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
PANEL = REPO / "research/n83_factor_base_100_results_20261014/materialize.py"
SOURCE_COMMIT = "42fc71945cd7087f78464abc65c4bae7c52a93e8"
SOURCE_SHA256 = "f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc"
LOCK_SHA256 = "ac5cd7a3459be776180afe1eb2975c76b1a88f84bd1e99aa20475927ec260a90"
POINT_SHA256 = "35e8b43966d15065cf0186ace762f823a602aa1be7d7c0a9039e020550d01491"
FIXTURE_SHA256 = "6e1181ba958708803316903a645463425bc79130074041d39c8094778e2f4379"
EXPECTED_BASE_HASH = "a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b"
EXPECTED_ROOT_ENTRIES = 903852075
PROFILES = {"standard": (3, 8, 536870912), "wide3": (3, 16, 1073741824),
            "wide5": (5, 16, 1073741824)}
RSS_CAP_BYTES = 14 * 1024 ** 3
ARMS = (("standard_a", "standard"), ("wide3_a", "wide3"),
        ("wide5_a", "wide5"), ("wide5_b", "wide5"),
        ("wide3_b", "wide3"), ("standard_b", "standard"))

spec = importlib.util.spec_from_file_location("n83_panel_materializer", PANEL)
panel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(panel)
panel.SOURCE_COMMIT = SOURCE_COMMIT
panel.SOURCE_SHA256 = SOURCE_SHA256


def write_gzip_once(path, value):
    raw = (json.dumps(value, indent=2, sort_keys=True,
                      ensure_ascii=False) + "\n").encode("utf-8")
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


def git_file(path):
    return subprocess.check_output(["git", "show", SOURCE_COMMIT + ":" + path], cwd=REPO)


def frozen_points():
    generator_path = "research/n83_k_abba_20261015/generate.py"
    exclusion_path = "research/n83_k_abba_20261015/excluded_bloom_targets.jsonl"
    for relative, expected in (
            (generator_path, "cc94a9f1dc971058acfc667ba19c4571076aa2b77cffc39553dba2122fca21be"),
            (exclusion_path, "ef616c5254a9d883e962024b9f4bd7f794cfe5576c73a93b099ac5be6fe45668")):
        if hashlib.sha256(git_file(relative)).hexdigest() != expected:
            raise ValueError(f"executed fixture source differs: {relative}")
    point_bytes = (REPO / "research/n83_k_abba_20261015/targets.jsonl").read_bytes()
    answer_bytes = (REPO / "research/n83_k_abba_20261015/answers.jsonl").read_bytes()
    if (hashlib.sha256(point_bytes).hexdigest() != POINT_SHA256
            or hashlib.sha256(answer_bytes).hexdigest() != FIXTURE_SHA256):
        raise ValueError("frozen fresh points or answer sidecar changed")
    points = [json.loads(line) for line in point_bytes.splitlines() if line]
    answers = [json.loads(line) for line in answer_bytes.splitlines() if line]
    if (len(points) != 4000 or len(answers) != 4000
            or len({tuple(point) for point in points}) != 4000
            or any(answer["fixture_index"] != index or answer["point"] != point
                   for index, (point, answer) in enumerate(zip(points, answers)))):
        raise ValueError("fresh point count, uniqueness or answer pairing differs")
    return points


def bloom_candidate(profile, summary, dump_path, field, curve):
    k, variant = 3300, "root11"
    record = panel.candidate(k, variant, summary, panel.read(dump_path),
                             dump_path, field, curve)
    if (summary["sparse_canonical_enabled"] is not True
            or summary["root_bloom_enabled"] is not True
            or summary["root_bloom_profile"] != profile
            or (summary["root_bloom_hashes"],
                summary["root_bloom_bits_per_request"],
                summary["root_bloom_bytes"]) != PROFILES[profile]
            or summary["base_hash"] != EXPECTED_BASE_HASH
            or summary["root_table_entries"] != EXPECTED_ROOT_ENTRIES):
        raise ValueError(f"{profile}: Bloom profile, root index or base differs")
    record["implementation"]["algorithm_flags"]["sparse_canonical"] = True
    record["implementation"]["algorithm_flags"]["root_bloom"] = True
    record["implementation"]["algorithm_flags"]["root_bloom_profile"] = profile
    record["point_decomposition"]["canonical_root_selector"] = "sparse_linear_map_exact"
    record["point_decomposition"]["root_bloom_filter"] = {
        "profile": profile, "hashes": PROFILES[profile][0],
        "bits_per_request": PROFILES[profile][1],
        "allocated_bytes": PROFILES[profile][2],
        "policy": "membership_prefilter_then_exact_root_lookup",
    }
    body = {key: value for key, value in record.items() if key != "candidate_id"}
    record["candidate_id"] = (
        f"IC1N83Ckb1a1fb{166*k}PDP4rootRCguidedLAgaussTDdirectISO0h"
        + panel.digest(body)[:12])
    return record


def seal(root, artifact_zip, artifact_sha256, analysis_sha256):
    if (panel.file_sha(artifact_zip) != artifact_sha256
            or panel.file_sha(root / "analysis.json") != analysis_sha256):
        raise ValueError("raw artifact or independent analysis digest mismatch")
    analysis = panel.read(root / "analysis.json")
    host = panel.read(root / "host.json")
    source = panel.read(root / "source.json")
    hashes = source["hashes"]
    if (host["git_sha"] != SOURCE_COMMIT
            or source["schema"] != "n83.bloom-profile-source/v1"
            or source["arm_order"] != [arm for arm, _ in ARMS]
            or hashes["examples/koblitz_orbit_dlp_fast.rs"] != SOURCE_SHA256
            or hashes["Cargo.lock"] != LOCK_SHA256
            or hashes["research/n83_k_abba_20261015/targets.jsonl"] != POINT_SHA256
            or hashes["research/n83_k_abba_20261015/answers.jsonl"] != FIXTURE_SHA256
            or source["panel_runner_sha256"] != panel.file_sha(
                REPO / "research/n83_bloom_profile_20261015/run.py")
            or analysis["schema"] != "n83.bloom-profile-analysis/v1"
            or analysis["hardware_tier"] != "cloud_instance_exploratory"
            or analysis["workload_role"] != "secondary_reused_4000_point_balanced_order_bloom_profile_comparison"
            or analysis["target_count"] != 4000
            or analysis["arm_order"] != [arm for arm, _ in ARMS]
            or len(analysis["arms"]) != len(ARMS)):
        raise ValueError("source, input or allocation identity mismatch")
    points = frozen_points()
    workload = {
        "schema_version": 1, "curve_id": panel.CURVE_ID,
        "target_count": 4000, "public_targets": points,
        "input_law": {"domain": "n83-k-abba-fresh-4000-20261015-v1",
                      "scalar_rule": "1+(little_endian(SHA256(domain||u32_le(index))) mod (r-1))",
                      "answer_sidecar_sha256": FIXTURE_SHA256},
        "target_sha256": POINT_SHA256,
        "cold_or_warm": "six_fresh_empty_index_and_logs_processes_in_standard_wide3_wide5_wide5_wide3_standard_order; 4000_public_points_reused_per_arm_after_each_setup",
    }
    workload_id = panel.digest(workload)[:12]
    workload_file = write_gzip_once(HERE / "workload.json.gz",
                                    {"workload_id": workload_id, **workload})
    field, curve = panel.curve_record()
    verified_by_arm = {row["arm"]: row for row in analysis["arms"]}
    if set(verified_by_arm) != {arm for arm, _ in ARMS}:
        raise ValueError("analysis arms differ from frozen protocol")
    rows = []
    for arm, profile in ARMS:
        verified = verified_by_arm[arm]
        row = {"arm": arm, "profile": profile,
               "candidate_id": None, "workload_id": workload_id,
               "run_id": None, "status": verified["status"],
               "hardware_tier": "cloud_instance_exploratory",
               "workload_role": "secondary_reused_4000_point_balanced_order_bloom_profile_comparison",
               "selection_status": analysis["status"],
               "controlled_speedup": None, "source_commit": SOURCE_COMMIT,
               "executed_binary_sha256": source["binary_sha256"],
               "host_receipt_sha256": panel.file_sha(root / "host.json"),
               "source_receipt_sha256": panel.file_sha(root / "source.json"),
               "correctness_receipt_sha256": analysis_sha256}
        arm_root = root / "runs" / arm
        usage_path = arm_root / "usage.json"
        if usage_path.exists():
            usage = panel.read(usage_path)
            row.update({"exit_status": usage["exit_status"],
                        "full_process_wall_seconds": usage["wall_seconds"],
                        "peak_rss_bytes": usage["peak_rss_bytes"],
                        "sampled_peak_rss_bytes": usage["sampled_peak_rss_bytes"],
                        "rss_guard_bytes": RSS_CAP_BYTES,
                        "stop_reason": usage["stop_reason"],
                        "major_faults": usage["major_faults"],
                        "involuntary_context_switches": usage["involuntary_context_switches"],
                        "resource_envelope_status": (
                            "OBSERVED_PEAK_ABOVE_14_GIB_GUARD"
                            if usage["peak_rss_bytes"] > RSS_CAP_BYTES else
                            "OBSERVED_PEAK_WITHIN_14_GIB_GUARD")})
        if verified["status"] == "PASS_INDEPENDENT_REPLAY":
            summary = panel.read(arm_root / "stdout.jsonl")
            targets = [json.loads(line) for line in
                       (arm_root / "targets_result.jsonl").read_text().splitlines() if line]
            if len(targets) != 4000:
                raise ValueError(f"{arm}: expected 4000 target rows")
            dump_path = arm_root / "base_dump.jsonl"
            manifest = bloom_candidate(profile, summary, dump_path, field, curve)
            candidate_id = manifest["candidate_id"]
            panel.write_once(HERE / f"candidate_K3300_root11_{profile}.json", manifest)
            run_number = 1 if arm.endswith("_a") else 2
            row.update({
                "candidate_id": candidate_id,
                "run_id": candidate_id + "W" + workload_id + f"R{run_number}",
                "actual_factor_base_points": verified["actual_factor_base_points"],
                "effective_columns": 3300, "final_rank": summary["rank"],
                "root_bloom_hashes": summary["root_bloom_hashes"],
                "root_bloom_bits_per_request": summary["root_bloom_bits_per_request"],
                "root_bloom_allocated_bytes": summary["root_bloom_bytes"],
                "root_table_allocated_bytes": verified["root_table_allocated_bytes"],
                "root_state_allocated_bytes": verified["root_state_allocated_bytes"],
                "regular_states": summary["regular_states"],
                "root_table_entries": verified["root_table_entries"],
                "target_probes": [target["probes"] for target in targets],
                "ic_online_ms_per_target_exploratory":
                    [float(target["target_ms"]) for target in targets],
                "ic_online_mean_ms_exploratory": verified["mean_online_ms_exploratory"],
                "ic_online_phases_ms_per_target": [
                    {key: target[key] for key in
                     ("target_query_ms", "target_pdp_ms", "target_relation_check_ms",
                      "target_descent_ms", "target_recovery_check_ms")}
                    for target in targets],
                "setup_ms": verified["setup_ms"],
                "rank_ms": summary["timing_ms"]["rank_stage"],
                "final_la_ms": summary["timing_ms"]["linear_algebra"],
                "cpu_throttled_periods_delta": verified["cpu_throttled_periods_delta"],
                "oom_kills_delta": verified["oom_kills_delta"],
            })
        rows.append(row)
    measurements_file = write_gzip_once(HERE / "measurements.json.gz", rows)
    eligible = [row for row in rows if row["status"] == "PASS_INDEPENDENT_REPLAY"
                and row.get("peak_rss_bytes", RSS_CAP_BYTES + 1) <= RSS_CAP_BYTES
                and row.get("cpu_throttled_periods_delta") == 0
                and row.get("oom_kills_delta") == 0]
    complete = len(eligible) == len(ARMS) and analysis["status"] == "COMPLETE"
    if complete != (analysis["balanced_comparison"] is not None):
        raise ValueError("balanced comparison eligibility changed")
    for row in eligible:
        if (verified_by_arm[row["arm"]]["base_hash"] != EXPECTED_BASE_HASH
                or verified_by_arm[row["arm"]]["root_table_entries"] != EXPECTED_ROOT_ENTRIES):
            raise ValueError(f"frozen {row['arm']} base hash changed")
    comparison = analysis["balanced_comparison"]
    selected_profile = analysis["selected_profile"]
    if complete:
        if (selected_profile != comparison["selected_profile"]
                or selected_profile not in (None, *PROFILES)):
            raise ValueError("frozen balanced selection changed")
        if set(comparison["profile_mean_ms"]) != set(PROFILES):
            raise ValueError("frozen profile comparison changed")
    elif selected_profile is not None:
        raise ValueError("incomplete result cannot select a Bloom profile")
    eligible_by_arm = {row["arm"]: row for row in eligible}
    selected = eligible_by_arm.get(f"{selected_profile}_a")
    selected_mean_ms = (sum(eligible_by_arm[f"{selected_profile}_{suffix}"]["ic_online_mean_ms_exploratory"]
                            for suffix in ("a", "b")) / 2
                        if selected is not None else None)
    selection = {
        "schema": "n83.bloom-profile-selection/v1",
        "status": (f"BALANCED_{selected_profile.upper()}_SELECTED" if selected_profile is not None else
                   "BOUNDED_UNKNOWN" if complete else "INCOMPLETE"),
        "selected_profile": selected_profile,
        "selected_arm": selected["arm"] if selected else None,
        "candidate_id": selected["candidate_id"] if selected else None,
        "candidate_manifest_sha256": None, "source_base_hash_blake3": None,
        "mean_online_ms_exploratory": selected_mean_ms,
        "balanced_comparison": comparison,
        "arm_order": analysis["arm_order"],
        "analysis_sha256": analysis_sha256,
        "artifact_sha256": artifact_sha256,
        "source_commit": SOURCE_COMMIT,
        "source_sha256": SOURCE_SHA256,
        "workload_id": workload_id,
    }
    if selected is not None:
        manifest_path = HERE / f"candidate_K3300_root11_{selected_profile}.json"
        manifest = panel.read(manifest_path)
        selection["candidate_manifest_sha256"] = panel.file_sha(manifest_path)
        selection["source_base_hash_blake3"] = manifest["factor_base"]["source_base_hash_blake3"]
    panel.write_once(HERE / "selection.json", selection)
    panel.write_once(HERE / "records_index.json", {
        "schema": "n83.bloom-profile-compressed-records/v1",
        "workload_id": workload_id,
        "measurement_rows": len(rows),
        "workload": workload_file,
        "measurements": measurements_file,
        "analysis_sha256": analysis_sha256,
        "artifact_sha256": artifact_sha256,
    })
    print(json.dumps({"workload_id": workload_id,
                      "candidate_ids": [row["candidate_id"] for row in rows],
                      "statuses": [row["status"] for row in rows],
                      "selected_arm": selection["selected_arm"],
                      "selection_status": selection["status"]}, sort_keys=True))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-zip", required=True, type=Path)
    parser.add_argument("--artifact-sha256", required=True)
    parser.add_argument("--analysis-sha256", required=True)
    parser.add_argument("artifact_results", type=Path)
    args = parser.parse_args()
    for name, value in (("artifact", args.artifact_sha256),
                        ("analysis", args.analysis_sha256)):
        if len(value) != 64 or any(ch not in "0123456789abcdef" for ch in value):
            parser.error(f"{name} SHA-256 must be 64 lowercase hex digits")
    seal(args.artifact_results.resolve(), args.artifact_zip.resolve(),
         args.artifact_sha256, args.analysis_sha256)


if __name__ == "__main__":
    main()

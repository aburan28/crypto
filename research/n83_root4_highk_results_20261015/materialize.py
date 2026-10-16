#!/usr/bin/env python3
"""Seal the balanced-order N83 root4 K3550/K3596 comparison."""

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
SOURCE_COMMIT = "8a37349f3d3cc2bc2960b03167adcf57da83ba77"
SOURCE_SHA256 = "c9824acebda0360416c3c507594bb07c17ca798c32cbf8f84bba95ca2c44d5b1"
LOCK_SHA256 = "ac5cd7a3459be776180afe1eb2975c76b1a88f84bd1e99aa20475927ec260a90"
POINT_SHA256 = "48a243619e270abc096e172fc669bb03494fd8c09d4e39c55935121408d44eea"
FIXTURE_SHA256 = "fa660a434e2618a6fb3a9c251843c6bc4c9cc596c612c244ca2bb2ff8d54262d"
EXPECTED_BASE_HASHES = {
    3596: "9bdadb868c18ccaad3446a62316356e3fe21ecc76512c0c27ccc1fa24d661bc1",
    3550: "1dfb31b3e371e6b846947fa28aa1ceb1ed440fd2c89135d9f8b5eb1505c32548",
}
EXPECTED_ROOT_ENTRIES = {3550: 1045986026, 3596: 1073268803}
RSS_CAP_BYTES = 14 * 1024 ** 3
ARMS = (("k3550_a", 3550, "root4"), ("k3596_a", 3596, "root4"),
        ("k3596_b", 3596, "root4"), ("k3550_b", 3550, "root4"))

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
    generator_path = "research/n83_root4_highk_20261015/generate.py"
    exclusion_path = "research/n83_root4_highk_20261015/excluded_bloom_targets.jsonl"
    for relative, expected in (
            (generator_path, "8bbe3128e25950ec188a6bde00b4651c4543e56f2f11b2b8be4c825adb370d98"),
            (exclusion_path, "ef616c5254a9d883e962024b9f4bd7f794cfe5576c73a93b099ac5be6fe45668")):
        if hashlib.sha256(git_file(relative)).hexdigest() != expected:
            raise ValueError(f"executed fixture source differs: {relative}")
    point_bytes = (REPO / "research/n83_root4_highk_20261015/targets.jsonl").read_bytes()
    answer_bytes = (REPO / "research/n83_root4_highk_20261015/answers.jsonl").read_bytes()
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


def sparse_candidate(k, variant, summary, dump_path, field, curve):
    record = panel.candidate(k, variant, summary, panel.read(dump_path),
                             dump_path, field, curve)
    arm = f"k{k}_{variant}"
    if (summary["sparse_canonical_enabled"] is not True
            or summary["root_bloom_enabled"] is not False
            or summary["root_bloom_bytes"] != 0
            or summary["root_table_indexed4_enabled"] is not True
            or summary["root_table_slot_bytes"] != 4
            or summary["root_state_slot_bytes"] != 11
            or summary["root_table_entries"] != EXPECTED_ROOT_ENTRIES[k]
            or summary["base_hash"] != EXPECTED_BASE_HASHES[k]):
        raise ValueError(f"{arm}: sparse selector, root index or base differs")
    record["implementation"]["algorithm_flags"]["sparse_canonical"] = True
    record["point_decomposition"]["canonical_root_selector"] = "sparse_linear_map_exact"
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
            or source["schema"] != "n83.root4-highk-source/v1"
            or source["arm_order"] != [arm for arm, _, _ in ARMS]
            or hashes["examples/koblitz_orbit_dlp_fast.rs"] != SOURCE_SHA256
            or hashes["Cargo.lock"] != LOCK_SHA256
            or hashes["research/n83_root4_highk_20261015/targets.jsonl"] != POINT_SHA256
            or hashes["research/n83_root4_highk_20261015/answers.jsonl"] != FIXTURE_SHA256
            or source["panel_runner_sha256"] != panel.file_sha(
                REPO / "research/n83_root4_highk_20261015/run.py")
            or analysis["schema"] != "n83.root4-highk-analysis/v1"
            or analysis["hardware_tier"] != "cloud_instance_exploratory"
            or analysis["workload_role"] != "secondary_fresh_4000_point_balanced_order_highk_root4_factor_base_comparison"
            or analysis["target_count"] != 4000
            or analysis["arm_order"] != [arm for arm, _, _ in ARMS]
            or len(analysis["arms"]) != len(ARMS)):
        raise ValueError("source, input or allocation identity mismatch")
    points = frozen_points()
    workload = {
        "schema_version": 1, "curve_id": panel.CURVE_ID,
        "target_count": 4000, "public_targets": points,
        "input_law": {"domain": "n83-root4-highk-fresh-4000-20261015-v1",
                      "scalar_rule": "1+(little_endian(SHA256(domain||u32_le(index))) mod (r-1))",
                      "answer_sidecar_sha256": FIXTURE_SHA256},
        "target_sha256": POINT_SHA256,
        "cold_or_warm": "four_fresh_empty_index_and_logs_processes_in_K3550_K3596_K3596_K3550_order; 4000_public_points_reused_per_arm_after_each_setup",
    }
    workload_id = panel.digest(workload)[:12]
    workload_file = write_gzip_once(HERE / "workload.json.gz",
                                    {"workload_id": workload_id, **workload})
    field, curve = panel.curve_record()
    verified_by_arm = {row["arm"]: row for row in analysis["arms"]}
    if set(verified_by_arm) != {arm for arm, _, _ in ARMS}:
        raise ValueError("analysis arms differ from frozen protocol")
    rows = []
    for arm, k, variant in ARMS:
        verified = verified_by_arm[arm]
        row = {"arm": arm, "candidate_id": None, "workload_id": workload_id,
               "run_id": None, "status": verified["status"],
               "hardware_tier": "cloud_instance_exploratory",
               "workload_role": "secondary_fresh_4000_point_balanced_order_highk_root4_factor_base_comparison",
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
            manifest = sparse_candidate(k, variant, summary, dump_path, field, curve)
            candidate_id = manifest["candidate_id"]
            panel.write_once(HERE / f"candidate_K{k}_{variant}.json", manifest)
            run_number = 1 if arm.endswith("_a") else 2
            row.update({
                "candidate_id": candidate_id,
                "run_id": candidate_id + "W" + workload_id + f"R{run_number}",
                "actual_factor_base_points": verified["actual_factor_base_points"],
                "effective_columns": k, "final_rank": summary["rank"],
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
        k = int(row["arm"][1:5])
        if verified_by_arm[row["arm"]]["base_hash"] != EXPECTED_BASE_HASHES[k]:
            raise ValueError(f"frozen {row['arm']} base hash changed")
    comparison = analysis["balanced_comparison"]
    selected_k = analysis["selected_k"]
    if complete:
        if selected_k != comparison["selected_k"] or selected_k not in (None, 3550, 3596):
            raise ValueError("frozen balanced selection changed")
        if comparison["contrast"] != "mean_k3596_a_b_minus_mean_k3550_a_b":
            raise ValueError("frozen comparison definition changed")
    elif selected_k is not None:
        raise ValueError("incomplete result cannot select a factor base")
    eligible_by_arm = {row["arm"]: row for row in eligible}
    selected = eligible_by_arm.get(f"k{selected_k}_a")
    selected_mean_ms = (sum(eligible_by_arm[f"k{selected_k}_{suffix}"]["ic_online_mean_ms_exploratory"]
                            for suffix in ("a", "b")) / 2
                        if selected is not None else None)
    selection = {
        "schema": "n83.root4-highk-selection/v1",
        "status": (f"BALANCED_K{selected_k}_SELECTED" if selected_k is not None else
                   "BOUNDED_UNKNOWN" if complete else "INCOMPLETE"),
        "selected_k": selected_k,
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
        manifest_path = HERE / f"candidate_K{selected_k}_root4.json"
        manifest = panel.read(manifest_path)
        selection["candidate_manifest_sha256"] = panel.file_sha(manifest_path)
        selection["source_base_hash_blake3"] = manifest["factor_base"]["source_base_hash_blake3"]
    panel.write_once(HERE / "selection.json", selection)
    panel.write_once(HERE / "records_index.json", {
        "schema": "n83.root4-highk-compressed-records/v1",
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

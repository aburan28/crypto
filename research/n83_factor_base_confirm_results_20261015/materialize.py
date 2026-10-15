#!/usr/bin/env python3
"""Seal replayed N83 confirmation receipts as canonical IC1 candidates and runs."""

import argparse
import hashlib
import importlib.util
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
PANEL = REPO / "research/n83_factor_base_100_results_20261014/materialize.py"
SOURCE_COMMIT = "e65f35916dfe4791a60b938b8ab0bdc875d0ae7b"
SOURCE_SHA256 = "5e9fd826d22816b429a82e658c517d9387916aabe7660f708e4c29901cce647e"
LOCK_SHA256 = "ac5cd7a3459be776180afe1eb2975c76b1a88f84bd1e99aa20475927ec260a90"
POINT_SHA256 = "451df8524b2c356dd4a523260620c4dff406975a122762b98ddddb7987b75c2c"
FIXTURE_SHA256 = "8e18c9c97ba9f57e5eb43f1b39b9e9575b5cfab87840c25c428dac150c868101"
EXPECTED_BASE_HASHES = {
    "k3300_root11": "a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b",
    "k3100_root22": "bea9c9cda2826a16c5e7309f2fd5267085190e91978a535eba34d254b00554e6",
}
RSS_CAP_BYTES = 14 * 1024 ** 3
ARMS = ((3300, "root11"), (3100, "root22"))

spec = importlib.util.spec_from_file_location("n83_panel_materializer", PANEL)
panel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(panel)


def git_file(path):
    return subprocess.check_output(["git", "show", SOURCE_COMMIT + ":" + path], cwd=REPO)


def frozen_points():
    point_bytes = git_file("research/n83_factor_base_confirm_20261015/targets.jsonl")
    answer_bytes = git_file("research/n83_factor_base_confirm_20261015/answers.jsonl")
    if (hashlib.sha256(point_bytes).hexdigest() != POINT_SHA256
            or hashlib.sha256(answer_bytes).hexdigest() != FIXTURE_SHA256):
        raise ValueError("committed confirmation points or answer sidecar changed")
    points = [json.loads(line) for line in point_bytes.splitlines() if line]
    answers = [json.loads(line) for line in answer_bytes.splitlines() if line]
    if (len(points) != 4000 or len(answers) != 4000
            or len({tuple(point) for point in points}) != 4000
            or any(answer["fixture_index"] != index or answer["point"] != point
                   for index, (point, answer) in enumerate(zip(points, answers)))):
        raise ValueError("confirmation point count, uniqueness or answer pairing differs")
    return points


def seal(root, artifact_zip, artifact_sha256, analysis_sha256):
    if (panel.file_sha(artifact_zip) != artifact_sha256
            or panel.file_sha(root / "analysis.json") != analysis_sha256):
        raise ValueError("raw artifact or independent analysis digest mismatch")
    analysis = panel.read(root / "analysis.json")
    host = panel.read(root / "host.json")
    source = panel.read(root / "source.json")
    hashes = source["hashes"]
    if (host["git_sha"] != SOURCE_COMMIT
            or source["schema"] != "n83.factor-base-confirm-source/v1"
            or hashes["examples/koblitz_orbit_dlp_fast.rs"] != SOURCE_SHA256
            or hashes["Cargo.lock"] != LOCK_SHA256
            or hashes["research/n83_factor_base_confirm_20261015/targets.jsonl"] != POINT_SHA256
            or analysis["schema"] != "n83.factor-base-confirm-analysis/v1"
            or analysis["hardware_tier"] != "cloud_instance_exploratory"
            or analysis["workload_role"] != "fresh_factor_base_confirmation"
            or analysis["target_count"] != 4000
            or len(analysis["arms"]) != len(ARMS)):
        raise ValueError("source, input or allocation identity mismatch")
    points = frozen_points()
    workload = {
        "schema_version": 1, "curve_id": panel.CURVE_ID,
        "target_count": 4000, "public_targets": points,
        "input_law": {"domain": "n83-factor-base-confirm-4000-20261015-v1",
                      "scalar_rule": "1+(little_endian(SHA256(domain||u32_le(index))) mod (r-1))",
                      "answer_sidecar_sha256": FIXTURE_SHA256},
        "target_sha256": POINT_SHA256,
        "cold_or_warm": "one_fresh_process_per_candidate_with_empty_index_and_logs; 4000_target_fresh_confirmation_after_shared_preparation",
    }
    workload_id = panel.digest(workload)[:12]
    panel.write_once(HERE / "workload.json", {"workload_id": workload_id, **workload})
    field, curve = panel.curve_record()
    verified_by_arm = {row["arm"]: row for row in analysis["arms"]}
    if set(verified_by_arm) != {f"k{k}_{variant}" for k, variant in ARMS}:
        raise ValueError("analysis arms differ from frozen protocol")
    rows = []
    for k, variant in ARMS:
        arm = f"k{k}_{variant}"
        verified = verified_by_arm[arm]
        row = {"arm": arm, "candidate_id": None, "workload_id": workload_id,
               "run_id": None, "status": verified["status"],
               "hardware_tier": "cloud_instance_exploratory",
               "workload_role": "fresh_factor_base_confirmation",
               "confirmation_status": analysis["confirmation_status"],
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
            manifest = panel.candidate(k, variant, summary, panel.read(dump_path),
                                       dump_path, field, curve)
            candidate_id = manifest["candidate_id"]
            panel.write_once(HERE / f"candidate_K{k}_{variant}.json", manifest)
            row.update({
                "candidate_id": candidate_id,
                "run_id": candidate_id + "W" + workload_id + "R1",
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
            })
        rows.append(row)
    panel.write_once(HERE / "measurements.json", rows)
    eligible = [row for row in rows if row["status"] == "PASS_INDEPENDENT_REPLAY"
                and row.get("peak_rss_bytes", RSS_CAP_BYTES + 1) <= RSS_CAP_BYTES]
    selected = (min(eligible, key=lambda row: row["ic_online_mean_ms_exploratory"])
                if eligible else None)
    if analysis["selected_arm_descriptive"] != (selected["arm"] if selected else None):
        raise ValueError("descriptive selection changed")
    complete = (len(eligible) == len(ARMS)
                and analysis["frozen_base_matches"] ==
                {arm: True for arm in EXPECTED_BASE_HASHES})
    for row in eligible:
        if row["arm"] in EXPECTED_BASE_HASHES and row["arm"] in verified_by_arm:
            if verified_by_arm[row["arm"]]["base_hash"] != EXPECTED_BASE_HASHES[row["arm"]]:
                raise ValueError("frozen confirmation base hash changed")
    confirmation = {
        "schema": "n83.factor-base-confirmation/v1",
        "status": "CONFIRMATION_RECORDED" if complete else "INCOMPLETE",
        "selected_arm_descriptive": selected["arm"] if selected else None,
        "candidate_id_descriptive": selected["candidate_id"] if selected else None,
        "candidate_manifest_sha256": None, "source_base_hash_blake3": None,
        "mean_online_ms_exploratory": (selected["ic_online_mean_ms_exploratory"]
                                        if selected else None),
        "confirmation_status": analysis["confirmation_status"],
        "paired_difference": analysis["paired_difference"],
        "frozen_base_matches": analysis["frozen_base_matches"],
        "analysis_sha256": analysis_sha256,
        "artifact_sha256": artifact_sha256,
        "source_commit": SOURCE_COMMIT,
        "source_sha256": SOURCE_SHA256,
        "workload_id": workload_id,
    }
    if complete:
        manifest_path = HERE / f"candidate_{selected['arm'].replace('k', 'K', 1)}.json"
        manifest = panel.read(manifest_path)
        confirmation["candidate_manifest_sha256"] = panel.file_sha(manifest_path)
        confirmation["source_base_hash_blake3"] = manifest["factor_base"]["source_base_hash_blake3"]
    panel.write_once(HERE / "confirmation.json", confirmation)
    print(json.dumps({"workload_id": workload_id,
                      "candidate_ids": [row["candidate_id"] for row in rows],
                      "statuses": [row["status"] for row in rows],
                      "selected_arm_descriptive": confirmation["selected_arm_descriptive"],
                      "confirmation_status": confirmation["confirmation_status"]}, sort_keys=True))


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

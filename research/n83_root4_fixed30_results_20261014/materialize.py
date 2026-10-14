#!/usr/bin/env python3
"""Seal independently replayed N83 fixed30 receipts as IC1 candidates and runs."""

import argparse
import copy
import hashlib
import importlib.util
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
PANEL = REPO / "research/n83_factor_base_100_results_20261014/materialize.py"
SOURCE_COMMIT = "7294b395961b8353c8be70e7873c5151c51c45be"
SOURCE_SHA256 = "68f9b22972ee09511e15ae74a998ba74b82e74eed36e3014bd0460c918b6cdcb"
LOCK_SHA256 = "ac5cd7a3459be776180afe1eb2975c76b1a88f84bd1e99aa20475927ec260a90"
POINT_SHA256 = "35e34151dcc1ad6048939560adf28fc7c1b1017833fe005a656c3f4e2bc07de8"
FIXTURE_SHA256 = "f8cbc90f7f592aae5014ca4f2693a3327a72e935970f22e8f02868e4aeed7f71"
ARMS = ((2940, "root4f30"), (2940, "root22"),
        (3200, "root4f30"), (3200, "root21"),
        (3300, "root4f30"), (3300, "root11"))
RSS_CAP_BYTES = 14 * 1024 ** 3

spec = importlib.util.spec_from_file_location("n83_panel_materializer", PANEL)
panel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(panel)


def frozen_points():
    points = subprocess.check_output(
        ["git", "show", SOURCE_COMMIT + ":research/n83_root4_fixed30_20261014/targets.jsonl"],
        cwd=REPO)
    answers = subprocess.check_output(
        ["git", "show", SOURCE_COMMIT + ":research/n83_root4_fixed30_20261014/answers.jsonl"],
        cwd=REPO)
    if (hashlib.sha256(points).hexdigest() != POINT_SHA256
            or hashlib.sha256(answers).hexdigest() != FIXTURE_SHA256):
        raise ValueError("frozen common points or answer sidecar changed")
    public = [json.loads(line) for line in points.splitlines() if line]
    sidecar = [json.loads(line) for line in answers.splitlines() if line]
    if (len(public) != 1000 or len(sidecar) != 1000
            or len({tuple(point) for point in public}) != 1000
            or any(row["fixture_index"] != index or row["point"] != point
                   for index, (point, row) in enumerate(zip(public, sidecar)))):
        raise ValueError("common panel count, uniqueness or pairing differs")
    return public


def candidate(k, variant, summary, dump, dump_path, field, curve):
    # The established builder checks every base and control layout field.
    # Supply only the new mode's corresponding five-byte layout for that
    # check, then replace all fields that describe its executed algorithm.
    checked = copy.deepcopy(summary)
    old_variant = variant
    if variant == "root4f30":
        old_variant = "root21"
        checked.update({"root_table_layout": "indexed83_5byte_exact",
                        "root_table_slot_bytes": 5,
                        "root_table_indexed4_enabled": False,
                        "root_table_indexed5_enabled": True,
                        "indexed5_fixed30_enabled": True})
    manifest = panel.candidate(k, old_variant, checked, dump, dump_path,
                               field, curve)
    record = {key: value for key, value in manifest.items()
              if key != "candidate_id"}
    pdp = record["point_decomposition"]
    flags = record["implementation"]["algorithm_flags"]
    for section in ("point_decomposition", "relation_collection",
                    "relation_linear_algebra", "target_descent"):
        for key in ("solver_source_sha256", "source_digest"):
            if key in record[section]:
                record[section][key] = SOURCE_SHA256
    record["implementation"]["source_commit"] = SOURCE_COMMIT
    record["implementation"]["source_sha256"]["ic"] = SOURCE_SHA256
    flags["indexed4_fixed30_request"] = variant == "root4f30"
    if variant == "root4f30":
        pdp["variant"] = ("S3_pair_root_index_four_summand_first_verified_witness_"
                          "indexed83_4byte_exact_slot_symmetric_state_quotient_"
                          "root_only21_ordinal_labels_fixed30")
        pdp["root_slot"].update({
            "bytes": 4, "fingerprint_bits": 2,
            "request_rule": "fixed_2pow29_request_2pow30_slots"})
        flags.update({"indexed_root4": True, "indexed_root5": False,
                      "indexed5_fixed30_request": False,
                      "indexed4_fixed31_request": False,
                      "root_only21": True})
    label = (f"IC1N83Ckb1a1fb{166 * k}PDP4rootRCguidedLAgaussTDdirectISO0h"
             + panel.digest(record)[:12])
    return {"candidate_id": label, **record}


def seal(root, artifact_zip, artifact_sha256, analysis_sha256):
    if (panel.file_sha(artifact_zip) != artifact_sha256
            or panel.file_sha(root / "analysis.json") != analysis_sha256):
        raise ValueError("raw artifact or independent analysis digest mismatch")
    analysis = panel.read(root / "analysis.json")
    host = panel.read(root / "host.json")
    source = panel.read(root / "source.json")
    hashes = source["hashes"]
    if (host["git_sha"] != SOURCE_COMMIT
            or source["schema"] != "n83.root4-fixed30-source/v1"
            or hashes["examples/koblitz_orbit_dlp_fast.rs"] != SOURCE_SHA256
            or hashes["Cargo.lock"] != LOCK_SHA256
            or hashes["research/n83_root4_fixed30_20261014/targets.jsonl"] != POINT_SHA256
            or analysis["schema"] != "n83.root4-fixed30-analysis/v1"
            or analysis["hardware_tier"] != "cloud_instance_exploratory"
            or analysis["workload_role"] != "secondary_factor_base_selection"
            or analysis["target_count"] != 1000
            or len(analysis["arms"]) != len(ARMS)):
        raise ValueError("executed source, inputs, allocation or analysis differ")
    public = frozen_points()
    workload = {
        "schema_version": 1, "curve_id": panel.CURVE_ID,
        "target_count": 1000, "public_targets": public,
        "input_law": {"domain": "n83-factor-base-refine-1000-20261014-v1",
                      "scalar_rule": "1+(little_endian(SHA256(domain||u32_le(index))) mod (r-1))",
                      "answer_sidecar_sha256": FIXTURE_SHA256},
        "target_sha256": POINT_SHA256,
        "cold_or_warm": "one_fresh_process_per_candidate_with_empty_index_and_logs; 1000_target_secondary_selection_after_shared_preparation",
    }
    workload_id = panel.digest(workload)[:12]
    panel.write_once(HERE / "workload.json", {"workload_id": workload_id,
                                             **workload})
    field, curve = panel.curve_record()
    verified_by_arm = {row["arm"]: row for row in analysis["arms"]}
    if set(verified_by_arm) != {f"k{k}_{variant}" for k, variant in ARMS}:
        raise ValueError("analysis arms differ from the frozen six-arm plan")
    rows = []
    for k, variant in ARMS:
        arm = f"k{k}_{variant}"
        verified = verified_by_arm[arm]
        row = {"arm": arm, "candidate_id": None, "workload_id": workload_id,
               "run_id": None, "status": verified["status"],
               "hardware_tier": "cloud_instance_exploratory",
               "workload_role": "secondary_factor_base_selection",
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
            if len(targets) != 1000:
                raise ValueError(f"{arm}: expected 1000 target rows")
            dump_path = arm_root / "base_dump.jsonl"
            manifest = candidate(k, variant, summary, panel.read(dump_path),
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
                and all(value is True for value in
                        analysis["same_k_base_key_probe_parity"].values()))
    selection = {
        "schema": "n83.root4-fixed30-selection/v1",
        "status": "PANEL_SELECTED" if complete else "PANEL_INCOMPLETE",
        "selected_arm": selected["arm"] if selected else None,
        "candidate_id": selected["candidate_id"] if selected else None,
        "candidate_manifest_sha256": None, "source_base_hash_blake3": None,
        "mean_online_ms_exploratory": (selected["ic_online_mean_ms_exploratory"]
                                        if selected else None),
        "analysis_sha256": analysis_sha256,
        "artifact_sha256": artifact_sha256,
        "source_commit": SOURCE_COMMIT, "source_sha256": SOURCE_SHA256,
        "workload_id": workload_id,
    }
    if complete:
        manifest_path = HERE / f"candidate_{selected['arm'].replace('k', 'K', 1)}.json"
        manifest = panel.read(manifest_path)
        selection["candidate_manifest_sha256"] = panel.file_sha(manifest_path)
        selection["source_base_hash_blake3"] = manifest["factor_base"]["source_base_hash_blake3"]
    panel.write_once(HERE / "selected_candidate.json", selection)
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

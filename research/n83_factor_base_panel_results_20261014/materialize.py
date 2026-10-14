#!/usr/bin/env python3
"""Bind replayed five-point N83 factor-base panel receipts to IC1 identities."""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
CURVE_ID = "EC1N83Ckb1a1haacb69ee9f37"
SOURCE_COMMIT = "6f070ffd3b65d11549afd8253aeb4d62efc193ce"
SOURCE_SHA256 = "07ce5094f0e561a75bd65ca143f2301830447742268936805c0b920f7a93f6ae"
LOCK_SHA256 = "ac5cd7a3459be776180afe1eb2975c76b1a88f84bd1e99aa20475927ec260a90"
POINT_SHA256 = "b96672b4d75df4b90a337c1e471d54dc147edbf236ef8d0054adb716d2314871"
FIXTURE_SHA256 = "c93edc32a77411443229191e2db4cc3eeb7e99976b7f79a8a753ba103bce9909"
R = 8569786107849059
ARMS = ((3400, "root11"), (3200, "root11"), (3240, "root21"),
        (3114, "root22"), (2940, "root22"), (2543, "root22"),
        (2543, "root11"))
ARTIFACT_SHA256 = "ef8b0221257789d8e45672c6125d102e6f11374947869ee03a8f4bd88c2fc469"
ANALYSIS_SHA256 = "e23c04abc39a2c517860c325ec66ab767f86e11b9ac41aa988a875d1d26b9663"
RSS_CAP_BYTES = 14 * 1024 ** 3


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode("utf-8")


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def file_sha(path):
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            sha.update(block)
    return sha.hexdigest()


def read(path):
    return json.loads(path.read_text())


def write_once(path, value):
    if path.exists():
        if read(path) != value:
            raise ValueError(f"existing sealed record differs: {path}")
        return
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def curve_record():
    field = {"p": 2, "n": 83, "basis": "polynomial",
             "defining_polynomial_int": 9671406556917033397649557,
             "element_encoding": "decimal strings of polynomial coefficient bitmasks"}
    curve = {"G": ["6614244651398184860292015", "8868260410029249391539529"],
             "coefficients": {"a1": 1, "a2": 1, "a3": 0, "a4": 0, "a6": 1},
             "cofactor": 1128547018, "curve_order": 9671406556910881928556062,
             "model": "y^2+x*y=x^3+x^2+1", "r": R,
             "target_group": "prime_order_r_subgroup", "trace": 6151469093347}
    if "EC1N83Ckb1a1h" + digest({"field": field, "curve": curve})[:12] != CURVE_ID:
        raise ValueError("curve identity changed")
    curve["curve_id"] = CURVE_ID
    return field, curve


def candidate(k, variant, summary, dump, dump_path, field, curve):
    points = dump["factor_base_point_coordinates"]
    if (len(points) != 166 * k or summary["factor_base_points"] != len(points)
            or summary["orbit_columns"] != k
            or summary["rank"] != k
            or summary["base_hash"] != dump["base_hash"]
            or summary["root_table_layout"] != "indexed83_5byte_exact"
            or summary["root_table_slot_bytes"] != 5
            or summary["root_state_layout"] != {
                "root11": "root_only83_11byte_states_reconstructed_second_root",
                "root21": "root_only83_21byte_states",
                "root22": "root_only83_22byte_states",
            }[variant]
            or summary["root_state_slot_bytes"] != {
                "root11": 11, "root21": 21, "root22": 22}[variant]
            or summary["symmetric_s3_states_enabled"] is not True
            or summary["indexed5_load75_enabled"] is not (variant == "root22")
            or summary["indexed5_fixed30_enabled"] is not (variant != "root22")
            or summary["root_only11_enabled"] is not (variant == "root11")):
        raise ValueError(f"K{k}: base identity mismatch")
    ordered = sorted(points, key=lambda item: (int(item[0]), int(item[1])))
    record = {
        "schema_version": 1, "field": field, "curve": curve,
        "isogeny": "none",
        "endomorphism": {"endomorphism_order_conductor": None,
                         "frobenius_order_conductor": None,
                         "status": "unproved_not_used", "volcano_levels": None},
        "factor_base": {
            "recipe": "scan polynomial x from one; rational lifts; cofactor projection; first K distinct signed-Frobenius orbits",
            "raw_x_start": 1,
            "subgroup_filter": "cofactor_1128547018_projection; discard_identity",
            "subspace": None, "nominal_dimension": None,
            "actual_usable_points": len(points), "effective_columns": k,
            "quotient_rule": "sign and all 83 Frobenius powers; 166 usable points per column",
            "enumerated_set_sha256": digest(ordered),
            "ordered_base_dump_sha256": file_sha(dump_path),
            "source_base_hash_blake3": summary["base_hash"],
        },
        "point_decomposition": {
            "m": 4, "family": "root",
            "variant": ("S3_pair_root_index_four_summand_first_verified_witness_"
                        "indexed83_5byte_exact_slot_symmetric_state_quotient_"
                        + {"root11": "root_only11_reconstruct_second_root_fixed30",
                           "root21": "root_only21_ordinal_labels_fixed30",
                           "root22": "root_only22_ordinal_labels_load75"}[variant]),
            "summation_chain": "two indexed leaves and two query leaves via S3 root lookup",
            "weil_descent_encoding": "none", "equation_order": "none",
            "monomial_order": "none", "internal_matrix_kernel": "none",
            "index_cache_policy": "reusable_target_independent",
            "query_stop": "scan_finite_root_index_once_until_verified_relation",
            "root_slot": {"bytes": 5, "state_id_bits": 29,
                          "root_selector_bits": 1, "fingerprint_bits": 10,
                          "empty_marker": "zero_slot",
                          "lookup_equality": "full_canonical_root_from_retained_state",
                          "duplicate_policy": "first_inserted_exact_root_wins",
                          "request_rule": (
                              "ceil_four_times_retained_states_over_three_then_round_twice_request_to_power_of_two"
                              if variant == "root22" else "fixed_2pow29_request_2pow30_slots")},
            "solver_source_sha256": SOURCE_SHA256,
        },
        "relation_collection": {
            "query_distribution": "guided_scalar_times_G_minus_first_pivotless_base_rep",
            "duplicate_policy": "echelon_rank_gain",
            "stop_criterion": f"rank_{k}",
            "verification": "curve_group_sum_and_folded_row",
            "source_digest": SOURCE_SHA256,
        },
        "relation_linear_algebra": {
            "modulus": R, "columns": k,
            "orbit_quotient": "sign_and_frobenius",
            "row_rule": "signed_frobenius_coefficients_mod_r",
            "rank_criterion": k, "solver": "gauss",
            "implementation": "incremental_echelon",
            "block_parameters": "none", "preconditioner": "none",
            "source_digest": SOURCE_SHA256,
        },
        "target_descent": {
            "policy": "direct_four_summand_target_relation_after_factor_logs",
            "success": "independent_scalar_replay",
            "source_digest": SOURCE_SHA256,
        },
        "implementation": {
            "source_commit": SOURCE_COMMIT,
            "source_sha256": {"ic": SOURCE_SHA256, "cargo_lock": LOCK_SHA256},
            "algorithm_flags": {"packed_root_table": True,
                                "compact_state28": False,
                                "compact_state26": False,
                                "root_only22": variant == "root22",
                                "root_only21": variant == "root21",
                                "root_only11": variant == "root11",
                                "indexed_root5": True,
                                "symmetric_s3_states": True,
                                "root_bloom": False, "dense_request": False,
                                "load75_request": False,
                                "indexed5_load75_request": variant == "root22",
                                "indexed5_fixed30_request": variant != "root22",
                                "sparse_canonical": False,
                                "rank_threads": 4, "target_threads": 1},
            "release_profile": "cargo_release_locked",
        },
    }
    label = (f"IC1N83Ckb1a1fb{166 * k}PDP4rootRCguidedLAgaussTDdirectISO0h"
             + digest(record)[:12])
    return {"candidate_id": label, **record}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-zip", required=True, type=Path,
                        help="downloaded ZIP for the exact five-point panel run")
    parser.add_argument("artifact_results", type=Path,
                        help="directory containing analysis.json, runs/, host.json, source.json")
    args = parser.parse_args()
    root = args.artifact_results.resolve()
    if (file_sha(args.artifact_zip) != ARTIFACT_SHA256
            or file_sha(root / "analysis.json") != ANALYSIS_SHA256):
        raise ValueError("raw artifact or independently replayed analysis digest mismatch")
    analysis = read(root / "analysis.json")
    host = read(root / "host.json")
    source = read(root / "source.json")
    hashes = source["hashes"]
    if (host["git_sha"] != SOURCE_COMMIT
            or hashes["examples/koblitz_orbit_dlp_fast.rs"] != SOURCE_SHA256
            or hashes["Cargo.lock"] != LOCK_SHA256
            or hashes["research/n83_evaluation_20261012/optimization/fixtures/holdout_5_targets.jsonl"] != POINT_SHA256
            or analysis["schema"] != "n83.factor-base-panel-analysis/v1"
            or analysis["hardware_tier"] != "cloud_instance_exploratory"
            or analysis["workload_role"] != "secondary_factor_base_selection"
            or analysis["target_count"] != 5):
        raise ValueError("source, input or hardware identity mismatch")
    point_bytes = subprocess.check_output(
        ["git", "show", SOURCE_COMMIT + ":research/n83_evaluation_20261012/optimization/fixtures/holdout_5_targets.jsonl"],
        cwd=REPO)
    if hashlib.sha256(point_bytes).hexdigest() != POINT_SHA256:
        raise ValueError("frozen five-point input differs from executed input")
    public_points = [json.loads(line) for line in point_bytes.splitlines() if line]
    receipt_bytes = subprocess.check_output(
        ["git", "show", SOURCE_COMMIT + ":research/n83_evaluation_20261012/optimization/fixtures/generator_receipt.json"],
        cwd=REPO)
    if hashlib.sha256(receipt_bytes).hexdigest() != FIXTURE_SHA256:
        raise ValueError("archived target-generation receipt changed")
    fixture = json.loads(receipt_bytes)
    if (len(public_points) != 5 or fixture["holdout_target_count"] != 5
            or fixture["holdout_combined_sha256"] != POINT_SHA256):
        raise ValueError("archived target panel size or receipt differs")
    workload = {
        "schema_version": 1, "curve_id": CURVE_ID,
        "target_count": 5, "public_targets": public_points,
        "input_law": {"seed": fixture["seed"],
                      "random_generator": fixture["random_generator"],
                      "sampling_law": fixture["sampling_law"],
                      "fixture_split": "holdout_1_through_holdout_5",
                      "generator_receipt_sha256": FIXTURE_SHA256},
        "target_sha256": POINT_SHA256,
        "cold_or_warm": "one_fresh_process_per_candidate_with_empty_index_and_logs; five_target_secondary_selection_after_shared_preparation",
    }
    workload_id = digest(workload)[:12]
    write_once(HERE / "workload.json", {"workload_id": workload_id, **workload})
    field, curve = curve_record()
    rows = []
    for k, variant in ARMS:
        arm = f"k{k}_{variant}"
        verified = next(row for row in analysis["arms"] if row["arm"] == arm)
        row = {"arm": arm, "candidate_id": None, "workload_id": workload_id,
               "run_id": None, "status": verified["status"],
               "hardware_tier": "cloud_instance_exploratory",
               "workload_role": "secondary_factor_base_selection",
               "controlled_speedup": None,
               "source_commit": SOURCE_COMMIT,
               "executed_binary_sha256": source["binary_sha256"],
               "host_receipt_sha256": file_sha(root / "host.json"),
               "source_receipt_sha256": file_sha(root / "source.json"),
               "correctness_receipt_sha256": file_sha(root / "analysis.json")}
        usage_path = root / "runs" / arm / "usage.json"
        if usage_path.exists():
            usage = read(usage_path)
            row.update({"exit_status": usage["exit_status"],
                        "full_process_wall_seconds": usage["wall_seconds"],
                        "peak_rss_bytes": usage["peak_rss_bytes"],
                        "sampled_peak_rss_bytes": usage["sampled_peak_rss_bytes"],
                        "rss_guard_bytes": RSS_CAP_BYTES,
                        "stop_reason": usage["stop_reason"],
                        "major_faults": usage["major_faults"],
                        "involuntary_context_switches":
                            usage["involuntary_context_switches"],
                        "resource_envelope_status": (
                            "OBSERVED_PEAK_ABOVE_14_GIB_GUARD"
                            if usage["peak_rss_bytes"] > RSS_CAP_BYTES else
                            "OBSERVED_PEAK_WITHIN_14_GIB_GUARD")})
        if verified["status"] == "PASS_INDEPENDENT_REPLAY":
            arm_root = root / "runs" / arm
            summary = read(arm_root / "stdout.jsonl")
            targets = [json.loads(line) for line in
                       (arm_root / "targets_result.jsonl").read_text().splitlines() if line]
            if len(targets) != 5:
                raise ValueError(f"{arm}: five target rows expected")
            dump_path = arm_root / "base_dump.jsonl"
            dump = read(dump_path)
            manifest = candidate(k, variant, summary, dump, dump_path, field, curve)
            candidate_id = manifest["candidate_id"]
            write_once(HERE / f"candidate_K{k}_{variant}.json", manifest)
            online_ms = [float(target["target_ms"]) for target in targets]
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
                "ic_online_ms_per_target_exploratory": online_ms,
                "ic_online_mean_ms_exploratory": verified["mean_online_ms_exploratory"],
                "ic_online_phases_ms_per_target": [
                    {key: target[key] for key in
                     ("target_query_ms", "target_pdp_ms",
                      "target_relation_check_ms", "target_descent_ms",
                      "target_recovery_check_ms")}
                    for target in targets],
                "setup_ms": verified["setup_ms"],
                "rank_ms": summary["timing_ms"]["rank_stage"],
                "final_la_ms": summary["timing_ms"]["linear_algebra"],
                "peak_rss_bytes": verified["peak_rss_bytes"],
            })
        rows.append(row)
    write_once(HERE / "measurements.json", rows)
    eligible = [row for row in rows if row["status"] == "PASS_INDEPENDENT_REPLAY"
                and row["peak_rss_bytes"] <= RSS_CAP_BYTES]
    selected = min(eligible, key=lambda row: row["ic_online_mean_ms_exploratory"]) if eligible else None
    if analysis["selected_arm_descriptive"] != (selected["arm"] if selected else None):
        raise ValueError("descriptive selection changed")
    selection = {
        "schema": "n83.factor-base-panel-selection/v1",
        "status": "PANEL_SELECTED" if selected and analysis["k2543_base_key_probe_parity"] is True
                  else "PANEL_INCOMPLETE",
        "selected_arm": selected["arm"] if selected else None,
        "candidate_id": selected["candidate_id"] if selected else None,
        "candidate_manifest_sha256": None,
        "source_base_hash_blake3": None,
        "mean_online_ms_exploratory": (selected["ic_online_mean_ms_exploratory"]
                                        if selected else None),
        "analysis_sha256": ANALYSIS_SHA256,
        "artifact_sha256": ARTIFACT_SHA256,
        "source_commit": SOURCE_COMMIT,
        "source_sha256": SOURCE_SHA256,
        "workload_id": workload_id,
    }
    if selection["status"] == "PANEL_SELECTED":
        manifest_path = HERE / f"candidate_{selected['arm'].replace('k', 'K', 1)}.json"
        manifest = read(manifest_path)
        selection["candidate_manifest_sha256"] = file_sha(manifest_path)
        selection["source_base_hash_blake3"] = manifest["factor_base"]["source_base_hash_blake3"]
    write_once(HERE / "selected_candidate.json", selection)
    print(json.dumps({"workload_id": workload_id,
                      "candidates": [row["candidate_id"] for row in rows],
                      "statuses": [row["status"] for row in rows],
                      "selected_arm": selection["selected_arm"],
                      "selection_status": selection["status"]}, sort_keys=True))


if __name__ == "__main__":
    main()

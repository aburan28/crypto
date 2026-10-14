#!/usr/bin/env python3
"""Bind replayed N83 21-byte fixed-slot receipts to immutable IC1 identities."""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
CURVE_ID = "EC1N83Ckb1a1haacb69ee9f37"
SOURCE_COMMIT = "d8e679f137dd13037cfe2b1471be3cfa76d52eb3"
SOURCE_SHA256 = "accd99313204d0e3327968230ed3c64df49d661ebddb466147a9a8fc844bd1c9"
LOCK_SHA256 = "ac5cd7a3459be776180afe1eb2975c76b1a88f84bd1e99aa20475927ec260a90"
POINT_SHA256 = "21c8dfc7d488909e5c4a42c839fcc01bdefa9821c11ba001626a1d15697240e8"
FIXTURE_SHA256 = "e4d74b8b35709c0e9f85b15af2aec13c289a945cf87b833fc32dfcd9a5f08d70"
R = 8569786107849059
ARMS = ((3240, "root21"), (3220, "root21"), (3190, "root21"),
        (2543, "root21"), (2543, "root22"))
WORKLOAD_ID = "64d976711e0e"
ARTIFACT_SHA256 = "7bca69a792f0f686aab8f3a5b2c3fec0546bb98643807f36056d23ce5fb8da2b"
ANALYSIS_SHA256 = "bd6c1d65343903f584e2336372424aad4e65cb0a0d2b5177572396db9802b4e2"
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
            or summary["root_state_layout"] != (
                "root_only83_21byte_states" if variant == "root21" else
                "root_only83_22byte_states")
            or summary["root_state_slot_bytes"] != (21 if variant == "root21" else 22)
            or summary["symmetric_s3_states_enabled"] is not True
            or summary["indexed5_load75_enabled"] is not (variant == "root22")
            or summary["indexed5_fixed30_enabled"] is not (variant == "root21")):
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
                        + ("root_only21_ordinal_labels_fixed30" if variant == "root21" else
                           "root_only22_ordinal_labels_load75")),
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
                              "fixed_2pow29_request_2pow30_slots" if variant == "root21" else
                              "ceil_four_times_retained_states_over_three_then_round_twice_request_to_power_of_two")},
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
                                "indexed_root5": True,
                                "symmetric_s3_states": True,
                                "root_bloom": False, "dense_request": False,
                                "load75_request": False,
                                "indexed5_load75_request": variant == "root22",
                                "indexed5_fixed30_request": variant == "root21",
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
                        help="downloaded GitHub artifact ZIP for run 38520546294")
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
            or hashes["research/n83_indexed5_20261013/primary_target.jsonl"] != POINT_SHA256
            or analysis["schema"] != "n83.root21-tier1-analysis/v1"
            or analysis["hardware_tier"] != "cloud_instance_exploratory"):
        raise ValueError("source, input or hardware identity mismatch")
    point_bytes = subprocess.check_output(
        ["git", "show", SOURCE_COMMIT + ":research/n83_indexed5_20261013/primary_target.jsonl"],
        cwd=REPO)
    if hashlib.sha256(point_bytes).hexdigest() != POINT_SHA256:
        raise ValueError("frozen public point differs from executed input")
    public_point = json.loads(point_bytes)
    workload = {
        "schema_version": 1, "curve_id": CURVE_ID,
        "target_count": 1, "public_targets": [public_point],
        "input_law": {"generator": "generic_curve_scalar_multiplication",
                      "digest": "blake3",
                      "domain": "n83-dense-root-held-out-20261013-v1",
                      "fixture_index": 0,
                      "scalar_rule": "1+low_u64_le(blake3(domain||u32_le(index)))%(r-1)",
                      "source_sha256": FIXTURE_SHA256},
        "target_sha256": POINT_SHA256,
        "cold_or_warm": "one_fresh_process_with_empty_index_and_logs; one_target_online_after_reusable_preparation",
    }
    workload_id = digest(workload)[:12]
    if workload_id != WORKLOAD_ID:
        raise ValueError("one-target workload identity changed")
    write_once(HERE / "workload.json", {"workload_id": workload_id, **workload})
    field, curve = curve_record()
    rho = analysis["rho"]
    rho_ms = (rho["single_target_online_ms_exploratory"]
              if rho["status"] == "PASS_INDEPENDENT_REPLAY" else None)
    rows = []
    for k, variant in ARMS:
        arm = f"k{k}_{variant}"
        verified = next(row for row in analysis["arms"] if row["arm"] == arm)
        row = {"arm": arm, "candidate_id": None, "workload_id": workload_id,
               "run_id": None, "status": verified["status"],
               "hardware_tier": "cloud_instance_exploratory",
               "controlled_speedup": None,
               "source_commit": SOURCE_COMMIT,
               "executed_binary_sha256": source["binary_sha256"],
               "host_receipt_sha256": file_sha(root / "host.json"),
               "source_receipt_sha256": file_sha(root / "source.json"),
               "correctness_receipt_sha256": file_sha(root / "analysis.json")}
        usage_path = root / "runs" / "primary" / arm / "usage.json"
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
            arm_root = root / "runs" / "primary" / arm
            summary = read(arm_root / "stdout.jsonl")
            target = read(arm_root / "targets_result.jsonl")
            dump_path = arm_root / "base_dump.jsonl"
            dump = read(dump_path)
            manifest = candidate(k, variant, summary, dump, dump_path, field, curve)
            candidate_id = manifest["candidate_id"]
            write_once(HERE / f"candidate_K{k}_{variant}.json", manifest)
            online_ms = float(target["target_ms"])
            row.update({
                "candidate_id": candidate_id,
                "run_id": candidate_id + "W" + workload_id + "R1",
                "actual_factor_base_points": verified["actual_factor_base_points"],
                "effective_columns": k, "final_rank": summary["rank"],
                "root_table_allocated_bytes": verified["root_table_allocated_bytes"],
                "root_state_allocated_bytes": verified["root_state_allocated_bytes"],
                "regular_states": verified["regular_states"],
                "root_table_entries": verified["root_table_entries"],
                "target_probes": verified["target_probes"],
                "ic_online_ms_exploratory": online_ms,
                "ic_online_phases_ms": {key: target[key] for key in
                    ("target_query_ms", "target_pdp_ms",
                     "target_relation_check_ms", "target_descent_ms",
                     "target_recovery_check_ms")},
                "rho_online_ms_exploratory": rho_ms,
                "rho_over_ic_online_exploratory": (rho_ms / online_ms
                    if rho_ms is not None else None),
                "setup_ms": verified["setup_ms"],
                "rank_ms": verified["rank_ms"],
                "final_la_ms": verified["final_la_ms"],
                "peak_rss_bytes": verified["peak_rss_bytes"],
            })
        rows.append(row)
    write_once(HERE / "measurements.json", rows)
    print(json.dumps({"workload_id": workload_id,
                      "candidates": [row["candidate_id"] for row in rows],
                      "statuses": [row["status"] for row in rows]}, sort_keys=True))


if __name__ == "__main__":
    main()

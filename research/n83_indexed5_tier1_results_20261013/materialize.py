#!/usr/bin/env python3
"""Bind frozen, independently replayed N83 indexed5 receipts to IC1 identities."""

import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
CURVE_ID = "EC1N83Ckb1a1haacb69ee9f37"
SOURCE_COMMIT = "8954fe0571bc5bb7c152c7f12d8ff387ae250ee0"
SOURCE_SHA256 = "64171575e58b97c4ad9041b4b49e5f2b5697d26f71be73861fc36ea879db9161"
LOCK_SHA256 = "ac5cd7a3459be776180afe1eb2975c76b1a88f84bd1e99aa20475927ec260a90"
POINT_SHA256 = "21c8dfc7d488909e5c4a42c839fcc01bdefa9821c11ba001626a1d15697240e8"
FIXTURE_SHA256 = "e4d74b8b35709c0e9f85b15af2aec13c289a945cf87b833fc32dfcd9a5f08d70"
R = 8569786107849059
ARMS = (1798, 1500, 1271)


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


def candidate(k, summary, dump, dump_path, field, curve):
    points = dump["factor_base_point_coordinates"]
    if (len(points) != 166 * k or summary["factor_base_points"] != len(points)
            or summary["orbit_columns"] != k
            or summary["rank"] != k
            or summary["base_hash"] != dump["base_hash"]
            or summary["root_table_layout"] != "indexed83_5byte_exact"
            or summary["root_table_slot_bytes"] != 5
            or summary["root_state_layout"] != "packed83_11byte_roots"
            or summary["root_state_slot_bytes"] != 28):
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
            "variant": "S3_pair_root_index_four_summand_first_verified_witness_indexed83_5byte_exact_slot",
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
                          "request_rule": "twofold_power_of_two_over_two_keys_per_state"},
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
                                "compact_state28": True, "indexed_root5": True,
                                "root_bloom": False, "dense_request": False,
                                "load75_request": False,
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
    parser.add_argument("artifact_results", type=Path,
                        help="directory containing analysis.json, runs/, host.json, source.json")
    args = parser.parse_args()
    root = args.artifact_results.resolve()
    analysis = read(root / "analysis.json")
    host = read(root / "host.json")
    source = read(root / "source.json")
    hashes = source["hashes"]
    if (host["git_sha"] != SOURCE_COMMIT
            or hashes["examples/koblitz_orbit_dlp_fast.rs"] != SOURCE_SHA256
            or hashes["Cargo.lock"] != LOCK_SHA256
            or hashes["research/n83_indexed5_20261013/primary_target.jsonl"] != POINT_SHA256
            or analysis["hardware_tier"] != "cloud_instance_exploratory"):
        raise ValueError("source, input or hardware identity mismatch")
    point_path = REPO / "research/n83_indexed5_20261013/primary_target.jsonl"
    if file_sha(point_path) != POINT_SHA256:
        raise ValueError("local public point differs from executed input")
    public_point = read(point_path)
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
    write_once(HERE / "workload.json", {"workload_id": workload_id, **workload})
    field, curve = curve_record()
    rho = analysis["rho"]
    rho_ms = (rho["single_target_online_ms_exploratory"]
              if rho["status"] == "PASS_INDEPENDENT_REPLAY" else None)
    rows = []
    for k in ARMS:
        arm = f"k{k}_indexed5"
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
                        "stop_reason": usage["stop_reason"],
                        "major_faults": usage["major_faults"],
                        "involuntary_context_switches":
                            usage["involuntary_context_switches"]})
        if verified["status"] == "PASS_INDEPENDENT_REPLAY":
            arm_root = root / "runs" / "primary" / arm
            summary = read(arm_root / "stdout.jsonl")
            target = read(arm_root / "targets_result.jsonl")
            dump_path = arm_root / "base_dump.jsonl"
            dump = read(dump_path)
            manifest = candidate(k, summary, dump, dump_path, field, curve)
            candidate_id = manifest["candidate_id"]
            write_once(HERE / f"candidate_K{k}.json", manifest)
            online_ms = float(target["target_ms"])
            row.update({
                "candidate_id": candidate_id,
                "run_id": candidate_id + "W" + workload_id + "R1",
                "actual_factor_base_points": verified["actual_factor_base_points"],
                "effective_columns": k, "final_rank": summary["rank"],
                "root_table_allocated_bytes": verified["root_table_allocated_bytes"],
                "root_state_allocated_bytes": verified["root_state_allocated_bytes"],
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

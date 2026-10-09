"""Freeze EC1/IC1/W/R identities for the prepared-S3 local panel.

The candidate and workload hash inputs contain no paths, timestamps, run
seeds, timings, or JSON floats. Run records retain those measurements.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CURVES = {
    41: {"r": 549756390943, "cofactor": 4, "order": 2199025563772,
         "trace": -2308219, "modulus_powers": [41, 3, 0],
         "generator": [2056947637384, 1635505394702]},
    53: {"r": 21044858204113, "cofactor": 428, "order": 9007199311360364,
         "trace": -56619371, "modulus_powers": [53, 6, 2, 1, 0],
         "generator": [198217578752339, 7929897206038174]},
}
SOURCES = {
    "candidate": ROOT / "experiments/koblitz-s3-pair-query-20261007-v4/candidate.rs",
    "prepared": ROOT / "examples/s3_prepared_targetspan.rs",
}


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()


def digest(value: object) -> str:
    return hashlib.sha256(canonical(value)).hexdigest()[:12]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        assert json.loads(path.read_text()) == value, f"refusing to replace {path}"
        return
    path.write_bytes(canonical(value) + b"\n")


def curve_record(n: int) -> dict:
    cfg = CURVES[n]
    field = {"characteristic_p": 2, "degree_n": n,
             "representation": "polynomial_basis",
             "irreducible_polynomial_powers": cfg["modulus_powers"],
             "element_encoding": "unsigned_integer_polynomial_coefficients_lsb_is_x0"}
    curve = {"model": "y2_plus_xy_equals_x3_plus_a2x2_plus_b",
             "coefficients": {"a1": 1, "a2": 0, "a3": 0, "a4": 0, "a6": 1},
             "group_order": cfg["order"], "trace_over_declared_field": cfg["trace"],
             "subgroup_order_r": cfg["r"], "cofactor": cfg["cofactor"],
             "generator_G": cfg["generator"],
             "point_encoding": "affine_x_y_as_unsigned_polynomial_basis_integers",
             "target_group": "cyclic_prime_order_subgroup_generated_by_G"}
    curve_id = f"EC1N{n}Ckb1h{digest({'field': field, 'curve': curve})}"
    curve["curve_id"] = curve_id
    return {"field": field, "curve": curve}


def candidate_record(n: int, arm: str, report: dict, curve: dict,
                     snapshot: str, source_tree: str, cargo_lock_sha: str) -> dict:
    cfg = CURVES[n]
    source_sha = sha(SOURCES[arm])
    prepared = arm == "prepared"
    return {
        "field": curve["field"], "curve": curve["curve"],
        "isogeny": "none",
        "endomorphism": {"order_conductor": None,
                         "frobenius_order_conductor": None,
                         "volcano_levels": None},
        "factor_base": {
            "construction": "natural_low_x_scan_then_subgroup_filter_then_signed_frobenius_orbits",
            "start_x": 1, "last_scanned_x": report["factor_base_scanned_x"],
            "enumerated_set_digest_algorithm": "blake3",
            "enumerated_set_digest": report["factor_base_digest"],
            "nominal_subspace_dimension": None,
            "requested_orbit_columns": 244,
            "geometric_point_count_before_subgroup_filter": None,
            "actual_usable_points_B": report["factor_base_points"],
            "quotient_rule": "signed_frobenius_orbit_labels",
            "effective_relation_columns": report["orbit_columns"],
            "source_sha256": source_sha,
        },
        "point_decomposition": {
            "summand_count_m": 4,
            "summation_polynomial": "S3_root_index_for_four_point_decompositions",
            "weil_descent_encoding": "normal_basis_s3_root_table",
            "equation_order": "left_right_relative_index_state_order",
            "solver_family": "compact_s3_root_index",
            "implementation": "parallel_sharded_paired_s3_inversion",
            "implementation_source_sha256": source_sha,
            "monomial_order": "none",
            "internal_matrix_kernel": "none",
            "root_table_shards": 64,
            "root_table_slots": 16777216,
            "orientation_policy": "one_target_seeded_pseudorandom_global_frobenius_shift_per_state",
            "limits": {"max_relation_attempts": "max(8*columns,columns+32)"},
            "cache_policy": "per_process_index_and_relations",
        },
        "relation_collection": {
            "query_distribution": "splitmix64(seed)_mod_(r-1)_plus_1_known_scalar_subgroup_points",
            "seed_policy": "run_supplied_u64",
            "batch_width": 14,
            "relation_verification": "four_group_points_sum_to_known_scalar_query",
            "duplicate_policy": "unique_relation_scalars; dependent_rows_retained_as_witnesses",
            "stop_criterion": ("full_factor_base_rank_before_target_input"
                               if prepared else "target_coefficient_row_in_span_after_batch"),
            "attempt_limit": "max(8*columns,columns+32)",
            "source_sha256": source_sha,
        },
        "relation_linear_algebra": {
            "modulus": cfg["r"],
            "matrix_row": "sum_signed_frobenius_labels_equals_known_scalar",
            "orbit_quotient": "signed_frobenius",
            "rank_criterion": "full_244" if prepared else "target_row_span_or_full_244_fallback",
            "solver": "incremental_modular_gaussian_elimination",
            "block_parameters": "none", "preconditioner": "none",
            "source_sha256": source_sha,
        },
        "target_descent": {
            "policy": "one_direct_four_point_s3_decomposition_then_linear_recovery",
            "recursive_solvers": "none",
            "success_rule": "recovered_scalar_times_G_equals_public_target",
            "stop_rule": "one_verified_target_or_attempt_limit",
            "source_sha256": source_sha,
        },
        "implementation": {
            "frozen_source_commit": snapshot,
            "rust_library_source_tree_git_oid": source_tree,
            "example_source_sha256": source_sha,
            "cargo_lock_sha256": cargo_lock_sha,
            "algorithm_flags": {"curve_a": 0, "columns": 244,
                                "relation_batch_width": 14,
                                "paired_s3_inversion": True},
        },
    }


def workload_record(n: int, fixture: dict, curve_id: str) -> dict:
    return {"curve_id": curve_id, "subgroup_order_r": CURVES[n]["r"],
            "target_count": 1, "target_point": fixture["public_point"],
            "target_label": fixture["label"],
            "target_input_law": "known_scalar_sha256_public_point_fixture_v1",
            "target_seed_label": "s3-pair-pilot-20261007-v1",
            "target_index": int(fixture["label"][1:]),
            "scalar_derivation_sha256": fixture["scalar_derivation_sha256"],
            "relation_seed": 20260928,
            "preparation_state_at_online_start": "reusable_factor_base_index_and_relations_ready"}


def main(ledger: Path, output: Path) -> None:
    entries = [json.loads(line) for line in ledger.read_text().splitlines()]
    assert len(entries) == 48
    snapshot = "833ef86191947ba8e4def56a78ee6c51a4ea0ff2"
    source_tree = subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", f"{snapshot}:src"], text=True
    ).strip()
    current_tree = subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True
    ).strip()
    assert source_tree == current_tree, "Rust library source differs from the frozen snapshot"
    cargo_lock_sha = sha(ROOT / "Cargo.lock")
    curves = {n: curve_record(n) for n in (41, 53)}
    candidates = {}
    workloads = {}
    for n in (41, 53):
        curve_id = curves[n]["curve"]["curve_id"]
        write(output / "curves" / f"{curve_id}.json", curves[n])
        for arm in ("candidate", "prepared"):
            report = next(row["report"] for row in entries if row["n"] == n and row["arm"] == arm)
            record = candidate_record(n, arm, report, curves[n], snapshot, source_tree, cargo_lock_sha)
            name = (f"IC1N{n}Ckb1fb{report['factor_base_points']}PDP4rootRCsample"
                    f"LAgaussTDdirectISO0h{digest(record)}")
            candidates[n, arm] = name
            write(output / "candidates" / f"{name}.json", {**record, "candidate_id": name})
        fixtures = json.loads((ROOT / f"experiments/koblitz-s3-fidelity-20261007/pilot/n{n}/fixtures.json").read_text())["fixtures"]
        for fixture in fixtures:
            record = workload_record(n, fixture, curve_id)
            workload_id = digest(record)
            workloads[n, fixture["label"]] = workload_id
            write(output / "workloads" / f"W{workload_id}.json", {**record, "workload_id": workload_id})
    run_rows = []
    for order, row in enumerate(entries, 1):
        report = row["report"]
        candidate_id = candidates[row["n"], row["arm"]]
        workload_id = workloads[row["n"], row["label"]]
        run_rows.append({"candidate_id": candidate_id, "workload_id": workload_id,
                         "run_id": f"{candidate_id}W{workload_id}R1",
                         "run_order": order, "arm": row["arm"], "n": row["n"],
                         "target_label": row["label"], "target": report["target"],
                         "status": "verified_exploratory_mac_contended" if row["verified_native_and_fixture"] else "failure",
                         "strict_host_isolation_receipt": None,
                         "verified_native_and_fixture": row["verified_native_and_fixture"],
                         "binary_sha256": row["binary_sha256"],
                         "target_sha256": row["target_sha256"],
                         "raw_sha256": row["raw_sha256"],
                         "resource_envelope": {"relation_workers": 14, "timeout_seconds": 180},
                         "rank_attempts": report["rank_attempts"],
                         "rank_failures": report["rank_failures"],
                         "verified_relations": report["rank_verified_relations"],
                         "novel_rows": report["rank_new_rows"],
                         "final_rank": report["rank"],
                         "actual_factor_base_points": report["factor_base_points"],
                         "folded_columns": report["orbit_columns"],
                         "target_relation_found": report["target_relation_found"],
                         "recovered_scalar": report["recovered_scalar"],
                         "peak_rss_bytes": report["peak_rss_bytes"],
                         "timing_ms": report["timing_ms"]})
    measurement_file = output / "measurements.jsonl"
    measurement_data = b"".join(canonical(row) + b"\n" for row in run_rows)
    if measurement_file.exists():
        assert measurement_file.read_bytes() == measurement_data
    else:
        measurement_file.write_bytes(measurement_data)
    print(json.dumps({"curves": [curves[n]["curve"]["curve_id"] for n in (41, 53)],
                      "candidates": list(candidates.values()),
                      "workload_count": len(workloads), "run_count": len(run_rows),
                      "measurements_sha256": sha(measurement_file)}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    main(args.ledger, args.output)

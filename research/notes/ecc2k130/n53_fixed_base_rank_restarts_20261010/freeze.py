#!/usr/bin/env python3
"""Freeze exact candidate and workload identities before n53 pilot runs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
INPUTS = HERE / "inputs"
CURVE_ID = "EC1N53Ce0hb097de99be9a"
BASE_HASH = "7af2460c8b5a2c29f9d1aa7fefecbcde3a6ce761293d0dab0cecfc3bc980b973"
PARENT_COMMIT = "1ac61222583ce43d5d986c925807e6d43070639c"
SEEDS = (530053, 530054, 530055)
CAPS = (None, 200_000, 400_000)


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()


def digest(value: object) -> str:
    return hashlib.sha256(canonical(value)).hexdigest()


def file_digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_frozen(path: Path, value: object) -> None:
    encoded = json.dumps(value, sort_keys=True, indent=2, ensure_ascii=False) + "\n"
    if path.exists():
        assert path.read_text() == encoded, f"frozen file changed: {path}"
    else:
        path.write_text(encoded)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--binary", type=Path,
        default=ROOT / "target/release/examples/koblitz_orbit_dlp_fast_online",
    )
    args = parser.parse_args()
    base_path = INPUTS / "base_n53_k220.jsonl"
    target_path = INPUTS / "development_q.jsonl"
    source_path = ROOT / "examples/koblitz_orbit_dlp_fast_online.rs"
    binary_path = args.binary
    base = json.loads(base_path.read_text())
    target = json.loads(target_path.read_text())
    assert base["base_hash"] == BASE_HASH
    assert (base["n"], base["a"], base["factor_base_points"], base["orbit_columns"]) == (
        53, 0, 23_320, 220
    )
    assert target == [7_960_849_849_661_793, 7_443_722_527_872_608]
    assert len(base["factor_base_point_coordinates"]) == 23_320
    assert len(base["factor_base_point_labels"]) == 23_320
    assert len(base["factor_base_representatives"]) == 220

    registry = json.loads((ROOT / "docs/curves/registry.json").read_text())
    matches = [
        (curve, rep)
        for curve in registry["curves"]
        for rep in curve.get("representations", [])
        if rep.get("ec1") == CURVE_ID
    ]
    assert len(matches) == 1
    registry_curve, rep = matches[0]
    assert digest({"field": rep["field"], "curve": rep["curve"]})[:12] == CURVE_ID[-12:]
    assert int(registry_curve["order"]) == 9_007_199_311_360_364
    assert int(registry_curve["trace"]) == -56_619_371
    assert int(rep["curve"]["subgroup_order"]) == base["subgroup_order"]
    assert rep["curve"]["cofactor"] == 428

    source_sha = file_digest(source_path)
    lock_sha = file_digest(ROOT / "Cargo.lock")
    base_sha = file_digest(base_path)
    target_sha = file_digest(target_path)
    protocol_sha = file_digest(HERE / "PROTOCOL.md")
    fixture_sha = file_digest(
        HERE / "inputs/n53_holdout_parent_fixture.jsonl"
    )
    candidate_ids: dict[str, str] = {}
    for cap in CAPS:
        label = "control" if cap is None else f"cap{cap}"
        record = {
            "field": rep["field"],
            "curve": {
                **rep["curve"],
                "curve_id": CURVE_ID,
                "order": int(registry_curve["order"]),
                "trace": int(registry_curve["trace"]),
            },
            "isogeny": "none",
            "endomorphism": {
                "order_conductor": "1",
                "order_status": "proved",
                "frobenius_order_conductor": "68476319",
                "frobenius_status": "proved",
                "volcano_levels": {"2": "not_applicable", "263": 0},
                "order_certificate_sha256": file_digest(ROOT / "docs/curves/traits.json"),
            },
            "factor_base": {
                "construction": "ascending nonzero x scan, lift, cofactor projection, signed Frobenius closure",
                "base_hash_blake3": BASE_HASH,
                "enumerated_set_sha256": base_sha,
                "nominal_dimension": None,
                "scan_stop_policy": "first 220 distinct usable signed Frobenius orbits",
                "geometric_point_count": 23_320,
                "B": 23_320,
                "orbit_quotient_rule": "signed Frobenius orbit",
                "effective_columns": 220,
                "subgroup_filter": "cofactor 428 projection and identity exclusion",
            },
            "point_decomposition": {
                "summands": 4,
                "solver_family": "root",
                "solver": "normal-basis S3 regular-root index plus exact point lift",
                "equation_order": "fixed left/right/relative-shift/root order",
                "weil_descent_encoding": "GF(2^53) normal-basis rotations with polynomial-basis point lift",
                "monomial_order": "none",
                "internal_matrix_kernel": "none",
                "support_probe_cap": cap,
                "cache_policy": "index constructed once per cold process",
                "source_sha256": source_sha,
            },
            "relation_collection": {
                "query_distribution": "LCG scalar stream modulo r-1 plus one",
                "first_rank_seed": "workload record",
                "pivot_rule": "first pivotless orbit column",
                "query_point": "[scalar]G - representative[pivotless column]",
                "restart_after_support_probes": cap,
                "partial_witness_reuse": False,
                "verification": "exact point lift and rank-row replay",
                "duplicate_policy": "record dependent rows; continue until rank 220",
                "stop_criterion": "full rank 220",
                "source_sha256": source_sha,
            },
            "relation_linear_algebra": {
                "modulus": str(base["subgroup_order"]),
                "matrix_row": "verified guided orbit-log relation",
                "orbit_quotient": "signed Frobenius coefficients",
                "rank_criterion": 220,
                "solver": "gauss",
                "block_parameters": "none",
                "preconditioner": "none",
                "source_sha256": source_sha,
            },
            "target_descent": {
                "method": "direct",
                "recursive_solvers": "none",
                "stop_rule": "first verified four-point target relation",
                "source_sha256": source_sha,
            },
            "implementation": {
                "example_source_sha256": source_sha,
                "cargo_lock_sha256": lock_sha,
                "library_parent_commit": PARENT_COMMIT,
                "flags": {
                    "KIC_RANK_PROBE_CAP": cap,
                    "KIC_DUMP_BASE": True,
                    "KIC_DUMP_RANK": True,
                    "rank_workers": 1,
                },
            },
        }
        candidate_id = (
            "IC1N53Ce0fb23320PDP4rootRCguidedLAgaussTDdirectISO0h"
            + digest(record)[:12]
        )
        candidate_ids[label] = candidate_id
        write_frozen(
            HERE / f"candidate_{label}.json",
            {"candidate_id": candidate_id, "record": record},
        )

    workload_ids: dict[str, str] = {}
    for seed in SEEDS:
        record = {
            "curve_id": CURVE_ID,
            "subgroup_order": str(base["subgroup_order"]),
            "target": target,
            "target_generation": "public_hash_to_curve_cofactor",
            "public_hash_seed": 53_261_207,
            "parent_fixture_sha256": fixture_sha,
            "input_law": "one previously published public point Q, scalar hidden from IC",
            "rank_seed": seed,
            "target_count": 1,
            "cache_state": "cold",
        }
        workload_id = "W" + digest(record)[:12]
        workload_ids[str(seed)] = workload_id
        write_frozen(
            HERE / f"workload_{seed}.json",
            {"workload_id": workload_id, "record": record},
        )

    write_frozen(
        HERE / "FROZEN.json",
        {
            "schema": "n53-fixed-base-rank-restart-freeze-v1",
            "candidate_ids": candidate_ids,
            "workload_ids": workload_ids,
            "source_sha256": source_sha,
            "cargo_lock_sha256": lock_sha,
            "binary_sha256": file_digest(binary_path),
            "base_file_sha256": base_sha,
            "base_hash_blake3": BASE_HASH,
            "expected_raw_x_scanned": 464,
            "development_q_sha256": target_sha,
            "parent_fixture_sha256": fixture_sha,
            "protocol_sha256": protocol_sha,
            "source_parent_commit": PARENT_COMMIT,
            "pilot_rank_seeds": list(SEEDS),
            "pilot_probe_caps": [None, 200_000, 400_000],
            "resource_envelope": {
                "threads": 1,
                "wall_cap_seconds": 60,
                "observed_rss_cap_bytes": 16 * 1024**3,
            },
        },
    )
    print(json.dumps({"candidates": candidate_ids, "workloads": workload_ids}, sort_keys=True))


if __name__ == "__main__":
    main()

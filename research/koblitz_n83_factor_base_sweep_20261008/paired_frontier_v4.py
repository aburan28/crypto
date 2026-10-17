#!/usr/bin/env python3
"""Lazy, source-gated cross-stream extension of the frozen N83 v3 grid."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

import paired_frontier_v3 as v3


ROOT = Path(__file__).resolve().parent
PARENT = ROOT / "paired-frontier-v3.json"
PILOT = ROOT / "verification/cross-stream-pilot-import-source-20261013/census.json"
OUTPUT = ROOT / "paired-frontier-v4.json"
PARENT_SHA256 = "f7b860c9a652154541c4e6cf7d52ad13bd5a7436a41d38855d206cde36d4b0ac"
PILOT_SHA256 = "d1a490e2e73552d96a5ee327e34367c5bb32fbc6f6d591eb9c0e5f029118a285"
STREAMS_PER_ARM = 27
OTHER_STREAMS = STREAMS_PER_ARM - 1
SIZE_COUNT = len(v3.EXPECTED_SIZES)
PAIRS_PER_ARM = STREAMS_PER_ARM * OTHER_STREAMS * SIZE_COUNT * SIZE_COUNT
PAIR_COUNT = 2 * PAIRS_PER_ARM
RESIDUAL_COUNT = len(v3.RESIDUAL_SUMMANDS)


def streams(base_specs: list[dict]) -> list[dict]:
    groups: dict[tuple, dict[str, int]] = defaultdict(dict)
    for index, base in enumerate(base_specs):
        if base["parameter_unit"] != "orbit_columns" or base["family"] not in v3.v2lib.POLICIES:
            continue
        key = (base["a"], base["family"], base["seed"], base["closure"])
        size = str(base["parameter"])
        v3.require(size not in groups[key], "duplicate public-x base size")
        groups[key][size] = index
    v3.require(len(groups) == 54, "public-x stream coverage changed")
    result = []
    for (arm, policy, seed, closure), indices in sorted(groups.items()):
        v3.require(set(indices) == {str(x) for x in v3.EXPECTED_SIZES},
                   "public-x size ladder changed")
        result.append({"curve_a": arm, "policy": policy, "seed": seed,
                       "closure": closure, "base_indices": indices})
    v3.require([row["curve_a"] for row in result] == [0] * 27 + [1] * 27,
               "arm stream partition changed")
    return result


def pair(design: dict, index: int) -> dict:
    v3.require(0 <= index < PAIR_COUNT, "cross-stream pair index outside grid")
    arm, local = divmod(index, PAIRS_PER_ARM)
    source_stream, local = divmod(local, OTHER_STREAMS * SIZE_COUNT * SIZE_COUNT)
    other_slot, local = divmod(local, SIZE_COUNT * SIZE_COUNT)
    pool_stream = other_slot + (other_slot >= source_stream)
    source_size, pool_size = divmod(local, SIZE_COUNT)
    base = design["streams"][arm * STREAMS_PER_ARM + source_stream]
    pool = design["streams"][arm * STREAMS_PER_ARM + pool_stream]
    base_k = v3.EXPECTED_SIZES[source_size]
    pool_k = v3.EXPECTED_SIZES[pool_size]
    return {
        "curve_a": arm,
        "base": {"policy": base["policy"], "seed": base["seed"],
                 "closure": base["closure"], "orbits": base_k,
                 "parent_base_index": base["base_indices"][str(base_k)]},
        "pool": {"policy": pool["policy"], "seed": pool["seed"],
                 "closure": pool["closure"], "orbits": pool_k,
                 "parent_base_index": pool["base_indices"][str(pool_k)]},
    }


def source_class(spec: dict, pilot_keys: set[tuple]) -> str:
    base, pool = spec["base"], spec["pool"]
    if base["closure"] != pool["closure"]:
        return "cross_closure_domain_unimplemented"
    if base["closure"] != "signed_frobenius":
        return "non_signed_closure_domain_unimplemented"
    if spec["curve_a"] != 0:
        return "primary_curve_only"
    if ((base["policy"], base["seed"], base["orbits"]) in pilot_keys
            and (pool["policy"], pool["seed"], pool["orbits"]) in pilot_keys):
        return "pilot_exact_import_gate_available_unexecuted"
    return "exact_archive_pair_unavailable"


def disposition(spec: dict, solver: dict, residual: int,
                pilot_keys: set[tuple]) -> str:
    """One source disposition per tuple; no empirical success is inferred."""
    if residual == 0:
        return "no_residual_pair_axis_redundant"
    if solver["max_large_primes"] != residual:
        return "large_prime_count_differs_from_residual_shape"
    if solver["summands"] not in (5, 6) or solver["summands"] - residual < 3:
        return "mixed_orbit_sum_shape_unavailable"
    if solver["backend"] != "sat_cdcl":
        return f"paired_backend_adapter_unavailable:{solver['backend']}"
    if solver["symmetry"] != v3.EXACT_PAIR_SYMMETRY:
        return "exact_first_pair_symmetry_adapter_unavailable"
    if solver["split_bits"] != 0:
        return "paired_split_adapter_unavailable"
    if solver["enumeration"] != "binary":
        return "paired_gray_enumeration_adapter_unavailable"
    if solver["threads"] != 1:
        return "paired_thread_setting_unmeasured"
    if spec["curve_a"] == 0 and solver["linear_algebra"] != "modular_biguint":
        return "paired_full_width_linear_algebra_required"
    return source_class(spec, pilot_keys)


def pilot_keys(receipt: dict) -> set[tuple]:
    v3.require(receipt["schema"] == "n83.cross-stream-pilot-source-census/v1"
               and receipt["ordered_cross_stream_pairs"] == 648
               and receipt["primary_rows"] == 27
               and receipt["archive_import_executed"] is False,
               "pilot source census changed")
    keys = {(row["policy"], row["seed"], row["columns"])
            for row in receipt["primary_objects"]}
    v3.require(len(keys) == 27, "pilot object identities changed")
    return keys


def case(design: dict, parent: dict, ordinal: int) -> dict:
    v3.require(0 <= ordinal < design["total_tuple_count"],
               "case ordinal outside v4 grid")
    if ordinal < design["v3_tuple_count_preserved"]:
        return v3.case(parent, json.loads(v3.PARENT.read_text()), ordinal)
    v3.require(v3.sha256(PILOT) == PILOT_SHA256,
               "pilot source census bytes changed")
    per_pair = RESIDUAL_COUNT * design["solver_tuple_count"]
    pair_index, inside = divmod(ordinal - design["v3_tuple_count_preserved"], per_pair)
    residual_index, solver_index = divmod(inside, design["solver_tuple_count"])
    spec = pair(design, pair_index)
    solver = v3.solver_from_index(solver_index, parent)
    residual = v3.RESIDUAL_SUMMANDS[residual_index]
    keys = pilot_keys(json.loads(PILOT.read_text()))
    return {
        "schema": "n83.factor-base-cross-stream-case/v4",
        "ordinal": ordinal, "pair_index": pair_index,
        "residual_index": residual_index, "solver_index": solver_index,
        "pair": spec, "residual_summands": residual,
        "base_summands": solver["summands"] - residual,
        "solver": solver,
        "conditional_unordered_point_multisets": None,
        "uniform_nonidentity_hit_ceiling": None,
        "exact_overlap_measured": False,
        "disposition": disposition(spec, solver, residual, keys),
        "executed_complete_pipeline": False,
        "selected_best_total_runtime": None,
    }


def build() -> dict:
    v3.require(v3.sha256(PARENT) == PARENT_SHA256, "v3 parent bytes changed")
    v3.require(v3.sha256(PILOT) == PILOT_SHA256, "pilot source census bytes changed")
    parent = json.loads(PARENT.read_text())
    v3.require(parent["schema"] == "n83.factor-base-design/v3-paired-frontier"
               and parent["total_tuple_count"] == 518918400
               and parent["solver_tuple_count"] == 51840,
               "v3 parent dimensions changed")
    base = json.loads(v3.PARENT.read_text())
    stream_rows = streams(base["base_specs"])
    keys = pilot_keys(json.loads(PILOT.read_text()))
    categories: Counter = Counter()
    representatives = {}
    for index in range(PAIR_COUNT):
        spec = pair({"streams": stream_rows}, index)
        category = source_class(spec, keys)
        categories[category] += 1
        representatives.setdefault(category, spec)
    v3.require(categories == {
        "cross_closure_domain_unimplemented": 117612,
        "non_signed_closure_domain_unimplemented": 34848,
        "primary_curve_only": 8712,
        "exact_archive_pair_unavailable": 8064,
        "pilot_exact_import_gate_available_unexecuted": 648,
    }, "cross-stream source classification changed")
    cross_tuples = PAIR_COUNT * RESIDUAL_COUNT * parent["solver_tuple_count"]
    dispositions: Counter = Counter()
    for category, multiplicity in categories.items():
        spec = representatives[category]
        for solver_index in range(parent["solver_tuple_count"]):
            solver = v3.solver_from_index(solver_index, parent)
            for residual in v3.RESIDUAL_SUMMANDS:
                dispositions[disposition(spec, solver, residual, keys)] += multiplicity
    v3.require(sum(dispositions.values()) == cross_tuples,
               "cross-stream disposition count mismatch")
    v3.require(dispositions["pilot_exact_import_gate_available_unexecuted"] == 2592,
               "pilot source-available solver count changed")
    return {
        "schema": "n83.factor-base-design/v4-cross-stream-frontier",
        "study": parent["study"],
        "parent_v3_sha256": PARENT_SHA256,
        "pilot_source_census_sha256": PILOT_SHA256,
        "v3_tuple_count_preserved": parent["total_tuple_count"],
        "streams": stream_rows,
        "stream_count": len(stream_rows),
        "cross_stream_pair_count": PAIR_COUNT,
        "source_class_pair_counts": dict(sorted(categories.items())),
        "disposition_counts": dict(sorted(dispositions.items())),
        "residual_summands_axis": list(v3.RESIDUAL_SUMMANDS),
        "solver_tuple_count": parent["solver_tuple_count"],
        "cross_stream_tuple_count": cross_tuples,
        "total_tuple_count": parent["total_tuple_count"] + cross_tuples,
        "coverage": "v3 ordinal prefix preserved; every ordered distinct-stream, same-arm public-x base/pool pair across all 11 sizes and all three closures, crossed with residual 0/1/2 and every v3 solver setting",
        "measurement_scope": "source availability only; cross-stream overlap and conditional support unmeasured",
        "selected_best_total_runtime": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    group.add_argument("--case", type=int)
    args = parser.parse_args()
    generated = v3.canonical_bytes(build())
    if args.case is not None:
        v3.require(OUTPUT.read_bytes() == generated, "saved v4 design differs")
        print(json.dumps(case(json.loads(OUTPUT.read_text()),
                              json.loads(PARENT.read_text()), args.case),
                         sort_keys=True, indent=2))
        return
    if args.check:
        v3.require(OUTPUT.read_bytes() == generated, "saved v4 design differs")
    elif OUTPUT.exists():
        v3.require(OUTPUT.read_bytes() == generated, "refusing changed v4 output")
    else:
        OUTPUT.write_bytes(generated)
    print(json.dumps({"status": "PASS", "sha256": hashlib.sha256(generated).hexdigest()}))


if __name__ == "__main__":
    main()

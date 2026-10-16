#!/usr/bin/env python3
"""Ordinal-addressable paired-base extension of the frozen N83 sweep.

This catalogs conditional recipes and source availability. It constructs no
factor base, runs no solver, and makes no total-runtime selection.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
from itertools import combinations, product
import json
import math
from pathlib import Path

import size_frontier_v2 as v2lib


ROOT = Path(__file__).resolve().parent
PARENT = ROOT / "size-frontier-v2.json"
SUPPORT = ROOT / "pilot-01/support-moments.json"
IMPORT = ROOT / "verification/k600-mixed-orbit-worker-source-20261012/import-receipts.json"
SEARCH_AUDIT = ROOT / "verification/mixed-orbit-search-15000-20261012/audit.json"
OUTPUT = ROOT / "paired-frontier-v3.json"
PARENT_SHA256 = "7526b6218dfe0b074368901f3869ddf4a978cb0bacffe8654056056700f8a972"
SUPPORT_SHA256 = "38abcfd36691b9565f97b6d86d4e89ce71ce0b18675586c65b40f945320b8c3f"
IMPORT_SHA256 = "9d1d8336ba3941ee03b37d67e4da8a447f2abe0c2613090f098d83ed4b936de0"
SEARCH_AUDIT_SHA256 = "32420a62bb6ccf8f9df7f03dea9040e518606740968110dc8c0ac47cb50cafd5"
RESIDUAL_SUMMANDS = (0, 1, 2)
EXACT_PAIR_SYMMETRY = "ordered_first_pair"
POINTS_PER_SIGNED_ORBIT = 166
EXPECTED_SIZES = (32, 64, 128, 256, 600, 900, 1182, 2048, 4096, 8192, 16627)
VERIFIED_PAIRS = {(600, 2048), (1182, 2048)}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_bytes(value: dict) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def paired_specs(parent: dict) -> list[dict]:
    groups: dict[tuple, dict[int, int]] = defaultdict(dict)
    for index, base in enumerate(parent["base_specs"]):
        if (base["parameter_unit"] != "orbit_columns"
                or base["family"] not in v2lib.POLICIES):
            continue
        group = (base["a"], base["family"], base["seed"], base["closure"])
        size = base["parameter"]
        require(size not in groups[group], "duplicate public-x base size")
        groups[group][size] = index
    require(len(groups) == 54, "public-x family/arm/seed/closure coverage changed")
    result = []
    for group, indices in sorted(groups.items()):
        sizes = tuple(sorted(indices))
        require(sizes == EXPECTED_SIZES, "public-x size ladder changed")
        arm, policy, seed, closure = group
        for small, envelope in combinations(sizes, 2):
            result.append({
                "curve_a": arm, "policy": policy, "seed": seed,
                "closure": closure, "small_orbits": small,
                "envelope_orbits": envelope,
                "small_base_index": indices[small],
                "envelope_base_index": indices[envelope],
            })
    require(len(result) == 2970, "paired specification count changed")
    return result


def pair_class(pair: dict) -> str:
    if pair["closure"] != "signed_frobenius":
        return "closure_incompatible"
    if pair["curve_a"] != 0:
        return "primary_curve_only"
    if (pair["policy"] != "public_x_hash" or pair["seed"] != 2026100801
            or (pair["small_orbits"], pair["envelope_orbits"]) not in VERIFIED_PAIRS):
        return "exact_nested_archives_unverified"
    return "source_available"


def disposition(pair: dict, solver: dict, residual: int) -> str:
    """First applicable source-availability reason; never a runtime result."""
    if residual == 0:
        return "no_residual_pair_axis_redundant"
    if solver["max_large_primes"] != residual:
        return "large_prime_count_differs_from_residual_shape"
    if solver["summands"] not in (5, 6) or solver["summands"] - residual < 3:
        return "mixed_orbit_sum_shape_unavailable"
    if solver["backend"] != "sat_cdcl":
        return f"paired_backend_adapter_unavailable:{solver['backend']}"
    if solver["symmetry"] != EXACT_PAIR_SYMMETRY:
        return "exact_first_pair_symmetry_adapter_unavailable"
    if solver["split_bits"] != 0:
        return "paired_split_adapter_unavailable"
    if solver["enumeration"] != "binary":
        return "paired_gray_enumeration_adapter_unavailable"
    if solver["threads"] != 1:
        return "paired_thread_setting_unmeasured"
    if pair["curve_a"] == 0 and solver["linear_algebra"] != "modular_biguint":
        return "paired_full_width_linear_algebra_required"
    return pair_class(pair) if pair_class(pair) != "source_available" else (
        "bounded_mixed_orbit_partial_worker_available_full_pipeline_unmeasured"
    )


def solver_axes(parent: dict) -> dict:
    axes = {key: list(values) for key, values in parent["solver_axes"].items()}
    require(EXACT_PAIR_SYMMETRY not in axes["symmetry"], "exact symmetry already in parent")
    axes["symmetry"].append(EXACT_PAIR_SYMMETRY)
    return axes


def solver_from_index(index: int, design: dict) -> dict:
    digits = v2lib.solver_digits(index, design["solver_radices"])
    return {key: design["solver_axes"][key][digit]
            for key, digit in zip(v2lib.AXIS_KEYS, digits)}


def support_count(pair: dict, summands: int, residual: int) -> int:
    small_summands = summands - residual
    if small_summands < 0:
        return 0
    small_points = POINTS_PER_SIGNED_ORBIT * pair["small_orbits"]
    residual_points = POINTS_PER_SIGNED_ORBIT * (
        pair["envelope_orbits"] - pair["small_orbits"])
    return (math.comb(small_points + small_summands - 1, small_summands)
            * math.comb(residual_points + residual - 1, residual))


def case(design: dict, parent: dict, ordinal: int) -> dict:
    require(0 <= ordinal < design["total_tuple_count"], "case ordinal outside v3 grid")
    require(design["support_receipt_sha256"] == sha256(SUPPORT),
            "support receipt changed")
    prior = design["v2_tuple_count_preserved"]
    if ordinal < prior:
        return v2lib.case(parent, ordinal)
    per_pair = len(RESIDUAL_SUMMANDS) * design["solver_tuple_count"]
    pair_index, within_pair = divmod(ordinal - prior, per_pair)
    residual_index, solver_index = divmod(within_pair, design["solver_tuple_count"])
    pair = design["pair_specs"][pair_index]
    solver = solver_from_index(solver_index, design)
    residual = RESIDUAL_SUMMANDS[residual_index]
    orders = {int(row["curve_a"]): int(row["subgroup_order"])
              for row in json.loads(SUPPORT.read_text())["cases"]}
    signed_closure = pair["closure"] == "signed_frobenius"
    combinations_count = (support_count(pair, solver["summands"], residual)
                          if signed_closure else None)
    order_minus_one = orders[pair["curve_a"]] - 1
    return {
        "schema": "n83.factor-base-paired-case/v3",
        "ordinal": ordinal, "pair_index": pair_index,
        "residual_index": residual_index, "solver_index": solver_index,
        "pair": pair, "residual_summands": residual,
        "small_summands": solver["summands"] - residual,
        "solver": solver,
        "conditional_unordered_point_multisets": (
            str(combinations_count) if signed_closure else None),
        "uniform_nonidentity_hit_ceiling": ({
            "numerator": str(min(combinations_count, order_minus_one)),
            "denominator": str(order_minus_one),
        } if signed_closure else None),
        "disposition": disposition(pair, solver, residual),
        "executed_complete_pipeline": False,
        "selected_best_total_runtime": None,
    }


def count_dispositions(pairs: list[dict], axes: dict) -> dict[str, int]:
    groups = Counter(pair_class(pair) for pair in pairs)
    counts: Counter = Counter()
    representative = {pair_class(pair): pair for pair in pairs}
    for category, multiplicity in groups.items():
        pair = representative[category]
        for residual in RESIDUAL_SUMMANDS:
            for choices in product(*(axes[key] for key in v2lib.AXIS_KEYS)):
                solver = dict(zip(v2lib.AXIS_KEYS, choices))
                counts[disposition(pair, solver, residual)] += multiplicity
    return dict(sorted(counts.items()))


def build() -> dict:
    require(sha256(PARENT) == PARENT_SHA256, "v2 parent bytes changed")
    require(sha256(SUPPORT) == SUPPORT_SHA256, "support receipt changed")
    require(sha256(IMPORT) == IMPORT_SHA256, "K600 import receipt changed")
    require(sha256(SEARCH_AUDIT) == SEARCH_AUDIT_SHA256, "prior search audit changed")
    parent = json.loads(PARENT.read_text())
    require(parent["schema"] == "n83.factor-base-design/v2-size-frontier",
            "v2 parent schema changed")
    require(parent["total_tuple_count"] == 57024000
            and parent["base_spec_count"] == 1320
            and parent["solver_tuple_count"] == 43200,
            "v2 parent dimensions changed")
    pairs = paired_specs(parent)
    axes = solver_axes(parent)
    radices = [len(axes[key]) for key in v2lib.AXIS_KEYS]
    solvers = math.prod(radices)
    require(solvers == 51840, "v3 solver-axis product changed")
    paired_tuples = len(pairs) * len(RESIDUAL_SUMMANDS) * solvers
    counts = count_dispositions(pairs, axes)
    require(sum(counts.values()) == paired_tuples, "paired disposition count mismatch")
    require(counts.get("bounded_mixed_orbit_partial_worker_available_full_pipeline_unmeasured") == 8,
            "source-available pair count changed")
    return {
        "schema": "n83.factor-base-design/v3-paired-frontier",
        "study": parent["study"],
        "parent_v2_sha256": PARENT_SHA256,
        "support_receipt_sha256": SUPPORT_SHA256,
        "k600_import_receipt_sha256": IMPORT_SHA256,
        "prior_k1182_search_audit_sha256": SEARCH_AUDIT_SHA256,
        "v2_tuple_count_preserved": parent["total_tuple_count"],
        "pair_specs": pairs, "pair_spec_count": len(pairs),
        "residual_summands_axis": list(RESIDUAL_SUMMANDS),
        "solver_axes": axes, "solver_radices": radices,
        "solver_tuple_count": solvers,
        "paired_tuple_count": paired_tuples,
        "total_tuple_count": parent["total_tuple_count"] + paired_tuples,
        "disposition_counts": counts,
        "coverage": "v2 ordinal prefix preserved; every same-arm/policy/seed/closure ordered public-x size pair crossed with 0/1/2 residuals and every v3 solver setting",
        "measurement_scope": "source availability and conditional point-multiset support only",
        "selected_best_total_runtime": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    group.add_argument("--case", type=int)
    args = parser.parse_args()
    if args.case is not None:
        require(OUTPUT.read_bytes() == canonical_bytes(build()),
                "saved v3 design differs")
        design = json.loads(OUTPUT.read_text())
        print(json.dumps(case(design, json.loads(PARENT.read_text()), args.case),
                         sort_keys=True, indent=2))
        return
    generated = canonical_bytes(build())
    if args.check:
        require(OUTPUT.read_bytes() == generated, "saved v3 design differs")
        print(json.dumps({"status": "PASS", "sha256": hashlib.sha256(generated).hexdigest()}))
        return
    if OUTPUT.exists():
        require(OUTPUT.read_bytes() == generated, "refusing changed v3 output")
    else:
        OUTPUT.write_bytes(generated)
    print(json.dumps({"status": "PASS", "sha256": hashlib.sha256(generated).hexdigest()}))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact, addressable size extension of the frozen N83 factor-base design.

This is a support and workload-size screen. It does not construct new bases,
measure a solver, or select a runtime winner.
"""

import argparse
import hashlib
import json
import math
from collections import Counter
from decimal import Decimal, localcontext
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parent
V1 = ROOT / "design.json"
SUPPORT = ROOT / "pilot-01" / "support-moments.json"
OUTPUT = ROOT / "size-frontier-v2.json"
NEW_SIZES = (1182, 2048, 4096, 8192, 16627)
POLICIES = ("public_x_sequential", "public_x_hash", "public_x_gray_prefix")
SEEDS = (2026100801, 2026100802, 2026100803)
CLOSURES = ("none", "sign", "signed_frobenius")
POINTS_PER_SIGNED_ORBIT = 166


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def canonical_bytes(value: dict) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def new_base_specs() -> list[dict]:
    return [
        {
            "a": arm,
            "family": policy,
            "parameter": columns,
            "parameter_unit": "orbit_columns",
            "seed": seed,
            "closure": closure,
        }
        for arm in (0, 1)
        for policy in POLICIES
        for columns in NEW_SIZES
        for seed in SEEDS
        for closure in CLOSURES
    ]


def support_count(points: int, summands: int) -> int:
    return math.comb(points + summands - 1, summands)


def threshold(order: int, summands: int) -> int:
    low, high = 0, 1
    while support_count(POINTS_PER_SIGNED_ORBIT * high, summands) < order - 1:
        high *= 2
    while high - low > 1:
        middle = (high + low) // 2
        if support_count(POINTS_PER_SIGNED_ORBIT * middle, summands) >= order - 1:
            high = middle
        else:
            low = middle
    return high


def decimal_ratio(numerator: int, denominator: int) -> str:
    with localcontext() as context:
        context.prec = 24
        return f"{Decimal(numerator) / Decimal(denominator):.6E}"


def solver_digits(ordinal: int, radices: list[int]) -> list[int]:
    count = math.prod(radices)
    require(0 <= ordinal < count, "solver ordinal outside grid")
    digits = [0] * len(radices)
    for index in reversed(range(len(radices))):
        digits[index] = ordinal % radices[index]
        ordinal //= radices[index]
    return digits


AXIS_KEYS = (
    "summands", "split_bits", "symmetry", "backend", "enumeration",
    "max_large_primes", "threads", "linear_algebra",
)


def disposition(base: dict, solver: dict) -> str:
    """One source-availability disposition, in priority order, per tuple."""
    backend = solver["backend"]
    if base["parameter_unit"] == "binary_dimension" and base["parameter"] > 20:
        return "resource_gate_dimension_above_20"
    if base["family"] == "frobenius_divisor":
        return "dimension_one_torsion_only"
    if solver["symmetry"] == "signed_frobenius" and base["closure"] != "signed_frobenius":
        return "incompatible_symmetry_and_domain"
    if (backend == "compact_s3_four_sum"
            and (solver["summands"] != 4
                 or not base["family"].startswith("public_x")
                 or base["closure"] != "signed_frobenius")):
        return "incompatible_compact_four_sum_domain"
    if solver["split_bits"] != 0 and backend not in ("f4", "wdsat"):
        return "split_adapter_not_implemented_for_backend"
    if solver["max_large_primes"] != 0:
        return "wide_large_prime_partial_producer_required"
    if backend.startswith("fes_"):
        return "wide_anf_filter_and_full_verifier_required"
    if base["a"] == 0 and solver["linear_algebra"] != "modular_biguint":
        return "subgroup_order_exceeds_u64"
    return "unexecuted_recipe_requires_n83_pipeline_adapter"


def disposition_class(base: dict) -> tuple:
    """Properties on which the ordered disposition rules depend."""
    return (
        base["a"],
        base["parameter_unit"] == "binary_dimension" and base["parameter"] > 20,
        base["family"] == "frobenius_divisor",
        base["family"].startswith("public_x"),
        base["closure"],
    )


def count_dispositions(base_specs: list[dict], axes: dict) -> dict[str, int]:
    representatives: dict[tuple, dict] = {}
    multiplicities: Counter = Counter()
    for base in base_specs:
        key = disposition_class(base)
        representatives[key] = base
        multiplicities[key] += 1
    counts: Counter = Counter()
    for key, multiplicity in multiplicities.items():
        base = representatives[key]
        for choices in product(*(axes[axis] for axis in AXIS_KEYS)):
            solver = dict(zip(AXIS_KEYS, choices))
            counts[disposition(base, solver)] += multiplicity
    return dict(sorted(counts.items()))


def case(design: dict, ordinal: int) -> dict:
    solver_count = design["solver_tuple_count"]
    require(0 <= ordinal < design["total_tuple_count"], "case ordinal outside grid")
    base_index, solver_index = divmod(ordinal, solver_count)
    digits = solver_digits(solver_index, design["solver_radices"])
    axes = design["solver_axes"]
    solver = {key: axes[key][digit] for key, digit in zip(AXIS_KEYS, digits)}
    base = design["base_specs"][base_index]
    return {
        "schema": "n83.factor-base-case/v2-size-frontier",
        "ordinal": ordinal,
        "base_index": base_index,
        "solver_index": solver_index,
        "base": base,
        "solver": solver,
        "disposition": disposition(base, solver),
        "executed": False,
        "runtime_rank_eligible": False,
    }


def build() -> dict:
    v1_bytes = V1.read_bytes()
    support_bytes = SUPPORT.read_bytes()
    v1 = json.loads(v1_bytes)
    support = json.loads(support_bytes)
    require(v1["schema"] == "n83.factor-base-design/v1", "frozen design schema")
    require(v1["base_spec_count"] == len(v1["base_specs"]) == 1050, "frozen base count")
    require(v1["solver_tuple_count"] == math.prod(v1["solver_radices"]) == 43200,
            "frozen solver count")
    require(v1["total_tuple_count"] == 45360000, "frozen tuple count")
    require(support["schema"] == "n83.support-moment-screen/v1", "support receipt schema")
    require(support["manifest_sha256"] == hashlib.sha256(
        (ROOT / "pilot-01" / "manifest.json").read_bytes()).hexdigest(),
        "support-to-panel binding")
    original = {
        (base["a"], base["family"], base["parameter"], base["seed"], base["closure"])
        for base in v1["base_specs"]
        if base["family"] in POLICIES
    }
    expected = {
        (arm, policy, columns, seed, closure)
        for arm in (0, 1)
        for policy in POLICIES
        for columns in (32, 64, 128, 256, 600, 900)
        for seed in SEEDS
        for closure in CLOSURES
    }
    require(original == expected, "frozen public-x grid differs from expected Cartesian set")
    additions = new_base_specs()
    require(len(additions) == 270, "new base count")
    require(not any(base in v1["base_specs"] for base in additions), "repeated base spec")

    orders = {}
    for row in support["cases"]:
        arm, order = row["curve_a"], int(row["subgroup_order"])
        require(arm not in orders or orders[arm] == order, "subgroup-order conflict")
        orders[arm] = order
    require(set(orders) == {0, 1}, "missing curve arm")
    thresholds = [
        {"curve_a": arm, "summands": m,
         "first_columns_with_vacuous_markov_ceiling": threshold(orders[arm], m)}
        for arm in (0, 1) for m in range(2, 7)
    ]
    prior_thresholds = {
        (row["curve_a"], row["summands"]): row["first_columns_with_vacuous_markov_ceiling"]
        for row in support["thresholds"]
    }
    require(all(prior_thresholds[row["curve_a"], row["summands"]]
                == row["first_columns_with_vacuous_markov_ceiling"] for row in thresholds),
            "threshold disagrees with frozen support receipt")

    sizes = sorted(set(NEW_SIZES) | {600, 900})
    rows = []
    for arm in (0, 1):
        order = orders[arm]
        for columns in sizes:
            points = POINTS_PER_SIGNED_ORBIT * columns
            for m in range(2, 7):
                multisets = support_count(points, m)
                capped = min(multisets, order - 1)
                # If each query yields at most one full-smooth row, this is a
                # necessary query count for the Markov rank ceiling to reach 1/2.
                half_rank_queries = max(
                    columns,
                    (columns * (order - 1) + 2 * capped - 1) // (2 * capped),
                )
                rows.append({
                    "curve_a": arm,
                    "orbit_columns": columns,
                    "conditional_signed_point_count": points,
                    "summands": m,
                    "unordered_multisets_with_repetition": str(multisets),
                    "uniform_nonidentity_hit_ceiling": {
                        "numerator": str(capped), "denominator": str(order - 1),
                        "decimal": decimal_ratio(capped, order - 1),
                    },
                    "minimum_queries_for_half_rank_ceiling": half_rank_queries,
                    "full_pair_multisets_if_materialized": points * (points + 1) // 2,
                })

    all_specs = v1["base_specs"] + additions
    disposition_counts = count_dispositions(all_specs, v1["solver_axes"])
    result = {
        "schema": "n83.factor-base-design/v2-size-frontier",
        "study": v1["study"],
        "parent_design_sha256": hashlib.sha256(v1_bytes).hexdigest(),
        "support_receipt_sha256": hashlib.sha256(support_bytes).hexdigest(),
        "v1_tuple_count_preserved": v1["total_tuple_count"],
        "new_sizes": list(NEW_SIZES),
        "new_base_spec_count": len(additions),
        "base_specs": all_specs,
        "base_spec_count": len(all_specs),
        "solver_axes": v1["solver_axes"],
        "solver_radices": v1["solver_radices"],
        "solver_tuple_count": v1["solver_tuple_count"],
        "total_tuple_count": len(all_specs) * v1["solver_tuple_count"],
        "disposition_counts": disposition_counts,
        "disposition_policy": "source-availability audit; wide large-prime row adapter exists but a partial producer is required",
        "subgroup_orders": {str(arm): str(order) for arm, order in sorted(orders.items())},
        "signed_orbit_point_count_conditional_on_successful_construction":
            POINTS_PER_SIGNED_ORBIT,
        "support_thresholds": thresholds,
        "size_cases": rows,
        "coverage": "all v1 ordinals preserved; new public-x size specs and every v1 solver-axis combination appended and addressable",
        "execution": "design and exact support counts only; new bases and solver tuples unexecuted",
        "s3_destination_for_future_new_objects":
            "s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/v2-size-frontier/",
        "selected_best_total_runtime": None,
    }
    require(result["total_tuple_count"] == 57024000, "v2 tuple count")
    require(sum(disposition_counts.values()) == result["total_tuple_count"],
            "disposition partition")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=int, help="print one addressable v2 tuple")
    args = parser.parse_args()
    design = build()
    encoded = canonical_bytes(design)
    if OUTPUT.exists():
        require(OUTPUT.read_bytes() == encoded, "immutable v2 design changed")
    else:
        OUTPUT.write_bytes(encoded)
    if args.case is not None:
        print(json.dumps(case(design, args.case), indent=2, sort_keys=True))
    else:
        print(json.dumps({"base_specs": design["base_spec_count"],
                          "tuples": design["total_tuple_count"]}, sort_keys=True))


if __name__ == "__main__":
    main()

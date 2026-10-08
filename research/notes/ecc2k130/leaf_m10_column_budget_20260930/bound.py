#!/usr/bin/env python3
"""Exact column-budget capacity bound for the measured disjoint m10 slots."""
from __future__ import annotations

from decimal import Decimal, getcontext
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
INPUTS = ROOT / "research/notes/ecc2k130/m10_export_capacity_20260925/inputs"
ANALYSIS = ROOT / "research/notes/ecc2k130/leaf_m10_support_outcome_20260930/ANALYSIS.json"
M = 10


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def max_product(columns: int) -> int:
    """Maximize ∏(1 + 2*c_i) over nonnegative integer c_i with sum≤columns."""
    assert columns >= 0
    low, high_slots = divmod(columns, M)
    return (1 + 2 * low) ** (M - high_slots) * (3 + 2 * low) ** high_slots


def minimum_columns_for_q(q: int) -> int:
    lo, hi = 0, 1
    while max_product(hi) < q:
        hi *= 2
    while lo < hi:
        mid = (lo + hi) // 2
        if max_product(mid) >= q:
            hi = mid
        else:
            lo = mid + 1
    assert max_product(lo - 1) < q <= max_product(lo)
    return lo


def decimal_ratio(numerator: int, denominator: int) -> str:
    return format(Decimal(numerator) / Decimal(denominator), ".15g")


def main() -> None:
    getcontext().prec = 80
    balanced_path = INPUTS / "balanced_m10_result.json"
    unequal_path = INPUTS / "unequal_m10_result.json"
    balanced = json.loads(balanced_path.read_text())
    unequal = json.loads(unequal_path.read_text())
    analysis = json.loads(ANALYSIS.read_text())
    assert analysis["status"] == "PASS_DIAGNOSTIC"
    assert balanced["q"] == unequal["q"]
    q = balanced["q"]
    assert balanced["nonzero_signed_columns"] == 3988
    assert unequal["normalized_global_nonzero_signed_columns"] == 8062
    assert all(all(arm["projected_columns_disjoint_across_slots"]
                   for arm in row["arms"].values()) for row in analysis["rows"])
    assert {row["curve"] for row in analysis["rows"]} == {
        "source", "leaf_1_0", "leaf_1_4"}

    q_threshold = minimum_columns_for_q(q)
    budgets = {"balanced": balanced["nonzero_signed_columns"],
               "unequal": unequal["normalized_global_nonzero_signed_columns"]}
    cases = {}
    for arm, budget in budgets.items():
        upper = max_product(budget)
        assert upper < q
        leaf_rows = {}
        for row in analysis["rows"]:
            if row["curve"] == "source":
                continue
            observed = row["arms"][arm]
            leaf_rows[row["curve"]] = {
                "full_domain_uncompressed_columns": observed[
                    "uncompressed_projected_sign_union"],
                "full_domain_physical_tuple_product": observed[
                    "physical_tuple_product"],
                "full_domain_product_over_q": observed[
                    "necessary_tuple_capacity_over_q"],
            }
            assert observed["uncompressed_projected_sign_union"] > q_threshold
        cases[arm] = {
            "source_normalized_nonzero_log_columns": budget,
            "max_physical_tuple_product_at_equal_columns": str(upper),
            "necessary_single_uniform_target_support_ceiling": decimal_ratio(upper, q),
            "max_product_slot_column_allocation": [budget // M +
                (i < budget % M) for i in range(M)],
            "minimum_disjoint_slot_columns_for_product_at_least_q": q_threshold,
            "minimum_columns_factor_vs_source": decimal_ratio(q_threshold, budget),
            "observed_full_leaf_domains": leaf_rows,
        }
    output = {
        "schema": "ecc2k130-leaf-m10-column-budget-bound-v1",
        "status": "PASS_CONDITIONAL_BOUND",
        "scope": "subsets of the archived collision-free, disjoint fixed-leaf m10 slots; no cross-slot log identifications",
        "field_degree": 131,
        "summands": M,
        "subgroup_order": str(q),
        "analysis_sha256": sha(ANALYSIS),
        "source_screen_sha256": {"balanced": sha(balanced_path),
                                 "unequal": sha(unequal_path)},
        "cases": cases,
        "PDP_yield": None,
        "full_ECDLP_cost": None,
        "method_crossover": None,
    }
    (HERE / "BOUND.json").write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": output["status"], "q_threshold": q_threshold,
                      "ceilings": {arm: row[
                          "necessary_single_uniform_target_support_ceiling"]
                          for arm, row in cases.items()}}, sort_keys=True))


if __name__ == "__main__":
    main()

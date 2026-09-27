#!/usr/bin/env python3
"""Exact-count transfer checks for the current compact S5 orbit extractor.

Integer counts and rational coverage/work bounds are exact architecture-specific
projections, not measured n=131 costs. Probability and log2 display fields
are floating-point approximations. The four-sum counting bound is universal:
at most F**4 subgroup elements are represented by ordered tuples.
"""

import json
import math


Q53 = 21044858204113
DEGREE53 = 53
ORBIT_COLUMNS53 = 220
POINTS53 = 23320
REGULAR_STATES53 = 2565200
UNIQUE_ROOTS53 = 5081560

Q131 = 680564733841876926932320129493409985129
DEGREE131 = 131
AUTOMORPHISMS131 = 2 * DEGREE131
PACKED_ROOT_BYTES = 24  # 17-byte field root + 7-byte pair/shift witness.


def ceil_div(a: int, b: int) -> int:
    return (a + b - 1) // b


def ceil_fourth_root_ratio(numerator: int, denominator: int) -> int:
    """Least nonnegative F with denominator*F**4 >= numerator."""
    lo, hi = 0, 1
    while denominator * hi**4 < numerator:
        hi *= 2
    while lo < hi:
        mid = (lo + hi) // 2
        if denominator * mid**4 >= numerator:
            hi = mid
        else:
            lo = mid + 1
    return lo


def orbit_scenario(name: str, columns: int) -> dict:
    n = DEGREE131
    points = 2 * n * columns
    slots = n * columns**2  # Every (left, right, relative shift) scan slot.
    root_candidates = 2 * slots  # Conditional: two regular S3 roots per slot.
    coverage_cap_numerator = min(points**4, Q131)
    probability_cap_approx = coverage_cap_numerator / Q131
    s3_calls_per_failed_query = 2 * n * slots
    failed_call_numerator = (
        columns * (Q131 - coverage_cap_numerator) * s3_calls_per_failed_query
    )
    failed_call_denominator = coverage_cap_numerator
    return {
        "scenario": name,
        "columns": columns,
        "points": points,
        "regular_scan_slots": slots,
        "root_candidates_if_two_per_slot": root_candidates,
        "packed_root_bytes_if_distinct": root_candidates * PACKED_ROOT_BYTES,
        "random_target_coverage_upper_bound_numerator": coverage_cap_numerator,
        "random_target_coverage_upper_bound_denominator": Q131,
        "random_target_coverage_upper_bound_approx": probability_cap_approx,
        "log2_random_target_coverage_upper_bound_approx": (
            math.log2(coverage_cap_numerator) - math.log2(Q131)
        ),
        "s3_calls_per_failed_query_if_all_regular": s3_calls_per_failed_query,
        "expected_failed_query_s3_calls_for_full_rank_numerator": failed_call_numerator,
        "expected_failed_query_s3_calls_for_full_rank_denominator": failed_call_denominator,
        "log2_expected_failed_query_s3_calls_for_full_rank_approx": (
            math.log2(failed_call_numerator) - math.log2(failed_call_denominator)
            if failed_call_numerator
            else None
        ),
    }


def main() -> None:
    assert POINTS53 == 2 * DEGREE53 * ORBIT_COLUMNS53
    assert REGULAR_STATES53 == DEGREE53 * ORBIT_COLUMNS53**2
    assert Q131.bit_length() == 130
    assert 2 * REGULAR_STATES53 > UNIQUE_ROOTS53

    one_percent_f = ceil_fourth_root_ratio(Q131, 100)
    half_f = ceil_fourth_root_ratio(Q131, 2)
    one_tib = 1 << 40
    one_tib_columns = math.isqrt(one_tib // (2 * DEGREE131 * PACKED_ROOT_BYTES))
    scenarios = [
        orbit_scenario("same_220_columns", ORBIT_COLUMNS53),
        orbit_scenario("one_TiB_packed_root_budget", one_tib_columns),
        orbit_scenario("one_percent_counting_threshold", ceil_div(one_percent_f, AUTOMORPHISMS131)),
        orbit_scenario("half_counting_threshold", ceil_div(half_f, AUTOMORPHISMS131)),
    ]
    assert scenarios[1]["packed_root_bytes_if_distinct"] <= one_tib
    assert orbit_scenario("next", one_tib_columns + 1)["packed_root_bytes_if_distinct"] > one_tib

    rho_single = math.sqrt(math.pi * Q131 / (2 * AUTOMORPHISMS131))
    rho_batch1024 = math.sqrt(2 * 1024 * Q131 / AUTOMORPHISMS131)
    print(json.dumps({
        "status": "conditional_projection_not_n131_measurement",
        "q53": Q53,
        "n53_recorded": {
            "degree": DEGREE53,
            "columns": ORBIT_COLUMNS53,
            "points": POINTS53,
            "regular_states": REGULAR_STATES53,
            "unique_root_entries": UNIQUE_ROOTS53,
            "root_occupancy_over_two_per_state": UNIQUE_ROOTS53 / (2 * REGULAR_STATES53),
            "four_sum_counting_ratio_F4_over_q": POINTS53**4 / Q53,
        },
        "q131": Q131,
        "degree131": DEGREE131,
        "signed_frobenius_order": AUTOMORPHISMS131,
        "packed_root_record_bytes_assumed": PACKED_ROOT_BYTES,
        "one_percent_minimum_F_before_orbit_rounding": one_percent_f,
        "half_minimum_F_before_orbit_rounding": half_f,
        "rho_single_nominal_additions": rho_single,
        "rho_batch1024_nominal_additions": rho_batch1024,
        "scenarios": scenarios,
    }, indent=2))


if __name__ == "__main__":
    main()

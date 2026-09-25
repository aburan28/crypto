#!/usr/bin/env python3
"""Exact multiset-count ceiling for the current four-factor-base-point design.

Each unordered 4-multiset of F base points has one group sum, so even an ideal
exhaustive extractor covers at most C(F+3,4) of the q subgroup elements.
This is a necessary bound, not a predicted hit rate or a measured cost.
"""

from __future__ import annotations

import json
from math import comb

Q131 = 680564733841876926932320129493409985129
N131 = 131
ROOT_BYTES = 24  # Same optimistic packed record as PR #739.
TIB = 1 << 40


def smallest_f_for_fraction(q: int, numerator: int, denominator: int,
                            summands: int = 4) -> int:
    lo, hi = 1, 1
    while comb(hi + summands - 1, summands) * denominator < q * numerator:
        hi *= 2
    while lo < hi:
        mid = (lo + hi) // 2
        if comb(mid + summands - 1, summands) * denominator >= q * numerator:
            hi = mid
        else:
            lo = mid + 1
    assert comb(lo + summands - 2, summands) * denominator < q * numerator
    return lo


def scenario(label: str, columns: int) -> dict:
    f = 2 * N131 * columns
    slots = N131 * columns * columns
    candidates = 2 * slots
    return {"label": label, "R": columns, "F": f,
            "ordered_tuple_count": f**4,
            "unordered_multiset_count": comb(f + 3, 4),
            "ordered_probability_upper_display": f**4 / Q131,
            "unordered_probability_upper_display": comb(f + 3, 4) / Q131,
            "regular_scan_slots_if_all_regular": slots,
            "two_root_candidates_if_all_regular": candidates,
            "optimistic_packed_root_bytes": candidates * ROOT_BYTES}


def main() -> None:
    one_tib_r = 0
    while (2*N131*(one_tib_r+1)**2)*ROOT_BYTES <= TIB:
        one_tib_r += 1
    assert one_tib_r == 13223
    f_one_pct = smallest_f_for_fraction(Q131, 1, 100)
    f_half = smallest_f_for_fraction(Q131, 1, 2)
    assert f_one_pct == 3574951633 and f_half == 9506325303
    alternative_thresholds = {str(m): smallest_f_for_fraction(Q131, 1, 100, m)
                              for m in (4, 5, 6)}
    assert alternative_thresholds == {"4": 3574951633, "5": 60591280, "6": 4121293}
    one_pct_r = (f_one_pct + 2*N131 - 1) // (2*N131)
    half_r = (f_half + 2*N131 - 1) // (2*N131)
    assert one_pct_r == 13644854 and half_r == 36283685
    scenarios = [scenario("one_TiB_optimistic_packed_roots", one_tib_r),
                 scenario("necessary_one_percent_coverage", one_pct_r),
                 scenario("necessary_half_coverage", half_r)]
    assert scenarios[0]["optimistic_packed_root_bytes"] <= TIB
    assert scenarios[1]["optimistic_packed_root_bytes"] == 1170712671804115008
    assert scenarios[2]["optimistic_packed_root_bytes"] == 8278188452662966800
    toy = []
    for n, q, bases in ((37, 230603167, (1, 2, 3)),
                        (41, 549756390943, (4, 8, 12))):
        for r in bases:
            f = 2*n*r
            toy.append({"n": n, "q": q, "R": r, "F": f,
                        "ordered_bound_full_group": f**4/q,
                        "unordered_bound_full_group": comb(f+3, 4)/q,
                        "unordered_bound_nonzero_targets": comb(f+3, 4)/(q-1)})
    print(json.dumps({"classification": "exact necessary counting bound, not a measurement",
                      "q131": Q131,
                      "smallest_unrounded_F_for_one_percent": f_one_pct,
                      "model_only_one_percent_minimum_F_by_summand_count": alternative_thresholds,
                      "smallest_unrounded_F_for_half": f_half,
                      "scenarios": scenarios, "toy_arms": toy},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

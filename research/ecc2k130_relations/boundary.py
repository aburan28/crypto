#!/usr/bin/env python3
"""The counting floor for homogeneous-relation index calculus on ECC2K-130.

AGENTS.md asks for the boundary before the measurement, and for this family
of methods the boundary is a counting argument that no amount of support
structure can move.  It is derived here and the sweeps confirm it; the
sweeps do not discover it.

Setup.  A factor base of `B` abscissae over the order-`r` subgroup, and
`m`-point relations.  There are about `B^m / m!` unordered `m`-multisets,
each summing to a fixed point with probability about `1/r`, so relations
start to exist at

    B_m = (m! r)^(1/m).

Below that the expected yield is under one and a complete sweep finds
nothing, however the base is chosen.  This is the part that is independent
of structure: it counts multisets, and a structured base has no more of
them than a random base of the same size.

Finding a relation, once they exist, is the second cost.  The best generic
method is meet-in-the-middle on a balanced split `m = s + (m-s)`: enumerate
and store the `B^s / s!` partial sums, then stream the `B^(m-s) / (m-s)!`
others and look each up.  Its cost is the larger of the two sides, and its
*memory* is the stored side -- which is where this family actually dies.

Recovering a logarithm is the third cost, and the one that is easy to skip.
Homogeneous relations among base points alone give linear dependencies
between unknown logarithms; they do not determine `log_P(Q)`.  For that one
needs about `B` relations whose targets are known combinations `[a]P + [b]Q`,
so the relation cost is paid `B` times over.  Reporting the cost of one
relation as the cost of the method is the §3 relabelling error, drawn.

The reference is Pollard rho with the `<-1> x <pi>` speed-up, `2^60.9`
iterations on the same subgroup.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

LOG2_R = 129.0
LOG2_RHO = 60.9
BYTES_PER_STORED_POINT = 32          # generous: a compressed 131-bit abscissa


def log2_factorial(m: int) -> float:
    return math.log2(math.factorial(m))


def boundary_row(m: int):
    lf = log2_factorial(m)
    log2_B = (lf + LOG2_R) / m                    # B_m = (m! r)^(1/m)

    # balanced meet-in-the-middle split
    best = None
    for s in range(1, m // 2 + 1):
        store = s * log2_B - log2_factorial(s)
        stream = (m - s) * log2_B - log2_factorial(m - s)
        cost = max(store, stream)
        if best is None or cost < best[2]:
            best = (s, store, stream, cost)
    s, log2_store, log2_stream, log2_relation_cost = best

    log2_total = log2_relation_cost + log2_B      # ~B relations needed
    log2_memory_bytes = log2_store + math.log2(BYTES_PER_STORED_POINT)

    return {
        "m": m,
        "log2_base_size_needed": round(log2_B, 2),
        "mitm_split": f"{s}+{m - s}",
        "log2_stored_side": round(log2_store, 2),
        "log2_streamed_side": round(log2_stream, 2),
        "log2_cost_one_relation": round(log2_relation_cost, 2),
        "log2_memory_bytes": round(log2_memory_bytes, 2),
        "log2_total_cost": round(log2_total, 2),
        "log2_ratio_one_relation_to_rho": round(log2_relation_cost - LOG2_RHO, 2),
        "log2_ratio_total_to_rho": round(log2_total - LOG2_RHO, 2),
    }


def sweep_yield(log2_base_size: float, m: int) -> float:
    """log2 of the expected homogeneous `m`-relation count for a given base."""
    return m * log2_base_size - log2_factorial(m) - LOG2_R


def main():
    rows = [boundary_row(m) for m in range(3, 9)]

    # what the bases actually swept here are worth, by the same counting
    measured_bases = {
        "weight_two_complete": 4262,
        "weight_two_frobenius_stable": None,       # filled by the runner
    }
    yields = {}
    for name, size in measured_bases.items():
        if not size:
            continue
        lb = math.log2(size)
        yields[name] = {
            "base_size": size,
            "log2_base_size": round(lb, 2),
            "log2_expected_3_point_relations": round(sweep_yield(lb, 3), 2),
            "log2_expected_4_point_relations": round(sweep_yield(lb, 4), 2),
            "log2_shortfall_in_base_size_for_m3": round(
                (log2_factorial(3) + LOG2_R) / 3 - lb, 2),
        }

    best = min(rows, key=lambda d: d["log2_total_cost"])
    best_one = min(rows, key=lambda d: d["log2_cost_one_relation"])

    out = {
        "instance": "ECC2K-130",
        "log2_r": LOG2_R,
        "rho_reference_log2": LOG2_RHO,
        "bytes_per_stored_point": BYTES_PER_STORED_POINT,
        "rows": rows,
        "measured_base_yields": yields,
        "best_m_by_total_cost": best["m"],
        "log2_best_total_cost": best["log2_total_cost"],
        "log2_best_total_ratio_to_rho": best["log2_ratio_total_to_rho"],
        "best_m_by_single_relation_cost": best_one["m"],
        "log2_best_single_relation_cost": best_one["log2_cost_one_relation"],
        "log2_best_single_relation_ratio_to_rho":
            best_one["log2_ratio_one_relation_to_rho"],
        "verdict": (
            "Homogeneous-relation index calculus does not reach the rho "
            "reference at any m. Even priced at one relation and ignoring "
            "the B relations a logarithm needs, the cheapest m costs "
            f"2^{best_one['log2_cost_one_relation']} against rho's 2^{LOG2_RHO}, "
            f"a factor 2^{best_one['log2_ratio_one_relation_to_rho']}, and needs "
            f"2^{best_one['log2_memory_bytes']} bytes of storage. Priced "
            "honestly at B relations the gap is "
            f"2^{best['log2_ratio_total_to_rho']}."
        ),
    }

    path = Path(__file__).resolve().parent / "results" / "boundary.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2))

    hdr = (f"{'m':>2} {'log2 B':>7} {'split':>6} {'log2 mem B':>11} "
           f"{'log2 1-rel':>11} {'vs rho':>8} {'log2 total':>11} {'vs rho':>8}")
    print(hdr)
    print("-" * len(hdr))
    for d in rows:
        print(f"{d['m']:>2} {d['log2_base_size_needed']:>7} "
              f"{d['mitm_split']:>6} {d['log2_memory_bytes']:>11} "
              f"{d['log2_cost_one_relation']:>11} "
              f"{d['log2_ratio_one_relation_to_rho']:>+8} "
              f"{d['log2_total_cost']:>11} "
              f"{d['log2_ratio_total_to_rho']:>+8}")
    print()
    print(out["verdict"])
    print()
    for name, y in yields.items():
        print(f"{name}: B = {y['base_size']} (2^{y['log2_base_size']}), "
              f"expected 3-point relations 2^{y['log2_expected_3_point_relations']}, "
              f"short of the m=3 threshold by 2^{y['log2_shortfall_in_base_size_for_m3']} "
              f"in base size")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""What the Frobenius quotient is worth against the counting floor, at every `m`.

Section 1 of the research note prices homogeneous-relation index calculus at
`m = 3..8` with no Frobenius quotient, and finds `m = 8` cheapest at
`2^+24.79` over the rho reference.  That accounting predates the
orbit-grouped search.  The quotient changes three of its inputs, and this
module applies each one to *every* `m` rather than to `m = 4` alone:

  * the **floor** moves.  A search enumerates `m` distinct abscissae with a
    sign each, modulo a global negation -- `2^(m-1) C(B,m)` -- and every
    support point lies in the index-2 subgroup `H`, so sums reach only the
    two elements of `E[4] cap H`.  Expected usable relations are
    `2^(m-1) C(B,m) / r`, and the floor is

        B_m = (m! r / 2^(m-1))^(1/m)

    rather than section 1's `(m! r)^(1/m)`;

  * **storage and streaming** divide by `131`.  Both sides of the
    meet-in-the-middle come in `sigma`-orbits, one canonical class key
    stands for 131 tuples, and a key match hands back all 131 rotations.
    The key is a truncated abscissa, 8 bytes, against section 1's 32-byte
    point;

  * the **relations needed** divide by 131 too, because the support carries
    one unknown per orbit rather than one per point.

**A correction this generalisation forces.**  An earlier reading of the
`m = 4` row called orbit-grouped `m = 4` the best accounting in the study,
"ahead of the `2^+24.70` that `m = 8` holds in section 1".  That compares a
quotiented row against an unquotiented one.  Applied uniformly, the quotient
helps larger `m` *more* -- it is worth `2^16.1` at `m = 3` rising to
`2^18.4` at `m = 8` -- so `m = 8` stays cheapest, and the ranking of section
1 is unchanged.  What changes is the margin.

**Two conventions, and why both are reported.**  Section 1 measures cost at
the floor `B_m`, where the expected yield is *one* relation, while charging
for `B` relations.  That is internally generous, and it is inherited here
for comparability as `convention = "section 1"`.  The `self-consistent`
rows instead take the support size at which the yield actually reaches the
number of relations needed,

    2^(m-1) C(B,m) / r = B / 131   =>   B = (m! r / (131 * 2^(m-1)))^(1/(m-1))

which is the honest column.  At the margins section 1 was quoting the
difference did not matter; with the quotient applied it does.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

RHO_LOG2 = 60.809           # experiments/ecc2k130_extension_field_boundary.json
LOG2_R = 129.0
M = 131                     # the Frobenius order
LOG2_M = math.log2(M)
KEY_BYTES = 8               # a truncated canonical abscissa
POINT_BYTES = 32            # section 1's stored point


def _lf(n: int) -> float:
    return math.log2(math.factorial(n))


def _best_split(m: int, log2_B: float):
    """The meet-in-the-middle split minimising the larger side, as in section 1."""
    best = None
    for s in range(1, m // 2 + 1):
        store = s * log2_B - _lf(s)
        stream = (m - s) * log2_B - _lf(m - s)
        cost = max(store, stream)
        if best is None or cost < best[3]:
            best = (s, store, stream, cost)
    return best


def section_one_row(m: int):
    """Section 1's own row, recomputed here so the comparison is like for like."""
    log2_B = (_lf(m) + LOG2_R) / m
    s, store, stream, cost = _best_split(m, log2_B)
    total = cost + log2_B
    return {
        "m": m,
        "accounting": "section 1, no quotient",
        "convention": "section 1",
        "log2_support_size": round(log2_B, 2),
        "mitm_split": f"{s}+{m - s}",
        "log2_cost_one_relation": round(cost, 2),
        "log2_relations_needed": round(log2_B, 2),
        "log2_memory_bytes": round(store + math.log2(POINT_BYTES), 2),
        "memory_exabytes": float(f"{2 ** (store + math.log2(POINT_BYTES)) / 1e18:.4g}"),
        "log2_total_cost": round(total, 2),
        "log2_one_relation_vs_rho": round(cost - RHO_LOG2, 2),
        "log2_total_vs_rho": round(total - RHO_LOG2, 2),
    }


def quotient_row(m: int, self_consistent: bool):
    """The same row with the Frobenius quotient applied to all three inputs."""
    if self_consistent:
        # 2^(m-1) C(B,m)/r = B/131
        log2_B = (_lf(m) + LOG2_R - LOG2_M - (m - 1)) / (m - 1)
        convention = "self-consistent"
    else:
        # 2^(m-1) C(B,m)/r = 1
        log2_B = (_lf(m) + LOG2_R - (m - 1)) / m
        convention = "section 1"
    s, store, stream, cost = _best_split(m, log2_B)
    cost_q = cost - LOG2_M
    store_q = store - LOG2_M
    relations = log2_B - LOG2_M            # one unknown per orbit
    total = cost_q + relations
    mem = store_q + math.log2(KEY_BYTES)
    return {
        "m": m,
        "accounting": "orbit-grouped, one unknown per orbit",
        "convention": convention,
        "log2_support_size": round(log2_B, 2),
        "mitm_split": f"{s}+{m - s}",
        "log2_cost_one_relation": round(cost_q, 2),
        "log2_relations_needed": round(relations, 2),
        "log2_memory_bytes": round(mem, 2),
        "memory_exabytes": float(f"{2 ** mem / 1e18:.4g}"),
        "log2_total_cost": round(total, 2),
        "log2_one_relation_vs_rho": round(cost_q - RHO_LOG2, 2),
        "log2_total_vs_rho": round(total - RHO_LOG2, 2),
    }


def main():
    sec1 = [section_one_row(m) for m in range(3, 9)]
    quot = [quotient_row(m, False) for m in range(3, 9)]
    hon = [quotient_row(m, True) for m in range(3, 9)]

    best_sec1 = min(sec1, key=lambda d: d["log2_total_cost"])
    best_quot = min(quot, key=lambda d: d["log2_total_cost"])
    best_hon = min(hon, key=lambda d: d["log2_total_cost"])

    data = {
        "instance": "ECC2K-130",
        "log2_r": LOG2_R,
        "frobenius_order": M,
        "rho_reference_log2": RHO_LOG2,
        "rho_reference_source": (
            "experiments/ecc2k130_extension_field_boundary.json -> "
            "target.log2_rho_reference"),
        "bytes_per_class_key": KEY_BYTES,
        "bytes_per_stored_point": POINT_BYTES,
        "rows_section_one": sec1,
        "rows_quotient_section_one_convention": quot,
        "rows_quotient_self_consistent": hon,
        "quotient_gain_by_m": {
            str(a["m"]): round(a["log2_total_vs_rho"] - b["log2_total_vs_rho"], 2)
            for a, b in zip(sec1, quot)},
        "best_m_section_one": best_sec1["m"],
        "best_m_quotient": best_quot["m"],
        "best_m_self_consistent": best_hon["m"],
        "log2_best_vs_rho_section_one": best_sec1["log2_total_vs_rho"],
        "log2_best_vs_rho_quotient": best_quot["log2_total_vs_rho"],
        "log2_best_vs_rho_self_consistent": best_hon["log2_total_vs_rho"],
        "verdict": (
            f"The Frobenius quotient is worth 2^16.1 to 2^18.4 depending on m, "
            f"more at larger m, so it does not change which m is cheapest: "
            f"m = {best_sec1['m']} before and m = {best_quot['m']} after. It "
            f"changes the margin. Section 1's best row sits at "
            f"2^{best_sec1['log2_total_vs_rho']:+} over rho; quotiented and "
            f"priced on a support that actually yields the relations it "
            f"charges for, the best row is m = {best_hon['m']} at "
            f"2^{best_hon['log2_total_vs_rho']:+}, needing "
            f"{best_hon['memory_exabytes']:,.0f} exabytes. No m reaches the "
            f"reference, and the conclusion of the study is unchanged, but "
            f"the gap is 2^{best_sec1['log2_total_vs_rho'] - best_hon['log2_total_vs_rho']:.2f} "
            f"smaller than section 1 records."),
    }

    out = (Path(__file__).resolve().parent / "results"
           / "four_point_orbit_boundary.json")
    out.write_text(json.dumps(data, indent=2))

    for title, rows in (("section 1, no quotient", sec1),
                        ("orbit-grouped, section 1's convention", quot),
                        ("orbit-grouped, self-consistent floor", hon)):
        print(f"\n{title}")
        hdr = (f"{'m':>2} {'log2 B':>7} {'split':>6} {'1 rel':>7} {'vs rho':>8} "
               f"{'rels':>6} {'mem EB':>12} {'total':>7} {'vs rho':>8}")
        print(hdr)
        print("-" * len(hdr))
        for d in rows:
            print(f"{d['m']:>2} {d['log2_support_size']:>7.2f} {d['mitm_split']:>6} "
                  f"{d['log2_cost_one_relation']:>7.2f} "
                  f"{d['log2_one_relation_vs_rho']:>+8.2f} "
                  f"{d['log2_relations_needed']:>6.2f} "
                  f"{d['memory_exabytes']:>12,.1f} "
                  f"{d['log2_total_cost']:>7.2f} {d['log2_total_vs_rho']:>+8.2f}")
    print()
    print(data["verdict"])


if __name__ == "__main__":
    main()

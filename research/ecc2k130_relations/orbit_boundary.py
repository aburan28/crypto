#!/usr/bin/env python3
"""What the Frobenius quotient is worth against the `m = 4` cost boundary.

Section 1 of the research note prices `m = 4` at `2^65.79` for one relation
and `2^99.19` for a logarithm, against a rho reference of `2^60.809`.  That
accounting predates the orbit-grouped search, and the quotient changes three
of its inputs at once:

  * the **floor** moves.  Section 1 counts multisets, `B^4 / 4!`.  A search
    enumerates four distinct abscissae with a sign each, modulo a global
    negation -- `8 C(B,4)` -- and accepts any sum in `E[4] cap H`, which is
    two of the `2r` elements the sums can reach.  Expected relations are
    `8 C(B,4) / r`, so the floor is `B_4 = (3r)^(1/4)`, a little below the
    `(4! r)^(1/4)` of section 1;

  * **storage and streaming** divide by `m = 131`.  Pair sums come in
    `sigma`-orbits, one canonical key stands for 131 pairs, and a key
    collision hands back all 131 rotations.  The key is a truncated
    abscissa, 8 bytes, rather than a 32-byte point;

  * the **number of relations needed** divides by 131 as well, because the
    support carries one unknown per orbit rather than one per point.

The first two make the single-relation row cheaper than rho.  That row is
exactly the one section 1 warns against quoting: the work has not gone away,
it has moved into the count of relations a logarithm actually needs.  Priced
there, `m = 4` remains above the reference.

Section 1's convention is kept throughout -- cost is measured at the floor
`B_4`, and a logarithm is priced at one unknown's worth of relations -- so
that the rows can be read against each other.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

RHO_LOG2 = 60.809           # experiments/ecc2k130_extension_field_boundary.json
LOG2_R = 129.0
M = 131                     # the Frobenius order


def _row(name, note, log2_B, log2_stored, bytes_per_entry, log2_relations):
    log2_mem = log2_stored + math.log2(bytes_per_entry)
    log2_total = log2_stored + log2_relations
    return {
        "accounting": name,
        "note": note,
        "log2_support_size": round(log2_B, 3),
        "log2_cost_one_relation": round(log2_stored, 3),
        "bytes_per_stored_entry": bytes_per_entry,
        "log2_memory_bytes": round(log2_mem, 3),
        "memory_exabytes": round(2 ** log2_mem / 1e18, 1),
        "log2_relations_needed": round(log2_relations, 3),
        "log2_total_cost": round(log2_total, 3),
        "log2_one_relation_vs_rho": round(log2_stored - RHO_LOG2, 3),
        "log2_total_vs_rho": round(log2_total - RHO_LOG2, 3),
    }


def rows():
    out = []

    # (a) Section 1, as published: multiset floor, 2+2 split, points stored.
    b_a = math.log2((math.factorial(4) * 2 ** LOG2_R) ** 0.25)
    out.append(_row(
        "section 1, as published",
        "multiset floor (4! r)^(1/4); B^2/2 pair sums stored as 32-byte points; "
        "B relations for a logarithm",
        b_a, 2 * b_a - 1, 32, b_a))

    # (b) Corrected floor: signed, distinct abscissae, E[4] cap H accepted.
    b_b = math.log2((3 * 2 ** LOG2_R) ** 0.25)
    out.append(_row(
        "corrected floor",
        "8 C(B,4)/r = 1 gives B_4 = (3r)^(1/4); search space and cofactor "
        "corrected, no Frobenius quotient yet",
        b_b, 2 * b_b - 1, 32, b_b))

    # (c) Frobenius quotient on storage and streaming only.
    out.append(_row(
        "+ orbit-grouped pair sums",
        "one canonical class key per sigma-orbit of pairs: 131x fewer entries, "
        "8-byte keys; unknowns still counted per point",
        b_b, 2 * b_b - 1 - math.log2(M), 8, b_b))

    # (d) And the collapse of the unknowns.
    out.append(_row(
        "+ one unknown per orbit",
        "log(sigma^k R) = s^k log(R), so a B-point support carries B/131 "
        "unknowns and needs B/131 relations",
        b_b, 2 * b_b - 1 - math.log2(M), 8, b_b - math.log2(M)))

    return out


def main():
    data = {
        "instance": "ECC2K-130",
        "m": 4,
        "log2_r": LOG2_R,
        "frobenius_order": M,
        "rho_reference_log2": RHO_LOG2,
        "rho_reference_source": (
            "experiments/ecc2k130_extension_field_boundary.json -> "
            "target.log2_rho_reference"),
        "convention": (
            "cost measured at the floor B_4, as in section 1 of "
            "RESEARCH_ECC2K130_RELATION_SWEEPS.md; a logarithm priced at one "
            "relation per unknown"),
        "rows": rows(),
    }
    out = Path(__file__).resolve().parent / "results" / "four_point_orbit_boundary.json"
    out.write_text(json.dumps(data, indent=2))
    hdr = f"{'accounting':<28} {'log2 B':>7} {'1 rel':>8} {'vs rho':>8} {'mem EB':>10} {'total':>8} {'vs rho':>8}"
    print(hdr)
    print("-" * len(hdr))
    for r in data["rows"]:
        print(f"{r['accounting']:<28} {r['log2_support_size']:>7.2f} "
              f"{r['log2_cost_one_relation']:>8.2f} "
              f"{r['log2_one_relation_vs_rho']:>+8.2f} "
              f"{r['memory_exabytes']:>10.1f} "
              f"{r['log2_total_cost']:>8.2f} {r['log2_total_vs_rho']:>+8.2f}")


if __name__ == "__main__":
    main()

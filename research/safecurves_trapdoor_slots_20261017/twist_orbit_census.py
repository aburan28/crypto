#!/usr/bin/env sage
"""Round 0 census: the full sextic-twist orbit of the SECG j = 0 curves.

For j = 0 curves y^2 = x^3 + b over F_p with p = 1 mod 6 there are six
F_p-isomorphism classes (the sextic twists), with six *different* group
orders.  SafeCurves-style twist criteria certify exactly one of the six: the
quadratic twist of the published curve.  This census factors all six orders
for secp192k1, secp224k1 and secp256k1 and records the rho cost of each.

Outputs census.md and census.json in the directory of this script.
"""
import json
import os
import sys

from sage.all import GF, EllipticCurve, factor, is_prime, next_prime, proof

proof.arithmetic(False)  # BPSW probable primes only; recorded as such.

HERE = os.path.dirname(os.path.abspath(__file__))

CURVES = {
    "secp192k1": (2**192 - 2**32 - 2**12 - 2**8 - 2**7 - 2**6 - 2**3 - 1, 3),
    "secp224k1": (2**224 - 2**32 - 2**12 - 2**11 - 2**9 - 2**7 - 2**4 - 2 - 1, 5),
    "secp256k1": (2**256 - 2**32 - 977, 7),
}


def largest_factor(n):
    f = factor(n)
    return int(f[-1][0]), f


def orbit(p_int, b_int):
    """Return the six sextic-twist classes of y^2 = x^3 + b over F_p.

    Enumerates the six distinct orders reachable as b * d, d in F_p*,
    tagged with the smallest witness d found.
    """
    F = GF(p_int)
    p = F.characteristic()
    assert p % 6 == 1, "p must be 1 mod 6 for a six-class orbit"
    found = {}
    d = 1
    while len(found) < 6:
        E = EllipticCurve(F, [0, F(b_int) * F(d)])
        n = int(E.order())
        if n not in found:
            found[n] = d
        d += 1
    return found


def census():
    rows = []
    for name, (p_int, b_int) in CURVES.items():
        found = orbit(p_int, b_int)
        F = GF(p_int)
        n_pub = int(EllipticCurve(F, [0, F(b_int)]).order())
        n_qtwist = 2 * (p_int + 1) - n_pub  # order of the quadratic twist
        sum_check = sum(found) == 6 * (p_int + 1)
        row = {
            "curve": name,
            "p_hex": hex(p_int),
            "b": b_int,
            "p_mod_6": int(p_int % 6),
            "published_order": n_pub,
            "published_order_prime": bool(is_prime(n_pub)),
            "quadratic_twist_order": n_qtwist,
            "orbit_sums_to_6p_plus_6": sum_check,
            "classes": [],
        }
        for n, d in sorted(found.items()):
            lf, f = largest_factor(n)
            role = "published" if n == n_pub else "quadratic twist" if n == n_qtwist else "other sextic class"
            row["classes"].append(
                {
                    "witness_d": d,
                    "order": n,
                    "factorization": str(f),
                    "largest_factor_bits": int(lf.nbits() if hasattr(lf, "nbits") else 0) or len(bin(lf)) - 2,
                    "rho_cost_bits": (len(bin(lf)) - 2) / 2.0,
                    "role": role,
                }
            )
        rows.append(row)
    return rows


def main():
    rows = census()
    with open(os.path.join(HERE, "census.json"), "w") as fh:
        json.dump(rows, fh, indent=1)
    lines = [
        "# Sextic-twist orbit census of the SECG j = 0 curves (round 0)",
        "",
        "Executed by `twist_orbit_census.py` (SageMath; probable-prime factorisations,",
        "BPSW, `proof.Arithmetic(False)`). Nothing is pinned by hand; every order is",
        "recomputed by SEA point counting.",
        "",
        "For each curve: six F_p-isomorphism classes y^2 = x^3 + b*d, tagged with the",
        "smallest witness d reaching that order. `rho_cost_bits` is half the bit length",
        "of the largest prime factor, the SafeCurves-style rho cost *on that twist*.",
        "",
    ]
    for row in rows:
        lines.append(f"## {row['curve']}  (p mod 6 = {row['p_mod_6']}, orbit sums to 6(p+1): {row['orbit_sums_to_6p_plus_6']})")
        lines.append("")
        lines.append(f"- published order prime: {row['published_order_prime']}")
        lines.append("")
        lines.append("| role | witness d | order | largest factor (bits) | rho cost (bits) | factorization |")
        lines.append("|---|---:|---|---:|---:|---|")
        for c in row["classes"]:
            lines.append(
                f"| {c['role']} | {c['witness_d']} | {c['order']} | {c['largest_factor_bits']} | {c['rho_cost_bits']:.1f} | {c['factorization']} |"
            )
        lines.append("")
    with open(os.path.join(HERE, "census.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("wrote census.md / census.json")


if __name__ == "__main__":
    sys.exit(main())

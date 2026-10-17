#!/usr/bin/env sage
"""Round 0 grind demonstration: a j = 0 curve family where the *quadratic*
twist passes a SafeCurves-strength bar while a *sextic* sibling class is weak.

Model: the designer chooses the prime p = 1 mod 6 (free in SECG-style
generation) and the sextic class of b.  Published-class checks mirror
SafeCurves, scaled to toy size:

  - published order has a large prime factor (>= 0.79 * bitlen(p)),
  - the quadratic twist of the published class (its +/- trace pair) also has
    a large prime factor,
  - but some *other* sextic class has largest prime factor <= 0.40 * bitlen(p),
    i.e. its rho cost is ~ (0.40/2) * bitlen(p) bits -- the withheld trapdoor.

Outputs grind_demo.md / grind_demo.json.
"""
import json
import os
import sys

from sage.all import GF, EllipticCurve, ZZ, factor, is_prime, next_prime, proof

proof.arithmetic(False)

HERE = os.path.dirname(os.path.abspath(__file__))

SIZES = [
    # (target bit length, strong fraction, weak fraction, max primes to scan)
    (48, 0.79, 0.40, 4000),
    (96, 0.79, 0.40, 3000),
]


def bits(x):
    return len(bin(int(x))) - 2


def largest_factor_bits(n):
    f = factor(n)
    return bits(f[-1][0])


def scan(size, strong_frac, weak_frac, max_primes):
    p = next_prime(2**size)
    scanned = 0
    strong_bits = strong_frac * size
    weak_bits = weak_frac * size
    while scanned < max_primes:
        if p % 6 == 1:
            scanned += 1
            F = GF(p)
            orders = {}
            d = 1
            while len(orders) < 6 and d < 40:
                E = EllipticCurve(F, [0, F(d)])
                n = int(E.order())
                if n not in orders:
                    orders[n] = d
                d += 1
            if len(orders) < 6:
                p = next_prime(p)
                continue
            ns = sorted(orders)
            lfs = {n: largest_factor_bits(n) for n in ns}
            # classes pair as (p+1-t, p+1+t): i and 5-i in the sorted list,
            # since the six orders are p+1 +/- t_i.
            for i, n_pub in enumerate(ns[:3]):
                n_twist = ns[5 - i]
                if lfs[n_pub] >= strong_bits and lfs[n_twist] >= strong_bits:
                    weak = [(n, lf) for n, lf in lfs.items() if lf <= weak_bits and n not in (n_pub, n_twist)]
                    if weak:
                        n_weak, lf_weak = min(weak, key=lambda kv: kv[1])
                        return {
                            "bit_length": size,
                            "p": int(p),
                            "scanned_primes": scanned,
                            "published_d": orders[n_pub],
                            "published_order": n_pub,
                            "published_largest_factor_bits": lfs[n_pub],
                            "quadratic_twist_order": n_twist,
                            "quadratic_twist_largest_factor_bits": lfs[n_twist],
                            "weak_d": orders[n_weak],
                            "weak_order": n_weak,
                            "weak_largest_factor_bits": lf_weak,
                            "weak_rho_cost_bits": lf_weak / 2.0,
                            "all_classes": [
                                {"d": orders[n], "order": n, "largest_factor_bits": lfs[n]}
                                for n in ns
                            ],
                        }
        p = next_prime(p)
    return None


def main():
    out = []
    for size, sf, wf, mp in SIZES:
        r = scan(size, sf, wf, mp)
        if r is None:
            r = {"bit_length": size, "found": False, "scanned_primes": mp}
        else:
            r["found"] = True
        out.append(r)
    with open(os.path.join(HERE, "grind_demo.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    lines = [
        "# Grind demonstration: strong quadratic twist, weak sextic sibling",
        "",
        "Executed by `sextic_grind_demo.py`. For each size, primes p = 1 mod 6 were",
        "scanned (all six sextic classes of y^2 = x^3 + b over F_p computed by SEA);",
        "the first (p, class) pair is reported where the published class and its",
        "quadratic twist both have largest prime factor >= 0.79 * bitlen(p) bits while",
        "some other sextic class has largest prime factor <= 0.40 * bitlen(p).",
        "",
    ]
    for r in out:
        if not r.get("found"):
            lines.append(f"## {r['bit_length']} bits: not found in {r['scanned_primes']} primes")
            lines.append("")
            continue
        lines.append(f"## {r['bit_length']} bits (found after {r['scanned_primes']} primes = 1 mod 6)")
        lines.append("")
        lines.append(f"- p = {r['p']}")
        lines.append(f"- publish b in class d = {r['published_d']}: order {r['published_order']}, largest factor {r['published_largest_factor_bits']} bits")
        lines.append(f"- its quadratic twist (class pair): order {r['quadratic_twist_order']}, largest factor {r['quadratic_twist_largest_factor_bits']} bits")
        lines.append(f"- weak sextic sibling d = {r['weak_d']}: order {r['weak_order']}, largest factor {r['weak_largest_factor_bits']} bits -> rho cost ~ 2^{r['weak_rho_cost_bits']:.1f}")
        lines.append("")
        lines.append("| class d | order | largest factor (bits) |")
        lines.append("|---:|---|---:|")
        for c in r["all_classes"]:
            lines.append(f"| {c['d']} | {c['order']} | {c['largest_factor_bits']} |")
        lines.append("")
    with open(os.path.join(HERE, "grind_demo.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("wrote grind_demo.md / grind_demo.json")


if __name__ == "__main__":
    sys.exit(main())

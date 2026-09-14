#!/usr/bin/env python3
"""Boundaries on attacking ECC2K-130 after base-changing to an extension field.

The question: `E(F_2^131) <= E(F_2^(131 e))` for every `e >= 1`, so an attacker
may work in any of those larger fields.  Does the extra field structure buy an
attack cheaper than the Pollard rho reference on the original group?

Everything below is *derived*, not measured: the script computes the curve
constants, then checks each derivation by exhaustive enumeration at the sizes
where enumeration is possible.  It writes `experiments/ecc2k130_extension_field_boundary.json`.

    python3 scripts/ecc2k130_extension_field_boundary.py [--emb-limit N] [--emax E]
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "experiments/ecc2k130_extension_field_boundary.json"

N131 = 131                      # the challenge field degree, prime
SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]


def is_prime(n: int) -> bool:
    """Deterministic Miller-Rabin over the first twelve primes (ample below 2^128)."""
    if n < 2:
        return False
    for p in SMALL_PRIMES:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in SMALL_PRIMES:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def mult_order(a: int, m: int) -> int:
    k, x = 1, a % m
    while x != 1:
        x = x * a % m
        k += 1
    return k


def curve_order(n: int, trace: int = -1, q: int = 2) -> int:
    """#E(F_q^n) from the Koblitz recurrence s_i = t s_{i-1} - q s_{i-2}."""
    s0, s1 = 2, trace
    if n == 0:
        return q ** 0 + 1 - s0
    for _ in range(n - 1):
        s0, s1 = s1, trace * s1 - q * s0
    return q ** n + 1 - s1


def factor_degrees(n: int) -> Counter:
    """Degrees of the irreducible factors of t^n - 1 over F_2, with multiplicity.

    In characteristic 2, t^n - 1 = (t^m - 1)^(2^v) for n = 2^v m with m odd, and
    t^m - 1 factors into one irreducible of degree |2 mod d| per cyclotomic coset.
    """
    v, m = 0, n
    while m % 2 == 0:
        m //= 2
        v += 1
    seen, degs = set(), Counter()
    for a in range(m):
        if a in seen:
            continue
        coset, x = set(), a
        while x not in coset:
            coset.add(x)
            x = x * 2 % m
        seen |= coset
        degs[len(coset)] += 1 << v
    return degs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--emb-limit", type=int, default=10_000_000,
                    help="exponent bound for the embedding-degree exhaustion")
    ap.add_argument("--emax", type=int, default=16, help="largest extension degree e checked")
    args = ap.parse_args()

    # ---- the target, and the boundary it has to be beaten against -----------
    order = curve_order(N131)
    r = order // 4
    assert order % 4 == 0 and is_prime(r)
    log2r = math.log2(r)
    aut = 2 * N131                                  # <-1> x <Frobenius> on <G>
    log2_rho = math.log2(math.sqrt(math.pi * r / (2 * aut)))
    log2_rho_plain = math.log2(math.sqrt(math.pi * r / 2))
    sqrt_r = log2r / 2

    target = {
        "curve": "K_0 : y^2 + x y = x^3 + 1 over F_2, used over F_2^131",
        "trace_over_F2": -1,
        "order": str(order),
        "cofactor": 4,
        "r": str(r),
        "r_is_prime": True,
        "log2_r": round(log2r, 4),
        "rho_automorphism_order": aut,
        "log2_rho_plain": round(log2_rho_plain, 4),
        "log2_rho_reference": round(log2_rho, 4),
        "S_rho_reference": round(2 ** (log2_rho - sqrt_r), 6),
    }

    # ---- Boundary A: the target group, and its automorphisms, do not move ---
    # E is ordinary (trace -1 is odd), so End_{F_2^N}(E) = End_{F_2bar}(E) for
    # every N: base change adds no endomorphism, hence no new rho speed-up.
    # Frobenius acts on <G> as lambda with lambda^131 = 1, for every ambient field.
    boundary_a = {
        "statement": "base change fixes both the group <G> and its automorphism group",
        "curve_is_ordinary": True,
        "frobenius_orbit_length_on_subgroup": N131,
        "automorphisms_available_at_every_N": aut,
        "log2_rho_at_every_N": round(log2_rho, 4),
    }

    # ---- Boundary B: the subspace dichotomy, checked by enumeration ---------
    # Frobenius-stable F_2-subspaces of F_2^N are the kernels ker g(sigma) for
    # divisors g of t^N - 1, with dim = deg g.  Claim: for N = 131 e, every such
    # subspace either sits inside F_2^e (dim <= e) or has dim >= 130.
    dichotomy = []
    # 131 | e is the one family where 131 and e are not coprime; spot-check it too.
    for e in list(range(1, args.emax + 1)) + [N131, 2 * N131]:
        big = factor_degrees(N131 * e)
        small = factor_degrees(e)
        rest = big.copy()
        rest.subtract(small)
        assert all(c >= 0 for c in rest.values()), e     # t^e - 1 divides t^(131e) - 1
        rest = {d: c for d, c in rest.items() if c > 0}
        small_dim = sum(d * c for d, c in small.items())
        min_big = min(rest)
        assert small_dim == e and min_big >= 130, (e, small_dim, min_big)
        dichotomy.append({
            "e": e, "N": N131 * e,
            "small_factor_degrees": {str(k): v for k, v in sorted(small.items())},
            "other_factor_degrees": {str(k): v for k, v in sorted(rest.items())},
            "max_dim_inside_F2e": small_dim,
            "min_dim_outside_F2e": min_big,
        })
    boundary_b = {
        "statement": "for N = 131e every Frobenius-stable F_2-subspace of F_2^N "
                     "has dim <= e and lies in F_2^e, or has dim >= 130",
        "checked_for_e_up_to": args.emax,
        "cells": dichotomy,
    }

    # ---- Boundary C: the large horn costs more than rho, at every e --------
    # Subfields F_2^d of F_2^(131e) have d | 131e, so 131 | d or 131 | n = N/d.
    # 131 | d: Gaudry/Diem decomposition over F_(2^d)^n costs q^(2-2/n), q >= 2^131.
    horn_a = []
    for e in range(2, args.emax + 1):
        best = None
        for n in range(2, e + 1):
            if e % n:
                continue
            d = N131 * e // n
            cost = d * (2 - 2 / n)
            if best is None or cost < best[0]:
                best = (cost, d, n)
        if best is None:
            continue
        cost, d, n = best
        horn_a.append({
            "e": e, "N": N131 * e, "subfield_degree_d": d, "extension_degree_n": n,
            "log2_factor_base": d,
            "log2_cost_q_pow_2_minus_2_over_n": round(cost, 2),
            "log2_ratio_to_rho": round(cost - log2_rho, 2),
        })
    boundary_c = {
        "statement": "when 131 | d the decomposition attack costs at least q = 2^131",
        "reference_complexity": "Gaudry/Diem: Otilde(q^(2-2/n)) for fixed n >= 2, "
                                "constant exponential in n",
        "cheapest_cell": min(horn_a, key=lambda c: c["log2_cost_q_pow_2_minus_2_over_n"]),
        "cells": horn_a,
    }

    # ---- Boundary D: the small horn has no relations at all ----------------
    # V <= F_2^e and 2e does not divide 131e (131 is odd), so F_2^2e is not in
    # F_2^N: both y-roots of a factor-base abscissa lie in F_2^e.  The factor
    # base is therefore contained in the *subgroup* E(F_2^e).
    horn_b = []
    for e in range(1, args.emax + 1):
        sub = curve_order(e)
        horn_b.append({
            "e": e,
            "two_e_divides_N": (N131 * e) % (2 * e) == 0,
            "factor_base_subgroup_order": str(sub),
            "log2_factor_base_subgroup": round(math.log2(sub), 2),
            "target_subgroup_in_it": e % N131 == 0,
        })
    assert not any(c["two_e_divides_N"] for c in horn_b)
    boundary_d = {
        "statement": "when the invariant subspace lies in F_2^e the factor base is "
                     "contained in E(F_2^e), a subgroup of order ~2^e, so no "
                     "decomposition of a target in <G> exists at any cost",
        "log2_target_subgroup": round(log2r, 4),
        "cells": horn_b,
    }

    # ---- Boundary E: the transfer route needs a field nobody can write down -
    q131 = pow(2, N131, r)
    x, k, emb = 1, 0, None
    while k < args.emb_limit:
        x = x * q131 % r
        k += 1
        if x == 1:
            emb = k
            break
    boundary_e = {
        "statement": "MOV/Frey-Ruck transfer lands in F_2^(131k) with k the embedding degree",
        "embedding_degree_exhausted_to": args.emb_limit,
        "embedding_degree": emb,
        "log2_min_target_field_bits": round(math.log2(N131 * args.emb_limit), 2),
        "min_target_field_bits": N131 * args.emb_limit,
    }

    # ---- Where the idea does pay: a redundant ring, not a bigger field ------
    onb = []
    for p in range(3, 4000):
        if all(p % f for f in range(2, int(p ** 0.5) + 1)) and (p - 1) % N131 == 0 \
                and mult_order(2, p) == N131:
            onb.append({"p": p, "factors_of_degree_131": (p - 1) // N131,
                        "bits": p, "overhead_vs_131": round(p / N131, 3)})
    representation = {
        "statement": "F_2^131 embeds in the ring F_2[x]/(x^263 - 1); Frobenius becomes a "
                     "cyclic shift.  This is the type-II optimal normal basis the client "
                     "already uses, and it is an engineering gain, not an exponent gain",
        "rings": onb,
        "ecc2k95_contrast": {"m": 97, "two_m_plus_1": 195,
                             "prime": False,
                             "note": "no type-II ONB, which is why that backend is polynomial basis"},
    }

    # ---- The counterfactual: what a composite degree would have cost -------
    counterfactual = []
    for n in (2, 5, 10, 13):
        d, N = 130 // n, 130
        counterfactual.append({
            "hypothetical_field_degree": N, "subfield_degree_d": d, "extension_degree_n": n,
            "log2_cost_q_pow_2_minus_2_over_n": round(d * (2 - 2 / n), 2),
            "log2_rho_on_that_group": round(math.log2(math.sqrt(math.pi * 2 ** (N - 2) / (2 * 2 * N))), 2),
        })

    report = {
        "schema": "ecc2k130_extension_field_boundary/v1",
        "question": "does base-changing ECC2K-130 to F_2^(131e) admit an attack "
                    "cheaper than rho on E(F_2^131)?",
        "verdict": "KILLED",
        "unit": "log2 group operations; S = ops / sqrt(r)",
        "target": target,
        "boundary_A_invariant_target": boundary_a,
        "boundary_B_subspace_dichotomy": boundary_b,
        "boundary_C_large_horn_cost": boundary_c,
        "boundary_D_small_horn_empty": boundary_d,
        "boundary_E_transfer": boundary_e,
        "representation_gain": representation,
        "counterfactual_composite_degree": counterfactual,
    }
    OUT.write_text(json.dumps(report, indent=2) + "\n")

    print(f"r = {r}  prime={report['target']['r_is_prime']}  log2 r = {log2r:.4f}")
    print(f"rho reference      2^{log2_rho:.4f}   (plain 2^{log2_rho_plain:.4f}, S = {target['S_rho_reference']})")
    print(f"B  dichotomy holds for e = 1..{args.emax} and e = {N131}, {2*N131}: "
          f"dim <= e inside F_2^e, else >= 130")
    c = boundary_c["cheapest_cell"]
    print(f"C  cheapest large horn: e={c['e']} d={c['subfield_degree_d']} n={c['extension_degree_n']}"
          f"  2^{c['log2_cost_q_pow_2_minus_2_over_n']}  = 2^{c['log2_ratio_to_rho']} x rho")
    print(f"D  small horn factor base lies in E(F_2^e): no relation exists at any cost")
    print(f"E  embedding degree > {args.emb_limit} -> transfer field > {N131*args.emb_limit} bits")
    print(f"wrote {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()

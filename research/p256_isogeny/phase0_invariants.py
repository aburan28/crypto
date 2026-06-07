#!/usr/bin/env python3
"""
Phase 0 - fix invariants & sanity checks for P-256.

Per the search plan (Phase 0): compute and store the same-order invariants
(trace, Frobenius discriminant, small-prime valuations, embedding-degree
checks, twist order, CM field) and run the immediate-rejection classical
weak-curve checks.  All of these are *same-order invariants*: if P-256 does
not have them, no same-order neighbour over F_p can acquire them.

Runs in pure Python (no Sage). Self-validates arithmetic via [n]G == O.
"""

from __future__ import annotations
import json
import math
import sys

from ecfp import Curve, INF, is_square, sqrt_mod
import params as PP


def factor_small(m: int, bound: int = 100000):
    """Trial-divide out primes <= bound. Returns (small_factors, cofactor)."""
    m = abs(m)
    facs = {}
    d = 2
    while d <= bound and d * d <= m:
        while m % d == 0:
            facs[d] = facs.get(d, 0) + 1
            m //= d
        d += 1 if d == 2 else 2
    return facs, m  # cofactor m may still be composite/large


def val(m: int, ell: int) -> int:
    v = 0
    m = abs(m)
    while m % ell == 0:
        m //= ell
        v += 1
    return v


def is_probable_prime(n: int) -> bool:
    if n < 2:
        return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % sp == 0:
            return n == sp
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def main():
    p, a, b, n = PP.P, PP.A, PP.B, PP.N
    E = Curve(a, b, p)

    # ---- self-validation: G on curve and [n]G == O -----------------------
    G = (PP.GX, PP.GY)
    assert E.is_on_curve(G), "base point not on curve"
    assert E.mul(n, G) is INF, "[n]G != O -- constants wrong!"

    t = p + 1 - n                      # Frobenius trace
    disc = t * t - 4 * p               # Frobenius discriminant (< 0)
    L = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]

    # twist order
    twist_n = 2 * p + 2 - n            # = p + 1 + t
    Etw = E.quadratic_twist()

    # CM field K = Q(sqrt(disc)); fundamental disc via small-square stripping
    sq = 1
    d2 = -disc
    for q in range(2, 100000):
        while d2 % (q * q) == 0:
            d2 //= q * q
            sq *= q
    # disc = -(sq^2) * d2  (partial: only square factors with root < 1e5 removed)
    cm_disc_sqfree_part = -d2

    # embedding degree: smallest k with n | p^k - 1  (MOV/Frey-Ruck)
    emb_k = None
    acc = p % n
    for k in range(1, 200):
        if acc == 1:
            emb_k = k
            break
        acc = acc * (p % n) % n

    small_facs_disc = {ell: val(disc, ell) for ell in L}
    disc_small, disc_cofactor = factor_small(disc)

    rep = {
        "curve": "NIST P-256 / secp256r1",
        "p_bits": p.bit_length(),
        "a": a,
        "b": b,
        "n": n,
        "n_is_prime": is_probable_prime(n),
        "cofactor": 1,
        "j_invariant": E.j_invariant(),
        "discriminant_curve": E.disc(),
        "trace_t": t,
        "trace_bits": t.bit_length(),
        "frobenius_disc": disc,
        "frobenius_disc_bits": disc.bit_length(),
        "frobenius_disc_small_factors": {str(k): v for k, v in disc_small.items()},
        "frobenius_disc_large_cofactor_is_prime": is_probable_prime(disc_cofactor),
        "v_ell_disc": small_facs_disc,
        "cm_field": f"Q(sqrt({cm_disc_sqfree_part}))  [square part partially stripped]",
        "twist_order": twist_n,
        "twist_order_prime": is_probable_prime(twist_n),
        "twist_a": Etw.a,
        "twist_b": Etw.b,
        "embedding_degree_k": emb_k,
        # ---- classical weak-curve rejection checks (same-order invariants) --
        "reject_anomalous_n_eq_p": (n == p),
        "reject_supersingular_t_eq_0_mod_p": (t % p == 0),
        "reject_smooth_order": (not is_probable_prime(n)),
        "reject_small_embedding_degree": (emb_k is not None and emb_k < 100),
        "self_check_nG_is_identity": True,
    }

    # twist smoothness summary
    tw_small, tw_cof = factor_small(twist_n)
    rep["twist_largest_known_prime_factor_is_cofactor_prime"] = is_probable_prime(tw_cof)
    rep["twist_small_factors"] = {str(k): v for k, v in tw_small.items()}

    print(json.dumps(rep, indent=2))

    # human-readable verdict
    print("\n--- Phase 0 verdict (classical weak-curve buckets) ---", file=sys.stderr)
    verdicts = [
        ("anomalous (n == p)", rep["reject_anomalous_n_eq_p"]),
        ("supersingular (t == 0 mod p)", rep["reject_supersingular_t_eq_0_mod_p"]),
        ("smooth order (n composite)", rep["reject_smooth_order"]),
        ("small embedding degree (MOV/FR, k<100)", rep["reject_small_embedding_degree"]),
    ]
    for name, flag in verdicts:
        print(f"  {name:42s}: {'WEAK' if flag else 'OK (ruled out)'}", file=sys.stderr)
    print("  -> all classical buckets are same-order invariants; ruled out for", file=sys.stderr)
    print("     every same-order neighbour of P-256 over F_p.", file=sys.stderr)


if __name__ == "__main__":
    main()

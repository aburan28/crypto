#!/usr/bin/env python3
"""
Self-tests for the pure-Python P-256 isogeny toolkit.

Includes an *independent* validation of the Velu codomain formulas against
brute-force point counting at toy scale, plus the P-256 Phase 0 / Phase 1
invariants.  Run:  python3 test_p256_isogeny.py
"""

from __future__ import annotations
import sys

from ecfp import (Curve, INF, two_torsion_x, velu_2_isogeny,
                  velu_from_kernel_xs, three_isogenous_neighbors,
                  isogenous_neighbors, poly_roots_mod_p, order_is)
import params as PP
import semaev


def brute_order(E: Curve) -> int:
    """#E(F_p) by exhaustive count (toy p only)."""
    p, a, b = E.p, E.a, E.b
    cnt = 1  # point at infinity
    for x in range(p):
        rhs = (x * x * x + a * x + b) % p
        if rhs == 0:
            cnt += 1
        elif pow(rhs, (p - 1) // 2, p) == 1:
            cnt += 2
    return cnt


def test_poly_roots():
    p = 1009
    # (x-3)(x-7)(x-100) = roots {3,7,100}
    import random
    rng = random.Random(1)

    def polymul(A, B):
        r = [0] * (len(A) + len(B) - 1)
        for i, ai in enumerate(A):
            for j, bj in enumerate(B):
                r[i + j] = (r[i + j] + ai * bj) % p
        return r
    f = [1]
    for r in (3, 7, 100):
        f = polymul(f, [(-r) % p, 1])
    assert sorted(poly_roots_mod_p(f, p)) == [3, 7, 100]
    print("  [ok] poly_roots_mod_p")


def test_velu_brute():
    """Velu codomains preserve order, checked against brute force at toy scale."""
    checked = 0
    for p in (101, 1009, 2003):
        for a in range(-3, 4):
            for b in range(1, 6):
                if (4 * a ** 3 + 27 * b ** 2) % p == 0:
                    continue
                E = Curve(a, b, p)
                nE = brute_order(E)
                # 2-isogenies
                for x0 in two_torsion_x(E):
                    E2 = velu_2_isogeny(E, x0)
                    assert brute_order(E2) == nE, (p, a, b, "2-iso order mismatch")
                    checked += 1
                # 3-isogenies
                nbrs, roots = three_isogenous_neighbors(E)
                for E3 in nbrs:
                    assert brute_order(E3) == nE, (p, a, b, "3-iso order mismatch")
                    checked += 1
    print(f"  [ok] velu_*_isogeny preserve order (brute-force, {checked} codomains)")


def test_velu_oddl_brute():
    """Schoof-kernel ell-isogenies (ell=5,7,11,13) preserve order, brute-checked."""
    checked = 0
    for p in (101, 1009, 2003):
        for a in range(0, 4):
            for b in range(1, 4):
                if (4 * a ** 3 + 27 * b ** 2) % p == 0:
                    continue
                E = Curve(a, b, p)
                nE = brute_order(E)
                t = p + 1 - nE
                for ell in (5, 7, 11, 13):
                    if ell * ell > 4 * p:
                        continue
                    for lam, h, E2 in isogenous_neighbors(E, ell, t):
                        assert len(h) - 1 == (ell - 1) // 2, "kernel poly degree"
                        assert brute_order(E2) == nE, (p, a, b, ell, "order")
                        checked += 1
    print(f"  [ok] odd-ell Schoof isogenies preserve order ({checked} codomains)")


def test_semaev3_vanishes():
    """S_3(x1,x2,x3) = 0 whenever P1+P2+P3 = O on a toy curve."""
    import random
    p, a, b = 1009, 2, 3
    E = Curve(a, b, p)
    S3 = semaev.semaev3(a, b, p)
    rng = random.Random(5)
    checked = 0
    for _ in range(300):
        P1 = E.random_point(rng)
        P2 = E.random_point(rng)
        S = E.add(P1, P2)
        if S is INF:
            continue
        x1, x2, x3 = P1[0], P2[0], S[0]
        val = sum(c * pow(x1, e1, p) * pow(x2, e2, p) * pow(x3, e3, p)
                  for (e1, e2, e3), c in S3.items()) % p
        assert val == 0, "S_3 did not vanish on a real collinear triple"
        checked += 1
    print(f"  [ok] Semaev S_3 vanishes on {checked} real triples")


def test_p256_phase0():
    p, a, b, n = PP.P, PP.A, PP.B, PP.N
    E = Curve(a, b, p)
    G = (PP.GX, PP.GY)
    assert E.is_on_curve(G)
    assert E.mul(n, G) is INF, "[n]G != O"
    t = p + 1 - n
    assert n != p                      # not anomalous
    assert t % p != 0                  # not supersingular
    print("  [ok] P-256 Phase 0 (on-curve, [n]G=O, not anomalous/supersingular)")


def test_p256_phase1():
    p, n = PP.P, PP.N
    E = Curve(PP.A, PP.B, p)
    assert len(two_torsion_x(E)) == 0, "P-256 should have no rational 2-torsion"
    nbrs, roots = three_isogenous_neighbors(E)
    assert len(roots) == 1, "P-256 has exactly one rational 3-isogeny (ramified)"
    E3 = nbrs[0]
    assert E3.j_invariant() != E.j_invariant(), "neighbour must differ from root"
    assert order_is(E3, n, p, trials=10, seed=7), "3-isogenous neighbour same order"
    print("  [ok] P-256 Phase 1 (0 rational 2-iso, 1 verified same-order 3-iso)")


def main():
    print("running P-256 isogeny toolkit self-tests...")
    test_poly_roots()
    test_velu_brute()
    test_velu_oddl_brute()
    test_semaev3_vanishes()
    test_p256_phase0()
    test_p256_phase1()
    print("ALL TESTS PASSED")


if __name__ == "__main__":
    sys.exit(main())

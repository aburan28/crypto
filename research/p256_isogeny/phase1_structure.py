#!/usr/bin/env python3
"""
Phase 1 - structural map of the small-degree isogeny neighbourhood of P-256.

For each prime ell in L, the rational ell-isogeny structure of an ordinary
curve is determined *exactly* by the Frobenius discriminant Delta_pi = t^2-4p
(a same-order invariant), via the Deuring/Kohel volcano picture:

    * ell ramified  (ell | Delta_pi):   ell | conductor or ell | D_K.
          - if v_ell(Delta_pi) is odd -> ell | D_K (ramified in K): 1 rational
            ell-isogeny (vertical, to/from the rest of the volcano).
    * ell split     ((Delta_pi/ell)=+1, ell ∤ Delta_pi):  ell splits in K.
          crater is a cycle -> exactly 2 horizontal rational ell-isogenies.
    * ell inert     ((Delta_pi/ell)=-1):  ell inert in K.
          0 horizontal isogenies; on a height-0 volcano -> 0 rational
          ell-isogenies at all.

The number of roots of the modular polynomial Phi_ell(j_E, Y) over F_p equals
the number of rational ell-isogenies, so this predicts the crawl degree without
needing Phi_ell.  We then *confirm* the ell=2 and ell=3 predictions by direct
construction (2-division roots; psi_3 roots + Velu).
"""

from __future__ import annotations
import json

from ecfp import Curve, two_torsion_x, three_isogenous_neighbors, order_is
import params as PP


def legendre(a: int, ell: int) -> int:
    a %= ell
    if a == 0:
        return 0
    ls = pow(a, (ell - 1) // 2, ell)
    return -1 if ls == ell - 1 else 1


def v_ell(m: int, ell: int) -> int:
    v, m = 0, abs(m)
    while m and m % ell == 0:
        m //= ell
        v += 1
    return v


def cheap_score(E: Curve):
    """Phase 4 cheap-invariant features for a candidate curve model."""
    p = E.p
    j = E.j_invariant()
    # automorphism group size: generic 2, j=1728 -> 4, j=0 -> 6
    if j == 0:
        aut = 6
    elif j == 1728 % p:
        aut = 4
    else:
        aut = 2
    return {
        "j_is_0": j == 0,
        "j_is_1728": j == 1728 % p,
        "a_is_minus3": E.a == (-3) % p,
        "a_hamming_weight": bin(E.a).count("1"),
        "b_hamming_weight": bin(E.b).count("1"),
        "automorphism_group_size": aut,
        "special": (j in (0, 1728 % p)) or E.a == (-3) % p,
    }


def classify(ell: int, disc: int):
    v = v_ell(disc, ell)
    if v > 0:
        kind = "ramified"
        # v odd  -> ell ramifies in K (1 vertical isogeny)
        # v even -> ell | conductor only; on height v/2 volcano
        predicted = 1 if v % 2 == 1 else (ell + 1)
        height = v // 2
    else:
        leg = legendre(disc, ell)
        if leg == 1:
            kind = "split"
            predicted = 2          # crater cycle: 2 horizontal edges
            height = 0
        else:
            kind = "inert"
            predicted = 0
            height = 0
    return {"ell": ell, "kind": kind, "v_ell_disc": v,
            "legendre_disc_mod_ell": legendre(disc, ell) if v == 0 else 0,
            "predicted_rational_isogenies": predicted, "volcano_height": height}


def main():
    p, n = PP.P, PP.N
    t = p + 1 - n
    disc = t * t - 4 * p
    E0 = Curve(PP.A, PP.B, p)

    L = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    rows = [classify(ell, disc) for ell in L]

    # ---- empirical confirmation for ell = 2 and ell = 3 ------------------
    roots2 = two_torsion_x(E0)
    nbr3, roots3 = three_isogenous_neighbors(E0)

    confirm = {
        "ell2_rational_2torsion_roots": len(roots2),
        "ell2_matches_prediction": len(roots2) == rows[0]["predicted_rational_isogenies"],
        "ell3_psi3_rational_roots": len(roots3),
        "ell3_matches_prediction": len(roots3) == rows[1]["predicted_rational_isogenies"],
    }

    # construct + verify the real 3-isogenous neighbour(s)
    j0 = E0.j_invariant()
    neighbors3 = []
    for E3, xk in zip(nbr3, roots3):
        same = order_is(E3, n, p, trials=8, seed=999)
        j3 = E3.j_invariant()
        neighbors3.append({
            "kernel_x": xk,
            "neighbor_a": E3.a,
            "neighbor_b": E3.b,
            "neighbor_j": j3,
            "same_order_as_p256": same,
            "distinct_from_root": j3 != j0,
            "cheap_invariants": cheap_score(E3),
        })

    out = {
        "frobenius_disc": disc,
        "structural_prediction": rows,
        "split_primes": [r["ell"] for r in rows if r["kind"] == "split"],
        "inert_primes": [r["ell"] for r in rows if r["kind"] == "inert"],
        "ramified_primes": [r["ell"] for r in rows if r["kind"] == "ramified"],
        "empirical_confirmation": confirm,
        "constructed_3_isogenous_neighbors": neighbors3,
    }
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()

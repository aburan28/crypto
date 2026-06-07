#!/usr/bin/env python3
"""
Phase 5 - endomorphism / GLV audit for P-256 (hypothesis H1).

A GLV speedup needs an *efficiently computable* endomorphism psi of small norm
whose eigenvalue lambda mod n yields a short scalar decomposition.  We quantify
whether P-256 (or any of its same-order neighbours) has one.

Endomorphism lattice.  For an ordinary curve, End(E) is an order in the CM field
K = Q(sqrt(Delta_pi)).  Every endomorphism a + b*pi (pi = Frobenius, pi^2 - t pi
+ p = 0) has norm

    N(a + b pi) = a^2 + a b t + b^2 p   = the binary quadratic form [1, t, p].

Reducing this form gives the minimum nonzero norm achievable with b != 0, i.e.
the degree of the smallest *non-scalar* endomorphism.  If that degree is ~p
(2^256-ish) the only cheap endomorphism is Frobenius itself, which acts as the
identity on the rational order-n subgroup (lambda = 1, useless for GLV).

Outputs the minimal non-scalar endomorphism norm, automorphism sizes for the
root and every catalogued neighbour, a self-isogeny-loop check, and the GLV
bits-saved estimate against the acceptance thresholds in the plan.
"""

from __future__ import annotations
import json
import math
import os

from ecfp import Curve
import params as PP


def reduce_form_exact(A, B, disc):
    """
    Gauss reduction of the positive-definite form [A, B, *] of discriminant
    disc < 0.  C is always derived from the invariant disc = B^2 - 4AC, which
    avoids arithmetic drift and the B = -A boundary loop.
    """
    C = (B * B - disc) // (4 * A)
    while True:
        # normalize B into (-A, A]
        B %= (2 * A)
        if B > A:
            B -= 2 * A
        C = (B * B - disc) // (4 * A)
        if A <= C:
            break
        A, B = C, -B            # swap (A,B,C) -> (C,-B,A)
    if A == C and B < 0:        # canonical sign on the boundary
        B = -B
    C = (B * B - disc) // (4 * A)
    return A, B, C


def main():
    p, n = PP.P, PP.N
    t = p + 1 - n
    disc = t * t - 4 * p

    # reduce the norm form [1, t, p]  (discriminant = disc)
    A, B, C = reduce_form_exact(1, t, disc)
    assert B * B - 4 * A * C == disc, "form reduction changed discriminant"

    # minimal non-scalar endomorphism norm = smallest value with y != 0.
    # In a reduced form [A,B,C] (A<=C), values at (0,1)->C, (1,1)->A+B+C, etc.;
    # the smallest with y!=0 is C (since A is the (1,0) scalar value = 1 here).
    min_nonscalar_norm = C
    frob_norm = p  # N(pi) = p

    E0 = Curve(PP.A, PP.B, p)
    j0 = E0.j_invariant()

    def aut_size(j):
        if j == 0:
            return 6
        if j == 1728 % p:
            return 4
        return 2

    # automorphism audit of root + catalogued neighbours
    neigh_path = os.path.join(os.path.dirname(__file__), "neighbors.jsonl")
    neighbor_auts = []
    self_loop = False
    if os.path.exists(neigh_path):
        with open(neigh_path) as f:
            for line in f:
                r = json.loads(line)
                neighbor_auts.append({"ell": r["ell"], "j_mod_1e9": r["j"] % 10**9,
                                      "aut": aut_size(r["j"]),
                                      "is_self_loop": r["j"] == j0})
                if r["j"] == j0:
                    self_loop = True

    # GLV bits saved:
    #  - automorphism group beyond the trivial {+-1} (size 2) would add
    #    0.5*log2(|Aut|/2) bits; root aut = 2 -> 0 extra bits.
    #  - the only small-norm (efficiently computable) endomorphism is Frobenius,
    #    which acts as identity on <G> -> lambda = 1 -> 0 bits.
    root_aut = aut_size(j0)
    aut_bits_saved = 0.5 * math.log2(root_aut / 2) if root_aut > 2 else 0.0

    out = {
        "curve": "NIST P-256",
        "frobenius_disc_bits": disc.bit_length(),
        "reduced_norm_form": [A, B, C],
        "frobenius_norm_log2": round(math.log2(frob_norm), 2),
        "min_nonscalar_endomorphism_norm_log2": round(math.log2(min_nonscalar_norm), 2),
        "interpretation": (
            "smallest non-scalar endomorphism has degree ~2^%d; only Frobenius "
            "(deg p) is efficiently computable and it acts as identity on <G>"
            % round(math.log2(min_nonscalar_norm))
        ),
        "root_automorphism_group_size": root_aut,
        "any_neighbor_special_automorphism": any(a["aut"] > 2 for a in neighbor_auts),
        "self_isogeny_loop_found": self_loop,
        "neighbor_automorphism_audit": neighbor_auts,
        "glv": {
            "efficient_nontrivial_endomorphism": False,
            "frobenius_eigenvalue_on_subgroup": 1,
            "automorphism_bits_saved": aut_bits_saved,
            "total_rho_bits_saved_estimate": aut_bits_saved,
        },
        "thresholds": {
            "interesting_bits": 1, "serious_bits": 8, "breakthrough_bits": 32,
            "meets_interesting": aut_bits_saved >= 1,
        },
        "verdict": (
            "H1 NOT supported: no efficiently computable non-trivial endomorphism; "
            "no special-j neighbour; rho bits saved ~ 0 (< interesting threshold 1)."
        ),
    }
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()

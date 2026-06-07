#!/usr/bin/env python3
"""
Phase 1 (Mode C/D) - split-prime crater walk around P-256.

Constructs honest same-order neighbours of P-256 via rational ell-isogenies for
odd ell, using the Schoof-style kernel (Frobenius eigenvalue -> kernel polynomial
-> Velu codomain from symmetric functions).  Two products:

  1. CATALOG: one step along every ramified/split prime ell <= 29, each codomain
     order-verified and cheap-scored.
  2. HORIZONTAL WALK: for the small split primes, walk the crater cycle several
     steps without backtracking, confirming every node is a distinct same-order
     neighbour.

Frobenius trace t = p+1-n is a same-order invariant, so it is constant along the
whole walk (every node has order n); we pass it once.

Writes neighbors.jsonl (catalog) and walk_<ell>.jsonl (walks).
"""

from __future__ import annotations
import json
import sys
import time

from ecfp import Curve, isogenous_neighbors, order_is
from analysis import cheap_score
import params as PP

P, N = PP.P, PP.N
T = P + 1 - N


def record(E: Curve, extra=None):
    r = {
        "a": E.a, "b": E.b, "j": E.j_invariant(),
        "cheap_invariants": cheap_score(E),
    }
    if extra:
        r.update(extra)
    return r


def build_catalog(primes):
    E0 = Curve(PP.A, PP.B, P)
    j0 = E0.j_invariant()
    rows = []
    for ell in primes:
        t0 = time.time()
        nbrs = isogenous_neighbors(E0, ell, T)
        dt = time.time() - t0
        for lam, h, E2 in nbrs:
            same = order_is(E2, N, P, trials=8, seed=2024)
            rows.append(record(E2, {
                "ell": ell, "frobenius_eigenvalue_mod_ell": lam,
                "kernel_poly_degree": len(h) - 1,
                "same_order_as_p256": same,
                "distinct_from_root": E2.j_invariant() != j0,
                "construct_seconds": round(dt, 3),
            }))
        print(f"  ell={ell:2d}: {len(nbrs)} neighbour(s)  ({dt:.2f}s)", file=sys.stderr)
    return rows


def horizontal_walk(ell, steps):
    """Walk the ell-crater forward (no backtracking) from P-256."""
    E = Curve(PP.A, PP.B, P)
    prev_j = None
    path = [record(E, {"step": 0})]
    for s in range(1, steps + 1):
        nbrs = isogenous_neighbors(E, ell, T)
        forward = [E2 for _, _, E2 in nbrs if E2.j_invariant() != prev_j]
        if not forward:
            break
        prev_j = E.j_invariant()
        E = forward[0]
        same = order_is(E, N, P, trials=6, seed=100 + s)
        path.append(record(E, {"step": s, "same_order_as_p256": same}))
    return path


def main():
    catalog_primes = [3, 5, 11, 13, 17, 23, 29]
    print("building catalog (one step per ramified/split prime)...", file=sys.stderr)
    catalog = build_catalog(catalog_primes)
    with open("neighbors.jsonl", "w") as f:
        for r in catalog:
            f.write(json.dumps(r) + "\n")

    walks = {}
    for ell, steps in [(11, 6), (13, 4)]:
        print(f"horizontal walk ell={ell}, {steps} steps...", file=sys.stderr)
        w = horizontal_walk(ell, steps)
        walks[ell] = w
        with open(f"walk_{ell}.jsonl", "w") as f:
            for r in w:
                f.write(json.dumps(r) + "\n")

    summary = {
        "root_j": Curve(PP.A, PP.B, P).j_invariant(),
        "catalog_size": len(catalog),
        "catalog_all_same_order": all(r["same_order_as_p256"] for r in catalog),
        "catalog_all_distinct_from_root": all(r["distinct_from_root"] for r in catalog),
        "catalog_any_special_cheap_score": any(r["cheap_invariants"]["special"] for r in catalog),
        "neighbor_count_by_ell": {
            str(ell): sum(1 for r in catalog if r["ell"] == ell) for ell in catalog_primes
        },
        "walks": {
            str(ell): {
                "length": len(w),
                "distinct_j": len({r["j"] for r in w}),
                "all_same_order": all(r.get("same_order_as_p256", True) for r in w),
                "any_special_cheap_score": any(r["cheap_invariants"]["special"] for r in w),
                "j_sequence_tail6": [r["j"] % 10**12 for r in w][:6],
            } for ell, w in walks.items()
        },
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

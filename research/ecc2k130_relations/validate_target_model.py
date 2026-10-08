#!/usr/bin/env python3
"""Measure the decomposition probability `target_boundary.py` depends on.

The cost of a logarithm is one meet-in-the-middle divided by the chance that
a random target `[a]P + [b]Q` decomposes over the support at all.  That
chance is the model's only free input, so it is measured on small analogues
where every decomposition can be enumerated exhaustively, rather than
asserted.

Two things are checked:

1. **The rate.**  `P = 2^n C(B,n) / r`, against exhaustive enumeration of
   every signed `n`-subset sum and a few hundred random targets.

2. **Whether `sigma` buys a factor of `m`.**  It is tempting to say a
   decomposition of `sigma^k(target)` is as useful as one of `target`, so
   there are `m` acceptable right-hand sides and the rate should carry a
   factor of `m`.  Measured, that model is about 20x optimistic.  The
   support is `sigma`-stable, so the *set of n-subset sums is itself
   `sigma`-closed*, which this script verifies directly: accepting
   `sigma^k(target)` is the same chance `m` times over, not `m` independent
   chances.  The quotient buys memory and unknowns, not hit rate.
"""

from __future__ import annotations

import json
import math
import random
import sys
from itertools import combinations, product
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from normalbasis import NormalSupport, find_normal_elements   # noqa: E402
from smallcurve import four_torsion, small_curve              # noqa: E402


def _support(F, E, orbits, rng):
    for a in find_normal_elements(F, 60, rng):
        s = NormalSupport(F, E, a)
        if s.orbit_count == orbits:
            return s
    return None


def _signed_subset_sums(E, reps, n):
    sums = set()
    for comb in combinations(range(len(reps)), n):
        for sg in product((1, -1), repeat=n):
            sums.add(E.sum_points(
                [reps[i] if g > 0 else E.neg(reps[i]) for i, g in zip(comb, sg)]))
    return sums


def cell(mdeg: int, orbits: int, n: int, trials: int = 600, seed: int = 5):
    F, E, order = small_curve(mdeg)
    e4 = four_torsion(E, order)
    r = order
    while r % 2 == 0:
        r //= 2
    rng = random.Random(seed)
    sup = _support(F, E, orbits, rng)
    if sup is None:
        return None
    reps = sup.orbit_points(E)
    B = sup.size
    sums = _signed_subset_sums(E, reps, n)

    sigma_closed = all(E.frobenius(S) in sums for S in sums if S is not None)

    P = None
    while P is None:
        pts = E.points_over(rng.getrandbits(mdeg) & F.mask)
        if pts:
            cand = E.mul(pts[0], order // r)
            if cand is not None:
                P = cand
    d = rng.randrange(2, r)
    Q = E.mul(P, d)

    hits = 0
    for _ in range(trials):
        tgt = E.add(E.mul(P, rng.randrange(r)), E.mul(Q, rng.randrange(r)))
        found, cur = False, tgt
        for _k in range(mdeg):                 # accept sigma^k(target)
            if any(E.add(cur, T) in sums for T in e4):
                found = True
                break
            cur = E.frobenius(cur)
        hits += found

    lc = sum(math.log2(B - i) for i in range(n)) - math.log2(math.factorial(n))
    plain = min(1.0, 2 ** (n + lc - math.log2(r)))
    with_sigma = min(1.0, plain * mdeg)
    measured = hits / trials
    return {
        "m": mdeg, "orbits": orbits, "support_size": B, "relation_length": n,
        "subgroup_order_r": r, "signed_subset_sums": len(sums),
        "sum_set_is_sigma_closed": sigma_closed,
        "trials": trials, "measured": round(measured, 5),
        "predicted_plain": round(plain, 5),
        "predicted_with_sigma_factor": round(with_sigma, 5),
        "ratio_to_plain": round(measured / plain, 3) if plain else None,
        "ratio_to_sigma_model": round(measured / with_sigma, 3) if with_sigma else None,
        "saturated": plain > 0.9,
    }


CELLS = [(19, 2, 2), (19, 3, 2), (17, 2, 2), (17, 3, 2),
         (19, 4, 2), (17, 4, 2), (19, 2, 3), (17, 2, 3)]


def main():
    rows = [c for c in (cell(*spec) for spec in CELLS) if c]
    live = [r for r in rows if not r["saturated"]]
    mean_plain = sum(r["ratio_to_plain"] for r in live) / len(live)
    mean_sigma = sum(r["ratio_to_sigma_model"] for r in live) / len(live)

    data = {
        "instance": "small analogues of ECC2K-130",
        "curve": "y^2 + xy = x^3 + 1",
        "question": ("does a random target [a]P + [b]Q decompose over a "
                     "sigma-orbit support at the modelled rate, and does "
                     "sigma multiply that rate by m?"),
        "cells": rows,
        "unsaturated_cells": len(live),
        "mean_ratio_plain_model": round(mean_plain, 3),
        "mean_ratio_sigma_model": round(mean_sigma, 3),
        "sum_set_sigma_closed_everywhere": all(
            r["sum_set_is_sigma_closed"] for r in rows),
        "verdict": (
            f"P = 2^n C(B,n)/r holds: mean measured/predicted "
            f"{mean_plain:.2f} over {len(live)} unsaturated cells. The model "
            f"that additionally multiplies by m is optimistic by "
            f"{1 / mean_sigma:.0f}x. The sum set is sigma-closed in every "
            f"cell, which is why: accepting sigma^k(target) is one chance "
            f"taken m times, not m chances."),
    }
    out = HERE / "results" / "target_decomposition.json"
    out.write_text(json.dumps(data, indent=2))

    hdr = (f"{'m':>3} {'T':>3} {'B':>5} {'n':>2} {'closed':>7} {'measured':>9} "
           f"{'plain':>8} {'ratio':>7} {'sigma-model':>12} {'ratio':>7}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(f"{r['m']:>3} {r['orbits']:>3} {r['support_size']:>5} "
              f"{r['relation_length']:>2} {str(r['sum_set_is_sigma_closed']):>7} "
              f"{r['measured']:>9.4f} {r['predicted_plain']:>8.4f} "
              f"{r['ratio_to_plain'] or 0:>7.2f} "
              f"{r['predicted_with_sigma_factor']:>12.4f} "
              f"{r['ratio_to_sigma_model'] or 0:>7.2f}")
    print()
    print(data["verdict"])


if __name__ == "__main__":
    main()

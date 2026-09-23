#!/usr/bin/env python3
"""The triple collector's base size, from a count of four-sums rather than W₃.

    python3 research/ic_triple_counted_20260923/model.py

The triple arm (`research/ic_triple_table_20260923`) sizes its base with W₃, the
pair model's units carried over.  W₃'s hit rate is `λ₃·cov₃` a rest,
`λ₃ = C(|F|+2, 3)/r`, so a target of `|F|` rests is modelled as meeting
`|F|·C(|F|+2, 3)/r` witnesses -- about four times what exists: a target is found
exactly when it is a sum of four base points with at least one summand in a
built row, and there are only `C(|F|+3, 4) − C(|F|−2nt+3, 4)` such multisets.

This file writes the count down, checks it against the triple arm's committed
trial counts (it was not fitted to them), and prints what it and W₃ choose.
Nothing here is a measurement of the counted arm.
"""
import json
import math
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
TRIPLE = HERE.parent / "ic_triple_table_20260923"
CELLS = {"n23a1": (23, 4_196_903), "n37a0": (37, 230_603_167), "n43a1": (43, 4_644_189_029),
         "n59a0": (59, 10_063_074_221), "n61a1": (61, 11_514_943_771)}
PROBES = 3.0  # TABLE_PROBES in the solver


def limit(n):
    """The packed entry's ceiling on K, as `triple_orbits` has it."""
    return max(min(65535 // (2 * n), 255, 64), 2)


# ── W₃, exactly as the committed triple arm evaluates it ────────────────────
def w3(n, r, K, t):
    s = 2 * n * K
    lam = s * (s + 1) * (s + 2) / 6 / r
    return t * s * (s + 1) / 2 + (K + PROBES) / max(lam * (1 - ((K - t) / K) ** 3), 1e-12)


def w3_choice(n, r):
    rows = lambda K: min(range(1, K + 1), key=lambda t: w3(n, r, K, t))
    K = min(range(2, limit(n) + 1), key=lambda K: w3(n, r, K, rows(K)))
    return K, rows(K)


# ── the count ───────────────────────────────────────────────────────────────
def p_target(n, r, K, t):
    """A fresh target is a four-sum of base points with a summand in the first t orbits."""
    s = 2 * n * K
    return -math.expm1(-(comb(s + 3, 4) - comb(s - 2 * n * t + 3, 4)) / r)


def units(n, r, K, t, hits=None):
    """Table entries plus rests scanned: K relations and one descent, |F| rests a target."""
    s = 2 * n * K
    hits = K + 1 if hits is None else hits
    return t * s * (s + 1) / 2 + hits * s / p_target(n, r, K, t)


def counted_choice(n, r):
    rows = lambda K: min(range(1, K + 1), key=lambda t: units(n, r, K, t))
    K = min(range(2, limit(n) + 1), key=lambda K: units(n, r, K, rows(K)))
    return K, rows(K)


def main():
    confirm = json.loads((TRIPLE / "confirm.json").read_text())
    rho = json.loads((TRIPLE / "rho-same-fixtures.json").read_text())
    gm = lambda v: math.exp(sum(map(math.log, v)) / len(v))

    print("1. The count against the triple arm's committed trials (32 fixtures a cell, not fitted):")
    for cell in ("n23a1", "n37a0", "n43a1"):
        n, r = CELLS[cell]
        K, t = w3_choice(n, r)
        rows = [q for q in confirm if q["cell"] == cell]
        mean = sum(q["trials_triple"] for q in rows) / len(rows)
        print(f"   {cell}: K={K} t={t}  predicted mean trials {(K + 1) / p_target(n, r, K, t):6.1f}   measured {mean:6.1f}")

    print("\n2. What each sizing chooses, in the count's units:")
    for cell, (n, r) in CELLS.items():
        (K3, t3), (Kc, tc) = w3_choice(n, r), counted_choice(n, r)
        u3, uc = units(n, r, K3, t3), units(n, r, Kc, tc)
        extra = units(n, r, Kc, tc, hits=Kc + 2)
        print(f"   {cell}: W3 -> K={K3} t={t3} ({u3:>7.0f})   count -> K={Kc} t={tc} ({uc:>7.0f})   "
              f"ratio {uc / u3:.3f}, with one extra hit {extra / u3:.3f}")

    # The whole job at n43a1: the fixed phases (curve and targets, final verification,
    # certification, start-up; phase-profile-probe.json) are ~2.03M instructions and do
    # not depend on the base; everything else is priced at the triple arm's own
    # measured instructions per unit.
    n, r = CELLS["n43a1"]
    fixed = 2.03e6
    rows = [q for q in confirm if q["cell"] == "n43a1"]
    mean_ir = sum(q["ir_triple"] for q in rows) / len(rows)
    per_unit = (mean_ir - fixed) / units(n, r, *w3_choice(n, r))
    lo = (fixed + per_unit * units(n, r, 3, 1)) / mean_ir
    hi = (fixed + per_unit * units(n, r, 3, 1, hits=5)) / mean_ir
    print(f"\n3. n43a1 whole job, counted/triple: {lo:.3f} (K+1 hits) .. {hi:.3f} (one extra hit);"
          f" {per_unit:.0f} instructions a unit")

    # Holdouts, against rho: the triple arm's n43a1 ratio carried by the count's units
    # and by rho's expected steps, sqrt(r / n) under the order-2n automorphism group.
    base_ratio = gm([q["ir_triple"] / q["ir_rho"] for q in rho if q["cell"] == "n43a1"])
    u43 = units(n, r, *w3_choice(n, r))
    print(f"\n4. counted/rho, carried from n43a1's triple/rho = {base_ratio:.3f} (secondary in the last study):")
    print(f"   n43a1: {base_ratio * lo:.3f} .. {base_ratio * hi:.3f}")
    for cell in ("n59a0", "n61a1"):
        nn, rr = CELLS[cell]
        steps = math.sqrt(rr / nn) / math.sqrt(r / n)
        u = units(nn, rr, *counted_choice(nn, rr))
        print(f"   {cell}: units x{u / u43:.3f}, rho steps x{steps:.3f} -> {base_ratio * (u / u43) / steps:.3f}")


if __name__ == "__main__":
    main()

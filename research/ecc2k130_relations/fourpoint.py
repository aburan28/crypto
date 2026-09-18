#!/usr/bin/env python3
"""How many four-point relations a support of size `B` actually holds.

The counting floor of the research note is about *multisets*.  A search
does not enumerate multisets: it enumerates four distinct support abscissae
with a sign on each, modulo a global negation.  Those two counts differ,
and the difference is not a constant -- it is the whole cancelling family.

**What a search space contains.**  Four distinct abscissae out of `B`, each
contributing a point up to sign, with a global negation that changes
nothing a collision search can see:

    configurations = C(B, 4) * 2^4 / 2 = 8 * C(B, 4)

**Why the abscissae must be distinct.**  Allowing a repeat admits
`R + (-R) + T + (-T) = O`, which is a relation, holds for every `R` and
`T`, and says nothing whatever about any logarithm.  It is the
construction equation of the research note in its purest form.  A supply
model that counts these overstates the yield by orders of magnitude at
small `B` -- that is the 2x-43x miss this module exists to have caught --
and it overstates it *in the validation code*, where an inflated model
makes a correct search look broken.

**What counts as a relation.**  The cofactor is 4.  A sum landing anywhere
in `E[4]`, not only on the identity, gives a usable equation once the
cofactor is cleared.

**Where the sums actually live.**  Not in `E`.  Every point of a
normal-basis weight-two support lies in the index-2 subgroup `H` of order
`2r` -- measured on the challenge curve and on the small analogues, with no
exceptions -- so a four-point sum lies in `H` always, and

    E[4] cap H = {O, the order-2 point},    |E[4] cap H| = 2.

The two order-4 points are unreachable, and exhaustive counts at `m = 11,
13, 19` put exactly zero relations on them across every support tried. So:

    expected_usable = 8 * C(B, 4) * 2 / (2r) = 8 * C(B, 4) / r
    expected_exact  = 8 * C(B, 4) * 1 / (2r)

and the coset factor -- what accepting `E[4]` buys over demanding `O` -- is
exactly **2**.  Assuming the four cosets are hit uniformly gives 4 and is
wrong; `2.33` was an earlier measurement of the same quantity taken before
the counter below was checked against brute force.

Note the sign of this.  Exact relations are *twice* as common as a uniform
model over `E` predicts, which looks like structure beating the independence
assumption -- and the usable rate, the one that matters, lands exactly on
`8 C(B,4) / r`.  The bias moves relations between cosets; it does not make
any more of them.
"""

from __future__ import annotations

import math
from itertools import product


def configurations(B: int) -> int:
    """`8 * C(B, 4)`: four distinct abscissae, signed, modulo global negation."""
    if B < 4:
        return 0
    return 8 * math.comb(B, 4)


def supply_prediction(B: int, group_order: int, n_four_torsion: int = 4):
    """Expected four-point relations from a size-`B` support, both rates.

    `group_order` is `#E`; the sums live in its index-2 subgroup, and the
    reachable part of `E[4]` is the half of it that lies there.
    """
    cfg = configurations(B)
    host = group_order // 2                 # the index-2 subgroup
    reachable = max(1, n_four_torsion // 2)  # E[4] cap H
    usable = cfg * reachable / host
    exact = cfg / host
    return {
        "support_size": B,
        "configurations": cfg,
        "log2_configurations": round(math.log2(cfg), 3) if cfg else None,
        "group_order": group_order,
        "host_subgroup_order": host,
        "four_torsion_size": n_four_torsion,
        "reachable_four_torsion": reachable,
        "coset_factor": reachable,
        "expected_usable": usable,
        "expected_exact": exact,
        "log2_expected_usable": round(math.log2(usable), 3) if cfg else None,
        "log2_expected_exact": round(math.log2(exact), 3) if cfg else None,
    }


def exhaustive_four_point(E, abscissae, four_torsion):
    """Ground truth: every four-point relation the support actually holds.

    Counts configurations `{i<j<k<l}` with signs, modulo global negation,
    whose four points sum into `E[4]` (usable) or to `O` (exact).

    Enumeration is by pair sums rather than by quadruples.  Every pair
    `{k, l}` with both signs free goes into a table keyed by its sum; then
    for each pair `{i, j}` containing the least index and each target `T`
    in `E[4]`, the partners are `table[T - S]`.  A configuration is reached
    once for each of the three ways its least element can be paired off, so
    the raw tally is divided by three -- checked against a direct quadruple
    enumeration on the smallest supports in `test_four_point_orbit.py`.
    """
    reps = []
    for x in abscissae:
        pts = E.points_over(x)
        if pts:
            reps.append(min(pts, key=lambda p: p[1]))
    B = len(reps)

    table: dict = {}
    for k in range(B):
        for l in range(k + 1, B):
            for sk, sl in product((1, -1), repeat=2):
                S = E.add(reps[k] if sk > 0 else E.neg(reps[k]),
                          reps[l] if sl > 0 else E.neg(reps[l]))
                table.setdefault(S, []).append((k, l))

    usable = exact = 0
    for i in range(B):
        for j in range(i + 1, B):
            for sj in (1, -1):
                S1 = E.add(reps[i], reps[j] if sj > 0 else E.neg(reps[j]))
                for T in four_torsion:
                    # partners must sum to T - S1 = T + (-S1)
                    want = E.add(T, E.neg(S1)) if S1 is not None else T
                    for (k, l) in table.get(want, ()):
                        if k <= i or k == j or l == j:
                            continue
                        usable += 1
                        if T is None:
                            exact += 1
    assert usable % 3 == 0 and exact % 3 == 0, "triple-count invariant broken"
    return {
        "support_size": B,
        "usable_relations": usable // 3,
        "exact_relations": exact // 3,
    }

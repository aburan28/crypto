"""Pollard rho reference: the boundary every row is measured against.

An r-adding walk (Teske, r = 16) with the negation map, distinguished
points and a 2-cycle escape by doubling, run on the very same subgroup,
target and accounting as the auxiliary-input algorithms.  This is the
repository's reference algorithm for the plain ECDLP, `S ≈ 1.3`.

Nothing about the auxiliary inputs is visible to it: it sees `G` and
`Q = [α]G` only.  Each step is one affine addition (one group operation);
the r table points, restarts, the escape doublings and the final
verification are all charged.
"""

from __future__ import annotations

import random
from typing import Dict, Optional, Tuple

from ec import Curve, Point


def _canon(E: Curve, P: Point) -> Tuple[Point, int]:
    """Negation-map representative (smaller y) and the sign applied."""
    x, y = P
    if y > E.q - y:
        return (x, E.q - y), -1
    return P, 1


def pollard_rho(E: Curve, G: Point, Q: Point, p: int, seed: int, r: int = 32,
                dp_bits: Optional[int] = None, max_ops: Optional[int] = None,
                cycle_window: int = 64) -> Dict:
    """Solve [alpha]G = Q.  Returns a dict with 'alpha' (None if `max_ops`
    was exhausted), 'ops' (group operations, all phases) and diagnostics.

    Fruitless cycles (Bernstein–Lange–Schwabe): a 2-cycle is escaped by
    doubling; longer short cycles are caught by a window of the last
    `cycle_window` points and escaped the same way, from the point where
    the cycle was detected (a deterministic function of the walk, so two
    merged walks stay merged)."""
    rng = random.Random(seed)
    ctr = E.ctr
    start = ctr.ops
    if dp_bits is None:
        dp_bits = max(2, (p.bit_length() // 2) - 7)
    dp_mask = (1 << dp_bits) - 1
    walk_cap = 20 * (1 << dp_bits)
    escapes = 0

    ctr.begin("rho:precompute")
    M = []
    for _ in range(r):
        a, b = rng.randrange(1, p), rng.randrange(1, p)
        M.append((E.add(E.mul(a, G), E.mul(b, Q)), a, b))
    ctr.end()

    ctr.begin("rho:walk")
    table: Dict[Point, Tuple[int, int]] = {}
    steps = walks = 0
    while True:
        walks += 1
        a, b = rng.randrange(1, p), rng.randrange(1, p)
        W = E.add(E.mul(a, G), E.mul(b, Q))
        if W is None:
            continue
        W, s = _canon(E, W)
        a, b = a * s % p, b * s % p
        recent: Dict[Point, int] = {}   # point -> step index, last `cycle_window` steps
        order: list = []
        for _ in range(walk_cap):
            steps += 1
            if max_ops is not None and ctr.ops - start > max_ops:
                ctr.end()
                return {"alpha": None, "ops": ctr.ops - start, "steps": steps, "walks": walks}
            if W[0] & dp_mask == 0:
                hit = table.get(W)
                if hit is not None:
                    a2, b2 = hit
                    if (b - b2) % p == 0:
                        break  # useless collision: restart the walk
                    alpha = (a2 - a) * pow(b - b2, -1, p) % p
                    ctr.end()
                    ctr.begin("rho:verify")
                    ok = E.mul(alpha, G) == Q
                    ctr.end()
                    return {"alpha": alpha if ok else None, "ops": ctr.ops - start,
                            "steps": steps, "walks": walks, "dp_bits": dp_bits,
                            "dps": len(table), "verified": ok, "cycle_escapes": escapes}
                table[W] = (a, b)
            Mj, aj, bj = M[W[0] % r]
            W2 = E.add(W, Mj)
            if W2 is None:
                break
            W2, s = _canon(E, W2)
            a2, b2 = (a + aj) * s % p, (b + bj) * s % p
            if W2 in recent:  # fruitless cycle: escape by doubling the current point
                escapes += 1
                W2 = E.dbl(W)
                if W2 is None:
                    break
                W2, s = _canon(E, W2)
                a2, b2 = 2 * a * s % p, 2 * b * s % p
                recent.clear()
                order.clear()
            recent[W] = steps
            order.append(W)
            if len(order) > cycle_window:
                old = order.pop(0)
                if recent.get(old) is not None and recent[old] <= steps - cycle_window:
                    del recent[old]
            W, a, b = W2, a2, b2

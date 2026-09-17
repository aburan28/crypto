#!/usr/bin/env python3
"""Four-point relations on ECC2K-130, by pair-sum collision.

`P1 + P2 + P3 + P4 = O` says `P1 + P2 = -(P3 + P4)`, and a point and its
negative share an abscissa.  So a four-point relation forces two distinct
pairs whose sums have the *same abscissa*, and the sweep is a collision
search over pair sums rather than an enumeration of quadruples: `B^2/2`
work instead of `B^4/24`.

Three corrections to the plan this experiment started from.

**The signed-Frobenius-orbit grouping does not apply.**  The intent was to
key pair sums by their orbit under `<-1> x <sigma>` and so store one entry
per 262 pairs.  That needs the base to be a union of Frobenius orbits, and
the weight-two base is not one: `z^i + z^j` squares out of weight two as
soon as `2i >= 131`, and only the `i < j <= 65` part survives (2145 of
8515).  Keying by orbit would therefore group sums whose partners are not
in the base, and the saving is not available on this support.  What *is*
available is the `<-1>` half of it, for free: keying on the abscissa of the
pair sum already identifies `S` with `-S`.

**Storing every pair is avoided by storing a digest, not by grouping.**
Each pair sum contributes one 64-bit truncation of `x(S)` to a flat array,
8 bytes rather than a point.  A truncation collision is a candidate, and is
then recomputed and checked in full, so truncation costs false positives
(about `2^-18.8` of them expected here) and never false negatives.

**The honest purpose is testing the boundary's independence assumption.**
`results/boundary.json` shows a base of this size expects `2^-85` four-point
relations, so this sweep is not going to find one and is not run in hope of
it.  What it measures is whether a structured support makes pair sums
collide more often than the birthday prediction -- that is, whether the
`1/r` heuristic the counting floor rests on can be broken by structure.  A
collision count matching prediction is a null result *about the assumption*,
which is the thing worth knowing.

    python3 run_four_point.py [--budget-seconds N] [--base NAME] [--out DIR]
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import random
import sys
import time
from array import array
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import relations as R                                  # noqa: E402
from fastfield import FastGF2m                         # noqa: E402

# A packed sort key: the top DIGEST_BITS are a truncation of x(S), the low
# INDEX_BITS are the pair's index, so the whole key fits one unsigned 64-bit
# word and the collision search is a sort of a flat array.
DIGEST_BITS = 38
INDEX_BITS = 25
DIGEST_MASK = (1 << DIGEST_BITS) - 1
INDEX_MASK = (1 << INDEX_BITS) - 1

# The digest must not be the low bits of x(S). For the two sums over the same
# pair, x(A+B) + x(A-B) = d^2 + d with d = x_B / (x_A + x_B), and squaring in a
# polynomial basis sends bit k to bit 2k, so the low bits of d^2 + d vanish far
# more often than a random element's would. Keying on low bits therefore
# manufactures collisions between exactly the two entries that share a pair.
# A multiplicative mix spreads those bits before truncation. This only ever
# affected the false-positive rate -- a full abscissa collision collides under
# any injective-on-x digest -- but it made the collision count untestable
# against the birthday prediction, which is the measurement this sweep exists
# to make.
MIX = 0x9E3779B97F4A7C15F39CC0605CEDC835
MIX_SHIFT = 67


def digest_of(x: int) -> int:
    return ((x * MIX) >> MIX_SHIFT) & DIGEST_MASK


def representative_points(E, abscissae):
    """One point per abscissa, and the abscissae that carry one.

    The full base is closed under negation, and that is a trap for a
    collision search: the pair `{-P_i, -P_j}` sums to `-(P_i + P_j)`, which
    shares an abscissa with `P_i + P_j`, so every pair would arrive with a
    guaranteed twin and the sweep would report ~18M collisions that are all
    the trivial relation `P + Q - P - Q = O`.

    Enumerating one representative per abscissa and carrying an explicit
    sign on the second summand covers every pair exactly once up to a global
    negation, which is the symmetry a relation is defined modulo anyway. It
    also halves the work.
    """
    reps, kept = [], []
    for x in abscissae:
        pts = E.points_over(x)
        if not pts:
            continue
        reps.append(min(pts, key=lambda p: p[1]))
        kept.append(x)
    return reps, kept


def pair_sum_sweep(E, reps, *, budget_seconds=None, chunk=20000, log=print,
                   name="sweep"):
    """Stream `P_i + eps P_j` over `i < j` and both signs, keeping one key each.

    Yields a flat list of packed keys.  Nothing about the pair is stored
    beyond 64 bits: the pair is recovered from the index when a collision
    needs checking, and recomputed from scratch then.
    """
    n = len(reps)
    total = n * (n - 1)                       # both signs
    keys = []
    started = time.time()
    last = [started]
    consumed = 0
    complete = True
    idx = 0

    buf_a, buf_b = [], []

    def flush():
        nonlocal consumed, idx
        if not buf_a:
            return
        sums = R.batch_add(E, buf_a, buf_b)
        for S in sums:
            if S is not None:
                keys.append((digest_of(S[0]) << INDEX_BITS) | idx)
            idx += 1
        consumed += len(buf_a)
        buf_a.clear(); buf_b.clear()

    stop = False
    for i in range(n):
        Pi = reps[i]
        for j in range(i + 1, n):
            Pj = reps[j]
            buf_a.append(Pi); buf_b.append(Pj)
            buf_a.append(Pi); buf_b.append(E.neg(Pj))
            if len(buf_a) >= chunk:
                flush()
                now = time.time()
                if now - last[0] >= 30:
                    last[0] = now
                    frac = consumed / total
                    el = now - started
                    log(f"[{name}]   {consumed:,}/{total:,} ({frac:6.2%}) "
                        f"kept={len(keys):,} elapsed={el:6.0f}s "
                        f"eta={el / frac - el if frac else 0:6.0f}s")
                if budget_seconds and now - started > budget_seconds:
                    complete = False
                    stop = True
                    break
        if stop:
            break
    flush()
    if idx > INDEX_MASK:
        raise RuntimeError(
            f"pair index {idx} overflows {INDEX_BITS} bits; widen INDEX_BITS")
    return keys, consumed, complete, time.time() - started


def decode_pair(index, n):
    """The `(i, j, sign)` that produced pair number `index`."""
    half, sign = divmod(index, 2)
    i = 0
    rem = half
    # pairs for a given i number (n - 1 - i)
    while rem >= (n - 1 - i):
        rem -= (n - 1 - i)
        i += 1
    j = i + 1 + rem
    return i, j, (1 if sign == 0 else -1)


def find_collisions(keys):
    """Groups of pair indices sharing a digest, found by sorting the keys."""
    keys.sort()
    groups = []
    k = 0
    m = len(keys)
    while k < m:
        d = keys[k] >> INDEX_BITS
        j = k + 1
        while j < m and (keys[j] >> INDEX_BITS) == d:
            j += 1
        if j - k > 1:
            groups.append([keys[t] & INDEX_MASK for t in range(k, j)])
        k = j
    return groups


def frobenius_pattern(E, points):
    """Express a relation as signed Frobenius powers of one of its points.

    Returns `[(k_i, eps_i)]` with `points[i] = eps_i * sigma^{k_i}(points[0])`,
    or `None` if some point is not in the Frobenius orbit of the first.
    """
    R0 = points[0]
    orbit = {}
    cur = R0
    for k in range(E.F.deg):
        orbit.setdefault(cur, (k, 1))
        orbit.setdefault(E.neg(cur), (k, -1))
        cur = E.frobenius(cur)
    pattern = []
    for pt in points:
        hit = orbit.get(pt)
        if hit is None:
            return None
        pattern.append(hit)
    return pattern


def is_frobenius_identity(E, points, trials=6, rng=None):
    """Does this relation's signed exponent pattern annihilate *every* point?

    The curve is a Koblitz curve with `#E(F_2) = 4`, so its trace is `t = -1`
    and Frobenius satisfies its characteristic equation

        sigma^2 + sigma + 2 = 0

    as an endomorphism: `sigma^2(R) + sigma(R) + [2]R = O` for every point `R`
    of the curve, the challenge `P` and `Q` included. Every consequence of it
    is free in the same way -- from `sigma^2 = -sigma - 2` one gets
    `sigma^3 = 2 - sigma`, hence `2R - sigma(R) - sigma^3(R) = O`, a different
    shape and equally empty. So the free relations are not one pattern to
    match against but the whole ideal generated by the characteristic
    equation, and matching a single shape undercounts them.

    The test used here is the definition rather than a shape: read the
    relation as a signed sum of Frobenius powers, then apply that same
    pattern to independent random points of the curve. A pattern that
    annihilates every point is an identity of the endomorphism ring and tells
    nobody anything about a logarithm; a pattern that annihilates only these
    points is a relation about these points.
    """
    pattern = frobenius_pattern(E, points)
    if pattern is None:
        return False
    rng = rng or random.Random(20260917)
    checked = 0
    guard = 0
    while checked < trials and guard < 200:
        guard += 1
        pts = E.points_over(rng.getrandbits(E.F.deg))
        if not pts:
            continue
        S = pts[0]
        acc = None
        for k, eps in pattern:
            T = S
            for _ in range(k):
                T = E.frobenius(T)
            acc = E.add(acc, T if eps > 0 else E.neg(T))
        if acc is not None:
            return False
        checked += 1
    return checked == trials


def verify_quadruples(E, reps, n, groups):
    """Recompute every candidate in full and classify what it actually is.

    Four outcomes, kept apart because collapsing them is how a structural
    identity gets reported as a discovery:

      * `distinct`   -- a relation on four distinct base points, the thing
                        four-point index calculus actually needs;
      * `repeated`   -- a genuine relation with a repeated summand, e.g.
                        `2 R + S + T = O`;
      * `frobenius`  -- a `repeated` relation that is an instance of
                        `sigma^2 + sigma + 2 = 0`, so it holds for every point
                        of the curve and says nothing about any logarithm;
      * `false`      -- a digest collision that is not an abscissa collision.
    """
    out = {"distinct": [], "repeated": [], "frobenius": [], "false": 0}
    for grp in groups:
        decoded = [decode_pair(g, n) for g in grp]
        for a in range(len(decoded)):
            for b in range(a + 1, len(decoded)):
                i, j, si = decoded[a]
                k, l, sk = decoded[b]
                quad = [reps[i], reps[j] if si > 0 else E.neg(reps[j]),
                        reps[k], reps[l] if sk > 0 else E.neg(reps[l])]
                hit = None
                if E.sum_points(quad) is None:
                    hit = quad
                else:
                    alt = quad[:2] + [E.neg(p) for p in quad[2:]]
                    if E.sum_points(alt) is None:
                        hit = alt
                if hit is None:
                    out["false"] += 1
                    continue
                record = [[hex(p[0]), hex(p[1])] for p in hit]
                if len({i, j, k, l}) == 4:
                    if is_frobenius_identity(E, hit):
                        out["frobenius"].append(record)
                    else:
                        out["distinct"].append(record)
                else:
                    out["repeated"].append(record)
                    if is_frobenius_identity(E, hit):
                        out["frobenius"].append(record)
    return out


def run(name, description, E, abscissae, r, *, budget_seconds, log=print):
    reps, kept = representative_points(E, abscissae)
    n = len(reps)
    total_pairs = n * (n - 1)
    log(f"[{name}] {len(abscissae)} abscissae -> {n} representatives, "
        f"{total_pairs:,} signed pairs")

    keys, consumed, complete, elapsed = pair_sum_sweep(
        E, reps, budget_seconds=budget_seconds, log=log, name=name)
    stored = len(keys)
    groups = find_collisions(keys)
    cls = verify_quadruples(E, reps, n, groups)

    candidates = sum(len(g) * (len(g) - 1) // 2 for g in groups)
    log2_pred_digest = (2 * math.log2(stored) - 1 - DIGEST_BITS
                        if stored > 1 else None)
    log2_pred_real = (2 * math.log2(stored) - 1 - 130.0
                      if stored > 1 else None)

    log(f"[{name}] done: {candidates} candidates | "
        f"{len(cls['distinct'])} four-distinct | {len(cls['repeated'])} repeated "
        f"({len(cls['frobenius'])} Frobenius) | {cls['false']} false positives")

    return {
        "name": name,
        "description": description,
        "abscissae_supplied": len(abscissae),
        "representatives": n,
        "signed_pairs_total": total_pairs,
        "signed_pairs_consumed": consumed,
        "complete": complete,
        "elapsed_seconds": round(elapsed, 3),
        "pair_rate_per_second": round(consumed / elapsed, 1) if elapsed else None,
        "keys_stored": stored,
        "key_bytes": stored * 8,
        "digest_bits": DIGEST_BITS,
        "digest_collision_groups": len(groups),
        "candidates_checked": candidates,
        "false_positive_candidates": cls["false"],
        "relations_four_distinct": len(cls["distinct"]),
        "relations_repeated_summand": len(cls["repeated"]),
        "relations_frobenius_characteristic": len(cls["frobenius"]),
        "relations_repeated_not_frobenius":
            len(cls["repeated"]) - len(cls["frobenius"]),
        "relations": {
            "distinct": cls["distinct"][:16],
            "repeated": cls["repeated"][:16],
            "frobenius": cls["frobenius"][:16],
        },
        "log2_predicted_digest_collisions": (
            round(log2_pred_digest, 3) if log2_pred_digest is not None else None),
        "log2_predicted_genuine_collisions": (
            round(log2_pred_real, 3) if log2_pred_real is not None else None),
        "log2_expected_relations_by_counting": round(
            4 * math.log2(2 * n) - math.log2(24) - math.log2(r), 3),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--budget-seconds", type=float, default=7200.0)
    ap.add_argument("--out", default=str(HERE / "results"))
    ap.add_argument("--only", default=None)
    ap.add_argument("--control-size", type=int, default=140)
    args = ap.parse_args()

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    F, E, P, Q, r = R.challenge_curve()
    rng = random.Random(20260917)

    w2 = R.base_from_abscissae(E, R.weight_two_abscissae(F))
    stable = R.base_from_abscissae(
        E, [(1 << i) ^ (1 << j) for i in range(66) for j in range(i + 1, 66)])

    specs = [
        ("four_weight_two", "the complete weight-two support", w2),
        ("four_weight_two_frobenius_stable",
         "the sigma-stable i<j<=65 part of the weight-two support", stable),
        ("four_random_matched", "uniform random abscissae, size matched",
         R.random_matched_base(E, len(w2), rng)),
    ]

    summary = {
        "instance": "ECC2K-130",
        "subgroup_order_r": str(r),
        "digest_bits": DIGEST_BITS,
        "budget_seconds": args.budget_seconds,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "experiments": [],
    }

    for name, desc, xs in specs:
        if args.only and name != args.only:
            continue
        res = run(name, desc, E, xs, r, budget_seconds=args.budget_seconds)
        (outdir / f"{name}.json").write_text(json.dumps(res, indent=2))
        summary["experiments"].append(
            {k: v for k, v in res.items() if k != "relations"})
        (outdir / "four_point_summary.json").write_text(
            json.dumps(summary, indent=2))
        print(json.dumps(summary["experiments"][-1], indent=2))


if __name__ == "__main__":
    main()

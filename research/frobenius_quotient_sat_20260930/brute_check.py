#!/usr/bin/env python3
"""Exhaustively re-check every UNSAT / verified answer in a results file.

For each target, decide by pair matching whether it is a signed sum of four
points with DISTINCT x-coordinates drawn from the admissible set: the
factor-space representatives themselves, or (for phase runs) all their
Frobenius conjugates.  Any disagreement with the solver is a mismatch."""
import collections, json, sys

from encode import Curve


def admissible(E, basis, phases):
    span = {0}
    for b in basis:
        span |= {w ^ b for w in span}
    reps = sorted(x for x in span if x and (P := E.lift(x)) is not None
                  and E.mul(E.h, P) is None)
    xs = set(reps)
    if phases:
        for x in reps:
            for _ in range(E.f.n):
                x = E.f.sq(x)
                xs.add(x)
    return reps, sorted(xs)


def decomposable(E, xs, T):
    pts = [(x, P) for x in xs for P in (E.lift(x), E.neg(E.lift(x)))]
    left = collections.defaultdict(list)
    for i, (xa, A) in enumerate(pts):
        for xb, B in pts[i + 1:]:
            if xa != xb:
                left[E.add(A, B)].append((xa, xb))
    for S, prs in left.items():
        need = E.add(T, E.neg(S))
        for xa, xb in left.get(need, ()):
            for xc, xd in prs:
                if len({xa, xb, xc, xd}) == 4:
                    return True
    return False


def main(path):
    lines = [json.loads(l) for l in open(path)]
    meta = lines[0]["meta"]
    E = Curve(meta["n"])
    basis = [int(b, 16) for b in meta["basis"]]
    phases = any(r.get("phases") for r in lines[1:])
    reps, xs = admissible(E, basis, phases)
    cache, bad = {}, 0
    for r in lines[1:]:
        T = (int(r["target"][0], 16), int(r["target"][1], 16))
        if T not in cache:
            cache[T] = decomposable(E, xs, T)
        if (r["status"] == "unsat" and cache[T]) or \
           (r["status"] == "verified" and not cache[T]):
            bad += 1
            print("MISMATCH", r["formulation"], r["status"], cache[T])
    print(json.dumps({"file": path, "phases": phases, "valid_reps": len(reps),
                      "admissible_x": len(xs), "targets": len(cache),
                      "decomposable_targets": sum(cache.values()),
                      "attempts": len(lines) - 1, "mismatches": bad}))


if __name__ == "__main__":
    for p in sys.argv[1:]:
        main(p)

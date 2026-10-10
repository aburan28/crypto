#!/usr/bin/env sage -python
"""Small-field cross-isogeny relation-yield pilot. Requires SageMath.

Usage:
  sage -python experiments/cross_isogeny_symmetry/pilot.py --prime 7 --extension 2 --ell 3 --output /tmp/isogeny-pilot.json

This is a *toy* experiment, not a practical ECDLP attack. It compares factor-base
membership and two-term relation counts; it never claims a complexity speedup.
"""
import argparse
import json
import random
import time
from sage.all import GF, EllipticCurve, ZZ


def key(P):
    return str(P)


def frob(P, p):
    if P.is_zero():
        return P
    E = P.curve()
    return E([P[0] ** p, P[1] ** p, 1])


def orbit(P, p, max_steps=10000):
    out = []
    seen = set()
    while key(P) not in seen:
        if len(out) >= max_steps:
            raise RuntimeError("Frobenius orbit bound exceeded")
        seen.add(key(P))
        out.append(P)
        P = frob(P, p)
    return out


def relation_stats(points, factor_base, targets):
    # Exact exhaustive two-term decompositions, with ordered pairs.
    start = time.perf_counter()
    fb = list(factor_base)
    sums = {}
    for P in fb:
        for R in fb:
            S = key(P + R)
            sums[S] = sums.get(S, 0) + 1
    hits = [sums.get(key(Q), 0) for Q in targets]
    return {
        "factor_base_size": len(fb),
        "factor_base_density": len(fb) / len(points),
        "targets": len(targets),
        "targets_with_relation": sum(h > 0 for h in hits),
        "total_ordered_decompositions": sum(hits),
        "enumeration_seconds": time.perf_counter() - start,
    }


def run(p, extension, ell, seed, max_points):
    if p < 5 or not ZZ(p).is_prime():
        raise ValueError("prime must be prime >=5")
    if extension < 1:
        raise ValueError("extension must be positive")
    F = GF(p ** extension, name="z")
    rng = random.Random(seed)
    # Seek an ordinary nonsingular curve with a rational degree-ell isogeny.
    attempts = 0
    while attempts < 250:
        attempts += 1
        a = F(rng.randrange(p))
        b = F(rng.randrange(p))
        try:
            E = EllipticCurve(F, [a, b])
            if E.is_supersingular():
                continue
            isogs = list(E.isogenies_prime_degree(ell))
            if isogs:
                break
        except (ArithmeticError, ValueError, NotImplementedError):
            continue
    else:
        raise RuntimeError("No suitable ordinary curve/isogeny found; try another seed or ell")
    phi = isogs[0]
    E2 = phi.codomain()
    points = list(E.points())
    if len(points) > max_points:
        raise ValueError("Curve too large for exhaustive toy enumeration")
    images = [phi(P) for P in points]
    # Full q-Frobenius commutation is required for maps defined over F_q.
    q = p ** extension
    assert all(phi(frob(P, q)) == frob(phi(P), q) for P in points)
    # p-Frobenius is generally NOT a self-map of an arbitrary curve over F_(p^n);
    # here a,b are deliberately chosen in F_p to make it one.
    assert all(phi(frob(P, p)) == frob(phi(P), p) for P in points)
    orbit_keys = {tuple(sorted(key(R) for R in orbit(P, p))) for P in points}
    # Baseline: x trace zero. Treatment: trace zero on image under phi.
    def trace_zero(P):
        return P.is_zero() or P[0].trace() == 0
    fb_base = [P for P in points if trace_zero(P)]
    fb_pullback = [P for P, image in zip(points, images) if trace_zero(image)]
    # Fair-size control: equal cardinality sampled from the baseline if possible.
    n = min(len(fb_base), len(fb_pullback))
    rng.shuffle(fb_base)
    rng.shuffle(fb_pullback)
    fb_base = fb_base[:n]
    fb_pullback = fb_pullback[:n]
    targets = list(points)
    rng.shuffle(targets)
    targets = targets[:min(50, len(targets))]
    base = relation_stats(points, fb_base, targets)
    treatment = relation_stats(points, fb_pullback, targets)
    return {
        "schema_version": 1, "status": "toy_exhaustive_not_cryptanalytic",
        "seed": seed, "p": p, "extension": extension, "q": q,
        "ell": ell, "curve_a": str(E), "curve_b": str(E2),
        "curve_order": len(points), "isogeny_degree": int(phi.degree()),
        "frobenius_orbits": len(orbit_keys),
        "factor_base_rule_baseline": "trace(x(P)) == 0",
        "factor_base_rule_treatment": "trace(x(phi(P))) == 0",
        "baseline": base, "treatment": treatment,
        "warning": "Two-term relation yield is not independent verified index-calculus relations; no solver or linear algebra included."
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prime", type=int, default=7)
    ap.add_argument("--extension", type=int, default=2)
    ap.add_argument("--ell", type=int, default=3)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--max-points", type=int, default=2000)
    ap.add_argument("--output", default="-")
    args = ap.parse_args()
    result = run(args.prime, args.extension, args.ell, args.seed, args.max_points)
    payload = json.dumps(result, indent=2, sort_keys=True)
    if args.output == "-":
        print(payload)
    else:
        with open(args.output, "w", encoding="utf-8") as f:
            f.write(payload + "\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env sage
# Diagnostic: why does the standard base yield zero point decompositions in
# every symmetric-S4 cell while geometric/random yield ~2^-5..2^-8?
# Replicates the engine's own x-level S3-chain labeling (label_targets_m3)
# to enumerate the full reachable target set S per structure, then inspects
# |S|, its linear span, and trace obstructions. Recomputes from the
# engine's functions; reads no files.
# Usage: sage analyze_standard.sage
import sys
load("../../sage/factor_base_dor.sage")

def reachable(field, spaces, b_int):
    K = field.K
    T = PolynomialRing(K, "T").gen()
    m1 = [field.from_int(v) for v in span_members(spaces[0])]
    m2 = m1
    m3 = m1
    b = field.from_int(b_int)
    roots_u = []
    for x1 in m1:
        for x2 in m2:
            poly = (x1 + x2)**2 * T**2 + x1 * x2 * T + (x1 * x2)**2 + b
            roots_u.extend(r for r, _ in poly.roots())
    roots_u = list(set(roots_u))
    S = set()
    for u in roots_u:
        for x3 in m3:
            p = (u + x3)**2 * T**2 + u * x3 * T + (u * x3)**2 + b
            for r, _ in p.roots():
                S.add(field.to_int(r))
    return S

def diag(n, l, seed=1):
    field = BinaryField(n)
    b_int = field.one_int
    print("n=%d l=%d seed=%d" % (n, l, seed))
    for name in ("geometric", "standard", "random"):
        spaces, err = symmetric_spaces(field, name, l, seed)
        if spaces is None:
            print("  %-10s skipped (%s)" % (name, err)); continue
        S = reachable(field, spaces, b_int)
        span = subspace_from_generators(field, [field.from_int(t) for t in S]) if S else []
        tr = {}
        for t in S:
            tv = int(field.from_int(t).trace())
            tr[tv] = tr.get(tv, 0) + 1
        print("  %-10s |S|=%d span_dim=%d pfrac=%.5f traces=%s"
              % (name, len(S), len(span), len(S) / float(2**n), tr))

diag(15, 2)
diag(15, 3)
diag(21, 3)

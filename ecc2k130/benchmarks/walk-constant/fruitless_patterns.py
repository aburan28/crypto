"""Fruitless cycles of the table walk that the device's cycle rule lets through.

ecc2k130/WALK-CONSTANT.md section 5.  A table step adds eps * sigma^k(T_h), and
Frobenius sigma satisfies tau^2 + tau + 2 = 0 on these curves, so a run of
steps returns a walk to a point it has already left, without any collision,
whenever for every branch h the signed sum of tau^k over that branch's steps
is 0 in Z[tau].  Tags that cancel in pairs are one way; tau-relations such as
tau^2 + tau + 1 + 1 = 0 are another.  The rule (eccTagFruitless in
include/tablewalk.h) refuses a step that undoes the last one or that closes a
pairwise 4-cycle with the three before it; this counts what it does not
refuse.

A pattern is a cyclic run of L tags.  The first tag of each branch-group is
free; every other tag is determined relative to it (same h, a given k offset
and sign).  A walk enters a given pattern at a given step with probability
(product over determined tags of the chance it comes up), which to leading
order is (sum p^2 / 2m) per determined tag when the tag pairs with one other
tag of its group, and sum p^4 / (2m)^3 for a four-step group.  Counted: the
patterns with at most MAXDET determined tags that are primitive (no shorter
return inside, at any rotation) and that the rule passes at every position
of the steady state (the cycle repeats forever).

Run: python3 fruitless_patterns.py [LMAX=6] [D=6] [MAXDET=3]
"""
import sys


def tau_pow(k):
    """tau^k = a + b tau for k >= 0, from tau^(j+1) = -2b + (a - b) tau."""
    a, b = 1, 0
    for _ in range(k):
        a, b = -2 * b, a - b
    return a, b


def zero_sum(group):
    """Whether sum eps tau^d over (d, eps) is 0 in Z[tau]."""
    lo = min(d for d, _ in group)
    A = B = 0
    for d, e in group:
        a, b = tau_pow(d - lo)
        A += e * a
        B += e * b
    return A == 0 and B == 0


def returns(tags):
    groups = set(g for g, _, _ in tags)
    return all(zero_sum([(d, e) for g, d, e in tags if g == gg]) for gg in groups)


def negation(x, y):
    return x[0] == y[0] and x[1] == y[1] and x[2] == -y[2]


def passes_rule(tags):
    """The rule never fires in the steady state: at each position, with the
    three before it (cyclically) as history."""
    L = len(tags)
    for i in range(L):
        t, t1, t2, t3 = tags[i], tags[(i - 1) % L], tags[(i - 2) % L], tags[(i - 3) % L]
        if negation(t, t1) or (negation(t, t2) and negation(t1, t3)):
            return False
    return True


def primitive(tags):
    L = len(tags)
    for i in range(L):
        for l in range(1, L):
            if returns([tags[(i + j) % L] for j in range(l)]):
                return False
    return True


def count(L, D, maxdet):
    found = {}

    def rec(seq, ngroups):
        if len(seq) - ngroups > maxdet:
            return
        if len(seq) == L:
            if returns(seq) and passes_rule(seq) and primitive(seq):
                key = (L - ngroups, tuple(sorted(sum(1 for g, _, _ in seq if g == gg) for gg in range(ngroups))))
                found[key] = found.get(key, 0) + 1
            return
        rec(seq + [(ngroups, 0, 1)], ngroups + 1)
        for g in range(ngroups):
            for d in range(-D, D + 1):
                for e in (1, -1):
                    rec(seq + [(g, d, e)], ngroups)

    rec([], 0)
    return found


if __name__ == "__main__":
    lmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    D = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    maxdet = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    for L in range(2, lmax + 1):
        found = count(L, D, maxdet)
        rows = ", ".join(f"{n} with {det} determined tags in groups of {sizes}" for (det, sizes), n in sorted(found.items()))
        print(f"L = {L}: {rows or 'none'}", flush=True)

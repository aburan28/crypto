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

Run: python3 fruitless_patterns.py [LMAX=6] [D=6] [MAXDET=3] [--rule v1|v2]
(v1, the default, is the rule before WALK-CONSTANT.md section 11 and
reproduces fruitless_patterns.txt; v2 is eccTagFruitless today.)
     python3 fruitless_patterns.py --residual 8
counts what v2 passes by group signature, out to six determined tags, which
the brute force above cannot reach (fruitless_patterns_v2.txt).
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


def tau_relation(t, t1, t2, t3):
    """eccTauRelation on pattern tags (group, offset, sign): one group, one
    repeated tag, its partners at +1, +2 (same sign) or +1, +3 (opposite)."""
    four = [t, t1, t2, t3]
    if len(set(g for g, _, _ in four)) != 1:
        return False
    for i in range(4):
        for j in range(i + 1, 4):
            if four[i] == four[j]:
                g, d, e = four[i]
                rest = sorted((x[1], x[2]) for k, x in enumerate(four) if k not in (i, j))
                return rest in (sorted([(d + 1, e), (d + 2, e)]), sorted([(d + 1, -e), (d + 3, -e)]))
    return False


def refused(rule, t, back):
    """Whether the rule refuses tag t after `back` (most recent first)."""
    if rule == "v1":
        t1, t2, t3 = back[0], back[1], back[2]
        return negation(t, t1) or (negation(t, t2) and negation(t1, t3))
    if any(negation(t, u) for u in back[:4]):
        return True
    return tau_relation(t, back[0], back[1], back[2])


def passes_rule(tags, rule):
    """The rule never fires in the steady state: at each position, with the
    four before it (cyclically) as history."""
    L = len(tags)
    return not any(refused(rule, tags[i], [tags[(i - j) % L] for j in range(1, 5)]) for i in range(L))


def linear_ok(seq, rule):
    """The rule does not fire on the last tag of a partial run, looking back
    only inside the run: a necessary condition for any completion."""
    i = len(seq) - 1
    back = [seq[i - j] if i - j >= 0 else ("none", 0, 0) for j in range(1, 5)]
    return not refused(rule, seq[i], back)


def primitive(tags):
    L = len(tags)
    for i in range(L):
        for l in range(1, L):
            if returns([tags[(i + j) % L] for j in range(l)]):
                return False
    return True


def count(L, D, maxdet, rule):
    found = {}

    def rec(seq, ngroups):
        if len(seq) - ngroups > maxdet:
            return
        if seq and not linear_ok(seq, rule):
            return
        # A contiguous run that already returns makes any completion non-primitive.
        for l in range(1, len(seq)):
            if returns(seq[len(seq) - l - 1:len(seq) - 1] + [seq[-1]]) and l + 1 < L:
                return
        if len(seq) == L:
            if returns(seq) and passes_rule(seq, rule) and primitive(seq):
                sizes = tuple(sorted(sum(1 for g, _, _ in seq if g == gg) for gg in range(ngroups)))
                key = (L - ngroups, sizes)
                found[key] = found.get(key, 0) + 1
            return
        rec(seq + [(ngroups, 0, 1)], ngroups + 1)
        for g in range(ngroups):
            for d in range(-D, D + 1):
                for e in (1, -1):
                    rec(seq + [(g, d, e)], ngroups)

    rec([], 0)
    return found


def group_relations(size, span):
    """Irreducible zero-sum multisets of `size` signed powers tau^d,
    0 <= d <= span, least power 0, one of each +- pair: the ways one branch's
    steps can cancel.  Irreducible (no proper part sums to zero) so that each
    pattern has one signature: two pairs on one branch are two groups of 2."""
    from itertools import combinations, combinations_with_replacement
    vals = [(d, e) for d in range(span + 1) for e in (1, -1)]
    last = {}
    for d, e in vals:
        a, b = tau_pow(d)
        last[(e * a, e * b)] = (d, e)
    out = set()
    for combo in combinations_with_replacement(vals, size - 1):
        A = sum(e * tau_pow(d)[0] for d, e in combo)
        B = sum(e * tau_pow(d)[1] for d, e in combo)
        t = last.get((-A, -B))
        if t is None:
            continue
        ms = tuple(sorted(combo + (t,)))
        if ms[0][0] != 0 or any(zero_sum(list(part)) for r in range(2, size - 1)
                                for part in combinations(ms, r)):
            continue
        out.add(min(ms, tuple(sorted((d, -e) for d, e in ms))))
    return sorted(out)


def structured(L, sizes, rule, span):
    """Patterns of L steps in groups of the given sizes, each group one of its
    zero-sum multisets in some order, counted as `count` counts them (each
    group's first tag free, the rest relative to it).  This reaches the
    patterns with many determined tags that `count` cannot enumerate in
    reasonable time; it agrees with `count` where both run."""
    from itertools import combinations, permutations, product
    orders = {}
    for s in set(sizes):
        seqs = set()
        for rel in group_relations(s, span):
            for p in set(permutations(rel)):
                seqs.add(tuple((d - p[0][0], e * p[0][1]) for d, e in p))
        orders[s] = sorted(seqs)
    total = 0

    def assign(slots, left, groups):
        nonlocal total
        free = [i for i in range(L) if slots[i] is None]
        if not free:
            for choice in product(*(orders[len(g)] for g in groups)):
                tags = [None] * L
                for gi, (g, o) in enumerate(zip(groups, choice)):
                    for pos, (d, e) in zip(g, o):
                        tags[pos] = (gi, d, e)
                total += passes_rule(tags, rule) and primitive(tags)
            return
        first = free[0]
        for s in sorted(set(left)):
            rest = list(left)
            rest.remove(s)
            for others in combinations(free[1:], s - 1):
                g = (first,) + others
                for pos in g:
                    slots[pos] = len(groups)
                assign(slots, rest, groups + [g])
                for pos in g:
                    slots[pos] = None

    assign([None] * L, list(sizes), [])
    return total


# Group signatures by determined tags.  A group has even size (every power of
# tau is 1 mod tau-bar, of norm 2), so these are all the signatures with at
# most four determined tags, which rule v2 must pass none of: it refuses a tag
# that negates any of the four before it, so a pair in an L-step cycle needs
# cyclic distance 5 or more both ways (L >= 10), and the tau test takes the
# four-step groups.
LOW_ORDER = [(2, (2,)), (4, (4,)), (4, (2, 2)), (6, (4, 2)), (6, (2, 2, 2)), (8, (2, 2, 2, 2))]
# And the signatures with five or six determined tags that rule v2 can pass
# (the others with five or six have a pair in fewer than ten steps).
RESIDUAL_V2 = [(6, (6,)), (8, (4, 4)), (10, (2, 2, 2, 2, 2)), (10, (4, 2, 2, 2)), (12, (2, 2, 2, 2, 2, 2))]


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("lmax", nargs="?", type=int, default=6)
    ap.add_argument("D", nargs="?", type=int, default=6)
    ap.add_argument("maxdet", nargs="?", type=int, default=3)
    ap.add_argument("--rule", choices=("v1", "v2"), default="v1",
                    help="v1: the rule before WALK-CONSTANT.md section 11; v2: eccTagFruitless today")
    ap.add_argument("--residual", type=int, metavar="SPAN",
                    help="count the v2 rule's survivors with five and six determined tags, relation spans up to SPAN")
    args = ap.parse_args()
    if args.residual is not None:
        print(f"rule v2 residual, relation spans <= {args.residual}", flush=True)
        for L, sizes in LOW_ORDER:
            n = structured(L, sizes, "v2", args.residual)
            print(f"at most four determined: L = {L}, groups {sizes}, {L - len(sizes)} determined tags: {n} patterns",
                  flush=True)
        for L, sizes in RESIDUAL_V2:
            n = structured(L, sizes, "v2", args.residual)
            print(f"residual v2: L = {L}, groups {sizes}, {L - len(sizes)} determined tags: {n} patterns", flush=True)
        sys.exit(0)
    print(f"rule {args.rule}, offsets |d| <= {args.D}, at most {args.maxdet} determined tags", flush=True)
    for L in range(2, args.lmax + 1):
        found = count(L, args.D, args.maxdet, args.rule)
        rows = ", ".join(f"{n} with {det} determined tags in groups of {sizes}" for (det, sizes), n in sorted(found.items()))
        print(f"L = {L}: {rows or 'none'}", flush=True)

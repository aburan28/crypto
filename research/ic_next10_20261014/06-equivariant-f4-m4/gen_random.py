#!/usr/bin/env python3
"""Arm D: random Boolean systems with a template's shape.

usage: python3 -I gen_random.py TEMPLATE.ms SEED OUT.ms
Keeps the template's variable count, equation count, and for every equation
its number of terms of each degree; draws the monomials of each degree
uniformly at random (distinct) and appends the field equations.  Prints the
semi-regular D_reg of the shape (first non-positive coefficient of
(1+z)^N / prod (1+z^{d_i})).
"""
import random, sys
from math import comb

tpl, seed, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
lines = open(tpl).read().split("\n")
vars_ = lines[0].split(",")
N = len(vars_)
polys = [p.strip() for p in "\n".join(lines[2:]).split(",") if p.strip()]
polys = [p for p in polys if "^2" not in p]  # drop the template's field equations
rng = random.Random(seed)
eqs, degs = [], []
for p in polys:
    by_deg = {}
    for t in p.split("+"):
        d = 0 if t == "1" else len(t.split("*"))
        by_deg[d] = by_deg.get(d, 0) + 1
    terms = []
    for d, cnt in sorted(by_deg.items()):
        if d == 0:
            terms.append("1")
            continue
        chosen = set()
        while len(chosen) < min(cnt, comb(N, d)):
            chosen.add(tuple(sorted(rng.sample(range(N), d))))
        terms += ["*".join(f"x{i}" for i in mono) for mono in sorted(chosen)]
    eqs.append("+".join(terms))
    degs.append(max(by_deg))
eqs += [f"x{i}^2+x{i}" for i in range(N)]
with open(out, "w") as f:
    f.write(",".join(vars_) + "\n2\n" + ",\n".join(eqs) + "\n")
cap = 2 * N + 8
c = [comb(N, i) if i <= N else 0 for i in range(cap + 1)]
for d in degs:
    for i in range(d, cap + 1):
        c[i] -= c[i - d]
dreg = next(i for i, x in enumerate(c) if x <= 0)
print(f"N={N} eqs={len(degs)} max_deg={max(degs)} semireg_dreg={dreg}")

#!/usr/bin/env python3
"""Tabulate the JSONL cells written by lfd.sage / lfd_fast.sage.

    python3 research/fall_degree_bounds_20261008/table.py

Columns: n, n', family, target (sat = decomposable over V, unsat = no
decomposition), draws, first fall degree, plain-XL full-ideal degree
(None where the fast script stopped early), HKY last fall degree, and the
semi-regular D_reg of n quadrics in n Boolean unknowns as a reference.
"""
import glob, json, os
from collections import defaultdict
from math import comb

def dreg(n_vars, degrees, cap=60):
    c = [comb(n_vars, i) for i in range(cap + 1)]
    for d in degrees:
        out = [0] * (cap + 1)
        for i in range(cap + 1):
            out[i] = c[i] - (out[i - d] if i >= d else 0)
        c = out
    return next(i for i, x in enumerate(c) if x <= 0)

here = os.path.dirname(os.path.abspath(__file__))
cells = defaultdict(list)
for f in sorted(glob.glob(os.path.join(here, "runs", "*.jsonl"))):
    for line in open(f):
        r = json.loads(line)
        cells[(r["n"], r["nprime"], r["family"], "unsat" if r["unsat"] else "sat")].append(r)

def fmt(vals):
    vals = [("≥%d" % (v["cap"] + 1)) if v is None else str(v) for v in vals]
    return " ".join(vals)

print("| n | n' | family | target | draws | d_ff | plain-XL full-ideal degree | d_F (HKY last fall) | semi-regular D_reg |")
print("|--:|--:|---|---|--:|---|---|---|--:|")
for key in sorted(cells):
    n, np_, fam, tgt = key
    rs = cells[key]
    ffd = " ".join(str(r["ffd"]) for r in rs)
    xl = " ".join(str(r["xl_solving"]) if r["xl_solving"] is not None else "–" for r in rs)
    dl = " ".join(str(r["d_last"]) if r["d_last"] is not None else ("≥%d" % (max(int(k) for k in r["closure_dims"]) + 1)) for r in rs)
    print(f"| {n} | {np_} | {fam} | {tgt} | {len(rs)} | {ffd} | {xl} | {dl} | {dreg(n, [2]*n)} |")

#!/usr/bin/env python3
"""Semi-regular prediction for each exported Boolean system.

For a system in v Boolean unknowns with equations of degrees d_1..d_k, the semi-regular
(generic) Hilbert series of the top-degree forms is

    H(z) = (1 + z)^v / prod_i (1 + z^{d_i}),

and the degree of regularity d_reg is the index of its first non-positive coefficient.
The registered window for the refutation degree is [d_reg, d_reg + 1].

    predict.py DUMP_DIR > predictions.jsonl     (one line per <cell>-d<draw>-<arm>.sing)

Reads only the ring header and the polynomial degrees; nothing is measured here.
"""
import json
import math
import os
import re
import sys


def system_shape(path):
    s = open(path).read()
    v = int(re.search(r"x\(1\.\.(\d+)\)", s).group(1))
    body = s.split("ideal I =", 1)[1].strip().rstrip(";")
    degs = []
    for poly in body.split(","):
        terms = [t.strip() for t in poly.split("+")]
        degs.append(max(0 if t == "1" else len(t.split("*")) for t in terms))
    return v, degs


def d_reg(v, degs):
    h = [math.comb(v, k) for k in range(v + 2)]
    for d in degs:
        if d == 0:  # a constant equation: the ideal is the whole ring
            return 0
        out = [0] * len(h)
        for k in range(len(h)):
            out[k] = h[k] - (out[k - d] if k >= d else 0)
        h = out
    return next((k for k, c in enumerate(h) if c <= 0), None)


def main():
    dump = sys.argv[1]
    for name in sorted(os.listdir(dump)):
        m = re.fullmatch(r"(K[01]n\d+l\d+)-d(\d+)-(rr|x4|ctrl)\.sing", name)
        if not m:
            continue
        cell, draw, arm = m.group(1), int(m.group(2)), m.group(3)
        v, degs = system_shape(os.path.join(dump, name))
        dr = d_reg(v, degs)
        hist = {}
        for d in degs:
            hist[str(d)] = hist.get(str(d), 0) + 1
        print(json.dumps({"cell": cell, "draw": draw, "arm": arm, "n_vars": v,
                          "n_eqs": len(degs), "degrees": hist, "d_reg": dr,
                          "window": None if dr is None else [dr, dr + 1]}))


if __name__ == "__main__":
    main()

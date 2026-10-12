#!/usr/bin/env python3
"""Post-hoc diagnostic (not registered; not evidence for the verdict).

For each exported system: the number of equations of each degree, the rank of their
top-degree forms, the degree profile after graded Gaussian elimination of the
equations themselves (duplicates and combinations removed), and the semi-regular
d_reg of that reduced profile.  Run after the round to explain readings outside the
registered windows.

    topforms.py DUMP_DIR [CELL ...] > topforms.jsonl
"""
import json
import math
import os
import re
import sys


def parse(path):
    s = open(path).read()
    v = int(re.search(r"x\(1\.\.(\d+)\)", s).group(1))
    polys = []
    for p in s.split("ideal I =", 1)[1].strip().rstrip(";").split(","):
        x = 0
        for t in (t.strip() for t in p.split("+")):
            m = 0 if t == "1" else sum(1 << (int(y.strip()[2:-1]) - 1) for y in t.split("*"))
            x ^= 1 << m
        polys.append(x)
    return v, polys


def monos(x):
    while x:
        hb = x.bit_length() - 1
        yield hb
        x &= ~(1 << hb)


def deg(x):
    return max((bin(m).count("1") for m in monos(x)), default=-1)


def lead(x):
    return max((bin(m).count("1"), m) for m in monos(x))


def reduced(polys):
    basis = {}
    for p in polys:
        x = p
        while x:
            ld = lead(x)
            if ld in basis:
                x ^= basis[ld]
            else:
                basis[ld] = x
                break
    return sorted((deg(b) for b in basis.values()), reverse=True)


def top_rank(polys, d):
    basis = {}
    r = 0
    for p in polys:
        if deg(p) != d:
            continue
        x = 0
        for m in monos(p):
            if bin(m).count("1") == d:
                x ^= 1 << m
        while x:
            hb = x.bit_length() - 1
            if hb in basis:
                x ^= basis[hb]
            else:
                basis[hb] = x
                r += 1
                break
    return r


def d_reg(v, degs):
    h = [math.comb(v, k) for k in range(v + 2)]
    for d in degs:
        if d <= 0:
            return 0
        o = [0] * len(h)
        for k in range(len(h)):
            o[k] = h[k] - (o[k - d] if k >= d else 0)
        h = o
    return next((k for k, c in enumerate(h) if c <= 0), None)


def main():
    dump, cells = sys.argv[1], set(sys.argv[2:])
    for name in sorted(os.listdir(dump)):
        m = re.fullmatch(r"(K[01]n\d+l\d+)-d(\d+)-(rr|x4)\.sing", name)
        if not m or (cells and m.group(1) not in cells):
            continue
        v, polys = parse(os.path.join(dump, name))
        degs = [deg(p) for p in polys]
        prof = reduced(polys)
        hist = {}
        for d in degs:
            hist[str(d)] = hist.get(str(d), 0) + 1
        rhist = {}
        for d in prof:
            rhist[str(d)] = rhist.get(str(d), 0) + 1
        print(json.dumps({"cell": m.group(1), "draw": int(m.group(2)), "arm": m.group(3), "n_vars": v,
                          "degrees": hist, "top_form_rank": {str(d): top_rank(polys, d) for d in sorted(set(degs))},
                          "reduced_degrees": rhist, "d_reg_face": d_reg(v, degs), "d_reg_reduced": d_reg(v, prof)}))


if __name__ == "__main__":
    main()

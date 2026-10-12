#!/usr/bin/env python3
"""Round 2 predictions, per exported system, written before any measurement.

x4: d_aff = first k at which the running sum of the semi-regular series
    (1+z)^v / prod(1+z^{d_i}) (face-value degrees) is <= 0.  Prediction: the
    refutation degree equals d_aff (secondary: it is at most d_aff).
rr: the equations are first reduced by graded Gaussian elimination (topforms.py);
    d_red = d_reg of the reduced degree profile.  Prediction: d_red - 1 for n >= 29,
    d_red for n < 29 (the rule round 1's diagnostic suggested; fitted on round 1).

    predict2.py DUMP_DIR > predictions.jsonl
"""
import json, math, os, re, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from topforms import parse, deg, reduced, d_reg  # noqa: E402


def series(v, degs, K):
    h = [math.comb(v, k) for k in range(K)]
    for d in degs:
        o = [0] * K
        for k in range(K):
            o[k] = h[k] - (o[k - d] if k >= d else 0)
        h = o
    return h


def d_aff(v, degs):
    if min(degs) == 0:
        return 0
    s = 0
    for k, c in enumerate(series(v, degs, v + 2)):
        s += c
        if s <= 0:
            return k
    return None


for name in sorted(os.listdir(sys.argv[1])):
    m = re.fullmatch(r"(K[01]n(\d+)l\d+)-d(\d+)-(rr|x4)\.sing", name)
    if not m:
        continue
    cell, n, draw, arm = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)
    v, polys = parse(os.path.join(sys.argv[1], name))
    degs = [deg(p) for p in polys]
    out = {"cell": cell, "draw": draw, "arm": arm, "n_vars": v, "n_eqs": len(degs)}
    if arm == "x4":
        out.update(d_reg=d_reg(v, degs), d_aff=d_aff(v, degs))
        out["predicted"] = out["d_aff"]
    else:
        prof = reduced(polys)
        dr = d_reg(v, prof)
        out.update(reduced_eqs=len(prof), d_red=dr)
        out["predicted"] = 0 if dr == 0 else (dr - 1 if n >= 29 else dr)
    print(json.dumps(out), flush=True)

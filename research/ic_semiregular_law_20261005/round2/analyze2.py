#!/usr/bin/env python3
"""Round 2 readout (PREREGISTRATION.md §5), applied mechanically.

    analyze2.py > runs/registered/readout.md
"""
import json, os, statistics
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
EXCLUDED = {("K1n59l4", 0, "rr"), ("K1n59l4", 3, "rr")}
NOT_RUN = {"rr": {"K0n53l3", "K1n57l3"}, "x4": {"K0n37l6"}}
FIT = {3: 8, 4: 9, 5: 10, 6: 11}

pred = {}
for l in open(os.path.join(HERE, "runs/registered/predictions.jsonl")):
    r = json.loads(l)
    pred[(r["cell"], r["draw"], r["arm"])] = r["predicted"]
recs = defaultdict(dict)
d = os.path.join(HERE, "runs/registered")
for f in os.listdir(d):
    if f.startswith("K") and f.endswith(".jsonl"):
        for l in open(os.path.join(d, f)):
            r = json.loads(l)
            recs[(r["cell"], r["draw"], r["arm"])].setdefault(r["dmax"], r)


def reading(k):
    by = recs.get(k, {})
    for D in sorted(by):
        r = by[D]
        if r["kind"] == "killed":
            return ("censored", D)
        if r.get("refuted_at", 0) > 0:
            return ("exact", r["refuted_at"])
    if by and max(by) >= pred[k] + 2:
        return ("censored", pred[k] + 3)
    return None


units = sorted(k for k in pred if k not in EXCLUDED and k[0] not in NOT_RUN[k[2]])
print("# Round 2 readout\n\n| cell | arm | readings | predicted | fit (`x4`) | median = predicted |\n|:--|:--|:--|:--|:--|:--|")
summary = {"x4": [0, 0, 0, 0], "rr": [0, 0, 0, 0]}  # resolved, exact, above, below
cellok = {"x4": True, "rr": True}
pending = []
for cell, arm in sorted({(c, a) for c, _, a in units}):
    vals, txt, p = [], [], None
    for (c, dr, a) in units:
        if (c, a) != (cell, arm):
            continue
        p = pred[(c, dr, a)]
        rd = reading((c, dr, a))
        if rd is None:
            txt.append("…"); pending.append(f"{c}:{dr}:{a}"); continue
        kind, v = rd
        if kind == "censored":
            txt.append(f"≥{v}"); continue
        txt.append(str(v)); vals.append(v)
        s = summary[arm]; s[0] += 1
        s[1] += v == p; s[2] += v > p; s[3] += v < p
    med = statistics.median(vals) if vals else None
    ok = med is not None and med == p
    cellok[arm] &= ok
    ell = int(cell.split("l")[1])
    print(f"| {cell} | {arm} | {' '.join(txt)} | {p} | {FIT[ell] if arm == 'x4' else ''} | {'yes' if ok else ('…' if med is None else 'no')} |")
for arm, name in (("x4", "A"), ("rr", "B")):
    n, ex, ab, be = summary[arm]
    rate = ex / n if n else 0
    print(f"\n**{name} ({arm})**: resolved {n}, equal to prediction {ex} ({rate:.0%}), above {ab}, below {be}. "
          f"Every cell's median equals the prediction: {cellok[arm]}. "
          f"Verdict: {'holds' if cellok[arm] and rate >= 0.9 else 'fails'}.")
    if arm == "x4":
        print(f"\n**A′**: no resolved `x4` draw above `d_aff`: {'holds' if ab == 0 else 'fails'}.")
if pending:
    print(f"\nUnfinished units ({len(pending)}): " + ", ".join(pending))

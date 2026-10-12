#!/usr/bin/env python3
"""Readout for round 1 (PREREGISTRATION.md §5), applied mechanically.

Sources: runs/registered (local lanes) and the completed remote lanes listed in
REMOTE.  runs/registered-remote-b is partial and is not read.  K1n61l4 is excluded
(seen before registration, §2).  One reading per (cell, draw, arm): the first refuting
degree, or censored ">=D" (scanned to the cap, killed, or refused).

    analyze.py > runs/registered/readout.md
"""
import json
import os
import statistics
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCES = ["registered", "registered-remote-c", "registered-remote-d", "registered-remote-e"]
EXCLUDED_CELLS = {"K1n61l4"}
FIT = {("x4", 4): 9, ("x4", 5): 10, ("x4", 6): 11, ("rr", 4): 6, ("rr", 5): 8, ("rr", 6): 8}
DISCRIMINATING = [("K0n61l6", "x4"), ("K1n61l6", "x4"), ("K1n29l5", "rr"), ("K1n45l5", "rr"),
                  ("K1n61l5", "rr"), ("K0n61l6", "rr"), ("K1n61l6", "rr")]
PLANNED = {("K1n23l6", "x4")}  # registered as not run


def load():
    pred = {}
    for line in open(os.path.join(HERE, "runs/registered/predictions.jsonl")):
        r = json.loads(line)
        pred[(r["cell"], r["draw"], r["arm"])] = r["d_reg"]
    lines = defaultdict(list)
    for src in SOURCES:
        d = os.path.join(HERE, "runs", src)
        if not os.path.isdir(d):
            continue
        for name in os.listdir(d):
            if name.startswith("K") and name.endswith(".jsonl"):
                for line in open(os.path.join(d, name)):
                    r = json.loads(line)
                    lines[(r["cell"], r["draw"], r["arm"])].append((src, r))
    return pred, lines


def reading(cell, draw, arm, recs, cap):
    """('exact', D) | ('censored', D) | None (unfinished)."""
    by_d = {}
    for _, r in recs:
        by_d.setdefault(r["dmax"], r)
    for d in sorted(by_d):
        r = by_d[d]
        if r["kind"] == "killed":
            return ("censored", d)
        if r.get("refuted_at", 0) > 0:
            return ("exact", r["refuted_at"])
    if by_d and max(by_d) >= cap and by_d[max(by_d)].get("refuted_at", 0) == 0:
        return ("censored", cap + 1)
    return None


def main():
    pred, lines = load()
    units = sorted({k for k in pred if k[0] not in EXCLUDED_CELLS and k[2] != "ctrl"
                    and (k[0], k[2]) not in PLANNED})
    res = {}
    for k in units:
        res[k] = reading(*k, lines.get(k, []), pred[k] + 2)
    print("# Round 1 readout\n")
    print("| cell | arm | readings | window | fit | in window |")
    print("|:--|:--|:--|:--|:--|:--|")
    cells = sorted({(c, a) for c, _, a in units})
    resolved = inside = below = above = 0
    medians = {}
    for c, a in cells:
        ell = int(c.split("l")[1])
        rs, ws = [], set()
        vals = []
        for (cc, d, aa) in units:
            if (cc, aa) != (c, a):
                continue
            dr = pred[(cc, d, aa)]
            ws.add((dr, dr + 1))
            r = res[(cc, d, aa)]
            if r is None:
                rs.append("…")
                continue
            kind, v = r
            if kind == "censored":
                rs.append(f"≥{v}")
                continue
            resolved += 1
            vals.append(v)
            if dr == 0:  # trace-constant draw: excluded from medians, as in earlier rounds
                rs.append("triv")
                inside += 1
                vals.pop()
                continue
            rs.append(str(v))
            if v < dr:
                below += 1
            elif v > dr + 1:
                above += 1
            else:
                inside += 1
        main_w = sorted(w for w in ws if w[0] > 0)
        med = statistics.median(vals) if vals else None
        medians[(c, a)] = (med, main_w[0] if main_w else None, len(vals))
        ok = "" if med is None else ("yes" if main_w and main_w[0][0] <= med <= main_w[0][1] else ("below" if med < main_w[0][0] else "above"))
        print(f"| {c} | {a} | {' '.join(rs)} | {main_w[0] if main_w else ''} | {FIT.get((a, ell), '')} | {ok} |")
    print(f"\nResolved draws: {resolved}; inside the window: {inside}; below: {below}; above: {above}.")
    p1 = inside / resolved if resolved else 0
    print(f"P1 (all resolved draws inside): {'holds' if below == above == 0 else 'fails'} ({p1:.0%} inside).")
    print("\n## Discriminating cells (P2)\n")
    p2 = True
    fitwins = 0
    for c, a in DISCRIMINATING:
        med, w, n = medians.get((c, a), (None, None, 0))
        ell = int(c.split("l")[1])
        fit = FIT[(a, ell)]
        inw = med is not None and n >= 3 and w[0] <= med <= w[1]
        p2 &= inw
        fitwins += med == fit
        print(f"- {c} {a}: median {med} over {n} draws, window {w}, fit {fit}: {'inside' if inside_ok(med, w) else 'outside'}")
    verdict = ("generic law holds out of sample" if p2 and p1 >= 0.9 else
               "fit holds, generic law refuted" if fitwins > len(DISCRIMINATING) / 2 else "neither")
    print(f"\nP2: {'holds' if p2 else 'fails'}.  Discriminating cells whose median equals the fit: {fitwins} of {len(DISCRIMINATING)}.")
    print(f"\n**Registered verdict: {verdict}.**")
    unfinished = [k for k, v in res.items() if v is None]
    if unfinished:
        print(f"\nUnfinished units ({len(unfinished)}): " + ", ".join(f"{c}:{d}:{a}" for c, d, a in unfinished))


def inside_ok(med, w):
    return med is not None and w is not None and w[0] <= med <= w[1]


if __name__ == "__main__":
    main()

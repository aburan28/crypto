#!/usr/bin/env python3
"""Score PREREGISTRATION.md hypotheses H1, H1b, H2a, H2b from committed JSONL.

    python3 research/fall_degree_bounds_20261008/score.py
"""
import glob, json, os
from collections import defaultdict
from statistics import median

HERE = os.path.dirname(os.path.abspath(__file__))
EXCLUDED_SYZ_SEEDS = {60000, 60001, 60002}  # seen before registration


def load(pattern):
    out = []
    for f in sorted(glob.glob(os.path.join(HERE, pattern))):
        out += [json.loads(l) for l in open(f) if l.strip()]
    return out


def verdict(ok, testable=True):
    return "not testable" if not testable else ("held" if ok else "FALSIFIED")


# ── E1 ────────────────────────────────────────────────────────────────
e1 = [r for r in load("exp/e1_*.jsonl") if "predicted_L" in r]  # as-run records
print(f"E1 as run (constant Tr(x3), see Amendment 1): {len(e1)} draws")
bad_span = [r for r in e1 if not r["L_in_span"]]
c1 = [r for r in e1 if not r["trace_on_V_zero"]]
bad_c1 = [r for r in c1 if r["ffd"] != 2]
c2 = [r for r in e1 if r["trace_on_V_zero"] and r["predicted_L"] == "1"]
bad_c2 = [r for r in c2 if not (r["unsat"] and r["d_last"] is not None and r["d_last"] <= 2)]
c3 = [r for r in e1 if r["trace_on_V_zero"] and r["predicted_L"] == "0"]
side = [r for r in c1 if r["n"] >= 12]
bad_side = [r for r in side if r["deg_le1_dim"] != 1]
print(f"  H1  L in span(F):                {len(e1)-len(bad_span)}/{len(e1)}  -> {verdict(not bad_span)}")
print(f"  C1  Tr|V != 0 => d_ff = 2:        {len(c1)-len(bad_c1)}/{len(c1)}  -> {verdict(not bad_c1, bool(c1))}")
print(f"  C2  Tr|V = 0, L = 1 => refuted@2: {len(c2)-len(bad_c2)}/{len(c2)}  -> {verdict(not bad_c2, bool(c2))}")
print(f"  C3  Tr|V = 0, L = 0 cells: {len(c3)}; d_ff there: "
      f"{sorted(r['ffd'] for r in c3)}; sat {sum(not r['unsat'] for r in c3)}")
print(f"  side dim(span∩R<=1) = 1, n>=12:  {len(side)-len(bad_side)}/{len(side)}  -> {verdict(not bad_side, bool(side))}")
for r in bad_span + bad_c1 + bad_c2 + bad_side:
    print("    decisive:", r["n"], r["family"], r["seed"], r["deg_le1"], r["ffd"], r["d_last"])

# ── E1, Amendment 1: registered constant Tr(b/x3^2) ──────────────────
a1 = load("exp/e1_rescored_amendment1.jsonl")
bad_a = [r for r in a1 if not r["L_in_span"]]
c2a = [r for r in a1 if r["trace_on_V_zero"] and r["c_registered"] == 1]
bad_c2a = [r for r in c2a if not (r["one_in_span"] and r["unsat"] and r["d_last"] is not None and r["d_last"] <= 2)]
c3a = [r for r in a1 if r["family"] == "random" and r["trace_on_V_zero"] and r["c_registered"] == 0]
print(f"E1 Amendment 1 (registered constant): {len(a1)} draws")
print(f"  H1  L in span(F):                {len(a1)-len(bad_a)}/{len(a1)}  -> {verdict(not bad_a)}")
print(f"  C2  Tr|V = 0, Tr(b/x3^2) = 1 => refuted at degree 2: {len(c2a)-len(bad_c2a)}/{len(c2a)}  -> {verdict(not bad_c2a, bool(c2a))}")
print(f"  C3  random V, Tr|V = 0, constant 0: d_ff {sorted(r['ffd'] for r in c3a)}, degree<=1 dims {sorted(r['deg_le1_dim'] for r in c3a)}")
for r in bad_a + bad_c2a:
    print("    decisive:", r)

# ── E1b ───────────────────────────────────────────────────────────────
e1b = [r for r in load("exp/e1b_*.jsonl") if r["seed"] not in EXCLUDED_SYZ_SEEDS]
scored = [r for r in e1b if r["family"] == "random" and not r["trace_on_V_zero"] and r["n"] >= 12]
bad = [r for r in scored if r["syz_dim"] != 1 or not r["syzygies"][0]["multiplier_equals_L_plus_1"]
       or not r["syzygies"][0]["F_S_equals_L"]]
print(f"E1b: {len(e1b)} draws, {len(scored)} scored")
print(f"  H1b unique syzygy = (L+1)·L:      {len(scored)-len(bad)}/{len(scored)}  -> {verdict(not bad, bool(scored))}")
for r in bad:
    print("    decisive:", r["n"], r["seed"], r["syz_dim"])
other = defaultdict(list)
for r in e1b:
    if r not in scored:
        other[(r["family"], r["trace_on_V_zero"], r["predicted_L"])].append(r["syz_dim"])
for k, v in sorted(other.items()):
    print("  unscored", k, "syz_dim:", sorted(v))

# ── E2 ────────────────────────────────────────────────────────────────
ctrl = load("exp/e2_*.jsonl")
sem = load("runs/*.jsonl")
def med(rows):
    vals = [r["d_last"] for r in rows if r["d_last"] is not None]
    return median(vals) if vals else None
print(f"E2: {len(ctrl)} control draws")
print("  n | Semaev sat | Semaev unsat | plain sat | plain unsat | planted sat | planted unsat")
gaps = {"sat": [], "unsat": []}
h2a_ok = True
for n in sorted({r["n"] for r in ctrl}):
    row = [n]
    s_cells = {}
    for st in ("sat", "unsat"):
        s_cells[st] = med([r for r in sem if r["n"] == n and r["family"] == "random" and (not r["unsat"]) == (st == "sat")])
    c_cells = {}
    for kind in ("plain", "planted_linear"):
        for st in ("sat", "unsat"):
            c_cells[(kind, st)] = med([r for r in ctrl if r["n"] == n and r["kind"] == kind and (not r["unsat"]) == (st == "sat")])
    print(f"  {n} | {s_cells['sat']} | {s_cells['unsat']} | {c_cells[('plain','sat')]} | {c_cells[('plain','unsat')]} "
          f"| {c_cells[('planted_linear','sat')]} | {c_cells[('planted_linear','unsat')]}")
    for st in ("sat", "unsat"):
        a, b = s_cells[st], c_cells[("planted_linear", st)]
        if a is not None and b is not None:
            gaps[st].append((n, b - a))
            if a > b:
                h2a_ok = False
print(f"  H2a Semaev d_F <= planted control (medians): -> {verdict(h2a_ok, any(gaps.values()))}")
h2b_ok = all(g2 >= g1 for st in gaps for (_, g1), (_, g2) in zip(gaps[st], gaps[st][1:]))
print(f"  H2b gap nondecreasing in n: gaps {dict(gaps)} -> {verdict(h2b_ok, any(len(g) > 1 for g in gaps.values()))}")

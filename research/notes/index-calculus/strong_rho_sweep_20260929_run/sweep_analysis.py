#!/usr/bin/env python3
"""Apply the pre-registered sweep analysis (RESEARCH_STRONG_RHO_SWEEP_PROTOCOL_20260929.md).

Reads every cell_n*/callgrind_{rho,ic}_result.json (+ native_result.json), the
ladder cell (n=53, L=1024) from ../strong_rho_ladder_20260929_run, and prints:
  * the per-cell table (Ir of both arms, IC/rho, native wall for reference),
  * OLS of log2(Ir) on log2(r) for each arm over the size sweep (L=1024, n in
    {37,41,53,61}) and of log2(IC/rho) on log2(r) (= the slope difference),
    with standard errors and the Student-t 95% interval (df = cells - 2),
  * the pre-registered decision.
The subgroup order r of each cell is the retained value (base headers of the
PR #830 evidence); it is cross-checked against r = ideal^2 * 2A / pi computed from
the rho records themselves.
"""
import glob
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LADDER = os.path.join(HERE, "..", "strong_rho_ladder_20260929_run")

# subgroup orders from the retained base headers (n=37,41,53) and PR #830's
# RESULT.md (n=61); each is re-derived below from the rho records.
ORDER = {37: 230603167, 41: 549756390943, 53: 21044858204113, 61: 162888033982417}
T975 = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365, 8: 2.306}


def ir_from_stderr(path):
    for line in open(path):
        if "I   refs:" in line:
            return int(line.split("refs:")[1].strip().replace(",", ""))
    return None


def order_from_records(path):
    with open(path) as fh:
        rec = json.loads(fh.readline())
    ideal, a = rec["ideal_independent_steps"], rec["automorphism_size"]
    return ideal * ideal * 2 * a / math.pi


def ols(x, y):
    n = len(x)
    xm, ym = sum(x) / n, sum(y) / n
    sxx = sum((xi - xm) ** 2 for xi in x)
    sxy = sum((xi - xm) * (yi - ym) for xi, yi in zip(x, y))
    slope = sxy / sxx
    icpt = ym - slope * xm
    ssr = sum((yi - (icpt + slope * xi)) ** 2 for xi, yi in zip(x, y))
    df = n - 2
    se = math.sqrt(ssr / df / sxx) if df > 0 else float("nan")
    t = T975.get(df, 1.96)
    return {"slope": slope, "se": se, "ci95": [slope - t * se, slope + t * se], "df": df, "resid_rms": math.sqrt(ssr / n)}


def cell_record(n, L, K, dirpath, ladder=False):
    if ladder:
        ic_ir = ir_from_stderr(os.path.join(dirpath, "callgrind_ic.stderr.log"))
        rho_ir = ir_from_stderr(os.path.join(dirpath, "callgrind_rung3.stderr.log"))
        rho_native = os.path.join(dirpath, "rung3_native.jsonl")
        rho_wall = None
        for line in open(rho_native):
            o = json.loads(line)
            if o["kind"] == "rho_ks_batch_summary":
                rho_wall = o["in_process_ms"] / 1000
        ic_wall = None
        p = os.path.join(dirpath, "ic_native.summary.json")
        if os.path.exists(p):
            ic_wall = json.load(open(p))["timing_ms"]["process_total"] / 1000
    else:
        ic_ir = rho_ir = None
        p = os.path.join(dirpath, "callgrind_ic_result.json")
        if os.path.exists(p):
            r = json.load(open(p))["ic"]
            ic_ir = r["ir"] if r["gates_G1"] else None
        p = os.path.join(dirpath, "callgrind_rho_result.json")
        if os.path.exists(p):
            r = json.load(open(p))["rho"]
            rho_ir = r["ir"] if r["gates_G1"] else None
        rho_native = os.path.join(dirpath, "rho_native.jsonl")
        nat = json.load(open(os.path.join(dirpath, "native_result.json")))
        rho_wall, ic_wall = nat["rho"]["in_process_ms"] / 1000, nat["ic"]["timing_ms"]["process_total"] / 1000
    order = ORDER[n]
    derived = order_from_records(rho_native)
    assert abs(derived - order) / order < 1e-6, f"order mismatch n={n}: {derived} vs {order}"
    return {"n": n, "L": L, "K": K, "r": order, "ic_ir": ic_ir, "rho_ir": rho_ir,
            "ic_over_rho": (ic_ir / rho_ir) if ic_ir and rho_ir else None,
            "ic_wall_s": ic_wall, "rho_wall_s": rho_wall}


def main():
    cells = [cell_record(53, 1024, 440, LADDER, ladder=True)]
    for line in open(os.path.join(HERE, "cells.txt")):
        n, L, K, _ = line.split()
        cells.append(cell_record(int(n), int(L), int(K), os.path.join(HERE, f"cell_n{n}_L{L}_K{K}")))
    cells.sort(key=lambda c: (c["L"], c["n"]))
    print("| n | L | K | log2 r | IC Ir | rho R3 Ir | IC / rho | IC wall s | rho wall s |")
    print("|--:|--:|--:|--:|--:|--:|--:|--:|--:|")
    for c in cells:
        f = lambda v, fmt: (fmt.format(v) if v is not None else "n/a")
        print(f"| {c['n']} | {c['L']:,} | {c['K']} | {math.log2(c['r']):.2f} | {f(c['ic_ir'], '{:,}')} | "
              f"{f(c['rho_ir'], '{:,}')} | {f(c['ic_over_rho'], '{:.3f}')} | {f(c['ic_wall_s'], '{:.1f}')} | {f(c['rho_wall_s'], '{:.1f}')} |")

    sizes = [c for c in cells if c["L"] == 1024 and c["ic_over_rho"] is not None]
    out = {"cells": cells}
    if len(sizes) >= 4:
        x = [math.log2(c["r"]) for c in sizes]
        fit_ic = ols(x, [math.log2(c["ic_ir"]) for c in sizes])
        fit_rho = ols(x, [math.log2(c["rho_ir"]) for c in sizes])
        fit_diff = ols(x, [math.log2(c["ic_over_rho"]) for c in sizes])
        out.update(fit_ic=fit_ic, fit_rho=fit_rho, fit_diff=fit_diff)
        print()
        for name, fit in (("IC", fit_ic), ("rho R3", fit_rho), ("log2(IC/rho) = slope difference", fit_diff)):
            print(f"slope {name}: {fit['slope']:.3f} +/- {fit['se']:.3f} (95% CI {fit['ci95'][0]:.3f}..{fit['ci95'][1]:.3f}, df={fit['df']})")
        d = fit_diff
        cond_i = any(c["ic_over_rho"] < 0.8 for c in cells if c["ic_over_rho"] is not None)
        cond_ii = d["slope"] <= -0.05 and d["ci95"][1] < 0
        out.update(cond_i_any_cell_below_0p8=cond_i, cond_ii_slope_gap=cond_ii)
        if cond_i or cond_ii:
            verdict = "SURVIVES as a scaling result"
        elif d["ci95"][0] <= 0 <= d["ci95"][1]:
            verdict = "DIES (no cell below 0.8, no resolved favourable slope gap); slope difference not resolved from 0"
        else:
            verdict = "DIES (no cell below 0.8); slope difference resolved but not favourable to IC"
        out["verdict"] = verdict
        print(f"\n(i) any cell IC/rho < 0.8: {cond_i};  (ii) slope gap <= -0.05 with CI excluding 0: {cond_ii}\nDecision: {verdict}")
    else:
        print(f"\nonly {len(sizes)} size cells complete; fit not yet available")
    json.dump(out, open(os.path.join(HERE, "sweep_analysis.json"), "w"), indent=2)


if __name__ == "__main__":
    sys.exit(main())

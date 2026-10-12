"""Build the tables in RESULTS.md from the frozen results/*.json files.

    python3 report.py > RESULTS.md

Every number printed here is read from a results file; nothing is computed
from formulas except the per-row ratios S = ops/√p and ops/√(p/d) that the
runs themselves already stored, and the medians over seeds.
"""

from __future__ import annotations

import glob
import json
import math
import os
import statistics

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")

VARIANT_ORDER = [
    "rho_reference",
    "plain_bsgs_as_in_rust",
    "cheon_p-1_bsgs_naive",
    "cheon_p-1_bsgs_comb",
    "cheon_p-1_kangaroo_comb",
    "cheon_p+1_bsgs_comb",
]

LABEL = {
    "rho_reference": "Pollard rho (reference, sees G and [α]G only)",
    "plain_bsgs_as_in_rust": "plain BSGS, d baby steps (old `cheon_attack.rs`)",
    "cheon_p-1_bsgs_naive": "Cheon p−1, BSGS, naive scalar mults",
    "cheon_p-1_bsgs_comb": "Cheon p−1, BSGS, fixed-base comb",
    "cheon_p-1_kangaroo_comb": "Cheon p−1, kangaroo, fixed-base comb",
    "cheon_p+1_bsgs_comb": "Cheon p+1, BSGS, fixed-base comb",
}

CLASS = {
    "rho_reference": "reference",
    "plain_bsgs_as_in_rust": "accounting",
    "cheon_p-1_bsgs_naive": "advance (DLPwAI floor, not ECDLP)",
    "cheon_p-1_bsgs_comb": "engineering",
    "cheon_p-1_kangaroo_comb": "engineering (memory)",
    "cheon_p+1_bsgs_comb": "advance (DLPwAI floor, not ECDLP)",
}


def load(case: str):
    files = sorted(glob.glob(os.path.join(RES, f"{case}_*.json")))
    return [json.load(open(f)) for f in files]


def med(xs):
    return statistics.median(xs) if xs else float("nan")


def fmt(x, nd=2):
    if isinstance(x, float) and (math.isnan(x)):
        return "—"
    if isinstance(x, float):
        return f"{x:,.{nd}f}"
    return f"{x:,}"


def main_table(case: str):
    data = load(case)
    if not data:
        return
    print(f"\n### {case} case: every variant, one unit, medians over seeds\n")
    print("| bits | p | d | log_p d | variant | seeds | correct | ops (median) | S = ops/√p | ops/√(p/d) | S / S_rho | class |")
    print("|---:|---:|---:|---:|:--|---:|:--|---:|---:|---:|---:|:--|")
    for run in data:
        p, d = run["curve"]["p"], run["d"]
        rows = run["rows"]
        by = {}
        for r in rows:
            by.setdefault(r["variant"], []).append(r)
        s_rho = med([r["S"] for r in by.get("rho_reference", [])])
        for v in VARIANT_ORDER:
            if v not in by:
                continue
            rs = by[v]
            ok = sum(1 for r in rs if r["correct"])
            S = med([r["S"] for r in rs])
            fl = med([r["ops_over_floor"] for r in rs])
            ops = med([r["ops"] for r in rs])
            ratio = S / s_rho if s_rho == s_rho and s_rho > 0 else float("nan")
            print(f"| {run['bits']} | {p} | {d} | {run['log_p_d']:.3f} | {LABEL[v]} | {len(rs)} | "
                  f"{ok}/{len(rs)} | {fmt(int(ops))} | {fmt(S, 3)} | {fmt(fl, 1)} | "
                  f"{fmt(ratio, 3)} | {CLASS[v]} |")


def phase_table(case: str, variant: str):
    data = load(case)
    if not data:
        return
    print(f"\n### {case}: where `{variant}` spends its operations (median over seeds)\n")
    phases = []
    for run in data:
        for r in run["rows"]:
            if r["variant"] == variant:
                for k in r["phases"]:
                    if k not in phases:
                        phases.append(k)
    print("| bits | " + " | ".join(phases) + " | total |")
    print("|---:|" + "---:|" * (len(phases) + 1))
    for run in data:
        rs = [r for r in run["rows"] if r["variant"] == variant]
        if not rs:
            continue
        cells = []
        for k in phases:
            cells.append(fmt(int(med([r["phases"].get(k, 0) for r in rs]))))
        print(f"| {run['bits']} | " + " | ".join(cells) + f" | {fmt(int(med([r['ops'] for r in rs])))} |")


def exponent_fit(case: str, variant: str):
    """Least-squares slope of log2(ops) against log2(p): the measured exponent."""
    data = load(case)
    xs, ys = [], []
    for run in data:
        rs = [r for r in run["rows"] if r["variant"] == variant and r["correct"]]
        if len(rs) == 0:
            continue
        xs.append(math.log2(run["curve"]["p"]))
        ys.append(math.log2(med([r["ops"] for r in rs])))
    if len(xs) < 3:
        return None
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    slope = sxy / sxx
    resid = [y - (my + slope * (x - mx)) for x, y in zip(xs, ys)]
    se = math.sqrt(sum(r * r for r in resid) / max(1, n - 2) / sxx)
    return slope, se, n


def torsion_table():
    files = sorted(glob.glob(os.path.join(RES, "torsion_*.json")))
    if not files:
        return
    print("\n### Torsion-point route (Kim–Cheon elliptic embedding), honest and oracle versions\n")
    print("| bits | p | d | δ = #E'/d | aux inputs needed (d²) | m | variant | correct / collision | ops | S = ops/√p | ops/√(p/d) | pairs tested |")
    print("|---:|---:|---:|---:|---:|---:|:--|:--|---:|---:|---:|---:|")
    for f in files:
        for run in json.load(open(f)):
            for r in run["rows"]:
                flag = ("yes" if r.get("correct") else "no") if "correct" in r else (
                    "collision" if r.get("collision") else "none")
                print(f"| {run['bits']} | {run['p']} | {run['d']} | {run['delta']} | {run['aux_inputs_needed']} | {run['m']} | "
                      f"{r['variant']} | {flag} | {fmt(r['ops'])} | {fmt(r['S'], 1)} | {fmt(r['ops_over_floor'], 1)} | "
                      f"{fmt(r.get('pairs_tested', r.get('samples_to_collision', 0)))} |")


def main():
    print("# Results: auxiliary inputs and torsion points against rho\n")
    print("Generated by `report.py` from `results/*.json`; do not edit by hand.\n")
    for case in ("pminus1", "pplus1"):
        main_table(case)
    print("\n### Fitted exponents, log2(ops) against log2(p), median ops per size\n")
    print("| case | variant | slope | s.e. | sizes | prediction |")
    print("|:--|:--|---:|---:|---:|:--|")
    preds = {
        ("pminus1", "rho_reference"): "1/2",
        ("pminus1", "cheon_p-1_bsgs_naive"): "1/4 + log factor (d ≈ p^{1/2})",
        ("pminus1", "cheon_p-1_bsgs_comb"): "1/4 + smaller log factor",
        ("pminus1", "cheon_p-1_kangaroo_comb"): "1/4 + log factor",
        ("pplus1", "rho_reference"): "1/2",
        ("pplus1", "cheon_p+1_bsgs_comb"): "1/3 + log factor (d ≈ p^{1/3})",
    }
    for (case, v), pred in preds.items():
        fit = exponent_fit(case, v)
        if fit:
            print(f"| {case} | {LABEL[v]} | {fit[0]:.3f} | {fit[1]:.3f} | {fit[2]} | {pred} |")
    phase_table("pminus1", "cheon_p-1_bsgs_comb")
    phase_table("pplus1", "cheon_p+1_bsgs_comb")
    torsion_table()


if __name__ == "__main__":
    main()

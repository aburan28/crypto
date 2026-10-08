"""Round eight analysis: reduced factor base with single large primes.

Usage: python3 analyse.py DIR   (reads DIR/L*.txt, DIR/bench.jsonl)
Writes DIR/summary.json and prints the table and the verdict on the
four hypotheses registered in the note before measuring.
"""
import glob
import json
import math
import os
import statistics as st
import sys

# Solve growth exponents s from round seven (registered in the protocol).
S_EXP = {2: 1.37, 3: 1.70, 4: 1.64}
IDEAL_WALK = 1.25331


def load(d):
    sweeps = {}
    for f in sorted(glob.glob(os.path.join(d, "L[0-9].txt"))):
        rows = [json.loads(l[5:]) for l in open(f) if l.startswith("JSON ")]
        sweeps[os.path.basename(f)[:-4]] = {
            (r["genus"], r["p"], r["theta"], r["large_primes"]): r for r in rows
        }
    return sweeps


def comb(P, B):
    return P - B * (1 - math.exp(-P / B)) if B > 0 else 0.0


def trials_needed(theta, g, m, sig, lp):
    R = theta * m + 1 + 8
    full = sig * theta ** g
    part = sig * g * theta ** (g - 1) * (1 - theta)
    lo, hi = 0.0, 1e15
    for _ in range(300):
        T = (lo + hi) / 2
        got = full * T + (comb(part * T, (1 - theta) * m) if lp else 0.0)
        lo, hi = (T, hi) if got < R else (lo, T)
    return hi


def model_ratio(base, theta, lp):
    """S(theta)/S(1) predicted from the theta = 1 row's own measurements."""
    g, m = base["genus"], base["m_full"]
    rel, pre = base["relation_ops"], base["precompute_ops"]
    orc, la = base["oracle_equiv"], base["la_equiv"]
    total = rel + orc + la + base["partial_equiv"]
    walked = max(rel - pre, 1.0)
    per_trial = (walked + orc) / walked
    t1 = trials_needed(1.0, g, m, base["smoothness"], True)
    tt = trials_needed(theta, g, m, base["smoothness"], lp)
    cost = per_trial * walked * tt / t1 + pre + la * theta ** S_EXP[g]
    return cost / total


def main(d):
    sweeps = load(d)
    labels = sorted(sweeps)
    bench = {}
    if os.path.exists(os.path.join(d, "bench.jsonl")):
        for line in open(os.path.join(d, "bench.jsonl")):
            rec = json.loads(line)
            bench[rec["label"]] = rec["run"]
    keys = sorted(set().union(*[set(s) for s in sweeps.values()]))
    complete = [k for k in keys if all(k in sweeps[l] for l in labels)]

    # Deterministic fields must agree across sweeps (only wall-derived
    # quantities may differ).
    det = ("N", "c", "m_full", "fb", "relation_ops", "precompute_ops", "oracle_modmuls",
           "la_modmuls", "partials", "combined", "rho_ops", "ic_correct", "rho_correct")
    drift = [(k, f) for k in complete for f in det
             if len({json.dumps(sweeps[l][k][f]) for l in labels}) > 1]

    med = {}
    for k in complete:
        rs = [sweeps[l][k] for l in labels]
        r = dict(rs[0])
        for f in ("S_ic", "S_rho", "S_walk", "ratio", "oracle_equiv", "la_equiv",
                  "partial_equiv", "wall_ic_ms", "wall_rho_ms", "floor_S"):
            r[f] = st.median(x[f] for x in rs)
            r[f + "_all"] = [x[f] for x in rs]
        r["verified_all"] = all(x["ic_correct"] and x["rho_correct"] for x in rs)
        r["ratio_corr"] = r["S_ic"] / (r["S_rho"] - r["S_walk"] + IDEAL_WALK)
        med[k] = r

    rows, cells = [], []
    instances = sorted({(k[0], k[1]) for k in med})
    print(f"sweeps {labels}; contended: {[l for l in labels if bench.get(l, {}).get('contended')]}; "
          f"deterministic drift across sweeps: {len(drift)}")
    print(f"{'g':>2}{'p':>6}{'N':>11}{'m':>6} {'config':>9} {'fb':>5} {'S_ic':>7} {'S/S1':>6} "
          f"{'model':>6} {'ratio':>6} {'corr':>6} {'solve':>6} {'part':>6} {'comb':>6} {'ok':>3}")
    for (g, p) in instances:
        base = med.get((g, p, 1.0, True))
        if base is None:
            continue
        for k in sorted((k for k in med if k[:2] == (g, p)), key=lambda k: (-k[2], not k[3])):
            r = med[k]
            rel = r["S_ic"] / base["S_ic"]
            pred = model_ratio(base, k[2], k[3]) if k[2] < 1 else 1.0
            row = dict(genus=g, p=p, N=r["N"], m_full=r["m_full"], theta=k[2], large_primes=k[3],
                       fb=r["fb"], S_ic=r["S_ic"], S_ic_all=r["S_ic_all"], S_rho=r["S_rho"],
                       S_walk=r["S_walk"], ratio=r["ratio"], ratio_corr=r["ratio_corr"],
                       rel_to_theta1=rel, model_pred=pred, la_equiv=r["la_equiv"],
                       oracle_equiv=r["oracle_equiv"], partial_equiv=r["partial_equiv"],
                       relation_ops=r["relation_ops"], partials=r["partials"],
                       combined=r["combined"], smoothness=r["smoothness"],
                       wall_ic_ms=r["wall_ic_ms"], wall_rho_ms=r["wall_rho_ms"],
                       verified_all=r["verified_all"])
            rows.append(row)
            if k[2] < 1:
                cells.append(row)
            cfg = f"{k[2]:.2f}{'+LP' if k[3] else '   '}"
            print(f"{g:>2}{p:>6}{r['N']:>11}{r['m_full']:>6} {cfg:>9} {r['fb']:>5} {r['S_ic']:>7.3f} "
                  f"{rel:>6.3f} {pred:>6.3f} {r['ratio']:>6.3f} {r['ratio_corr']:>6.3f} "
                  f"{r['la_equiv']:>6.0f} {r['partials']:>6.0f} {r['combined']:>6.0f} "
                  f"{'y' if r['verified_all'] else 'N':>3}")

    def best_lt1(g, p):
        c = [x for x in cells if (x["genus"], x["p"]) == (g, p)]
        return min(c, key=lambda x: x["rel_to_theta1"]) if c else None

    h1_rows = [best_lt1(g, p) for (g, p) in instances if g in (3, 4)]
    h1_viol = [(x["genus"], x["p"], x["theta"], x["large_primes"], round(x["rel_to_theta1"], 3))
               for x in h1_rows if x and x["rel_to_theta1"] < 0.90]
    h2_small = [best_lt1(2, p) for (g, p) in instances if g == 2 and p <= 251]
    h2_large = [best_lt1(2, p) for (g, p) in instances if g == 2 and p >= 2003]
    h2_small_viol = [(x["p"], round(x["rel_to_theta1"], 3)) for x in h2_small if x and x["rel_to_theta1"] < 0.95]
    h2_large_fail = [(x["p"], round(x["rel_to_theta1"], 3)) for x in h2_large if x and x["rel_to_theta1"] > 0.95]
    within = [abs(x["rel_to_theta1"] / x["model_pred"] - 1) <= 0.25 for x in cells]
    within_lp = [abs(x["rel_to_theta1"] / x["model_pred"] - 1) <= 0.25 for x in cells if x["large_primes"]]
    h4 = []
    for (g, p) in instances:
        a, b = med.get((g, p, 0.5, True)), med.get((g, p, 0.5, False))
        if a and b:
            h4.append(((g, p), a["S_ic"] < b["S_ic"], round(a["S_ic"] / b["S_ic"], 3)))
    verdict = {
        "H1_no_win_at_genus_3_4": {"holds": not h1_viol, "violations": h1_viol,
                                   "best_per_row": [(x["genus"], x["p"], x["theta"], x["large_primes"], round(x["rel_to_theta1"], 3)) for x in h1_rows if x]},
        "H2_genus2_crossover": {"holds": not h2_small_viol and not h2_large_fail and bool(h2_large),
                                "small_p_violations": h2_small_viol, "large_p_failures": h2_large_fail,
                                "large_p_best": [(x["p"], x["theta"], x["large_primes"], round(x["rel_to_theta1"], 3)) for x in h2_large if x]},
        "H3_model_fidelity": {"holds": sum(within) >= 0.8 * len(within), "cells": len(within),
                              "within_25pct": sum(within), "fraction": round(sum(within) / max(len(within), 1), 3),
                              "lp_only_fraction": round(sum(within_lp) / max(len(within_lp), 1), 3)},
        "H4_large_primes_beat_control": {"holds": all(h for _, h, _ in h4), "rows": h4},
        "all_verified": all(r["verified_all"] for r in rows),
    }
    print(json.dumps(verdict, indent=1))
    json.dump({"sweeps": {l: {k: bench.get(l, {}).get(k) for k in ("wall_seconds", "contended", "other_cpu_seconds")} for l in labels},
               "deterministic_drift": drift, "rows": rows, "verdict": verdict},
              open(os.path.join(d, "summary.json"), "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1])

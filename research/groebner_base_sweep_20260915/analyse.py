#!/usr/bin/env python3
"""Rank invariant factor bases by measured stage cost, over independent
target samples.

    python3 analyse.py sample-a sample-b [--output summary.json]

Per curve and summand count, one row per candidate invariant subspace:
dimension, size, relation columns, whether it lies inside `ker Tr`,
measured coverage with a Wilson 95% interval, the expected trials
`T = (U+1)/p` that the factor-base search minimises today, and the
expected stage word XORs `C = T · ops/target` that collection actually
spends.  The samples are different target draws and are pooled into one
binomial, so a coverage difference inside the interval is reported as
what it is — sampling error, not a property of the base.
"""
import argparse, json, pathlib, sys

INF = float("inf")

def wilson(successes, trials, z=1.96):
    if trials == 0:
        return 0.0, 0.0, 0.0
    p = successes / trials
    d = 1 + z * z / trials
    centre = (p + z * z / (2 * trials)) / d
    half = z * ((p * (1 - p) / trials + z * z / (4 * trials ** 2)) ** 0.5) / d
    return p, max(centre - half, 0.0), min(centre + half, 1.0)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("samples", nargs="+")
    ap.add_argument("--output")
    args = ap.parse_args()
    docs = [json.loads(pathlib.Path(s, "base_sweep.json").read_text()) for s in args.samples]
    if len({d["target_seed"] for d in docs}) != len(docs):
        sys.exit("samples must use distinct target seeds")

    out_rows = []
    for index, rung in enumerate(docs[0]["rungs"]):
        rungs = [d["rungs"][index] for d in docs]
        assert len({(r["curve"], r["m"]) for r in rungs}) == 1
        by_key = [{tuple(c["divisor"]): c for c in r["candidates"]} for r in rungs]
        shared = set(by_key[0])
        for m in by_key[1:]:
            shared &= set(m)

        rows = []
        for k in sorted(shared):
            cands = [m[k] for m in by_key]
            if not all(c.get("admissible") for c in cands):
                continue
            decided = sum(c["decided"] for c in cands)
            decomposed = sum(c["decomposed"] for c in cands)
            ops = sum(c["ops_per_target"] * c["decided"] for c in cands) / max(decided, 1)
            p, lo, hi = wilson(decomposed, decided)
            u = cands[0]["unknowns"]
            cost = lambda q: (u + 1) / q * ops if q > 0 else INF
            rows.append({
                "divisor": list(k), "ell": cands[0]["ell"], "points": cands[0]["points"],
                "unknowns": u,
                # V lies inside ker Tr exactly when (x+1) does not divide
                # the divisor: Tr = h(sigma) with h = (t^n+1)/(t+1), so
                # ker g(sigma) is inside ker Tr iff g divides h.
                "trace_zero": 0 not in k,
                "is_pipeline_default": cands[0].get("is_pipeline_default", False),
                "decided": decided, "decomposed": decomposed,
                "coverage": p, "coverage_ci": [lo, hi],
                "per_sample_coverage": [c["coverage"] for c in cands],
                "ops_per_target": ops,
                "T": (u + 1) / p if p > 0 else INF,
                "C": cost(p), "C_range": [cost(hi), cost(lo)],
            })
        if not rows:
            continue
        finite = [r for r in rows if r["C"] < INF]
        if not finite:
            continue
        bestC = min(r["C"] for r in finite)
        bestT = min(r["T"] for r in finite)
        rows.sort(key=lambda r: r["C"])
        print(f"\n=== {rungs[0]['curve']}  m = {rungs[0]['m']} "
              f"({sum(r['decided'] for r in rows[:1])} targets per base, pooled) ===")
        print(f"{'divisor':<10} {'l':>3} {'|F|':>6} {'U':>5} {'Tr0':>4} {'decomposed':>11} "
              f"{'coverage 95% CI':>19} {'ops/target':>11} {'T':>8} {'C':>10} "
              f"{'C 95% band':>23} {'C/best':>7} flags")
        for r in rows:
            flags = ("D" if r["is_pipeline_default"] else "-")
            flags += "T" if r["T"] <= bestT * 1.000001 else "-"
            flags += "C" if r["C"] <= bestC * 1.000001 else "-"
            ratio = f"{r['C']/bestC:7.2f}" if r["C"] < INF else "      -"
            band = (f"{r['C_range'][0]:10.3g}-{r['C_range'][1]:10.3g}"
                    if r["C"] < INF else "                      -")
            print(f"{'.'.join(map(str, r['divisor'])):<10} {r['ell']:>3} {r['points']:>6} "
                  f"{r['unknowns']:>5} {str(r['trace_zero'])[0]:>4} "
                  f"{r['decomposed']:>5}/{r['decided']:<5} "
                  f"{r['coverage']*100:>6.0f}% ({r['coverage_ci'][0]*100:3.0f}-{r['coverage_ci'][1]*100:3.0f}%) "
                  f"{r['ops_per_target']:>11.3g} {r['T']:>8.4g} {r['C']:>10.3g} "
                  f"{band:>23} {ratio} {flags}")
        t_choice = min(finite, key=lambda r: r["T"])
        c_choice = min(finite, key=lambda r: r["C"])
        penalty = t_choice["C"] / c_choice["C"]
        # Conservative: the trials pick at the cheap end of its interval
        # against the cost pick at the dear end of its own.
        worst = t_choice["C_range"][0] / c_choice["C_range"][1]
        print(f"  selecting on trials picks {'.'.join(map(str, t_choice['divisor']))}: "
              f"{penalty:.2f}x the cost-optimal {'.'.join(map(str, c_choice['divisor']))}"
              f" ({worst:.2f}x at the worst end of both coverage intervals)")
        out_rows.append({
            "curve": rungs[0]["curve"], "m": rungs[0]["m"], "candidates": rows,
            "trials_optimal": t_choice["divisor"], "cost_optimal": c_choice["divisor"],
            "penalty_of_selecting_on_trials": penalty, "penalty_lower_bound": worst,
        })

    if args.output:
        pathlib.Path(args.output).write_text(json.dumps({
            "unit": "64-bit word XORs in the Macaulay elimination",
            "scope": "decomposition-oracle stage only; not an end-to-end ECDLP cost",
            "metric": "C = (U+1)/p * E[word XORs per target], coverage pooled over samples",
            "samples": [{"seed": d["target_seed"]} for d in docs],
            "census": docs[0]["census"],
            "rungs": out_rows,
        }, indent=2) + "\n")
        print(f"\nWrote {args.output}")

if __name__ == "__main__":
    main()

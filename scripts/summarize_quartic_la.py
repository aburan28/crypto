#!/usr/bin/env python3
"""Section 11.19's table: r = LA / rho at k = 4, measured, against the
outcomes section 11.18 registered.

`research/notes/index-calculus/RESEARCH_RESIDUAL_WALKS.md` §11.16 derived that a
plain k = 4 method tends to r_inf = c_LA / c_rho = 0.34-0.63x rho; §11.18
registered measuring it.  Each run of `examples/gaudry_quartic_la` collects
relations (meet in the middle), filters and solves them as §11.7 does,
verifies the logarithm, and runs rho 16 times on the same group.  This prints
one row per run, the inputs §11.16 derived beside their measured values, the
linear algebra's exponent in n, and r both per curve and with rho's S pooled
over every curve (rho's S does not depend on n, and one curve's 16 runs leave
a standard error near 15 %).

With the frozen C4 measurement as a second argument it also re-derives the
crossover section 11.17 extrapolated at the derived r_inf range, now at the
measured r_inf: p* = 12 C4 / (S_rho c_add (1 - r_inf)), n* = p*^4.  That is
still an extrapolation, on the exponents n^(1/4) (relations) and n^(1/2)
(linear algebra and rho), and is printed as one.

Usage:
    python3 scripts/summarize_quartic_la.py LA.json [C4.json]
"""

from __future__ import annotations

import json
import math
import statistics
import sys

MODN = 16.0  # F_p multiplications per multiplication mod n, as section 11.16 charges


def slope(xs, ys):
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)


def main() -> None:
    rows = json.load(open(sys.argv[1]))
    for r in rows:
        assert r["solved"] and r["correct"], (r["p"], r["seed"])
        assert all(x["correct"] for x in r["rho"]), (r["p"], r["seed"])
    all_s = [x["s"] for r in rows for x in r["rho"]]
    s_pool = statistics.mean(all_s)
    s_se = statistics.stdev(all_s) / math.sqrt(len(all_s))
    print("| p | log2 n | |F| | rate (1/24 = 0.042) | phi (0.73) | weight (5) | LA mod-n muls | LA / N^2 (20) "
          "| rho S, 16 runs | r (this curve's rho) | r (pooled rho) | MITM relation phase / rho (not in r) |")
    print("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    r_pool = []
    for r in rows:
        la_fp = r["la_ops"] * MODN
        rho_unit = math.sqrt(r["n"]) * r["fp_muls_per_add"]
        rp = la_fp / (s_pool * rho_unit)
        r_pool.append(rp)
        print(f"| {r['p']} | {r['bits']:.1f} | {r['base']} | {r['decomposition_rate']:.4f} | {r['phi']:.3f} "
              f"| {r['row_weight']:.2f} | {r['la_ops']:.3e} | {r['wiedemann_constant']:.2f} "
              f"| {r['rho_s_mean']:.3f} ± {r['rho_s_sd'] / math.sqrt(len(r['rho'])):.3f} "
              f"| {r['r']:.3f} | {rp:.3f} | {r['relation_group_ops'] / (s_pool * math.sqrt(r['n'])):.1f} |")
    print()
    print(f"rho S pooled over {len(all_s)} runs on {len(rows)} curves: {s_pool:.3f} ± {s_se:.3f} (s.e.)")
    xs = [math.log(r["n"]) for r in rows]
    la_exp = slope(xs, [math.log(r["la_ops"]) for r in rows])
    r_exp = slope(xs, [math.log(v) for v in r_pool])
    print(f"linear algebra ~ n^{la_exp:.3f} (section 11.16: 1/2); r (pooled) ~ n^{r_exp:+.3f} (flat: 0)")
    by_p: dict = {}
    for r, rp in zip(rows, r_pool):
        by_p.setdefault(r["p"], []).append(rp)
    means = {p: statistics.mean(v) for p, v in by_p.items()}
    print("r (pooled rho) by size: " + ", ".join(f"p={p}: {m:.3f}" for p, m in means.items()))
    r_all = statistics.mean(r_pool)
    r_sd = statistics.stdev(r_pool)
    # The statistic's own uncertainty: seed spread of the LA term, plus rho's s.e.
    r_unc = math.hypot(r_sd / math.sqrt(len(r_pool)), r_all * s_se / s_pool)
    print(f"r_inf = {r_all:.3f} ± {r_unc:.3f}  (cap on the plain k = 4 method's gain over rho: 1/r = {1 / r_all:.2f}x)")
    rate = statistics.mean(r["decomposition_rate"] for r in rows)
    phi = statistics.mean(r["phi"] for r in rows)
    cw = statistics.mean(r["wiedemann_constant"] for r in rows)
    print(f"inputs: rate {rate:.4f} (derived 0.0417), phi {phi:.3f} (borrowed 0.73), "
          f"Wiedemann {cw:.2f} N^2 (derived 20), rho S {s_pool:.3f} (convention 1.3)")
    rising = r_exp > 0 and abs(r_exp) > 2 * (r_sd / r_all) / (max(xs) - min(xs))
    if r_all >= 1 or rising:
        verdict = "WITHDRAWN: a plain k = 4 method does not beat rho"
    elif 0.34 <= r_all <= 0.63:
        verdict = "CONFIRMED: inside the derived 0.34-0.63"
    else:
        verdict = "CORRECTED: below 1 but outside 0.34-0.63"
    print(f"registered outcome: {verdict}")
    if len(sys.argv) > 2:
        c4_rep = json.load(open(sys.argv[2]))
        c4 = statistics.mean(r["solve"]["total_muls"] for r in c4_rep["residuals"])
        c_add = statistics.mean(r["fp_muls_per_add"] for r in rows)
        print()
        print(f"crossover at the measured r_inf (extrapolated; C4 = {c4:.4e}, rho S {s_pool:.3f}, "
              f"c_add {c_add:.1f}):")
        for name, ri in [("r_inf - unc", r_all - r_unc), ("r_inf", r_all), ("r_inf + unc", r_all + r_unc)]:
            pstar = 12.0 * c4 / (s_pool * c_add * (1.0 - ri))
            print(f"  {name} = {ri:.3f}: p* = 2^{math.log2(pstar):.2f}, n* = 2^{4 * math.log2(pstar):.1f}")


if __name__ == "__main__":
    main()

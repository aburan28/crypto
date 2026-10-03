#!/usr/bin/env python3
"""The rho-parity ledger: every measured index-calculus route on the
E(F_{p^k}) subspace harness, its distance to rho at the largest size it was
run at, the exponent that distance moves with, and the numeric condition
parity would need.  Every figure is read from the frozen experiment files;
the two pair-only numbers come from RESEARCH_RESIDUAL_WALKS.md section 11.9
and are marked as such.

Usage:
    python3 scripts/parity_ledger.py [experiments/]

Prints the tables of RESEARCH_RHO_PARITY_PROGRAMME.md.  Extrapolations are
printed with the exponents they rest on and are marked "extrapolated".
"""

from __future__ import annotations

import json
import math
import os
import statistics
import sys
from collections import defaultdict

EXP = sys.argv[1] if len(sys.argv) > 1 else "experiments"


def load(name):
    path = os.path.join(EXP, name)
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return json.load(fh)


def fit(xs, ys):
    """Least-squares exponent of y in x (log-log), with its standard error."""
    lx = [math.log2(x) for x in xs]
    ly = [math.log2(y) for y in ys]
    n = len(lx)
    mx, my = sum(lx) / n, sum(ly) / n
    sxx = sum((a - mx) ** 2 for a in lx)
    sxy = sum((a - mx) * (b - my) for a, b in zip(lx, ly))
    slope = sxy / sxx
    if n > 2:
        resid = [b - my - slope * (a - mx) for a, b in zip(lx, ly)]
        se = math.sqrt(sum(r * r for r in resid) / (n - 2) / sxx)
    else:
        se = float("nan")
    return slope, se


def by_size(rows, key):
    d = defaultdict(list)
    for r in rows:
        d[key(r)].append(r)
    return dict(sorted(d.items()))


def mean(v):
    return sum(v) / len(v)


def doublings_to_parity(ratio, closing):
    """Doublings of n until ratio·2^(closing·k) = 1, closing < 0."""
    if closing >= 0:
        return None
    return math.log2(ratio) / (-closing)


# ── k = 3, plain and the GLV quotient (22_glv_quotient_seeds6.json) ──────

def k3_stores():
    rows = load("22_glv_quotient_seeds6.json")
    if rows is None:
        return None
    out = {}
    groups = by_size([r["report"] for r in rows], lambda r: r["p"])
    for label in ("control", "quotient"):
        sizes, S, rho, rel, la, ratio = [], [], [], [], [], []
        for p, rs in groups.items():
            n = rs[0]["n"]
            c = rs[0]["fp_muls_per_add"]
            sizes.append(n)
            S.append(mean([r[label]["s"] for r in rs]))
            rho.append(mean([r["rho_s_mean"] for r in rs]))
            # relation phase: group ops + oracle, in F_p multiplications
            rel.append(mean([r[label]["group_ops_at_solve"] * c + r[label]["oracle_fp_muls_at_solve"] for r in rs]))
            la.append(mean([r[label]["la_ops"] * 9 for r in rs]))
            ratio.append(S[-1] / rho[-1])
        out[label] = dict(sizes=sizes, S=S, rho=rho, rel=rel, la=la, ratio=ratio)
    return out


# ── k = 3, the linear-algebra variants (21_gaudry_cubic_la.json) ─────────

def k3_la_variants():
    rows = load("21_gaudry_cubic_la.json")
    if rows is None:
        return None
    out = {}
    for r in rows:
        g = r["gaudry"]
        key = "double large primes" if g["max_large_primes"] > 0 else g["la_mode"]
        out.setdefault(key, []).append(r)
    res = {}
    for key, rs in out.items():
        groups = by_size(rs, lambda r: r["gaudry"]["p"])
        sizes, S, rho, rel, la, ratio = [], [], [], [], [], []
        for p, g in groups.items():
            n = g[0]["gaudry"]["n"]
            c = g[0]["gaudry"]["fp_muls_per_add"]
            sizes.append(n)
            S.append(mean([x["gaudry"]["s"] for x in g]))
            rho.append(mean([x["rho"]["s"] for x in g]))
            rel.append(mean([x["gaudry"]["group_ops"] * c + x["gaudry"]["oracle_fp_muls"] for x in g]))
            la.append(mean([x["gaudry"]["la_ops"] * x["gaudry"]["la_fp_equiv"] for x in g]))
            ratio.append(S[-1] / rho[-1])
        res[key] = dict(sizes=sizes, S=S, rho=rho, rel=rel, la=la, ratio=ratio)
    return res


# ── k = 4, full decompositions: C4 and r_inf ─────────────────────────────

def k4_full():
    c4 = load("24_gaudry_quartic_c4.json")
    la = load("25_gaudry_quartic_la.json")
    if c4 is None or la is None:
        return None
    C4 = mean([r["solve"]["total_muls"] for r in c4["residuals"] if r["solve"]["outcome"] == "solved"])
    c_add = c4["fp_muls_per_add"]
    rho_all = [x["s"] for r in la for x in r["rho"]]
    S_rho = mean(rho_all)
    pooled_r = [r["la_ops"] * 16 / (S_rho * math.sqrt(r["n"]) * r["fp_muls_per_add"]) for r in la]
    r_inf = mean(pooled_r)
    r_se = statistics.stdev(pooled_r) / math.sqrt(len(pooled_r))
    return dict(C4=C4, c_add=c_add, S_rho=S_rho, r_inf=r_inf, r_se=r_se,
                sizes=[r["n"] for r in la], rho_runs=len(rho_all))


# ── k = 4, Joux–Vitse three-point decompositions ─────────────────────────

def k4_jv():
    cp = load("26_jv_quartic_cprime.json")
    dl = load("26_jv_quartic_dlp.json")
    out = {}
    if cp is not None:
        groups = by_size(cp, lambda r: r["p"])
        out["cprime"] = [
            dict(p=p, n=g[0]["n"], base=mean([r["base"] for r in g]),
                 c_prime=mean([r["c_prime"]["mean"] for r in g]),
                 weil=mean([r["weil_muls"]["mean"] for r in g]),
                 f4=mean([r["f4_muls"]["mean"] for r in g]),
                 rows=mean([r["max_rows"]["mean"] for r in g]),
                 cols=mean([r["max_cols"]["mean"] for r in g]),
                 c_second_at_y=mean([r["c_second_at_y"]["mean"] for r in g]) if "c_second_at_y" in g[0] else None,
                 c_second_at_x=mean([r["c_second_at_x"]["mean"] for r in g]) if "c_second_at_x" in g[0] else None,
                 ms=mean([r["f4_ms"]["mean"] for r in g]),
                 planted=sum(r["planted_found"] for r in g),
                 constructed=sum(r["constructed_residuals"] for r in g),
                 random=sum(r["random_residuals"] for r in g),
                 mismatches=sum(r["mismatches"] for r in g),
                 undetermined=sum(r["undetermined"] for r in g),
                 mitm=mean([r["mitm3_group_ops"] for r in g]),
                 c_add=g[0]["fp_muls_per_add"])
            for p, g in groups.items()
        ]
    if dl is not None and dl:
        groups = by_size(dl, lambda r: r["p"])
        out["dlp"] = [
            dict(p=p, n=g[0]["n"], seeds=len(g), base=mean([r["base"] for r in g]),
                 residuals=mean([r["residuals"] for r in g]),
                 relations=mean([r["relations"] for r in g]),
                 rate=mean([r["decomposition_rate"] for r in g]),
                 expected_rate=g[0]["expected_rate"],
                 res_floor=mean([r["residual_ratio"] for r in g]),
                 c_prime=mean([r["c_prime"] for r in g]),
                 S=mean([r["s"] for r in g]),
                 walk=mean([r["walk_muls"] / r["total_muls"] for r in g]),
                 oracle=mean([r["oracle_muls"] / r["total_muls"] for r in g]),
                 la=mean([r["la_muls"] / r["total_muls"] for r in g]),
                 rho=mean([r["rho_s_mean"] for r in g]),
                 rho_runs=sum(len(r["rho"]) for r in g),
                 ratio=mean([r["s_over_rho"] for r in g]),
                 rel_ratio=mean([r["relation_over_rho"] for r in g]),
                 r=mean([r["r"] for r in g]),
                 r_last=mean([r["r_last"] for r in g]),
                 la_attempts=mean([r["la_attempts"] for r in g]),
                 formula=mean([3.0 * r["c_prime"] / (r["rho_s_mean"] * r["fp_muls_per_add"]) + r["r"] for r in g]),
                 checked=sum(r["cross_checked"] for r in g),
                 mismatches=sum(r["mismatches"] for r in g),
                 correct=all(r["correct"] for r in g) and all(x["correct"] for r in g for x in r["rho"]),
                 phi=mean([r["phi"] for r in g]),
                 wall_s=mean([r["wall_ms"] / 1e3 for r in g]),
                 c_add=g[0]["fp_muls_per_add"])
            for p, g in groups.items()
        ]
        rho_all = [x["s"] for r in dl for x in r["rho"]]
        out["dlp_rho_pooled"] = (mean(rho_all), statistics.stdev(rho_all) / math.sqrt(len(rho_all)), len(rho_all))
    return out


# ── k = 5, Joux–Vitse four-point decompositions: C″ (27_jv_quintic_*.json) ─

def k5():
    """Per size: C″ for the plain (Weierstrass) and the 2-torsion (Edwards)
    symmetrisation, and the crossover each implies on the derived exponents
    (residuals ∝ n^{2/5}, rho ∝ n^{1/2})."""
    out = {}
    for name, key in (("27_jv_quintic_csecond.json", "weierstrass"), ("27_jv_quintic_edwards_csecond.json", "edwards"), ("28_jv_quintic_edwards4_csecond.json", "edwards4")):
        rows = load(name)
        if not rows:
            continue
        groups = by_size(rows, lambda r: r["p"])
        out[key] = [
            dict(p=p, n=g[0]["n"], seeds=len(g), base=mean([r["base"] for r in g]),
                 c_add=g[0]["fp_muls_per_add"],
                 c_second=mean([r["c_second"]["mean"] for r in g]),
                 c_second_dec=mean([r["c_second_decomposable"]["mean"] for r in g]),
                 weil=mean([r["weil_muls"]["mean"] for r in g]),
                 degree=mean([r["degree_reached"]["mean"] for r in g]),
                 rows=mean([r["max_rows"]["mean"] for r in g]),
                 cols=mean([r["max_cols"]["mean"] for r in g]),
                 c_second_at_y=mean([r["c_second_at_y"]["mean"] for r in g]) if "c_second_at_y" in g[0] else None,
                 c_second_at_x=mean([r["c_second_at_x"]["mean"] for r in g]) if "c_second_at_x" in g[0] else None,
                 ms=mean([r["f4_ms"]["mean"] for r in g]),
                 random=sum(r["random_residuals"] for r in g),
                 decomposable=sum(r["random_decomposable"] for r in g),
                 planted=sum(r["planted_found"] for r in g),
                 constructed=sum(r["constructed_residuals"] for r in g),
                 mismatches=sum(r["mismatches"] for r in g),
                 unverified=sum(r.get("unverified", 0) for r in g),
                 undetermined=sum(r["undetermined"] for r in g),
                 timed_out=sum(r["timed_out"] for r in g),
                 mitm=mean([r["mitm4_group_ops"] for r in g]))
            for p, g in groups.items()
        ]
    census = load("28_jv_quintic_edwards_rate_census.json")
    if census:
        out["census"] = census
    return out


# (M, columns per p, rho reference in the subgroup, cofactor) per symmetrisation
K5_ROUTES = {
    "weierstrass": dict(label="plain", m=16, cpp=0.5, rho_s=1.3, cof=1, rate="1/(24p)"),
    "edwards": dict(label="2-torsion", m=32, cpp=0.25, rho_s=1.3, cof=4, rate="1/(192p)"),
    "edwards4": dict(label="2-torsion, residual saturated by Q₄", m=64, cpp=0.25, rho_s=1.3, cof=4, rate="1/(96p)"),
}


def k5_residuals_per_p2(m, columns_per_p):
    """Residuals per relation set, in units of p².  A class-quadruple of the
    base gives 16 signed sums; M of them (with the free torsion translates)
    are distinct points of the subgroup the residuals live in, so the rate
    of decomposable residuals is M · columns⁴ / (24 · p⁵) and residuals =
    columns / rate = 24 p² / (M · columns_per_p³).  Weierstrass: M = 16 on
    p/2 columns, 12p² (rate 1/(24p)).  Edwards with T free: M = 32 on p/4
    columns, 48p² (rate 1/(192p); the run's `expected_rate` field carried
    1/(24p), corrected by the census 28_jv_quintic_edwards_rate_census.json).
    Edwards with Q₄ free: M = 64, 24p² (rate 1/(96p))."""
    return 24.0 / (m * columns_per_p ** 3)


def k5_crossover(c_second, c_add, m, columns_per_p, rho_s=1.3, cofactor=1):
    """S/rho = residuals · C″ / (rho_s · √n · c_add) with n = p⁵ / cofactor
    the order of the subgroup the logarithm lives in (prime order: cofactor
    1; Edwards: cofactor 4).  With residuals = k5_residuals_per_p2 · p²:
    S/rho = k · C″ · √cofactor / (rho_s · c_add · √p), parity at
    √p* = k · C″ · √cofactor / (rho_s · c_add).  `n*` is the subgroup order
    at `p*`, the convention of the tables' size column."""
    k = k5_residuals_per_p2(m, columns_per_p)
    sqrt_p = k * c_second * math.sqrt(cofactor) / (rho_s * c_add)
    p_star = sqrt_p ** 2
    return p_star, 5 * math.log2(p_star) - math.log2(cofactor)


def k5_c_needed(bits, c_add, m, columns_per_p, rho_s=1.3, cofactor=1):
    """C″ for parity at subgroup order 2^bits: p = (cofactor · 2^bits)^{1/5}."""
    p = (cofactor * 2.0 ** bits) ** 0.2
    k = k5_residuals_per_p2(m, columns_per_p)
    return math.sqrt(p) * rho_s * c_add / (k * math.sqrt(cofactor))


# ── The cover-and-decomposition route on E(F_{p^6}) (30_jv_cover_*.json) ──

def cover(ccov_prefix="30_jv_cover_ccov", dlp_prefix="30_jv_cover_dlp"):
    """Per size: C_cov, the cost of one Nagao test on the genus-3 Jacobian,
    from the cost runs (the oracle run and the sizes run pooled), and the
    end-to-end rows; the registered formula S/rho = 720 C_cov / (rho_S c_add p^2)
    + r, the pooled decomposition rate against 1/720, the fitted exponent of
    S/rho in p, and the crossover the measured constants give."""
    ccov = (load(ccov_prefix + "_oracle.json") or []) + (load(ccov_prefix + ".json") or [])
    dlp = []
    # every end-to-end file of record; the files of defective or superseded runs
    # stay in the repository under their own names and are not read
    names = sorted(
        n for n in os.listdir(EXP)
        if n.startswith(dlp_prefix) and n.endswith(".json")
        and "superseded" not in n and "defective" not in n
    )
    for name in names:
        dlp += load(name) or []
    if not ccov and not dlp:
        return None
    sizes = {}
    for p, g in by_size(ccov, lambda r: r["p"]).items():
        sizes[p] = dict(
            p=p, seeds=len(g), l=g[0]["l"], base=mean([r["base"] for r in g]),
            c_add_e=mean([r["c_add_e"] for r in g]), c_add_j=mean([r["c_add_j"] for r in g]),
            c_cov=mean([r["c_cov"]["mean"] for r in g]), c_cov_dec=mean([r["c_cov_decomposable"]["mean"] for r in g if r["c_cov_decomposable"]["mean"] > 0] or [0]),
            weil=mean([r["weil_muls"]["mean"] for r in g]), f4=mean([r["f4_muls"]["mean"] for r in g]),
            lin=mean([r["lin_muls"]["mean"] for r in g]), post=mean([r["post_muls"]["mean"] for r in g]),
            delta=mean([r["delta"]["mean"] for r in g]), degree=mean([r["degree_reached"]["mean"] for r in g]),
            rows=mean([r["max_rows"]["mean"] for r in g]), cols=mean([r["max_cols"]["mean"] for r in g]),
            f4_ms=mean([r["f4_ms"]["mean"] for r in g]), total_ms=mean([r["total_ms"]["mean"] for r in g]),
            random=sum(r["random_residuals"] for r in g), decomposable=sum(r["random_decomposable"] for r in g),
            constructed=sum(r["constructed_residuals"] for r in g), planted=sum(r["planted_found"] for r in g),
            oracle=sum(r["oracle_checked"] for r in g), mismatches=sum(r["mismatches"] for r in g),
            unverified=sum(r["unverified"] for r in g), incomplete=sum(r["incomplete"] for r in g),
            timed_out=sum(r["timed_out"] for r in g), stopped=sum(r.get("stopped", 0) for r in g),
            stop=g[0].get("stop_staircase"))
    return dict(sizes=sizes, dlp=dlp)


def cover_crossover(c_cov, c_add, rho_s=1.3):
    """S/rho = 360 p C / (rho_S (p^3/2) c_add) = 720 C / (rho_S c_add p^2):
    parity at p*^2 = 720 C / (rho_S c_add); the subgroup there has order
    ~ p*^6 / 4."""
    p2 = 720.0 * c_cov / (rho_s * c_add)
    p_star = math.sqrt(p2)
    return p_star, 6 * math.log2(p_star) - 2


def two_term_minimum(sizes, rel, la, rho_s, c_add):
    """S(n) = A n^(a-1/2) + B n^(b-1/2) fitted through the top size; its
    minimum over n and where.  Returns (min S/rho, log2 n at the minimum,
    a, b) or None when both terms fall."""
    a, _ = fit(sizes, rel)
    b, _ = fit(sizes, la)
    n0 = sizes[-1]
    A = rel[-1] / n0 ** a
    B = la[-1] / n0 ** b
    al, be = a - 0.5, b - 0.5
    if al >= 0 or be <= 0:
        return None
    n_min = (-al * A / (be * B)) ** (1.0 / (be - al))
    s_min = (A * n_min ** al + B * n_min ** be) / c_add
    return s_min / rho_s, math.log2(n_min), a, b


def print_cover(cv, title, label):
    print(f"\n## {title}\n")
    if True:
        if cv["sizes"]:
            print(label + "\n")
            print("| p | ℓ | seeds | columns | c_add E(F_{p⁶}) | c_add Jac | C_cov | Weil | F4 | solver | roots etc. | ideal degree | F4 degree | F4 matrix | F4 ms | test ms | random / decomposable | planted found | oracle checked / mismatches | unverified / incomplete / timed out |")
            print("|---:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|--:|--:|:--|:--|:--|:--|")
            for p, x in cv["sizes"].items():
                print(f"| {p} | 2^{math.log2(x['l']):.1f} | {x['seeds']} | {x['base']:.0f} | {x['c_add_e']:.0f} | {x['c_add_j']:.0f} | {x['c_cov']:.3e} | {x['weil']:,.0f} | {x['f4']:.3e} | {x['lin']:.3e} | {x['post']:,.0f} | {x['delta']:.0f} | {x['degree']:.1f} | {x['rows']:,.0f} × {x['cols']:,.0f} | {x['f4_ms']:.1f} | {x['total_ms']:.1f} | {x['random']} / {x['decomposable']} | {x['planted']}/{x['constructed']} | {x['oracle']} / {x['mismatches']} | {x['unverified']} / {x['incomplete']} / {x['timed_out']} |")
        if cv["dlp"]:
            print("\nThe method end to end (every phase in F_p multiplications, rho on the same group; `*` = rho's S taken from the pooled smaller sizes, its table not fitting in memory).  `exact rate` is 16·C(|F|,6)/ℓ, the number of six-sums over the |F| classes landing in the subgroup, divided by its order: it equals 1/720 only as |F| = p/2 grows (C(|F|,6) = |F|⁶/720 · (1 − 15/|F| + …)); `S/rho` is against the pooled rho S of every walk (rho S scatters ±0.5 over 16 walks at ℓ = 2^32), `own rho` against the row's own walks:\n")
            print("| p | seed | ℓ | columns | residuals | relations | rate | exact rate | residuals / (unknowns / exact rate) | c_add E | C_cov as paid | LA ops (·u²) | S | S / rho (pooled rho S) | own rho S / S/rho | relation + LA | predicted (exact rate) | solved / correct | mismatches |")
            print("|---:|--:|:--|--:|--:|--:|:--|--:|--:|--:|--:|:--|--:|--:|:--|:--|--:|:--|--:|")
            rows = sorted(cv["dlp"], key=lambda r: (r["p"], r["seed"]))
            walks = [w["s"] for r in rows for w in r["rho"]]
            rho_pool = mean(walks) if walks else 1.3
            for r in rows:
                star = "*" if not r["rho"] else ""
                ex = 16.0 * math.comb(r["base"], 6) / r["l"]
                need = (r["base"] + 1) / ex
                sqrt_l = math.sqrt(r["l"])
                s_pool = r["total_muls"] / (rho_pool * sqrt_l * r["c_add_e"]) if "total_muls" in r else r["s_over_rho"]
                pred = need * r["c_cov"] / (rho_pool * sqrt_l * r["c_add_e"]) + r["r"] * r["rho_s_mean"] / rho_pool
                print(f"| {r['p']} | {r['seed']} | 2^{r['bits']:.1f} | {r['base']} | {r['residuals']:,} | {r['relations']} | {r['decomposition_rate']:.5f} | {ex:.5f} | {r['residuals']/need:.2f} | {r['c_add_e']:.0f} | {r['c_cov']:.3e} | {r['la_ops']:,} ({r['wiedemann_constant']:.1f}) | {r['s']:.3e} | {s_pool:.3f} | {r['rho_s_mean']:.3f}{star} / {r['s_over_rho']:.3f} | {r['relation_over_rho']:.3f} + {r['r']:.4f} | {pred:.3f} | {r['solved']} / {r['correct']} | {r['mismatches']} |")
            print(f"\nPooled rho S over {len(walks)} walks: {rho_pool:.3f} ± {statistics.stdev(walks) / math.sqrt(len(walks)):.3f}." if len(walks) > 1 else "")
            tot_res = sum(r["residuals"] for r in rows)
            tot_rel = sum(r["relations"] for r in rows)
            exp = tot_res / 720.0
            exact = sum(r["residuals"] * 16.0 * math.comb(r["base"], 6) / r["l"] for r in rows)
            print(f"\nPooled decomposition rate: {tot_rel} relations in {tot_res:,} residuals = {tot_rel/tot_res*720:.3f}/720 (1/720 predicts {exp:.0f} ± {math.sqrt(exp):.0f}); the exact count 16·C(|F|,6)/ℓ predicts {exact:.1f} ± {math.sqrt(exact):.1f}.")
            g = by_size([r for r in rows if r["solved"]], lambda r: r["p"])
            if len(g) >= 3:
                ps = list(g)
                ys = [mean([r["s_over_rho"] for r in v]) for v in g.values()]
                ys = [mean([r["total_muls"] / (rho_pool * math.sqrt(r["l"]) * r["c_add_e"]) for r in v]) for v in g.values()]
                slope, se = fit(ps, ys)
                print(f"Fitted exponent of S/rho in p over {len(ps)} sizes: {slope:.3f} ± {se:.3f} (registered P4: −2 ± 0.2; derivation: −2).")
        if cv["sizes"]:
            top = list(cv["sizes"].values())[-1]
            cs = mean([x["c_cov"] for x in cv["sizes"].values()])
            ce = mean([x["c_add_e"] for x in cv["sizes"].values()])
            p_star, bits = cover_crossover(cs, ce)
            print(f"\nCrossover implied by the measured constants (mean C_cov = {cs:.3e}, c_add = {ce:.0f}, rho S = 1.3, extrapolated on S/rho ∝ p^{{−2}}): p* = {p_star:,.0f}, subgroup order n* ≈ p*⁶/4 = 2^{bits:.0f}; C_cov/c_add = {cs/ce:,.0f} (registered P5: ≥ 10³).")
            if cv["dlp"]:
                last = sorted(cv["dlp"], key=lambda r: (r["p"], r["seed"]))[-1]
                print(f"At the largest end-to-end size (p = {last['p']}, ℓ = 2^{last['bits']:.1f}): S/rho = {last['s_over_rho']:.3f}.")



def sieve(prefix="32_jv_cover_sieve"):
    """The sieving variant (note §11): every end-to-end row of record, the
    registered P7–P11 checks, and the measured crossover."""
    rows = []
    names = sorted(
        n for n in os.listdir(EXP)
        if n.startswith(prefix) and n.endswith(".json")
        and "superseded" not in n and "defective" not in n
    )
    for name in names:
        rows += load(name) or []
    if not rows:
        return None
    return sorted(rows, key=lambda r: (r["p"], r["seed"]))


def print_sieve(rows, title):
    print(f"\n## {title}\n")
    print("Relations among factor-base points by the sieve of the note's §11 (sums of m points, F(x, s) quadratic in s along a line), the two residuals' descent by the six-point test of §10, Wiedemann mod ℓ, rho on the same group (`*` = pooled reference).  `S` counts F_p multiplications (the ledger's unit); `S⁺` also charges every addition and table lookup of the sieve as a multiplication.  `rate` is relations per line against the registered p/m! (P7); `C_rel` the multiplications per relation, enumeration (the B's factored) + sieve + extraction + verification (P8).\n")
    print("| p | seed | ℓ | columns | m (rule) | available | B's | lines | relations | rate / (p/m!) | C_rel | enum / sieve / extract / verify per relation | mults per base step | descent tests (per success) | C_cov | S | S⁺ | rho S | S / rho | S⁺ / rho | relation + descent + LA | solved / correct | failed verify |")
    print("|---:|--:|:--|--:|:--|--:|--:|--:|--:|--:|--:|:--|--:|:--|--:|--:|--:|:--|--:|--:|:--|:--|--:|")
    for r in rows:
        star = "*" if not r["rho"] else ""
        n = max(r["relations"], 1)
        per = lambda k: r[k] / n
        print(f"| {r['p']} | {r['seed']} | 2^{r['bits']:.1f} | {r['base']} | {r['m_final']} ({r['m_rule']}) | {r['relations_available']:,.0f} | {r['bs']:,} | {r['lines']:,} | {r['relations']} | {r['rate_ratio']:.2f} | {r['c_rel']:.3e} | {per('enum_muls'):.2e} / {per('sieve_muls'):.2e} / {per('extract_muls'):.2e} / {per('verify_muls'):.2e} | {r['sieve_muls'] / max(r['base_steps'], 1):.2f} | {r['descent_residuals']:,} ({r['descent_tests_per_success']:.0f}) | {r['descent_c_cov']:.3e} | {r['s']:.3e} | {r['s_plus']:.3e} | {r['rho_s_mean']:.3f}{star} | {r['s_over_rho']:.4f} | {r['s_plus_over_rho']:.4f} | {r['relation_over_rho']:.4f} + {r['descent_over_rho']:.4f} + {r['la_over_rho']:.5f} | {r['solved']} / {r['correct']} | {r['rels_failed_verify']} |")
    solved = [r for r in rows if r["solved"]]
    g = by_size(solved, lambda r: r["p"])
    if len(g) >= 2:
        ps = list(g)
        ys = [mean([r["s_over_rho"] for r in v]) for v in g.values()]
        yp = [mean([r["s_plus_over_rho"] for r in v]) for v in g.values()]
        below = [p for p, y in zip(ps, ys) if y < 1]
        above = [p for p, y in zip(ps, ys) if y >= 1]
        if below and above and max(above) < min(below):
            p0, p1 = max(above), min(below)
            y0, y1 = ys[ps.index(p0)], ys[ps.index(p1)]
            # log-log interpolation between the two sizes that bracket parity
            t = math.log(y0) / (math.log(y0) - math.log(y1))
            p_star = math.exp(math.log(p0) + t * (math.log(p1) - math.log(p0)))
            y0p, y1p = yp[ps.index(p0)], yp[ps.index(p1)]
            if y0p >= 1 > y1p:
                tp = math.log(y0p) / (math.log(y0p) - math.log(y1p))
                p_plus = math.exp(math.log(p0) + tp * (math.log(p1) - math.log(p0)))
                plus = f"; on S⁺: p* ≈ {p_plus:,.0f} (2^{6 * math.log2(p_plus) - 2:.0f})"
            else:
                plus = ""
            print(f"\nMeasured crossover: S / rho ≥ 1 at p = {p0} ({y0:.3f}) and < 1 at p = {p1} ({y1:.3f}); log-log interpolation between them gives p* ≈ {p_star:,.0f}, subgroup order ≈ 2^{6 * math.log2(p_star) - 2:.0f} (registered P9: p* = 340 ± 60, [280, 420]){plus}.")
        elif below and not above:
            print(f"\nS / rho < 1 at every measured size ({ps[0]}–{ps[-1]}).")
        else:
            print(f"\nNo measured size below parity; the smallest S / rho is {min(ys):.3f} at p = {ps[ys.index(min(ys))]}.")
        m9 = [r for r in solved if r["m_final"] == 9]
        if m9:
            c = mean([r["c_rel"] for r in m9])
            cc = mean([r["descent_c_cov"] for r in m9]) * 720
            print(f"At m = 9 (p = {sorted(set(r['p'] for r in m9))}): C_rel = {c:.3e} per relation against 720 · C_cov = {cc:.3e} per Nagao relation, a ratio of {cc / c:,.0f} (registered P11: 10³ within 3×; P8: C_rel in [5·10⁵, 10⁷]).")
        ratios = [r["rate_ratio"] for r in solved]
        print(f"Rate against p/m! over every solved row: {min(ratios):.2f}–{max(ratios):.2f} (registered P7: [0.5, 2]); relations failing the group check: {sum(r['rels_failed_verify'] for r in rows)}; false hits: {sum(r['false_hits'] for r in rows)}.")
        desc = [r["descent_over_rho"] / r["s_over_rho"] for r in solved if r["s_over_rho"] < 1]
        if desc:
            print(f"Share of the descent in S below parity: {min(desc):.2f}–{max(desc):.2f} (registered P9: more than half).")



def main():
    print("# The parity ledger (printed by scripts/parity_ledger.py)\n")
    k3 = k3_stores()
    la3 = k3_la_variants()
    k4 = k4_full()
    jv = k4_jv()

    print("## A. Where every measured route stands\n")
    print("| route | sizes | S at top | rho S | S / rho | relation phase ∝ n^a | linear algebra ∝ n^b | S / rho ∝ n^c (measured) | asymptote | parity |")
    print("|:--|:--|--:|--:|--:|--:|--:|--:|:--|:--|")

    def row(name, d, rho_ref=None, asym="", parity=""):
        sizes = d["sizes"]
        a, ase = fit(sizes, d["rel"])
        b, bse = fit(sizes, d["la"])
        c, cse = fit(sizes, d["S"])
        rho = rho_ref if rho_ref is not None else d["rho"][-1]
        ratio = d["S"][-1] / rho
        span = f"2^{math.log2(sizes[0]):.1f}–2^{math.log2(sizes[-1]):.1f}"
        print(f"| {name} | {span} | {d['S'][-1]:,.0f} | {rho:.2f} | {ratio:,.0f}× | {a:.2f} ± {ase:.2f} | {b:.2f} ± {bse:.2f} | {c:+.2f} ± {cse:.2f} | {asym} | {parity} |")
        return ratio, c

    if k3:
        for label, name in (("control", "k = 3, plain, ⟨−1⟩ base, O(1) S₄ solve"), ("quotient", "k = 3, ⟨ψ⟩ quotient (j = 0)")):
            d = k3[label]
            tt = two_term_minimum(d["sizes"], d["rel"], d["la"], d["rho"][-1], 63.0)
            asym = "relations n^{1/3}, LA n^{2/3}: S rises after its minimum"
            par = "never" if tt is None else f"never; minimum ≈ {tt[0]:,.0f}× rho near 2^{tt[1]:.0f} (extrapolated on a = {tt[2]:.2f}, b = {tt[3]:.2f})"
            row(name, d, rho_ref=d["rho"][-1], asym=asym, parity=par)
    if la3:
        for key, name in (("wiedemann", "k = 3, plain, Wiedemann (§11.7)"), ("double large primes", "k = 3, double large primes (§11.7)")):
            if key not in la3:
                continue
            d = la3[key]
            rho_ref = 1.3
            if key == "double large primes":
                c, _ = fit(d["sizes"], d["S"])
                ratio = d["S"][-1] / rho_ref
                k_meas = doublings_to_parity(ratio, c)
                k_theory = doublings_to_parity(ratio, -1.0 / 18)
                asym = "n^{4/9} end to end: closes as n^{-1/18}"
                par = f"2^{math.log2(d['sizes'][-1]) + k_meas:.0f} at the measured n^{{{c:+.3f}}}; 2^{math.log2(d['sizes'][-1]) + k_theory:.0f} at n^{{-1/18}} (extrapolated)"
            else:
                tt = two_term_minimum(d["sizes"], d["rel"], d["la"], rho_ref, 63.0)
                asym = "relations n^{1/3}, LA n^{0.68}"
                par = "never" if tt is None else f"never; minimum ≈ {tt[0]:,.0f}× rho near 2^{tt[1]:.0f} (extrapolated)"
            row(name, d, rho_ref=rho_ref, asym=asym, parity=par)
    print("| k = 3, pair-only (k − 1) decompositions (§11.9, from the note) | 2^24.2–2^33.1 | 1,659 | 1.6 | 1,037× | 0.72 | — | +0.22 | residuals n^{2/3}: S rises | never |")
    if k4:
        p_star = 12 * k4["C4"] / (k4["S_rho"] * k4["c_add"] * (1 - k4["r_inf"]))
        n_star = 4 * math.log2(p_star)
        print(f"| k = 4, full decompositions (S₅ solve, §11.17 C₄ + §11.19 r∞) | 2^32.3–2^40.1 | — | {k4['S_rho']:.2f} | relation phase ≫ 10⁸× at 2^32 | 1/4 (derived) | 0.48 ± 0.02 measured (1/2 derived) | → r∞ = {k4['r_inf']:.3f} ± {k4['r_se']:.3f} | tends to r∞ below one | 2^{n_star:.0f} (extrapolated on C₄ = {k4['C4']:.2e}, n^{{1/4}}, n^{{1/2}}) |")
    if jv and "dlp" in jv:
        d = jv["dlp"]
        sizes = [x["n"] for x in d]
        S = [x["S"] for x in d]
        c, cse = fit(sizes, S) if len(sizes) > 1 else (float("nan"), float("nan"))
        rho_p, rho_se, nr = jv["dlp_rho_pooled"]
        ratio = S[-1] / rho_p
        span = f"2^{math.log2(sizes[0]):.1f}–2^{math.log2(sizes[-1]):.1f}"
        print(f"| **k = 4, Joux–Vitse three-point decompositions (this note)** | {span} | {S[-1]:,.0f} | {rho_p:.2f} (pooled, {nr} runs) | **{ratio:,.0f}×** | 1/2 (derived; measured below) | 1/2 (derived) | {c:+.2f} ± {cse:.2f} | **constant**: 3C′/(S_rho·c_add) + r | needs C′ < {(1 - sum(x['r_last'] for x in d) / len(d)) * rho_p * d[-1]['c_add'] / 3:.0f} F_p multiplications; measured C′ = {d[-1]['c_prime']:,.0f} |")

    print("\n## B. The k = 4 Joux–Vitse route, measured\n")
    if jv and "cprime" in jv:
        print("### B.1 C′, the cost of one three-point test (26_jv_quartic_cprime.json)\n")
        print("| p | n | \\|F\\| | random residuals | constructed, planted found | mismatches vs oracle | undetermined | C′ (F_p muls, mean over random) | Weil | F4 | F4 matrix (rows × cols) | F4 ms | oracle (group ops) |")
        print("|---:|:--|--:|--:|:--|--:|--:|--:|--:|--:|:--|--:|--:|")
        for x in jv["cprime"]:
            print(f"| {x['p']} | 2^{math.log2(x['n']):.1f} | {x['base']:.0f} | {x['random']} | {x['planted']}/{x['constructed']} | {x['mismatches']} | {x['undetermined']} | {x['c_prime']:,.0f} | {x['weil']:,.0f} | {x['f4']:,.0f} | {x['rows']:.0f} × {x['cols']:.0f} | {x['ms']:.2f} | {x['mitm']:,.0f} |")
    if jv and "dlp" in jv:
        print("\n### B.2 The method end to end (26_jv_quartic_dlp.json)\n")
        print("| p | n | seeds | \\|F\\| | residuals | relations | rate (1/6p) | residuals / floor | C′ paid | S | walk | oracle | LA | rho S (16 runs per seed) | S / rho | relation phase / rho | r (LA / rho, every attempt) | attempts | r, last attempt | formula 3C′/(S_rho c_add) + r | cross-checked, mismatches | correct |")
        print("|---:|:--|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|")
        for x in jv["dlp"]:
            print(f"| {x['p']} | 2^{math.log2(x['n']):.1f} | {x['seeds']} | {x['base']:.0f} | {x['residuals']:,.0f} | {x['relations']:.0f} | {x['rate']:.5f} ({x['expected_rate']:.5f}) | {x['res_floor']:.2f} | {x['c_prime']:,.0f} | {x['S']:,.0f} | {100*x['walk']:.1f} % | {100*x['oracle']:.1f} % | {100*x['la']:.1f} % | {x['rho']:.2f} | {x['ratio']:,.0f}× | {x['rel_ratio']:,.0f}× | {x['r']:.3f} | {x['la_attempts']:.0f} | {x['r_last']:.3f} | {x['formula']:,.0f}× | {x['checked']}, {x['mismatches']} | {'yes' if x['correct'] else 'NO'} |")
        rho_p, rho_se, nr = jv["dlp_rho_pooled"]
        sizes = [x["n"] for x in jv["dlp"]]
        if len(sizes) > 1:
            cS, cSe = fit(sizes, [x["S"] for x in jv["dlp"]])
            cR, cRe = fit(sizes, [x["residuals"] for x in jv["dlp"]])
            cC, cCe = fit(sizes, [x["c_prime"] for x in jv["dlp"]])
            print(f"\nPooled rho S = {rho_p:.3f} ± {rho_se:.3f} over {nr} runs.  Fitted exponents over {len(sizes)} sizes: S ∝ n^{{{cS:+.3f} ± {cSe:.3f}}} (rho: 0), residuals ∝ n^{{{cR:.3f} ± {cRe:.3f}}} (derived 1/2), C′ ∝ n^{{{cC:+.3f} ± {cCe:.3f}}} (derived 0).")
            print(f"S / rho against the pooled reference: " + ", ".join(f"{x['S']/rho_p:,.0f}× at 2^{math.log2(x['n']):.1f}" for x in jv["dlp"]) + ".")

    k5d = k5()
    if k5d:
        print("\n## D. The k = 5 stage: C″ measured (27_jv_quintic_*.json, 28_jv_quintic_edwards4_csecond.json)\n")
        print("| symmetrisation | p | n | seeds | columns | c_add | C″ (random residuals) | C″ (decomposable) | Weil | F4 degree reached | F4 matrix (rows × cols) | F4 s | random / decomposable | planted found | mismatches | unverified | undetermined / timed out | oracle (group ops) |")
        print("|:--|---:|:--|--:|--:|--:|--:|--:|--:|--:|:--|--:|:--|:--|--:|--:|:--|--:|")
        for key, label in (("weierstrass", "plain, Weierstrass `x`, e(x)"), ("edwards", "2-torsion, Edwards `y`, e(y²) + π"), ("edwards4", "2-torsion + Q₄ saturation (F4 at y_R and at x_R)")):
            for x in k5d.get(key, []):
                print(f"| {label} | {x['p']} | 2^{math.log2(x['n']):.1f} | {x['seeds']} | {x['base']:.0f} | {x['c_add']:.0f} | {x['c_second']:.3e} | {x['c_second_dec']:.3e} | {x['weil']:,.0f} | {x['degree']:.1f} | {x['rows']:,.0f} × {x['cols']:,.0f} | {x['ms']/1e3:.2f} | {x['random']} / {x['decomposable']} | {x['planted']}/{x['constructed']} | {x['mismatches']} | {x['unverified']} | {x['undetermined']} / {x['timed_out']} | {x['mitm']:,.0f} |")
        if "census" in k5d:
            print("\nThe rate of decomposable residuals, counted exactly over every signed four-point sum of the base (28_jv_quintic_edwards_rate_census.json).  `predicted` is 8 and 16 subgroup points per class-quadruple, C(|F|, 4) quadruples; `1/(192p)` and `1/(96p)` are these with |F| = p/4 and C(|F|, 4) = |F|⁴/24, so the measured multiple of 1/(192p) is compared with the same multiple predicted from the actual |F|:\n")
            print("| p | seed | n | columns | class-quadruples | 2-torsion: subgroup points (distinct / predicted) | rate | rate · 192p (measured / predicted) | Q₄-saturated: subgroup points (distinct / predicted) | rate | rate · 96p (measured / predicted) |")
            print("|---:|--:|:--|--:|--:|:--|--:|:--|:--|--:|:--|")
            for c in k5d["census"]:
                pf2 = c['two_torsion_predicted'] / c['n'] * 192 * c['p']; pf4 = c['saturated_predicted'] / c['n'] * 96 * c['p']
                print(f"| {c['p']} | {c['seed']} | 2^{math.log2(c['n']):.1f} | {c['base']} | {c['quadruples']:,} | {c['two_torsion_distinct']:,} / {c['two_torsion_predicted']:,} | {c['two_torsion_rate']:.3e} | {c['two_torsion_rate_x_192p']:.3f} / {pf2:.3f} | {c['saturated_distinct']:,} / {c['saturated_predicted']:,} | {c['saturated_rate']:.3e} | {c['saturated_rate_x_96p']:.3f} / {pf4:.3f} |")
        print("\nCrossovers implied by the measured C″, extrapolated on residuals ∝ n^{2/5} and rho ∝ n^{1/2} (a stage diagnostic: no end-to-end k = 5 S exists).  Residuals = 24p²/(M · (columns/p)³) with M the subgroup points per class-quadruple (16 signs × the free translates that land in the subgroup):\n")
        print("| symmetrisation | C″ (top size) | columns per p | M | rate | residuals | rho reference | parity at p* | n* |")
        print("|:--|--:|--:|--:|--:|--:|--:|--:|--:|")
        for key, r in K5_ROUTES.items():
            xs = k5d.get(key)
            if not xs:
                continue
            x = xs[-1]
            p_star, bits = k5_crossover(x["c_second"], x["c_add"], r["m"], r["cpp"], r["rho_s"], r["cof"])
            print(f"| {r['label']} | {x['c_second']:.2e} | {r['cpp']} | {r['m']} | {r['rate']} | {k5_residuals_per_p2(r['m'], r['cpp']):.0f}p² | {r['rho_s']}·√n, n = p⁵{'/4' if r['cof'] == 4 else ''} | {p_star:.2e} | 2^{bits:.0f} |")
        if "edwards" in k5d and "weierstrass" in k5d:
            e = k5d["edwards"][-1]["c_second"]; w = k5d["weierstrass"][-1]["c_second"]
            ca = k5d["edwards"][-1]["c_add"]
            r = K5_ROUTES["edwards"]
            n128 = k5_c_needed(128, ca, r["m"], r["cpp"], r["rho_s"], r["cof"]); n160 = k5_c_needed(160, ca, r["m"], r["cpp"], r["rho_s"], r["cof"]); n100 = k5_c_needed(100, ca, r["m"], r["cpp"], r["rho_s"], r["cof"])
            print(f"\nThe 2-torsion symmetry cuts C″ by {w/e:,.0f}× at the top size.  For parity at a subgroup order of 2^128 on the 2-torsion route (p = (4·2^128)^{{1/5}} = 2^{(130/5):.1f}): C″ < {n128:.2e}, {e/n128:,.0f}× below the measurement; at 2^160: C″ < {n160:.2e} ({e/n160:,.0f}×); at 2^100: C″ < {n100:.2e} ({e/n100:,.0f}×).")
        if "edwards4" in k5d and "edwards" in k5d:
            e = k5d["edwards"][-1]; e4 = k5d["edwards4"][-1]
            pe = k5_crossover(e["c_second"], e["c_add"], 32, 0.25, 1.3, 4)[1]
            pe4 = k5_crossover(e4["c_second"], e4["c_add"], 64, 0.25, 1.3, 4)[1]
            print(f"\nSaturating the residual by Q₄ doubles the rate (1/(96p)) and costs {e4['c_second']/e['c_second']:.2f}× the test ({e4['c_second_at_y']:.2e} at y_R + {e4['c_second_at_x']:.2e} at x_R against {e['c_second']:.2e}): the crossover moves from 2^{pe:.0f} to 2^{pe4:.0f}, a wash (engineering, {2 * e['c_second'] / e4['c_second']:.2f}× on S).")

    cv = cover()
    if cv:
        print_cover(cv, "E. The cover-and-decomposition route on E(F_{p⁶}) (30_jv_cover_*.json)",
                    "One Nagao test on Jac_H(F_{p²}), genus 3: six quadrics in six unknowns over F_p (Weil restriction of a monic sextic over F_{p²}), F4 and a zero-dimensional solver, every hit verified in the group; oracle = meet in the middle over the three-point sums.")
    cvs = cover("31_jv_cover_stop_ccov", "31_jv_cover_stop_dlp")
    if cvs:
        print_cover(cvs, "E.2 The same route with F4 stopped at the Bézout staircase (31_jv_cover_stop_*.json)",
                    "Identical to E except that F4 stops as soon as the leading monomials found leave 64 standard monomials (the Bézout count of six quadrics in six unknowns), at which point the partial basis is already a Gröbner basis and the remaining steps only certify it; `stopped` counts the tests that stopped there.")
        if cv and cv["sizes"] and cvs["sizes"]:
            common = [p for p in cv["sizes"] if p in cvs["sizes"]]
            if common:
                a = mean([cv["sizes"][p]["c_cov"] for p in common]); b = mean([cvs["sizes"][p]["c_cov"] for p in common])
                fa = mean([cv["sizes"][p]["f4"] for p in common]); fb = mean([cvs["sizes"][p]["f4"] for p in common])
                stopped = sum(cvs["sizes"][p]["stopped"] for p in common); tests = sum(cvs["sizes"][p]["random"] + cvs["sizes"][p]["constructed"] for p in common)
                print(f"\nOver the {len(common)} common sizes: F4 {fa:.3e} → {fb:.3e} ({fa/fb:.2f}×), C_cov {a:.3e} → {b:.3e} ({a/b:.2f}×); {stopped} of {tests} tests stopped at the staircase.  Class: engineering (a constant on the test; the exponent of S/rho is the route's).")
    sv = sieve()
    if sv:
        print_sieve(sv, "F. The sieving variant of the cover route (32_jv_cover_sieve_*.json; note §11)")
    print("\n## C. What parity needs, per route\n")
    if k4:
        for bits in (80, 128, 160):
            p = 2 ** (bits / 4)
            c4_needed = p * k4["S_rho"] * k4["c_add"] * (1 - k4["r_inf"]) / 12
            if c4_needed < k4["C4"]:
                print(f"- k = 4, full decompositions, parity at 2^{bits}: C₄ < {c4_needed:.2e} F_p multiplications ({k4['C4']/c4_needed:,.0f}× below the measured {k4['C4']:.2e}); the S₅ solve's own floor is 7.1·10¹⁰ (§11.16).")
            else:
                print(f"- k = 4, full decompositions, parity at 2^{bits}: C₄ < {c4_needed:.2e} F_p multiplications — met by the measured {k4['C4']:.2e} (the handover 2^{4*math.log2(12*k4['C4']/(k4['S_rho']*k4['c_add']*(1-k4['r_inf']))):.0f} is below 2^{bits}); an extrapolation on n^{{1/4}} and n^{{1/2}}.")
    if jv and "cprime" in jv:
        x = jv["cprime"][-1]
        r = (sum(x["r_last"] for x in jv["dlp"]) / len(jv["dlp"])) if "dlp" in jv else 0.518
        rho_ref = jv["dlp_rho_pooled"][0] if "dlp" in jv else 1.319
        need = (1 - r) * rho_ref * x["c_add"] / 3
        print(f"- k = 4, Joux–Vitse three-point: parity needs C′ < {need:.1f} F_p multiplications at r = {r:.3f} (last attempt); the Weil restriction alone costs {x['weil']:,.0f} and the measured C′ is {x['c_prime']:,.0f} ({x['c_prime']/need:,.0f}× over).  A test {x['c_prime']/1513:.0f}× cheaper would reach §11.16's assumed 36×; no test reaches one.")
    if la3 and "double large primes" in la3:
        d = la3["double large primes"]
        print(f"- k = 3, double large primes: at n^{{-1/18}} the constant must fall by the whole gap for parity at any size that fits: {d['S'][-1]/1.3:,.0f}× at 2^{math.log2(d['sizes'][-1]):.1f}; each 2× on C₃ buys 18 doublings of n.")
    print("- k = 3, plain: no constant reaches parity (the linear algebra's n^{2/3} decides it).")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Summarise results.json (item 06) into summary.json and the figures. Deterministic bootstrap."""
import json
import math
import random

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

R = json.load(open("results.json"))
rng = random.Random(6)
B = 2000
ARMS = ["A1_pkm_mult_chain", "A2_subgroup_flat", "A3_pkm_isogeny_chain", "A4_dickson_chain",
        "A5_dickson_flat", "A6_small_x", "A7_random_x"]
PKM_T4 = {1: 4.00, 2: 4.00, 3: 4.00, 4: 4.00, 5: 4.00, 6: 5.00, 7: 5.00, 8: 5.00}  # D_av^corr


def pois_cdf(k, mu):
    s, t = 0.0, math.exp(-mu)
    for i in range(k + 1):
        s += t
        t *= mu / (i + 1)
    return s


def garwood(k, a=0.05):
    if k > 150:  # exp(-mu) underflows; normal approximation on sqrt scale
        return (max(0.0, math.sqrt(k) - 0.98) ** 2, (math.sqrt(k + 1) + 0.98) ** 2)
    def solve(f, lo, hi):
        for _ in range(200):
            mid = (lo + hi) / 2
            if f(mid):
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2
    lo = 0.0 if k == 0 else solve(lambda mu: 1 - pois_cdf(k - 1, mu) < a / 2, 0, 10 * k + 10)
    hi = solve(lambda mu: pois_cdf(k, mu) > a / 2, 0, 10 * k + 50)
    return lo, hi


def mean(x):
    return sum(x) / len(x)


summary = {"levels": [], "slopes": {}, "pair46": R["pair46"], "stopped": R["stopped"], "cut": R["cut"]}
series = {a: [] for a in ARMS}
for lv in R["levels"]:
    n1 = lv["n1"]
    rho = lv["rho"]
    rho_t = [r["time"] for r in rho]
    row = {"n1": n1, "p": lv["p"], "pD": lv["pD"], "rho_ops_mean": mean([r["ops"] for r in rho]),
           "rho_time_mean": mean(rho_t), "rho_ops_pred": math.sqrt(math.pi * lv["N"] / 2), "arms": {}}
    for arm in ARMS + ["A8_planted"]:
        a = lv["arms"][arm]
        k = sum(a["valid"])
        nt = len(a["valid"])
        pred = a["pred_valid"] * nt
        lo, hi = garwood(k)
        ent = {"F": a["F"], "valid_total": k, "nbsols_total": sum(a["nbsols"]), "targets": nt,
               "pred_total": pred, "ratio": k / pred if pred else None,
               "ratio_ci": [lo / pred, hi / pred] if pred else None}
        g = a.get("gb", [])
        if g:
            T = [s["time"] for s in g]
            ent.update({"T_mean": mean(T), "n_gb": len(g), "timeouts": sum(s["timeout"] for s in g),
                        "D_corr_mean": mean([s["D_corr"] for s in g]),
                        "Dl_corr_mean": mean([s["Dl_corr"] for s in g]),
                        "D_max": max(s["D_max"] for s in g), "Dl_max": max(s["Dl_max"] for s in g),
                        "all_ok": all(s["ok"] for s in g),
                        "sympy_checked": sum(1 for s in g if s.get("sympy_match") is not None),
                        "sympy_match": all(s.get("sympy_match", True) for s in g)})
            # projected IC relation-stage time and ratio to rho, bootstrap over targets
            need = a["F"] + 10
            vals = a["valid"]

            def ic_over_rho(Tb, vb, rb):
                yv = mean(vb)
                if yv == 0:
                    return float("inf")
                return (need / yv) * mean(Tb) / mean(rb)
            pt = ic_over_rho(T, vals, rho_t)
            bs = []
            for _ in range(B):
                idx = [rng.randrange(len(T)) for _ in T]
                idv = [rng.randrange(len(vals)) for _ in vals]
                idr = [rng.randrange(len(rho_t)) for _ in rho_t]
                bs.append(ic_over_rho([T[i] for i in idx], [vals[i] for i in idv], [rho_t[i] for i in idr]))
            bs.sort()
            ent["ic_over_rho"] = pt
            ent["ic_over_rho_ci"] = [bs[int(0.025 * B)], bs[int(0.975 * B) - 1]]
            if not any(s["timeout"] for s in g) and len(g) == lv["arms"][arm].get("gb") and True:
                pass
            series[arm].append((lv["pD"] if "dickson" in arm else lv["p"], T, n1))
        row["arms"][arm] = ent
    summary["levels"].append(row)

# growth slopes log2 T vs log2 p, bootstrap over targets within level
for arm in ARMS:
    pts = [(p, T) for p, T, n1 in series[arm] if len(T) == 12]
    if len(pts) < 3:
        continue

    def slope(xy):
        xs = [x for x, _ in xy]
        ys = [y for _, y in xy]
        mx, my = mean(xs), mean(ys)
        return sum((x - mx) * (y - my) for x, y in xy) / sum((x - mx) ** 2 for x in xs)
    base = [(math.log2(p), math.log2(mean(T))) for p, T in pts]
    s0 = slope(base)
    bs = []
    for _ in range(B):
        xy = []
        for p, T in pts:
            Tb = [T[rng.randrange(len(T))] for _ in T]
            xy.append((math.log2(p), math.log2(mean(Tb))))
        bs.append(slope(xy))
    bs.sort()
    # slope over the top three completed levels too (less overhead-dominated)
    s3 = slope(base[-3:])
    summary["slopes"][arm] = {"slope_all": s0, "ci": [bs[int(0.025 * B)], bs[int(0.975 * B) - 1]],
                              "slope_last3": s3, "levels": [n1 for _, _, n1 in series[arm] if True][:len(pts)],
                              "ic_exponent": s0 + 0.5}
# rho slope
rp = [(math.log2(lv["p"]), math.log2(mean([r["ops"] for r in lv["rho"]]))) for lv in R["levels"]]
mx, my = mean([x for x, _ in rp]), mean([y for _, y in rp])
summary["rho_ops_slope"] = sum((x - mx) * (y - my) for x, y in rp) / sum((x - mx) ** 2 for x, _ in rp)

# published PKM Table 4 / Table 5 Magma times (not our data): slope of log2 T on log2 p
def _sl(ps, ts):
    xs = [math.log2(v) for v in ps]
    ys = [math.log2(v) for v in ts]
    mx, my = mean(xs), mean(ys)
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
summary["pkm_published_slopes"] = {
    "table4_n1_7_11": _sl([17921, 65537, 262657, 1051649, 4206593], [1.14, 9.09, 59.87, 454.79, 4975.07]),
    "table5_last3": _sl([967, 3467, 13619], [0.03, 0.15, 1.16])}

# positive controls
pc = []
for row in summary["levels"]:
    a1 = row["arms"]["A1_pkm_mult_chain"]
    if "D_corr_mean" in a1:
        pc.append({"n1": row["n1"], "pkm_table4": PKM_T4.get(row["n1"]), "D_corr_sugar": a1["D_corr_mean"],
                   "Dl_corr_lcm": a1["Dl_corr_mean"]})
summary["positive_control_degree"] = pc
summary["positive_control_planted"] = [{"n1": r["n1"], "ratio": r["arms"]["A8_planted"]["ratio"]}
                                       for r in summary["levels"]]
json.dump(summary, open("summary.json", "w"), indent=1)

# ---------------- figures ----------------
COL = {"A1_pkm_mult_chain": "#1f5fbf", "A2_subgroup_flat": "#7aa6e8", "A3_pkm_isogeny_chain": "#2a9d6f",
       "A4_dickson_chain": "#c0392b", "A5_dickson_flat": "#e8908a", "A6_small_x": "#8e6bbf",
       "A7_random_x": "#888888"}
LAB = {"A1_pkm_mult_chain": "A1 PKM p-1 subgroup, chained", "A2_subgroup_flat": "A2 same V, flat x^d=1",
       "A3_pkm_isogeny_chain": "A3 PKM isogeny kernel, chained", "A4_dickson_chain": "A4 Dickson image, chained",
       "A5_dickson_flat": "A5 Dickson image, flat", "A6_small_x": "A6 small-x (Semaev), flat",
       "A7_random_x": "A7 random-x (null), flat"}
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
ax = axs[0]
for arm in ARMS:
    xs, ys, lo, hi = [], [], [], []
    for row in summary["levels"]:
        e = row["arms"][arm]
        if "T_mean" in e:
            xs.append(row["n1"])
            ys.append(e["T_mean"])
    if xs:
        ax.plot(xs, ys, "o-", color=COL[arm], label=LAB[arm], lw=1.6, ms=4)
xs = [r["n1"] for r in summary["levels"]]
ax.plot(xs, [r["rho_time_mean"] for r in summary["levels"]], "k--", label="Pollard rho, full DLP (same E)")
ax.set_yscale("log")
ax.set_xlabel("n1  (|V| = 2^n1, p ~ 2^(2 n1))")
ax.set_ylabel("CPU seconds (mean of 12)")
ax.set_title("M2: Groebner time per relation system vs full rho")
ax.legend(fontsize=7)
ax = axs[1]
for arm in ARMS:
    xs, ys = [], []
    for row in summary["levels"]:
        e = row["arms"][arm]
        if "Dl_corr_mean" in e:
            xs.append(row["n1"])
            ys.append(e["Dl_corr_mean"])
    if xs:
        ax.plot(xs, ys, "o-", color=COL[arm], label=LAB[arm], lw=1.6, ms=4)
ks = sorted(k for k in PKM_T4 if k in [r["n1"] for r in summary["levels"]])
ax.plot(ks, [PKM_T4[k] for k in ks], "ks", mfc="none", ms=9, label="PKM Table 4 (Magma, A1)")
ax.set_xlabel("n1")
ax.set_ylabel("max lcm degree of nonzero-reducing pair")
ax.set_yscale("log", base=2)
ax.set_title("M2: solving degree (amendment metric)")
ax.legend(fontsize=7)
ax = axs[2]
for arm in ARMS:
    xs, ys = [], []
    for row in summary["levels"]:
        e = row["arms"][arm]
        if e.get("ic_over_rho") not in (None, float("inf")):
            xs.append(row["n1"])
            ys.append(e["ic_over_rho"])
            ax.plot([row["n1"]] * 2, e["ic_over_rho_ci"], color=COL[arm], lw=0.8, alpha=0.6)
    if xs:
        ax.plot(xs, ys, "o-", color=COL[arm], label=LAB[arm], lw=1.6, ms=4)
ax.axhline(1, color="k", lw=0.8)
ax.set_yscale("log")
ax.set_xlabel("n1")
ax.set_ylabel("projected IC relation stage / rho (time)")
ax.set_title("M3: IC/rho, 95% bootstrap CI")
ax.legend(fontsize=7)
fig.suptitle("Item 06 PKM rational-map factor bases, m=2, 12 targets/level, pure-Python Buchberger, 1 core", fontsize=10)
fig.tight_layout()
fig.savefig("fig_pkm.svg")
fig.savefig("fig_pkm.png", dpi=150)

fig, ax = plt.subplots(figsize=(7, 4.2))
off = {a: (i - 3.5) * 0.08 for i, a in enumerate(ARMS + ["A8_planted"])}
for arm in ARMS + ["A8_planted"]:
    for row in summary["levels"]:
        e = row["arms"][arm]
        if e["ratio"] is None:
            continue
        x = row["n1"] + off[arm]
        ax.plot([x, x], e["ratio_ci"], color=COL.get(arm, "#e67e22"), lw=1)
        ax.plot(x, e["ratio"], "o", color=COL.get(arm, "#e67e22"), ms=4,
                label=(LAB.get(arm, "A8 planted known-log (positive control)") if row["n1"] == 2 else None))
ax.axhline(1, color="k", lw=0.8)
ax.set_yscale("log")
ax.set_xlabel("n1")
ax.set_ylabel("valid relations / heuristic 2|F|^2/N")
ax.set_title("M1: relation yield vs heuristic (95% Garwood CI, pooled 12 targets)")
ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig("fig_yield.svg")
fig.savefig("fig_yield.png", dpi=150)
print(json.dumps({k: summary[k] for k in ("slopes", "pkm_published_slopes", "rho_ops_slope", "positive_control_degree",
                                          "positive_control_planted", "stopped", "cut")}, indent=1))

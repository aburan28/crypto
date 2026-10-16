"""Lays out recorded values only (AGENTS.md: Python never computes a result here).

Reads yield/yield-a*-n*.json (written by examples/cyclotomic_falsifier_yield.rs)
and tables/table-*.md (ecbench table) of research/cyclotomic_falsifier_20261015
and draws SVG graphs, each also rendered to a 150 dpi PNG.
The only arithmetic is the display transform value / null median - 1.
Usage: python3 figures/plot.py <study dir>
"""
import json, os, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = sys.argv[1]
OUT = os.path.join(H, "figures")
plt.rcParams.update({"font.size": 8, "svg.fonttype": "none"})
COL = {"cyclotomic": "#c0392b", "compact_orbit_scan": "#1b6ca8", "onb_weight": "#7d3c98",
       "positive_lambda_progression": "#117a65"}
LABEL = {"cyclotomic": "cyclotomic C_w", "compact_orbit_scan": "compact-orbit-scan",
         "onb_weight": "ONB weight (A2's comparator)", "positive_lambda_progression": "positive control (λ-progression)"}
MARK = {"cyclotomic": "D", "compact_orbit_scan": "o", "onb_weight": "s", "positive_lambda_progression": "^"}


def kind(base):
    return "cyclotomic" if base.startswith("cyclotomic") else base


rungs = [("a1-n23", "K_1, n=23"), ("a0-n23", "K_0, n=23"), ("a1-n11", "K_1, n=11")]
data = {}
for tag, _ in rungs:
    p = os.path.join(H, "yield", f"yield-{tag}.json")
    if os.path.exists(p):
        data[tag] = json.load(open(p))

# Figure 1: every pre-registered cell, null spread vs the arms (relative to the null median).
rows = []
for tag, title in rungs[:2]:  # n = 23 only; n = 11 (r = 991) is tabulated in the README
    if tag not in data:
        continue
    d = data[tag]
    for cell in d["cells"]:
        w, k = cell["w"], cell["k"]
        if k == 1:
            continue  # degenerate by identity (amendment A1): every 1-column base is equal
        for metric in cell["metrics"]:
            nulls = [float(v) for v in d["null"][str(k)][metric]]
            tests = [t for t in d["tests"] if t["w_cell"] == w and t["k"] == k and t["metric"] == metric]
            rows.append((f"{title}\nw={w}, K={k}, {metric}", nulls, tests))

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.6, 0.42 * len(rows) + 1.6), gridspec_kw={"width_ratios": [3, 1.3]})
for i, (label, nulls, tests) in enumerate(rows):
    y = len(rows) - 1 - i
    med = sorted(nulls)[len(nulls) // 2]
    rel = [100 * (v / med - 1) for v in nulls]
    srel = sorted(rel)
    ax.plot([srel[0], srel[-1]], [y, y], color="#bbbbbb", lw=1)
    q = lambda f: srel[min(len(srel) - 1, round(f * (len(srel) - 1)))]
    ax.plot([q(0.025), q(0.975)], [y, y], color="#555", lw=3, alpha=0.5, solid_capstyle="butt")
    for t in tests:
        kd = kind(t["base"])
        v = 100 * (t["value"] / med - 1)
        target = ax2 if kd == "positive_lambda_progression" else ax
        off = {"cyclotomic": 0.22, "compact_orbit_scan": -0.22, "onb_weight": -0.34}.get(kd, 0)
        target.scatter([v], [y + off], marker=MARK[kd], color=COL[kd], s=22 if kd == "cyclotomic" else 14,
                       edgecolor="black" if t["significant"] else "none", lw=0.8, zorder=3)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([r[0] for r in rows][::-1], fontsize=6.5)
ax.axvline(0, color="#999", lw=0.6)
ax.set_xlabel("% difference from the null median (1999 random Frobenius-stable bases: thin = min-max, bar = 2.5-97.5%)")
ax.set_title("Exact yield over all r targets, matched columns", fontsize=8.5)
ax2.set_yticks(range(len(rows)))
ax2.set_yticklabels([])
ax2.axvline(0, color="#999", lw=0.6)
ax2.set_xlabel("% vs null median")
ax2.set_title("positive control", fontsize=8.5)
for kd in ["cyclotomic", "compact_orbit_scan", "onb_weight", "positive_lambda_progression"]:
    ax.scatter([], [], marker=MARK[kd], color=COL[kd], label=LABEL[kd])
ax.legend(fontsize=6.5, loc="lower left", ncol=2)
for a in (ax, ax2):
    a.grid(alpha=0.3, axis="x")
fig.text(0.01, 0.005, "Black edge = significant at p < 0.05/15 (two-sided rank). Source: yield/yield-*.json "
         "(examples/cyclotomic_falsifier_yield.rs, seed 20261015). K=1 cells omitted: identical for every base; n=11 in the README table.",
         fontsize=6, color="#555")
fig.tight_layout(rect=(0, 0.02, 1, 1))
fig.savefig(os.path.join(OUT, "yield-vs-null.svg"))
fig.savefig(os.path.join(OUT, "yield-vs-null.png"), dpi=150)

# Figure 2: ecbench mean S per arm, both curves.
def table(tag):
    rows = {}
    p = os.path.join(H, "tables", f"table-{tag}.md")
    if not os.path.exists(p):
        return rows
    for line in open(p):
        c = [x.strip() for x in line.split("|")]
        if len(c) > 9 and c[1].startswith("ECBS"):
            lo, hi = re.findall(r"[\d.]+", c[8])
            rows[c[2]] = dict(sid=c[1], s=float(c[7]), lo=float(lo), hi=float(hi))
    return rows

fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.4))
for ax, (tag, title) in zip(axes, rungs[:2]):
    t = table(tag)
    names = list(t.keys())
    for i, nm in enumerate(names):
        v = t[nm]
        c = "#c0392b" if "cyc" in nm else ("#1b6ca8" if "cos" in nm else "#555")
        ax.errorbar([i], [v["s"]], yerr=[[v["s"] - v["lo"]], [v["hi"] - v["s"]]], fmt="o", color=c, capsize=3)
        ax.annotate(f"{v['s']:.2f}", (i, v["hi"]), textcoords="offset points", xytext=(0, 3), ha="center", fontsize=6.5)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=30, ha="right", fontsize=7)
    ax.set_yscale("log")
    ax.set_ylabel("mean S (group operations / sqrt r)")
    sid = next(iter(t.values()))["sid"] if t else "-"
    ax.set_title(f"{title}: 16 targets x 2 rounds ({sid})", fontsize=8)
    ax.grid(alpha=0.3, which="both")
fig.text(0.01, 0.005, "ecbench table, operations only (ic.pipeline, rho.signed_frobenius), L0; bars = recorded 95% interval. "
         "Red: cyclotomic; blue: compact-orbit-scan at equal columns.", fontsize=6, color="#555")
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(os.path.join(OUT, "cost-s-by-arm.svg"))
fig.savefig(os.path.join(OUT, "cost-s-by-arm.png"), dpi=150)
print("ok")

"""Lays out recorded values only (AGENTS.md: Python never computes a result here).

Reads ecbench tables, comparison files and the `ecbench trend` output of
research/ecbench_frobenius_orbit_vs_rho_20261014 and draws three SVG graphs,
each also rendered to a 150 dpi PNG. Usage: python3 figures/plot.py <study dir>
"""
import json, re, sys, glob, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = sys.argv[1]
OUT = os.path.join(H, "figures")
CURVES = ["a0-n13", "a1-n17", "a0-n19", "a0-n23", "a1-n29", "a0-n31", "a0-n37", "a0-n41"]
plt.rcParams.update({"font.size": 9, "svg.fonttype": "none"})

def table(s):
    rows = {}
    for line in open(os.path.join(H, "tables", f"table-{s}.md")):
        c = [x.strip() for x in line.split("|")]
        if len(c) > 9 and c[1].startswith("ECBS"):
            lo, hi = re.findall(r"[\d.]+", c[8])
            rows[c[2]] = dict(sid=c[1], log2r=float(c[5]), s=float(c[7]), lo=float(lo), hi=float(hi))
    return rows

def arm(rows, pat):
    for k, v in rows.items():
        if re.match(pat, k):
            return k, v
    return None, None

T = {s: table(s) for s in CURVES}
SERIES = [("rho-frob", r"^rho-frob$", "#1b6ca8", "o"), ("rho-strong", r"^rho-strong$", "#7a7a7a", "s"),
          ("IC orbit, stride-renew", r"^ic-orbit-.*-renew$", "#c0392b", "D"),
          ("IC orbit, walk targets", r"^ic-orbit-.*-walk$", "#e59866", "v"),
          ("IC subspace negfold, stride-renew", r"^ic-negfold", "#6c3483", "^")]

# Figure 1: mean S of every arm against log2 r.
fig, ax = plt.subplots(figsize=(7.2, 4.2))
for label, pat, col, mk in SERIES:
    xs, ys, lo, hi = [], [], [], []
    for s in CURVES:
        _, v = arm(T[s], pat)
        if v:
            xs.append(v["log2r"]); ys.append(v["s"]); lo.append(v["s"] - v["lo"]); hi.append(v["hi"] - v["s"])
    ax.errorbar(xs, ys, yerr=[lo, hi], fmt=mk, color=col, label=label, capsize=2, ms=5, lw=1)
for s in CURVES:
    _, v = arm(T[s], r"^rho-frob$")
    ax.annotate(s.replace("a", "a=").replace("-n", ", n="), (v["log2r"], 0.07), fontsize=6.5, rotation=90, ha="center", va="bottom", color="#555")
ax.set_yscale("log"); ax.set_ylim(0.05, 120)
ax.set_xlabel("log2 r (prime subgroup order)"); ax.set_ylabel("mean S = GAE / sqrt(r), fully priced")
ax.set_title("Every arm, every curve: 16 targets x 3 rounds; bars = recorded 95% interval", fontsize=9)
ax.legend(fontsize=7, loc="upper right"); ax.grid(alpha=0.3, which="both")
fig.text(0.01, 0.005, "Sessions: research/ecbench_frobenius_orbit_vs_rho_20261014/sessions/*; ecbench table", fontsize=6, color="#555")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "s-by-arm.svg")); fig.savefig(os.path.join(OUT, "s-by-arm.png"), dpi=150)

# Figure 2: paired ratios and the trend fits.
def comps(pattern):
    pts = []
    for f in sorted(glob.glob(os.path.join(H, "sessions", "*", "comparisons", pattern))):
        d = json.load(open(f))
        c = d["curves"][0]
        pts.append((c["slug"], c["log2_r"], c["ratio_b_over_a"], c["ci95"]))
    return pts
trend = json.load(open(os.path.join(H, "trend", "trend-orbit-vs-rho-frob.json")))
fig, ax = plt.subplots(figsize=(7.2, 4.2))
for pat, col, mk, label in [("rho-frob__ic-orbit-*-renew.json", "#c0392b", "D", "IC orbit ÷ rho-frob"),
                            ("rho-strong__ic-orbit-*-renew.json", "#7a7a7a", "s", "IC orbit ÷ rho-strong")]:
    p = comps(pat)
    ax.errorbar([x[1] for x in p], [x[2] for x in p], yerr=[[x[2] - x[3][0] for x in p], [x[3][1] - x[2] for x in p]],
                fmt=mk, color=col, capsize=2, ms=5, lw=1, label=label)
    if "frob" in pat:
        for x in p:
            ax.annotate("n=" + re.search(r"f2m(\d+)", x[0]).group(1), (x[1], x[3][1]), fontsize=6.5, xytext=(3, 2), textcoords="offset points")
pf = {f["name"]: f for f in trend["fits"]}
xs = [pt["log2_r"] for pt in trend["points"]]
for name, style, lab in [("free_r_n", "x", "fit ln ratio = c + β ln r + γ ln n (per curve)"), ("model", "+", "model β = 1/6, γ = -1/2, c fitted")]:
    ax.scatter(xs, pf[name]["fitted"], marker=style, color="#1b6ca8" if name == "model" else "#117a65", s=30, label=lab, zorder=5)
b = pf["free_r"]["coefficients"][1]
ax.axhline(1, color="k", lw=0.8)
ax.set_yscale("log"); ax.set_xlabel("log2 r"); ax.set_ylabel("paired ratio of S (B ÷ A), cluster-bootstrap 95%")
ax.set_title(f"IC ÷ rho grows with r: free_r slope β = {b[1]:.3f} [{b[2][0]:.3f}, {b[2][1]:.3f}] per ln r", fontsize=9)
ax.legend(fontsize=7, loc="upper left"); ax.grid(alpha=0.3, which="both")
fig.text(0.01, 0.005, "ecbench compare --save; ecbench trend (trend/trend-orbit-vs-rho-frob.json)", fontsize=6, color="#555")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "ratio-trend.svg")); fig.savefig(os.path.join(OUT, "ratio-trend.png"), dpi=150)

# Figure 3: phase split of the orbit arm (derived means of recorded phase GAE).
rows = [json.loads(l) for l in open(os.path.join(H, "tables", "derived-phase-means.jsonl"))]
fig, ax = plt.subplots(figsize=(9.0, 3.8))
labels = [r["s"] for r in rows]
import math
parts = [("factor base (construction)", "ic_fb_phase", "#a9cce3"), ("oracle set-up (pair table)", "ic_setup", "#5499c7"),
         ("relations (trials)", "ic_rel", "#c0392b")]
bottom = [0] * len(rows)
for lab, key, col in parts:
    vals = [r[key] / r["ic_total"] for r in rows]
    ax.bar(labels, vals, bottom=bottom, color=col, label=lab)
    bottom = [a + b for a, b in zip(bottom, vals)]
rest = [1 - b for b in bottom]
ax.bar(labels, rest, bottom=bottom, color="#d5d8dc", label="linear algebra + verify")
ax.plot(labels, [r["ic_fb_pre_uncharged"] / r["ic_total"] for r in rows], "k.-", lw=0.8, label="share newly priced (was *_uncharged)")
ax.set_ylabel("share of IC orbit total GAE (mean over 48 runs)"); ax.set_ylim(0, 1)
ax.set_title("Where the fully priced orbit IC spends its cost", fontsize=9); ax.legend(fontsize=6.5, loc="upper left", bbox_to_anchor=(1.01, 1.0))
fig.text(0.01, 0.005, "tables/derived-phase-means.jsonl (jq means of recorded phase GAE)", fontsize=6, color="#555")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "orbit-phase-split.svg")); fig.savefig(os.path.join(OUT, "orbit-phase-split.png"), dpi=150)
print("ok")

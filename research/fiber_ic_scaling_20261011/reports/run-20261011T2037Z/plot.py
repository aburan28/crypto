"""Plot recorded mean S per arm against log2 #E for ecbench session scaling_h4 (recorded values only)."""
import json, os, collections
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt
here = os.path.dirname(os.path.abspath(__file__))
recs = [json.loads(l) for l in open(os.path.join(here, "../../sessions/scaling_h4/records.jsonl"))]
by = collections.defaultdict(list)
order = {}
for r in recs:
    if r["warmup"] or r["outcome"]["status"] != "verified":
        continue
    c = r["workload"]["curve"]
    by[(r["arm"], c["slug"])].append(r["cost"]["s"])
    order[c["slug"]] = c["group_order"]
import math
arms = sorted({a for a, _ in by}, key=lambda a: (a != "rho-neg", a.startswith("blind"), int(''.join(ch for ch in a if ch.isdigit()) or 0)))
fig, ax = plt.subplots(figsize=(7.5, 4.4))
cmap = plt.get_cmap("viridis")
ic = [a for a in arms if a != "rho-neg"]
for a in arms:
    pts = sorted((math.log2(order[s]), sum(v) / len(v), min(v), max(v)) for (aa, s), v in by.items() if aa == a)
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    lo = [p[1] - p[2] for p in pts]; hi = [p[3] - p[1] for p in pts]
    if a == "rho-neg":
        kw = dict(color="black", marker="o", lw=2, label="rho-neg (reference)")
    elif a.startswith("blind"):
        kw = dict(color="#c0632b", marker="s", ls="--", label=a)
    else:
        kw = dict(color=cmap(ic.index(a) / max(1, len(ic) - 1)), marker="^", label=a)
    ax.errorbar(xs, ys, yerr=[lo, hi], capsize=2, ms=4, **kw)
ax.set_yscale("log"); ax.set_xlabel("log2 #E"); ax.set_ylabel("S (cost ÷ floor), mean over verified runs; bars min–max")
ax.set_title("scaling_h4 (ECBS1h4a5009a8f716), 2 curves per size, 3 targets × 2 rounds per arm", fontsize=9)
ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=7, ncol=2)
fig.tight_layout(); fig.savefig(os.path.join(here, "s_by_arm.svg")); print("ok")

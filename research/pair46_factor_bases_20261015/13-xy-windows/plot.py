"""Layout only: plots values recorded in results.json (no result is computed here)."""
import json, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = sys.argv[1] if len(sys.argv) > 1 else "."
r = json.load(open(f"{d}/results.json"))
arms = r["arms"]
g = r["toy"]["groups"]

# Figure 1: matched-size yield ratios vs x-only window X
fig, axs = plt.subplots(1, 2, figsize=(11, 4.3), sharey=False)
show = [a for a in range(len(arms)) if arms[a] != "X"]
for ax, key, title in [(axs[0], "ratio_D", "decompositions per target D / D(X)"), (axs[1], "ratio_M2", "targets decomposed M2 / M2(X)")]:
    for gi, (gname, mk, off) in enumerate([("all", "o", 0.0), ("F", "s", -0.2), ("N2", "^", 0.2)]):
        for j, a in enumerate(show):
            v = g[gname][key][a]
            col = "#c0392b" if arms[a].startswith("XY") else ("#2e86c1" if arms[a] == "Xbits4" else "#7f8c8d")
            ax.errorbar(j + off, v["point"], yerr=[[v["point"] - v["ci95"][0]], [v["ci95"][1] - v["point"]]], fmt=mk, color=col, ms=4, capsize=2,
                        label=(f"{gname} curves" if j == 0 else None), mec="k", mew=0.3)
    ax.axhspan(0.85, 1.15, color="#2ecc71", alpha=0.12, label="null band [0.85,1.15]")
    ax.axhline(1, color="k", lw=0.6)
    ax.axhline(0.25, color="#999", lw=0.6, ls=":")
    ax.set_xticks(range(len(show))); ax.set_xticklabels([arms[a] for a in show], rotation=40, ha="right", fontsize=8)
    ax.set_title(title, fontsize=10); ax.set_ylim(0, 1.35)
axs[0].legend(fontsize=7, loc="lower left")
fig.suptitle(f"13-xy-windows: matched size K=ceil(sqrt(N/2)), {g['all']['n_curves']} toy curves x {r['toy']['T']} targets (95% curve bootstrap)", fontsize=10)
fig.tight_layout()
for ext in ("svg", "png"): fig.savefig(f"{d}/fig_yield.{ext}", dpi=150)

# Figure 2: pair-46 dispersion and corner density
fig, axs = plt.subplots(1, 2, figsize=(11, 4.3))
cs = r["pair46"]["curves"]
cols = {"crater": "#c0392b", "floor": "#e67e22"}
for ci, c in enumerate(cs):
    col = cols.get(c["name"], "#7f8c8d")
    for si, dd in enumerate(c["dispersion"]):
        x = si + (ci - 2.5) * 0.1
        axs[0].errorbar(x, dd["index"], yerr=[[dd["index"] - dd["ci95"][0]], [dd["ci95"][1] - dd["index"]]], fmt="o", color=col, ms=4, capsize=2, label=c["name"] if si == 0 else None)
        if "planted" in dd:
            q = dd["planted"]
            axs[0].errorbar(si + 0.35, q["index"], yerr=[[q["index"] - q["ci95"][0]], [q["ci95"][1] - q["index"]]], fmt="D", color="#27ae60", ms=4, capsize=2, label="planted (crater +8 in 10%)" if si == 0 else None)
axs[0].axhline(1, color="k", lw=0.6); axs[0].axhspan(0.97, 1.03, color="#2ecc71", alpha=0.12)
axs[0].set_xticks(range(len(cs[0]["dispersion"]))); axs[0].set_xticklabels([f"B1=2^{dd['B1'].bit_length()-1}" for dd in cs[0]["dispersion"]])
axs[0].set_ylabel("variance / mean of box counts (lambda~16)"); axs[0].set_title("pair 46: box-count dispersion, x in [0,2^23)", fontsize=10)
p = r["pair46"]["p"]
for si, dd in enumerate(cs[0]["dispersion"]):
    nl = 1 - 0.5 * dd["H"] / ((p + 1) // 2)
    axs[0].hlines(nl, si - 0.35, si + 0.25, color="#2e86c1", lw=1.2, ls="--", label="exact binomial null 1-q (post hoc)" if si == 0 else None)
axs[0].legend(fontsize=7)
for ci, c in enumerate(cs):
    col = cols.get(c["name"], "#7f8c8d")
    ks = [e["k"] for e in c["corner"]]
    axs[1].errorbar([k + (ci - 2.5) * 0.08 for k in ks], [e["obs_over_exp"] for e in c["corner"]],
                    yerr=[[max(0, e["obs_over_exp"] - e["ci95"][0]) for e in c["corner"]], [max(0, e["ci95"][1] - e["obs_over_exp"]) for e in c["corner"]]],
                    fmt=".", color=col, ms=4, capsize=1.5, label=c["name"])
axs[1].axhline(1, color="k", lw=0.6)
axs[1].set_xlabel("k  (|y| < p/2^(k+1))"); axs[1].set_ylabel("observed / expected (2^-k)")
axs[1].set_title("pair 46: small-x, small-|y| corner density", fontsize=10); axs[1].legend(fontsize=7)
fig.tight_layout()
for ext in ("svg", "png"): fig.savefig(f"{d}/fig_pair46.{ext}", dpi=150)
print("plots written")

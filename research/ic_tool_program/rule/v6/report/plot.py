"""Rule v6's graph: the cold ratio, index calculus over rho on the same point,
at rule v3 and rule v6, by size, with 95% intervals; and the online speedup.
Lays out the two analysis.json files; computes nothing new."""
import json, sys
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

v6 = json.load(open(sys.argv[1]))["sizes"]
v3 = json.load(open(sys.argv[2]))["sizes"]
out = sys.argv[3]
fig, (a, b) = plt.subplots(1, 2, figsize=(10.0, 4.2))
for ax, key, title, ylab in (
    (a, "cold_ratio_ic_over_rho", "cold: IC over rho (above 1, rho faster)", "cold time, IC / rho"),
    (b, "online_speedup", "online: rho over IC (above 1, IC faster)", "online speedup"),
):
    for sizes, label, colour, dx in ((v3, "rule v3", "#9a9a9a", -0.1), (v6, "rule v6", "#1f5fa8", 0.1)):
        xs = [s["log2_r"] + dx for s in sizes]
        val = [s[key].get("value", s[key].get("mean_ratio")) for s in sizes]
        lo = [v - s[key]["ci95"][0] for v, s in zip(val, sizes)]
        hi = [s[key]["ci95"][1] - v for v, s in zip(val, sizes)]
        ax.errorbar(xs, val, yerr=[lo, hi], fmt="o", color=colour, capsize=3, ms=4, label=label)
    ax.axhline(1.0, color="black", lw=0.8)
    ax.set_yscale("log")
    ax.set_xlabel("log₂ r (§23's six sizes)")
    ax.set_ylabel(ylab)
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8)
fig.text(0.01, 0.01, "research/ic_tool_program/rule/v6/analysis.json and rule/v3/analysis.json; 408 checked rows each; 95% bootstrap intervals over 64 targets", fontsize=7)
fig.tight_layout(rect=(0, 0.04, 1, 1))
fig.savefig(out + ".svg")
fig.savefig(out + ".png", dpi=150)

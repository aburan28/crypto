"""R10's graph: cold-time ratio v3/main per size with 95% intervals, suite and
holdouts, beside the round's own A/A band. Lays out analysis.json; computes
nothing new."""
import json, sys
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

a = json.load(open(sys.argv[1]))
out = sys.argv[2]
suite = a["suite"]
hold = {h["slug"]: h for h in a["holdouts"]}
fig, ax = plt.subplots(figsize=(8.0, 4.4))
xs = range(len(suite))
for i, s in enumerate(suite):
    aa = s["aa_cold"]
    ax.add_patch(plt.Rectangle((i - 0.35, aa["lo"]), 0.7, aa["hi"] - aa["lo"], color="#d9d9d9", zorder=0, lw=0))
    c = s["cold"]
    ax.errorbar(i - 0.1, c["geomean"], yerr=[[c["geomean"] - c["lo"]], [c["hi"] - c["geomean"]]], fmt="o", color="#1f5fa8", capsize=3, ms=4, label="suite rows" if i == 0 else None)
    if s["slug"] in hold:
        h = hold[s["slug"]]["cold"]
        ax.errorbar(i + 0.1, h["geomean"], yerr=[[h["geomean"] - h["lo"]], [h["hi"] - h["geomean"]]], fmt="s", color="#c2571a", capsize=3, ms=4, label="fresh holdouts" if s["slug"] == suite[-3]["slug"] else None)
ax.axhline(1.0, color="black", lw=0.8)
ax.set_xticks(list(xs))
ax.set_xticklabels([f"{s['log2_r']:.1f}" for s in suite])
ax.set_xlabel("log₂ r of the suite curve (ICV1 slugs in the report's table)")
ax.set_ylabel("cold time, v3 (rebuilt) / main 92246530")
ax.set_title("R10: main's head against v3, paired, 95% t intervals; grey: the round's own A/A band")
ax.legend(loc="upper left", fontsize=8)
fig.text(0.01, 0.01, "research/ic_tool_program/rounds/R10-main-head/analysis.json; 943 isolated processes, 3 contended (excluded)", fontsize=7)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(out + ".svg")
fig.savefig(out + ".png", dpi=150)

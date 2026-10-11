"""Plot the recorded E4 values (e4/records.jsonl). Lays out recorded numbers only."""
import json, sys, os
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
recs = [json.loads(l) for l in open(os.path.join(here, "../../e4/records.jsonl"))]
sys_ = [r for r in recs if r["kind"] == "system"]

# Figure 1: row-reduction multiplications against the base-polynomial degree.
fig, ax = plt.subplots(figsize=(7, 4.2))
style = {1: ("o", "#3b6ea8", "direct on E (|H| = 1)"), 2: ("s", "#c0632b", "quotient E/H, |H| = 2"),
         4: ("^", "#3f8f4f", "quotient E/H, |H| = 4")}
for h, (m, c, lab) in style.items():
    for kind, fill in (("random", c), ("planted", "none")):
        pts = [(r["base_degree"], r["field_ops"]) for r in sys_ if r["subgroup_order"] == h and r["target_kind"] == kind]
        if pts:
            xs, ys = zip(*pts)
            ax.scatter(xs, ys, marker=m, facecolors=fill, edgecolors=c, s=46,
                       label=f"{lab}, {kind} target", zorder=3)
ax.set_xscale("log", base=2); ax.set_yscale("log")
ax.set_xticks([8, 16, 32, 64]); ax.set_xticklabels(["8", "16", "32", "64"])
ax.set_xlabel("degree of the factor-base polynomial in each unknown (B / |H|)")
ax.set_ylabel("F_p multiplications in F4 row reduction")
ax.set_title("E4: S3 systems on E and on E/H, 1128 systems, 72 curves", fontsize=10)
ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=7, loc="upper left")
fig.tight_layout(); fig.savefig(os.path.join(here, "e4_ops_vs_degree.svg"))

# Figure 2: crater E[2] against floor and floor-model C4 at B = 64, random targets.
pairs = sorted({r["pair"] for r in sys_ if r["role"] == "crater" and r["two_rank"] == 2})
def mean_ops(pair, role, sub):
    v = [r["field_ops"] for r in sys_ if r["pair"] == pair and r["role"] == role and r["subgroup"] == sub
         and r["base_size"] == 64 and r["target_kind"] == "random"]
    return sum(v) / len(v) if v else None
fig, ax = plt.subplots(figsize=(7, 4.2))
w = 0.27
for k, (role, sub, c) in enumerate((("crater", "E[2]", "#3b6ea8"), ("floor", "C4", "#c0632b"), ("floor-iso", "C4", "#3f8f4f"))):
    ys = [mean_ops(p, role, sub) or 0 for p in pairs]
    ax.bar([i + (k - 1) * w for i in range(len(pairs))], ys, w, color=c, label=f"{role}, H = {sub}")
ax.set_xticks(range(len(pairs))); ax.set_xticklabels(pairs, rotation=40, ha="right", fontsize=7)
ax.set_ylabel("F_p row-reduction mults\n(mean, 2 random targets)")
ax.set_title("E4 at B = 64, |H| = 4: crater with rational E[2] vs floor and floor model with C4", fontsize=9)
ax.set_ylim(0, 45000); ax.legend(fontsize=7, ncol=3, loc="upper center"); ax.grid(True, axis="y", alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(here, "e4_gap_b64.svg"))
print("ok")

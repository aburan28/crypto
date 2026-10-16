"""Layout only: plots values recorded in data/summary_p*.json.

python research/smith_iso_granularity_20261015/figures/plot.py
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
PRIOR = HERE.parent.parent / "smith_prop24_native_20261014" / "data"
plt.rcParams.update({"font.size": 8, "svg.fonttype": "none"})

ps = sorted(int(f.stem.split("_p")[1]) for f in DATA.glob("summary_p*.json"))
S = {p: json.loads((DATA / f"summary_p{p}.json").read_text()) for p in ps}
P = {p: json.loads((PRIOR / f"census_p{p}_summary.json").read_text()) for p in ps}

fig, axes = plt.subplots(1, 2, figsize=(7.4, 2.9), gridspec_kw={"width_ratios": [1, 1.25]})
ax = axes[0]
for j, p in enumerate(ps):
    s, x = S[p], j
    ax.hlines([s["D"], s["U"]], x - 0.32, x + 0.32, colors=["#999", "#999"], lw=1, ls=":")
    ax.text(x + 0.34, s["U"], "U", va="center", fontsize=7, color="#666")
    ax.text(x + 0.34, s["D"], "D", va="center", fontsize=7, color="#666")
    pm = P[p]["null"]
    ax.errorbar(x - 0.15, pm["class_fraction_median"], fmt="s", color="#c77",
                label="model-level null (2026-10-14)" if j == 0 else None)
    n = s["null_primary"]
    ax.errorbar(x, n["median"], yerr=[[n["median"] - n["lo"]], [n["hi"] - n["median"]]],
                fmt="o", color="#3a6fb0", capsize=3,
                label="iso-level null, 95% (999 perms)" if j == 0 else None)
    ax.plot(x + 0.15, s["C"], "D", color="#e08a1e", label="observed C" if j == 0 else None)
ax.set_xticks(range(len(ps)), [f"p = {p}" for p in ps])
ax.set_xlim(-0.6, len(ps) - 0.3)
ax.set_ylabel("share of models in an L-class\nwith a full_success member")
ax.set_title("(a) class-level share C against the nulls", fontsize=8)
ax.legend(fontsize=6.5, loc="lower center", bbox_to_anchor=(0.5, 0.12))

ax = axes[1]
w = 0.38
for j, p in enumerate(ps):
    bk = S[p]["by_k"]
    labs = sorted(bk, key=lambda k: int(k.strip().rstrip("+")))
    xs = [i + (j - 0.5) * w for i in range(len(labs))]
    obs = [bk[k]["no_full"] / bk[k]["classes"] for k in labs]
    exp = [bk[k]["expected_no_full"] / bk[k]["classes"] for k in labs]
    ax.bar(xs, obs, w * 0.95, color=["#e08a1e", "#b5651d"][j], label=f"observed, p = {p}")
    ax.plot(xs, exp, "_", ms=9, mew=2, color="#1d3f6e",
            label="independent iso classes" if j == 0 else None)
ax.set_xticks(range(len(labs)), [k.strip() for k in labs])
ax.set_xlabel("k = tractable isomorphism classes in the L-class")
ax.set_ylabel("share of tractable L-classes\nwith no full_success member")
ax.set_title("(b) empty classes by k: observed vs independence", fontsize=8)
ax.legend(fontsize=6.5)
fig.tight_layout()
fig.savefig(HERE / "iso_granularity.svg")
fig.savefig(HERE / "iso_granularity.png", dpi=150)

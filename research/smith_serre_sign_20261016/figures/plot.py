"""Layout only: plots values recorded in data/serre_p*.json.

python research/smith_serre_sign_20261016/figures/plot.py
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
plt.rcParams.update({"font.size": 8, "svg.fonttype": "none"})
ps = sorted(int(f.stem.split("_p")[1]) for f in DATA.glob("serre_p*.json"))
S = {p: json.loads((DATA / f"serre_p{p}.json").read_text()) for p in ps}
col = {7: "#e08a1e", 11: "#3a6fb0"}

fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.0), gridspec_kw={"width_ratios": [1.2, 1]})
ax = axes[0]
ax.plot([0, 1], [0, 1], color="#999", lw=1, ls=":", label="N1 / (2p+2)")
for p in ps:
    b = S[p]["h1"]["bins"]
    x = [d["N1"] / (2 * p + 2) for d in b]
    y = [d["r"] for d in b]
    ax.plot(x, y, "-o", ms=3.5, color=col[p], label=f"p = {p}" + (" (exploratory)" if p == 7 else ""))
    for d in b:
        if not d["in_test"]:
            ax.plot(d["N1"] / (2 * p + 2), d["r"], "o", ms=5, mfc="white", mec=col[p])
ax.set_xlabel("N1 / (2p + 2),  N1 = #H(F_p)")
ax.set_ylabel("Prop. 4 pass rate among non-degenerate\nProp. 2 passes, r(N1)")
ax.set_title("(a) H1: the Prop. 4 sign follows #H(F_p)", fontsize=8)
ax.legend(fontsize=6.5, loc="upper left")

ax = axes[1]
lab = ["stratum\n(2026-10-15)", "random label\n+ Prop. 2 (control)", "N1 + Prop. 2\n(H2)"]
keys = ["gap_closed_stratum", "gap_closed_random", "gap_closed_n1"]
w = 0.38
for j, p in enumerate(ps):
    h = S[p]["h2"]
    xs = [i + (j - 0.5) * w for i in range(3)]
    ax.bar(xs, [h[k] for k in keys], w * 0.95, color=col[p], label=f"p = {p}")
ax.axhline(1, color="#999", lw=1, ls=":")
ax.set_xticks(range(3), lab, fontsize=6.5)
ax.set_ylim(0, 1.05)
ax.set_ylabel("share of the gap U − C closed\nby the null's median (999 perms)")
ax.set_title("(b) H2: how much clustering each null absorbs", fontsize=8)
ax.legend(fontsize=6.5, loc="upper left")
fig.tight_layout()
fig.savefig(HERE / "serre_sign.svg")
fig.savefig(HERE / "serre_sign.png", dpi=150)

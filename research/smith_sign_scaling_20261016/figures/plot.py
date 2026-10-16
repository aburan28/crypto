"""Layout only: plots values recorded in data/scaling_p*.json and the
2026-10-16 exhaustive bins (research/smith_serre_sign_20261016/data).

python research/smith_sign_scaling_20261016/figures/plot.py
"""

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
plt.rcParams.update({"font.size": 8, "svg.fonttype": "none"})
R = {r["p"]: r for r in (json.loads(f.read_text()) for f in DATA.glob("scaling_p*.json"))}
V = json.loads((DATA / "verdict.json").read_text())
ps = sorted(R)
col = {101: "#e08a1e", 1009: "#3a6fb0", 10007: "#2b8a5a"}

fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.0))
ax = axes[0]
for p in ps:
    b = [d for d in R[p]["bins_t"] if d["n"] >= 2000]
    t = [(d["t_lo"] + d["t_hi"]) / 2 for d in b]
    se = [math.sqrt(d["rate"] * (1 - d["rate"]) / d["n"]) for d in b]
    ax.errorbar(t, [d["rate"] for d in b], yerr=[1.96 * s for s in se], fmt="-o", ms=3,
                capsize=2, color=col[p], label=f"p = {p}")
ax.axhline(0.5, color="#999", lw=1, ls=":")
ax.set_xlabel("normalised trace t = a1 / (2√p)")
ax.set_ylabel("Prop. 4 pass rate, non-degenerate\nProp. 2 passes (95% binomial)")
ax.set_title("(a) the rate flattens toward 1/2 as p grows", fontsize=8)
ax.legend(fontsize=6.5)

ax = axes[1]
x = [p for p in ps]
y = [R[p]["beta"] for p in ps]
e = [1.96 * R[p]["beta_se"] for p in ps]
ax.errorbar(x, y, yerr=e, fmt="o", color="#333", capsize=3, label="β_p, 95% CI (CR1)")
ax.axhline(1, color="#999", lw=1, ls=":", label="β = 1  (N1 / (2p+2))")
x0 = ps[0]
ax.plot(x, [y[0] * math.sqrt(p / x0) for p in x], "--", color="#c44", lw=1,
        label="persist: β ∝ √p from p = 101")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("p")
ax.set_ylabel("slope β of pass rate on a1 / (2p + 2)")
ax.set_title(f"(b) γ = {V['gamma']:.2f} [{V['gamma_ci'][0]:.2f}, {V['gamma_ci'][1]:.2f}]: {V['h3']}",
             fontsize=7)
ax.legend(fontsize=6.5)
fig.tight_layout()
fig.savefig(HERE / "sign_scaling.svg")
fig.savefig(HERE / "sign_scaling.png", dpi=150)

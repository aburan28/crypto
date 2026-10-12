"""Plot recorded values only (run_r*.json, fit.json); no new computation.

Writes fig_total_vs_L.svg, fig_per_target.svg, fig_class_tree.svg.
"""
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

D = Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).parent)
RUNS = ["run_r28.json", "run_r32.json", "run_r36.json"]
COL = {"A": "#c0392b", "B": "#1f6fb2", "C": "#2a9d55"}
MARK = {0: "o", 1: "s", 2: "^"}
runs = [json.loads((D / f).read_text()) for f in RUNS]
fit = json.loads((D / "fit.json").read_text())

# Figure 1: total S vs L.
fig, ax = plt.subplots(figsize=(7.2, 4.8))
for ri, run in enumerate(runs):
    lr = run["base_curve"]["log2_r"]
    summ = run["summary"]
    ls = [s["L"] for s in summ]
    n = summ[0]["trials"]
    for arm in "ABC":
        m = [s[f"{arm}_S"] for s in summ]
        e = [s[f"{arm}_S_sd"] / math.sqrt(s["trials"]) for s in summ]
        x = [l * (1 + 0.04 * (ri - 1)) for l in ls]
        ax.errorbar(x, m, yerr=e, fmt=MARK[ri], color=COL[arm], ms=5, capsize=2,
                    mfc="white" if arm == "C" else COL[arm], lw=1,
                    label=f"arm {arm}, r=2^{lr:.1f} (n={n}/L)")
summ = runs[0]["summary"]
ls = [s["L"] for s in summ]
ax.plot(ls, [s["independent_reference_S_derived"] for s in summ], "--", color="#888",
        lw=1, label=r"derived: $L\sqrt{\pi/2}$ (independent)")
ax.plot(ls, [s["ks_reference_S_derived"] for s in summ], ":", color="#222", lw=1.2,
        label=r"derived: Kuhn–Struik $\sum_{i<L}\binom{2i}{i}4^{-i}\sqrt{\pi/2}$")
ax.set_xscale("log", base=2)
ax.set_yscale("log", base=2)
ax.set_xlabel("number of targets L (curves E_1..E_L of one isogeny class)")
ax.set_ylabel(r"S = total group operations / $\sqrt{r}$")
ax.set_title("Lead 16.5: A independent per curve, B transported to E_0 + multi-target,\n"
             "C L targets planted on E_0 (control). Mean ± s.e., measured.", fontsize=9)
ax.legend(fontsize=6.5, ncol=2, loc="upper left")
ax.grid(alpha=0.3, which="both")
fig.tight_layout()
fig.savefig(D / "fig_total_vs_L.svg")

# Figure 2: per-target cost for L=64, B vs C.
fig, ax = plt.subplots(figsize=(7.2, 4.2))
ks = fit["ks_per_target_S_derived"]
for ri, run in enumerate(runs):
    sr = math.sqrt(run["base_curve"]["r"])
    s64 = [s for s in run["summary"] if s["L"] == 64][0]
    lr = run["base_curve"]["log2_r"]
    for arm, key in (("B", "per_target_steps_B_mean"), ("C", "per_target_steps_C_mean")):
        y = [v / sr for v in s64[key]]
        ax.plot(range(1, 65), y, MARK[ri], color=COL[arm], ms=3,
                mfc="white" if arm == "C" else COL[arm],
                label=f"arm {arm}, r=2^{lr:.1f} (n={s64['trials']})")
ax.plot(range(1, 65), ks, "-", color="#222", lw=1.2,
        label=r"derived: $\sqrt{\pi/2}\binom{2i}{i}4^{-i}$")
ax.set_xscale("log", base=2)
ax.set_yscale("log", base=2)
ax.set_xlabel("target index i in the L=64 multi-target run")
ax.set_ylabel(r"mean walk steps for target i / $\sqrt{r}$")
ax.set_title("Per-target rho steps in the shared-table run (L=64), measured", fontsize=9)
ax.legend(fontsize=6.5, ncol=2)
ax.grid(alpha=0.3, which="both")
fig.tight_layout()
fig.savefig(D / "fig_per_target.svg")

# Figure 3: the class tree used at r=2^28.8 (recorded parents and degrees).
run = runs[0]
cs = run["class"]["curves"]
depth = {c["i"]: c["depth"] for c in cs}
bylevel = {}
for c in cs:
    bylevel.setdefault(c["depth"], []).append(c["i"])
pos = {}
for d, ids in bylevel.items():
    for k, i in enumerate(ids):
        pos[i] = ((k + 1) / (len(ids) + 1), -d)
fig, ax = plt.subplots(figsize=(7.2, 3.6))
for c in cs[1:]:
    (x1, y1), (x2, y2) = pos[c["parent"]], pos[c["i"]]
    ax.plot([x1, x2], [y1, y2], "-", color="#1f6fb2" if c["ell"] == 2 else "#e08a1e", lw=0.8)
for i, (x, y) in pos.items():
    ax.plot(x, y, "o", color="#222" if i == 0 else "#888", ms=5 if i == 0 else 3)
    if i <= 64:
        ax.annotate(str(i), (x, y), fontsize=4.5, xytext=(2, 2), textcoords="offset points")
ax.plot([], [], "-", color="#1f6fb2", label="2-isogeny (Vélu, rational 2-torsion kernel)")
ax.plot([], [], "-", color="#e08a1e", label="3-isogeny (Vélu, rational root of ψ3)")
ax.set_yticks([-d for d in bylevel])
ax.set_yticklabels([f"depth {d}" for d in bylevel])
ax.set_xticks([])
bc = run["base_curve"]
ax.set_title(f"BFS tree of the toy class at r=2^{bc['log2_r']:.1f} (p={bc['p']}); "
             "every edge has a verified dual (dual∘φ=[ℓ])", fontsize=8)
ax.legend(fontsize=6.5, loc="lower right")
fig.tight_layout()
fig.savefig(D / "fig_class_tree.svg")
print("ok")

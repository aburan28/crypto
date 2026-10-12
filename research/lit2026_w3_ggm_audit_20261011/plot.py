"""Plot the W3 audit's recorded values.  Lays out what
`examples/lit2026_w3_ggm_audit.rs` wrote; computes nothing new (the
ratios plotted are fields of table.json / fit.json / measure.json).

usage: python plot.py <run-dir>   (reads table.json, fit.json, measure.json;
writes graph.svg)
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402

run = Path(sys.argv[1])
table = json.loads((run / "table.json").read_text())
fit = json.loads((run / "fit.json").read_text())
meas = json.loads((run / "measure.json").read_text())

INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7"]
plt.rcParams.update({
    "font.size": 8.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "svg.fonttype": "none",
})

fig, axes = plt.subplots(2, 2, figsize=(11, 8.6))
(a, b), (c, d) = axes

# (a) measured vs bound, group-operation families.
fams = ["prime", "char2", "koblitz_round5", "koblitz", "koblitz_batch"]
labels = {"prime": "prime field (Table B)", "char2": "random binary (Table B)",
          "koblitz_round5": "Koblitz round 5 (n ≤ 41)", "koblitz": "Koblitz one target (Table A)",
          "koblitz_batch": "Koblitz 32-target batch (Table C)"}
lo, hi = 1e-1, 1e12
for i, fam in enumerate(fams):
    rows = [r for r in table["rows"] if r["family"] == fam]
    nv = [r for r in rows if r["delta_counting_upper"] < 1]
    vac = [r for r in rows if r["delta_counting_upper"] >= 1]
    a.scatter([r["bound_theta"] for r in nv], [r["measured"] for r in nv], s=22,
              color=SERIES[i], label=f"{labels[fam]} ({len(rows)})", edgecolor="white", linewidth=0.6, zorder=3)
    a.scatter([r["bound_theta"] for r in vac], [r["measured"] for r in vac], s=22,
              facecolor="none", edgecolor=SERIES[i], linewidth=1.0, zorder=3)
a.plot([lo, hi], [lo, hi], color=INK, linewidth=1.2)
a.text(40, 70, "measured = bound", color=INK, rotation=22)
a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(0.5, 1e5); a.set_ylim(1e1, 1e12)
a.set_xlabel("bound min(√(r/A), 1/δ), group operations  (δ = max(counting ceiling, relations/trials))")
a.set_ylabel("measured whole-process cost, GAE (= S·√r·k)")
a.set_title("(a) SGGM bound vs measured cost; hollow = δ_upper = 1 (bound vacuous)", loc="left", fontsize=9)
a.legend(fontsize=7, loc="upper right", frameon=False)

# (b) compact-orbit ladder: probes vs √(r/A), and S·T/(r/A).
rows = [r for r in table["rows"] if r["family"] == "koblitz_compact_orbit"]
tgt = [r for r in rows if r["id"].endswith("target_probes_online")]
tot = [r for r in rows if r["id"].endswith("rank_plus_target_probes")]
ns = [int(r["id"].split("/")[0].split("-n")[1]) for r in tgt]
b.plot(ns, [r["probes_over_sqrt_r_over_aut"] for r in tot], marker="o", color=SERIES[0], linewidth=2,
       markersize=6, label="(rank + target) probes / √(r/A)")
b.plot(ns, [r["probes_over_sqrt_r_over_aut"] for r in tgt], marker="s", color=SERIES[1], linewidth=2,
       markersize=6, label="online target probes / √(r/A)")
b.plot(ns, [r["s_t_over_r_over_a_at_rank_mean"] for r in tgt], marker="^", color=SERIES[2], linewidth=2,
       markersize=6, label="S·T/(r/A), S = K²n index entries, T = rank mean")
b.plot(ns, [r["s_t_over_r_over_a"] for r in tgt], marker="v", color=SERIES[3], linewidth=2,
       markersize=6, label="S·T/(r/A), T = frozen target probes")
b.axhline(1, color=INK, linewidth=1.2)
b.set_yscale("log"); b.set_xlabel("n (Koblitz curve over GF(2^n), A = 2n)")
b.set_ylabel("ratio (probes are S3 solves + lookups, not group ops)")
b.set_title("(b) compact-orbit ladder n37–n83: δ_upper = 1 on every rung", loc="left", fontsize=9)
b.legend(fontsize=7, loc="upper left", frameon=False)

def short(cv):
    name = cv["curve"]
    if name.startswith("prime-floor"):
        return f"prime field, log2 r = {cv['log2_r']:.1f}, A = {cv['aut']:g}"
    if name.startswith("icv1"):
        name = name.split("-tm")[0]
    return f"{name}, A = {cv['aut']:g}"


# (c) multi-target: per-target cost over k = 1 vs k.
mk = ["o", "s", "^", "D", "v"]
i = 0
for src, tag in ((fit, "records"), (meas, "fresh")):
    for cv in src["fits"]:
        ks = [p["k"] for p in cv["per_k"]]
        y = [p["per_target_over_k1"] for p in cv["per_k"]]
        c.plot(ks, y, marker=mk[i % 5], color=SERIES[i], linewidth=2, markersize=6,
               label=f"{tag}: {short(cv)} (slope {cv['slope']:.3f} [{cv['slope_ci95_bootstrap'][0]:.3f}, {cv['slope_ci95_bootstrap'][1]:.3f}])")
        i += 1
law = {p["k"]: p["batch_law"] for src in (fit, meas) for cv in src["fits"] for p in cv["per_k"]}
law_k = sorted(law)
c.plot(law_k, [law[k] for k in law_k], color=INK, linewidth=1.2,
       linestyle="--", label="Kuhn–Struik batch law")
c.set_xscale("log", base=2); c.set_yscale("log")
c.set_xlabel("targets per batch k"); c.set_ylabel("per-target cost / cost at k = 1")
c.set_title("(c) multi-target rho: total ∝ k^slope (generic √(k·r) ⇒ ½)", loc="left", fontsize=9)
c.legend(fontsize=6.5, loc="lower left", frameon=False)

# (d) multiplier closure on the non-adaptive line.
pre = table["precompute_multiplier_closure"]
xs = list(range(len(pre)))
cell_color = {"c23": SERIES[0], "c37": SERIES[1], "c37b": SERIES[2]}
d.bar(xs, [p["canonical_times_t_over_r_over_a"] for p in pre], width=0.75,
      color=[cell_color[p["id"].split("/")[1]] for p in pre], edgecolor="white", linewidth=0.8)
d.axhline(1, color=INK, linewidth=1.2)
d.set_xticks(xs); d.set_xticklabels([p["id"].split("/", 1)[1] for p in pre], rotation=90, fontsize=6)
d.set_ylabel("canonical entries · mean first-hit probes / (r/A)")
d.set_title("(d) multiplier closure: S·T ≈ r/A (non-adaptive line, DKM Thm 1.1)", loc="left", fontsize=9)
from matplotlib.patches import Patch  # noqa: E402

d.legend(handles=[Patch(color=col, label=name.upper()) for name, col in cell_color.items()],
         fontsize=7, frameon=False, loc="upper left")

fig.suptitle("lit2026 W3: IC and precompute configurations against generic lower bounds "
             "(ePrint 2026/384, arXiv 2505.00894, ePrint 2026/2253)", fontsize=10, x=0.01, ha="left")
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(run / "graph.svg")
print("wrote", run / "graph.svg")

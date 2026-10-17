# Layout only: draws raw/summary.json (written by
# examples/gaudry_quartic_d4_exponent_fit.rs). Computes nothing that is
# reported; the fitted lines are drawn from the summary's own coefficients.
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
S = json.loads((HERE.parent / "raw/summary.json").read_text())
C = {"rel": "#2a78d6", "la": "#eb6834", "gauss": "#6250d6", "rho": "#1baf7a", "lp2": "#eb6834", "lp1": "#6250d6", "ink": "#0b0b0b", "muted": "#52514e", "grid": "#e4e3df"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": C["muted"], "axes.labelcolor": C["ink"], "xtick.color": C["muted"], "ytick.color": C["muted"], "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none"})


def save(fig, name):
    fig.tight_layout()
    fig.savefig(HERE / f"{name}.svg")
    fig.savefig(HERE / f"{name}.png", dpi=150)
    plt.close(fig)


def grid(ax):
    ax.grid(True, which="major", color=C["grid"], linewidth=0.6)
    ax.set_axisbelow(True)


def line(ax, fit, ps, color, style="-"):
    xs = [ps[0] * (ps[1] / ps[0]) ** (i / 50) for i in range(51)]
    ax.plot(xs, [math.exp(fit["const_ln"] + fit["slope"] * math.log(x)) for x in xs], style, color=color, linewidth=1.2)


def fmt_ci(f):
    a, b = f["slope_ci95_union"]
    return f"β = {f['slope']:.2f} [{a:.2f}, {b:.2f}]"


# 1. Phases and rho against p (main arm).
w = S["arms"]["ic-d4-w"]
pr = w["per_rung"]
ps = [r["p"] for r in pr]
fig, ax = plt.subplots(figsize=(7.2, 5.4))
series = [
    ("relations", C["rel"], "o", "relations (D₄ solves)"),
    ("linear_algebra", C["la"], "s", "linear algebra, Wiedemann (priced)"),
    ("shadow_gauss", C["gauss"], "^", "dense elimination, same rows (derived, not priced)"),
]
for key, col, mk, lab in series:
    for r in pr:
        ax.scatter([r["p"]] * len(r[f"{key}_values"]), r[f"{key}_values"], s=9, color=col, alpha=0.35, linewidths=0)
    ax.scatter(ps, [r[key] for r in pr], s=26, color=col, marker=mk, label=f"{lab}: {fmt_ci(w['fits'][key])}", zorder=3)
    line(ax, w["fits"][key], (ps[0], ps[-1]), col)
rho = S["rho"]
rps = [r["p"] for r in rho["per_rung"]]
ax.scatter(rps, [r["gae"] for r in rho["per_rung"]], s=26, color=C["rho"], marker="D", zorder=3,
           label=f"rho (negation), mean of {rho['records']} runs: β = {rho['fit_gae']['slope']:.2f}")
line(ax, rho["fit_gae"], (rps[0], rps[-1]), C["rho"], "--")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("p  (curve over F_{p⁴})")
ax.set_ylabel("cost (GAE, group-addition equivalents)")
ax.set_title("k = 4 index calculus by phase, and rho, against p (means per rung; dots: targets)", fontsize=9, color=C["ink"])
grid(ax)
ax.legend(fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
save(fig, "phases_vs_p")

# 2. Extrapolation to the crossover.
fig, ax = plt.subplots(figsize=(7.2, 4.4))
pmax = 2 ** 34
xs = [29 * (pmax / 29) ** (i / 200) for i in range(201)]
f = w["fits"]
tot = [sum(math.exp(f[k]["const_ln"] + f[k]["slope"] * math.log(x)) for k in ("relations", "linear_algebra", "other")) for x in xs]
gtot = [sum(math.exp(f[k]["const_ln"] + f[k]["slope"] * math.log(x)) for k in ("relations", "shadow_gauss", "other")) for x in xs]
ax.plot(xs, tot, color=C["rel"], linewidth=1.6, label="IC, sum of fitted phases (Wiedemann)")
ax.plot(xs, gtot, color=C["gauss"], linewidth=1.2, linestyle=":", label="IC with the dense elimination priced instead")
for h, st in ((4, "--"), (16, "-.")):
    ax.plot(xs, [rho["s_rho_mean"] * math.sqrt(x ** 4 / h) for x in xs], st, color=C["rho"], linewidth=1.2, label=f"rho, S̄ = {rho['s_rho_mean']:.2f}, h = {h}")
ax.axvspan(29, ps[-1], color=C["grid"], alpha=0.6, linewidth=0)
ax.text(math.sqrt(29 * ps[-1]), 1e6, "measured", ha="center", fontsize=8, color=C["muted"])
for c in w["crossovers"]:
    if c["variant"] == "wiedemann" and c["rho"] == "h=4" and c["log2_p_star"]:
        x = 2 ** c["log2_p_star"]
        ax.axvline(x, color=C["muted"], linewidth=0.8)
        lo, hi = c["log2_p_star_boot95"]
        ax.axvspan(2 ** lo, 2 ** hi, color=C["rel"], alpha=0.12, linewidth=0)
        ax.text(x, 1e5, f" measured fit, h = 4:\n p* = 2^{c['log2_p_star']:.1f}, group 2^{c['log2_group_star']:.0f}", fontsize=7, color=C["ink"])
sens = [c["log2_group_star"] for c in w["crossover_relations_slope_sensitivity"] if c["h"] == 4.0 and c["log2_group_star"]]
if len(sens) == 2:
    ax.axvspan(2 ** (min(sens) / 4), 2 ** (max(sens) / 4), facecolor="none", alpha=0.5, linewidth=0, hatch="////", edgecolor=C["rel"])
    ax.text(2 ** (min(sens) / 4), 1e20, f" relations slope at its t-interval ends:\n group 2^{min(sens):.0f}…2^{max(sens):.0f}", fontsize=7, color=C["muted"])
m = S["model_11_16_from_pr1863"]["log2_group_star"]["h=4"]
ax.axvspan(2 ** (m[0] / 4), 2 ** (m[1] / 4), color=C["la"], alpha=0.25, linewidth=0)
ax.text(2 ** (m[0] / 4), 3e1, f"§11.16 model, h = 4:  \n group 2^{m[0]:.1f}…2^{m[1]:.1f}  ", ha="right", va="bottom", fontsize=7, color=C["ink"])
ax.set_xscale("log", base=2)
ax.set_yscale("log")
ax.set_xlabel("p (log₂ scale)")
ax.set_ylabel("cost (GAE)")
ax.set_title("Crossover with rho from the measured fits (extrapolated beyond the shaded range)", fontsize=9, color=C["ink"])
grid(ax)
ax.legend(fontsize=7, frameon=False, loc="upper left")
save(fig, "crossover")

# 3. Large primes against the plain arm.
fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.6))
for arm, col, mk, lab in (("ic-d4-w", C["rel"], "o", "no large prime"), ("ic-d4-lp2", C["lp2"], "s", "two large primes"), ("ic-d4-lp1", C["lp1"], "^", "one large prime")):
    if arm not in S["arms"]:
        continue
    a = S["arms"][arm]
    q = [r["p"] for r in a["per_rung"]]
    for ax, key in zip(axs, ("linear_algebra", "total")):
        ax.scatter(q, [r[key] for r in a["per_rung"]], s=22, color=col, marker=mk, zorder=3, label=f"{lab}: {fmt_ci(a['fits'][key])}")
        line(ax, a["fits"][key], (q[0], q[-1]), col)
for ax, t in zip(axs, ("linear algebra (Wiedemann)", "total cost")):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("p")
    ax.set_ylabel("GAE")
    ax.set_title(t, fontsize=9, color=C["ink"])
    grid(ax)
    ax.legend(fontsize=6.5, frameon=False, loc="upper left")
save(fig, "large_primes")

# 4. Slopes with intervals.
rows = []
for arm, lab in (("ic-d4-w", "plain"), ("ic-d4-lp2", "LP2"), ("ic-d4-lp1", "LP1")):
    if arm not in S["arms"]:
        continue
    for key, kl in (("relations", "relations"), ("linear_algebra", "lin. alg. (Wiedemann)"), ("shadow_gauss", "dense elimination"), ("total", "total")):
        fk = S["arms"][arm]["fits"].get(key)
        if fk:
            rows.append((f"{lab}: {kl}", fk["slope"], fk["slope_ci95_union"]))
rows.append(("rho", rho["fit_gae"]["slope"], [min(rho["fit_gae"]["slope_t95"][0], rho["fit_gae"]["slope_boot95"][0]), max(rho["fit_gae"]["slope_t95"][1], rho["fit_gae"]["slope_boot95"][1])]))
fig, ax = plt.subplots(figsize=(7.0, 0.32 * len(rows) + 1.0))
for i, (lab, b, (lo, hi)) in enumerate(rows):
    ax.plot([lo, hi], [i, i], color=C["muted"], linewidth=1.2)
    ax.scatter([b], [i], s=24, color=C["rel"], zorder=3)
    ax.text(hi + 0.05, i, f"{b:.2f}", va="center", fontsize=7, color=C["ink"])
for v in (1, 1.5, 2, 3):
    ax.axvline(v, color=C["grid"], linewidth=0.8, zorder=0)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([r[0] for r in rows], fontsize=7.5)
ax.invert_yaxis()
ax.set_xlabel("fitted exponent β in cost ∝ p^β (95 % interval: union of OLS-t and bootstrap)")
save(fig, "slopes")
print("ok", file=sys.stderr)

# 5. Run report: every recorded S, by arm, against p.
RUN = HERE.parent / "reports/run-20261015T0457Z"
RUN.mkdir(parents=True, exist_ok=True)
fig, ax = plt.subplots(figsize=(7.2, 4.2))
for arm, col, mk, lab in (("ic-d4-w", C["rel"], "o", "ic-d4-w"), ("ic-d4-lp2", C["lp2"], "s", "ic-d4-lp2"), ("ic-d4-lp1", C["lp1"], "^", "ic-d4-lp1")):
    if arm not in S["arms"]:
        continue
    xs, ys = [], []
    for r in S["arms"][arm]["per_rung"]:
        xs += [r["p"]] * len(r["s_values"])
        ys += r["s_values"]
    ax.scatter(xs, ys, s=16, color=col, marker=mk, label=f"{lab} ({len(ys)} runs)", zorder=3)
ax.scatter([r["p"] for r in rho["per_rung"]], [r["s"] for r in rho["per_rung"]], s=16, color=C["rho"], marker="D", label=f"rho-neg, rung mean ({rho['records']} distinct runs)")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("p")
ax.set_ylabel("S (GAE / √r), counted")
ax.set_title("Every verified run's S (sessions under sessions/; summary raw/summary.json)", fontsize=9, color=C["ink"])
grid(ax)
ax.legend(fontsize=7, frameon=False, loc="center right")
fig.tight_layout()
fig.savefig(RUN / "s_by_run.svg")
fig.savefig(RUN / "s_by_run.png", dpi=150)
plt.close(fig)

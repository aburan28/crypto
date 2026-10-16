"""Figures for research/bitop_unit_20261016: recorded values only.

Every plotted number is read from a file ecbench wrote (bound records,
session records with their bitop/v1 work block, audit-ir output, the
archive summary); nothing here computes a result.  Run from the repository
root:

    python research/bitop_unit_20261016/figures/plots.py
"""

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path("research/bitop_unit_20261016")
OUT = ROOT / "figures"
# Reference palette, light mode, first three slots (validated all-pairs).
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update(
    {
        "font.size": 9,
        "axes.edgecolor": INK2,
        "axes.labelcolor": INK,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "svg.fonttype": "none",
    }
)


def save(fig, name):
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def bound(name):
    return json.loads((ROOT / "bounds" / f"{name}.json").read_text())


# 1. Exponents in additions, with their intervals.
rows = []
for study, arms in (
    ("prime ladder", ["rho-neg", "bsgs-neg", "kangaroo", "ic-mitm-b32", "ic-mitm-b128", "ic-rootfind-m3-b32"]),
    ("fiber ladder", ["rho-neg", "closed-b16", "closed-b32", "closed-b64", "closed-b128", "closed-b256", "blind-b128"]),
):
    prefix = "prime-ladder" if study == "prime ladder" else "fiber"
    for arm in arms:
        b = bound(f"{prefix}-{arm}")
        f = b["fit"]
        rows.append((f"{study}: {arm}", f["alpha_in_additions"], f["alpha_in_additions_ci95"], f["sizes"]))
fig, ax = plt.subplots(figsize=(6.4, 0.32 * len(rows) + 0.9))
for i, (label, a, ci, sizes) in enumerate(reversed(rows)):
    ax.plot(ci, [i, i], color=BLUE, linewidth=2, solid_capstyle="round")
    ax.plot([a], [i], "o", color=BLUE, markersize=5, markeredgecolor="white", markeredgewidth=1)
ax.axvline(0.5, color=INK2, linewidth=1, linestyle=(0, (3, 3)))
ax.text(0.505, -0.9, "generic theory: ½", color=INK2, va="top", fontsize=8)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([r[0] for r in reversed(rows)])
ax.set_xlabel("exponent α of bitops / gates(one addition) = C·r^α, with 95% interval")
ax.set_title("Fitted exponents in bitop/v1, in additions\n(sessions ECBS1h3674c4a0b139, ECBS1hef866ef12457)", loc="left", fontsize=9)
ax.set_ylim(-1.6, len(rows) - 0.5)
ax.grid(axis="y", visible=False)
save(fig, "exponents")

# 2. Every run of the prime ladder: gates against r.
recs = [json.loads(l) for l in (ROOT / "sessions/prime-ladder/records.jsonl").read_text().splitlines()]
fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), sharey=True)
for ax, arms, title in (
    (axes[0], ["rho-neg", "bsgs-neg", "kangaroo"], "generic methods"),
    (axes[1], ["ic-mitm-b32", "ic-mitm-b128", "ic-rootfind-m3-b32"], "index calculus (fixed base sizes)"),
):
    for arm, color in zip(arms, (BLUE, ORANGE, AQUA)):
        pts = [
            (math.log2(r["workload"]["curve"]["r"]), r["work"]["log2_bitops"])
            for r in recs
            if r["arm"] == arm and r["outcome"]["status"] == "verified" and r.get("work")
        ]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o", color=color, markersize=4, alpha=0.75,
                markeredgecolor="white", markeredgewidth=0.6, label=arm)
    ax.set_title(title, loc="left", fontsize=9)
    ax.set_xlabel("log₂ r")
    ax.legend(frameon=False, fontsize=8)
axes[0].set_ylabel("log₂ bitops of one verified run")
fig.suptitle("Prime ladder, 16–26 bits, 8 targets a size, every run (session ECBS1h3674c4a0b139)", x=0.01, ha="left", fontsize=9)
save(fig, "ladder_runs")

# 3. Instruction audit: each arm against rho on the same curve.
audit = [json.loads(l) for l in (ROOT / "audit/final-panel.jsonl").read_text().splitlines()]
by_curve = {}
for a in audit:
    by_curve.setdefault(a["curve"], []).append(a)
labels, arm_v, rho_v = [], [], []
for curve, items in by_curve.items():
    rhos = [x for x in items if x["method"].startswith("rho.")]
    if not rhos:
        continue
    ref = min(x["ir_per_kilo_bitop"] for x in rhos)
    for x in items:
        if x["method"].startswith("rho."):
            continue
        what = x["params"].get("oracle", "").split(":")[0]
        solver = x["params"].get("solver", "").split(":")[0]
        engine = [p for p in x["params"].get("oracle", "").split(",") if p.startswith("engine=")]
        tag = solver or (engine[0][7:] if engine else what)
        labels.append(f"{tag} · {curve.split('-')[1]} · seq {x['seq']}")
        arm_v.append(x["ir_per_kilo_bitop"])
        rho_v.append(ref)
fig, ax = plt.subplots(figsize=(6.4, 0.34 * len(labels) + 0.9))
ys = list(range(len(labels)))
ax.plot(rho_v, ys, "s", color=ORANGE, markersize=6, label="rho, same curve", markeredgecolor="white")
ax.plot(arm_v, ys, "o", color=BLUE, markersize=6, label="index-calculus arm", markeredgecolor="white")
ax.set_xscale("log")
ax.set_yticks(ys)
ax.set_yticklabels(labels)
ax.set_xlabel("Callgrind instructions per thousand counted bit operations (log scale)")
ax.set_title("Instruction audit: a ratio near rho's means the counters see the work", loc="left", fontsize=9)
ax.legend(frameon=False, fontsize=8, loc="lower right")
ax.grid(axis="y", visible=False)
save(fig, "audit")

# 4. Old unit against new, per audited arm (same execution): the committed
# record's GAE ratio to the floor and the audit replay's bitop/v1 ratio.
old_new = []
for x in audit:
    if x["method"].startswith("rho."):
        continue
    old_new.append((x, x.get("committed_ratio_to_floor")))
pairs = [(x, o) for x, o in old_new if o]
if pairs:
    fig, ax = plt.subplots(figsize=(6.4, 0.34 * len(pairs) + 0.9))
    ys = list(range(len(pairs)))
    ax.plot([o for _, o in pairs], ys, "s", color=ORANGE, markersize=6, label="committed record, old unit (GAE)", markeredgecolor="white")
    ax.plot([x["work"]["ratio_to_floor"] for x, _ in pairs], ys, "o", color=BLUE, markersize=6, label="bitop/v1, complete counters", markeredgecolor="white")
    ax.set_xscale("log")
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{(x['params'].get('solver') or x['params'].get('oracle','')).split(':')[0]} · seq {x['seq']}" for x, _ in pairs])
    ax.set_xlabel("ratio to the generic floor √(πr/2A) (log scale)")
    ax.set_title("What the old unit left out: the same executions, priced both ways\n(table oracles: old and new ratios within 15% of each other, so their marks nearly coincide)", loc="left", fontsize=9)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.grid(axis="y", visible=False)
    save(fig, "old_vs_new")
print("figures written to", OUT)

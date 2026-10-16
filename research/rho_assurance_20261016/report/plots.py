"""Lay out recorded values from ../summary.json and the Callgrind function
summaries as SVG (and 150 dpi PNG) graphs.  Plotting only: every number is
read from a committed record; nothing is computed beyond per-step division of
recorded instruction counts by recorded step counts."""
import json
import pathlib
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

here = pathlib.Path(__file__).resolve().parent
study = here.parent
summary = json.loads((study / "summary.json").read_text())


def save(fig, name):
    fig.savefig(here / f"{name}.svg")
    fig.savefig(here / f"{name}.png", dpi=150)
    plt.close(fig)


# 1. Q2: steps / ideal per fixture, mean and 95% interval, per size.
q2 = summary["q2_step_counts"]
fig, ax = plt.subplots(figsize=(7, 4))
for i, cell in enumerate(q2):
    r = cell["ratios"]
    xs = [i + (j % 40 - 20) * 0.006 for j in range(len(r))]
    ax.scatter(xs, r, s=4, alpha=0.25, color="#4c72b0")
    lo, hi = cell["ratio"]["ci95"]
    ax.errorbar(i + 0.3, cell["ratio"]["mean"], yerr=[[cell["ratio"]["mean"] - lo], [hi - cell["ratio"]["mean"]]],
                fmt="o", color="black", capsize=4)
    ax.errorbar(i + 0.4, cell["visited_ratio"]["mean"],
                yerr=[[cell["visited_ratio"]["mean"] - cell["visited_ratio"]["ci95"][0]],
                      [cell["visited_ratio"]["ci95"][1] - cell["visited_ratio"]["mean"]]],
                fmt="s", color="#dd8452", capsize=4)
ax.axhline(1.0, color="grey", lw=1)
ax.axhline(1.10, color="#c44e52", lw=1, ls="--")
ax.set_xticks(range(len(q2)))
ax.set_xticklabels([f"n = {c['n']}\n{c['slug']}\n({c['fixtures']} fixtures)" for c in q2], fontsize=7)
ax.set_ylabel("walk steps / √(πr / 2A)")
ax.set_ylim(0, 3)
ax.set_title("Strong rho step counts against the ideal walk (Q2)", fontsize=10)
ax.text(0.01, 0.97, "dots: fixtures · black: mean, 95% CI of walk steps · "
        "orange: mean, 95% CI of steps + walk starts\n"
        "grey: ideal (1.0) · red dashed: declared 1.10 limit",
        transform=ax.transAxes, fontsize=7, va="top")
fig.tight_layout()
save(fig, "q2-step-counts")

# 2. Q3/Q4: instructions per walk step by function, per arm (n53, 4 fixtures).
groups = [
    ("canonicalize", r"::canonicalize"),
    ("solve body (hashes, lanes, point add)", r"::solve"),
    ("batch inversion", r"batch_inv"),
    ("field multiply (calls)", r"mul_fold"),
    ("inversion (walk starts)", r"(::inv \[|::inv$|inv_clmul)"),
    ("table hashing (SipHash)", r"(sip::Hasher|DefaultHasher|RandomState)"),
    ("128-bit modulo", r"__umodti3"),
]
arms = [("baseline", "baseline"), ("apply2 only", "candidate"), ("final", "candidate2")]
rows = {}
for label, arm in arms:
    text = (study / f"runs/q4/callgrind/fixture-n53-4-{arm}.annotate.txt").read_text()
    total = int(re.search(r"([\d,]+) \(100\.0%\)\s+PROGRAM TOTALS", text).group(1).replace(",", ""))
    steps = next(a["walk_steps"] for a in summary["q4_instructions"]["arms"] if a["arm"] == arm)
    parts = {g: 0 for g, _ in groups}
    for line in text.splitlines():
        m = re.match(r"\s*([\d,]+) \(\s*[\d.]+%\)\s+(\S.*)", line)
        if not m or "PROGRAM TOTALS" in line:
            continue
        ir = int(m.group(1).replace(",", ""))
        for g, pat in groups:
            if re.search(pat, m.group(2)):
                parts[g] += ir
                break
    parts["other"] = total - sum(parts.values())
    rows[label] = {k: v / steps for k, v in parts.items()}
fig, ax = plt.subplots(figsize=(7, 4.8))
names = [g for g, _ in groups] + ["other"]
colors = ["#4c72b0", "#55a868", "#c44e52", "#8172b2", "#ccb974", "#64b5cd", "#937860", "#bbbbbb"]
for i, (label, _) in enumerate(arms):
    left = 0
    for name, color in zip(names, colors):
        v = rows[label][name]
        ax.barh(i, v, left=left, color=color, label=name if i == 0 else None)
        left += v
    ax.text(left + 10, i, f"{left:.0f}", va="center", fontsize=8)
ax.set_yticks(range(len(arms)))
ax.set_yticklabels([a for a, _ in arms])
ax.invert_yaxis()
ax.set_xlabel("instructions per walk step (Callgrind Ir / walk_steps)")
ax.set_title("Strong rho, n = 53, 4 fixtures: where a step's instructions go (Q3, Q4)", fontsize=10)
ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3)
fig.tight_layout()
save(fig, "q4-instructions-per-step")

# 3. Q4 native time per step, when recorded.
t = summary.get("q4_native_time")
if t:
    fig, ax = plt.subplots(figsize=(7, 3.6))
    for i, cell in enumerate(t["cells"]):
        for j, (arm, color) in enumerate([("baseline_ns", "#4c72b0"), ("candidate_ns", "#55a868")]):
            s = cell[arm]
            ax.bar(i * 3 + j, s["median"], color=color, label=arm.split("_")[0] if i == 0 else None)
            ax.errorbar(i * 3 + j, s["median"], yerr=[[s["median"] - s["min"]], [s["max"] - s["median"]]],
                        color="black", capsize=4)
        ax.text(i * 3 + 0.5, max(cell["baseline_ns"]["max"], cell["candidate_ns"]["max"]) * 1.03,
                f"×{cell['ratio']:.3f} (A/A {100 * cell['aa_spread']:.1f}%)", ha="center", fontsize=8)
    ax.set_xticks([i * 3 + 0.5 for i in range(len(t["cells"]))])
    ax.set_xticklabels([c["slug"] for c in t["cells"]], fontsize=7)
    ax.set_ylim(0, max(c["baseline_ns"]["max"] for c in t["cells"]) * 1.15)
    ax.set_ylabel("ns per walk step (median, range)")
    ax.set_title("Strong rho native time per step, 5 A/B rounds under isolated_bench, uncontended runs (Q4)", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    save(fig, "q4-native-time")

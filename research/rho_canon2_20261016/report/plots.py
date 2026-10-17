"""Lay out recorded values from ../summary.json and the Callgrind function
summaries as SVG (and 150 dpi PNG) graphs.  Plotting only: every number is
read from a committed record; the only arithmetic is recorded instruction
counts divided by recorded walk-step counts."""
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


# 1. Instructions per walk step by arm.
arms = summary["instructions"]["arms"]
labels = {"baseline": "baseline (main)", "candidate1-4": "changes 1-4",
          "candidate": "changes 1-5 (Shoup)", "final": "changes 1-6 (inlined)"}
fig, ax = plt.subplots(figsize=(7, 3.4))
ys = range(len(arms))
ax.barh(list(ys), [a["ir_per_step"] for a in arms], color=["#4c72b0", "#8fb0d9", "#8fb0d9", "#55a868"])
for i, a in enumerate(arms):
    ax.text(a["ir_per_step"] + 8, i, f"{a['ir_per_step']:.0f}", va="center", fontsize=8)
ax.set_yticks(list(ys))
ax.set_yticklabels([labels[a["arm"]] for a in arms])
ax.invert_yaxis()
ax.set_xlim(0, 1050)
ax.set_xlabel("instructions per walk step (Callgrind Ir / walk_steps)")
ax.set_title("Strong rho, n = 53, 4 fixtures: instructions per step, outputs identical", fontsize=10)
fig.tight_layout()
save(fig, "instructions-per-step")

# 2. Where the step goes, baseline vs final.
groups = [
    ("canonicalize", r"::canonicalize"),
    ("solve body (incl. inlined canonicalize in final)", r"::solve"),
    ("batch inversion", r"batch_inv"),
    ("field multiply (calls)", r"mul_fold"),
    ("inversion (walk starts)", r"(::inv \[|::inv\s*$|inv_clmul)"),
]
rows = {}
for arm in ["baseline", "final"]:
    text = (study / f"runs/callgrind/fixture-n53-4-{arm}.annotate.txt").read_text()
    total = int(re.search(r"([\d,]+) \(100\.0%\)\s+PROGRAM TOTALS", text).group(1).replace(",", ""))
    steps = next(a["walk_steps"] for a in arms if a["arm"] == arm)
    parts = {g: 0 for g, _ in groups}
    for line in text.splitlines():
        m = re.match(r"\s*([\d,]+) \(\s*[\d.]+%\)\s+(\S.*)", line)
        if not m or "PROGRAM TOTALS" in line:
            continue
        for g, pat in groups:
            if re.search(pat, m.group(2)):
                parts[g] += int(m.group(1).replace(",", ""))
                break
    parts["other"] = total - sum(parts.values())
    rows[arm] = {k: v / steps for k, v in parts.items()}
fig, ax = plt.subplots(figsize=(7, 3.8))
names = [g for g, _ in groups] + ["other"]
colors = ["#4c72b0", "#55a868", "#c44e52", "#8172b2", "#ccb974", "#bbbbbb"]
for i, arm in enumerate(["baseline", "final"]):
    left = 0
    for name, color in zip(names, colors):
        ax.barh(i, rows[arm][name], left=left, color=color, label=name if i == 0 else None)
        left += rows[arm][name]
    ax.text(left + 8, i, f"{left:.0f}", va="center", fontsize=8)
ax.set_yticks([0, 1])
ax.set_yticklabels(["baseline", "final"])
ax.invert_yaxis()
ax.set_xlabel("instructions per walk step")
ax.set_title("Where a step's instructions go (n = 53, 4 fixtures)", fontsize=10)
ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
fig.tight_layout()
save(fig, "instructions-by-function")

# 3. Native time per step.
t = summary["native_time"]
fig, ax = plt.subplots(figsize=(7, 3.6))
for i, cell in enumerate(t["cells"]):
    for j, (arm, color) in enumerate([("baseline_ns", "#4c72b0"), ("candidate_ns", "#55a868")]):
        s = cell[arm]
        ax.bar(i * 3 + j, s["median"], color=color, label=("baseline" if j == 0 else "final") if i == 0 else None)
        ax.errorbar(i * 3 + j, s["median"], yerr=[[s["median"] - s["min"]], [s["max"] - s["median"]]], color="black", capsize=4)
    ax.text(i * 3 + 0.5, max(cell["baseline_ns"]["max"], cell["candidate_ns"]["max"]) * 1.03,
            f"x{cell['ratio']:.3f} (A/A {100 * cell['aa_spread']:.1f}%)", ha="center", fontsize=8)
ax.set_ylim(0, max(c["baseline_ns"]["max"] for c in t["cells"]) * 1.18)
ax.set_xticks([i * 3 + 0.5 for i in range(len(t["cells"]))])
ax.set_xticklabels([c["slug"] for c in t["cells"]], fontsize=7)
ax.set_ylabel("ns per walk step (median, range)")
ax.set_title(f"Native time, 5 A/B rounds, isolated_bench ({t['isolation']['runs']} runs, "
             f"{len(t['isolation']['contended'])} contended)", fontsize=10)
ax.legend(fontsize=8, loc="lower right")
fig.tight_layout()
save(fig, "native-time")

"""Draw the recorded values of the F4/F5/F6 cross-family round as SVG.

Plots only numbers already written to ../../ANALYSIS.json by the native
analyser (examples/f4f5f6_cross_family_analyze.rs): per-arm medians,
minima and maxima, and the per-target values it lists.  Nothing is
computed here beyond choosing what to draw.
Usage: python plot.py   (from this directory)
"""
import json
import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt

A = json.load(open("../../ANALYSIS.json"))
ENGINES = [("matrix-f4", "matrix F4", "#7a7a7a"),
           ("matrix-f5", "matrix F5", "#b04a8a"),
           ("inherited-f4", "inherited F4", "#2a6fb0"),
           ("f6-ic", "F6-IC", "#d0702a")]
RHOS = [("rho-strong", "strong rho (Frobenius)", "#3a9a5a"),
        ("rho-negation", "rho, negation only", "#9ab83a")]
CELLS = [("K13", "Koblitz m=13"), ("R13", "random binary m=13"),
         ("K19", "Koblitz m=19"), ("R19", "random binary m=19")]


def cell(session, family=None):
    s = A["sessions"].get(session)
    if not s:
        return None
    for c in s["cells"]:
        if family is None or c["family"] == family:
            return c
    return None


def x19(family):
    return cell("X19", family)


def stat(arms, arm, field):
    a = arms.get(arm)
    if not a or not a.get(field):
        return None
    return a[field]


# 1. Solver word XORs per engine and cell (median, min-max bars).
fig, ax = plt.subplots(figsize=(9, 4.4))
x, ticks, labels = 0, [], []
for key, title in CELLS:
    c = cell(key)
    if c is None:
        continue
    arms = dict(c["arms"])
    extra = x19("koblitz" if key == "K19" else "binary") if key.endswith("19") else None
    start = x
    for arm, lab, col in ENGINES:
        src = arms
        if arm not in arms and extra is not None:
            src = extra["arms"]
        v = stat(src, arm, "word_xors")
        if v is not None:
            ax.bar(x, v["median"], color=col, width=0.8,
                   hatch="//" if src is not arms else None, edgecolor="white")
            ax.errorbar(x, v["median"], yerr=[[v["median"] - v["min"]], [v["max"] - v["median"]]],
                        color="black", capsize=2, linewidth=0.8)
            ax.text(x, v["max"] * 1.3, f"n={v['n']}", ha="center", fontsize=6)
        elif src.get(arm):
            st = src[arm]["statuses"]
            ax.text(x, 1e6, "\n".join(f"{k} {n}" for k, n in st.items()), ha="center", fontsize=6, rotation=90)
        x += 1
    ticks.append((start + x - 1) / 2)
    labels.append(title)
    x += 1.2
ax.set_yscale("log")
ax.set_xticks(ticks, labels, fontsize=8)
ax.set_ylabel("solver word XORs per run (counted; log scale)")
ax.set_title("Counted solver work per one-target run, median and range\n"
             "8 public targets per cell; hatched = SPEC-X19 (1 target, 1800 s cap)", fontsize=10)
for _, lab, col in ENGINES:
    ax.bar(0, 0, color=col, label=lab)
ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=4)
fig.tight_layout()
fig.savefig("solver_word_xors.svg", bbox_inches="tight")

# 2. F6-IC against inherited F4: per-target log2(W_inherited / W_F6).
fig, ax = plt.subplots(figsize=(7.5, 3.8))
for i, (key, title) in enumerate(CELLS):
    c = cell(key)
    q = c and c["questions"].get("q3_log2_inherited_over_f6_word_xors")
    if not q:
        continue
    vals = q["log2_ratio"]["values"]
    ax.scatter([i + (j - len(vals) / 2) * 0.04 for j in range(len(vals))], vals,
               color="#d0702a", s=16, zorder=3)
    ax.hlines(q["log2_ratio"]["median"], i - 0.25, i + 0.25, color="black", linewidth=1.5)
    ax.text(i, min(vals) - 0.02, q["verdict"], ha="center", va="top", fontsize=7)
ax.axhline(0.1, color="#444", linestyle="--", linewidth=0.8)
ax.axhspan(-0.05, 0.05, color="#dddddd", zorder=0)
ax.text(3.45, 0.105, "rule: median >= 0.1 and all > 0", ha="right", va="bottom", fontsize=7)
ax.text(3.45, 0.0, "'no saving' band +/-0.05", ha="right", va="center", fontsize=7)
ax.set_xticks(range(len(CELLS)), [t for _, t in CELLS], fontsize=8)
ax.set_ylabel("log2( W inherited F4 / W F6-IC )")
ax.set_title("F6-IC against inherited F4, per target (dots) and median (bar)\n"
             "W = counted solver word XORs; same public workloads", fontsize=10)
fig.tight_layout()
fig.savefig("f6_vs_inherited_log2.svg", bbox_inches="tight")

# 3. Lower-bound S per arm against rho and the floor.
fig, ax = plt.subplots(figsize=(9, 4.4))
x, ticks, labels = 0, [], []
for key, title in CELLS:
    c = cell(key)
    if c is None:
        continue
    start = x
    for arm, lab, col in RHOS + ENGINES:
        v = stat(c["arms"], arm, "s")
        if v is None:
            continue
        ax.bar(x, v["median"], color=col, width=0.8)
        ax.errorbar(x, v["median"], yerr=[[v["median"] - v["min"]], [v["max"] - v["median"]]],
                    color="black", capsize=2, linewidth=0.8)
        x += 1
    ax.hlines(c["floor_s"], start - 0.5, x - 0.5, color="red", linestyle=":", linewidth=1.2)
    ticks.append((start + x - 1) / 2)
    labels.append(title)
    x += 1.2
ax.set_yscale("log")
ax.set_xticks(ticks, labels, fontsize=8)
ax.set_ylabel("S = charged GAE / sqrt(r) (log scale)")
ax.set_title("Cost per one-target run in ecbench's unit, median and range (8 targets)\n"
             "IC S is a lower bound (word XORs unpriced); red dotted line = floor of the family", fontsize=10)
for _, lab, col in RHOS + ENGINES:
    ax.bar(0, 0, color=col, label=lab)
ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=6)
fig.tight_layout()
fig.savefig("s_per_arm.svg", bbox_inches="tight")
print("wrote solver_word_xors.svg f6_vs_inherited_log2.svg s_per_arm.svg")

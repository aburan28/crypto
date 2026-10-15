"""Layout only: draws raw/summary.json (written by the Rust example
gaudry_k4_crossover_summary) and #1883's raw/summary.json. No result is
computed here."""
import json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
S = json.load(open(os.path.join(STUDY, "raw", "summary.json")))
OLD = json.load(open(os.path.join(STUDY, "..", "gaudry_quartic_d4_char2_20261015", "raw", "summary.json")))
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3})
BLUE, ORANGE, GREEN, PURPLE, GREY = "#2a78d6", "#e8673a", "#1fa37a", "#8b5cf6", "#555555"


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, name + ".svg"))
    fig.savefig(os.path.join(HERE, name + ".png"), dpi=150)
    plt.close(fig)


# 1. IC / rho against log2 n: #1883 (dense matrix, m = 6..14) and this study.
fig, ax = plt.subplots(figsize=(7.6, 4.6))
old = [r for r in OLD["sessions"]["ic"] if r["arm"] == "ic-frobenius" and r["session"] != "char2-forms-m6"]
ax.scatter([r["log2_n"] for r in old], [r["ratio_to_rho_theory"] for r in old], s=28, color=GREY,
           label="#1883: dense matrix, m = 6–14 (h = 2–14)")
E = S["part_e"]["targets"]
ok = [t for t in E if t["status"] == "verified"]
ax.scatter([t["log2_n"] for t in ok], [t["ic_over_rho_theory"] for t in ok], s=46, color=BLUE, zorder=3,
           label="this study: wiedemann, h = 2, ÷ rho theory (0.886)")
meas = [t for t in ok if t.get("ic_over_rho_measured")]
ax.scatter([t["log2_n"] for t in meas], [t["ic_over_rho_measured"] for t in meas], s=46, marker="x",
           color=ORANGE, zorder=3, label="this study, ÷ the measured rho run on the same target")
for t in ok:
    ax.annotate(f"m={t['m']}", (t["log2_n"], t["ic_over_rho_theory"]), textcoords="offset points",
                xytext=(5, 4), fontsize=8)
ax.axhline(1.0, color="black", lw=1)
ax.text(ax.get_xlim()[0] + 0.5, 1.08, "rho", fontsize=9)
c = S["part_e"].get("crossover")
if c:
    ax.axvline(c["log2_n_star"], color=GREEN, ls="--", lw=1)
    ax.text(c["log2_n_star"] + 0.3, 6, f"measured crossover\nlog2 n* = {c['log2_n_star']:.1f}", color=GREEN, fontsize=8)
ax.set_yscale("log")
ax.set_xlabel("log2 group order n")
ax.set_ylabel("index calculus ÷ rho (S)")
ax.set_title("k = 4 in characteristic 2: IC against rho, measured")
ax.legend(fontsize=8, loc="upper right")
save(fig, "crossover")

# 2. Linear-algebra row_ops against columns.
fig, ax = plt.subplots(figsize=(6.8, 4.4))
dense, wied = [], []
for r in S["runs"]:
    if r["matrix_work"] is None or not str(r["arm"]).startswith("ic-"):
        continue
    cols = 2 ** r["m"] / 4
    (dense if r["arm"] == "ic-dense" else wied).append((cols, r["matrix_work"]))
for pts, col, lab in ((dense, GREY, "dense incremental elimination"), (wied, BLUE, "wiedemann (sparse)")):
    if pts:
        pts.sort()
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o", color=col, label=lab)
ax.set_xscale("log", base=2)
ax.set_yscale("log")
ax.set_xlabel("factor-base columns (≈ q/4)")
ax.set_ylabel("row_ops (multiplications mod r)")
ax.set_title("The relation matrix: dense ∝ |F|³, wiedemann ∝ |F|²")
ax.legend(fontsize=8)
save(fig, "la_scaling")

# 3. Large primes: S ratio over the plain arm.
F = S["part_f"]
if F:
    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    arms = ["lp1-700", "lp2-700", "lp2-500"]
    ms = sorted({f["m"] for f in F})
    w = 0.25
    for k, (arm, col) in enumerate(zip(arms, (ORANGE, PURPLE, GREEN))):
        ys = [next((f["s_ratio_over_base"] for f in F if f["arm"] == arm and f["m"] == m), float("nan")) for m in ms]
        ax.bar([i + (k - 1) * w for i in range(len(ms))], ys, w, color=col, label=arm)
    ax.axhline(1.0, color="black", lw=1)
    ax.set_xticks(range(len(ms)))
    ax.set_xticklabels([f"m = {m}" for m in ms])
    ax.set_ylabel("S ÷ S without large primes")
    ax.set_title("One and two large primes against the full base")
    ax.legend(fontsize=8)
    save(fig, "large_primes")
print("figures written")

"""Plot recorded values from raw/summary.json and raw/frozen.json (both
written by the Rust `analyze` / `tune` subcommands).  Computes nothing new.
python3 -I plot.py"""
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
S = json.loads((HERE / "raw" / "summary.json").read_text())
F = json.loads((HERE / "raw" / "frozen.json").read_text())
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

ARMS = [("d4-2lp", "D₄ + 2LP, tuned (α = 0.95)", "#2a6fdb"),
        ("d4-2lp-gaudry", "D₄ + 2LP, Gaudry α = 0.75", "#2b8a3e"),
        ("d4-half", "positive control (D₄, odd trials discarded)", "#999999")]


def save(fig, name):
    fig.tight_layout()
    for ext, kw in (("svg", {}), ("png", {"dpi": 150})):
        fig.savefig(FIG / f"{name}.{ext}", **kw)


# Figure 1: paired whole-cost ratio against the D₄ incumbent, per curve and per size.
fig, axs = plt.subplots(1, 3, figsize=(12.5, 4.0), sharex=True)
for ax, ph, title in zip(axs, ("whole", "relation", "la"),
                         ("whole charged cost", "relation phase", "linear algebra")):
    for arm, label, c in ARMS:
        P = S["paired_vs_d4"][arm][ph]
        for s in P["per_size"]:
            ax.scatter([s["log2_E"]] * len(s["per_curve"]), s["per_curve"], s=10, color=c, alpha=0.45)
        xs = [s["log2_E"] for s in P["per_size"]]
        ys = [s["ratio"] for s in P["per_size"]]
        lo = [s["ratio"] - s["ci95"][0] for s in P["per_size"]]
        hi = [s["ci95"][1] - s["ratio"] for s in P["per_size"]]
        pooled = P["pooled"]
        ax.errorbar(xs, ys, yerr=[lo, hi], color=c, lw=1.2, capsize=2,
                    label=f"{label}: pooled {pooled['ratio']:.3f} [{pooled['ci95'][0]:.3f}, {pooled['ci95'][1]:.3f}]")
    ax.axhline(1.0, color="k", lw=0.8, ls="--")
    ax.set_yscale("log", base=2)
    ax.set_xlabel("log₂ #E(F_{p⁴})")
    ax.set_title(title)
axs[0].set_ylabel("candidate / D₄ (paired, same curve and stream)")
axs[2].set_ylim(bottom=2 ** -3.4)
for ax, loc in zip(axs, ("upper left", "upper left", "lower left")):
    ax.legend(fontsize=6.5, frameon=False, loc=loc)
fig.suptitle("Item 01: D₄ + double large primes vs D₄, held-out test curves (24 curves, p = 101…1301); "
             "dashed = D₄ incumbent; bars = 95 % bootstrap over 3 curves per size", fontsize=9.5)
save(fig, "paired_ratio")

# Figure 2: phase costs against log2 #E with fitted exponents.
COST = [("d4", "D₄ (incumbent)", "#d9480f"), ("d4-2lp", "D₄ + 2LP tuned", "#2a6fdb"),
        ("d4-2lp-gaudry", "D₄ + 2LP Gaudry", "#2b8a3e")]
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), sharex=True)
for ax, ph, key, title in ((axs[0], "relation", "relation_gae", "relation phase (trials × C₄/97 + adds)"),
                           (axs[1], "la", "la_gae", "linear algebra (mod-r mults × 16/97)")):
    for arm, label, c in COST:
        xs = [r["log2_E"] for r in S["rows"]]
        ys = [r["arms"][arm][key] for r in S["rows"]]
        f = S["fits"][arm][ph]
        ax.scatter(xs, ys, s=14, color=c,
                   label=f"{label}: {f['slope']:.3f} [{f['ci95'][0]:.3f}, {f['ci95'][1]:.3f}]")
        x0, x1 = min(xs), max(xs)
        ax.plot([x0, x1], [2 ** (f["intercept_log2"] + f["slope"] * x0), 2 ** (f["intercept_log2"] + f["slope"] * x1)],
                color=c, lw=1)
    ax.set_yscale("log", base=2)
    ax.set_xlabel("log₂ #E")
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=7.5, frameon=False, title="exponent in #E [95 % CI]", title_fontsize=7.5)
axs[0].set_ylabel("GAE")
save(fig, "phases")

# Figure 3: training scores of the 60-config grid.
fig, ax = plt.subplots(figsize=(7.2, 4.2))
styles = {}
for t in F["ranking"]:
    c = t["cfg"]
    key = (c["merge_cap"], c["pick"], c["omega"])
    styles.setdefault(key, []).append((c["alpha"], t["score_gm_ratio_p_ge_min"]))
for (cap, pick, om), pts in sorted(styles.items()):
    pts.sort()
    ax.plot([p[0] for p in pts], [p[1] for p in pts], marker="o", ms=3, lw=0.9,
            ls="-" if pick == "uniform" else ":",
            label=f"cap {cap if cap else '∞'}, {pick}, ω {om}")
ax.axhline(1.0, color="k", lw=0.8, ls="--")
ax.set_yscale("log", base=2)
ax.set_xlabel("large-prime bound α  (|F'| = ⌈columns^α⌉)")
ax.set_ylabel("geo-mean whole-cost ratio vs D₄ (train, p ≥ 521, 6 curves)")
ax.set_title(f"Training grid (14 curves, seeds 4–5; dotted minlp coincides with uniform)\nfrozen: {F['frozen']['name']} = {F['frozen']['score']:.4f}", fontsize=9)
ax.legend(fontsize=6.5, frameon=False, ncol=2)
save(fig, "training_grid")
print("figures written")

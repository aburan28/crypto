# Lays out recorded values only: raw/summary.json (written by
# examples/gaudry_quartic_d4_summary.rs from the ladder rows) and the ecbench
# session's records. No result is computed here.
# Figures: c4_ablation (C4 by configuration, p = 269), c4_ladder (C4 against
# p), crossover (the summary's n* against §11.16's model curves, drawn from
# the summary's inputs), s_by_arm (ecbench S per arm and target).
import json, os, math
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.join(here, "..")
S = json.load(open(os.path.join(root, "raw", "summary.json")))
cfg = {c["name"]: c for c in S["configs"]}
C = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "plum": "#8e5bd1", "ink": "#0b0b0b",
     "ink2": "#52514e", "grid": "#e6e5e1", "surface": "#fcfcfb", "wash": "#9aa3ad"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": C["grid"], "axes.labelcolor": C["ink2"],
                     "xtick.color": C["ink2"], "ytick.color": C["ink2"], "axes.spines.top": False,
                     "axes.spines.right": False,
                     "figure.facecolor": C["surface"], "axes.facecolor": C["surface"]})


def save(fig, name):
    fig.tight_layout()
    for ext, kw in (("svg", {}), ("png", {"dpi": 150})):
        fig.savefig(os.path.join(here, f"{name}.{ext}"), **kw)
    plt.close(fig)


# 1. Ablation at p = 269.
ABL = [("abl_committed", "committed\n(plain, deg 14)"), ("abl_s", "+ sparsest"), ("abl_sc", "+ sparsest\n+ chains"),
       ("abl_scr", "+ sparsest\n+ chains\n+ replay"), ("abl_w", "weighted\n(deg 26)"), ("abl_ws", "weighted\n+ sparsest"),
       ("abl_wc", "weighted\n+ chains"), ("abl_wcs", "weighted\n+ chains\n+ sparsest"),
       ("abl_wcsr", "weighted\n+ chains\n+ sparsest\n+ replay")]
ABL = [a for a in ABL if a[0] in cfg]
PH = [("echelon_muls", "elimination", C["blue"]), ("nf_muls", "normal forms", C["orange"]),
      ("charpoly_muls", "char. polynomial / det A", C["aqua"]), ("eigenvector_muls", "eigenvectors", C["plum"])]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.6), gridspec_kw={"width_ratios": [1.3, 1]})
for i, (k, lab) in enumerate(ABL):
    c = cfg[k]
    m, (lo, hi) = c["mean_muls"], c["ci95_muls"]
    col = C["wash"] if i == 0 else (C["blue"] if k.startswith("abl_w") else C["plum"])
    a1.bar(i, m, color=col, width=0.62)
    a1.errorbar(i, m, yerr=[[m - lo], [hi - m]], color=C["ink"], lw=0.8, capsize=2)
    a1.text(i, m * 1.25, f"{m:.2g}", ha="center", fontsize=7.5, color=C["ink"])
a1.set_yscale("log")
a1.set_ylim(1e7, 1e11)
a1.set_xticks(range(len(ABL)), [a[1] for a in ABL], fontsize=6.5)
a1.set_ylabel("F_p multiplications per decomposition (C₄)")
a1.grid(axis="y", color=C["grid"], lw=0.6)
a1.set_title("C₄ by option, p = 269, 10 planted targets (mean, 95 % bootstrap)", fontsize=9)
bottom = [0.0] * len(ABL)
for key, lab, col in PH:
    fr = [cfg[k]["phase_means"][key] / sum(cfg[k]["phase_means"].values()) for k, _ in ABL]
    a2.bar(range(len(ABL)), fr, bottom=bottom, color=col, width=0.62, label=lab, edgecolor=C["surface"], lw=1)
    bottom = [b + f for b, f in zip(bottom, fr)]
a2.set_xticks(range(len(ABL)), [str(i + 1) for i in range(len(ABL))], fontsize=8)
a2.set_xlabel("configuration (left panel, in order)")
a2.set_ylabel("share of C₄")
a2.set_ylim(0, 1)
a2.legend(frameon=False, fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2)
a2.set_title("where the multiplications go", fontsize=9)
save(fig, "c4_ablation")

# 2. Ladder.
fig, ax = plt.subplots(figsize=(6.4, 4))
for kind, col, mk, lab in [("new_planted", C["blue"], "o", "new, planted"), ("new_uniform", C["aqua"], "s", "new, uniform"),
                           ("committed_planted", C["wash"], "D", "committed, planted"),
                           ("committed_uniform", C["ink2"], "v", "committed, uniform")]:
    pts = [(cfg[f"p{p}_{kind}"]["p"], cfg[f"p{p}_{kind}"]) for p in (29, 269, 1009, 4049) if f"p{p}_{kind}" in cfg]
    if not pts:
        continue
    xs = [p for p, _ in pts]
    ys = [c["mean_muls"] for _, c in pts]
    lo = [c["mean_muls"] - c["ci95_muls"][0] for _, c in pts]
    hi = [c["ci95_muls"][1] - c["mean_muls"] for _, c in pts]
    ax.errorbar(xs, ys, yerr=[lo, hi], color=col, marker=mk, ms=6, lw=1.6, capsize=2, label=lab)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylim(1e7, 1e11)
ax.set_xlabel("p (curves over F_p⁴)")
ax.set_ylabel("C₄ (F_p multiplications, mean, 95 % bootstrap)")
ax.grid(color=C["grid"], lw=0.6)
ax.legend(frameon=False, fontsize=8)
ax.set_title("one decomposition against p: the committed and the new oracle", fontsize=9)
save(fig, "c4_ladder")

# 3. Crossover: the model curves with the summary's inputs.
X = [x for x in S["crossover"] if "model" in x]
if X:
    inp = X[0]
    coef = 181 / (4 * 126)
    l1, r1, l2, r2 = math.log2(1009), inp["roots_muls_at_1009"], math.log2(4049), inp["roots_muls_at_4049"]
    base = inp["c4_new_at_4049"] - r2
    cnew = lambda p: base + r1 + (r2 - r1) / (l2 - l1) * (math.log2(p) - l1)
    lg = [80 + 0.5 * k for k in range(161)]
    fig, ax = plt.subplots(figsize=(6.4, 4))
    for lab, c4, col in [("committed C₄ = 2.54·10¹⁰", lambda p: 2.54e10, C["wash"]), ("new C₄ (measured, p = 4049)", cnew, C["blue"])]:
        for h, ls in [(1, "-"), (2, "--")]:
            ys = [math.sqrt(h) * (coef * c4(2 ** (x / 4)) / 2 ** (x / 4) + 0.34 / 4) for x in lg]
            ax.plot(lg, ys, color=col, ls=ls, lw=1.6, label=f"{lab}, h = {h}")
    for x in S["crossover"]:
        if x.get("r_inf") == "phi=0.73" and x.get("h") == 2:
            for key, col in [("log2_n_star_committed_c4_2.54e10", C["wash"]), ("log2_n_star_new", C["blue"])]:
                ax.axvline(x[key], color=col, lw=0.8, ls=":")
                ax.text(x[key] + 0.6, 0.065, f"n* ≈ 2^{x[key]:.1f} (h = 2)", color=C["ink"], fontsize=7.5)
    ax.axhline(1, color=C["ink2"], lw=0.8)
    ax.text(81, 1.2, "rho", fontsize=8, color=C["ink2"])
    ax.set_yscale("log")
    ax.set_ylim(0.05, 1e4)
    ax.set_xlabel("log₂ group order n")
    ax.set_ylabel("index calculus ÷ rho (model, r∞ at φ = 0.73)")
    ax.grid(color=C["grid"], lw=0.6)
    ax.legend(frameon=False, fontsize=7.5)
    ax.set_title("§11.16 crossover model (extrapolated, not measured)", fontsize=9)
    save(fig, "crossover")

# 4. ecbench S per arm and target.
rec = os.path.join(root, "sessions", "k4-p29-d4cost", "records.jsonl")
if os.path.exists(rec):
    runs = [json.loads(l) for l in open(rec)]
    arms = ["rho-neg", "ic-d4-committed", "ic-d4-chains", "rho-neg-aa"]
    wl = sorted({r["workload"]["workload_id"] for r in runs}, key=lambda w: [r["workload"]["target_index"] for r in runs if r["workload"]["workload_id"] == w][0])
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    w = 0.38
    for j, wid in enumerate(wl):
        for i, a in enumerate(arms):
            rs = [r for r in runs if r["arm"] == a and r["workload"]["workload_id"] == wid and r["outcome"]["status"] == "verified"]
            if not rs:
                continue
            s = rs[0]["cost"]["s"]
            ax.bar(i + (j - 0.5) * w, s, width=w * 0.92, color=[C["blue"], C["orange"]][j],
                   label=f"target {j} ({wid})" if i == 0 else None)
            ax.text(i + (j - 0.5) * w, s * 1.4, f"{s:.3g}", ha="center", fontsize=7, color=C["ink"])
    ax.set_yscale("log")
    ax.set_ylim(0.3, 1e10)
    ax.set_xticks(range(len(arms)), arms, fontsize=8)
    ax.set_ylabel("S (group-addition equivalents / √r)")
    ax.grid(axis="y", color=C["grid"], lw=0.6)
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("ecbench session k4-p29-d4cost, E(F_{29⁴}), one round", fontsize=9)
    save(fig, "s_by_arm")

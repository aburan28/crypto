# Lays out recorded values only: raw/summary.json, written by
# examples/gaudry_quartic_bin_summary.rs from the ladder rows and the
# ecbench sessions. No result is computed here; the trend line and the
# model curve are drawn from the summary's fitted coefficients and inputs.
import json, os, math
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
S = json.load(open(os.path.join(here, "..", "raw", "summary.json")))
cfg = {c["name"]: c for c in S["configs"]}
C = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "plum": "#8e5bd1", "ink": "#0b0b0b",
     "ink2": "#52514e", "grid": "#e6e5e1", "surface": "#fcfcfb", "wash": "#9aa3ad"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": C["grid"], "axes.labelcolor": C["ink2"],
                     "xtick.color": C["ink2"], "ytick.color": C["ink2"], "axes.spines.top": False,
                     "axes.spines.right": False, "figure.facecolor": C["surface"], "axes.facecolor": C["surface"]})


def save(fig, name):
    fig.tight_layout()
    for ext, kw in (("svg", {}), ("png", {"dpi": 150})):
        fig.savefig(os.path.join(here, f"{name}.{ext}"), **kw)
    plt.close(fig)


# 1. Ablation at m = 8.
ABL = [("abl_plain_s", "plain\n+ sparsest", C["wash"]), ("abl_plain_sr", "plain\n+ s + replay", C["wash"]),
       ("abl_weighted", "weighted\n(2,2,2,1)", C["plum"]), ("abl_weighted_s", "weighted\n+ sparsest", C["plum"]),
       ("abl_weighted_sr", "weighted\n+ s + replay", C["plum"]), ("abl_frob", "Frobenius\n(2,2,1,1)", C["blue"]),
       ("abl_frob_s", "Frobenius\n+ sparsest", C["blue"]), ("abl_frob_sr", "Frobenius\n+ s + replay", C["blue"])]
ABL = [a for a in ABL if a[0] in cfg]
PH = [("echelon_muls", "elimination", C["blue"]), ("nf_muls", "normal forms", C["orange"]),
      ("charpoly_muls", "det A", C["aqua"]), ("roots_muls", "roots", C["plum"]),
      ("split_muls", "split and lifts", C["wash"])]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.6), gridspec_kw={"width_ratios": [1.3, 1]})
for i, (k, lab, col) in enumerate(ABL):
    c = cfg[k]
    m, (lo, hi) = c["mean_muls"], c["ci95_muls"]
    a1.bar(i, m, color=col, width=0.62)
    a1.errorbar(i, m, yerr=[[m - lo], [hi - m]], color=C["ink"], lw=0.8, capsize=2)
    a1.text(i, m * 1.25, f"{m:.2g}", ha="center", fontsize=7.5, color=C["ink"])
a1.set_yscale("log")
a1.set_ylim(1e4, 1e9)
a1.set_xticks(range(len(ABL)), [a[1] for a in ABL], fontsize=6.5)
a1.set_ylabel("F_q multiplications per decomposition (C₄)")
a1.grid(axis="y", color=C["grid"], lw=0.6)
a1.set_title("C₄ by form, m = 8 (q = 256), 10 planted targets (mean, 95 % bootstrap)", fontsize=9)
bottom = [0.0] * len(ABL)
for key, lab, col in PH:
    fr = [cfg[k]["phase_means"][key] / sum(cfg[k]["phase_means"].values()) for k, _, _ in ABL]
    a2.bar(range(len(ABL)), fr, bottom=bottom, color=col, width=0.62, label=lab, edgecolor=C["surface"], lw=1)
    bottom = [b + f for b, f in zip(bottom, fr)]
a2.set_xticks(range(len(ABL)), [str(i + 1) for i in range(len(ABL))], fontsize=8)
a2.set_xlabel("configuration (left panel, in order)")
a2.set_ylabel("share of C₄")
a2.set_ylim(0, 1)
a2.legend(frameon=False, fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3)
a2.set_title("where the multiplications go", fontsize=9)
save(fig, "c4_ablation")

# 2. C4 against m.
fig, ax = plt.subplots(figsize=(6.4, 4))
for kind, col, mk, lab in [("planted", C["blue"], "o", "planted (exact against MITM)"),
                           ("uniform", C["aqua"], "s", "uniform (exact against MITM)")]:
    pts = [(c["m"], c) for c in S["configs"] if c["name"] == f"m{c['m']}_{kind}"]
    if pts:
        ax.errorbar([p for p, _ in pts], [c["mean_muls"] for _, c in pts],
                    yerr=[[c["mean_muls"] - c["ci95_muls"][0] for _, c in pts], [c["ci95_muls"][1] - c["mean_muls"] for _, c in pts]],
                    color=col, marker=mk, ms=6, lw=1.6, capsize=2, label=lab)
if "m15_uniform_cost_only" in cfg:
    c = cfg["m15_uniform_cost_only"]
    ax.plot([15], [c["mean_muls"]], marker="D", ms=7, color=C["plum"], ls="none", label="uniform, cost only (no MITM)")
ax.set_yscale("log")
ax.set_ylim(5e4, 1e6)
ax.set_xlabel("m  (q = 2^m, curves over F_(q⁴))")
ax.set_ylabel("C₄ (F_q multiplications, Frobenius form)")
ax.grid(color=C["grid"], lw=0.6)
ax.legend(frameon=False, fontsize=8)
ax.set_title("one decomposition against m (mean, 95 % bootstrap)", fontsize=9)
save(fig, "c4_ladder")

# 3. The measured trend and the model.
T = S.get("trend", {})
ic = [p for p in S["sessions"]["ic"] if p["arm"] == "ic-frobenius"]
if ic:
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    xs = [p["log2_n"] for p in ic if p["status"] == "verified"]
    ys = [p["ratio_to_rho_theory"] for p in ic if p["status"] == "verified"]
    ax.plot(xs, ys, "o", color=C["blue"], ms=6, label="ecbench, IC (Frobenius) ÷ rho's theory S (measured)")
    if "slope" in T:
        a, b = T["intercept"], T["slope"]
        lg = [20 + 0.5 * k for k in range(121)]
        ax.plot(lg, [2 ** (a + b * x) for x in lg], color=C["blue"], lw=1.2, ls="--",
                label=f"fit, m ≥ 8: slope {b:.3f}, n* ≈ 2^{T['log2_n_star']:.1f} (extrapolated)")
    M = S.get("model", {})
    if "crossover" in M:
        spr, cadd = M["spr_trials_per_relation"], M["c_add_fq_muls"]
        r12, r15, c15 = M["roots_at_m12"], M["roots_at_m15"], M["c4_at_m15"]
        base = c15 - r15
        lg = [20 + 0.5 * k for k in range(121)]
        def model(x, h=2.0):
            mm = x / 4
            q = 2 ** mm
            c4 = base + r15 + (r15 - r12) / 3 * (mm - 15)
            return math.sqrt(h) * (spr * (q / 4) * c4 / (1.3 * q * q * cadd) + 0.34 / 4)
        ax.plot(lg, [model(x) for x in lg], color=C["orange"], lw=1.4,
                label="§11.16 model, measured inputs, h = 2 (r∞ assumed)")
    ax.axhline(1, color=C["ink2"], lw=0.8)
    ax.text(21, 1.25, "rho", fontsize=8, color=C["ink2"])
    ax.set_yscale("log")
    ax.set_ylim(0.05, 1e4)
    ax.set_xlim(20, 80)
    ax.set_xlabel("log₂ group order n")
    ax.set_ylabel("index calculus ÷ rho")
    ax.grid(color=C["grid"], lw=0.6)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")
    ax.set_title("k = 4 in characteristic 2: measured points, fitted trend, model", fontsize=9)
    save(fig, "crossover")

# 4. The m = 6 forms session.
f6 = [p for p in S["sessions"]["ic"] if p["session"] == "char2-forms-m6"]
if f6:
    arms = ["ic-plain", "ic-weighted", "ic-frobenius"]
    wl = sorted({p["workload"] for p in f6})
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    w = 0.38
    for j, wid in enumerate(wl):
        for i, a in enumerate(arms):
            ps = [p for p in f6 if p["arm"] == a and p["workload"] == wid and p["status"] == "verified"]
            if not ps:
                continue
            s = ps[0]["s"]
            ax.bar(i + (j - 0.5) * w, s, width=w * 0.92, color=[C["blue"], C["orange"]][j],
                   label=f"target {j} ({wid})" if i == 0 else None)
            ax.text(i + (j - 0.5) * w, s * 1.3, f"{s:.3g}", ha="center", fontsize=7, color=C["ink"])
    rhos = [r["s"] for r in S["sessions"]["rho"] if r["session"] == "char2-forms-m6"]
    if rhos:
        ax.axhline(sum(rhos) / len(rhos), color=C["ink2"], lw=0.8)
        ax.text(2.3, sum(rhos) / len(rhos) * 1.3, "rho (mean of 4 runs)", fontsize=7.5, color=C["ink2"])
    ax.set_yscale("log")
    ax.set_xticks(range(len(arms)), ["plain (deg 14)", "weighted (deg 26)", "Frobenius (deg 11)"], fontsize=8)
    ax.set_ylabel("S (group-addition equivalents / √r)")
    ax.grid(axis="y", color=C["grid"], lw=0.6)
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("ecbench session char2-forms-m6, E(F_(2^6)⁴), one round", fontsize=9)
    save(fig, "s_forms_m6")

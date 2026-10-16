"""Lay out the recorded hwcost round-2 results as graphs (SVG + 150 dpi PNG).

Plotting only: every value drawn is read from the files `hwcost measure`,
`fit`, `compare`, `gpu-check` and `ic-price` wrote; nothing is computed here
beyond unit conversions for display (AGENTS.md "Implementation language").

    python plot.py RUNS_DIR OUT_DIR
"""

import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

RUNS = Path(sys.argv[1]).resolve()
OUT = Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 9, "svg.fonttype": "none"})
FIELDS = ["fp64", "fp64bf", "w1", "w2", "w3", "w4"]
ISA = ["isa-base-1", "isa-native-1", "isa-base-2", "isa-native-2"]
ISA_STYLE = {
    "isa-base-1": ("#1b6ca8", "o", "-"),
    "isa-base-2": ("#1b6ca8", "s", "--"),
    "isa-native-1": ("#d1495b", "o", "-"),
    "isa-native-2": ("#d1495b", "s", "--"),
}


def save(fig, name):
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def read_csv(p):
    with open(p) as f:
        return list(csv.DictReader(f))


def started(run):
    p = RUNS / run / "started_utc.txt"
    return p.read_text().strip() if p.exists() else "?"


# A. Kernel time per field, baseline vs native, all four runs.
val = {r: read_csv(RUNS / r / "validation.csv") for r in ISA if (RUNS / r / "validation.csv").exists()}
if val:
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.6))
    axes = axes.flatten()
    for ax, f in zip(axes, FIELDS):
        for r, rows in val.items():
            rs = sorted((x for x in rows if x["field"] == f), key=lambda x: int(x["batch"]))
            c, mk, ls = ISA_STYLE[r]
            ax.plot([int(x["batch"]) for x in rs], [float(x["measured_ns"]) for x in rs], ls, marker=mk, ms=3, color=c, label=r)
        ax.set_xscale("log", base=2)
        ax.set_ylim(bottom=0)
        ax.set_title(f, fontsize=9)
        ax.set_xlabel("walks per batch B")
        ax.grid(alpha=0.3)
    for ax in axes[::3]:
        ax.set_ylabel("ns per rho iteration (median of 7)")
    axes[0].legend(fontsize=7)
    fig.suptitle("Rho kernel, baseline (blue) vs -C target-cpu=native (red) builds, runs A1 B1 A2 B2:\n" + ", ".join(f"{r} {started(r)}" for r in val), fontsize=8)
    fig.tight_layout()
    save(fig, "fig_isa_kernel")

# A. Ratio native/baseline per arm (r1, r2) with the drift margin.
cmp_path = RUNS / "compare" / "compare.csv"
if cmp_path.exists():
    rows = read_csv(cmp_path)
    groups = [("field ops", [r for r in rows if r["arm"].split(".")[0] in FIELDS]), ("kernel", [r for r in rows if r["arm"].startswith("kernel.")])]
    fig, axes = plt.subplots(2, 1, figsize=(14, 6.5))
    colors = {"candidate faster": "#2e8540", "candidate slower": "#d1495b", "no detected effect": "#888888"}
    for ax, (title, rs) in zip(axes, groups):
        for i, r in enumerate(rs):
            m = float(r["margin"])
            ax.plot([i, i], [1 - m, 1 + m], color="#dddddd", lw=4, zorder=0)
            ax.plot(i, float(r["r1"]), "o", color=colors[r["class"]], ms=4)
            ax.plot(i, float(r["r2"]), "s", mfc="none", color=colors[r["class"]], ms=4)
        ax.axhline(1, color="black", lw=0.6)
        ax.set_xticks(range(len(rs)))
        ax.set_xticklabels([r["arm"].replace("kernel.", "") for r in rs], rotation=90, fontsize=5.5)
        ax.set_ylabel("native / baseline")
        ax.set_title(f"{title}: r1 (filled) and r2 (open); grey = ±max(2·drift, 2%); green faster, red slower, grey no detected effect", fontsize=8)
        ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    save(fig, "fig_isa_ratios")

# B. fp64 vs fp64bf kernel in every run.
pairs_path = RUNS / "compare" / "fp64_pairs.csv"
if pairs_path.exists():
    pr = read_csv(pairs_path)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 3.6))
    for r in ISA:
        rs = sorted((x for x in pr if x["run"] == r), key=lambda x: int(x["batch"]))
        if not rs:
            continue
        c, mk, ls = ISA_STYLE[r]
        b = [int(x["batch"]) for x in rs]
        a1.plot(b, [float(x["fp64_ns"]) for x in rs], ls, marker=mk, ms=3, color=c, label=f"fp64, {r}")
        a1.plot(b, [float(x["fp64bf_ns"]) for x in rs], ls, marker=mk, ms=3, color=c, alpha=0.45, label=f"fp64bf, {r}")
        a2.plot(b, [float(x["ratio"]) for x in rs], ls, marker=mk, ms=3, color=c, label=r)
    a2.axhline(0.85, color="black", lw=0.8, ls=":", label="F1 threshold 0.85")
    a2.axvspan(4096, 2**22 * 1.3, color="#eeeeee", zorder=0)
    for ax in (a1, a2):
        ax.set_xscale("log", base=2)
        ax.set_xlabel("walks per batch B")
        ax.grid(alpha=0.3)
    a1.set_ylabel("ns per iteration")
    a1.set_ylim(bottom=0)
    a2.set_ylabel("fp64bf / fp64")
    a1.legend(fontsize=6, ncol=2)
    a2.legend(fontsize=7)
    a2.set_title("grey: B ≥ 4096 (F1 cells)", fontsize=8)
    fig.suptitle("Same field, same values: Fp64 vs branch-free Fp64Bf in the rho kernel", fontsize=9)
    save(fig, "fig_fp64_branches")

# B, addendum 1: fp64 vs fp64cm (select_unpredictable) in cm-1 and cm-2, with
# the F3 instruction counts of the three kernels and the a priori model.
cm_path = RUNS / "cm-pair" / "pairs.csv"
if cm_path.exists():
    pr = read_csv(cm_path)
    fig, ((a1, a2), (a3, a4)) = plt.subplots(2, 2, figsize=(12, 8))
    CM_STYLE = {"cm-1": ("o", "-"), "cm-2": ("s", "--")}
    for r, (mk, ls) in CM_STYLE.items():
        rs = sorted((x for x in pr if x["run"] == r), key=lambda x: int(x["batch"]))
        b = [int(x["batch"]) for x in rs]
        a1.plot(b, [float(x["fp64_ns"]) for x in rs], ls, marker=mk, ms=3, color="#1b6ca8", label=f"fp64, {r}")
        a1.plot(b, [float(x["fp64cm_ns"]) for x in rs], ls, marker=mk, ms=3, color="#2e8540", label=f"fp64cm, {r}")
        a2.plot(b, [float(x["ratio"]) for x in rs], ls, marker=mk, ms=3, color="#2e8540", label=f"fp64cm / fp64, {r}")
    if pairs_path.exists():
        for r in ["isa-base-1", "isa-base-2"]:
            rs = sorted((x for x in read_csv(pairs_path) if x["run"] == r), key=lambda x: int(x["batch"]))
            a2.plot([int(x["batch"]) for x in rs], [float(x["ratio"]) for x in rs], ":", color="#888888", lw=1, label=f"fp64bf / fp64, {r}")
    a2.axhline(0.85, color="black", lw=0.8, ls=":", label="F1′ threshold 0.85")
    a2.axvspan(4096, 2**22 * 1.3, color="#eeeeee", zorder=0)
    for ax in (a1, a2):
        ax.set_xscale("log", base=2)
        ax.set_xlabel("walks per batch B")
        ax.grid(alpha=0.3)
    a1.set_ylabel("ns per rho iteration (median of 7)")
    a1.set_ylim(bottom=0)
    a1.legend(fontsize=6)
    a1.set_title("kernel time, baseline build", fontsize=8)
    a2.set_ylabel("kernel ratio")
    a2.set_ylim(0, 1.2)
    a2.legend(fontsize=6)
    a2.set_title("grey band: B ≥ 4096 (F1′ cells)", fontsize=8)
    # a priori (V2) predictions against measured, both cm runs.
    for r, (mk, ls) in CM_STYLE.items():
        v = read_csv(RUNS / r / "validation.csv")
        for f, c in [("fp64", "#1b6ca8"), ("fp64cm", "#2e8540")]:
            rs = sorted((x for x in v if x["field"] == f), key=lambda x: int(x["batch"]))
            b = [int(x["batch"]) for x in rs]
            a3.plot(b, [float(x["measured_ns"]) for x in rs], ls, marker=mk, ms=3, color=c, label=f"{f} measured, {r}")
            a3.plot(b, [float(x["a_priori_ns"]) for x in rs], ls, color=c, alpha=0.4, lw=2.5, label=f"{f} a priori, {r}")
    a3.set_xscale("log", base=2)
    a3.set_ylim(bottom=0)
    a3.set_xlabel("walks per batch B")
    a3.set_ylabel("ns per rho iteration")
    a3.grid(alpha=0.3)
    a3.legend(fontsize=5.5, ncol=2)
    a3.set_title("measured vs a priori model (no fitted parameter)", fontsize=8)
    f3 = RUNS / "f3" / "counts.csv"
    if f3.exists():
        cnt = read_csv(f3)
        labels = [f"{x['label']}\n{x['function'].replace('step_', '')}" for x in cnt]
        xs = range(len(cnt))
        w = 0.2
        for k, (col, c) in enumerate([("jcc", "#d1495b"), ("cmov", "#2e8540"), ("mulx", "#e1a100"), ("vector", "#7a5195")]):
            a4.bar([i + (k - 1.5) * w for i in xs], [int(x[col]) for x in cnt], w, color=c, label=col)
        for i, x in enumerate(cnt):
            a4.text(i, max(int(x[k]) for k in ("jcc", "cmov", "mulx", "vector")) + 1, f"{x['instructions']} ins", ha="center", fontsize=6)
        a4.set_ylim(0, 1.15 * max(int(x[k]) for x in cnt for k in ("jcc", "cmov", "mulx", "vector")) + 3)
        a4.set_xticks(list(xs))
        a4.set_xticklabels(labels, fontsize=6.5)
        a4.set_ylabel("instructions in the kernel function")
        a4.legend(fontsize=6)
        a4.set_title("F3: objdump -d of step_* (runs/f3)", fontsize=8)
        a4.grid(axis="y", alpha=0.3)
    fig.suptitle("Addendum 1: Fp64 vs Fp64Cm (select_unpredictable), same field and values\nruns cm-1 " + started("cm-1") + ", cm-2 " + started("cm-2") + "; host Xeon family 6 model 207, L0", fontsize=9)
    fig.tight_layout()
    save(fig, "fig_fp64_cm")

# C. GPU: predicted vs measured rate.
gpu_path = RUNS / "gpu" / "gpu_check.csv"
if gpu_path.exists():
    g = read_csv(gpu_path)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
    panels = [
        ("A priori (no fitted parameter)", lambda r: r["set"] in ("fit", "same-gpu"), "a_priori_bps"),
        ("η from each GPU's base build", lambda r: r["set"] in ("fit", "same-gpu"), "calibrated_bps"),
        ("η from the RTX PRO 6000 base build", lambda r: r["set"] == "cross-gpu" or (r["gpu"] == "rtx-pro-6000" and r["set"] in ("fit", "same-gpu")), "calibrated_bps"),
    ]
    style = {("rtx-pro-6000", "separate"): ("#1b6ca8", "o"), ("rtx-pro-6000", "shared-int"): ("#1b6ca8", "^"), ("b200", "separate"): ("#d1495b", "o"), ("b200", "shared-int"): ("#d1495b", "^")}
    for ax, (title, sel, key) in zip(axes, panels):
        for (gpu, model), (c, mk) in style.items():
            rs = [r for r in g if r["gpu"] == gpu and r["model"] == model and sel(r)]
            ax.plot([float(r["measured_bps"]) for r in rs], [float(r[key]) for r in rs], mk, color=c, mfc="none" if model == "shared-int" else c, ms=5, label=f"{gpu}, {model} pipes")
        lim = [12, 26]
        ax.plot(lim, lim, color="black", lw=0.7)
        for k, ls in [(0.1, ":"), (0.3, "--")]:
            ax.plot(lim, [x * (1 + k) for x in lim], color="grey", lw=0.6, ls=ls)
            ax.plot(lim, [x * (1 - k) for x in lim], color="grey", lw=0.6, ls=ls)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_xlabel("measured, B updates/s")
        ax.set_title(title, fontsize=9)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("predicted, B updates/s")
    axes[0].legend(fontsize=7)
    fig.suptitle("ECC2K-130 GPU builds (ecc2k130/benchmarks/roofline/calibration.json, 46 builds): dotted ±10%, dashed ±30%", fontsize=9)
    save(fig, "fig_gpu")

# D. IC vs rho in wall time: measured, band-predicted, and projections.
for dom in ["koblitz", "prime"]:
    d = RUNS / f"ic-{dom}"
    if not (d / "pairs.csv").exists():
        continue
    p = read_csv(d / "pairs.csv")
    proj = read_csv(d / "projections.csv")
    fig, ((a1, a2), (a3, a4)) = plt.subplots(2, 2, figsize=(12, 9))
    x = [float(r["log2_r"]) for r in p]
    a1.plot(x, [float(r["ic_wall_ms"]) for r in p], "o-", color="black", label=f"{p[0]['ic_arm']} measured")
    a1.plot(x, [float(r["rho_wall_ms"]) for r in p], "s-", color="#1b6ca8", label=f"{p[0]['rho_arm']} measured")
    for k, (lab, c) in enumerate([("w0", "#2e8540"), ("w01", "#e1a100"), ("w1", "#d1495b")]):
        a1.plot(x, [float(r[f"ic_pred_ms_{lab}"]) for r in p], "--", color=c, label=f"IC predicted, w={['0', '0.1', '1'][k]}")
    a1.set_yscale("log")
    a1.set_xlabel("log2 r")
    a1.set_ylabel("ms per solve (mean of 32)")
    a1.legend(fontsize=7)
    a1.set_title("measured vs band-predicted wall time", fontsize=9)
    a2.plot(x, [float(r["w_star"]) for r in p], "o-", color="black")
    a2.axhspan(0, 1, color="#e8f4e8", zorder=0, label="§5a band [0, 1]")
    a2.set_xlabel("log2 r")
    a2.set_ylabel("implied weight w* (GAE per unpriced op)")
    a2.legend(fontsize=7)
    a2.set_title("weight that wall time charged", fontsize=9)
    xp = [float(r["log2_r"]) for r in proj]
    a3.plot(xp, [float(r["ic_core_years"]) for r in proj], "o-", color="black", label="IC (extrapolated)")
    a3.plot(xp, [float(r["rho_core_years"]) for r in proj], "s-", color="#1b6ca8", label="rho (extrapolated)")
    a3.set_yscale("log")
    a3.set_xlabel("log2 r")
    a3.set_ylabel("core-years per solve")
    a3.legend(fontsize=7)
    a3.set_title("pre-registered: fitted wall-time exponents (not measured;\nrho's ~3 ms per-solve floor flattens its exponent)", fontsize=8)
    gae_path = d / "projections_gae.csv"
    if gae_path.exists():
        g = read_csv(gae_path)
        for w, c in [("0", "#2e8540"), ("0.1", "#e1a100"), ("1", "#d1495b")]:
            rs = [r for r in g if r["w"] == w]
            xg = [float(r["log2_r"]) for r in rs]
            a4.plot(xg, [float(r["ic_core_years"]) for r in rs], "o-", color=c, ms=3, label=f"IC, w={w}")
            a4.plot(xg, [float(r["rho_core_years"]) for r in rs], "s--", color=c, ms=3, mfc="none", label=f"rho, w={w}")
        a4.set_yscale("log")
        a4.set_xlabel("log2 r")
        a4.set_ylabel("core-years per solve")
        a4.legend(fontsize=6, ncol=2)
        a4.set_title(f"post hoc: fitted P(w) exponents × {g[0]['ns_per_gae_arm']}'s ns per GAE\nat log2 r {g[0]['ns_per_gae_log2_r']} (not measured)", fontsize=8)
    for ax in (a1, a2, a3, a4):
        ax.grid(alpha=0.3)
    fig.suptitle(f"{dom}: IC vs matched rho in wall time, research/ic_frontier_campaign_20261015 Phase 1\n(32 verified runs per arm and curve; host Xeon 2.30 GHz, L0)", fontsize=9)
    fig.tight_layout()
    save(fig, f"fig_ic_{dom}")

print("wrote", sorted(q.name for q in OUT.glob("fig_*")))

"""Plot recorded values only (AGENTS.md: Python may lay out recorded results, never compute them).

Reads the JSON records the native Rust examples wrote into ../ and draws SVGs into this
directory. Ratios shown are the ratios of recorded medians; nothing is re-measured.
"""

import glob
import json
import os

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
plt.rcParams.update({"font.size": 8.5, "svg.fonttype": "none", "font.family": "DejaVu Sans"})

INK, MUTED, GRID = "#0b0b0b", "#52514e", "#c3c2b7"
C = {"wnaf": "#8a8a85", "gls2": "#2a6fb0", "affine": "#d08a2c", "proj": "#b4442c", "aa": "#c3c2b7"}


def load(name):
    with open(os.path.join(STUDY, name)) as fh:
        return json.load(fh)


def arm(rec, name):
    return next(a for a in rec["arms"] if a["name"] == name)


# ---- 1. D = -619: per-round ns per arm, every benchmark run on D = -619
def fig_bench_arms():
    runs = sorted(glob.glob(os.path.join(STUDY, "bench_run*.json")))
    names = ["wnaf", "wnaf (A/A)", "gls2", "glvgls4 (affine loop)", "glvgls4 (projective loop)"]
    cols = [C["wnaf"], C["aa"], C["gls2"], C["affine"], C["proj"]]
    fig, axes = plt.subplots(1, len(runs), figsize=(3.6 * len(runs), 2.9), squeeze=False, sharey=True)
    for ax, path in zip(axes[0], runs):
        rec = json.load(open(path))
        data = [[x / 1000 for x in arm(rec, n)["rounds_ns"]] for n in names]
        bp = ax.boxplot(data, widths=0.55, patch_artist=True, showfliers=True, flierprops={"markersize": 2.5})
        for patch, c in zip(bp["boxes"], cols):
            patch.set_facecolor(c)
            patch.set_alpha(0.75)
        for med in bp["medians"]:
            med.set_color(INK)
        ax.set_xticks(range(1, len(names) + 1), ["wNAF", "wNAF\n(A/A)", "GLS-2", "GLV-GLS-4\naffine\nloop", "GLV-GLS-4\nproj.\nloop"], fontsize=7)
        ax.grid(axis="y", color=GRID, lw=0.5)
        ax.set_title(f"{os.path.basename(path)[:-5]}: D = {rec['d']}, {rec['rounds']} rounds x {rec['pairs']} scalars", fontsize=8)
    axes[0][0].set_ylabel("microseconds per scalar multiplication\n(per-round mean over the scalars)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "bench_arms.svg"))
    plt.close(fig)


# ---- 2. measured GLS-2 / GLV-GLS-4 ratio across D, against round 2's derived model
# Round-2 model ratios are the recorded values in research/isogeny_loop_glv_20261011/README.md section 3
# (embedding s; derived operation counts, not measured).
MODEL_R2 = {-23: 1.316, -71: 1.292, -431: 1.248, -619: 1.220, -3019: 1.158}


def fig_bench_trend():
    recs = []
    for path in sorted(glob.glob(os.path.join(STUDY, "bench_*.json"))):
        rec = json.load(open(path))
        recs.append((os.path.basename(path)[:-5], rec))
    fig, ax = plt.subplots(figsize=(6.2, 3.1))
    for tag, rec in recs:
        g = arm(rec, "gls2")["median_ns"]
        cost = rec["loop"]["cost_M"]
        rp = g / arm(rec, "glvgls4 (projective loop)")["median_ns"]
        ra = g / arm(rec, "glvgls4 (affine loop)")["median_ns"]
        ax.scatter([cost], [rp], color=C["proj"], s=22, zorder=3)
        ax.scatter([cost], [ra], color=C["affine"], s=18, marker="s", zorder=3)
        ax.annotate(f"D={rec['d']}\n{tag}", (cost, rp), textcoords="offset points", xytext=(5, 4), fontsize=6.5, color=MUTED)
    xs = sorted(MODEL_R2)
    costs = {-23: 37.5, -71: 60, -431: 97.5, -619: 127.5, -3019: 210}
    ax.plot([costs[d] for d in xs], [MODEL_R2[d] for d in xs], "--", color=INK, lw=0.8, label="round-2 derived op-count model")
    ax.scatter([], [], color=C["proj"], s=22, label="measured, projective loop (median / median)")
    ax.scatter([], [], color=C["affine"], s=18, marker="s", label="measured, affine loop (median / median)")
    ax.axhline(1.0, color=GRID, lw=0.8)
    ax.set_xlabel("certified cheapest loop cost (M, projective model 7.5 l per step)")
    ax.set_ylabel("GLS-2 time / GLV-GLS-4 time")
    ax.grid(color=GRID, lw=0.4)
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "bench_trend.svg"))
    plt.close(fig)


# ---- 3. binary scans: funnel per m, and the loops found
def fig_binary():
    recs = [load(f"binary_scan_m{m}.json") for m in (83, 127, 131)]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.4, 2.9), gridspec_kw={"width_ratios": [1.1, 1]})
    labels = ["scanned D = 1 mod 8", "ord[p2] = m", "prime r", "loop, l <= 61"]
    for i, rec in enumerate(recs):
        vals = [rec["scanned"], rec["order_m_hits"], rec["prime_r_hits"], len(rec["kept"])]
        xs = [j + (i - 1) * 0.26 for j in range(4)]
        bars = a1.bar(xs, vals, width=0.25, label=f"m = {rec['m']} (|D| <= {-rec['d_reached']:,})", color=["#2a6fb0", "#b4442c", "#d08a2c"][i])
        for b, v in zip(bars, vals):
            a1.text(b.get_x() + b.get_width() / 2, v * 1.15, str(v), ha="center", fontsize=6, rotation=90)
    a1.set_yscale("log")
    a1.set_ylim(0.5, 5e7)
    a1.set_xticks(range(4), ["scanned\nD = 1 mod 8", "ord[p2]\n= m", "prime r", "loop with\nl <= 61"], fontsize=7)
    a1.set_ylabel("discriminants (log scale)")
    a1.legend(fontsize=6.5, frameon=False, loc="upper right")
    a1.grid(axis="y", color=GRID, lw=0.4)
    rows = []
    for i, rec in enumerate(recs):
        for k in rec["kept"]:
            rows.append((i, rec["m"], k))
    for y, (i, m, k) in enumerate(rows):
        bl = k["best_loop"]
        a2.barh([y], [bl["model_cost_M"]], color=["#2a6fb0", "#b4442c", "#d08a2c"][i], height=0.6)
        a2.text(bl["model_cost_M"] + 6, y, f"l={bl['l']}, {bl['frobenius_steps']} F", va="center", fontsize=6.5, color=MUTED)
    a2.set_yticks(range(len(rows)), [f"m={m}, D={k['D']}" for _, m, k in rows], fontsize=6.5)
    a2.invert_yaxis()
    a2.set_xlim(0, 480)
    a2.set_xlabel("cheapest loop, model cost (M)")
    a2.grid(axis="x", color=GRID, lw=0.4)
    a2.set_title("m = 83: none of 4 prime-r D has a loop with l <= 61", fontsize=7, color=MUTED)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "binary_scan.svg"))
    plt.close(fig)


# ---- 4. what is computed, validated, timed (status diagram; no data)
def fig_status():
    from matplotlib.patches import FancyBboxPatch

    fig, ax = plt.subplots(figsize=(7.6, 3.0))
    ax.set_xlim(0, 4)
    ax.set_ylim(0, 2.1)
    ax.axis("off")
    W, H = 0.88, 0.62
    boxes = {
        "loop": (0.04, 1.35, "certified cheapest loop\n(exhaustive to the\nnorm bound), native", "#e7eef6"),
        "hd": (1.04, 1.35, "H_D: analytic q-series,\nfixed point, rounding\n< 2^-16 asserted", "#e7eef6"),
        "walk": (2.04, 1.35, "CM curve / F_p, crater\nwalk (Velu), closing\nisomorphism", "#e7eef6"),
        "val": (3.04, 1.35, "twist E'/F_(p^2): psi = +mu,\nloop = \u00b1lambda(\u00b1s)\n24 pts x 5 D: PASS", "#e3f1e3"),
        "bench": (3.04, 0.15, "GLV-GLS-4 vs GLS-2\nvs wNAF timed\n(5 D, pinned core)", "#e3f1e3"),
        "bscan": (1.04, 0.15, "binary GLS: Frobenius-\nassisted loops scanned\nm = 83, 127, 131", "#fbf0dc"),
        "bval": (2.04, 0.15, "binary, x-only check\nm = 31, GF(2^62):\n2 D, PASS", "#fbf0dc"),
    }
    for x, y, text, fc in boxes.values():
        ax.add_patch(FancyBboxPatch((x, y), W, H, boxstyle="round,pad=0.02", fc=fc, ec=GRID))
        ax.text(x + W / 2, y + H / 2, text, ha="center", va="center", fontsize=6.8)

    def arrow(a, b, side=("r", "l"), ls="-"):
        (xa, ya, *_), (xb, yb, *_) = boxes[a], boxes[b]
        pa = {"r": (xa + W, ya + H / 2), "b": (xa + W / 2, ya), "l": (xa, ya + H / 2)}[side[0]]
        pb = {"l": (xb, yb + H / 2), "t": (xb + W / 2, yb + H), "r": (xb + W, yb + H / 2)}[side[1]]
        ax.annotate("", xy=pb, xytext=pa, arrowprops={"arrowstyle": "->", "color": MUTED, "ls": ls})

    arrow("loop", "hd")
    arrow("hd", "walk")
    arrow("walk", "val")
    arrow("val", "bench", ("b", "t"))
    arrow("hd", "bscan", ("b", "t"), "--")
    arrow("bscan", "bval")
    ax.text(0.04, 0.55, "blue: computed natively\ngreen: validated / measured\namber: binary field (char 2)\ndashed: H_D reduced mod 2", fontsize=6.5, color=MUTED, va="center")
    fig.savefig(os.path.join(HERE, "status.svg"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig_bench_arms()
    fig_bench_trend()
    fig_binary()
    fig_status()

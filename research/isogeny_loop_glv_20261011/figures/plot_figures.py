#!/usr/bin/env python3
"""Lay out the recorded values of this study as SVG figures.

Plotting only (AGENTS.md "Implementation language: no Python"; research-visuals
skill): every number drawn here is read from a JSON record written by a native
Rust example, and nothing is computed beyond unit-free ratios already present
in those records.

    /tmp/report-venv/bin/python research/isogeny_loop_glv_20261011/figures/plot_figures.py

Inputs:  ../degree_bound.json   (examples/isogeny_loop_degree_bound.rs)
         ../glvgls_lattice.json (examples/isogeny_loop_glv_gls_lattice.rs)
Outputs: degree_bound.svg, glvgls_bits.svg, glvgls_model.svg, construction.svg
"""
import json
import os

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
# reference palette (dataviz skill, light mode), fixed slot order; validated
SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
MARKERS = ["o", "s", "D", "^"]

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2,
    "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none",
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
})


def load(name):
    with open(os.path.join(STUDY, name)) as fh:
        return json.load(fh)


def short(name):
    return {"id-GostR3410-2001-CryptoPro-B-ParamSet": "CryptoPro-B",
            "id-GostR3410-2001-CryptoPro-A-ParamSet": "CryptoPro-A",
            "icv1-f2m83-tm6151469093347-debefd74": "icv1-f2m83-tm…debefd74",
            "icv1-f2m83-t6151469093347-cdcc5432": "icv1-f2m83-t…cdcc5432"}.get(name, name)


def fig_degree_bound():
    d = load("degree_bound.json")
    rows = [r for r in d["rows"] if r["log2_dmin"] > 2]  # curves with a non-trivial bound
    rows.sort(key=lambda r: (r["family"], r["log2_r"]))
    fig, ax = plt.subplots(figsize=(7.2, 0.34 * len(rows) + 2.2))
    ys = range(len(rows))
    series = [("bound, projective (paper)", lambda r: r["bounds"]["projective"]["ratio"]),
              ("bound, affine", lambda r: r["bounds"]["affine"]["ratio"]),
              ("bound, optimistic x-only", lambda r: r["bounds"]["xonly"]["ratio"]),
              ("cheapest loop found (projective)", lambda r: (r["cheapest_found_M"] or float("nan")) / r["saving_M"])]
    for i, (label, f) in enumerate(series):
        ax.scatter([f(r) for r in rows], list(ys), s=46, marker=MARKERS[i], color=SLOT[i], label=label,
                   edgecolors=SURFACE, linewidths=1.5, zorder=3)
    ax.axvline(1.0, color=INK, lw=1, zorder=2)
    ax.text(1.06, -0.9, "loop cost = doublings saved", color=INK, fontsize=8, va="center")
    ax.set_ylim(-1.4, len(rows) - 0.4)
    ax.set_xscale("log")
    ax.set_yticks(list(ys), [f"{short(r['curve'])} ({r['family']}, {r['log2_r']}-bit r)" for r in rows])
    ax.set_xlabel("loop cost / doublings a 2-dimensional GLV split saves (log scale)")
    ax.grid(axis="x", color=GRID, lw=0.6, zorder=0)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=2, frameon=False)
    ax.set_title("Class-group-free lower bound vs cheapest loop found\n"
                 "source: research/isogeny_loop_glv_20261011/degree_bound.json", fontsize=9, loc="left", color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "degree_bound.svg"))
    plt.close(fig)


def fig_glvgls():
    d = load("glvgls_lattice.json")
    rows = [r for r in d["rows"] if r.get("instance")]  # D = -4 has none (see report)
    discs = []
    for r in rows:
        if r["D"] not in discs:
            discs.append(r["D"])
    by = {(r["D"], r["embedding"]): r for r in rows}
    # bits figure
    fig, ax = plt.subplots(figsize=(7.2, 0.46 * len(discs) + 2.3))
    ys = list(range(len(discs)))
    pts = [("GLS-2 (psi only)", lambda D: by[(D, "s")]["gls2"]["mean_bits"]),
           ("GLV-GLS-4, embedding s", lambda D: by[(D, "s")]["glvgls4"]["mean_bits"]),
           ("GLV-GLS-4, embedding -s", lambda D: by[(D, "-s")]["glvgls4"]["mean_bits"])]
    offs = (0.0, 0.17, -0.17)  # vertical offsets so the two embeddings stay visible
    for i, (label, f) in enumerate(pts):
        ax.scatter([f(D) for D in discs], [y + offs[i] for y in ys], s=40, marker=MARKERS[i], color=SLOT[i],
                   label=label, edgecolors=SURFACE, linewidths=1.5, zorder=3)
    lr = by[(discs[0], "s")]["instance"]["log2_r"]
    for x, lab in ((lr / 2, "log2 r / 2"), (lr / 4, "log2 r / 4")):
        ax.axvline(x, color=INK2, lw=0.8, ls="--", zorder=1)
        ax.text(x - 0.8, -0.85, lab, color=INK2, fontsize=8, ha="right", va="center")
    ax.set_ylim(-1.2, len(discs) - 0.5)
    ax.set_yticks(ys, [f"D = {D}" for D in discs])
    ax.set_xlabel(f"mean over {d['samples']} seeded scalars of max_i bits(|k_i|)  (r of about {lr:.0f} bits, p of {d['bits']} bits)")
    ax.grid(axis="x", color=GRID, lw=0.6, zorder=0)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3, frameon=False)
    ax.set_title("Decomposition coefficient size on the GLS twist E'(F_(p^2))\n"
                 "source: research/isogeny_loop_glv_20261011/glvgls_lattice.json", fontsize=9, loc="left", color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "glvgls_bits.svg"))
    plt.close(fig)
    # model figure: derived operation counts, one point per D (embedding s)
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    xs = [by[(D, "s")]["loop"]["cost_M"] for D in discs]
    ys2 = [by[(D, "s")]["model"]["ratio"] for D in discs]
    ax.scatter(xs, ys2, s=40, marker="o", color=SLOT[0], edgecolors=SURFACE, linewidths=1.5, zorder=3)
    for x, y, D in zip(xs, ys2, discs):
        off = (5, -12) if D == -8 else (5, 4)  # D = -7 and -8 share a cost
        ax.annotate(f"D = {D}", (x, y), textcoords="offset points", xytext=off, fontsize=8, color=INK)
    ax.axhline(1.0, color=INK, lw=1)
    ax.set_xscale("log")
    ax.set_xlabel("cost of the cheapest loop endomorphism, F_(p^2) multiplications (log scale)")
    ax.set_ylabel("model cost GLS-2 / GLV-GLS-4")
    ax.grid(color=GRID, lw=0.6, zorder=0)
    ax.set_title("Derived operation-count model (not a measurement): w = 4, c_dbl = 8, c_add = 12, c_psi = 4\n"
                 "source: research/isogeny_loop_glv_20261011/glvgls_lattice.json", fontsize=9, loc="left", color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "glvgls_model.svg"))
    plt.close(fig)


def fig_construction():
    fig, ax = plt.subplots(figsize=(7.4, 4.0))
    ax.set_xlim(0, 10.6)
    ax.set_ylim(0, 5.6)
    ax.axis("off")
    W, H = 2.5, 1.25

    def box(x, y, text):
        ax.add_patch(FancyBboxPatch((x, y), W, H, boxstyle="round,pad=0.05,rounding_size=0.12",
                                    fc=SURFACE, ec=INK2, lw=1.2))
        ax.text(x + W / 2, y + H / 2, text, ha="center", va="center", fontsize=7.6, color=INK)

    def harrow(x0, x1, y, above, below=""):
        ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>", mutation_scale=10, color=INK2, lw=1.1))
        ax.text((x0 + x1) / 2, y + 0.12, above, ha="center", va="bottom", fontsize=7, color=INK2)
        if below:
            ax.text((x0 + x1) / 2, y - 0.12, below, ha="center", va="top", fontsize=6.6, color=INK2, style="italic")

    def varrow(x, y0, y1, label):
        ax.add_patch(FancyArrowPatch((x, y0), (x, y1), arrowstyle="-|>", mutation_scale=10, color=INK2, lw=1.1))
        ax.text(x + 0.1, (y0 + y1) / 2, label, ha="left", va="center", fontsize=7, color=INK2)

    xs = (0.1, 4.05, 8.0)
    top, bot = 3.85, 1.35
    box(xs[0], top, "E / F_p\nCM by O_K, trace t\n4p = t^2 + |D| v^2")
    box(xs[1], top, "loop of horizontal\nl-isogenies on the crater\nalpha in O_K, deg N(alpha)")
    box(xs[2], top, "eigenvalue on E'[r]\nlambda = (A + B s)/2\ns = (2 mu - t)/v")
    box(xs[0], bot, "E' / F_(p^2)\nquadratic twist\n#E' = (p-1)^2 + t^2 = h r")
    box(xs[1], bot, "psi = phi o pi_p o phi^-1\npsi^2 = -1\nmu = (p-1)/t mod r")
    box(xs[2], bot, "4-dim lattice\nk0 + k1 mu + k2 lambda\n+ k3 lambda mu = 0 mod r")
    yt, yb = top + H / 2, bot + H / 2
    harrow(xs[0] + W + 0.08, xs[1] - 0.08, yt, "cheapest loop", "certified, native")
    harrow(xs[1] + W + 0.08, xs[2] - 0.08, yt, "ring map", "derived, asserted")
    harrow(xs[0] + W + 0.08, xs[1] - 0.08, yb, "GLS map")
    harrow(xs[1] + W + 0.08, xs[2] - 0.08, yb, "with lambda", "LLL, native")
    varrow(xs[0] + W / 2, top - 0.08, bot + H + 0.08, "twist phi")
    varrow(xs[2] + W / 2, top - 0.08, bot + H + 0.08, "")
    ax.text(5.3, 0.55, "alpha' = phi o alpha o phi^-1 is F_(p^2)-rational. Its action on points of a constructed E' is not\n"
            "validated here (no curve is built): the eigenvalue lambda is derived from the ring map, status 'not point-validated'.",
            ha="center", va="center", fontsize=7.2, color=INK)
    ax.set_title("GLV-GLS with an isogeny-loop endomorphism: what is computed and what is only derived",
                 fontsize=9, loc="left", color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "construction.svg"))
    plt.close(fig)


if __name__ == "__main__":
    fig_construction()
    fig_degree_bound()
    if os.path.exists(os.path.join(STUDY, "glvgls_lattice.json")):
        fig_glvgls()

#!/usr/bin/env python3
"""Plot recorded PS-4 values (raw/plot_data.json, emitted by `pdp_ps4 summary --data`).

Python here only lays out what the Rust binary recorded and derived; it computes nothing.

  plot_ps4.py raw/plot_data.json OUT_DIR
writes OUT_DIR/ps4_candidate_factor.svg, ps4_decision.svg, ps4_leaf.svg, ps4_mechanism.svg
"""
import json
import sys

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "sans-serif"

# Categorical slots 1-3 of the validated reference palette (light surface), fixed order.
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#dcdcd8"


def style(ax):
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=7)
    ax.yaxis.label.set_color(INK)
    ax.xaxis.label.set_color(INK)


def label(s):
    return f"n={s['n']}\nw={s['w']}"


def main():
    data = json.load(open(sys.argv[1]))
    out = sys.argv[2]
    S = sorted(data["settings"], key=lambda s: (s["w"], s["n"]))
    L = sorted(data["leaf"], key=lambda s: (s["w"], s["n"]))
    T = sorted(data["torsion"], key=lambda s: (s["w"], s["n"]))

    # ---- Figure 1: candidate factor against n --------------------------------
    fig, ax = plt.subplots(figsize=(7.2, 3.9))
    ns = [5 + i for i in range(0, 40)]
    ax.plot([10, 43], [10, 43], color=INK2, lw=0.9, ls="--")
    ax.plot([10, 43], [5, 21.5], color=INK2, lw=0.9, ls=":")
    ax.text(43.3, 43, "n", fontsize=7, color=INK2, va="center")
    ax.text(43.3, 21.5, "n/2", fontsize=7, color=INK2, va="center")
    for w, mk in ((2, "o"), (3, "s")):
        rows = [s for s in S if s["w"] == w]
        x = [s["n"] for s in rows]
        ax.plot(x, [s["pair_factor_ord"] for s in rows], marker=mk, color=BLUE, lw=1.2, ms=5,
                label=f"w={w}: pairs, ordered unreduced / reduced")
        ax.plot(x, [s["pair_factor_unord"] for s in rows], marker=mk, color=ORANGE, lw=1.2, ms=5,
                label=f"w={w}: pairs, unordered unreduced / reduced")
        ax.plot(x, [s["M_factor_ord"] for s in rows], marker=mk, color=BLUE, lw=0, ms=9, mfc="none", mew=0.8)
        ax.plot(x, [s["M_factor_unord"] for s in rows], marker=mk, color=ORANGE, lw=0, ms=9, mfc="none", mew=0.8)
    tx = [t["n"] for t in T if t["w"] == 3]
    ax.plot([t["n"] for t in T], [t["pair_factor"] for t in T], marker="^", color=AQUA, lw=0, ms=5,
            label="2-torsion on the T-closed base B+-, pairs (ordered)")
    ax.set_xlabel("field degree n (type-II ONB)")
    ax.set_ylabel("factor (unreduced / reduced), per target")
    ax.set_xlim(9, 45.5)
    ax.set_ylim(0, 46)
    style(ax)
    ax.legend(fontsize=6.3, loc="upper left", frameon=False)
    ax.set_title("Candidate pairs fall by about n (ordered) or n/2 (unordered); the 2-torsion factor is 2 and only on B+-\n"
                 "(filled markers: candidate pairs; hollow rings: counted field multiplications M, same settings; circles w=2, squares w=3)",
                 fontsize=7.5, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(f"{out}/ps4_candidate_factor.svg")
    plt.close(fig)

    # ---- Figure 2: what the reduction decides ----------------------------------
    fig, axs = plt.subplots(3, 1, figsize=(7.2, 7.0), sharex=True)
    x = list(range(len(S)))
    xl = [label(s) for s in S]
    ax = axs[0]
    w_ = 0.38
    ax.bar([i - w_ / 2 for i in x], [s["verdict_kept_frac"] for s in S], w_, color=BLUE,
           label="reduced loop on the target alone: verdict equals the unreduced one")
    ax.bar([i + w_ / 2 for i in x], [s["closure_eq_frac"] for s in S], w_, color=ORANGE,
           label="reduced loop over all conjugates of the target: solution set equals the unreduced one")
    ax.set_ylabel("fraction of the\nsetting's targets")
    ax.set_ylim(0, 1.18)
    style(ax)
    ax.legend(fontsize=6.3, loc="upper center", ncol=1, frameon=False)
    ax.set_title("(a) verdicts and solution sets, reduced against unreduced (8 targets per setting)", fontsize=7.5, loc="left")
    ax = axs[1]
    ax.bar(x, [s["kept_sols_frac"] for s in S], 0.55, color=BLUE, label="solutions the reduced loop returns / solutions of the oracle (target alone)")
    ax.set_ylabel("fraction of\nsolutions")
    ax.set_ylim(0, 1.0)
    style(ax)
    ax.legend(fontsize=6.3, loc="upper right", frameon=False)
    ax.set_title("(b) the solutions the reduced loop keeps on a single target", fontsize=7.5, loc="left")
    ax = axs[2]
    ax.bar([i - w_ / 2 for i in x], [s["dec_fw_ratio"] for s in S], w_, color=BLUE,
           label="decomposable targets, to first witness: M reduced / M unreduced (unordered)")
    ax.bar([i + w_ / 2 for i in x], [s["closure_M_over_ord"] for s in S], w_, color=ORANGE,
           label="all conjugates of the target, full scan: M reduced / M unreduced (ordered)")
    ax.axhline(1.0, color=INK2, lw=0.9, ls="--")
    ax.set_ylabel("ratio of counted\nmultiplications M")
    ax.set_ylim(0, 1.5)
    ax.set_xticks(x)
    ax.set_xticklabels(xl, fontsize=6.5)
    style(ax)
    ax.legend(fontsize=6.3, loc="upper left", frameon=False)
    ax.set_title("(c) cost of deciding: ratio 1 means no saving", fontsize=7.5, loc="left")
    fig.tight_layout()
    fig.savefig(f"{out}/ps4_decision.svg")
    plt.close(fig)

    # ---- Figure 3: leaf cost per pair, weight base against the subspace control ---
    fig, axs = plt.subplots(2, 1, figsize=(7.2, 5.0), sharex=True)
    x = list(range(len(L)))
    xl = [f"n={s['n']}\nw={s['w']}\n(l={s['near_l']})" for s in L]
    ax = axs[0]
    ax.bar([i - w_ / 2 for i in x], [s["M_pair_onb"] for s in L], w_, color=BLUE, label="weight base, ONB (this task)")
    ax.bar([i + w_ / 2 for i in x], [s["M_pair_ctl"] for s in L], w_, color=ORANGE, label="subspace control, nearest-n rung l, polynomial basis (PS-2 kernel)")
    ax.set_ylabel("counted field\nmultiplications M\nper pair")
    ax.set_ylim(0, 13)
    style(ax)
    ax.legend(fontsize=6.3, loc="lower center", frameon=False, ncol=1)
    ax.set_title("(a) counted multiplications per pair (inversion products and batch products included)", fontsize=7.5, loc="left")
    ax = axs[1]
    ax.bar([i - w_ / 2 for i in x], [s["bits_ratio"] for s in L], w_, color=BLUE, label="model bit operations, weight / control, each at its own degree")
    ax.bar([i + w_ / 2 for i in x], [s["bits_ratio_scaled"] for s in L], w_, color=AQUA,
           label="same, control scaled to the weight base's degree by (n/n_l)^2")
    ax.axhline(2.0, color=INK, lw=1.0, ls="--")
    ax.text(len(L) - 0.5, 2.05, "stop condition: ratio > 2", fontsize=6.5, ha="right", va="bottom", color=INK)
    ax.axhline(1.0, color=INK2, lw=0.7, ls=":")
    ax.set_ylabel("weight leaf / subspace leaf\nper pair (model)")
    ax.set_xticks(x)
    ax.set_xticklabels(xl, fontsize=6.3)
    style(ax)
    ax.legend(fontsize=6.3, loc="upper left", frameon=False)
    ax.set_title("(b) the same pairs priced in the bit-operation model (a convention, model.rs)", fontsize=7.5, loc="left")
    fig.tight_layout()
    fig.savefig(f"{out}/ps4_leaf.svg")
    plt.close(fig)

    # ---- Figure 4: the accounting, as a diagram ----------------------------------
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 52)
    ax.axis("off")

    def box(x0, y0, w, h, text, fc, fs=6.6, ec=INK2):
        ax.add_patch(FancyBboxPatch((x0, y0), w, h, boxstyle="round,pad=0.3,rounding_size=1.2", fc=fc, ec=ec, lw=0.8))
        ax.text(x0 + w / 2, y0 + h / 2, text, ha="center", va="center", fontsize=fs, color=INK)

    box(1, 34, 30, 15, "unreduced\none target R, every base point as\nfirst summand P1\npairs m(m+1)/2 (or m^2)\nfinds all of S(R)", "#e3eefb")
    box(35, 34, 30, 15, "reduced, target alone\none representative per orbit as P1\npairs reps x m  (about m^2/n)\nfinds S(R) restricted to triples\nwith a representative member", "#fdebe3")
    box(69, 34, 30, 15, "reduced, all conjugates of R\nR, R^2, R^4, ... each with\nrepresentatives\npairs about m^2: equals the\nunreduced ordered loop", "#e1f5ec")
    ax.annotate("", xy=(35, 41.5), xytext=(31.5, 41.5), arrowprops=dict(arrowstyle="->", color=INK2))
    ax.annotate("", xy=(69, 41.5), xytext=(65.5, 41.5), arrowprops=dict(arrowstyle="->", color=INK2))
    ax.text(33.3, 44.2, "n fewer\npairs", ha="center", fontsize=6, color=INK2)
    ax.text(67.3, 44.2, "x n\ntargets", ha="center", fontsize=6, color=INK2)
    ax.text(50, 28.5, "Frobenius maps a decomposition of R to a decomposition of its conjugate, not of R itself:\n"
            "S(R^2) = sigma S(R).   The target is fixed, so the symmetry trades targets for first summands one for one.",
            ha="center", va="top", fontsize=7, color=INK)
    box(1, 3, 47, 14, "2-torsion T = (0,1):  x -> 1/x,  u = 1/(x+1): u -> u + 1\n"
        "(P1+T) + (P2+T) + P3 = R  for the same R:\na symmetry of a fixed target, on a T-closed base only\n"
        "weight base: u + 1 has weight n - wt(u) > w: not closed", "#fdebe3", fs=6.4)
    box(52, 3, 47, 14, "on B+- (weight base and its complements, in u):\n"
        "first summand one per class {x, 1/x}  (+ T alone)\n"
        "pairs fall by 2 exactly; closure under translating two\nsummands by T equals the unreduced set", "#e1f5ec", fs=6.4)
    ax.set_title("What each reduction decides (mechanism; the measured values are in the other figures)", fontsize=8, loc="left", color=INK)
    fig.tight_layout()
    fig.savefig(f"{out}/ps4_mechanism.svg")
    plt.close(fig)


if __name__ == "__main__":
    main()

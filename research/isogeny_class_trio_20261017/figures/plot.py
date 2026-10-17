"""Lay out recorded values from ../data as SVG + PNG figures. No computation
beyond reading the JSON written by the native selector (AGENTS.md: Python only
plots recorded values)."""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
COL = {"K": "#2a78d6", "B": "#eb6834", "P": "#1baf7a"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "svg.fonttype": "none",
                     "axes.spines.top": False, "axes.spines.right": False})


def load(name):
    return json.loads((DATA / name).read_text())


def save(fig, stem):
    fig.savefig(HERE / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(HERE / f"{stem}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def strata():
    picks = [("K61", "K-61  (f = 1951·587551)"), ("B61", "B-61  (f = 7²)"), ("P61", "P-61  (f = 3²·5³)"),
             ("K41", "K-41  (f = 409·1721)"), ("B47", "B-47  (f = 31²)"), ("P47", "P-47  (f = 3⁵)")]
    fig, axes = plt.subplots(2, 3, figsize=(10, 5.6))
    for ax, (tag, title) in zip(axes.flat, picks):
        r = load(f"pick-{tag}.json")
        s = r["strata"]
        hs = [int(x["h"]) for x in s]
        labels = [x["conductor_str"] for x in s]
        # the root's stratum: crater (c = 1) for Koblitz, floor (c = f) for the seeded roots
        root = 0 if tag.startswith("K") else len(s) - 1
        colors = [COL[tag[0]] if i == root else "#b9c7d8" if tag[0] == "K" else
                  "#f3c3ad" if tag[0] == "B" else "#a9dcc6" for i in range(len(s))]
        ax.bar(range(len(s)), hs, color=colors, width=0.7, edgecolor="white", linewidth=2)
        ax.set_yscale("log")
        ax.set_xticks(range(len(s)))
        ax.set_xticklabels(labels, rotation=60 if len(s) > 6 else 0, fontsize=7)
        ax.set_title(f"{title}\nH = {int(r['class_size']):,} curves", fontsize=9, color=INK)
        ax.grid(axis="y", color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.annotate("root", (root, hs[root]), textcoords="offset points", xytext=(0, 4),
                    ha="center", fontsize=7, color=INK)
    for ax in axes[:, 0]:
        ax.set_ylabel("curves in stratum  h(c² d_K)")
    for ax in axes[1]:
        ax.set_xlabel("stratum: conductor c of End(E), c | f")
    fig.suptitle("Strata of each isogeny class (exact class numbers; top row = the 61-bit picks, bottom = the 41/47-bit trio)",
                 fontsize=10, color=INK)
    fig.tight_layout()
    save(fig, "strata")


def search():
    runs = [("B47", "B-47, 131072 seeds"), ("P47", "P-47, 131072 seeds"),
            ("B61", "B-61, 32768 seeds"), ("P61", "P-61, 32768 seeds")]
    fig, ax = plt.subplots(figsize=(7.5, 3.4))
    width = 0.2
    shades = ["#eb6834", "#1baf7a", "#f3a27f", "#7fd1b0"]
    for i, (tag, label) in enumerate(runs):
        t = load(f"tally-{tag}.json")
        h = {int(k): v for k, v in t["max_depth_histogram"].items()}
        xs = sorted(h)
        ax.bar([x + (i - 1.5) * width for x in xs], [h[x] for x in xs], width=width * 0.9,
               color=shades[i], label=f"{label} ({t['prime_subgroup_curves']} prime-subgroup curves)")
        for x in xs:
            ax.annotate(str(h[x]), (x + (i - 1.5) * width, h[x]), textcoords="offset points",
                        xytext=(0, 2), ha="center", fontsize=6, color=INK2)
    ax.set_yscale("log")
    ax.set_xlabel("max over l of volcano height v_l(f)")
    ax.set_ylabel("curves (log scale)")
    ax.set_title("Seeded searches: volcano height among standard-cofactor prime-subgroup curves", fontsize=9)
    ax.legend(fontsize=7, frameon=False)
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    save(fig, "search-depth")


def koblitz():
    rows = [json.loads(line) for line in (DATA / "koblitz.jsonl").read_text().splitlines()]
    rows = [r for r in rows if r["m_prime"]]
    fig, ax = plt.subplots(figsize=(7, 3.3))
    for a, mk in [(0, "o"), (1, "s")]:
        sub = [r for r in rows if r["a"] == a]
        ax.scatter([r["m"] for r in sub], [r["r_bits"] for r in sub], marker=mk, s=40,
                   facecolors=[COL["K"] if r["near_prime"] else "white" for r in sub],
                   edgecolors=COL["K"], label=f"a = {a}")
        for r in sub:
            ax.annotate(f"h={r['cofactor']}", (r["m"], r["r_bits"]), textcoords="offset points",
                        xytext=((0, -12) if r["r_bits"] < 40 else (-6, 5)) if r["m"] >= 59 else (6, -3),
                        ha=(("right" if r["m"] == 59 else "center") if r["r_bits"] < 40 else "right") if r["m"] >= 59 else "left",
                        fontsize=6, color=INK2)
    ax.plot([41, 61], [41, 61], color=GRID, linewidth=1)
    ax.text(55, 56.5, "r = m bits", fontsize=6, color=INK2, rotation=14)
    ax.set_xlabel("field degree m (prime)")
    ax.set_ylabel("bits of largest prime factor r")
    ax.set_title("Every prime-degree Koblitz curve, m in [41, 66]: filled = standard cofactor (4 or 2)", fontsize=9)
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    ax.grid(color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    save(fig, "koblitz-table")


def volcanoes():
    """Schematic of one volcano component per primary class, from recorded strata."""
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.6))
    specs = [
        ("K61", "K-61: the single crater curve (End = Z[τ])\n1951- and 587551-steps (inert craters)",
         [("crater  O_K", 1), ("O_1951", 1952), ("O_587551", 587552), ("floor  O_f", 1146901504)]),
        ("B61", "B-61: height-2 7-volcano, split crater\nper crater vertex (crater cycles not drawn)",
         [("crater", 1), ("level 1", 6), ("floor", 42)]),
        ("P61", "P-61: 5-volcano (height 3, split) per crater vertex;\neach vertex also heads a 3-tree 1 → 4 → 12 (inert)",
         [("crater", 1), ("level 1", 4), ("level 2", 20), ("floor", 100)]),
    ]
    for ax, (tag, title, tiers) in zip(axes, specs):
        n = len(tiers)
        for i, (name, count) in enumerate(tiers):
            y = n - 1 - i
            w = 0.45 + 0.5 * (i / max(n - 1, 1))
            ax.barh(y, w, left=(1 - w) / 2, height=0.7, color=COL[tag[0]], alpha=0.25 + 0.75 * (1 - i / n),
                    edgecolor="white", linewidth=2)
            ax.text(0.5, y, f"{name}: {count:,}", ha="center", va="center", fontsize=7.5, color=INK)
        ax.annotate("root (seeded)" if tag != "K61" else "root (Koblitz curve)",
                    xy=(0.95 if tag != "K61" else 0.73, 0 if tag != "K61" else n - 1), xytext=(1.05, (n - 1) / 2),
                    fontsize=7, color=INK, ha="left", va="center",
                    arrowprops=dict(arrowstyle="->", color=INK2, linewidth=0.8))
        ax.set_xlim(0, 1.45)
        ax.set_ylim(-0.6, n - 0.4)
        ax.axis("off")
        ax.set_title(title, fontsize=8.5, color=INK)
    fig.suptitle("Volcano strata (vertex counts per tier, from the recorded class numbers)", fontsize=10, color=INK)
    fig.tight_layout()
    save(fig, "volcanoes")


if __name__ == "__main__":
    strata()
    search()
    koblitz()
    volcanoes()

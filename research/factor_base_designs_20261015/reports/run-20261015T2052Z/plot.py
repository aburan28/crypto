"""Plot recorded D and max_cols against b for the Amendment 2 series (F1, F2, F4),
from runs/*.jsonl only (no computation beyond reading the rows and taking log
scale). Rendering aid only (AGENTS.md "no Python")."""
import glob, json, os
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
runs = os.path.join(here, "..", "..", "runs")


def load(pattern):
    out = []
    for f in sorted(glob.glob(os.path.join(runs, pattern))):
        out += [json.loads(l) for l in open(f) if l.strip()]
    return out


# (panel title, [(label, rows, colour, marker, x offset)])
PANELS = [
    ("F1  m = 2,  p2 = 2013265921", [
        ("A0", load("F1-A0-m2-p2.jsonl"), "#2a78d6", "o", -0.06),
        ("A11", load("F1-A11-m2-p2.jsonl"), "#008300", "P", 0.06)]),
    ("F1  m = 3,  p2 = 2013265921", [
        ("A0", load("F1-A0-m3-p2.jsonl"), "#2a78d6", "o", -0.06),
        ("A11", load("F1-A11-m3-p2.jsonl"), "#008300", "P", 0.06)]),
    ("F4  m = 4,  p1 = 786433", [
        ("A0", load("F4-A0-m4.jsonl"), "#2a78d6", "o", -0.06),
        ("A11", load("F4-A11-m4.jsonl"), "#008300", "P", 0.06)]),
    ("F2  m = 2,  p1 = 786433", [
        ("N11", load("F2-N11-m2.jsonl"), "#eb6834", "s", -0.06),
        ("A11 (committed p1 rows)", load("A11-m2.jsonl"), "#008300", "P", 0.06)]),
    ("F2  m = 3,  p1 = 786433", [
        ("N11", load("F2-N11-m3.jsonl"), "#eb6834", "s", -0.06),
        ("A11 (committed p1 rows)", load("A11-m3.jsonl"), "#008300", "P", 0.06)]),
]

plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": "#888", "axes.labelcolor": "#222",
                     "xtick.color": "#444", "ytick.color": "#444"})


def figure(key, ylabel, fname, title, log):
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.6))
    flat = axes.flatten()
    for ax, (ptitle, series) in zip(flat, PANELS):
        for label, rows, col, mk, off in series:
            ok = sorted((r["b"], r[key]) for r in rows if r["status"] == "ok" and r.get(key) is not None)
            ax.plot([b + off for b, _ in ok], [v for _, v in ok], color=col, lw=1.4, marker=mk, ms=5.5,
                    mec="white", mew=0.8, label=label)
            for r in rows:
                if r["status"] == "ok":
                    continue
                if r.get(key) is not None:  # F4 budget stop: value reached (a lower bound)
                    ax.plot([r["b"] + off], [r[key]], marker=mk, ms=7, mfc="none", mec=col, mew=1.4, ls="none")
                else:  # hard kill or memory abort: nothing recorded; a cross on the axis floor
                    ax.plot([r["b"] + off], [0.02], marker="x", ms=7, color=col, mew=1.6, ls="none",
                            transform=ax.get_xaxis_transform())
        ax.set_title(ptitle, color="#222", fontsize=9)
        ax.set_xlabel("b = log2|V| per summand (bits)")
        ax.set_ylabel(ylabel)
        if log:
            ax.set_yscale("log", base=2)
        ax.grid(color="#e5e5e5", lw=0.6)
        ax.set_axisbelow(True)
        ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    note = flat[5]
    note.axis("off")
    note.text(0.0, 0.95,
              "Rows from runs/F1-*, F2-*, F4-* and the committed\n"
              "runs/A11-m2.jsonl, A11-m3.jsonl (F2's reference).\n\n"
              "Filled markers: completed systems, both targets\n"
              "plotted (n = 2 per cell; coincident points overlap).\n"
              "Open markers: systems stopped by the 300 s F4\n"
              "budget, at the value reached (a lower bound).\n"
              "Crosses on the floor: hard kill or memory abort,\n"
              "no value recorded.\n\n"
              "Markers are offset by ±0.06 in b so that\n"
              "overlapping series stay visible.\n\n"
              "run-20261015T2052Z, plotted by plot.py",
              va="top", fontsize=8, color="#444", transform=note.transAxes)
    fig.suptitle(title, fontsize=10, color="#222")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(os.path.join(here, fname + ".svg"))
    fig.savefig(os.path.join(here, fname + ".png"), dpi=150)


figure("D", "D = solving_degree_max", "d_vs_b",
       "Solving degree vs factor-base bits, Amendment 2 series (run-20261015T2052Z)", False)
figure("max_cols", "largest F4 matrix width (max_cols)", "cols_vs_b",
       "Largest F4 matrix width vs factor-base bits, Amendment 2 series (run-20261015T2052Z)", True)
print("wrote d_vs_b and cols_vs_b (svg, png)")

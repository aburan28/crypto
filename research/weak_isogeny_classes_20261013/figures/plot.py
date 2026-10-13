"""Lay out recorded values from ../data/*.json as SVG + 150 dpi PNG.

Plotting only: every number comes from the native run's JSON. No values are
computed here beyond ratios of recorded counts.
"""
import json, pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE.parent / "data"
plt.rcParams.update({"font.size": 9, "svg.fonttype": "none"})
C = {"direct": "#4C6EF5", "class": "#E8590C", "walk": "#2F9E44", "nec": "#868E96"}


def save(fig, name):
    fig.savefig(HERE / f"{name}.svg", bbox_inches="tight")
    fig.savefig(HERE / f"{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def ghs(n_field):
    return json.loads((DATA / f"ghs_N{n_field}.json").read_text())


# 1. Coverage against the magic-number threshold for selected (N, n, l).
picks = [(16, 4, 4), (20, 4, 5), (24, 4, 6), (18, 6, 3), (21, 7, 3), (24, 8, 3)]
fig, axes = plt.subplots(2, 3, figsize=(10, 6), sharey=True)
for ax, (N, n, l) in zip(axes.flat, picks):
    rep = ghs(N)
    row = next(f for f in rep["factorisations"] if f["n"] == n)
    th = [t for t in row["thresholds"] if t["max_m"] >= 2]
    M = [t["max_m"] for t in th]
    tot = th[0]["curves_total"]
    ax.plot(M, [t["direct_weak_curves"] / tot for t in th], "o-", c=C["direct"], label="direct: 2 ≤ m ≤ M")
    ax.plot(M, [t["class_covered_curves"] / tot for t in th], "s-", c=C["class"], label="class (same trace) has 2 ≤ m ≤ M")
    if th[0]["walk3_covered_curves"] is not None:
        ax.plot(M, [t["walk3_covered_curves"] / tot for t in th], "^-", c=C["walk"], label="joined by 3-isogeny/Frobenius walk")
    ax.set_yscale("log")
    ax.set_title(f"N={N}, n={n}, l={l}  ({tot} curves)")
    ax.set_xlabel("magic-number threshold M")
    ax.grid(alpha=.3)
for ax in axes[:, 0]:
    ax.set_ylabel("fraction of curves over F_2^N")
axes.flat[0].legend(fontsize=7, loc="lower right")
fig.suptitle("GHS weak set: curves vs isogeny classes (exhaustive, a=0 family + twist under MMT (8))", fontsize=10)
fig.tight_layout()
save(fig, "ghs_coverage_vs_threshold")

# 2. Enlargement factor class/direct at the smallest non-subfield threshold
#    m <= 3 (genus <= 4) for every (N, n) with n >= 3.
fig, ax = plt.subplots(figsize=(7, 4))
for n, mk in [(4, "s"), (6, "^"), (8, "v")]:
    xs, ys, ws = [], [], []
    for N in [8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 22, 24]:
        p = DATA / f"ghs_N{N}.json"
        if not p.exists():
            continue
        rep = json.loads(p.read_text())
        row = next((f for f in rep["factorisations"] if f["n"] == n), None)
        if row is None:
            continue
        t = next(t for t in row["thresholds"] if t["max_m"] == 3) if n >= 3 else None
        if t and t["direct_weak_curves"]:
            xs.append(N)
            ys.append(t["class_covered_curves"] / t["direct_weak_curves"])
            if t["walk3_covered_curves"] is not None:
                ws.append((N, t["walk3_covered_curves"] / t["direct_weak_curves"]))
    ax.plot(xs, ys, mk + "-", label=f"n={n}: class / direct")
    if ws:
        ax.plot([w[0] for w in ws], [w[1] for w in ws], mk + ":", alpha=.6, label=f"n={n}: 3-walk / direct")
ax.set_yscale("log")
ax.set_xlabel("N (field F_2^N)")
ax.set_ylabel("curves covered ÷ directly weak curves")
ax.set_title("Enlargement of the GHS weak set at 2 ≤ m ≤ 3 (genus ≤ 4)\n(n = 3, 5 omitted: no b has m = 2, and m = 3 = n is every non-subfield curve)", fontsize=9)
ax.grid(alpha=.3)
ax.legend(fontsize=7, ncol=2)
save(fig, "ghs_enlargement_m3")

# 3. Smith: direct vs class closure per p, with S_8 Chebotarev constants.
cheb = json.loads((DATA / "smith_chebotarev.json").read_text())
ps = [p for p in (3, 5, 7, 11, 13) if (DATA / f"smith_p{p}.json").exists()]
reps = [json.loads((DATA / f"smith_p{p}.json").read_text()) for p in ps]
fig, ax = plt.subplots(figsize=(7, 4))
w = 0.2
x = range(len(ps))
ser = [
    ("rational tractable subgroup (direct)", "direct_pair_kernel_curves", C["direct"]),
    ("class contains one (hyperelliptic)", "class_pair_kernel_curves", C["class"]),
    ("any stable Lagrangian (135)", "direct_any_lagrangian_curves", C["walk"]),
    ("χ mod 2 = g·g* (necessary)", "necessary_condition_curves", C["nec"]),
]
for k, (lab, key, col) in enumerate(ser):
    ax.bar([i + (k - 1.5) * w for i in x], [r[key] / r["curves"] for r in reps], w, label=lab, color=col)
ax.axhline(cheb["s8_pair_fraction"], ls="--", c="k", lw=.8)
ax.text(len(ps) - .5, cheb["s8_pair_fraction"] + .01, f"S8 limit {cheb['s8_pair_fraction']:.4f}", ha="right", fontsize=7)
ax.axhline(0.18574, ls=":", c="k", lw=.8)
ax.text(len(ps) - .5, 0.19, "Smith success 0.1857 (Hyp. 2–3)", ha="right", fontsize=7)
ax.set_xticks(list(x))
ax.set_xticklabels([f"p={p}\n{r['curves']} curves\n{r['classes']} classes" for p, r in zip(ps, reps)])
ax.set_ylabel("fraction of genus-3 models y² = f(x), deg f = 8")
ax.set_ylim(0, .95)
ax.set_title("Smith kernels: per curve vs per isogeny class (exhaustive)")
ax.legend(fontsize=7, loc="upper center", ncol=2)
save(fig, "smith_coverage")

# 4. Pipeline diagram with status per edge.
fig, ax = plt.subplots(figsize=(9, 3.2))
ax.axis("off")
boxes = [
    (0.06, "E / F_2^N\n(any curve)"),
    (0.28, "E′ ~ E, m(b′) ≤ M\n(same trace)"),
    (0.50, "C / F_2^l\nGHS, genus 2^(m−1)"),
    (0.72, "H genus 3 hyp.\n(char ≥ 5 needed)"),
    (0.94, "plane quartic\nDiem Õ(q)"),
]
for xc, t in boxes:
    ax.text(xc, .55, t, ha="center", va="center", fontsize=8,
            bbox=dict(boxstyle="round", fc="#F1F3F5", ec="#495057"))
edges = [
    (0.06, 0.28, "isogeny walk\nmeasured: class & 3-walk\ncoverage (this study)", C["walk"]),
    (0.28, 0.50, "GHS descent\nexisting: m=1 executable,\nm≥2 symbolic", C["nec"]),
    (0.50, 0.72, "BLOCKED: char 2\n(problem record)", "#C92A2A"),
    (0.72, 0.94, "Smith (2,2,2)\nexistence measured,\nsquare tests not built", C["class"]),
]
for a, b, t, col in edges:
    ax.annotate("", xy=(b - .085, .55), xytext=(a + .085, .55),
                arrowprops=dict(arrowstyle="->", color=col, lw=1.6))
    ax.text((a + b) / 2, .18, t, ha="center", va="center", fontsize=7, color=col)
ax.set_title("Where the two leads meet: pipeline and status of each step", fontsize=10)
save(fig, "pipeline_status")
print("ok")

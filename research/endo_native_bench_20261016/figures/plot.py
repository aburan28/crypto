"""Lay out the committed records of research/endo_native_bench_20261016
(AGENTS.md: Python only plots; every number comes from the native
`endo-chain-bench` JSON records beside this directory, nothing is re-run).

Reads only `../<curve>.counts.json`, `../<curve>.timing.json` and
`../<curve>.isolated_bench_record.jsonl`.  Derived values are ratios of
recorded values, computed exactly as the README defines them:

- counted ratio = min counted M_eq over baseline arms / min counted M_eq over
  GLV arms (the README's "best baseline / best GLV-2");
- model ratio   = the same with the frozen model's `model_m_eq` per arm;
- timed ratios vs the reference are the recorded `ratio_reference_over_arm`.

Writes, beside this file (SVG and a 150-dpi PNG of each):
- fig_counted_vs_model  - counted vs model ratio per curve, and every total
                          arm's counted deviation from the model (%);
- fig_timed             - declared GLV arm vs the reference: counted, median
                          and min ratio with the A/A floor, and per-round
                          wall times normalised by the reference median;
- fig_grid              - per-width / per-table grid of the declared chain
                          (counted M_eq and median us);
- fig_phi               - phi(P) stage diagnostics: counted M_eq vs median ns.
"""
import json
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

here = Path(__file__).resolve().parent
study = here.parent
CURVES = [("bandersnatch", "Bandersnatch (SW)"), ("tom384", "Tom-384"), ("cp6_782", "CP6-782 G1")]

INK, MUTED, GRID = "#1f2328", "#57606a", "#d0d7de"
C_GLV, C_GLV2, C_BASE, C_BASE2, C_AA = "#0072b2", "#56b4e9", "#d55e00", "#e69f00", "#6e7781"
plt.rcParams.update({
    "font.size": 8.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "svg.fonttype": "none", "svg.hashsalt": "endo-native-bench",
    "legend.frameon": False,
})
CAPTION = ("endo-chain-bench, commit 33e2dd10b; 256 (k, P) pairs, seed 20261016; "
           "timing 20 rounds on CPU 3 via isolated_bench (L0 cloud VM)")


def load(curve):
    counts = json.loads((study / f"{curve}.counts.json").read_text())
    timing = json.loads((study / f"{curve}.timing.json").read_text())
    iso = json.loads((study / f"{curve}.isolated_bench_record.jsonl").read_text().splitlines()[0])
    return counts, timing, iso


def by_key(arms):
    return {a["key"]: a for a in arms}


def save(fig, name):
    fig.savefig(here / f"{name}.svg", metadata={"Date": None})
    fig.savefig(here / f"{name}.png", dpi=150)
    plt.close(fig)


def caption(fig, text=CAPTION):
    fig.text(0.01, 0.005, textwrap.fill(text, 160), fontsize=6.5, color=MUTED, ha="left", va="bottom")


DATA = {c: load(c) for c, _ in CURVES}


def declared(timing):
    return f"glv/{timing['inputs']['grid_chain']}/optimised/jacobian/w5"


def aa_key(timing):
    return next(a["key"] for a in timing["arms"] if a["key"].startswith("aa/"))


def reference(timing):
    tab, w = timing["inputs"]["reference"].split("/")
    return f"baseline/{tab}/w{w}"


# -- fig_counted_vs_model -----------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.0, 3.6), gridspec_kw={"width_ratios": [1, 1.5]})
xs = range(len(CURVES))
counted_r, model_r = [], []
for curve, _ in CURVES:
    arms = DATA[curve][0]["arms"]
    base = [a for a in arms if a["key"].startswith("baseline/")]
    glv = [a for a in arms if a["key"].startswith("glv/")]
    counted_r.append(min(a["m_eq"] for a in base) / min(a["m_eq"] for a in glv))
    model_r.append(min(a["model_m_eq"] for a in base) / min(a["model_m_eq"] for a in glv))
wbar = 0.36
ax1.bar([x - wbar / 2 for x in xs], counted_r, wbar, color=C_GLV, label="counted (this run)")
ax1.bar([x + wbar / 2 for x in xs], model_r, wbar, color=C_GLV2, label="frozen model")
for x, c, m in zip(xs, counted_r, model_r):
    ax1.text(x - wbar / 2, c + 0.004, f"{c:.3f}", ha="center", va="bottom", fontsize=7)
    ax1.text(x + wbar / 2, m + 0.004, f"{m:.3f}", ha="center", va="bottom", fontsize=7)
    ax1.plot([x - 0.45, x + 0.45], [0.99 * m] * 2, color=INK, lw=0.8, ls=":")
ax1.set_xticks(list(xs), [lbl for _, lbl in CURVES])
ax1.set_ylim(1.40, 1.58)
ax1.set_ylabel("best baseline / best GLV-2 (counted M_eq)")
ax1.set_title("(a) counted vs model ratio; dotted: 99 % of model (cond. 3)", fontsize=8.5, loc="left")
ax1.legend(loc="upper center", fontsize=7, ncol=2)

for i, (curve, _) in enumerate(CURVES):
    arms = [a for a in DATA[curve][0]["arms"] if a["key"].split("/")[0] in ("baseline", "glv", "aa")]
    for n, a in enumerate(arms):
        dev = 100.0 * (a["m_eq"] - a["model_m_eq"]) / a["model_m_eq"]
        kind = a["key"].split("/")[0]
        col = C_GLV if kind == "glv" else (C_AA if kind == "aa" else C_BASE)
        jitter = ((n * 0.618034) % 1.0 - 0.5) * 0.5  # deterministic spread
        ax2.scatter(i + jitter, dev, s=9, color=col, alpha=0.75, linewidths=0)
ax2.axhspan(-1, 1, color=GRID, alpha=0.35, lw=0)
ax2.axhline(0, color=MUTED, lw=0.8)
ax2.set_xticks(list(xs), [f"{lbl}\n{sum(1 for a in DATA[c][0]['arms'] if a['key'].split('/')[0] in ('baseline','glv','aa'))} total arms"
                          for c, lbl in CURVES])
ax2.set_ylim(-1.1, 1.1)
ax2.set_ylabel("(counted − model) / model, %")
ax2.set_title("(b) every total arm, counted vs model (band: ±1 %, cond. 2)", fontsize=8.5, loc="left")
for col, lbl in ((C_BASE, "baseline"), (C_GLV, "GLV-2"), (C_AA, "A/A")):
    ax2.scatter([], [], s=12, color=col, label=lbl)
ax2.legend(loc="upper right", fontsize=7, ncol=3)
fig.tight_layout(rect=(0, 0.06, 1, 1))
caption(fig, "Counting runs (" + ", ".join(f"{c}.counts.json" for c, _ in CURVES)
        + "), opcount build at commit 33e2dd10b; counted M_eq = M + S, mean over 256 pairs; model averages 400 other scalars")
save(fig, "fig_counted_vs_model")

# -- fig_timed ----------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.0, 3.7), gridspec_kw={"width_ratios": [1.2, 1]})
labels = ["counted", "median", "min"]
cols = [C_GLV2, C_GLV, "#003f63"]
wbar = 0.25
aa_text = []
for i, (curve, lbl) in enumerate(CURVES):
    t = by_key(DATA[curve][1]["arms"])
    r = t[declared(DATA[curve][1])]["ratio_reference_over_arm"]
    aa = t[aa_key(DATA[curve][1])]["ratio_reference_over_arm"]["median"]
    floor = abs(1 - aa)
    for j, k in enumerate(labels):
        x = i + (j - 1) * wbar
        ax1.bar(x, r[k], wbar, color=cols[j], label=k if i == 0 else None)
        ax1.text(x, r[k] + 0.004, f"{r[k]:.3f}", ha="center", va="bottom", fontsize=6.5)
    rc = r["counted"]
    ax1.fill_between([i - 0.42, i + 0.42], rc * (1 - floor), rc * (1 + floor), facecolor="none",
                     edgecolor=MUTED, hatch="///", lw=0.5, zorder=3,
                     label="counted ± A/A floor" if i == 0 else None)
    aa_text.append(f"{lbl}\n(A/A median {aa:.3f})")
ax1.set_xticks(range(len(CURVES)), aa_text, fontsize=7.5)
ax1.set_ylim(1.40, 1.63)
ax1.set_ylabel("reference / declared GLV arm")
ax1.set_title("(a) declared GLV (opt/jac/w5) vs the reference baseline", fontsize=8.5, loc="left")
ax1.legend(loc="upper left", fontsize=6.5, ncol=2)

pos = 0
ticks, ticklabels = [], []
for i, (curve, lbl) in enumerate(CURVES):
    t = by_key(DATA[curve][1]["arms"])
    ref = t[reference(DATA[curve][1])]
    med = ref["median_ns_per_op"]
    for key, col, short in ((reference(DATA[curve][1]), C_BASE, "ref"), (aa_key(DATA[curve][1]), C_AA, "A/A"),
                            (declared(DATA[curve][1]), C_GLV, "GLV")):
        vals = [v / med for v in t[key]["ns_per_op_per_round"]]
        bp = ax2.boxplot([vals], positions=[pos], widths=0.6, patch_artist=True, showfliers=True,
                         flierprops={"markersize": 2, "markeredgecolor": col},
                         medianprops={"color": INK, "lw": 1})
        for b in bp["boxes"]:
            b.set_facecolor(col)
            b.set_alpha(0.55)
        ticks.append(pos)
        ticklabels.append(short)
        pos += 1
    ax2.text(pos - 2, 1.62, lbl, ha="center", fontsize=7, color=INK)
    pos += 0.8
ax2.set_xticks(ticks, ticklabels, fontsize=7)
ax2.set_ylim(0.5, 1.7)
ax2.set_ylabel("round ns/op ÷ reference median")
ax2.set_title("(b) per-round wall times, 20 rounds per arm", fontsize=8.5, loc="left")
fig.tight_layout(rect=(0, 0.06, 1, 1))
caption(fig, CAPTION + "; ratios = recorded ratio_reference_over_arm against the reference baseline (Bandersnatch jacobian/w5, "
        "Tom-384 and CP6-782 jacobian/w6); CP6-782's fastest baseline by median is jacobian/w5 (README condition 4: 1.570)")
save(fig, "fig_timed")

# -- fig_grid -----------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(9.2, 5.4), sharex=True)
series = [("glv/{chain}/optimised/jacobian", C_GLV, "-", "GLV-2 opt, Jacobian table"),
          ("glv/{chain}/optimised/affine", C_GLV2, "--", "GLV-2 opt, affine table"),
          ("baseline/jacobian", C_BASE, "-", "baseline, Jacobian table"),
          ("baseline/affine", C_BASE2, "--", "baseline, affine table")]
for i, (curve, lbl) in enumerate(CURVES):
    counts, timing, _ = DATA[curve]
    c, t = by_key(counts["arms"]), by_key(timing["arms"])
    chain = timing["inputs"]["grid_chain"]
    widths = timing["inputs"]["widths"]
    for pat, col, ls, name in series:
        keys = [pat.format(chain=chain) + f"/w{w}" for w in widths]
        axes[0, i].plot(widths, [c[k]["m_eq"] for k in keys], color=col, ls=ls, marker="o", ms=3, label=name)
        axes[1, i].plot(widths, [t[k]["median_ns_per_op"] / 1000 for k in keys], color=col, ls=ls, marker="o", ms=3)
    axes[0, i].set_title(f"{lbl}, chain {chain}", fontsize=8.5)
    axes[1, i].set_xlabel("NAF window width w")
    axes[1, i].set_xticks(widths)
axes[0, 0].set_ylabel("counted M_eq per k·P (mean, 256 pairs)")
axes[1, 0].set_ylabel("median µs per k·P (20 rounds)")
handles, names = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, names, loc="lower center", ncol=4, fontsize=7.5, bbox_to_anchor=(0.5, 0.025))
fig.tight_layout(rect=(0, 0.075, 1, 1))
caption(fig, CAPTION + "; top: counting runs, bottom: timing runs")
save(fig, "fig_grid")

# -- fig_phi ------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.3))
for i, (curve, lbl) in enumerate(CURVES):
    counts, timing, _ = DATA[curve]
    c, t = by_key(counts["arms"]), by_key(timing["arms"])
    for ev, col, mk in (("generic", C_BASE, "s"), ("optimised", C_GLV, "o")):
        keys = [k for k in t if k.startswith("phi/") and k.endswith("/" + ev)]
        axes[i].scatter([c[k]["m_eq"] for k in keys], [t[k]["median_ns_per_op"] / 1000 for k in keys],
                        s=12, color=col, marker=mk, alpha=0.8, linewidths=0, label=f"{ev} evaluator ({len(keys)})")
    axes[i].set_title(f"{lbl}: φ(P) per chain", fontsize=8.5)
    axes[i].set_xlabel("counted M_eq per φ(P) call")
    axes[i].legend(fontsize=6.5, loc="upper left")
axes[0].set_ylabel("median µs per φ(P) call (20 rounds)")
fig.tight_layout(rect=(0, 0.07, 1, 1))
caption(fig, CAPTION + "; stage diagnostics, not a speed of k·P; counts identical on every input (min = max)")
save(fig, "fig_phi")

print("wrote", ", ".join(sorted(p.name for p in here.glob("fig_*.svg"))))
for (curve, _), cr, mr in zip(CURVES, counted_r, model_r):
    print(f"{curve}: counted ratio {cr:.3f}, model ratio {mr:.3f}")

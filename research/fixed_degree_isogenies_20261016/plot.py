"""Lay out the recorded values in data/records.jsonl (written by the native harness
`isogeny_algos/examples/fixed_degree_report.rs`). Plotting only: no value is computed here
beyond selecting and arranging the records.

    python3 plot.py            # writes figures/*.svg, *.pdf and *.png (150 dpi)
"""
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)
RECS = [json.loads(l) for l in open(os.path.join(HERE, "data", "records.jsonl")) if l.strip()]
RUN = open(os.path.join(HERE, "data", "run_id.txt")).read().strip()

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "svg.fonttype": "none"})
C = {"Cornacchia": "#2f6db3", "CoppersmithBivariate": "#d97b29",
     "CoppersmithTrivariate": "#4a9a5b", "naive": "#888888"}


def save(fig, name):
    for ext in ("svg", "pdf"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, f"{name}.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


# 1. pipeline: per-step wall time per case
pl = [r for r in RECS if r["kind"] == "pipeline"]
fig, ax = plt.subplots(figsize=(7.2, 3.4))
steps = [("end_ring_ms", "1 End(E1), End(E2)", "#9bb7d4"), ("connecting_ms", "2-3 connecting ideal, norm form", "#555555"),
         ("solve_ms", "4-5 norm equation, J", "#d97b29"), ("isogeny_ms", "6 ideal -> isogeny (+ certificate)", "#4a9a5b")]
labels = []
for i, r in enumerate(pl):
    bottom = 0.0
    if not r.get("ok"):
        ax.text(i, 30, "no solution\nin box", rotation=90, ha="center", va="bottom", fontsize=6.5)
        s = {"Cornacchia": "Cor", "CoppersmithBivariate": "Coron"}.get(r["solver"], r["solver"])
        labels.append(f"{r['log_p_d']:.2f}\n{s}\n{'planted' if r['planted'] else 'random'}")
        continue
    for key, lab, col in steps:
        v = r[key]
        ax.bar(i, v, bottom=bottom, color=col, label=lab if i == 0 else None, width=0.7)
        bottom += v
    mark = "K" if r["certificate"] == "kani2_certified" else ("p" if r.get("pairing_ok") else "")
    ax.text(i, bottom * 1.02, mark, ha="center", va="bottom", fontsize=7)
    s = {"Cornacchia": "Cor", "CoppersmithBivariate": "Coron"}.get(r["solver"], r["solver"])
    labels.append(f"{r['log_p_d']:.2f}\n{s}\n{'planted' if r['planted'] else 'random'}")
ax.set_xticks(range(len(pl)))
ax.set_xticklabels(labels, fontsize=6.5)
ax.set_ylabel("wall time per step (ms)")
ax.set_xlabel("instance: log_p d / solver / target")
ax.set_title(f"End-to-end fixed-degree isogenies, p = 1605631 (K = Kani-certified, p = pairing only)\nrun {RUN}, one run per instance, release build", fontsize=8)
ax.legend(fontsize=7, frameon=False, loc="upper left")
save(fig, "fig1_pipeline_steps")

# 2. steps 2-5 at cryptographic size
st = [r for r in RECS if r["kind"] == "steps25"]
fig, ax = plt.subplots(figsize=(6.4, 3.4))
for (bits, solver), mk in [((61, "Cornacchia"), "o"), ((61, "CoppersmithBivariate"), "s"),
                           ((61, "CoppersmithTrivariate"), "^"), ((127, "CoppersmithBivariate"), "D")]:
    xs = [r["log_p_d"] for r in st if r["p_bits"] == bits and r["solver"] == solver and r["ok"]]
    ys = [r["ms"] for r in st if r["p_bits"] == bits and r["solver"] == solver and r["ok"]]
    fx = [r["log_p_d"] for r in st if r["p_bits"] == bits and r["solver"] == solver and not r["ok"]]
    if xs:
        ax.scatter(xs, ys, marker=mk, color=C[solver], label=f"{solver}, {bits}-bit p", s=22,
                   facecolors="none" if bits == 127 else C[solver])
    for x in fx:
        ax.axvline(x, color=C[solver], ls=":", lw=0.8)
ax.set_yscale("log")
ax.set_xlabel("log_p d   (d = p^(1/2 + eps), planted)")
ax.set_ylabel("steps 2-5 wall time (ms)")
ax.set_title(f"Norm-equation solvers inside steps 2-5 (KV ideal, norm form, solve, J)\nrun {RUN}; 2 instances per setting; dotted = no solution found in the box", fontsize=8)
ax.legend(fontsize=7, frameon=False)
save(fig, "fig2_steps25_solvers")

# 3. Strassen
sr = [r for r in RECS if r["kind"] == "strassen"]
fig, ax = plt.subplots(figsize=(4.8, 3.2))
ax.loglog([r["l"] for r in sr], [r["naive_ms"] for r in sr], "o-", color=C["naive"], label="naive, l - 1 products")
ax.loglog([r["l"] for r in sr], [r["strassen_ms"] for r in sr], "s-", color=C["Cornacchia"], label="Strassen (Sec. 2.3), Karatsuba trees")
ax.set_xlabel("l")
ax.set_ylabel("l! mod (2^61 - 1), ms (mean of reps)")
ax.set_title(f"Strassen's modular factorial, run {RUN}", fontsize=8)
ax.legend(fontsize=7, frameon=False)
save(fig, "fig3_strassen")

# 4. baselines
bl = [r for r in RECS if r["kind"] == "baseline"]
fig, ax = plt.subplots(figsize=(4.8, 3.2))
es = sorted({r["e"] for r in bl})
ax.semilogy([r["e"] for r in bl], [r["mitm_leaves"] for r in bl], "o", color=C["Cornacchia"], label="MITM leaves (measured)")
ax.semilogy(es, [2 * 3 * 2 ** (e // 2 - 1) for e in es], "-", color=C["Cornacchia"], lw=0.8, label="MITM 2 * 3 * 2^(e/2-1)")
ax.semilogy([r["e"] for r in bl], [r["vow_f_evals"] for r in bl], "s", color=C["CoppersmithBivariate"], label="vOW f-evaluations (measured, w = 16)")
vexp = {r["e"]: r["vow_expected"] for r in bl}
ax.semilogy(es, [vexp[e] for e in es], "-", color=C["CoppersmithBivariate"], lw=0.8, label="vOW 2.5 sqrt(|S|^3 / w)")
ax.set_xlabel("e  (d = 2^e)")
ax.set_ylabel("count")
ax.set_title(f"Baselines for d = 2^e, p = {bl[0]['p'] if bl else '?'}; run {RUN}", fontsize=8)
ax.legend(fontsize=6.5, frameon=False)
save(fig, "fig4_baselines")

# 5. diagram of the pipeline (what runs where)
fig, ax = plt.subplots(figsize=(8.2, 2.8))
ax.axis("off")
boxes = ["E1, E2, d", "1  End(E1), End(E2)\nBFS on the Phi_2\ngraph + sub-ideal\nwalk",
         "2  I = N O1 O2\n(Kirschmer-\nVoight)", "3  Q = 2 Nrd/N(I)\non an LLL basis",
         "4  Q(x) = 2d:\nCornacchia |\nCoron | trivariate", "5  J = I conj(a)/N(I)\nN(J) = d",
         "6  kernel (smooth d)\nor torsion images\n+ Kani certificate"]
W, STEP = 0.122, 0.1425
for i, text in enumerate(boxes):
    x = 0.004 + i * STEP
    ax.add_patch(FancyBboxPatch((x, 0.33), W, 0.44, boxstyle="round,pad=0.008", fc="#eef3f9", ec="#2f6db3"))
    ax.text(x + W / 2, 0.55, text, ha="center", va="center", fontsize=6.2)
    if i:
        ax.annotate("", xy=(x - 0.002, 0.55), xytext=(x - STEP + W + 0.012, 0.55), arrowprops=dict(arrowstyle="->", lw=0.8))
ax.text(0.15, 0.2, "steps 1 and 6 use the Deuring correspondence on E0 over F_(p^4): p < 2^31", fontsize=6.5, color="#4a9a5b")
ax.text(0.29, 0.1, "steps 2-5 are pure quaternion arithmetic, any p (tested at 61 and 127 bits)", fontsize=6.5, color="#d97b29")
ax.text(0.15, 0.88, "hybrid (smooth d = l^e): guess an l^e1 prefix as sub-ideals after step 1, run 2-5 on the rest, compose before 6", fontsize=6.5)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
save(fig, "fig5_pipeline_diagram")
print("figures written to", OUT)

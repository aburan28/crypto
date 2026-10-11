"""Plot recorded ratios from the r4_16_4 JSON outputs (no new computation
beyond reading fields and taking reciprocals for orientation)."""
import json
import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt

def rows(name):
    return json.load(open(name))["rows"]

cen = rows("census_torsion.json")
pt = rows("pipeline_torsion.json")
pp = rows("pipeline_prime.json")
sp = rows("stream_prime.json")
sw = rows("sweep_torsion_10.json") + rows("sweep_torsion_12.json")

yield_cats = [
    ("distinct rows\nm=2, (c)/(b)", [r["size_matched_distinct_rows_c_over_b"][0] for r in cen], 1 / 3),
    ("distinct rows\nm=3, (c)/(b)", [r["size_matched_distinct_rows_c_over_b"][1] for r in cen], 1 / 9),
    ("rows per unknown\nm=2, (c)/(b)", [r["size_matched_rows_per_unknown_c_over_b"][0] for r in cen], 1.0),
    ("rows per unknown\nm=3, (c)/(b)", [r["size_matched_rows_per_unknown_c_over_b"][1] for r in cen], 1 / 3),
]
cost_cats = [
    ("U_b/U_a\n(all curves)", [r["U_b"] / r["U_a"] for r in pp + pt], 1 / 3),
    ("total (b)/(a)\nm=2 line oracle", [r["total_b_over_a"] for r in pp + pt], 1 / 3),
    ("trials to full rank\n(b)/(a), prime family", [r["ab_stream"]["b"]["square_trials"] / r["ab_stream"]["a"]["square_trials"] for r in sp], 1 / 3),
    ("total (c)/(b)\nm=2 line oracle", [r["size_matched_total_c_over_b"] for r in pt], 1.0),
    ("total (c)/(b) m=3\nquotient sweep", [r["size_matched_total_c_quotient_over_b"] for r in sw], 1.0),
    ("total (c)/(b) m=3\nsweep on E (derived x3)", [3 * r["size_matched_total_c_quotient_over_b"] for r in sw], 3.0),
]

fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw={"width_ratios": [4, 6]})
for ax, cats, title in [
    (axes[0], yield_cats, "Yield (census, exact counts, 16 torsion curves, 1000 targets each)"),
    (axes[1], cost_cats, "Cost ratios (F_p multiplications; columns; relations)"),
]:
    for i, (label, vals, law) in enumerate(cats):
        xs = [i + (k - (len(vals) - 1) / 2) * 0.03 for k in range(len(vals))]
        ax.scatter(xs, vals, s=14, color="#2f6f9f", alpha=0.75, zorder=3,
                   label="measured, one curve" if i == 0 else None)
        ax.hlines(law, i - 0.35, i + 0.35, colors="#c0392b", linewidth=2, zorder=2,
                  label="Prop. 15.2 prediction" if i == 0 else None)
    ax.set_yscale("log")
    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels([c[0] for c in cats], fontsize=7.5)
    ax.set_title(title, fontsize=9)
    ax.grid(axis="y", which="both", alpha=0.3)
    ax.legend(fontsize=7.5, loc="upper left")
axes[0].set_ylabel("ratio (size-matched to |F_c| = |F|)")
axes[1].set_ylabel("ratio")
fig.suptitle("Lead 16.4: j = 0 curves over F_{p^2}, base x in F_p; (a) <-1>, (b) <-1, zeta>, (c) T-saturated phi^{-1}(x' in F_p)", fontsize=10)
fig.tight_layout()
fig.savefig("ratios.svg")

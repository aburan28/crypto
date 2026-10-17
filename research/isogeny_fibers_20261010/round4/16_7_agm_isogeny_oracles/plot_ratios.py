"""Plot the class-factor ratios recorded by r4_16_7_agm_pushforward (check 5).

Reads out_seed7.json and draws recorded values only; no computation beyond
selecting fields.  Run with the shared venv python.
"""
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

here = pathlib.Path(__file__).resolve().parent
d = json.loads((here / "out_seed7.json").read_text())
rows = {r["m"]: r for r in d["check5_ratios_derived"]}
labels = [
    ("image_law_phi_plain_over_plain_E", "image_law_predicted", "phi(F) on E' / F on E\n(image law, i = 2)"),
    ("image_law_phi_sat_over_plain_E2", "image_law_predicted", "phi(F_sat) / plain' on E'\n(image law, i = 2)"),
    ("saturation_law_sat_over_plain_E", "saturation_law_predicted", "F_sat / F on E\n(saturation, 1/2^(m-1))"),
]
fig, ax = plt.subplots(figsize=(8.2, 4.4))
w = 0.18
colors = {2: "#2f5d9b", 3: "#b7791f"}
for gi, (key, pkey, lab) in enumerate(labels):
    for k, m in enumerate((2, 3)):
        x = gi + (k - 0.5) * 2 * w
        v = rows[m][key]
        pred = rows[m][pkey]
        ax.bar(x, v, width=w * 1.8, color=colors[m], label=f"m = {m} (measured ratio)" if gi == 0 else None)
        ax.hlines(pred, x - w, x + w, colors="black", linestyles="dashed",
                  label="heuristic prediction" if (gi == 0 and m == 2) else None)
        ax.text(x, v + 0.04, f"{v:.3f}", ha="center", fontsize=9)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels([lab for _, _, lab in labels], fontsize=9)
ax.set_ylabel("ratio of class factors c = D_m(G) m! / |F|^m")
ax.set_ylim(0, 2.4)
c = d["curve_E"]
ax.set_title(
    f"Lead 16.7 check 5: class-factor ratios on the toy pair E, E' = E/<(0,0)>\n"
    f"p = {c['p']}, N = 4 r, r = {c['r']}; exact census of all m-sums, seed {d['seed']}; one curve pair, no error bars",
    fontsize=9,
)
ax.legend(fontsize=8, loc="upper right")
fig.tight_layout()
fig.savefig(here / "class_factor_ratios.svg")

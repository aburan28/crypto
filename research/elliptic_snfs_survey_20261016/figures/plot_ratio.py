# Layout only: plots recorded values from ../data/ratio_to_rho.tsv (copied
# from the 2026-10-15 frontier campaign). Computes nothing beyond min/max of
# the three recorded band columns for the error bars.
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
rows = []
with open(os.path.join(here, "..", "data", "ratio_to_rho.tsv")) as fh:
    lines = [l for l in fh if not l.startswith("#")]
for r in csv.DictReader(lines, delimiter="\t"):
    rows.append(r)

style = {
    "ic-pkm": dict(color="#c0392b", marker="o",
                   label="ic-pkm: PKM Kummer tower, 2-adicity >= 12 primes (\"special\" p)"),
    "ic-pair": dict(color="#2c6fbb", marker="s",
                    label="ic-pair: pair-table IC, generic primes"),
}

fig, ax = plt.subplots(figsize=(7.2, 4.4))
for arm, st in style.items():
    pts = sorted((float(r["log2_r"]), r) for r in rows if r["arm"] == arm)
    xs = [x for x, _ in pts]
    w0 = [float(r["ratio_w0"]) for _, r in pts]
    band = [[float(r[k]) for k in ("ratio_w0", "ratio_w0.1", "ratio_w1")] for _, r in pts]
    lo = [w - min(b) for w, b in zip(w0, band)]
    hi = [max(b) - w for w, b in zip(w0, band)]
    ax.errorbar(xs, w0, yerr=[lo, hi], fmt=st["marker"], color=st["color"],
                capsize=3, lw=1.4, ms=5, label=st["label"])
    for x, (_, r), y in zip(xs, pts, w0):
        # label each point with its curve's field size (the "fpNN" in its id)
        ax.annotate(r["curve"].split("-")[1], (x, y), textcoords="offset points",
                    xytext=(6, -3), fontsize=6.5, color=st["color"])

ax.axhline(1.0, color="#333", ls="--", lw=1)
ax.text(18.3, 1.15, "matched rho (rho-neg) = 1", fontsize=8, color="#333")
ax.set_yscale("log")
ax.set_ylim(0.7, 400)
ax.set_xlim(18, 41)
ax.set_xlabel("log2 r (prime subgroup order, bits)")
ax.set_ylabel("mean S(IC arm) / mean S(rho-neg), same curve")
ax.set_title("Prime-field index calculus vs matched rho (L0 operation counts)", fontsize=10)
ax.legend(fontsize=7.5, loc="upper left", frameon=False)
ax.grid(True, which="both", alpha=0.25)
fig.text(0.01, 0.01,
         "Points: w = 0 GAE/native op; bars: range over w in {0, 0.1, 1}. n = 32 verified solves per arm per curve "
         "(16 targets x 2 rounds).\nSource: research/ic_frontier_campaign_20261015/analysis/ratios_to_rho.tsv "
         "(commit 5b2ef3ebc). Arms ran on different curve ladders: not a paired comparison.",
         fontsize=6.3, color="#444")
fig.tight_layout(rect=(0, 0.07, 1, 1))
fig.savefig(os.path.join(here, "ratio_to_rho.svg"))
fig.savefig(os.path.join(here, "ratio_to_rho.png"), dpi=150)

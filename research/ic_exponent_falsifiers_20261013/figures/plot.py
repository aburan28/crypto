# Lays out the C4 falsifier's recorded output (raw/c4_exponent_sum_p20_t8.json,
# written by examples/kummer_exponent_bias.rs) as one figure: largest
# residue-class share of the exponent sum per k, curve against the null and
# the positive control.  Recorded values only.
import json, os, sys
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
raw = os.path.join(here, "..", "raw", "c4_exponent_sum_p20_t8.json")
d = json.load(open(raw))
C = {"blue": "#2a78d6", "orange": "#eb6834", "ink": "#0b0b0b", "ink2": "#52514e",
     "grid": "#e6e5e1", "surface": "#fcfcfb", "wash": "#9aa3ad", "aqua": "#1baf7a"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": C["grid"], "axes.labelcolor": C["ink2"],
                     "xtick.color": C["ink2"], "ytick.color": C["ink2"], "axes.titlecolor": C["ink"],
                     "figure.facecolor": C["surface"], "axes.facecolor": C["surface"]})
k = [r["k"] for r in d["rows"]]
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.plot(k, [r["positive_control"] for r in d["rows"]], "o-", color=C["aqua"], label="positive control (exponent sum fixed)")
ax.plot(k, [r["null_same_size"] for r in d["rows"]], "s--", color=C["wash"], label="null: random triples, same count")
ax.plot(k, [r["curve"] for r in d["rows"]], "^-", color=C["blue"], label="curve decompositions")
ax.plot(k, [r["uniform"] for r in d["rows"]], ":", color=C["ink2"], label="1/2^k")
ax.set_xlabel("k (residue of e_a + e_b + e_c mod 2^k)")
ax.set_ylabel("mean largest class share")
ax.set_ylim(0, 1.05)
ax.grid(color=C["grid"], lw=0.6)
ax.set_title(f"p = {d['p']}, t = {d['t']}: {d['targets']} targets, {d['decompositions']} decompositions", fontsize=9)
ax.legend(frameon=False, fontsize=8, loc="center right")
fig.tight_layout()
for ext, kw in (("svg", {}), ("png", {"dpi": 150})):
    fig.savefig(os.path.join(here, f"c4_class_share.{ext}"), **kw)

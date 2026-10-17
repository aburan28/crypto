"""Lay out the recorded census (census/summary.json, written by classify.gp) as
one lane per srand48 run. Plotting only: every value comes from the file or
from the published seed table; nothing is computed here."""
import json, pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).parent
runs = json.loads((HERE / "census" / "summary.json").read_text())
PUBLISHED = {  # srand48 value -> {candidate index: curve name}; FIPS 186-2 / X9.62 seeds
    "863711343": {2185: "prime239v1", 2320: "prime239v2", 3048: "prime239v3"},
    "863966734": {139: "P-192", 5228: "prime192v3", 6013: "prime192v2"},
    "871002724": {1136: "P-521"}, "871696354": {1842: "P-256"},
    "872952608": {2688: "P-384"}, "891607927": {7026: "P-224"},
}
DATE = {"863711343": "1997-05-15", "863966734": "1997-05-18", "871002724": "1997-08-08",
        "871696354": "1997-08-16", "872952608": "1997-08-30", "891607927": "1998-04-03"}

fig, ax = plt.subplots(figsize=(11, 4.6))
for row, r in enumerate(runs):
    ts = r["run"][8:-1]; y = len(runs) - row
    pub = PUBLISHED[ts]; last = max(r["accepted"])
    ax.hlines(y, 0, last, color="#c9ccd3", lw=1, zorder=0)
    tw = [i for i in r["only_twist_prime"] if i <= last]
    bl = [i for i in r["n_longer_than_p"] if i <= last]
    ax.scatter(tw, [y + 0.18] * len(tw), marker="|", s=40, color="#9aa3b5",
               label="prime only for the quadratic twist (the other root -b when p = 3 mod 4)" if row == 0 else None)
    ax.scatter(bl, [y - 0.18] * len(bl), marker="x", s=22, color="#b0702a",
               label="prime order, but bitlen(n) > bitlen(p)" if row == 0 else None)
    skp = [i for i in r["accepted"] if i not in pub]; pp = [i for i in r["accepted"] if i in pub]
    ax.scatter(skp, [y] * len(skp), s=70, facecolors="white", edgecolors="#c0392b", linewidths=1.6, zorder=3,
               label="accepted, not published (skipped)" if row == 0 else None)
    ax.scatter(pp, [y] * len(pp), s=70, color="#1f4e9c", zorder=4, label="published" if row == 0 else None)
    for i in pp:
        dx, ha = {"prime239v1": (-4, "right"), "prime239v2": (4, "left"), "prime239v3": (6, "left")}.get(pub[i], (0, "center"))
        ax.annotate(pub[i], (i, y), xytext=(dx, 9), textcoords="offset points", ha=ha, fontsize=8, color="#1f4e9c")
ax.set_yticks(range(len(runs), 0, -1))
ax.set_yticklabels([f"srand48({r['run'][8:-1]})\n{DATE[r['run'][8:-1]]}  p{r['p_bits']}" for r in runs], fontsize=8)
ax.set_xlabel("candidate index (consecutive 20-byte blocks of drand48 output)")
ax.set_ylim(0.4, len(runs) + 0.7); ax.set_xlim(-100, 7300)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=8, frameon=False)
ax.set_title("Every prime-order candidate up to each run's last published seed "
             f"({sum(r['candidates'] for r in runs):,} point-counted candidates, PARI/GP 2.15.4 ellsea)", fontsize=9, pad=36)
fig.tight_layout()
fig.savefig(HERE / "streams.svg"); fig.savefig(HERE / "streams.png", dpi=150)

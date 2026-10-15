"""Plot the Phase 1 frontier from `ecbench band --json` output (layout only).

Every number comes from the native band report; this script draws it.
Marker = mean S at w = 0.1; vertical bar = the band (w = 0 to w = 1).
Open marker: some run of that arm on that curve did not verify (incomplete).

    python3 plot_frontier.py analysis/band-koblitz.json analysis/band-prime.json
    # writes analysis/frontier.svg and analysis/frontier.png
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PAL = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#999999", "#000000"]
ORDER = ["rho-strong", "rho-neg", "rho-frob", "bsgs-neg", "kangaroo", "ic-compact", "ic-negfold",
         "ic-f4", "ic-shared", "ic-pair", "ic-pkm"]

paths = sys.argv[1:]
fig, axes = plt.subplots(1, len(paths), figsize=(5.2 * len(paths), 4.0), squeeze=False)
for ax, path in zip(axes[0], paths):
    rep = json.load(open(path))
    w = rep["weights"]
    i_mid, i_lo, i_hi = w.index(0.1), w.index(0.0), w.index(1.0)
    arms = sorted({r["arm"] for r in rep["rows"] if not r["arm"].endswith("-aa")},
                  key=lambda a: ORDER.index(a) if a in ORDER else 99)
    for k, arm in enumerate(arms):
        rows = sorted((r for r in rep["rows"] if r["arm"] == arm), key=lambda r: r["log2_r"])
        xs, ys, lo, hi, full = [], [], [], [], []
        for r in rows:
            m = r["mean_s"]
            if m[i_mid] is None:
                continue
            xs.append(r["log2_r"]); ys.append(m[i_mid])
            lo.append(m[i_mid] - (m[i_lo] if m[i_lo] is not None else m[i_mid]))
            hi.append((m[i_hi] if m[i_hi] is not None else m[i_mid]) - m[i_mid])
            full.append(r["verified"] == r["runs"])
        c = PAL[k % len(PAL)]
        ax.errorbar(xs, ys, yerr=[lo, hi], color=c, lw=1.0, capsize=2, fmt="-", alpha=0.8, label=arm)
        for x, y, f in zip(xs, ys, full):
            ax.plot(x, y, "o", ms=4.5, color=c, mfc=c if f else "white")
    ax.set_yscale("log")
    ax.set_xlabel("log2 r")
    ax.set_ylabel("mean S priced at w = 0.1 (bar: w = 0 to 1)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7, frameon=False, ncol=2)
    ax.set_title(os.path.basename(path).replace(".json", ""), fontsize=9)
fig.suptitle("IC frontier campaign, Phase 1: one public target per run, 16 targets x 2 rounds per curve,"
             " counts (L0), AGENTS.md §5a band", fontsize=8.5)
fig.tight_layout()
out = os.path.join(os.path.dirname(paths[0]), "frontier")
fig.savefig(out + ".svg")
fig.savefig(out + ".png", dpi=150)

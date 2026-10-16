#!/usr/bin/env python3
"""Render the frozen ladder results (l=9, l=10) of the binary MITM oracle.

Plots recorded per-target medians only. Outputs SVG and 150-dpi PNG.
"""
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path("research/binary_mitm_oracle_l9_l10_20261015")
OUT = ROOT / "reports" / "run-20261015T000000Z"


def med(v):
    s = sorted(v)
    return s[len(s) // 2]


fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2), sharey=False)
rungs = [(9, "l9n27"), (10, "l10n31")]
for ax, (l, sub) in zip(axes, rungs):
    rows = json.loads((ROOT / sub / "results.json").read_text())["decomposition"]["per_target"]
    mitm = sorted(r["mitm_ms"] for r in rows)
    pairs = sorted(r["pairs_ms"] for r in rows)
    ax.plot(range(1, len(mitm) + 1), mitm, "o-", ms=3, lw=1,
           label=f"MITM, median {med(mitm):.2f} ms")
    ax.plot(range(1, len(pairs) + 1), pairs, "s-", ms=3, lw=1,
           label=f"pairs-and-solve, median {med(pairs):.0f} ms")
    ax.axhline(0.5 * med(pairs), color="gray", ls="--", lw=1,
               label=f"gate: half of pairs median ({0.5 * med(pairs):.0f} ms)")
    ax.set_yscale("log")
    ax.set_xlabel("target rank (sorted per-target wall time, 64 targets)")
    ax.set_title(f"$\\ell={l}$, $n={27 if l == 9 else 31}$")
    ax.legend(fontsize=7)
    ax.grid(True, which="both", alpha=0.3)
axes[0].set_ylabel("per-target wall time (ms, host-contended)")
fig.suptitle("Binary MITM oracle ladder: sub-$2^{2\\ell}$ holds at $\\ell=9$ and $\\ell=10$", fontsize=10)
fig.text(0.5, 0.01,
        "64 targets per rung, 0 disagreements vs pairs-and-solve, max 2,076 / 4,212 counted point-adds "
        "vs 134,912 / 525,312 quartic solves; run binary_mitm_oracle_l9_l10_20261015",
        ha="center", fontsize=7)
fig.tight_layout(rect=(0, 0.04, 1, 1))
OUT.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT / "mitm_ladder.svg")
fig.savefig(OUT / "mitm_ladder.png", dpi=150)
print(f"wrote {OUT / 'mitm_ladder.svg'} and .png")

#!/usr/bin/env python3
"""Render the frozen results.json of the binary MITM oracle beat.

Plots recorded per-target medians only (no new computation beyond sorting
and medians, which the report states). Outputs SVG and 150-dpi PNG.
"""
import json
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULTS = pathlib.Path("research/binary_mitm_oracle_l8_20261015/results.json")
OUT = pathlib.Path("research/binary_mitm_oracle_l8_20261015/reports/run-20261015T000000Z")

data = json.loads(RESULTS.read_text())
rows = data["decomposition"]["per_target"]
mitm = sorted(r["mitm_ms"] for r in rows)
pairs = sorted(r["pairs_ms"] for r in rows)


def med(v):
    return v[len(v) // 2]


fig, ax = plt.subplots(figsize=(7.2, 4.2))
ax.plot(range(1, len(mitm) + 1), mitm, "o-", ms=3, lw=1,
       label=f"MITM pair table, median {med(mitm):.3f} ms")
ax.plot(range(1, len(pairs) + 1), pairs, "s-", ms=3, lw=1,
       label=f"pairs-and-solve, median {med(pairs):.3f} ms")
ax.set_xlabel("target rank (sorted per-target wall time, 64 targets)")
ax.set_ylabel("per-target wall time (ms, host-contended)")
ax.set_title(
    "First binary-regime sub-$2^{2\\ell}$ decomposition oracle, "
    "$\\ell=8$, $n=25$, non-Koblitz $a=z$"
)
ax.set_yscale("log")
ax.axhline(0.5 * med(pairs), color="gray", ls="--", lw=1,
           label=f"gate: half of pairs-and-solve median ({0.5 * med(pairs):.3f} ms)")
ax.legend(fontsize=8)
ax.grid(True, which="both", alpha=0.3)
caption = (
    f"n=64 targets (20 found / 44 refuted), 0 disagreements, max {max(r['mitm_ops'] for r in rows)} "
    f"counted point-adds vs 2^15=32768 pairs; run research/binary_mitm_oracle_l8_20261015"
)
fig.text(0.5, 0.01, caption, ha="center", fontsize=7)
fig.tight_layout(rect=(0, 0.04, 1, 1))
OUT.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT / "mitm_vs_pairs.svg")
fig.savefig(OUT / "mitm_vs_pairs.png", dpi=150)
print(f"wrote {OUT / 'mitm_vs_pairs.svg'} and .png")
print(f"median mitm {med(mitm):.3f} ms, median pairs {med(pairs):.3f} ms, "
      f"ratio {med(pairs) / med(mitm):.1f}x")

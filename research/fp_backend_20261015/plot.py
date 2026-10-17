"""Lay out the recorded fp_backend benchmark values (no computation).

Reads runs/final/summary.csv and runs/final/claims.csv (written by the
Rust programs examples/fp_backend_bench.rs and fp_backend_report.rs) and
draws fig_ns_per_op.{svg,png} and fig_claims.{svg,png}.
"""
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "runs", "final")

INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3df"
BLUE = "#2a78d6"   # this module's backends
GRAY = "#9a9893"   # existing code / references
GOOD = "#0ca30c"
BAD = "#e34948"

plt.rcParams.update({
    "font.size": 8.5, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.titlesize": 9.5,
    "svg.fonttype": "none",
})

summary = {}
with open(os.path.join(RUN, "summary.csv")) as fh:
    for r in csv.DictReader(fh):
        summary[r["arm"]] = r

OURS = ("Fp32", "Fp64", "FpP192", "FpWide<3>", "FpWide<4>", "FpWide<4+spare>")
ROWS = [
    ("f4_fp_tower", "p31"), ("Fp32", "p31"), ("Fp32", "p32"),
    ("prime_fast::Fp64", "p61"), ("Fp64", "p61"), ("Fp64", "p63"),
    ("BigUint", "P-192"), ("ct_bignum", "P-192"), ("FpWide<3>", "P-192"), ("FpP192", "P-192"),
    ("BigUint", "P-224"), ("ct_bignum", "P-224"), ("FpWide<4>", "P-224"), ("FpWide<4+spare>", "P-224"),
    ("BigUint", "P-256"), ("ct_bignum", "P-256"), ("FpWide<4>", "P-256"),
]
PANELS = [("mul_tp", "multiplication (throughput)"), ("inv", "inversion"),
          ("sqrt", "square root (of a square)"), ("axpy_lazy_1e4", "F4 row update, lazy (per entry)")]

fig, axes = plt.subplots(1, 4, figsize=(11, 5.2), sharey=True)
labels = [f"{b} · {p}" for b, p in ROWS]
for ax, (op, title) in zip(axes, PANELS):
    for i, (b, p) in enumerate(ROWS):
        r = summary.get(f"{b}/{p}/{op}")
        if r is None:
            continue
        med, lo, hi = float(r["median_ns"]), float(r["min_ns"]), float(r["max_ns"])
        col = BLUE if b in OURS else GRAY
        ax.barh(i, med, color=col, height=0.62)
        ax.plot([lo, hi], [i, i], color=INK, lw=0.9)
        txt = f"{med:.3g}" if med < 1000 else f"{med/1000:.3g} µs"
        ax.text(med * 1.25, i, txt, va="center", fontsize=7, color=INK2)
    ax.set_xscale("log")
    ax.set_title(title, color=INK)
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_xlabel("ns per operation (log)")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    lo_, hi_ = ax.get_xlim()
    ax.set_xlim(lo_, hi_ * 6)
axes[0].set_yticks(range(len(ROWS)))
axes[0].set_yticklabels(labels)
axes[0].invert_yaxis()
from matplotlib.patches import Patch  # noqa: E402
fig.legend(handles=[Patch(color=BLUE, label="fp_backend (this module)"),
                    Patch(color=GRAY, label="existing code / reference")],
           loc="lower center", ncol=2, frameon=False)
fig.suptitle("fp_backend: median ns per operation, 9 interleaved repetitions "
             "(whiskers: min–max); one core, shared 4-vCPU host", color=INK, fontsize=10)
fig.tight_layout(rect=(0, 0.05, 1, 0.95))
for ext, kw in (("svg", {}), ("png", {"dpi": 150})):
    fig.savefig(os.path.join(HERE, f"fig_ns_per_op.{ext}"), **kw)

claims = list(csv.DictReader(open(os.path.join(RUN, "claims.csv"))))
claims = [c for c in claims if c["verdict"] in ("PASS", "FAIL")]
fig, ax = plt.subplots(figsize=(8.5, 0.24 * len(claims) + 1.2))
for i, c in enumerate(claims):
    med, lo, hi, thr = (float(c[k]) for k in ("ratio_median", "ratio_min", "ratio_max", "threshold"))
    ok = c["verdict"] == "PASS"
    ax.plot([lo, hi], [i, i], color=INK2, lw=1)
    ax.plot([med], [i], marker="o" if ok else "X", ms=6, color=GOOD if ok else BAD,
            mec="white", mew=0.8, zorder=3)
    ax.plot([thr, thr], [i - 0.38, i + 0.38], color=INK, lw=1.4)
ax.set_yticks(range(len(claims)))
ax.set_yticklabels([f"{c['id']}  {c['description'][:70]}" for c in claims], fontsize=7)
ax.invert_yaxis()
ax.set_xscale("log")
ax.axvline(1.0, color=GRID, lw=1, zorder=0)
ax.set_xlabel("paired ratio baseline / candidate (log; > 1 = candidate faster)")
ax.set_title("Pre-registered claims: median paired ratio (● PASS, ✕ FAIL), min–max line, "
             "threshold tick", color=INK, fontsize=9)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.grid(axis="x", color=GRID, lw=0.5)
fig.tight_layout()
for ext, kw in (("svg", {}), ("png", {"dpi": 150})):
    fig.savefig(os.path.join(HERE, f"fig_claims.{ext}"), **kw)
print("ok")

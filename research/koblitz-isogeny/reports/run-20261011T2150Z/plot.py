"""Lay out recorded values (experiments/koblitz_m3, experiments/koblitz_presentation) as SVG.
Plots recorded numbers only; computes no result."""
import json
from pathlib import Path
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parents[4]
here = Path(__file__).parent
rows = lambda p: [json.loads(l) for l in open(root / p)]
diag = rows("experiments/koblitz_m3/geom_cost_diag.jsonl")
pin = rows("experiments/koblitz_m3/pin_vs_e3.jsonl")
full = rows("experiments/koblitz_m3/full_solves_pin.jsonl")
icc = json.load(open(root / "experiments/koblitz_presentation/icc_n19_m3_geom.json"))["n=19 a2=1 l=6 m=3"]

# 1: ms per call by solver variant
fig, ax = plt.subplots(figsize=(7.2, 4.0))
labels, n17, n19 = [], [], []
def val(rs, n, pred):
    v = [r["us_per_call"] / 1e3 for r in rs if r["n"] == n and pred(r)]
    return v[0] if v else float("nan")
variants = [
    ("e3 serial, mono", diag, lambda r: r["wordfield_serial"] and not r["geometric"]),
    ("e3 serial, geom", diag, lambda r: r["wordfield_serial"] and r["geometric"]),
    ("e3 clmul, mono", diag, lambda r: not r["wordfield_serial"] and not r["geometric"]),
    ("e3 clmul, geom", diag, lambda r: not r["wordfield_serial"] and r["geometric"]),
    ("PinX3, geom", pin, lambda r: r["m3_solver"] == "PinX3"),
]
x = range(len(variants))
ax.bar([i - 0.2 for i in x], [val(rs, 17, p) for _, rs, p in variants], 0.4, label="n = 17, l = 5")
ax.bar([i + 0.2 for i in x], [val(rs, 19, p) for _, rs, p in variants], 0.4, label="n = 19, l = 6")
ax.set_xticks(list(x))
ax.set_xticklabels([v[0] for v in variants], rotation=15, fontsize=8)
ax.set_yscale("log")
ax.set_ylabel("ms per call (measured)")
ax.set_title("m = 3, a6 = 1: diag 20 probes/cell; PinX3 40 probes/cell", fontsize=9)
ax.legend(); ax.grid(alpha=0.3, axis="y")
fig.tight_layout(); fig.savefig(here / "per_call.svg")

# 2: ICC with CI, n = 19, m = 3
fig, ax = plt.subplots(figsize=(6.4, 3.4))
names = ["yield_per_Fm", "log_us_per_call", "work_per_call", "log_projected_solve"]
for i, k in enumerate(names):
    p, (lo, hi) = icc[k]["icc_curve"], icc[k]["ci95"]
    ax.errorbar([i], [p], yerr=[[p - lo], [hi - p]], fmt="o", capsize=4, color="C0")
ax.axhline(0.05, ls="--", color="green", lw=1, label="H0 bound 0.05")
ax.axhline(0.10, ls=":", color="red", lw=1, label="H1 threshold 0.10")
ax.set_xticks(range(len(names))); ax.set_xticklabels(names, fontsize=8)
ax.set_ylabel("ICC_curve (95 % bootstrap CI)")
ax.set_title("n = 19, m = 3, l = 6: 32 curves x 4 geometric V, 400 probes/cell", fontsize=9)
ax.legend(fontsize=8); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(here / "icc_n19_m3.svg")

# 3: PinX3 full solves
fig, ax = plt.subplots(figsize=(6.0, 3.2))
for geo, c in [(False, "C0"), (True, "C1")]:
    pts = sorted((r["n"], r["seconds"]) for r in full if r["geometric"] == geo)
    ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", color=c, label="geometric V" if geo else "monomial V")
ax.set_yscale("log"); ax.set_xlabel("n"); ax.set_ylabel("wall seconds per full solve")
ax.set_title("m = 3 full IC solves with PinX3, a6 = 1; all BSGS-verified", fontsize=9)
ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(here / "full_solves_pin.svg")
print("ok")

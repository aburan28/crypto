"""Lay out recorded ICC values (experiments/koblitz_presentation/icc_n31_m3_geom*.json) as SVG.
Plots recorded numbers only; computes no result."""
import json
from pathlib import Path
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt
root = Path(__file__).resolve().parents[4]
here = Path(__file__).parent
d = json.load(open(root / "experiments/koblitz_presentation/icc_n31_m3_geom_diag.json"))["n=31 a2=0 l=10 m=3"]
names = ["yield_per_Fm", "log_us_per_call", "work_per_call", "log_projected_solve",
         "diag_factor_base_points", "diag_pins_per_call", "diag_roots_per_pin", "diag_work_per_call_per_F"]
fig, ax = plt.subplots(figsize=(8.0, 3.8))
for i, k in enumerate(names):
    p, (lo, hi) = d[k]["icc_curve"], d[k]["ci95"]
    ax.errorbar([i], [p], yerr=[[p - lo], [hi - p]], fmt="o", capsize=4, color="C0" if i < 4 else "C2")
ax.axhline(0.05, ls="--", color="green", lw=1, label="H0 bound 0.05")
ax.axhline(0.10, ls=":", color="red", lw=1, label="H1 threshold 0.10")
ax.axvline(3.5, color="grey", lw=0.8)
ax.text(1.5, 0.47, "pre-registered", ha="center", fontsize=8)
ax.text(5.5, 0.47, "post-hoc diagnostics", ha="center", fontsize=8)
ax.set_xticks(range(len(names))); ax.set_xticklabels([n.replace("diag_", "") for n in names], rotation=20, fontsize=7)
ax.set_ylim(-0.02, 0.5)
ax.set_ylabel("ICC_curve (95 % bootstrap CI)")
ax.set_title("n = 31, a2 = 0, m = 3, l = 10, PinX3: 32 curves x 4 geometric V, 2000 probes/cell", fontsize=9)
ax.legend(fontsize=8, loc="upper right"); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(here / "icc_n31_m3.svg")
print("ok")

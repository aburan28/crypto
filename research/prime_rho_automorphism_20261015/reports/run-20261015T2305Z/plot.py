"""Plot the recorded ecbench comparison ratios (no new computation)."""
import json, sys, math
from pathlib import Path
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

here = Path(__file__).resolve().parent
study = here.parent.parent
panels = [("cubic", "rho.cubic / rho.negation  (j = 0, A = 6)", 1 / math.sqrt(3), "1/√3"),
          ("quartic", "rho.quartic / rho.negation  (j = 1728, A = 4)", 1 / math.sqrt(2), "1/√2")]
fig, axes = plt.subplots(2, 1, figsize=(8, 7.5))
for ax, (k, title, pred, plabel) in zip(axes, panels):
    c = json.loads(next((study / "sessions" / k / "comparisons").glob("*.json")).read_text())
    curves = sorted(c["curves"], key=lambda x: x["log2_r"])
    ys = range(len(curves))
    for y, cv in zip(ys, curves):
        lo, hi = cv["ci95"]
        ax.plot([lo, hi], [y, y], color="0.45", lw=1.5)
        ax.plot(cv["ratio_b_over_a"], y, "o", color="#1f5fa8")
    o = c["ops"]
    yp = len(curves) + 0.6
    ax.plot(o["ci95"], [yp, yp], color="black", lw=2.5)
    ax.plot(o["ratio_b_over_a"], yp, "D", color="black")
    ax.axvline(pred, color="#c0392b", ls="--", lw=1.2, label=f"predicted {plabel} = {pred:.3f}")
    ax.axvline(1.0, color="0.6", ls=":", lw=1, label="no change")
    ax.set_yticks(list(ys) + [yp])
    ax.set_ylim(-0.6, len(curves) + 1.4)
    ax.set_yticklabels([f"{cv['slug']}  (log₂ r = {cv['log2_r']:.1f})" for cv in curves] + [f"pooled, {o['pairs']} pairs"], fontsize=7)
    ax.set_xlabel("ratio of mean S (group-addition equivalents), 95% cluster interval")
    ax.set_title(f"{title}\n{c['comparison_id']}", fontsize=9)
    ax.set_xlim(0, 1.6)
    ax.legend(fontsize=7, loc="upper right", framealpha=0.9)
fig.suptitle("Rho on the full automorphism group vs negation-only rho, round-2 E2 crater curves\n"
             "sessions ECBS1h5e1b9912a46d (cubic), ECBS1hf8d1171784f1 (quartic); 4 targets × 3 rounds per curve", fontsize=9)
fig.tight_layout()
fig.savefig(here / "ratios.svg")
fig.savefig(here / "ratios.png", dpi=150)

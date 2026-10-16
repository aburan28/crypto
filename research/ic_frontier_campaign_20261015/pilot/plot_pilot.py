"""Plot the Phase 0 pilot's recorded S per arm against log2 r (recorded values only).

    python3 plot_pilot.py   # writes pilot_s.svg, pilot_s.png
"""
import json, math, glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rows = []
for path in sorted(glob.glob("sessions/*/records.jsonl")):
    for line in open(path):
        r = json.loads(line)
        if r.get("warmup") or r["cost"].get("s") is None:
            continue
        rows.append((r["workload"]["curve"]["family"], math.log2(int(r["workload"]["curve"]["r"])),
                     r["arm"], r["cost"]["s"], r["cost"]["lower_bound"], r["outcome"]["status"]))
fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6))
cols = {}
pal = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#999999", "#000000"]
for ax, fam, title in ((axes[0], "koblitz", "Koblitz E_a over F_2^n"), (axes[1], "prime", "prime a=-3 / PKM curves")):
    arms = sorted({a for f, _, a, *_ in rows if f == fam})
    for i, arm in enumerate(arms):
        pts = [(x, s, lb, st) for f, x, a, s, lb, st in rows if f == fam and a == arm]
        for x, s, lb, st in pts:
            mk = "o" if st == "verified" else "x"
            ax.plot(x, s, mk, ms=5, mfc=("white" if lb and mk == "o" else pal[i % 8]), mec=pal[i % 8], color=pal[i % 8])
        ax.plot([], [], "o", color=pal[i % 8], label=arm)
    ax.set_yscale("log"); ax.set_xlabel("log2 r"); ax.set_ylabel("charged S (single train targets)")
    ax.set_title(title, fontsize=9); ax.grid(alpha=0.25); ax.legend(fontsize=6.5, frameon=False, ncol=2)
fig.suptitle("Phase 0 pilot, 2026-10-15: feasibility probes, 1-2 public train targets per arm, L0 counts.\n"
             "open marker = lower bound (unpriced native work); x = exhausted. Not a measurement.", fontsize=8)
fig.tight_layout()
fig.savefig("pilot_s.svg"); fig.savefig("pilot_s.png", dpi=150)

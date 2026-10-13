"""Plot recorded gpu/nist rho rates (data.csv). Lays out recorded values only."""
import csv, sys, pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else ".")
rows = [r for r in csv.DictReader(l for l in open(here / "data.csv") if not l.startswith("#"))]

fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True)
for ax, key, title in ((axes[0], "P256_billion_it_s", "P-256"), (axes[1], "P384_billion_it_s", "P-384")):
    for th, mk in (("64", "^"), ("128", "s"), ("256", "o")):
        pts = sorted((int(r["W"]), float(r[key])) for r in rows
                     if r["round"] == "1" and r["threads"] == th)
        ax.plot([p[0] for p in pts], [p[1] for p in pts], marker=mk, label=f"round 1, {th} thr/block")
    pts = sorted((int(r["W"]), float(r[key])) for r in rows
                 if r["round"] == "2" and r["maxreg"] == "0" and r["threads"] == "256")
    ax.plot([p[0] for p in pts], [p[1] for p in pts], "o", mfc="white", mec="black", ms=8,
            label="round 2, 256 thr/block (PR 1778)")
    first = [r for r in rows if r["round"] == "first" and r["build"] == "portable"][0]
    ax.axhline(float(first[key]), ls=":", c="gray", label="first working build (portable, 128 thr, W=8)")
    ax.set_xscale("log", base=2)
    ax.set_xticks([2, 4, 8, 16, 32], ["2", "4", "8", "16", "32"])
    ax.set_title(f"{title} Pollard rho, RTX PRO 6000")
    ax.set_xlabel("walks per thread W (one shared inversion)")
    ax.set_ylabel("billion iterations/s (one cell, 8,192 iterations)")
    ax.grid(alpha=0.3)
axes[0].legend(fontsize=7, loc="lower right")
fig.suptitle("gpu/nist rho walk: batch width x launch geometry, PTX build, unrestricted registers", fontsize=10)
fig.tight_layout()
fig.savefig(here / "nist_rho_w_sweep.svg")
fig.savefig(here / "nist_rho_w_sweep.png", dpi=150)

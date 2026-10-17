"""Layout only: plots outcomes recorded in this study's runs and the scale
study's (research/dreg_closure_scale_20261016/runs).  Each cell is placed at
(l, N); filled = the closure refutes at degree 4 on every draw, hollow = it
saturates at degree 4 without refuting on every draw; circle m = 3, square
m = 4; this study's cells are outlined in black."""
import glob, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
studies = {"scale": os.path.join(here, "../../../dreg_closure_scale_20261016/runs"),
           "driver": os.path.join(here, "../../runs")}
cells = {}
for tag, d in studies.items():
    for p in sorted(glob.glob(os.path.join(d, "*.d*.jsonl"))):
        if os.path.basename(p).startswith("control-"):
            continue
        t = open(p).read().strip()
        if not t:
            continue
        r = json.loads(t.splitlines()[0])
        k = (r["m"], r["n"], r["ell"], r["n_vars"], tag)
        cells.setdefault(k, []).append(r["closure_degree"])

fig, ax = plt.subplots(figsize=(8.5, 5.2))
for (m, n, ell, N, tag), ds in sorted(cells.items()):
    refutes = all(d is not None and d <= 4 for d in ds)
    col = "tab:green" if m == 3 else "tab:purple"
    x = ell + (-0.08 if m == 3 else 0.08)
    ax.scatter(x, N, s=150, marker="o" if m == 3 else "s",
               facecolors=col if refutes else "white", edgecolors="black" if tag == "driver" else col,
               linewidths=2.0 if tag == "driver" else 1.5, zorder=3)
    ax.annotate(f"m{m} ({n},{ell}) m·ℓ={m*ell}\n{sum(d is not None and d<=4 for d in ds)}/{len(ds)} refute at 4",
                (x, N), xytext=(10, -4), textcoords="offset points", fontsize=6.5)
for lab, kw in [("m=3, refutes at 4", dict(marker="o", facecolors="tab:green", edgecolors="tab:green")),
                ("m=3, degree 4 saturates without 1", dict(marker="o", facecolors="white", edgecolors="tab:green")),
                ("m=4, refutes at 4", dict(marker="s", facecolors="tab:purple", edgecolors="tab:purple")),
                ("this study's cells (black outline)", dict(marker="o", facecolors="white", edgecolors="black"))]:
    ax.scatter([], [], s=80, label=lab, linewidths=1.5, **kw)
ax.set_xlabel("ℓ (subspace dimension)")
ax.set_ylabel("unknowns N")
ax.set_xlim(2.5, 8.2)
ax.grid(alpha=0.3)
ax.legend(fontsize=7, loc="lower right")
ax.set_title("Where the MutantXL closure stops refuting at degree 4 — 56 draws, two studies\n"
             "research/dreg_closure_driver_20261016 (run 2026-10-16T20:49Z) + dreg_closure_scale_20261016", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(here, "closure_driver.svg"))
fig.savefig(os.path.join(here, "closure_driver.png"), dpi=150)

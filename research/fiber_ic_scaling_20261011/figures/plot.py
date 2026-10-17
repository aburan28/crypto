"""Plot envelope S_IC / S_rho against log2 #E from tables.json (written by
examples/fiber_scaling_tables.rs), with the fits that tool computed."""
import json, os
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt
here = os.path.dirname(os.path.abspath(__file__))
t = json.load(open(os.path.join(here, "../tables.json")))
fig, ax = plt.subplots(figsize=(7.2, 4.4))
for name, label, color, marker in (("scaling_h4", "h = 4", "#3b6ea8", "o"),
                                   ("scaling_h72_combined", "h = 72 (with amendment 1)", "#c0632b", "s"),
                                   ("scaling_h72", "h = 72, pre-registered sweep only", "#c0632b", "x")):
    d = t[name]
    xs = [c["log2_order"] for c in d["curves"]]; ys = [c["ratio"] for c in d["curves"]]
    f = d["fit"]
    lab = f"{label}: slope {f['slope']:.3f} [{f['slope_low']:.3f}, {f['slope_high']:.3f}]"
    ax.scatter(xs, ys, color=color, marker=marker, s=40, zorder=3, label=lab,
               facecolors=color if marker != "x" else None)
    if name != "scaling_h72":
        lo, hi = min(xs), max(xs)
        ax.plot([lo, hi], [2 ** (f["intercept"] + f["slope"] * lo), 2 ** (f["intercept"] + f["slope"] * hi)],
                color=color, lw=1.5)
ax.axhline(1, color="black", lw=1, ls=":")
ax.text(17.3, 1.15, "rho (ratio 1)", fontsize=8)
ax.set_yscale("log", base=2)
ax.set_xlabel("log2 #E")
ax.set_ylabel("best fiber-closed IC S ÷ rho S (per curve)")
ax.set_title("Fiber-closed MITM IC against rho: the gap grows as #E^(1/6)  (prediction 0.167)", fontsize=9)
ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=7, loc="lower right")
fig.tight_layout(); fig.savefig(os.path.join(here, "envelope_vs_size.svg")); print("ok")

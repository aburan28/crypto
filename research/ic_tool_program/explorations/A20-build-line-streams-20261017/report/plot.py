"""A20's smoke: each process's build time, A19 against A20, per size.
Lays out smoke.tsv; computes nothing new."""
import csv, sys
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

rows = list(csv.DictReader((l for l in open(sys.argv[1]) if not l.startswith("#")), delimiter="\t"))
out = sys.argv[2]
curves = [("icv1-f2m53-tm56619371-dac20a85", "2^44.3"), ("icv1-f2m61-t158598901-ab42b6c5", "2^47.2")]
fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.4))
for ax, (c, label) in zip(axes, curves):
    for i, b in enumerate(["ic-a19-v8-6c998fb0", "ic-a20-v8-21564f85"]):
        v = [float(r["build_ms"]) for r in rows if r["curve"] == c and r["binary"] == b]
        ax.plot([i] * len(v), v, "o", ms=6, color=["#1f5fa8", "#c2571a"][i])
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["A19", "A20 (line streams)"], fontsize=8)
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylabel("build, ms (median repetition)")
    ax.set_title(f"{label}  {c}", fontsize=7.5)
fig.suptitle("A20 smoke: two isolated processes per arm and size; no visible build gain", fontsize=9)
fig.text(0.01, 0.01, "explorations/A20-build-line-streams-20261017/smoke.tsv", fontsize=6.5)
fig.tight_layout(rect=(0, 0.04, 1, 0.93))
fig.savefig(out + ".svg")
fig.savefig(out + ".png", dpi=150)

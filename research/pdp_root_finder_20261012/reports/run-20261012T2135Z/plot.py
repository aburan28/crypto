# Run report graph for run-20261012T2135Z: per-round wall time of every timed
# row, read from the raw records of the timing stage.  Lays out recorded
# values only.
import json, os, collections
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
raw = os.path.join(here, "..", "..", "raw")
C = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "yellow": "#eda100",
     "ink": "#0b0b0b", "ink2": "#52514e", "grid": "#e6e5e1", "surface": "#fcfcfb", "wash": "#9aa3ad"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": C["grid"], "axes.labelcolor": C["ink2"],
                     "xtick.color": C["ink2"], "ytick.color": C["ink2"], "axes.titlecolor": C["ink"],
                     "figure.facecolor": C["surface"], "axes.facecolor": C["surface"]})

def rows(name):
    out = []
    p = os.path.join(raw, name)
    if not os.path.exists(p):
        return out
    for line in open(p):
        line = line.strip()
        if line:
            try:
                out.append(json.loads(line))
            except ValueError:
                pass
    return out

# Which labels were contended, from the isolation tool's records.
contended = set()
for r in rows("isolated.jsonl"):
    if r.get("run", {}).get("contended"):
        contended.add(r.get("label", ""))

STYLE = {"baseline": (C["wash"], "o"), "group-scalar": (C["yellow"], "D"), "rf-pairs": (C["orange"], "s"), "rf-mitm": (C["blue"], "^"),
         "pairs3": (C["orange"], "s"), "mitm3": (C["blue"], "^"), "subtract2": (C["aqua"], "o"), "mitm4": (C["yellow"], "D")}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
# Binary: every round's ns per target, by arm and ℓ (n = 3ℓ rungs only).
bt = [r for r in rows("bin_time.jsonl") if r.get("label", "").startswith("round") and r["n"] == 3 * r["l"]]
by = collections.defaultdict(list)
for r in bt:
    by[(r["arm"], r["l"])].append((r["ns_per_target"], r["label"] in contended))
for arm in ["baseline", "group-scalar", "rf-pairs", "rf-mitm"]:
    ls = sorted({l for (a, l) in by if a == arm})
    col, mk = STYLE[arm]
    med = []
    for l in ls:
        vals = [v for v, c in by[(arm, l)] if not c]
        cont = [v for v, c in by[(arm, l)] if c]
        ax1.scatter([l] * len(vals), [v / 1e3 for v in vals], color=col, s=14, zorder=3)
        if cont:
            ax1.scatter([l] * len(cont), [v / 1e3 for v in cont], facecolors="none", edgecolors=C["ink"], s=30, zorder=4)
        if vals:
            med.append((l, sorted(vals)[len(vals) // 2]))
    if med:
        ax1.plot([m[0] for m in med], [m[1] / 1e3 for m in med], color=col, marker=mk, markersize=6, linewidth=1.5, label=f"{arm} (median of rounds)")
ax1.set_yscale("log")
ax1.set_xlabel("ℓ (n = 3ℓ); one dot per round; hollow rings: runs the isolation tool flagged as contended")
ax1.set_ylabel("µs per target")
ax1.set_title("Binary arms, every timed round", fontsize=9)
ax1.legend(fontsize=7, frameon=False)
# Prime: every round's ns per full scan by arm and t.
pt = [r for r in rows("prime_time.jsonl") if r.get("label", "").startswith("round")]
by = collections.defaultdict(list)
for r in pt:
    by[(r["arm"], r["t"])].append((r["ns_per_target"], r["label"] in contended))
for arm in ["pairs3", "mitm3", "subtract2", "mitm4"]:
    ts = sorted({t for (a, t) in by if a == arm})
    col, mk = STYLE[arm]
    med = []
    for t in ts:
        vals = [v for v, c in by[(arm, t)] if not c]
        cont = [v for v, c in by[(arm, t)] if c]
        ax2.scatter([t] * len(vals), [v / 1e3 for v in vals], color=col, s=14, zorder=3)
        if cont:
            ax2.scatter([t] * len(cont), [v / 1e3 for v in cont], facecolors="none", edgecolors=C["ink"], s=30, zorder=4)
        if vals:
            med.append((t, sorted(vals)[len(vals) // 2]))
    if med:
        ax2.plot([m[0] for m in med], [m[1] / 1e3 for m in med], color=col, marker=mk, markersize=6, linewidth=1.5, label=f"{arm} (median of rounds)")
ax2.set_yscale("log")
ax2.set_xlabel("tower height t (Kummer, p₁ = 2013265921); one dot per round")
ax2.set_ylabel("µs per full scan (8 uniform targets)")
ax2.set_title("Prime arms, every timed round", fontsize=9)
ax2.legend(fontsize=7, frameon=False)
for ax in (ax1, ax2):
    ax.grid(True, color=C["grid"], linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
fig.suptitle("run-20261012T2135Z: wall time of every timed row on the reserved core (tools/isolated_bench.py, CPU 3)", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(here, "rounds.svg"))
fig.savefig(os.path.join(here, "rounds.png"), dpi=150)
print("ok", len(bt), "binary rows,", len(pt), "prime rows,", len(contended), "contended labels")

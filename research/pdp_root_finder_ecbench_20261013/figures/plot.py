# Lays out the sessions' records (records.jsonl, written and sealed by
# ecbench) as the study's figures: S per arm and size with the floor, the
# phase split of the IC arms, and solve wall time.  Recorded values only;
# the means, ratios and intervals the text quotes come from `ecbench table`
# and `ecbench compare`, not from here.
import json, os, sys, collections
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
root = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "sessions")
C = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "yellow": "#eda100",
     "ink": "#0b0b0b", "ink2": "#52514e", "grid": "#e6e5e1", "surface": "#fcfcfb", "wash": "#9aa3ad", "plum": "#8e5bd1"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": C["grid"], "axes.labelcolor": C["ink2"],
                     "xtick.color": C["ink2"], "ytick.color": C["ink2"], "axes.titlecolor": C["ink"],
                     "figure.facecolor": C["surface"], "axes.facecolor": C["surface"]})
STYLE = {"rho-strong": (C["ink"], "o"), "rho": (C["wash"], "o"), "ic-f4": (C["plum"], "D"), "ic-mitm": (C["yellow"], "s"),
         "ic-rootfind": (C["blue"], "^"), "ic-rootfind-pairs": (C["orange"], "v"), "rho-strong-aa": (C["aqua"], "o")}
ORDER = ["rho-strong", "rho-strong-aa", "rho", "ic-f4", "ic-mitm", "ic-rootfind", "ic-rootfind-pairs"]

sessions = {}
for name in sorted(os.listdir(root)):
    p = os.path.join(root, name, "records.jsonl")
    if not os.path.exists(p):
        continue
    rows = [json.loads(l) for l in open(p) if l.strip()]
    rows = [r for r in rows if not r.get("warmup")]
    if not rows:
        continue
    n = rows[0]["workload"]["curve_spec"]["n"]
    sessions[n] = rows
ns = sorted(sessions)
if not ns:
    print("no sessions under", root)
    sys.exit(0)

def style(ax):
    ax.grid(True, color=C["grid"], linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(here, name + ".svg"))
    fig.savefig(os.path.join(here, name + ".png"), dpi=150)

# ── Figure 1: S per run, by arm and size; the floor; verified status ──
fig, ax = plt.subplots(figsize=(11, 4.6))
arms_seen = []
for n in ns:
    rows = sessions[n]
    floor = rows[0]["boundaries"]["floor_s"]
    ax.hlines(floor, n - 0.4, n + 0.4, color=C["ink2"], linestyle=":", linewidth=1.2, label="floor √(π/4n)" if n == ns[0] else None)
    for k, arm in enumerate(ORDER):
        vals = [r["cost"]["s"] for r in rows if r["arm"] == arm and r["outcome"]["status"] == "verified"]
        bad = [r["cost"]["s"] for r in rows if r["arm"] == arm and r["outcome"]["status"] != "verified" and r["cost"].get("s")]
        if not vals and not bad:
            continue
        col, mk = STYLE[arm]
        x = n + (k - 3) * 0.1
        ax.scatter([x] * len(vals), vals, color=col, marker=mk, s=16, zorder=3, label=arm if arm not in arms_seen else None)
        if bad:
            ax.scatter([x] * len(bad), bad, facecolors="none", edgecolors=col, marker=mk, s=40, zorder=4)
        arms_seen.append(arm)
ax.set_yscale("log")
ax.set_xticks(ns)
ax.set_xlabel("field degree n (E_0 over F_2^n; one marker per verified measured run, hollow: not verified)")
ax.set_ylabel("S = charged group-addition equivalents / √r")
ax.set_title("Every measured run: S by arm and size (IC with the F4 oracle is a lower bound: its solver's word XORs are uncharged)", fontsize=9)
ax.legend(fontsize=7, frameon=False, ncol=2)
style(ax)
save(fig, "s_by_size")

# ── Figure 2: phase split of the IC arms (mean GAE per phase over verified runs) ──
ic_arms = [a for a in ["ic-f4", "ic-mitm", "ic-rootfind", "ic-rootfind-pairs"] if any(r["arm"] == a for n in ns for r in sessions[n])]
fig, axes = plt.subplots(1, len(ns), figsize=(2.4 * len(ns) + 2, 4.2), sharey=False)
if len(ns) == 1:
    axes = [axes]
phase_cols = {"factor_base": C["wash"], "oracle_setup": C["yellow"], "relations": C["blue"], "linear_algebra": C["orange"], "verification": C["aqua"]}
for ax, n in zip(axes, ns):
    rows = sessions[n]
    xs = []
    for k, arm in enumerate(ic_arms):
        runs = [r for r in rows if r["arm"] == arm and r["outcome"]["status"] == "verified"]
        if not runs:
            continue
        bottom = 0.0
        sqrt_r = 2 ** (rows[0]["workload"]["curve"].get("log2_r", 0)) if False else None
        for ph in ["factor_base", "oracle_setup", "relations", "linear_algebra", "verification"]:
            vals = []
            for r in runs:
                g = sum(p["gae"] for p in r["phases"] if p["name"] == ph)
                vals.append(g / (r["cost"]["total_gae"] / r["cost"]["s"]) if r["cost"]["total_gae"] else 0.0)
            mean = sum(vals) / len(vals)
            ax.bar(k, mean, bottom=bottom, color=phase_cols[ph], width=0.7, label=ph if (n == ns[0] and k == 0) else None)
            bottom += mean
        xs.append(k)
    ax.set_xticks(range(len(ic_arms)))
    ax.set_xticklabels(ic_arms, rotation=30, ha="right", fontsize=7)
    ax.set_title(f"n = {n}", fontsize=9)
    ax.set_yscale("log")
    style(ax)
axes[0].set_ylabel("mean S by phase (stacked; log scale)")
axes[0].legend(fontsize=7, frameon=False)
fig.suptitle("IC arms: where the charged work sits, by phase (means over verified runs; the F4 arm's solver term is wall-priced in its relations phase)", fontsize=9)
save(fig, "phase_split")

# ── Figure 3: solve wall time per run ──
fig, ax = plt.subplots(figsize=(11, 4.4))
arms_seen = []
for n in ns:
    rows = sessions[n]
    for k, arm in enumerate(ORDER):
        vals = [r["time"]["solve_wall_ns"] / 1e6 for r in rows if r["arm"] == arm and r["outcome"]["status"] == "verified"]
        if not vals:
            continue
        col, mk = STYLE[arm]
        x = n + (k - 3) * 0.1
        ax.scatter([x] * len(vals), vals, color=col, marker=mk, s=16, zorder=3, label=arm if arm not in arms_seen else None)
        arms_seen.append(arm)
ax.set_yscale("log")
ax.set_xticks(ns)
ax.set_xlabel("field degree n")
ax.set_ylabel("solve wall time per run (ms; the record's solve_wall_ns)")
ax.set_title("Wall time of every verified measured run, by arm and size (descriptive: read the admitted comparisons for a wall result)", fontsize=9)
ax.legend(fontsize=7, frameon=False, ncol=2)
style(ax)
save(fig, "wall_by_size")
print("ok", {n: len(sessions[n]) for n in ns})

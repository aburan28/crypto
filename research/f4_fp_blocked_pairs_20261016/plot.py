"""Lays out the recorded A/B/msolve wall times (raw/ab_full.tsv); computes
only medians and ratios of recorded values for display."""
import sys, statistics, collections
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

raw, out = sys.argv[1], sys.argv[2]
t = collections.defaultdict(list)
order = []
for line in open(raw):
    r, case, arm, ms, fp = line.rstrip("\n").split("\t")
    if case not in order:
        order.append(case)
    t[(case, arm)].append(float(ms))
med = {k: statistics.median(v) for k, v in t.items()}
rows = []
for c in order:
    b, a, m = med[(c, "ours_B1")], med[(c, "ours_B0")], med[(c, "msolve")]
    rows.append((c, b, a, m, b / a, b / m, a / m,
                 min(t[(c, "ours_B1")]), max(t[(c, "ours_B1")]),
                 min(t[(c, "ours_B0")]), max(t[(c, "ours_B0")])))
with open(f"{out}/table.tsv", "w") as f:
    f.write("case\tblocked_ms\tcontrol_ms\tmsolve_ms\tblocked/control\tblocked/msolve\tcontrol/msolve\tblocked_min\tblocked_max\tcontrol_min\tcontrol_max\n")
    for r in rows:
        f.write("\t".join([r[0]] + [f"{x:.3f}" for x in r[1:]]) + "\n")

fig, ax = plt.subplots(figsize=(9, 4.6))
x = range(len(rows))
ax.axhline(1.0, color="#888", lw=0.8)
ax.plot(x, [r[4] for r in rows], "o", color="#1f6fb4", label="blocked / control (same binary)")
ax.plot(x, [r[5] for r in rows], "s", color="#d1495b", label="blocked / msolve 0.6.5")
ax.plot(x, [r[6] for r in rows], "^", color="#999", mfc="none", label="control / msolve 0.6.5")
ax.set_yscale("log")
ax.set_xticks(list(x))
ax.set_xticklabels([r[0] for r in rows], rotation=60, ha="right", fontsize=7)
ax.set_ylabel("ratio of median wall time (log; < 1 = blocked faster)")
ax.set_title("f4_fp, 1 thread: median of 3 interleaved rounds per arm (run-20261016T1720Z)", fontsize=9)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(f"{out}/fig_ratio.svg")
fig.savefig(f"{out}/fig_ratio.png", dpi=150)

big = [r for r in rows if r[3] >= 20]
fig, ax = plt.subplots(figsize=(8, 4.4))
w = 0.27
xs = range(len(big))
for k, (idx, lab, col) in enumerate([(1, "f4_fp blocked", "#1f6fb4"), (2, "f4_fp control", "#9ecae1"), (3, "msolve 0.6.5", "#d1495b")]):
    ax.bar([i + (k - 1) * w for i in xs], [r[idx] for r in big], w, label=lab, color=col)
ax.set_yscale("log")
ax.set_xticks(list(xs))
ax.set_xticklabels([r[0] for r in big], rotation=45, ha="right", fontsize=7)
ax.set_ylabel("median wall time, ms (log)")
ax.set_title("Cases where msolve takes ≥ 20 ms; 1 thread, n = 3 rounds (run-20261016T1720Z)", fontsize=9)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(f"{out}/fig_times.svg")
fig.savefig(f"{out}/fig_times.png", dpi=150)

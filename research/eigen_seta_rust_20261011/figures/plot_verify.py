"""Graphs for run-20261011T1955Z-verify, drawn from its records only
(plotting recorded values; nothing is computed beyond reading them).
Usage: python plot_verify.py <records_dir> <out_dir>
"""
import os, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as tk

rec, out = sys.argv[1:3]
RID = os.path.basename(os.path.dirname(os.path.abspath(rec)))

def summary(path):
    rows = [l.split() for l in open(path).read().splitlines()[1:] if l.strip()]
    return [(r[0].replace("cg.", ""), int(r[2]), int(r[3])) for r in rows]

panels = [("Toy16", summary(os.path.join(rec, "toyprobe/ctprobe_summary.txt"))),
          ("L1", summary(os.path.join(rec, "l1probe/ctprobe_summary.txt")))]
fig, axes = plt.subplots(2, 2, figsize=(8, 5), gridspec_kw={"width_ratios": [3, 1.5]})
fmt = tk.FuncFormatter(lambda v, _: f"{v:,.0f}")
for col, (name, rows) in enumerate(panels):
    labels = [r[0] for r in rows]; alloc = [r[1] for r in rows]; rest = [r[2] for r in rows]
    a1, a2 = axes[0][col], axes[1][col]
    a1.bar(range(len(rows)), rest, color="tab:blue")
    lo, hi = min(rest), max(rest)
    pad = max(1, lo // 20000)
    a1.set_ylim(lo - pad, hi + pad); a1.yaxis.set_major_formatter(fmt); a1.yaxis.set_major_locator(tk.MaxNLocator(3))
    a1.set_title(f"{name}: non-allocator {'identical' if lo == hi else 'DIFFERS'}\n{lo:,}", fontsize=8)
    a1.set_xticks([]); a1.tick_params(labelsize=7)
    a2.bar(range(len(rows)), alloc, color="tab:orange")
    span = max(alloc) - min(alloc)
    a2.set_ylim(min(alloc) - span - 1, max(alloc) + span + 1); a2.yaxis.set_major_formatter(fmt); a2.yaxis.set_major_locator(tk.MaxNLocator(3))
    a2.set_xticks(range(len(rows))); a2.set_xticklabels(labels, rotation=70, fontsize=6); a2.tick_params(axis="y", labelsize=7)
axes[0][0].set_ylabel("non-allocator\ninstructions", fontsize=8)
axes[1][0].set_ylabel("allocator + mem*\ninstructions", fontsize=8)
fig.text(0.01, 0.005, f"measured (callgrind, decaps_core only); records: {RID}/records/{{toyprobe,l1probe}}/ctprobe_summary.txt", fontsize=6)
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(os.path.join(out, "ct_probe.svg")); plt.close(fig)

txt = open(os.path.join(rec, "f_filter_rate.txt")).read()
rows = re.findall(r"^(\w+) drawn (\d+) accepted (\d+) rate ([0-9.]+) \(log2 (-?[0-9.]+)\)", txt, re.M)
fig, ax = plt.subplots(figsize=(6, 3.4))
names = [r[0] for r in rows]; lg = [float(r[4]) for r in rows]; acc = [int(r[2]) for r in rows]
err = [1.4427 / (a ** 0.5) for a in acc]  # binomial standard error of log2(rate), approx 1/(ln2 sqrt(k))
ax.bar(names, lg, yerr=err, capsize=3, color="tab:purple")
ax.set_xticks(range(len(names))); ax.set_xticklabels([f"{n}\n{v:.2f} (k={a})" for n, v, a in zip(names, lg, acc)], fontsize=7)
ax.set_ylabel("log2(acceptance rate of public f filters)", fontsize=8); ax.tick_params(labelsize=8)
ax.set_title("Fraction of odd f candidates accepted by key generation's public filters", fontsize=9)
fig.text(0.01, 0.005, f"measured; error bars approx. 1/(ln2*sqrt(k)), k accepted; records: {RID}/records/f_filter_rate.txt", fontsize=6)
fig.tight_layout(rect=(0, 0.04, 1, 1)); fig.savefig(os.path.join(out, "f_filter.svg")); plt.close(fig)
print("wrote", out)

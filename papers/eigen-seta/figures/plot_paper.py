"""Figures for the Eigen-Seta paper, drawn from committed run records only
(plotting recorded values; no result is computed here).

Usage: python figures/plot_paper.py <research_dir>
  <research_dir> = crypto:research/eigen_seta_rust_20261011
Writes figures/kem_times.pdf and figures/ct_probe.pdf.
"""
import os, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

R = sys.argv[1]
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)))
SETS = ["Toy16", "Toy32", "L1", "L3", "L5"]

txt = open(os.path.join(R, "reports/run-20261011T1813Z-bench/records/bench_output.txt")).read()
keygen = {m[0]: float(m[1]) for m in re.findall(r"^(\w+) keygen_s ([0-9.]+)", txt, re.M)}
enc, dec = {}, {}
for name, e, d in re.findall(r"^(\w+) trial \d+ encaps_s ([0-9.]+) decaps_s ([0-9.]+) ok true", txt, re.M):
    enc.setdefault(name, []).append(float(e)); dec.setdefault(name, []).append(float(d))
fig, ax = plt.subplots(figsize=(6.4, 3.6))
for off, (lab, data, mk, col) in zip([-0.2, 0, 0.2], [("keygen (1 run)", {s: [keygen[s]] for s in SETS}, "s", "tab:blue"), ("encaps (per trial)", enc, "o", "tab:orange"), ("decaps (per trial)", dec, "^", "tab:green")]):
    for i, s in enumerate(SETS):
        ys = data.get(s, [])
        ax.scatter([i + off] * len(ys), ys, marker=mk, color=col, label=lab if i == 0 else None, zorder=3)
ax.set_xticks(range(len(SETS))); ax.set_xticklabels([f"{s}\n(n={len(dec.get(s, []))})" for s in SETS])
ax.set_yscale("log"); ax.set_ylabel("wall time [s], log scale"); ax.grid(True, which="both", axis="y", alpha=0.3)
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(os.path.join(OUT, "kem_times.pdf")); plt.close(fig)

def summary(path):
    rows = [l.split() for l in open(path).read().splitlines()[1:] if l.strip()]
    return [(r[0].replace("cg.", ""), int(r[2]), int(r[3])) for r in rows]

V = os.path.join(R, "reports/run-20261011T1955Z-verify/records")
panels = [("Toy16", summary(os.path.join(V, "toyprobe/ctprobe_summary.txt"))), ("L1", summary(os.path.join(V, "l1probe/ctprobe_summary.txt")))]
fig, axes = plt.subplots(2, 2, figsize=(7.6, 4.6), gridspec_kw={"width_ratios": [3, 1.5]})
for col, (name, rows) in enumerate(panels):
    labels = [r[0] for r in rows]; alloc = [r[1] for r in rows]; rest = [r[2] for r in rows]
    a1, a2 = axes[0][col], axes[1][col]
    a1.bar(range(len(rows)), rest, color="tab:blue")
    lo, hi = min(rest), max(rest)
    pad = max(1, lo // 20000)
    a1.set_ylim(lo - pad, hi + pad)
    a1.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    a1.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
    a1.set_title(f"{name}: {'identical' if lo == hi else 'DIFFER'}\n{lo:,}", fontsize=8)
    a1.set_xticks([]); a1.tick_params(labelsize=7); a1.yaxis.get_offset_text().set_fontsize(7)
    a2.bar(range(len(rows)), alloc, color="tab:orange")
    span = max(alloc) - min(alloc)
    a2.set_ylim(min(alloc) - span - 1, max(alloc) + span + 1)
    a2.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    a2.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
    a2.set_xticks(range(len(rows))); a2.set_xticklabels(labels, rotation=70, fontsize=6)
    a2.tick_params(axis="y", labelsize=7); a2.yaxis.get_offset_text().set_fontsize(7)
axes[0][0].set_ylabel("non-allocator", fontsize=8)
axes[1][0].set_ylabel("allocator + mem*", fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "ct_probe.pdf")); plt.close(fig)
print("wrote", OUT)

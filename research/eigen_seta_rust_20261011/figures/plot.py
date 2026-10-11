"""Plot recorded values only (no computation of results).

Inputs are the raw records committed under ../reports/run-*/records/.
Usage: python plot.py <bench_records_dir> <ctprobe_records_dir> <out_dir>
"""
import re, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

bench, ctp, out = sys.argv[1:4]
BID = os.path.basename(os.path.dirname(os.path.abspath(bench)))
CID = os.path.basename(os.path.dirname(os.path.abspath(ctp)))
os.makedirs(out, exist_ok=True)
SETS = ["Toy16", "Toy32", "L1", "L3", "L5"]

# ---- KEM timings (recorded per trial) ----
txt = open(os.path.join(bench, "bench_output.txt")).read()
keygen = {m[0]: float(m[1]) for m in re.findall(r"^(\w+) keygen_s ([0-9.]+)", txt, re.M)}
enc, dec = {}, {}
for name, e, d in re.findall(r"^(\w+) trial \d+ encaps_s ([0-9.]+) decaps_s ([0-9.]+) ok true", txt, re.M):
    enc.setdefault(name, []).append(float(e)); dec.setdefault(name, []).append(float(d))
sets = [s for s in SETS if s in keygen]
fig, ax = plt.subplots(figsize=(7, 4))
x = range(len(sets))
for off, (lab, data, mk, col) in zip([-0.2, 0, 0.2], [("keygen (1 run)", {s: [keygen[s]] for s in sets}, "s", "tab:blue"), ("encaps (per trial)", enc, "o", "tab:orange"), ("decaps (per trial)", dec, "^", "tab:green")]):
    for i, s in enumerate(sets):
        ys = data.get(s, [])
        ax.scatter([i + off] * len(ys), ys, marker=mk, color=col, label=lab if i == 0 else None, zorder=3)
ax.set_xticks(list(x)); ax.set_xticklabels([f"{s}\n(n={len(dec.get(s, []))})" for s in sets])
ax.set_yscale("log"); ax.set_ylabel("wall time [s], log scale"); ax.grid(True, which="both", axis="y", alpha=0.3)
ax.set_title("Eigen-Séta v2 (Rust): KEM wall time per operation, one thread")
ax.legend(fontsize=8)
fig.text(0.01, 0.01, f"measured; records: {BID}/records/bench_output.txt", fontsize=7)
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(os.path.join(out, "kem_times.svg")); plt.close(fig)

# ---- field operations (recorded, same binary) ----
rows = re.findall(r"^(\w+)\s+limbs\s+(\d+)\s+mul\s+(\d+) ns\s+fermat\s+(\d+) ns\s+safegcd\s+(\d+) ns", txt, re.M)
limbs = [int(r[1]) for r in rows]
fig, ax = plt.subplots(figsize=(7, 4))
for idx, lab, mk in [(2, "F_p mul (Montgomery)", "o"), (3, "inversion, Fermat", "s"), (4, "inversion, safegcd", "^")]:
    ax.plot(limbs, [int(r[idx]) for r in rows], marker=mk, label=lab)
for r in rows:
    ax.annotate(r[0], (int(r[1]), int(r[3])), textcoords="offset points", xytext=(0, 6), fontsize=7, ha="center")
ax.set_xlabel("limbs (64-bit)"); ax.set_ylabel("time per operation [ns], log scale"); ax.set_yscale("log")
ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=8)
ax.set_title("F_p arithmetic: best of 15 runs (mul), mean of 200 (inversions)")
fig.text(0.01, 0.01, f"measured; records: {BID}/records/bench_output.txt", fontsize=7)
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(os.path.join(out, "field_ops.svg")); plt.close(fig)

# ---- constant-time probe (recorded instruction counts) ----
lines = [l.split() for l in open(os.path.join(ctp, "ctprobe_summary.txt")).read().splitlines()[1:] if l.strip()]
names = [l[0].replace("cg.", "") for l in lines]
alloc = [int(l[2]) for l in lines]; rest = [int(l[3]) for l in lines]
fig, (a1, a2) = plt.subplots(2, 1, figsize=(7, 5), sharex=True)
a1.bar(range(len(names)), rest, color="tab:blue")
lo, hi = min(rest), max(rest)
a1.set_ylim(lo - 2e6, hi + 2e6); a1.set_ylabel("non-allocator\ninstructions")
a1.set_title(f"decaps_core instruction counts (callgrind, Toy16)\nnon-allocator count {'identical in all runs' if lo == hi else 'DIFFERS between runs'}: {lo:,}" + ("" if lo == hi else f" .. {hi:,}"), fontsize=10)
a2.bar(range(len(names)), alloc, color="tab:orange"); a2.set_ylabel("glibc allocator\n+ mem* instructions")
a2.set_ylim(min(alloc) - 3e5, max(alloc) + 3e5)
a2.set_xticks(range(len(names))); a2.set_xticklabels(names, rotation=60, fontsize=7)
a2.set_xlabel("key seed . ciphertext seed . mode (valid = accepted, foreign = FO-rejected)")
fig.text(0.01, 0.01, f"measured; records: {CID}/records/ctprobe_summary.txt", fontsize=7)
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(os.path.join(out, "ct_probe.svg")); plt.close(fig)
print("wrote", out)

"""Lay out the recorded hwcost run as graphs (SVG + 150 dpi PNG).

Plotting only: every value drawn is read from the run directory written by
`hwcost measure` / `hwcost fit` / `hwcost rho`; nothing is computed here
beyond unit conversions for display (AGENTS.md "Implementation language").

    python plot.py RUN_DIR OUT_DIR
"""

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

RUN = Path(sys.argv[1]).resolve()
OUT = Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
_started = (RUN / "started_utc.txt").read_text().strip() if (RUN / "started_utc.txt").exists() else "?"
_commit = (RUN / "commit.txt").read_text().strip()[:9] if (RUN / "commit.txt").exists() else "?"
RUN_ID = f"{RUN.parent.parent.name}/{RUN.parent.name}/{RUN.name} @ {_started}, {_commit}"
FIELDS = ["fp64", "w1", "w2", "w3", "w4"]
LABEL = {
    "fp64": "fp64 (Fp64, 2^63-25)",
    "w1": "w1 (FpWide<1>, 2^64-59)",
    "w2": "w2 (FpWide<2>, 2^128-159)",
    "w3": "w3 (FpWide<3>, P-192)",
    "w4": "w4 (FpWide<4>, P-256)",
}
COLOR = {"fp64": "#1b6ca8", "w1": "#4d9de0", "w2": "#e1a100", "w3": "#d1495b", "w4": "#2e8540"}
plt.rcParams.update({"font.size": 9, "svg.fonttype": "none"})


def save(fig, name):
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def read_csv(p):
    with open(p) as f:
        return list(csv.DictReader(f))


rows = read_csv(RUN / "validation.csv")
val = json.loads((RUN / "validation.json").read_text())
med = {r["arm"]: float(r["median"]) for r in read_csv(RUN / "medians.csv")}
meas = json.loads((RUN / "measurements.json").read_text())
fit_max = val["fit_max_batch"]

# 1. Measured kernel time against batch size, with models A, C, S.
fig, axes = plt.subplots(1, 5, figsize=(15, 3.4), sharex=True)
for ax, f in zip(axes, FIELDS):
    rs = sorted((r for r in rows if r["field"] == f), key=lambda r: int(r["batch"]))
    if not rs:
        ax.set_visible(False)
        continue
    b = [int(r["batch"]) for r in rs]
    m = [float(r["measured_ns"]) for r in rs]
    lo = [m_ - float(r["measured_q1"]) for m_, r in zip(m, rs)]
    hi = [float(r["measured_q3"]) - m_ for m_, r in zip(m, rs)]
    ax.axvspan(min(b) / 1.5, fit_max * 1.4, color="#eeeeee", zorder=0)
    ax.errorbar(b, m, yerr=[lo, hi], fmt="o", color="black", ms=4, capsize=2, label="measured (median, IQR)")
    ax.plot(b, [float(r["a_priori_ns"]) for r in rs], "--", color=COLOR[f], label="A: a priori")
    ax.plot(b, [float(r["calibrated_ns"]) for r in rs], "-", color=COLOR[f], label="C: calibrated")
    s = [(bb, float(r["size_extrapolated_ns"])) for bb, r in zip(b, rs) if r["size_extrapolated_ns"]]
    if s:
        ax.plot([x for x, _ in s], [y for _, y in s], ":", color="purple", lw=2, label="S: size-extrapolated")
    ax.set_xscale("log", base=2)
    ax.set_title(LABEL[f], fontsize=8.5)
    ax.set_xlabel("walks per batch B")
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.3)
axes[0].set_ylabel("ns per rho iteration (1 thread)")
handles, labels = [], []
for ax in axes:
    for h, l in zip(*ax.get_legend_handles_labels()):
        if l not in labels:
            handles.append(h)
            labels.append(l)
fig.legend(handles, labels, loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.08), frameon=False)
fig.suptitle(
    f"Rho kernel: measured vs predicted (run {RUN_ID}; n = {meas['params']['reps']} reps per point; grey = fit set B ≤ {fit_max})",
    y=1.14,
    fontsize=9,
)
save(fig, "fig_validation")

# 2. Prediction error per field and batch size.
fig, ax = plt.subplots(figsize=(8, 3.6))
models = [("err_a_priori", "A: a priori", "s", 0.0), ("err_calibrated", "C: calibrated", "o", 0.0), ("err_size_extrapolated", "S: size-extrapolated", "^", 0.0)]
xt, xl = [], []
x = 0
for f in FIELDS:
    rs = sorted((r for r in rows if r["field"] == f), key=lambda r: int(r["batch"]))
    if rs:
        ax.text(x + (len(rs) - 1) / 2, -0.12, f, transform=ax.get_xaxis_transform(), ha="center", fontsize=8, fontweight="bold")
    for r in rs:
        for key, name, mk, _ in models:
            if r[key]:
                ax.plot(
                    x,
                    100 * float(r[key]),
                    mk,
                    color=COLOR[f],
                    mfc="white" if key == "err_a_priori" else (COLOR[f] if key == "err_calibrated" else "purple"),
                    ms=5,
                )
        xt.append(x)
        xl.append(str(int(r["batch"]).bit_length() - 1))
        x += 1
    x += 1
ax.axhspan(-30, 30, color="#e8f4e8", zorder=0, label="±30% (protocol threshold)")
ax.axhline(0, color="black", lw=0.6)
ax.set_xticks(xt)
ax.set_xticklabels(xl, fontsize=5.5)
ax.set_ylabel("predicted / measured − 1 (%)")
ax.set_xlabel("log2 B (walks per batch), grouped by field", labelpad=22)
ax.set_title(f"Prediction error; open squares A (a priori), filled circles C (calibrated), purple triangles S (size-extrapolated)\nrun {RUN_ID}", fontsize=8.5)
ax.grid(axis="y", alpha=0.3)
save(fig, "fig_errors")

# 3. Memory roofline: measured read bandwidth against working-set size, and
#    the bandwidth each kernel arm actually drew (bytes moved / measured time).
fig, ax = plt.subplots(figsize=(7.5, 4))
sizes = meas["params"]["mem_sizes"]
threads = meas["params"]["threads"]
for t, style in [(1, "-o"), (threads[-1], "--o")]:
    ys = [med.get(f"mem.{t}t.{s}") for s in sizes]
    ax.plot([s for s, y in zip(sizes, ys) if y], [y for y in ys if y], style, color="black", ms=3, label=f"sequential read, {t} thread(s) (aggregate)")
for c in meas["caches"]:
    if c["kind"] != "Instruction":
        ax.axvline(c["size_bytes"], color="grey", lw=0.6, ls=":")
        ax.text(c["size_bytes"], 0.97, f" L{c['level']}", fontsize=7, color="grey", va="top", transform=ax.get_xaxis_transform())
for f in FIELDS:
    rs = sorted((r for r in rows if r["field"] == f), key=lambda r: int(r["batch"]))
    elem = {"fp64": 8, "w1": 8, "w2": 16, "w3": 24, "w4": 32}[f]
    ax.plot(
        [float(r["working_set_bytes"]) for r in rs],
        [12 * elem / float(r["measured_ns"]) for r in rs],
        "x-",
        color=COLOR[f],
        ms=4,
        lw=0.8,
        label=f"kernel {f}: 12×{elem} B / measured ns",
    )
ax.set_xscale("log", base=2)
ax.set_yscale("log")
ax.set_xlabel("working set per thread (bytes)")
ax.set_ylabel("GB/s")
ax.set_title(f"Roofline check: bandwidth available vs drawn (run {RUN_ID})")
ax.grid(alpha=0.3, which="both")
ax.legend(fontsize=7, loc="lower left")
save(fig, "fig_roofline")

# 4. Op cost against limb count, with the size fits and their extrapolation.
fig, axes = plt.subplots(1, 4, figsize=(14, 3.2))
fig.subplots_adjust(wspace=0.32)
ghz = val["clock_ghz"]
for ax, op in zip(axes, ["mul", "sqr", "add", "inv"]):
    for mode, mk in [("throughput", "o"), ("latency", "s")]:
        xs = [1, 2, 3, 4]
        ys = [med.get(f"w{n}.{op}.{mode}") for n in xs]
        ax.plot(xs, [y * ghz for y in ys], mk, color="black" if mode == "latency" else "#1b6ca8", label=f"measured {mode}")
        sf = next(s for s in val["size_fits"] if s["op"] == f"{op}.{mode}")
        xf = [1 + i * 0.05 for i in range(61)]
        quad = sf["model"].endswith("n^2")
        ax.plot(xf, [(sf["alpha"] + sf["slope"] * (x * x if quad else x)) * ghz for x in xf], "-", lw=0.8, color="black" if mode == "latency" else "#1b6ca8")
        ax.plot([4], [sf["extrapolated_ns"] * ghz], "*", ms=10, color="purple")
    fp = [med.get(f"fp64.{op}.throughput") * ghz, med.get(f"fp64.{op}.latency") * ghz]
    ax.plot([1, 1], fp, "D", mfc="none", color="#d1495b", label="fp64 (dedicated)")
    ax.set_title(f"{op}: {next(s['model'] for s in val['size_fits'] if s['op'] == op + '.throughput')} fitted on n = 1..3", fontsize=8)
    ax.set_xlabel("limbs n")
    ax.set_xticks([1, 2, 3, 4])
    ax.grid(alpha=0.3)
axes[0].set_ylabel(f"cycles per op (at {ghz:.2f} GHz measured)")
axes[0].legend(fontsize=7)
fig.suptitle(f"Field-op costs vs limb count; purple stars = extrapolation to n = 4 (run {RUN_ID})", y=1.03, fontsize=9)
save(fig, "fig_opcost")

# 5. Application: predicted rho wall time and cost against group size.
rho = []
for name, neg in [("rho.csv", False), ("rho_negation.csv", True)]:
    p = RUN / name
    if p.exists():
        for r in read_csv(p):
            r["neg"] = neg
            rho.append(r)
if rho:
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 3.8))
    year = 365.25 * 86400
    for mach, color in [("1.0", "#1b6ca8"), ("1000.0", "#d1495b")]:
        for neg, ls in [(False, "-"), (True, "--")]:
            rs = [r for r in rho if r["machines"] in (mach, mach.rstrip("0").rstrip(".")) and r["neg"] == neg]
            rs.sort(key=lambda r: float(r["group_log2"]))
            if not rs:
                continue
            g = [float(r["group_log2"]) for r in rs]
            lab = f"{float(mach):g} machine(s){', negation map' if neg else ''}"
            a1.plot(g, [float(r["wall_s"]) / year for r in rs], ls + "o", ms=3, color=color, label=lab)
            a2.plot(g, [float(r["total_usd"]) for r in rs], ls + "o", ms=3, color=color, label=lab)
    for a, yl in [(a1, "predicted wall time (years)"), (a2, "predicted cost (USD, illustrative price)")]:
        a.set_yscale("log")
        a.set_xlabel("log2 of the prime group order n")
        a.set_ylabel(yl)
        a.grid(alpha=0.3, which="major")
    a1.axhline(1, color="grey", lw=0.6)
    a1.text(48, 2, "1 year", fontsize=7, color="grey")
    a1.axhline(13.8e9, color="grey", lw=0.6)
    a1.text(48, 13.8e9 * 2, "age of the universe", fontsize=7, color="grey")
    a1.legend(fontsize=7)
    fig.suptitle(
        f"Parallel rho on this profile (4 cores per machine), extrapolated from run {RUN_ID}; price and network are illustrative/assumed",
        fontsize=9,
    )
    save(fig, "fig_rho")

print("wrote", sorted(p.name for p in OUT.glob("fig_*")))

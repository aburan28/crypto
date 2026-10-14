"""Lay out recorded values from raw/*.jsonl as SVG + 150 dpi PNG.

Layout only (AGENTS.md "Implementation language"): every number plotted is
read from the run records; medians and min-max over seeds are the only
arithmetic, and summary.json records them for the report.
"""
import json, statistics, sys, pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = pathlib.Path(__file__).resolve().parent.parent
raw = here / "raw"
figs = here / "figures"
figs.mkdir(exist_ok=True)

def load(name):
    p = raw / name
    if not p.exists():
        return []
    return [json.loads(l) for l in p.read_text().splitlines() if l.strip()]

grid = load("grid.jsonl")
amend = load("grid_amend1.jsonl")
e2e = load("e2e.jsonl")

SERIES = {  # arm -> (label, color slot, marker)
    "igauss": ("IncrementalGauss (dense, ecbench default)", "#e34948", "v"),
    "sgauss": ("StructuredGauss (sparse RREF)", "#eda100", "^"),
    "core-bw": ("ours: block Wiedemann, unmerged core", "#2a78d6", "o"),
    "cado-bwc": ("CADO bwc, same unmerged core", "#eb6834", "s"),
    "full-bw": ("ours: filter+merge + block Wiedemann", "#1baf7a", "D"),
    "full-bl": ("ours: filter+merge + block Lanczos", "#6250d6", "P"),
    "full-bw-pos": ("ours: filter+merge + BW (+1 coeffs)", "#008300", "d"),
    "cado-full": ("CADO purge+merge-dl+bwc (+1 coeffs)", "#e87ba4", "X"),
}

def cells(rows, arm, key="solve_s", weight=3, pass_=None, threads=None):
    by = {}
    for r in rows:
        if r.get("arm") != arm or r.get("failed") or not r.get("verified"):
            continue
        if r.get("weight", 3) != weight:
            continue
        if pass_ and r.get("pass") != pass_:
            continue
        if threads is not None and r.get("threads_requested", 4) != threads:
            continue
        by.setdefault(r["cols"], []).append(r[key])
    return {d: (statistics.median(v), min(v), max(v), len(v)) for d, v in sorted(by.items())}

def style(ax, title, ylabel):
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("factor-base columns D (m = 3 summands per relation)")
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontsize=11)
    ax.grid(True, which="major", color="#e4e3df", linewidth=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

def save(fig, name):
    fig.savefig(figs / f"{name}.svg", bbox_inches="tight")
    fig.savefig(figs / f"{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

summary = {}

def plot_lines(ax, rows, arms, pass_=None, threads=None):
    for arm in arms:
        c = cells(rows, arm, pass_=pass_, threads=threads)
        if not c:
            continue
        label, color, marker = SERIES[arm]
        xs = list(c)
        med = [c[d][0] for d in xs]
        lo = [c[d][0] - c[d][1] for d in xs]
        hi = [c[d][2] - c[d][0] for d in xs]
        ax.errorbar(xs, med, yerr=[lo, hi], color=color, marker=marker, markersize=6,
                    linewidth=2, capsize=2, label=label)
        summary.setdefault(f"{arm}|{pass_ or 'all'}|t{threads or 4}", {})
        summary[f"{arm}|{pass_ or 'all'}|t{threads or 4}"] = {
            str(d): {"median_s": v[0], "min_s": v[1], "max_s": v[2], "n": v[3]} for d, v in c.items()}

# 1. Every solver, wall time against D (four threads).
fig, ax = plt.subplots(figsize=(9, 5.5))
plot_lines(ax, grid, ["igauss", "sgauss", "core-bw", "cado-bwc", "full-bw", "full-bl"], pass_="main")
plot_lines(ax, amend, ["full-bw-pos", "cado-full"])
style(ax, "Linear-algebra solve, wall time (median of 3 seeds, bars min-max; 4 threads)", "wall seconds")
ax.legend(fontsize=8, frameon=False, loc="upper left")
save(fig, "solve_time_vs_D")

# 2. Engine parity on the identical unmerged core, 4 and 1 threads.
fig, ax = plt.subplots(figsize=(9, 4.5))
for t, ls in ((4, "-"), (1, "--")):
    pass_ = "main" if t == 4 else "onethread"
    a = {}
    b = {}
    for r in grid:
        if r.get("pass") != pass_ or r.get("failed") or not r.get("verified"):
            continue
        key = (r["cols"], r["seed"])
        if r["arm"] == "core-bw":
            a[key] = r["solve_s"]
        elif r["arm"] == "cado-bwc":
            b[key] = r["solve_s"]
    ratios = {}
    for k in a:
        if k in b:
            ratios.setdefault(k[0], []).append(b[k] / a[k])
    if not ratios:
        continue
    xs = sorted(ratios)
    med = [statistics.median(ratios[d]) for d in xs]
    ax.errorbar(xs, med, yerr=[[m - min(ratios[d]) for m, d in zip(med, xs)],
                               [max(ratios[d]) - m for m, d in zip(med, xs)]],
                color="#2a78d6" if t == 4 else "#eb6834", marker="o" if t == 4 else "s",
                linestyle=ls, linewidth=2, capsize=2, label=f"{t} thread{'s' if t > 1 else ''}")
    summary[f"ratio_cado_over_ours_core_t{t}"] = {str(d): {"median": statistics.median(ratios[d]),
        "min": min(ratios[d]), "max": max(ratios[d]), "n": len(ratios[d])} for d in xs}
ax.axhline(1, color="#8a8984", linewidth=1)
style(ax, "CADO bwc / ours (block Wiedemann), same unmerged core: >1 means ours is faster", "time ratio")
ax.legend(fontsize=8, frameon=False)
save(fig, "engine_ratio")

# 3. Core size after filtering.
fig, ax = plt.subplots(figsize=(9, 4.5))
for arm, key, rows, label, color, marker in (
    ("core-bw", "core_cols", grid, "unmerged core (singleton + clique only)", "#2a78d6", "o"),
    ("full-bw", None, grid, "ours: after merge (core_dimension)", "#1baf7a", "D"),
    ("cado-full", "small_dim", amend, "CADO merge-dl (+1 coeffs)", "#e87ba4", "X"),
    ("full-bw-pos", None, amend, "ours: after merge (+1 coeffs)", "#008300", "d"),
):
    by = {}
    for r in rows:
        if r.get("arm") != arm or r.get("failed") or r.get("weight", 3) != 3:
            continue
        if rows is grid and r.get("pass") != "main":
            continue
        v = r[key] if key else r["report"]["core_dimension"]
        by.setdefault(r["cols"], []).append(v)
    if by:
        xs = sorted(by)
        ax.plot(xs, [statistics.median(by[d]) for d in xs], color=color, marker=marker,
                linewidth=2, label=label)
        summary[f"core_dim|{arm}"] = {str(d): statistics.median(by[d]) for d in xs}
ax.plot([1024, 262144], [1024, 262144], color="#8a8984", linewidth=1, linestyle=":", label="D (no filtering)")
style(ax, "Matrix handed to the Krylov solver (median over seeds)", "core dimension")
ax.legend(fontsize=8, frameon=False)
save(fig, "core_dimension")

# 4. End-to-end ic logs: collection against linear algebra.
if e2e:
    fig, ax = plt.subplots(figsize=(9, 4.5))
    labels, col, la = [], [], []
    for r in e2e:
        if r.get("failed"):
            continue
        labels.append(r["label"])
        col.append(r["collection_seconds"])
        la.append(r["linear_algebra_seconds"])
    y = range(len(labels))
    ax.barh(y, col, color="#2a78d6", label="relation collection", height=0.6)
    ax.barh(y, la, left=col, color="#eb6834", label="linear algebra", height=0.6)
    ax.set_yticks(list(y), labels, fontsize=7)
    ax.set_xscale("log")
    ax.set_xlabel("wall seconds (log)")
    ax.set_title("ic logs end to end: where the precompute time goes", loc="left", fontsize=11)
    ax.legend(fontsize=8, frameon=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    save(fig, "e2e_split")

(here / "summary.json").write_text(json.dumps(summary, indent=1, sort_keys=True))
print("figures:", sorted(p.name for p in figs.iterdir()))

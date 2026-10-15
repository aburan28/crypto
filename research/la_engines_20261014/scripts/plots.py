"""Lay out recorded values from raw/*.jsonl as SVG + 150 dpi PNG.

Layout only (AGENTS.md "Implementation language"): every plotted number is
read from the run records; medians, min-max over seeds and per-seed ratios
are the only arithmetic, and summary.json records them for the report.
"""
import json, statistics, pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = pathlib.Path(__file__).resolve().parent.parent
raw = here / "raw"
figs = here / "figures"
figs.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#8a8984",
                     "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e",
                     "ytick.color": "#52514e", "svg.fonttype": "none"})

def load(name):
    p = raw / name
    return [json.loads(l) for l in p.read_text().splitlines() if l.strip()] if p.exists() else []

grid, a1, a2, e2e = load("grid.jsonl"), load("grid_amend1.jsonl"), load("grid_amend2.jsonl"), load("e2e.jsonl")
a4 = load("grid_amend4.jsonl")
a5 = load("grid_amend5.jsonl")

# Amendment 3: a record whose run overlapped a build or diagnostic window is
# replaced by its idle-host re-run (matched by pass, arm, D, m, seed; the
# contended record is identified by its end time). Originals stay in raw/.
contended = load("contended_runs.jsonl")
reruns = load("rerun_contended.jsonl")
replaced = []
for c, rr in zip(contended, reruns):
    if rr.get("failed"):
        rr = dict(rr, cols=c["cols"], seed=c["seed"], weight=c["weight"])
    assert (c["pass"], c["arm"], c["cols"], c["seed"]) == (rr["pass"], rr["arm"], rr["cols"], rr["seed"]), (c, rr)
    for rows in (grid, a1, a2):
        for i, r in enumerate(rows):
            if r.get("utc") == c["utc"] and r.get("arm") == c["arm"] and r.get("cols") == c["cols"] \
                    and r.get("seed") == c["seed"] and r.get("pass") == c["pass"]:
                rows[i] = dict(rr, weight=rr.get("weight", c["weight"]), replaced_utc=c["utc"])
                replaced.append({"pass": c["pass"], "arm": c["arm"], "cols": c["cols"], "weight": c["weight"],
                                 "seed": c["seed"], "original_s": r.get("solve_s"), "rerun_s": rr.get("solve_s"),
                                 "rerun_verified": rr.get("verified")})

# The grid driver kept the bench's own arm name; the A/A repeat is the second
# `core-bw` record for (main, D, seed 1). Relabel it so it never enters a median.
seen = set()
for r in grid:
    if r.get("pass") == "main" and r.get("arm") == "core-bw" and r.get("seed") == 1:
        k = r["cols"]
        if k in seen:
            r["arm"] = "core-bw-aa"
        seen.add(k)

# Fixed categorical order (reference palette, light mode); marker shapes carry
# identity too, and every chart has a legend.
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SERIES = {
    "igauss": ("IncrementalGauss: dense, rows fed one by one (ecbench default)", C[7], "v"),
    "sgauss": ("StructuredGauss: sparse RREF", C[3], "^"),
    "full-bw": ("ours: filter + merge + block Wiedemann", C[0], "o"),
    "full-bl": ("ours: filter + merge + block Lanczos", C[6], "P"),
    "core-bw": ("ours: block Wiedemann on unmerged core", C[2], "D"),
    "cado-bwc": ("CADO bwc on the same unmerged core", C[1], "s"),
    "full-bw-pos": ("ours, default merge (8, 32)", C[0], "o"),
    "full-bw-m32r100": ("ours, merge (32, 100)", C[2], "D"),
    "full-bw-m128r200": ("ours, merge (128, 200)", C[6], "P"),
    "cado-full-fb": ("CADO purge + merge-dl + replay-dl + bwc", C[1], "s"),
    "legacy-bw": ("ours before: legacy filter + block Wiedemann", C[7], "v"),
    "cado-bw": ("ours after: CADO-style filter + block Wiedemann", C[0], "o"),
    "cado-bl": ("ours after: CADO-style filter + block Lanczos", C[6], "P"),
    "cado-nfs": ("CADO-NFS purge + merge-dl + bwc (+1 coeffs)", C[1], "s"),
    "cado-bw-pos": ("ours after (+1 coeffs)", C[2], "D"),
    "legacy-bw-pos": ("ours before (+1 coeffs)", C[3], "^"),
}
summary = {}

def cells(rows, arm, key="solve_s", pass_=None):
    by = {}
    for r in rows:
        if r.get("arm") != arm or r.get("failed") or r.get("verified") is not True:
            continue
        if r.get("weight", 3) != 3 or (pass_ and r.get("pass") != pass_):
            continue
        v = r[key] if not callable(key) else key(r)
        by.setdefault(r["cols"], []).append(v)
    return {d: (statistics.median(v), min(v), max(v), len(v)) for d, v in sorted(by.items())}

def style(ax, title, ylabel, ylog=True):
    ax.set_xscale("log", base=2)
    if ylog:
        ax.set_yscale("log")
    ax.set_xlabel("factor-base columns D (synthetic, m = 3 summands per relation)")
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontsize=10, color="#0b0b0b")
    ax.grid(True, which="major", color="#e4e3df", linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

def save(fig, name):
    fig.savefig(figs / f"{name}.svg", bbox_inches="tight")
    fig.savefig(figs / f"{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

def lines(ax, rows, arms, tag, pass_=None, key="solve_s"):
    for arm in arms:
        c = cells(rows, arm, key=key, pass_=pass_)
        if not c:
            continue
        label, color, marker = SERIES[arm]
        xs = list(c)
        med = [c[d][0] for d in xs]
        ax.errorbar(xs, med, yerr=[[c[d][0] - c[d][1] for d in xs], [c[d][2] - c[d][0] for d in xs]],
                    color=color, marker=marker, markersize=6, linewidth=2, capsize=2, label=label,
                    markeredgecolor="#fcfcfb", markeredgewidth=1)
        summary[f"{tag}|{arm}"] = {str(d): {"median": v[0], "min": v[1], "max": v[2], "n": v[3]} for d, v in c.items()}

# 1. Every solver on the frozen grid (±1 coefficients, four threads).
fig, ax = plt.subplots(figsize=(8.5, 5))
lines(ax, grid, ["igauss", "sgauss", "core-bw", "cado-bwc", "full-bw", "full-bl"], "main_solve_s", pass_="main")
ax.annotate("1,800 s timeout: both unmerged-core arms at D = 2^18", (262144, 1800), color="#52514e",
            fontsize=7, ha="right", va="bottom", xytext=(0, 3), textcoords="offset points")
ax.axhline(1800, color="#c3c2b7", linewidth=0.8, linestyle=":")
style(ax, "Solve wall time against D: median of 3 seeds, bars min-max, 4 threads", "wall seconds (log)")
ax.legend(fontsize=7.5, frameon=False, loc="upper left")
save(fig, "solve_time_vs_D")

# 2. Engine against engine on the identical unmerged core, per-seed ratio.
fig, ax = plt.subplots(figsize=(8.5, 4))
for t, pass_, color, marker in ((4, "main", C[0], "o"), (1, "onethread", C[1], "s")):
    a, b = {}, {}
    for r in grid:
        if r.get("pass") != pass_ or r.get("failed") or r.get("verified") is not True or r.get("weight", 3) != 3:
            continue
        (a if r["arm"] == "core-bw" else b if r["arm"] == "cado-bwc" else {})[(r["cols"], r["seed"])] = r["solve_s"]
    ratios = {}
    for k in a:
        if k in b:
            ratios.setdefault(k[0], []).append(b[k] / a[k])
    xs = sorted(ratios)
    med = [statistics.median(ratios[d]) for d in xs]
    ax.errorbar(xs, med, yerr=[[m - min(ratios[d]) for m, d in zip(med, xs)], [max(ratios[d]) - m for m, d in zip(med, xs)]],
                color=color, marker=marker, linewidth=2, capsize=2, label=f"{t} thread{'s' if t > 1 else ''} (CADO thr={'2x2' if t > 1 else '1x1'})",
                markeredgecolor="#fcfcfb", markeredgewidth=1)
    summary[f"ratio_cado_bwc_over_core_bw_t{t}"] = {str(d): {"median": statistics.median(ratios[d]), "min": min(ratios[d]), "max": max(ratios[d]), "n": len(ratios[d])} for d in xs}
ax.axhline(1, color="#8a8984", linewidth=1)
style(ax, "CADO bwc time / our block Wiedemann time, same matrix and rhs (above 1: ours faster)", "ratio (log)")
ax.legend(fontsize=8, frameon=False)
save(fig, "engine_ratio")

# 3. Pipeline against pipeline (+1 coefficients, amendments 1-2).
fig, ax = plt.subplots(figsize=(8.5, 4.5))
lines(ax, a1, ["full-bw-pos"], "pipeline_solve_s")
lines(ax, a2, ["full-bw-m32r100", "full-bw-m128r200", "cado-full-fb"], "pipeline_solve_s")
style(ax, "Whole linear-algebra pipeline, filter to verified logs (+1 coefficients, 4 threads)", "wall seconds (log)")
ax.legend(fontsize=8, frameon=False, loc="upper left")
save(fig, "pipeline_time")

# 4. What each filter hands its Krylov solver.
fig, ax = plt.subplots(figsize=(8.5, 4.5))
core = lambda r: r.get("small_dim") or r["report"]["core_dimension"]
for rows, arm in ((a1, "full-bw-pos"), (a2, "full-bw-m32r100"), (a2, "full-bw-m128r200"), (a2, "cado-full-fb")):
    c = cells(rows, arm, key=core)
    label, color, marker = SERIES[arm]
    xs = list(c)
    ax.plot(xs, [c[d][0] for d in xs], color=color, marker=marker, linewidth=2, label=label,
            markeredgecolor="#fcfcfb", markeredgewidth=1)
    summary[f"core_dim|{arm}"] = {str(d): c[d][0] for d in xs}
ax.plot([1024, 262144], [1024, 262144], color="#8a8984", linewidth=1, linestyle=":", label="D itself (no filtering)")
style(ax, "Matrix dimension after filtering and merging (median of 3 seeds)", "core dimension (log)")
ax.legend(fontsize=8, frameon=False, loc="upper left")
save(fig, "core_dimension")

# 5. Real ic logs: linear algebra against everything else (derived).
if e2e:
    rows = [r for r in e2e if not r.get("failed") and r.get("verified")]
    labels, other, la = [], [], []
    for r in rows:
        if r["seed"] != 1:
            continue
        mode = r["linear_algebra"] + ("" if r["linear_algebra"] == "dense" else f" {'BW' if 'wiedemann' in r['sparse_solver'] else 'BL'}")
        labels.append(f"{r['label']} ({r['columns']} cols), {mode}")
        other.append(r["other_seconds_derived"] * 1e3)
        la.append(r["linear_algebra_seconds"] * 1e3)
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    y = list(range(len(labels)))
    ax.barh([v - 0.2 for v in y], other, height=0.38, color=C[0], label="collection + certification (derived: elapsed − LA)")
    ax.barh([v + 0.2 for v in y], la, height=0.38, color=C[1], label="linear algebra (report's linear_algebra.seconds)")
    ax.set_yticks(y, labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("milliseconds (log), seed 1")
    ax.set_title("ic logs end to end on Koblitz curves: the linear algebra is 0.1–1.9% of the precompute", loc="left", fontsize=10)
    ax.grid(True, axis="x", color="#e4e3df", linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(fontsize=7.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
    save(fig, "e2e_split")
    summary["e2e"] = [{k: r[k] for k in ("label", "linear_algebra", "sparse_solver", "seed", "columns", "relations",
                       "elapsed_seconds", "linear_algebra_seconds", "other_seconds_derived")} for r in rows]

# 6. Amendments 4 and 5: before and after, and CADO-NFS, each in one session.
def before_after(rows, tag, title):
    fig, ax = plt.subplots(figsize=(8.5, 5))
    lines(ax, rows, ["legacy-bw", "cado-bw", "cado-bl", "legacy-bw-pos", "cado-bw-pos", "cado-nfs"], f"{tag}_solve_s")
    style(ax, title, "wall seconds (log)")
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")
    save(fig, f"{tag}_pipeline")

    fig, ax = plt.subplots(figsize=(8.5, 4))
    for num, den, label, color, marker in (("legacy-bw", "cado-bw", "legacy / CADO-style filter, ±1 coefficients", C[0], "o"),
                                           ("legacy-bw-pos", "cado-bw-pos", "legacy / CADO-style filter, +1 coefficients", C[2], "D"),
                                           ("cado-nfs", "cado-bw-pos", "CADO-NFS / ours after, +1 coefficients", C[1], "s")):
        a = {(r["cols"], r["seed"]): r["solve_s"] for r in rows if r.get("arm") == num and r.get("verified") is True}
        b = {(r["cols"], r["seed"]): r["solve_s"] for r in rows if r.get("arm") == den and r.get("verified") is True}
        ratios = {}
        for k in a:
            if k in b:
                ratios.setdefault(k[0], []).append(a[k] / b[k])
        xs = sorted(ratios)
        if not xs:
            continue
        med = [statistics.median(ratios[d]) for d in xs]
        ax.errorbar(xs, med, yerr=[[m - min(ratios[d]) for m, d in zip(med, xs)], [max(ratios[d]) - m for m, d in zip(med, xs)]],
                    color=color, marker=marker, linewidth=2, capsize=2, label=label, markeredgecolor="#fcfcfb", markeredgewidth=1)
        summary[f"{tag}_ratio|{num}/{den}"] = {str(d): {"median": statistics.median(ratios[d]), "min": min(ratios[d]),
                                                         "max": max(ratios[d]), "n": len(ratios[d])} for d in xs}
    ax.axhline(1, color="#8a8984", linewidth=1)
    style(ax, "Time ratio per (D, seed) within one session (above 1: the second-named pipeline is faster)", "time ratio (log)")
    ax.legend(fontsize=8, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3)
    save(fig, f"{tag}_speedup")

    core = lambda r: r.get("small_dim") or r["report"]["core_dimension"]
    for arm in ("legacy-bw", "cado-bw", "legacy-bw-pos", "cado-bw-pos", "cado-nfs"):
        c = cells(rows, arm, key=core)
        summary[f"{tag}_core_dim|{arm}"] = {str(d): c[d][0] for d in c}
    if not any(r.get("arm", "").startswith("cado-bw-pos-d") for r in rows):
        return
    fig, ax = plt.subplots(figsize=(8.5, 4))
    for i, (arm, label) in enumerate((("cado-bw-pos-d0", "target density 0 (stop on columns × nonzeros)"),
                                      ("cado-bw-pos-d50", "target density 50"),
                                      ("cado-bw-pos", "target density 100 (shipped)"),
                                      ("cado-bw-pos-d170", "target density 170"))):
        c = cells(rows, arm)
        xs = list(c)
        if not xs:
            continue
        ax.errorbar(xs, [c[d][0] for d in xs], yerr=[[c[d][0] - c[d][1] for d in xs], [c[d][2] - c[d][0] for d in xs]],
                    color=[C[3], C[2], C[0], C[6]][i], marker="oDsP"[i], linewidth=2, capsize=2, label=label,
                    markeredgecolor="#fcfcfb", markeredgewidth=1)
        summary[f"{tag}_density|{arm}"] = {str(d): {"median": v[0], "min": v[1], "max": v[2], "n": v[3]} for d, v in c.items()}
    style(ax, "CADO-style filter: sensitivity to the merge density target (+1 coefficients)", "wall seconds (log)")
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    save(fig, f"{tag}_density")



if a4:
    before_after(a4, "amend4", "After the port: whole linear-algebra pipeline, one session (median of 3 seeds, 4 threads)")
if a5:
    before_after(a5, "amend5", "After the filter speed-ups: whole pipeline, one session (median of 3 seeds, 4 threads)")

summary["amendment3_replaced"] = replaced
(here / "summary.json").write_text(json.dumps(summary, indent=1, sort_keys=True))
print("figures:", sorted(p.name for p in figs.iterdir()))

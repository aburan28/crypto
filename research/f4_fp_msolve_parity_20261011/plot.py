# Lays out the recorded values of one run directory (argument): medians of
# the recorded rounds per arm, and their ratios.  No new computation.
import collections, glob, json, os, re, statistics as st, sys
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

run = os.path.abspath(sys.argv[1])
out = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else run
R = lambda *p: os.path.join(run, *p)
v = collections.defaultdict(list)
fps = collections.defaultdict(set)
for f in glob.glob(R("raw", "ours_*_t*_*.round*.jsonl")):
    arm, th = re.search(r"ours_(main|new)_t(\d)_", f).groups()
    for l in open(f):
        d = json.loads(l)
        v[(d["case"], f"f4_fp {arm}, {th} thr")].append(d["ours_f4_ms_median"])
        fps[(d["case"], arm)].add(d["basis_fp"])
for f in glob.glob(R("raw", "external.round*.tsv")):
    for l in open(f):
        p = l.rstrip("\n").split("\t")
        if p[1].startswith("msolve"):
            v[(p[0], f"msolve, {p[1][-1]} thr")].append(float(p[3]))
        elif p[2] == "ok":
            m = re.search(r"std_ms (\d+)", p[3])
            if m:
                v[(p[0], "Singular std")].append(float(m.group(1)))
med = {k: st.median(x) for k, x in v.items()}
json.dump({f"{k[0]}|{k[1]}": [med[k], len(v[k])] for k in sorted(med)},
          open(os.path.join(out, "medians.json"), "w"), indent=1)

cases = ["quad_n8_p65521", "quad_n9_p65521", "cubic_n5_p31", "overdet_n12_m24_p65521",
         "quad_n10_p65521", "quad_n11_p65521", "cubic_n6_p31"]
arms = [("f4_fp main, 1 thr", "#9aa3ad"), ("f4_fp new, 1 thr", "#2a6fdb"), ("msolve, 1 thr", "#d9822b"),
        ("f4_fp new, 4 thr", "#7fa8ee"), ("msolve, 4 thr", "#f0b97a"), ("Singular std", "#5a9e5a")]
for name, scale in (("fp_vs_msolve", "log"),):
    fig, ax = plt.subplots(figsize=(11, 4.6))
    w = 0.13
    for j, (a, c) in enumerate(arms):
        xs, ys = [], []
        for i, cs in enumerate(cases):
            if (cs, a) in med:
                xs.append(i + (j - 2.5) * w)
                ys.append(med[(cs, a)])
                ax.scatter([i + (j - 2.5) * w] * len(v[(cs, a)]), v[(cs, a)], color="black", s=3, zorder=3)
        ax.bar(xs, ys, width=w, color=c, label=a)
    ax.set_yscale(scale)
    ax.set_xticks(range(len(cases)))
    ax.set_xticklabels([c.replace("_p65521", "\np=65521").replace("_p31", "\np=31") for c in cases], fontsize=8)
    ax.set_ylabel("wall ms (log scale)")
    ax.legend(fontsize=7.5, ncol=6, loc="upper left")
    ax.set_title(f"Reduced grevlex Gröbner basis over F_p, same systems; bars = median of rounds, dots = rounds\n"
                 f"{os.path.basename(run)}; Singular missing where it hit the time limit", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(out, f"{name}.svg"))
    fig.savefig(os.path.join(out, f"{name}.png"), dpi=150)

allcases = sorted({k[0] for k in med}, key=lambda c: med.get((c, "f4_fp main, 1 thr"), 0))
fig, ax = plt.subplots(figsize=(11, 4.2))
xs = range(len(allcases))
for a, b, lab, c, mk in (("f4_fp main, 1 thr", "f4_fp new, 1 thr", "new vs main (1 thr)", "#2a6fdb", "o"),
                         ("msolve, 1 thr", "f4_fp new, 1 thr", "msolve ÷ new (1 thr)", "#d9822b", "s"),
                         ("msolve, 4 thr", "f4_fp new, 4 thr", "msolve ÷ new (4 thr)", "#f0b97a", "^")):
    ys = [med[(cs, a)] / med[(cs, b)] if (cs, a) in med and (cs, b) in med else float("nan") for cs in allcases]
    ax.plot(list(xs), ys, marker=mk, linestyle="none", color=c, label=lab)
ax.axhline(1, color="#888", lw=0.8)
ax.set_yscale("log")
ax.set_xticks(list(xs))
ax.set_xticklabels([c.replace("_p65521", "").replace("_p31", "").replace("_p29", "").replace("_p101", "") for c in allcases],
                   rotation=40, fontsize=7.5, ha="right")
ax.set_ylabel("ratio of medians (> 1: f4_fp new faster)")
ax.legend(fontsize=8)
ax.set_title(f"Speed ratios per system, {os.path.basename(run)} (ratios of recorded medians)", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(out, "ratios.svg"))
fig.savefig(os.path.join(out, "ratios.png"), dpi=150)
print(json.dumps({c: sorted(fps[(c, 'main')] | fps[(c, 'new')]) for c in sorted({k[0] for k in fps})}))

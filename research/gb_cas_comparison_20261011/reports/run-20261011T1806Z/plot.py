# Lays out recorded values from raw/ (medians of the three recorded rounds
# for F_p; per-target values for the Boolean roots).  No new computation.
import collections, glob, json, os, re, statistics as st
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
R = lambda *p: os.path.join(here, *p)

v = collections.defaultdict(list)
for f in glob.glob(R("raw", "fp_ours_t*.round*.jsonl")):
    th = re.search(r"_t(\d)\.", f).group(1)
    for l in open(f):
        d = json.loads(l)
        v[(d["case"], f"ours f4_fp, {th} thr")].append(d["ours_f4_ms_median"])
for f in glob.glob(R("raw", "fp_external.round*.tsv")):
    for l in open(f):
        p = l.rstrip("\n").split("\t")
        if p[1].startswith("msolve"):
            v[(p[0], f"msolve, {p[1][-1]} thr")].append(float(p[3]))
        else:
            v[(p[0], "Singular std")].append(float(re.search(r"std_ms (\d+)", p[2]).group(1)))
cases = ["quad_n7_p65521", "quad_n8_p65521", "quad_n9_p65521", "cubic_n5_p31",
         "overdet_n10_m15_p65521", "overdet_n12_m24_p65521"]
arms = ["ours f4_fp, 1 thr", "ours f4_fp, 4 thr", "msolve, 1 thr", "msolve, 4 thr", "Singular std"]
colors = ["#2a6fdb", "#7fa8ee", "#d9822b", "#f0b97a", "#5a9e5a"]
fig, ax = plt.subplots(figsize=(10, 4.4))
w = 0.16
for j, (a, c) in enumerate(zip(arms, colors)):
    xs = [i + (j - 2) * w for i in range(len(cases))]
    meds = [st.median(v[(cs, a)]) for cs in cases]
    ax.bar(xs, meds, width=w, color=c, label=a)
    for x, cs in zip(xs, cases):
        ax.scatter([x] * len(v[(cs, a)]), v[(cs, a)], color="black", s=4, zorder=3)
ax.set_yscale("log")
ax.set_xticks(range(len(cases)))
ax.set_xticklabels([c.replace("_p65521", "\np=65521").replace("_p31", "\np=31") for c in cases], fontsize=8)
ax.set_ylabel("wall ms (log scale)")
ax.legend(fontsize=8, ncol=5, loc="upper left")
ax.set_title("Reduced grevlex Gröbner basis over F_p, same systems; bars = median of 3 rounds, dots = rounds\n"
             "run-20261011T1806Z; msolve wall includes ~6 ms process start-up", fontsize=9)
fig.tight_layout()
fig.savefig(R("fp.svg"))

ours = {}
for l in open(R("raw", "bool_n17_ours.jsonl")):
    if '"target"' in l:
        d = json.loads(l)
        ours[(d["arm"], d["target"])] = d["ms"]
ext = collections.defaultdict(dict)
for l in open(R("raw", "bool_n17_external.tsv")):
    p = l.rstrip("\n").split("\t")
    k = int(p[0][4:])
    if p[1] == "msolve_t1":
        ext["msolve, 1 thr (full GB)"][k] = float(p[3])
    elif p[1] == "polybori":
        ext["PolyBoRi (full GB)"][k] = float(re.search(r"polybori_ms ([0-9.]+)", p[4]).group(1))
    elif p[1] == "singular_std":
        ext["Singular std: no result in 120 s"][k] = float(p[3])
fig, ax = plt.subplots(figsize=(10, 4.2))
ts = list(range(1, 9))
series = [("inherited F4 (tree)", [ours[("f4", t)] for t in ts], "#2a6fdb", "o"),
          ("F6-IC + index (tree)", [ours[("f6p", t)] for t in ts], "#7fa8ee", "s"),
          ("msolve, 1 thr (full GB)", [ext["msolve, 1 thr (full GB)"][t] for t in ts], "#d9822b", "^"),
          ("PolyBoRi (full GB)", [ext["PolyBoRi (full GB)"][t] for t in ts], "#a05a9e", "D"),
          ("Singular std: no result in 120 s", [ext["Singular std: no result in 120 s"][t] for t in ts], "#5a9e5a", "x")]
for name, ys, c, mk in series:
    ax.plot(ts, ys, marker=mk, color=c, label=name, linestyle="none", markersize=7)
ax.set_yscale("log")
ax.set_xlabel("target k (decomposition system of [k]G); k = 3 decomposes, the others are the unit ideal")
ax.set_ylabel("wall ms per target (log scale)")
ax.set_xticks(ts)
ax.legend(fontsize=8, loc="center right")
ax.set_title("Deciding the m = 3 decomposition system, K(1,17), dim-6 factor base (35 Boolean vars, 34 cubics)\n"
             "run-20261011T1806Z; ours: best of 3 reps; external: one run each", fontsize=9)
fig.tight_layout()
fig.savefig(R("bool_n17.svg"))

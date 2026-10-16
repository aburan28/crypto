"""Layout only: plots values recorded in ../../runs/*.jsonl (AGENTS.md: no
Python computes a result).  Panel A: committed XL degree vs measured closure
degree per draw.  Panel B: recorded closure word ops and matrix columns at the
closure's refuting degree, against the recorded unknown count N."""
import glob, json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
rows = []
for p in sorted(glob.glob(os.path.join(here, "../../runs/*.jsonl"))):
    if "satisfiable" in p:
        continue
    for line in open(p):
        if line.strip():
            rows.append(json.loads(line))

cells = sorted({(r["m"], r["n"], r["ell"], r["n_vars"]) for r in rows},
               key=lambda c: (c[0], c[3], c[2]))
fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.8))
for i, (m, n, ell, nv) in enumerate(cells):
    draws = [r for r in rows if (r["m"], r["n"], r["ell"]) == (m, n, ell)]
    for j, r in enumerate(draws):
        o = r["committed_outcome"]
        x = i + (j - 1.5) * 0.12
        a.scatter(x, o["degree"], marker="^" if o["kind"] == "resolved" else "2",
                  color="tab:red", s=28, zorder=3,
                  label="XL (committed)" if i == 0 and j == 0 else None)
        a.scatter(x, r["closure_degree"], marker="o", color="tab:blue", s=22, zorder=3,
                  label="MutantXL closure (measured)" if i == 0 and j == 0 else None)
a.scatter([], [], marker="2", color="tab:red", label="XL lower bound (≥D, m=4 ℓ=3)")
a.set_xticks(range(len(cells)))
a.set_xticklabels([f"m{m} ({n},{ell})" for m, n, ell, _ in cells], rotation=75, fontsize=7)
a.set_ylabel("refutation degree")
a.set_ylim(2.5, 7.5)
a.grid(axis="y", alpha=0.3)
a.legend(fontsize=7, loc="upper left")
a.set_title("A. Refutation degree, 4 unsat draws per cell (92 draws)", fontsize=9)

for m, col in [(3, "tab:green"), (4, "tab:purple")]:
    pts = []
    for r in rows:
        if r["m"] != m:
            continue
        s = next(s for s in r["scan"] if s["degree"] == r["closure_degree"])
        pts.append((r["n_vars"], s["closure_word_ops"], s["columns"]))
    b.scatter([p[0] for p in pts], [p[1] for p in pts], color=col, s=16,
              label=f"m={m}: closure word ops (64-bit)")
    b.scatter([p[0] for p in pts], [p[2] for p in pts], color=col, s=16, marker="x",
              label=f"m={m}: columns at closure degree")
b.set_yscale("log", base=2)
b.set_xlabel("unknowns N")
b.set_ylabel("count (log2 scale)")
b.grid(alpha=0.3)
b.legend(fontsize=7)
b.set_title("B. Recorded cost at the closure's refuting degree", fontsize=9)
fig.suptitle("XL vs MutantXL closure on committed chained-S3 draws — "
             "research/dreg_closure_mutantxl_20261016, run 2026-10-16T02:35Z", fontsize=9)
fig.tight_layout()
for ext in ("svg",):
    fig.savefig(os.path.join(here, f"closure_vs_xl.{ext}"))
fig.savefig(os.path.join(here, "closure_vs_xl.png"), dpi=150)

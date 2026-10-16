"""Layout only: plots values recorded in ../../runs/*.d*.jsonl (AGENTS.md: no
Python computes a result).  A: closure refutation degree per draw against the
unknowns N (a draw whose degree-4 closure saturated without refuting is drawn
at '>4').  B: recorded closure word ops against recorded columns at the
refuting degree."""
import glob, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
rows = []
for p in sorted(glob.glob(os.path.join(here, "../../runs/*.d*.jsonl"))):
    if os.path.basename(p).startswith("control-"):
        continue
    t = open(p).read().strip()
    if t:
        rows.append(json.loads(t.splitlines()[0]))

style = {3: ("tab:green", "o"), 4: ("tab:purple", "s")}
fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.8))
seen = set()
for r in rows:
    col, mk = style[r["m"]]
    cd = r["closure_degree"]
    # m = 3 drawn just below the row, m = 4 just above, so shared N stay visible.
    y = (cd if cd is not None else 5) + (-0.08 if r["m"] == 3 else 0.08)
    jitter = (r["draw"] % 4 - 1.5) * 0.18
    lab = f"m={r['m']}" + (" (ℓ shown)" if False else "")
    a.scatter(r["n_vars"] + jitter, y, color=col, marker=mk if cd is not None else "^",
              s=40, label=lab if (r["m"], cd is None) not in seen and cd is not None else None)
    seen.add((r["m"], cd is None))
    if r["draw"] == min(q["draw"] for q in rows if (q["m"], q["n"], q["ell"]) == (r["m"], r["n"], r["ell"])):
        a.annotate(f"m{r['m']} ℓ={r['ell']}", (r["n_vars"], y + (-0.22 if r["m"] == 3 else 0.14)),
                   fontsize=7, ha="center", color=col)
a.scatter([], [], color="tab:green", marker="^", s=40, label="m=3: degree 4 saturated, no refutation (>4)")
a.set_yticks([3, 4, 5]); a.set_yticklabels(["3", "4", ">4"])
a.set_ylim(2.6, 5.5)
a.set_xlabel("unknowns N")
a.set_ylabel("closure refutation degree")
a.grid(alpha=0.3)
a.legend(fontsize=7, loc="upper left")
a.set_title("A. Closure degree per draw (40 unsat draws, 10 cells)", fontsize=9)

for m in (3, 4):
    col, mk = style[m]
    pts = []
    for r in rows:
        if r["m"] != m or r["closure_degree"] is None:
            continue
        s = next(s for s in r["scan"] if s["degree"] == r["closure_degree"])
        pts.append((s["columns"], s["closure_word_ops"]))
    b.scatter([p[0] for p in pts], [p[1] for p in pts], color=col, marker=mk, s=18, label=f"m={m}")
b.set_xscale("log", base=2); b.set_yscale("log", base=2)
b.set_xlabel("columns C(N, ≤4) at the refuting degree")
b.set_ylabel("counted 64-bit word ops (closure build)")
b.grid(alpha=0.3); b.legend(fontsize=7)
b.set_title("B. Recorded closure work at degree 4 (32 refuting draws)", fontsize=9)
fig.suptitle("MutantXL closure at scale — research/dreg_closure_scale_20261016, run 2026-10-16T17:05Z", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(here, "closure_scale.svg"))
fig.savefig(os.path.join(here, "closure_scale.png"), dpi=150)

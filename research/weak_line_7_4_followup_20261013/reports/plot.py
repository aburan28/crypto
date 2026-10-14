#!/usr/bin/env python3
"""Layout only: plots recorded values from results/f3.jsonl, f4.jsonl, f5.jsonl (no computation of results)."""
import json, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

res = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "weak_line_extra")
rd = lambda n: [json.loads(l) for l in open(os.path.join(res, n))]
f3 = [r for r in rd("f3.jsonl") if r["part"] == "field"]
f4 = {(r["q"], r["n"]): r for r in rd("f4.jsonl")}
f5 = rd("f5.jsonl")

def parse4(r):
    return {int(s.split(":")[0]): int(s.split(":")[1]) for s in r["hist_lines"]}

def parse3(r):
    h = {}
    for s in r["hist_census_lines_autsize"]:
        k, v = s.split(":"); lines = int(k.strip("()").split(",")[0]); h[lines] = h.get(lines, 0) + int(v)
    return h

fig, (ax, bx) = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"width_ratios": [1.6, 1]})
order = sorted([(r["n"] % 2, r["q"], r["n"]) for r in f3], key=lambda t: (t[0], t[2], t[1]))
labels = []; xs = []
cols = {2: "#1f77b4", 4: "#ff7f0e", 5: "#2ca02c"}
for i, (_, q, n) in enumerate(order):
    r3 = next(r for r in f3 if r["q"] == q and r["n"] == n); h3 = parse3(r3); tot3 = r3["curves"]
    labels.append(f"({q},{n})\nN={tot3}")
    if (q, n) in f4:
        h4 = parse4(f4[(q, n)]); tot4 = f4[(q, n)]["g_total"]; bottom = 0
        for k in (2, 4, 5):
            v = h4.get(k, 0) / tot4
            if v: ax.bar(i, v, bottom=bottom, color=cols[k], width=0.6, label=f"exact, {k} lines" if i == 0 or True else None); bottom += v
    extra = sum(v for k, v in h3.items() if k > 1) / tot3
    ax.plot(i, extra, "kD", ms=5, label="control sample (census)")
hand, lab = ax.get_legend_handles_labels(); seen = {}
for h, l in zip(hand, lab): seen.setdefault(l, h)
ax.legend(seen.values(), seen.keys(), fontsize=8, loc="center right")
ax.set_xticks(range(len(order))); ax.set_xticklabels(labels, fontsize=7)
ax.set_ylabel("fraction with more than one weak line")
ax.set_title("A. Extra weak lines: exact population over all g (bars)\nand census of fresh weak curves (diamonds); odd n on the right", fontsize=9)
ax.axvline(sum(1 for o in order if o[0] == 0) - 0.5, color="grey", lw=0.8, ls=":")
ax.text(0.01, 0.80, "bars: Lemma-A tally over all monic squarefree g, nonzero lower coeffs;\nno bar: field not enumerated, or fraction 0", transform=ax.transAxes, fontsize=7, va="top")

names = []; lines_obs = []
for r in f5:
    names.append(f"({r['q']},{r['n']})\n{r['arm'].replace('x5+c-supp','x5+c*')}")
    lines_obs.append(max(int(s.split(":")[0].strip("()").split(",")[0]) for s in r["hist_census_lines_autsize"]))
bx.bar(range(len(f5)), lines_obs, color=["#ff7f0e" if r["arm"] == "odd" else "#9467bd" for r in f5])
for i, r in enumerate(f5):
    bx.text(i, lines_obs[i] + 0.1, f"{r['lemmaA_equals_census']}/{r['curves']}", ha="center", fontsize=7)
bx.set_xticks(range(len(f5))); bx.set_xticklabels(names, fontsize=6, rotation=0)
bx.set_ylabel("weak lines found by census")
bx.set_title("B. Positive controls: odd g (orange), x^5+c over F_11 (purple)\nlabel = curves where Lemma-A count equals census", fontsize=9)
bx.set_ylim(0, 6.2)
fig.suptitle("Second weak lines follow affine symmetries of g with norm 1 (runs: results/f3, f4, f5.jsonl; seed 20261012 streams 7000+, 8000+)", fontsize=9)
fig.tight_layout()
fig.savefig(out + ".svg"); fig.savefig(out + ".png", dpi=150)
print("wrote", out)

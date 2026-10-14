"""Plot recorded values only: results.json -> cost_vs_genus.svg (no computation)."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = json.load(open("results.json"))
rho = d["log2_rho_n130"]
sw = d["sweep_n130_l1"]
g = [p["genus"] for p in sw]
c = [p["log2_cost"] for p in sw]
rows = {(t["field_degree"], t["rule"]): t["rows"] for t in d["tables"]}
old = next(r for r in rows[(130, "legacy")] if r["subfield_degree_l"] == 1)
new = next(r for r in rows[(130, "mmt")] if r["subfield_degree_l"] == 1)

fig, ax = plt.subplots(figsize=(7.2, 4.2))
ax.plot(g, c, "-", color="#1f4e79", lw=1.6, label="model cost, genus g over $F_2$ (random-polynomial counts)")
ax.axhline(rho, color="#b03a2e", ls="--", lw=1.3, label=f"Pollard rho, $N=130$: $2^{{{rho:.2f}}}$")
ax.axhline(rho - 0, color="none")
ax.plot([int(old["genus"])], [old["log2_cost"]], "o", color="#7f7f7f", ms=8, label=f"pre-fix rule: $m=8$, $g=128$: $2^{{{old['log2_cost']:.2f}}}$")
ax.plot([int(new["genus"])], [new["log2_cost"]], "s", color="#1e8449", ms=8, label=f"MMT rule: $m=9$ (type I), $g=255$: $2^{{{new['log2_cost']:.2f}}}$")
ax.annotate(f"margin {rho-old['log2_cost']:.2f} bits", (128, old["log2_cost"]), (135, 42), fontsize=8, color="#555",
            arrowprops=dict(arrowstyle="-", color="#777", lw=0.6))
ax.annotate(f"margin {rho-new['log2_cost']:.2f} bits", (255, new["log2_cost"]), (180, 66), fontsize=8, color="#1e8449",
            arrowprops=dict(arrowstyle="-", color="#1e8449", lw=0.6))
ax.set_xlabel("genus g of the descended curve")
ax.set_ylabel("log$_2$ cost (model units)")
ax.set_title("N = 130, l = 1: index-calculus model cost vs genus (n = %d points)" % len(sw), fontsize=10)
ax.set_ylim(30, 72)
ax.grid(alpha=0.25)
ax.legend(fontsize=7.5, loc="lower right")
fig.text(0.01, 0.005, "source: results.json (ghs_counterfactual/v1); model figures for a hypothetical curve, not measurements",
         fontsize=6.5, color="#555")
fig.tight_layout(rect=(0, 0.02, 1, 1))
fig.savefig("cost_vs_genus.svg")
fig.savefig("cost_vs_genus.png", dpi=150)

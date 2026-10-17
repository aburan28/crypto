# Joins Magma's RESULT lines (results/magma_output.txt, pasted from the
# calculator) with our recorded solves (ours_*.jsonl) and lays them out:
# results/table.json, results/table.md and figures/magma_vs_ours.{svg,png}.
# Recorded values only.  Exactness is checked here, not assumed: a Magma
# verdict that disagrees with ours is reported as a disagreement.
import json, os, re, sys
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "results", "magma_output.txt")
magma = {}
for line in open(src):
    m = re.match(r"\s*(RESULT-GB|RESULT-VAR)\s+(.*)", line)
    if not m:
        continue
    kv = dict(x.split("=", 1) for x in m.group(2).split())
    magma.setdefault(kv["name"], {}).update({m.group(1): kv})

ours = {}
for l in open(os.path.join(here, "ours_k4.jsonl")):
    r = json.loads(l)
    ours[r["name"]] = {"family": "k=4 S5" if r["name"].startswith("s5") else f"k=4 generic d={r['equation_degree']}",
                       "ours_s": r["ours_total_ms"] / 1e3, "ours_engine": "gaudry_quartic (Macaulay)",
                       "expect": {"solutions": r["ours_solutions"]}, "size": f"dim {r['ours_dim']}"}
for l in open(os.path.join(here, "ours_d4.jsonl")):
    r = json.loads(l)
    ours[r["name"]] = {"family": "k=4 D4", "ours_s": r["ours_total_ms"] / 1e3, "ours_engine": "gaudry_quartic_d4 (Macaulay, 5 unknowns)",
                       "expect": {"solutions": r["ours_solutions"]}, "size": f"dim {r['ours_dim']}"}
for l in open(os.path.join(here, "ours_pkm.jsonl")):
    r = json.loads(l)
    ours[r["name"]] = {"family": f"PKM m={r['m']}", "ours_s": r["ours_ms"] / 1e3, "ours_engine": r["ours_engine"],
                       "expect": {"unit": r["ours_inconsistent"]}, "size": f"N={r['N']}"}
for l in open(os.path.join(here, "ours_bool.jsonl")):
    r = json.loads(l)
    if r.get("arm") == "f4" and "ms" in r:
        ours[f"bool_n17_root{r['target']}"] = {"family": "binary n=17", "ours_s": r["ms"] / 1e3,
                                               "ours_engine": "inherited F4 (decision)", "expect": {}, "size": "35 vars"}

rows = []
for name, o in sorted(ours.items()):
    mg = magma.get(name, {})
    gb = mg.get("RESULT-GB")
    var = mg.get("RESULT-VAR")
    row = {"name": name, **{k: v for k, v in o.items() if k != "expect"}, "magma_gb_s": None,
           "magma_variety_s": None, "verdict": "not run"}
    if gb:
        row["magma_gb_s"] = float(gb["magma_gb_s"])
        ok = True
        if "unit" in o["expect"]:
            ok = (gb["unit"] == "true") == bool(o["expect"]["unit"])
        if var:
            row["magma_variety_s"] = float(var["magma_variety_s"])
            if "agree" in var:
                ok = ok and var["agree"] == "true"
        elif "solutions" in o["expect"]:
            ok = None  # the basis finished, the variety did not (cut-off)
        row["verdict"] = {True: "agree", False: "DISAGREE", None: "basis only"}[ok]
        t = row["magma_gb_s"] + (row["magma_variety_s"] or 0.0)
        row["magma_total_s"] = t
        row["ours_over_magma"] = o["ours_s"] / max(t, 1e-3)
    rows.append(row)
json.dump(rows, open(os.path.join(here, "results", "table.json"), "w"), indent=1)

with open(os.path.join(here, "results", "table.md"), "w") as f:
    f.write("| system | family | size | ours (s) | ours engine | Magma GB (s) | Magma variety (s) | ours ÷ Magma | verdict |\n")
    f.write("|:--|:--|:--|--:|:--|--:|--:|--:|:--|\n")
    for r in rows:
        g = "—" if r["magma_gb_s"] is None else f"{r['magma_gb_s']:.3f}"
        v = "—" if r["magma_variety_s"] is None else f"{r['magma_variety_s']:.3f}"
        q = "—" if "ours_over_magma" not in r else f"{r['ours_over_magma']:.3g}"
        f.write(f"| `{r['name']}` | {r['family']} | {r['size']} | {r['ours_s']:.3f} | {r['ours_engine']} | {g} | {v} | {q} | {r['verdict']} |\n")
print(open(os.path.join(here, "results", "table.md")).read())

done = [r for r in rows if "ours_over_magma" in r]
if done:
    C = {"k=4": "#2a78d6", "PKM": "#eb6834", "binary": "#1baf7a"}
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for fam, col in C.items():
        pts = [r for r in done if r["family"].startswith(fam)]
        if pts:
            ax.scatter([r["ours_s"] for r in pts], [r["magma_total_s"] for r in pts], color=col, label=fam, s=26)
    lo = min(min(r["ours_s"] for r in done), min(r["magma_total_s"] for r in done)) / 2
    hi = max(max(r["ours_s"] for r in done), max(r["magma_total_s"] for r in done)) * 2
    ax.plot([lo, hi], [lo, hi], color="#9aa3ad", lw=0.8, ls="--", label="equal time")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("ours, s (this repository's engine, its recorded run)")
    ax.set_ylabel("Magma F4 GroebnerBasis (+ Variety), s")
    ax.grid(color="#e6e5e1", lw=0.6)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    for ext, kw in (("svg", {}), ("png", {"dpi": 150})):
        fig.savefig(os.path.join(here, "figures", f"magma_vs_ours.{ext}"), **kw)

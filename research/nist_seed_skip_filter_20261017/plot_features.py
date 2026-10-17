"""Plotting only: lays out recorded values from out/features.jsonl and
out/analysis.json (written by features.gp / analysis.gp). Computes nothing."""
import json, pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = pathlib.Path(__file__).parent
rows = [json.loads(l) for l in (here / "out/features.jsonl").read_text().splitlines()]
an = json.loads((here / "out/analysis.json").read_text())
P_C, S_C = "#2a6fdb", "#e0782f"

fig, (a, b) = plt.subplots(1, 2, figsize=(12, 6.2), gridspec_kw={"width_ratios": [1, 1.25]})
keys = [("lpf_twist", "twist order"), ("lpf_nm1", "n − 1"), ("lpf_np1", "n + 1")]
for k, (key, name) in enumerate(keys):
    for r in rows:
        lo, hi, full = r[key]
        x = k + (-0.13 if r["label"] == "P" else 0.13)
        c = P_C if r["label"] == "P" else S_C
        if full:
            a.scatter(x, lo / r["p_bits"], color=c, s=28, zorder=3,
                      marker="o" if r["label"] == "P" else "s")
        else:  # not fully factored in 600 s: recorded interval
            a.plot([x, x], [lo / r["p_bits"], hi / r["p_bits"]], color=c, lw=1.2, alpha=.6)
a.set_xticks(range(3), [n for _, n in keys])
a.set_ylabel("bits of largest prime factor / bits of p")
a.set_ylim(0, 1.05)
a.scatter([], [], color=P_C, marker="o", label="published (10)")
a.scatter([], [], color=S_C, marker="s", label="skipped (11)")
a.plot([], [], color="grey", label="interval: not fully factored")
a.legend(loc="lower right", fontsize=8, frameon=False)
a.set_title("(a) H8/H9 largest prime factors (21 accepted candidates)", fontsize=10)
cont = {c["feature"]: c for c in an["continuous"]}
a.text(0.01, 0.99, "exact Mann–Whitney p (2-sided): " + ", ".join(
    f"{n} {cont[f]['p_exact']:.2f}" for n, f in
    [("twist", "H8 lpf(twist)/pbits"), ("n−1", "H9 lpf(n-1)/pbits"), ("n+1", "H9 lpf(n+1)/pbits")]),
    transform=a.transAxes, fontsize=7.5, va="top")

bins = an["binary"]
y = range(len(bins))
b.barh([i - 0.2 for i in y], [f["P_true"] / 10 for f in bins], height=0.4, color=P_C, label="published: share true (n = 10)")
b.barh([i + 0.2 for i in y], [f["S_true"] / 11 for f in bins], height=0.4, color=S_C, label="skipped: share true (n = 11)")
b.set_yticks(list(y), [f["feature"] + (" *" if f["undecided"] else "") for f in bins], fontsize=7.5)
b.invert_yaxis()
b.set_xlim(0, 1)
b.set_xlabel("share of group for which the feature is true")
b.set_title("(b) 21 binary features: a perfect split needs 1.0 vs 0.0 — none does", fontsize=10)
b.legend(loc="lower right", fontsize=7.5, frameon=False)
nc = an["negative_control"]
fig.text(0.5, 0.005, f"* some rows undecided (not fully factored). Positive control (planted label) detected: "
         f"{bool(an['positive_control_detected'])}. Negative control: {nc['any_perfect']} of {nc['labellings']:,} random "
         "labellings split perfectly. Data: out/features.jsonl, out/analysis.json.", ha="center", fontsize=7.5)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(here / "features.svg")
fig.savefig(here / "features.png", dpi=150)

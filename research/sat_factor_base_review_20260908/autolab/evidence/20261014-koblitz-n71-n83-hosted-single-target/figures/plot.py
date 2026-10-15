"""Plot recorded hosted single-target panel values (SVG + 150 dpi PNG).

Plots only values recorded in the frozen run records; no new computation
beyond medians the harness already reported.
"""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

E = "research/sat_factor_base_review_20260908/autolab/evidence/20261014-koblitz-n71-n83-hosted-single-target"
data = json.load(open(f"{E}/figures/data.json"))

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
for ax, tag, rid in zip(axes, ("n71", "n83"),
                        ("gha-38585913645-1", "gha-38597929540-1")):
    rows = [r for r in data[tag] if "seed" in r]
    seeds = [r["seed"] for r in rows]
    ic = [r["ic_online_ms"] for r in rows]
    rho = [r["rho_online_ms"] for r in rows]
    median = data[tag][-1]["median_speedup"]
    ax.plot(seeds, rho, "o-", label="rho (strong rung 3) online ms", color="#c44")
    ax.plot(seeds, ic, "s-", label="IC compact orbit online ms", color="#146")
    ax.set_yscale("log")
    ax.set_xlabel("one-target eval workload seed")
    ax.set_ylabel("online wall (ms, log scale)")
    ax.set_title(f"{tag}: paired single-target online, hosted EPYC\n"
                 f"run {rid}, K=600, W=16, median rho/IC = {median:.2f}x")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
fig.suptitle("Koblitz compact-orbit single-target panels, hosted isolated reruns (AMD EPYC 7763, 0 contended)", fontsize=9)
fig.tight_layout()
for out in (f"{E}/figures/n71-n83-hosted-online-ms.svg",):
    fig.savefig(out)
fig.savefig(f"{E}/figures/n71-n83-hosted-online-ms.png", dpi=150)
print("written")

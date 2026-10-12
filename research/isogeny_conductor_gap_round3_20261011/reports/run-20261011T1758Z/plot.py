"""Plot the recorded values of ecbench session e5_combined_n17_t1 (measured, non-warm-up records only)."""
import json, os
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt
here = os.path.dirname(os.path.abspath(__file__))
recs = [json.loads(l) for l in open(os.path.join(here, "../../sessions/e5_combined_n17_t1/records.jsonl"))]
arms = ["rho-neg", "f4-blind", "f4-lifts", "f4-combined"]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(8, 3.6))
for k, arm in enumerate(arms):
    rs = [r for r in recs if r["arm"] == arm and not r["warmup"]]
    walls = [max(r["online"]["wall_ns"] / 1e9, 1e-4) for r in rs]
    a1.scatter([k] * len(walls), walls, color="#3b6ea8", s=40, zorder=3)
    ops = [r["solver"]["ops"] for r in rs if r.get("solver")]
    if ops:
        a2.scatter([k] * len(ops), ops, color="#c0632b", s=40, zorder=3)
for ax, lab in ((a1, "online wall per target (s)"), (a2, "F4 solver operations")):
    ax.set_yscale("log"); ax.set_xticks(range(4)); ax.set_xticklabels(arms, rotation=20, fontsize=8)
    ax.set_ylabel(lab); ax.grid(True, which="both", alpha=0.3); ax.set_xlim(-0.5, 3.5)
fig.suptitle("e5_combined_n17_t1 (2 measured rounds per arm; rho has no solver)", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(here, "e5_arms.svg"))
print("ok")

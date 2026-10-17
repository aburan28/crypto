import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

recs = json.load(open("research/p256_trapdoor_round2_20261012/records.json"))["probes"]

# Graph 1: Hamming-weight z-scores (recorded values only).
zs = {"b": 0.81, "r": -0.31, "j": -0.12, "n": 4.88}
fig, ax = plt.subplots(figsize=(7, 4))
names = list(zs); vals = [zs[k] for k in names]
colors = ["#4477aa" if abs(v) < 5 else "#cc3311" for v in vals]
ax.bar(names, vals, color=colors)
ax.axhline(0, color="black", lw=0.8)
for y in (5, -5):
    ax.axhline(y, color="#cc3311", ls="--", lw=1, label="±5σ anomaly threshold" if y == 5 else None)
ax.set_ylabel("Hamming-weight z-score vs binomial(256, 1/2)")
ax.set_title("P-256 round-2 probes: recorded z-scores (run-20261012T2150Z, n=256 bits each)")
ax.legend()
fig.tight_layout()
fig.savefig("research/p256_trapdoor_round2_20261012/reports/run-20261012T2150Z/graph.svg")

# Graph 2: factor-size profiles of p-1 and p+1 (recorded).
fig, ax = plt.subplots(figsize=(7, 4))
labels = ["p-1", "p+1"]
sizes = [[1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 160], [96, 1, 19, 140]]
unfactored = [0, 140]
for i, (lab, s) in enumerate(zip(labels, sizes)):
    ax.bar([f"{lab} f{j}" for j in range(len(s))], s, color="#4477aa")
ax.set_ylabel("factor size (bits)")
ax.set_title("Recorded factorizations of p−1 and p+1 (run-20261012T2150Z); p+1 has a 140-bit unfactored residue")
plt.xticks(rotation=45, ha="right", fontsize=7)
fig.tight_layout()
fig.savefig("research/p256_trapdoor_round2_20261012/reports/run-20261012T2150Z/factors.svg")
print("graphs written")

"""Plot recorded V2 results (figures/zero_rows.csv); no computation beyond layout."""
import csv, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rows = list(csv.DictReader(open(sys.argv[1])))
style = {"seen": ("#8a8f98", "o", "seen during V1 diagnosis"),
         "pristine": ("#2a6fdb", "s", "pristine holdout"),
         "forward": ("#d9480f", "D", "forward prediction, scored"),
         "pending": ("#d9480f", "D", "forward prediction, census running")}
fig, ax = plt.subplots(figsize=(7.2, 4.2))
done = set()
for r in rows:
    p = int(r["p"]); c, m, lab = style[r["role"]]
    lo, mid, hi = float(r["lo95"]), float(r["predicted"]), float(r["hi95"])
    ax.vlines(p, lo, hi, color="#c9ced6", lw=6, zorder=1)
    ax.plot(p, mid, "_", color="#4a4f57", ms=14, mew=2, zorder=2)
    if r["observed_zero_rows"]:
        ax.plot(p, float(r["observed_zero_rows"]), m, color=c, ms=7, zorder=3,
                label=None if r["role"] in done else lab)
    else:
        ax.plot(p, mid, m, mfc="none", color=c, ms=7, zorder=3,
                label=None if r["role"] in done else lab)
    done.add(r["role"])
ax.plot([], [], "_", color="#4a4f57", ms=14, mew=2, label="predicted mean")
ax.plot([], [], color="#c9ced6", lw=6, label="95% interval (Poisson-binomial)")
ax.set_xlabel("base prime p  (curves over $F_{p^6}$)")
ax.set_ylabel("high-depth classes with no weak curve (t > 0)")
ax.set_title("ISO-1 residual zero classes: observed vs registered model V2", fontsize=10)
ax.legend(fontsize=7.5, loc="upper left", frameon=False)
ax.grid(axis="y", color="#eceef1")
fig.text(0.01, 0.005, "Certified censuses (certified_20261016); predictions from class data only. "
         "Source: residual_model_20261016/model_v2_results.txt, predictions_v2.txt", fontsize=6.5, color="#555")
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(sys.argv[2])
print("wrote", sys.argv[2])

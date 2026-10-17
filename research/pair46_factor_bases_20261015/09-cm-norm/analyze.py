"""Layout/derivation only: ratios of recorded coverages (raw.json from the Rust run) + plots. No new measurement."""
import json, random, statistics as st
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
raw = json.load(open("raw.json"))
arms = ["ideal_l1", "ideal_l2", "ap_pos_and_floor", "smallx"]
summ = {}
for m in (1, 2, 3):
    for B in (32, 64, 128, 256, 512, 1024):
        per = {a: [] for a in arms}; nullr = []
        for c in raw["toys"]:
            r = [x for x in c["rows"] if x["B"] == B and x["m"] == m]
            if not r: continue
            r = r[0]; mu = st.mean(r["null"]); ns = sorted(r["null"])
            for a in arms: per[a].append(r[a] / mu)
            nullr.append((ns[0] / mu, ns[-1] / mu))
        if not nullr: continue
        rng = random.Random(7); out = {}
        for a in arms:
            v = per[a]; bs = sorted(st.mean(rng.choice(v) for _ in v) for _ in range(2000))
            out[a] = {"mean_ratio": st.mean(v), "ci95_boot_curves": [bs[50], bs[1949]], "min": min(v), "max": max(v)}
        out["null_ratio_range_over_20_seeds"] = [min(x[0] for x in nullr), max(x[1] for x in nullr)]
        out["n_curves"] = len(nullr); summ[f"m{m}_B{B}"] = out
json.dump({"protocol_sha256": open("PROTOCOL.sha256").read().split()[0], "pair46_exact": raw["pair46_exact"], "toy_curves": [{k: c[k] for k in ("idx","p","t","u","N","a","b","lambda1","lambda2")} for c in raw["toys"]], "coverage_ratio_vs_null": summ}, open("results.json", "w"), indent=1)
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
cols = {"ideal_l1": "#c0392b", "ideal_l2": "#e67e22", "ap_pos_and_floor": "#2c3e50", "smallx": "#27ae60"}
lab = {"ideal_l1": "ideal/norm-form (lambda1)", "ideal_l2": "ideal/norm-form (lambda2)", "ap_pos_and_floor": "AP = pos. control = floor Z[pi]", "smallx": "small-x baseline"}
for ax, m, Bs in ((axs[0], 2, (32, 64, 128, 256, 512, 1024)), (axs[1], 3, (32, 64, 128))):
    lo = [summ[f"m{m}_B{B}"]["null_ratio_range_over_20_seeds"][0] for B in Bs]; hi = [summ[f"m{m}_B{B}"]["null_ratio_range_over_20_seeds"][1] for B in Bs]
    ax.fill_between(Bs, lo, hi, color="#999", alpha=.35, label="random null (min-max, 20 seeds x curves)")
    for a in arms:
        y = [summ[f"m{m}_B{B}"][a]["mean_ratio"] for B in Bs]
        e = [[y[i] - summ[f"m{m}_B{B}"][a]["ci95_boot_curves"][0] for i, B in enumerate(Bs)], [summ[f"m{m}_B{B}"][a]["ci95_boot_curves"][1] - y[i] for i, B in enumerate(Bs)]]
        ax.errorbar(Bs, y, yerr=e, marker="o", ms=4, capsize=2, color=cols[a], label=lab[a], ls="--" if a == "ideal_l2" else "-")
    ax.set_xscale("log", base=2); ax.set_yscale("log"); ax.set_xlabel("factor-base size B (pairs +-)"); ax.set_ylabel(f"coverage_{m} / mean random coverage_{m}")
    ax.set_title(f"m = {m}: sums of <= {m} factor-base points"); ax.grid(alpha=.3)
axs[0].legend(fontsize=7, loc="lower left")
fig.suptitle("09-cm-norm: relation coverage on 8 toy j=-96^3 prime-order curves (N~2^20, 2^22); 95% bootstrap CI over curves", fontsize=9)
fig.tight_layout(); fig.savefig("fig_coverage.svg"); fig.savefig("fig_coverage.png", dpi=150)
for k in ("m2_B256", "m2_B1024", "m3_B128"): print(k, json.dumps(summ[k]))

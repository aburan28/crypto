"""Merge the per-cell raw JSON (perm and random-base arms), apply the frozen
decision rule to the recorded values, write summary.csv / summary.json and the
figures.  Usage: python3 -I .../analysis/report.py
"""
import csv
import glob
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ITEM = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RAW = os.path.join(ITEM, "raw")
FIG = os.path.join(ITEM, "figures")
STATS = ["maxZ", "Ztmax", "Z_C4max", "Zt_C4max", "M1", "M10", "M100", "ZtM1", "ZtM10", "ZtM100", "marginal"]
MASS = ["M1", "M10", "M100", "ZtM1", "ZtM10", "ZtM100"]


def q(v, a):
    return float(np.quantile(np.asarray(v, float), a))


def load_cells():
    cells = {}
    for f in sorted(glob.glob(os.path.join(RAW, "cell_*_perm.json"))):
        perm = json.load(open(f))
        base = perm["name"][: -len("_perm")]
        rf = os.path.join(RAW, base + "_rand.json")
        rand = json.load(open(rf)) if os.path.exists(rf) else None
        for extra in sorted(glob.glob(os.path.join(RAW, base + "_rand[0-9].json"))):
            e = json.load(open(extra))
            if rand is None:
                rand = e
                continue
            assert abs(e["observed"]["Ztmax"] - rand["observed"]["Ztmax"]) < 1e-9
            for k in STATS:
                rand["rand_values"][k] += e["rand_values"][k]
            rand["seconds"] += e["seconds"]
        for extra in sorted(glob.glob(os.path.join(RAW, base + "_perm[0-9].json"))):
            e = json.load(open(extra))
            assert abs(e["observed"]["Ztmax"] - perm["observed"]["Ztmax"]) < 1e-9
            for k in STATS:
                perm["perm_values"][k] += e["perm_values"][k]
            perm["seconds"] += e["seconds"]
        cells[base] = (perm, rand)
    return cells


def decide(base, perm, rand):
    o = perm["observed"]
    pv = perm["perm_values"]
    rv = rand["rand_values"] if rand else {k: [] for k in STATS}
    assert rand is None or abs(rand["observed"]["Ztmax"] - o["Ztmax"]) < 1e-9
    row = dict(cell=base, p=perm["p"], a=perm["a"], b=perm["b"], t=perm["t"], m=perm["m"], N=perm["N"],
               P=o["P"], targets_ge2=perm["info"]["targets_ge2"], nperm=len(pv["maxZ"]), nrand=len(rv["maxZ"]))
    for k in STATS:
        row[k] = o[k]
        row[k + "_perm_q999"] = q(pv[k], 0.999)
        row[k + "_perm_q001"] = q(pv[k], 0.001)
        row[k + "_perm_med"] = q(pv[k], 0.5)
        row[k + "_perm_p"] = (1 + sum(v >= o[k] for v in pv[k])) / (1 + len(pv[k]))
        if rv[k]:
            row[k + "_rand_q999"] = q(rv[k], 0.999)
            row[k + "_rand_med"] = q(rv[k], 0.5)
            row[k + "_rand_p"] = (1 + sum(v >= o[k] for v in rv[k])) / (1 + len(rv[k]))
    row["exceed_Z"] = o["maxZ"] > row["maxZ_perm_q999"]
    row["exceed_Zt"] = o["Ztmax"] > row["Ztmax_perm_q999"]
    row["exceed_C4"] = o["Z_C4max"] > row["Z_C4max_perm_q999"] or o["Zt_C4max"] > row["Zt_C4max_perm_q999"]
    row["exceed_marginal"] = o["marginal"] > row["marginal_perm_q999"]
    row["mass_out"] = [k for k in MASS if not (row[k + "_perm_q001"] <= o[k] <= row[k + "_perm_q999"])]
    row["Zargmax"] = o["Zargmax"]
    row["Ztargmax"] = o["Ztargmax"]
    if "planted" in perm:
        row["planted"] = [dict(rho=r["rho"], maxZ=r["maxZ"], Ztmax=r["Ztmax"],
                               det_Z=r["maxZ"] > row["maxZ_perm_q999"], det_Zt=r["Ztmax"] > row["Ztmax_perm_q999"],
                               argmax_in_subgroup=r["argmax_in_cstar_subgroup"],
                               Zt_cstar=r["Z_cstar"][1]) for r in perm["planted"]["rows"]]
    row["sec"] = perm["seconds"] + (rand["seconds"] if rand else 0)
    return row


def label(row):
    if row["cell"].startswith("cell_c4"):
        return "C4 curve\nm=3 t=8"
    bits = {1073153: 20, 4206593: 22, 16801793: 24}[row["p"]]
    c = row["cell"].split("_c")[-1]
    return f"{bits}b m={row['m']} c{c}"


def fig_cells(rows, cells):
    fig, axes = plt.subplots(2, 1, figsize=(13, 8.5), sharex=True)
    xs = np.arange(len(rows))
    for ax, key, title in [(axes[0], "maxZ", "frozen normalisation Z (max |Z| over all character classes)"),
                           (axes[1], "Ztmax", "type-standardised Zt (max |Zt| over all character classes)")]:
        for i, r in enumerate(rows):
            perm, rand = cells[r["cell"]]
            pv = perm["perm_values"][key]
            ax.scatter(np.full(len(pv), i - 0.18) + np.random.default_rng(i).uniform(-0.07, 0.07, len(pv)), pv,
                       s=3, color="#9aa5b1", alpha=0.6, lw=0, label="permutation null (200)" if i == 0 else None)
            if rand:
                rv = rand["rand_values"][key]
                ax.scatter(np.full(len(rv), i + 0.18) + np.random.default_rng(100 + i).uniform(-0.07, 0.07, len(rv)), rv,
                           s=3, color="#e0a458", alpha=0.7, lw=0, label="random-factor-base null" if i == 0 else None)
            ax.plot([i - 0.32, i + 0.32], [r[key + "_perm_q999"]] * 2, color="#c0392b", lw=1.2,
                    label="perm-null 99.9th pct (kill threshold)" if i == 0 else None)
            ax.scatter([i], [r[key]], s=36, marker="D", color="#1f4e79", zorder=5,
                       label="Kummer tower (observed)" if i == 0 else None)
        ax.set_ylabel(key)
        ax.set_title(title, fontsize=10, loc="left")
        ax.set_yscale("log")
        ax.grid(axis="y", alpha=0.3)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, fontsize=8, loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.0), markerscale=1.5)
    axes[1].set_xticks(xs)
    axes[1].set_xticklabels([label(r) for r in rows], rotation=70, fontsize=7)
    fig.suptitle("Item 4: max Fourier coefficient of the target-conditioned decomposition spectrum, tower vs nulls", fontsize=11)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    for ext, kw in [("svg", {}), ("png", {"dpi": 150})]:
        fig.savefig(os.path.join(FIG, f"maxcoef_by_cell.{ext}"), **kw)
    plt.close(fig)


def fig_planted(rows):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, key in [(axes[0], "maxZ"), (axes[1], "Ztmax")]:
        for r in rows:
            if "planted" not in r:
                continue
            rho = [x["rho"] for x in r["planted"]]
            val = [x[key] for x in r["planted"]]
            l, = ax.plot(rho, [v / r[key + "_perm_q999"] for v in val], marker="o", label=label(r).replace("\n", " ").replace(" c0", ""))
        ax.axhline(1.0, color="#c0392b", lw=1, ls="--", label="kill threshold (perm 99.9th)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("planted fraction ρ of orbits on c*·e ≡ s(R)")
        ax.set_ylabel(f"planted {key} / perm-null 99.9th pct")
        ax.grid(alpha=0.3)
    axes[1].legend(fontsize=7)
    fig.suptitle("Positive control: planted bias in one weighted character c* = (1,3) / (1,3,5)", fontsize=11)
    fig.tight_layout()
    for ext, kw in [("svg", {}), ("png", {"dpi": 150})]:
        fig.savefig(os.path.join(FIG, f"planted_control.{ext}"), **kw)
    plt.close(fig)


def main():
    os.makedirs(FIG, exist_ok=True)
    cells = load_cells()
    rows = [decide(b, *cells[b]) for b in cells]
    order = lambda r: (r["cell"].startswith("cell_c4"), r["m"], r["p"], r["cell"])
    rows.sort(key=order)
    json.dump(rows, open(os.path.join(ITEM, "summary.json"), "w"), indent=1, default=float)
    cols = ["cell", "p", "a", "b", "t", "m", "P", "targets_ge2", "nperm", "nrand",
            "maxZ", "maxZ_perm_q999", "maxZ_perm_p", "maxZ_rand_p", "Ztmax", "Ztmax_perm_q999", "Ztmax_perm_p", "Ztmax_rand_p",
            "Z_C4max", "Z_C4max_perm_q999", "Zt_C4max", "Zt_C4max_perm_q999", "M100", "M100_perm_q001", "M100_perm_q999",
            "ZtM100", "ZtM100_perm_q001", "ZtM100_perm_q999", "marginal", "marginal_perm_q999",
            "exceed_Z", "exceed_Zt", "exceed_C4", "exceed_marginal", "mass_out"]
    with open(os.path.join(ITEM, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            w.writerow([r.get(c) if not isinstance(r.get(c), float) else f"{r[c]:.5g}" for c in cols])
    fig_cells(rows, cells)
    fig_planted(rows)
    n_ex = sum(r["exceed_Z"] or r["exceed_Zt"] for r in rows)
    print(f"cells={len(rows)} exceed(Z or Zt)={n_ex} C4={sum(r['exceed_C4'] for r in rows)} "
          f"marginal={sum(r['exceed_marginal'] for r in rows)} mass_out={sum(bool(r['mass_out']) for r in rows)}")
    for r in rows:
        print(f"{r['cell']:<22} P={r['P']:>8} Z={r['maxZ']:7.2f}/{r['maxZ_perm_q999']:7.2f} p={r['maxZ_perm_p']:.3f} "
              f"Zt={r['Ztmax']:5.2f}/{r['Ztmax_perm_q999']:5.2f} p={r['Ztmax_perm_p']:.3f} "
              f"rand_p={r.get('Ztmax_rand_p', float('nan')):.3f} C4={r['exceed_C4']} mass={r['mass_out']} "
              + (" planted=" + " ".join(f"{x['rho']}:{int(x['det_Z'])}{int(x['det_Zt'])}" for x in r["planted"]) if "planted" in r else ""))


if __name__ == "__main__":
    main()

"""Lay out the recorded values of benchmarks/best-of-two as SVG/PNG graphs.
Reads the JSON lines and TSVs the runs wrote; computes nothing beyond placing
them (and the c +- 2 se bars the records carry).

    python plot.py
"""
import csv
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
plt.rcParams.update({"font.size": 8.5, "svg.fonttype": "none"})

rows = []
for f in ("emulation-v1.jsonl", "emulation-v2.jsonl", "emulation-v3.jsonl"):
    rows += [json.loads(l) for l in (HERE / f).read_text().splitlines() if l.strip()]


def label(d):
    w, sc, sg = d["walk"], d.get("score"), d.get("sign")
    if sg is None:  # emulation-v2 predates --sign: bit 4 under native, a hash otherwise
        sg = "hw4" if d["dist"] == "native" else "hash"
    if w == "sigma":
        return "shipping σ walk"
    if w == "sigmapm":
        return {"sj": "± by parity of x∧σʲx (SIGMA_SIGN=2)", "hw4": "± by weight bit 4 (SIGMA_SIGN=1)",
                "hash": "± by class hash (fair sign)"}[sg]
    if w == "sigma2":
        return f"best of R±σʲR, score {sc}"
    if w == "sigma2t":
        return f"threshold on R+σʲR, score {sc}"
    return w


STYLE = {
    "shipping σ walk": ("#333333", "o", "-"),
    "± by parity of x∧σʲx (SIGMA_SIGN=2)": ("#c0392b", "s", "-"),
    "± by weight bit 4 (SIGMA_SIGN=1)": ("#e59866", "s", ":"),
    "± by class hash (fair sign)": ("#e59866", "s", "--"),
    "best of R±σʲR, score hw": ("#2f6db5", "^", "-"),
    "best of R±σʲR, score key": ("#2f6db5", "^", "--"),
    "threshold on R+σʲR, score hw": ("#27ae60", "v", "-"),
    "threshold on R+σʲR, score key": ("#27ae60", "v", "--"),
}

fig, axes = plt.subplots(1, 2, figsize=(7.4, 4.4), sharey=True)
for ax, dist, title in zip(axes, ("ecc2k130", "native"),
                           ("ECC2K-130 branch probabilities\n(hashed; averaged over mappings)",
                            "device branch choice\n(one fixed mapping per degree)")):
    series = {}
    for d in rows:
        if d["dist"] != dist:
            continue
        lab = label(d)
        if d["walk"] == "sigma" and lab in series and any(x[0] == d["n"] for x in series[lab]):
            continue  # the v3 repeat of the shipping walk; v1's row is plotted
        series.setdefault(lab, []).append((d["n"], d["c"], d["c_se"]))
    for lab, pts in sorted(series.items(), key=lambda kv: list(STYLE).index(kv[0]) if kv[0] in STYLE else 99):
        pts.sort()
        col, mk, ls = STYLE.get(lab, ("#999999", ".", "-"))
        ax.errorbar([p[0] for p in pts], [p[1] for p in pts], yerr=[2 * p[2] for p in pts],
                    color=col, marker=mk, ls=ls, ms=3.5, lw=1, capsize=1.5, label=lab)
    ax.axhline(1.0, color="#999999", lw=0.7)
    ax.set_title(title, fontsize=8, loc="left")
    ax.set_xlabel("degree n")
axes[0].set_ylabel("c = iterations / random mapping (±2 se)")
merged = {}
for ax in axes:
    for h, l in zip(*ax.get_legend_handles_labels()):
        merged.setdefault(l, h)
order = [l for l in STYLE if l in merged]
fig.legend([merged[l] for l in order], order, fontsize=6.5, frameon=False, loc="lower center", ncol=3)
fig.tight_layout(rect=(0, 0.13, 1, 1))
for ext in ("svg", "png"):
    fig.savefig(HERE / f"walk-constant.{ext}", dpi=160)

# GPU timing panels
fig, axes = plt.subplots(1, 2, figsize=(7.4, 2.6), sharey=True)
for ax, run, name in zip(axes, ("gpu-run1", "gpu-run2"), ("SIGMA_SIGN=1 (weight bit 4)", "SIGMA_SIGN=2 (parity of x∧σʲx)")):
    with open(HERE / run / "results" / "samples.tsv") as f:
        rs = list(csv.DictReader(f, delimiter="\t"))
    for phase, off in (("aa", 0), ("ab", 6)):
        sel = [r for r in rs if r["phase"] == phase]
        for r in sel:
            x = off + int(r["pair"]) + (0.15 if r["variant"] in ("candidate", "control_b") else -0.15)
            col = "#c0392b" if r["variant"] == "candidate" else "#333333"
            ax.plot(x, float(r["rateMps"]) / 1000, "o", color=col, ms=3.5)
    ax.set_xticks([3, 9], ["A/A pairs", "A/B pairs"])
    ax.set_title(name, fontsize=8, loc="left")
axes[0].set_ylabel("billion iterations / s")
fig.text(0.5, 0.01, "RTX PRO 6000, campaign knobs, 385,024 workers, 64 launches × 1024 steps per sample; red = candidate",
         ha="center", fontsize=7)
fig.tight_layout(rect=(0, 0.05, 1, 1))
for ext in ("svg", "png"):
    fig.savefig(HERE / f"gpu-rate.{ext}", dpi=160)
print("wrote walk-constant and gpu-rate graphs")

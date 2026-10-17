"""Figures for revision 2 of the Eigen-Seta paper, drawn from committed
records only (plotting recorded values; no result is computed here).

Usage: python figures/plot_assessment.py <assessment_dir> <out_dir>
  <assessment_dir> = crypto:research/eigen_seta_assessment_20261017
Writes pk_test.{pdf,svg,png} and attack_margins.{pdf,svg,png} into <out_dir>.
"""
import os, re, sys, statistics
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

A, OUT = sys.argv[1], sys.argv[2]
REC = os.path.join(A, "records")
RID = open(os.path.join(REC, "RUN_ID")).read().strip()
SETS = ["L1", "L3", "L5"]


def kv(line):
    return dict(f.split("=", 1) for f in line.rstrip("\n").split("\t") if "=" in f)


def save(fig, stem):
    for ext in ("pdf", "svg"):
        fig.savefig(os.path.join(OUT, f"{stem}.{ext}"))
    fig.savefig(os.path.join(OUT, f"{stem}.png"), dpi=150)
    plt.close(fig)


# --- public-key degree test vs decapsulation -------------------------------
rows = {}
for s in SETS:
    t_pk, t_dec, e_wrong = [], [], []
    e_true = big_e = None
    for line in open(os.path.join(REC, RID, f"{s}.tsv")):
        d = kv(line)
        if line.startswith("keygen"):
            big_e = int(d["E"])
        elif line.startswith("time"):
            t_pk.append(float(d["pk_test_s"]))
            t_dec.append(float(d["decaps_s"]))
            e_true = int(d["pk_test_e"])
        elif line.startswith("wrong") and "e" in d:
            e_wrong.append(int(d["e"]))
    rows[s] = dict(E=big_e, e=e_true, t_pk=t_pk, t_dec=t_dec, e_wrong=e_wrong)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.2))
x = range(len(SETS))
ax1.bar([i - 0.2 for i in x], [rows[s]["E"] for s in SETS], 0.4, label="decapsulation chain E", color="#9aa5b1")
ax1.bar([i + 0.2 for i in x], [rows[s]["e"] for s in SETS], 0.4, label="public-key test e (true f)", color="#2f6db5")
for i, s in enumerate(SETS):
    ys = rows[s]["e_wrong"]
    ax1.scatter([i + 0.2] * len(ys), ys, s=6, color="k", zorder=3, label="e for wrong candidates" if i == 0 else None)
ax1.set_xticks(list(x)); ax1.set_xticklabels(SETS)
ax1.set_ylabel("(2,2)-steps"); ax1.legend(fontsize=7); ax1.grid(True, axis="y", alpha=0.3)
for i, s in enumerate(SETS):
    n = len(rows[s]["t_pk"])
    ax2.scatter([i - 0.12] * n, rows[s]["t_dec"], marker="^", color="#9aa5b1", label="decapsulation" if i == 0 else None, zorder=3)
    ax2.scatter([i + 0.12] * n, rows[s]["t_pk"], marker="o", color="#2f6db5", label="public-key test" if i == 0 else None, zorder=3)
ax2.set_xticks(list(x)); ax2.set_xticklabels([f"{s}\n(n={len(rows[s]['t_pk'])})" for s in SETS])
ax2.set_yscale("log"); ax2.set_ylabel("wall time per test [s], log"); ax2.legend(fontsize=7); ax2.grid(True, which="both", axis="y", alpha=0.3)
fig.suptitle(f"run {RID}", fontsize=8)
fig.tight_layout()
save(fig, "pk_test")

# --- attack-cost margins from the cost model ------------------------------
txt = open(os.path.join(REC, "security_estimate.txt")).read()
blocks = re.findall(r"^(L1|L3|L5|Q3|Q5) \(lambda (\d+).*?floor: classical 2\^([0-9.]+), quantum queries 2\^([0-9.]+), quantum gates 2\^([0-9.]+)/2\^([0-9.]+)/2\^([0-9.]+)\s+reference AES-\d+ Grover gates: 2\^(\d+)/2\^(\d+)/2\^(\d+)", txt, re.M | re.S)
names = [b[0] for b in blocks]
lam = [int(b[1]) for b in blocks]
cl = [float(b[2]) - l for b, l in zip(blocks, lam)]
qq = [float(b[3]) - l for b, l in zip(blocks, lam)]
gm = [min(float(b[4 + k]) - float(b[7 + k]) for k in range(3)) for b in blocks]
fig, ax = plt.subplots(figsize=(6.4, 3.4))
xs = range(len(names))
w = 0.26
b1 = ax.bar([i - w for i in xs], cl, w, label="classical floor − λ", color="#2f6db5")
b2 = ax.bar([i for i in xs], qq, w, label="quantum-query floor − λ (Tier B)", color="#d9822b")
b3 = ax.bar([i + w for i in xs], gm, w, label="min over MAXDEPTH of gates − AES ref. (Tier A)", color="#3a9a5b")
for b in (b1, b2, b3):
    ax.bar_label(b, fmt="%.1f", fontsize=5.5, padding=1)
ax.axhline(0, color="k", lw=0.8)
ax.set_xticks(list(xs)); ax.set_xticklabels([f"{n}\nλ={l}" for n, l in zip(names, lam)])
ax.set_ylabel("margin [bits]; < 0 fails the bar")
ax.legend(fontsize=7); ax.grid(True, axis="y", alpha=0.3)
fig.suptitle("Derived by eigenx examples/security_estimate (cost model, attacker-favourable)", fontsize=8)
fig.tight_layout()
save(fig, "attack_margins")
print({s: (rows[s]["e"], rows[s]["E"], statistics.median(rows[s]["t_pk"]), statistics.median(rows[s]["t_dec"])) for s in SETS})

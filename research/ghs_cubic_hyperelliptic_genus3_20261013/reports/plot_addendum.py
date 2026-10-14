"""Plot recorded values of the addendum run (no computation beyond reading the records).
Usage (repo root): python3 research/ghs_cubic_hyperelliptic_genus3_20261013/reports/plot_addendum.py <run id>"""
import json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = "research/ghs_cubic_hyperelliptic_genus3_20261013/"
RUN = sys.argv[1]
inst = [json.loads(l) for l in open(D + "results/" + RUN + ".jsonl") if '"rec":"instance"' in l]
cen = []
for q in (11, 13):
    p = D + "results/" + RUN + f".census_q{q}.jsonl"
    cen += [json.loads(l) for l in open(p) if '"rec":"curve"' in l]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"width_ratios": [1.5, 1]})
lab, tot, pf, ker, col = [], [], [], [], []
for r in inst:
    d = r["P4_detail"] if r["arm"] == "split" else r["P4c_detail"]
    pfv = d["ext_pushforward_all_E0K_eq_4P"] if r["arm"] == "split" else d["pushforward_all_E0K_eq_4P"]
    lab.append(f"q={r['q']} {'S' if r['arm'] == 'split' else 'C'}({r['Q'][0]},{r['Q'][1]})")
    tot.append(d["E0_points"]); pf.append(pfv); ker.append(d["kernel_size"])
    col.append("#1f6fb4" if r["arm"] == "split" else "#d9822b")
x = range(len(inst)); w = 0.4
ax1.bar([i - w / 2 for i in x], tot, w, color="#bbb", label="|E0(K)|")
ax1.bar([i + w / 2 for i in x], pf, w, color=col, label="π0_*φ(P) = 4P verified, split arm (imaginary model)")
for i in x:
    ax1.text(i, tot[i] * 1.1, f"ker {ker[i]}", ha="center", fontsize=6, rotation=90)
ax1.set_yscale("log"); ax1.set_ylim(10, 3e4); ax1.set_xticks(list(x)); ax1.set_xticklabels(lab, fontsize=6.5, rotation=60, ha="right")
ax1.set_ylabel("points of E0(K)"); ax1.set_title("Transfer φ on all 20 instances (hom. 400/400, 0 failures each)", fontsize=9.5)
ax1.bar([0], [0], 0, color="#d9822b", label="π0_*φ(P) = 4P verified, conjugate arm (real model, addendum)"); ax1.legend(fontsize=6.5, loc="upper left")

names, g3hyp, g3quart = [], [], []
for c in cen:
    t = {(e["m"], e["r"], e["hyperelliptic"]): e["structures"] for e in c["types_at_min"]}
    names.append(f"q={c['q']} ({c['A_idx']},{c['B_idx']})")
    g3hyp.append(t.get((3, 5, True), 0)); g3quart.append(t.get((2, 6, None), 0))
y = range(len(cen))
ax2.barh([i + 0.2 for i in y], g3hyp, 0.4, color="#2a9d5c", label="genus 3, m=3, r=5 (hyperelliptic)")
ax2.barh([i - 0.2 for i in y], g3quart, 0.4, color="#7b5ea7", label="genus 3, m=2, r=6")
ax2.set_yticks(list(y)); ax2.set_yticklabels(names, fontsize=6.5); ax2.invert_yaxis()
ax2.set_xlabel("F_q-structures attaining genus 3\n(of 1,786,323 at q=11; 4,855,539 at q=13)", fontsize=8)
ax2.set_title("P5c: census over every structure, q = 11, 13", fontsize=9.5); ax2.legend(fontsize=6.5, loc="lower right")
fig.suptitle(f"Hyperelliptic genus-3 GHS cover, addendum run {RUN}", fontsize=10)
fig.tight_layout()
out = D + "reports/" + RUN + "/addendum"
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out + ".svg"); fig.savefig(out + ".png", dpi=150)

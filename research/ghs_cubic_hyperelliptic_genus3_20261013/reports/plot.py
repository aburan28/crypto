"""Plot recorded values of run 20261014T0029Z (no computation beyond reading the records).
Usage (repo root): python3 research/ghs_cubic_hyperelliptic_genus3_20261013/reports/plot.py"""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = "research/ghs_cubic_hyperelliptic_genus3_20261013/"
RUN = "run-20261014T0029Z"
inst = [json.loads(l) for l in open(D + "results/" + RUN + ".jsonl") if '"rec":"instance"' in l]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))
for arm, mk, col in [("split", "o", "#1f6fb4"), ("conjugate", "s", "#d9822b")]:
    xs = [r["order_E0"] for r in inst if r["arm"] == arm]
    ys = [r["L_of_1"] for r in inst if r["arm"] == arm]
    ax1.scatter(xs, ys, marker=mk, s=46, facecolors="none", edgecolors=col, label=f"{arm} arm: L(1) from N1,N2,N3 on w²=h(z)")
sp = [r for r in inst if r["arm"] == "split"]
ax1.scatter([r["order_E0"] for r in sp], [r["Jac_Fq_by_Mumford_enumeration"] for r in sp], marker="x", s=30, color="#1f6fb4", label="split arm: Mumford-pair enumeration")
lim = [80, 2400]
ax1.plot(lim, lim, color="#888", lw=0.8, ls="--", label="#Jac(C)(F_q) = #E0(K)")
ax1.set_xscale("log"); ax1.set_yscale("log"); ax1.set_xlim(lim); ax1.set_ylim(lim)
ax1.set_xlabel("#E0(K), K = F_{q³} (counted on y0² = F0(x))"); ax1.set_ylabel("#Jac(C)(F_q)")
ax1.set_title("P3: Jacobian order of the F_q model (20 instances)", fontsize=10); ax1.legend(fontsize=7, loc="upper left")
lab = [f"q={r['q']} ({r['Q'][0]},{r['Q'][1]})" for r in sp]
x = range(len(sp)); w = 0.27
ax2.bar([i - w for i in x], [r["P4_detail"]["E0_points"] for r in sp], w, color="#bbb", label="|E0(K)|")
ax2.bar([i for i in x], [r["P4_detail"]["pushforward_eq_4P"] for r in sp], w, color="#1f6fb4", label="π0_*φ(P)=4P, preregistered set (no special-fibre points)")
ax2.bar([i + w for i in x], [r["P4_detail"]["ext_pushforward_all_E0K_eq_4P"] for r in sp], w, color="#2a9d5c", label="π0_*φ(P)=4P, all of E0(K) (amendment)")
for i, r in enumerate(sp):
    ax2.text(i, r["P4_detail"]["E0_points"] * 1.08, f"ker {r['P4_detail']['kernel_size']}", ha="center", fontsize=6.5)
ax2.set_yscale("log"); ax2.set_xticks(list(x)); ax2.set_xticklabels(lab, fontsize=7, rotation=50, ha="right"); ax2.set_ylim(10, 2e4); ax2.set_xlabel("instance: q, (a, b)")
ax2.set_ylabel("points of E0(K)"); ax2.set_title("P4: transfer φ, split arm (0 failures; hom. 400/400 each)", fontsize=10)
ax2.legend(fontsize=6.5, loc="upper left")
fig.suptitle(f"Hyperelliptic genus-3 GHS cover (m=3, r=5), run {RUN}", fontsize=10)
fig.tight_layout()
out = D + "reports/" + RUN + "/checks"
import os; os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out + ".svg"); fig.savefig(out + ".png", dpi=150)

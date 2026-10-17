"""Structure diagram for the Eigen-Seta assessment (no data; layout only).
Usage: python diagram.py <out_dir>  -> attack_map.{svg,pdf,png}"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

OUT = sys.argv[1]
fig, ax = plt.subplots(figsize=(10.0, 5.6))
ax.set_xlim(0, 12); ax.set_ylim(0, 6.4); ax.axis("off")

def node(x, y, title, sub, col):
    ax.add_patch(FancyBboxPatch((x - 1.4, y - 0.55), 2.8, 1.1, boxstyle="round,pad=0.05", fc=col, ec="k", lw=0.8))
    ax.text(x, y + 0.25, title, ha="center", va="center", fontsize=9, weight="bold")
    ax.text(x, y - 0.18, sub, ha="center", va="center", fontsize=6.5)

def arrow(x0, y0, x1, y1, label, style="-|>", ls="-", col="k", dy=0.22):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle=style, lw=1.0, ls=ls, color=col))
    ax.text((x0 + x1) / 2, (y0 + y1) / 2 + dy, label, ha="center", va="center", fontsize=6.5, color=col)

node(1.8, 5.2, "E_D  (j = -3375)", "public; CM by Z[(1+sqrt-7)/2]\nomega-eigenbasis on E_D[2^E]", "#e8eef7")
node(6.0, 5.2, "E_A  (public key)", "lambda1*eps(P_D), lambda2*eps(Q_D)\ntau = eps omega eps^ hidden (f unknown)", "#e8eef7")
node(10.2, 5.2, "E_C  (ciphertext)", "a*phi(P_A), b*phi(Q_A)\nab public, b/a hidden", "#e8eef7")
arrow(3.25, 5.2, 4.55, 5.2, "eps\ndegree f (secret)", dy=0.3)
arrow(7.45, 5.2, 8.75, 5.2, "phi\ndegree d = 3^F", dy=0.3)

node(6.0, 2.6, "E_A x E_A", "(2^e,2^e) Kani chain\ne ~ 2 log2 f + 3 (new test)", "#dff0e3")
arrow(6.0, 4.6, 6.0, 3.2, "", col="#1d6b35")
ax.text(6.1, 3.9, "candidate f:\ngamma = u psi_A + c", fontsize=6.5, color="#1d6b35", ha="left", va="center")
ax.text(6.0, 1.6, "true f: splits to E_A x E_A\nwrong f: full-length walk, no split", fontsize=6.5, ha="center", va="center", color="#1d6b35")

node(10.2, 2.6, "E_C x E_C", "decapsulation\n(2^E,2^E) chain, E ~ 4 log2 f", "#f2f2f2")
arrow(10.2, 4.6, 10.2, 3.2, "", )
ax.text(10.3, 3.9, "secret key\n(c, c', lambda_S, lambda_T)", fontsize=6.5, ha="left", va="center")
ax.text(10.2, 1.3, "binding attack: guess b/a mod 2^n, 2^(2n) > d\n2^lambda classical; 2^(lambda/2) Grover on L1, L3\nmask-free embedding: open (risk 2)", ha="center", fontsize=6.5)

node(1.8, 2.6, "Z[f omega]-orientation", "hidden; public: 2E steps\nof the p-cycle through E_A", "#fbe9e0")
arrow(4.7, 4.65, 2.9, 3.2, "", col="#b5482f", ls="--")
ax.text(3.0, 4.15, "eigenlines are horizontal\n2-walks (risk 1)", fontsize=6.5, color="#b5482f", ha="center")
ax.text(1.8, 1.4, "if made effective: Kuperberg,\nsubexponential quantum key\nrecovery at every set", ha="center", fontsize=6.5, color="#b5482f")
ax.text(6.0, 0.35, "Solid: computed map. Dashed: structural relation, not an algorithm. Status: diagram of the analysis, not a measurement.", ha="center", fontsize=6)
fig.tight_layout()
for ext in ("svg", "pdf"):
    fig.savefig(os.path.join(OUT, f"attack_map.{ext}"))
fig.savefig(os.path.join(OUT, "attack_map.png"), dpi=150)

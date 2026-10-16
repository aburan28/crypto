# Layout only: draws the NFS -> ECDLP dictionary of README section 3 as a
# diagram. Every cell's text and status comes from the report; nothing is
# computed here.
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

here = os.path.dirname(os.path.abspath(__file__))

STATUS = {
    "exists": ("#1e8449", "exists (literature)"),
    "ext": ("#b9770e", "extension fields only"),
    "open": ("#7d3c98", "open / unproven"),
    "none": ("#922b21", "no analogue known"),
}

rows = [
    ("1  Global ring with unique\n   factorisation: Z[a], ideals",
     [("E(K) over a number field: Mordell-Weil rank tiny, heights quadratic\n"
       "(Miller 1986; xedni fails, Jacobson et al. 2000)", "none"),
      ("Weil restriction + summation polynomials (Semaev, Gaudry, Diem)", "ext"),
      ("Pairing to F_{p^k}*: the field itself (MOV / Frey-Rueck)", "exists")]),
    ("2  Size function with a\n   smoothness probability (norm)",
     [("decomposition probability ~ |V|^m / (m! r): combinatorial, no Dickman rho,\n"
       "no L(1/3)", "open")]),
    ("3  'Special' knob that shrinks\n   the size: p = f(m), small f",
     [("subfield F_q in F_{q^n}: Diem, subexponential for suitable (q, n)", "ext"),
      ("smooth subgroup of F_p* or smooth auxiliary curve (Petit-Kosters-Messeng 2016)", "open"),
      ("p = p(u) for pairing curves (BN, BLS) -> SexTNFS on F_{p^12} (Kim-Barbulescu 2016)", "exists"),
      ("hidden isogeny to a GHS-weak curve (Teske trapdoor, F_{2^161})", "ext")]),
    ("4  The sieve: divisibility on a\n   public lattice finds smooth (a, b)",
     [("no divisibility on E(F_p): each decomposition test is a polynomial-system\n"
       "solve (Groebner / root finding), the dominant cost", "none")]),
    ("5  Sparse linear algebra mod r", [("same: sparse LA mod r over factor-base columns", "exists")]),
    ("6  Descent for one target", [("target point enters relation search directly", "exists")]),
]

fig_h = 10.4
fig, ax = plt.subplots(figsize=(10.5, fig_h))
ax.set_xlim(0, 10.5)
ax.set_ylim(0, fig_h)
ax.axis("off")
ax.text(0.15, fig_h - 0.35, "NFS for F_p* (what SNFS relies on)", fontsize=10, weight="bold")
ax.text(3.75, fig_h - 0.35, "Elliptic-curve analogue (ECDLP)", fontsize=10, weight="bold")

y = fig_h - 0.75
for left, cells in rows:
    h = 0.62 * len(cells) + 0.1
    ax.add_patch(FancyBboxPatch((0.1, y - h), 3.3, h, boxstyle="round,pad=0.02",
                                fc="#eef2f7", ec="#5d6d7e", lw=0.8))
    ax.text(0.22, y - h / 2, left, va="center", fontsize=8.3)
    cy = y - 0.08
    for text, st in cells:
        col, _ = STATUS[st]
        ax.add_patch(FancyBboxPatch((3.7, cy - 0.54), 6.65, 0.52, boxstyle="round,pad=0.02",
                                    fc="white", ec=col, lw=1.4))
        ax.add_patch(FancyBboxPatch((3.7, cy - 0.54), 0.12, 0.52, boxstyle="square,pad=0",
                                    fc=col, ec=col))
        ax.text(3.95, cy - 0.28, text, va="center", fontsize=7.4)
        ax.annotate("", xy=(3.68, cy - 0.28), xytext=(3.42, y - h / 2),
                    arrowprops=dict(arrowstyle="->", color="#5d6d7e", lw=0.7))
        cy -= 0.62
    y -= h + 0.18

lx = 0.15
for st, (col, lab) in STATUS.items():
    ax.add_patch(FancyBboxPatch((lx, 0.18), 0.25, 0.2, boxstyle="square,pad=0", fc=col, ec=col))
    ax.text(lx + 0.33, 0.28, lab, va="center", fontsize=8)
    lx += 2.55
ax.text(0.15, 0.62, "Status is of the analogue as a route to ECDLP on a prime-order subgroup; "
        "proposed or derived links are not measured speedups. See README sections 3-4.",
        fontsize=7.2, color="#444")
fig.savefig(os.path.join(here, "dictionary.svg"), bbox_inches="tight")
fig.savefig(os.path.join(here, "dictionary.png"), dpi=150, bbox_inches="tight")

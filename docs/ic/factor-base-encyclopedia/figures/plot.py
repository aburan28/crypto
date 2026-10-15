"""Lay out the factor-base encyclopedia figures from data.json.

Plots quoted values only; computes nothing. Run with any Python that has
matplotlib:  python3 plot.py   (writes taxonomy.svg, cost_vs_rho.svg,
pkm_solving_degree.svg next to this file).
"""
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
DATA = json.loads((HERE / "data.json").read_text())

VERDICT_COLOR = {
    "refuted": "#d1495b",
    "engineering": "#edae49",
    "open": "#00798c",
    "not-ECDLP": "#8d99ae",
    "below-rho-weak-class": "#66a182",
}
VERDICT_LABEL = {
    "refuted": "refuted / negative",
    "engineering": "engineering / constant factor",
    "open": "open / inconclusive",
    "not-ECDLP": "genus >= 2 Jacobian (not an ECDLP route alone)",
    "below-rho-weak-class": "below rho, weak class only",
}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "svg.fonttype": "none"})


def taxonomy():
    columns = [
        ("Prime fields E(F_p)", [
            ("0.1 Kummer g*mu_{2^t} [P-224]", "open", True),
            ("0.2 Isogeny x(R'+<Q>) [P-192]", "open", True),
            ("A1 small abscissae (baseline)", "refuted", False),
            ("A2 geometry family (FB3)", "refuted", False),
            ("A3 automorphism orbits", "engineering", False),
            ("A4 subgroup/coset, one binomial", "refuted", False),
            ("A5 Dickson; P-256 torus (24 rds)", "refuted", False),
            ("A6 P-256 affine bit-box", "refuted", False),
            ("A7 P-256 PGL2 root orbit", "refuted", False),
            ("A8 Eisenstein-smooth (j=0)", "refuted", False),
            ("A9 arith. progression; mod-l", "refuted", False),
            ("A10 fixed-degree / map images", "refuted", False),
            ("A11 isogeny-augmented, class nbrs", "refuted", False),
            ("A12 height box, degree-d maps", "open", False),
        ]),
        ("Extension fields E(F_{p^k})", [
            ("B1 Gaudry-Diem x in F_p, k=2..5", "refuted", False),
            ("B2 psi / GLS / eigenline folds", "engineering", False),
            ("B3 exotic coordinates", "engineering", False),
            ("# Binary fields E(F_{2^n})", None, False),
            ("C1 coordinate window (FPPR)", "refuted", False),
            ("C2 random subspace", "refuted", False),
            ("C3 subfield F_{2^l}", "refuted", False),
            ("C4 EXP-H shapes (G1-G6 open)", "open", False),
            ("C5 trace-zero", "engineering", False),
            ("C6 Hamming ball", "refuted", False),
            ("C7 quasi-subfield", "refuted", False),
            ("C8 GG / Shantz-Teske", "engineering", False),
            ("C9-C11 affine, adapted, quotient", "refuted", False),
        ]),
        ("Koblitz curves", [
            ("D1 Frobenius-invariant subspace", "engineering", False),
            ("D2 Frobenius-orbit unions", "engineering", False),
            ("D3 Z[tau] multiplier closure", "engineering", False),
            ("D4 cyclotomic gamma-weight", "open", False),
            ("D5 u-frame symmetrised", "refuted", False),
            ("D6 prefix Frobenius closure", "engineering", False),
            ("D7 isogeny-transported", "refuted", False),
            ("# Jacobians and covers", None, False),
            ("E1 degree-1 places, g=2..4", "not-ECDLP", False),
            ("E2 Joux-Vitse g3 cover", "below-rho-weak-class", False),
            ("E3 P-192 CM class ideals", "open", False),
            ("E4 GHS covers ECC2K-130", "refuted", False),
            ("E5 Smith quartic codomain", "refuted", False),
        ]),
    ]
    fig, ax = plt.subplots(figsize=(10, 9.2))
    ax.set_xlim(0, 3)
    ax.set_ylim(-17.2, 1.2)
    ax.axis("off")
    for ci, (title, items) in enumerate(columns):
        x = ci + 0.03
        ax.text(ci + 0.5, 0.75, title, ha="center", va="center", fontsize=10, fontweight="bold")
        for ri, (label, verdict, flagship) in enumerate(items):
            y = -ri * 1.08
            if verdict is None:
                ax.text(ci + 0.5, y, label[2:], ha="center", va="center", fontsize=10, fontweight="bold")
                continue
            box = FancyBboxPatch((x, y - 0.42), 0.94, 0.84, boxstyle="round,pad=0.02,rounding_size=0.08",
                                 fc=VERDICT_COLOR[verdict], ec="black" if flagship else "white",
                                 lw=2.2 if flagship else 0.6, alpha=0.9)
            ax.add_patch(box)
            ax.text(x + 0.47, y, label, ha="center", va="center", fontsize=8.6,
                    color="white" if verdict in ("refuted", "open") else "black",
                    fontweight="bold" if flagship else "normal", wrap=True)
    handles = [plt.Rectangle((0, 0), 1, 1, fc=VERDICT_COLOR[k]) for k in VERDICT_LABEL]
    ax.legend(handles, list(VERDICT_LABEL.values()), loc="lower left", bbox_to_anchor=(0.0, 0.0),
              fontsize=9, frameon=False, title="recorded verdict (black border = flagship families, §0)")
    ax.set_title("Factor-base families tried, by field regime\n(verdicts as recorded by each source; 2026-10-15)",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(HERE / "taxonomy.svg")
    plt.close(fig)


def cost_vs_rho():
    rows = DATA["cost_vs_rho"]
    fig, ax = plt.subplots(figsize=(8.2, 7))
    for i, r in enumerate(rows):
        y = len(rows) - 1 - i
        c = VERDICT_COLOR[r["verdict"]]
        if r["hi"] > r["lo"]:
            ax.plot([r["lo"], r["hi"]], [y, y], color=c, lw=6, solid_capstyle="butt")
        else:
            ax.plot([r["lo"]], [y], "o", color=c, ms=8)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r["label"] for r in reversed(rows)], fontsize=8)
    ax.set_xscale("log")
    ax.axvline(1.0, color="black", lw=1, ls="--")
    ax.text(1.05, len(rows) - 0.4, "rho parity", fontsize=8)
    ax.set_xlabel("recorded cost / matched rho (log scale; bar = range over the sizes the source ran)")
    ax.set_title("Best recorded row per family, as each source reported it\n"
                 "(different curve sizes and accounting per row: not a ranking)", fontsize=10)
    handles = [plt.Line2D([], [], color=VERDICT_COLOR[k], lw=6) for k in VERDICT_LABEL]
    ax.legend(handles, list(VERDICT_LABEL.values()), fontsize=7, loc="lower right")
    ax.grid(axis="x", which="major", alpha=0.3)
    fig.tight_layout()
    fig.savefig(HERE / "cost_vs_rho.svg")
    plt.close(fig)


def pkm_degree():
    d = DATA["pkm_solving_degree"]
    fig, ax = plt.subplots(figsize=(8, 5))
    markers = {2: "o", 3: "s", 4: "^"}
    for s in d["series"]:
        m = int(s["label"].split("m=")[1][0])
        color = "#00798c" if s["family"] == "kummer" else "#d1495b"
        ls = "-" if "p1" in s["label"] else "--"
        ax.plot(s["N"], s["D"], ls, marker=markers[m], color=color, label=s["label"], ms=6, alpha=0.85)
    ax.axhline(d["admissible_m4_deployed"], color="black", lw=1, ls=":")
    ax.text(12.5, d["admissible_m4_deployed"] + 0.1,
            "admissible D at a 256-bit p, m=4 (requirements README §6.4)", fontsize=8)
    ax.set_xlabel("N = m*t  (tower variables; deployed size needs N ~ log2 r ~ 192-256)")
    ax.set_ylabel("F4 solving degree D (measured)")
    ax.set_ylim(2.5, 7.6)
    ax.set_title("PKM tower bases: measured solving degree, toy primes only\n"
                 "p0 = 786433 (20-bit), p1 = 2013265921 (31-bit); one or two targets per cell", fontsize=10)
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(HERE / "pkm_solving_degree.svg")
    plt.close(fig)


if __name__ == "__main__":
    taxonomy()
    cost_vs_rho()
    pkm_degree()

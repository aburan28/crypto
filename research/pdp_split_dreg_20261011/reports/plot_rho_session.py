"""Lay out the ecbench table of session ECBS1he23809badf5b as an SVG.

Layout only: the means and 95 % intervals are ecbench's own
(`ecbench table`, saved as audits/rho-01.table.md); nothing is computed here.

    python plot_rho_session.py ../audits/rho-01.table.md OUT.svg
"""

import re
import sys

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402

ARM_STYLE = {
    "rho-strong": ("#2a78d6", "o", "rho.signed_frobenius_strong (reference)"),
    "rho-frob": ("#eb6834", "s", "rho.signed_frobenius (baseline)"),
}


def main():
    rows = []
    for line in open(sys.argv[1], encoding="utf-8"):
        cols = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cols) >= 9 and cols[0].startswith("ECBS1"):
            lo, hi = map(float, re.findall(r"[\d.]+", cols[7]))
            rows.append((cols[1], cols[3], float(cols[4]), float(cols[6]), lo, hi, cols[5]))
    curves = sorted({(r[2], r[1]) for r in rows})
    fig, ax = plt.subplots(figsize=(8, 3.8))
    for arm, (color, marker, label) in ARM_STYLE.items():
        xs, ys, los, his = [], [], [], []
        for xi, (_, slug) in enumerate(curves):
            for r in rows:
                if r[0] == arm and r[1] == slug:
                    xs.append(xi + (0.12 if arm == "rho-frob" else -0.12))
                    ys.append(r[3])
                    los.append(r[3] - r[4])
                    his.append(r[5] - r[3])
        ax.errorbar(xs, ys, yerr=[los, his], fmt=marker, color=color, ecolor=color,
                    elinewidth=1.5, capsize=3, markersize=6, label=label)
    ax.set_yscale("log")
    ax.set_xticks(range(len(curves)))
    ax.set_xticklabels([f"{slug}\nlog2 r = {lr:.2f}" for lr, slug in curves], fontsize=7)
    ax.set_ylabel("mean S per run, 95 % interval (log)")
    ax.grid(True, color="#e4e3df", linewidth=0.6)
    ax.legend(fontsize=8, frameon=False)
    ax.set_title("Session ECBS1he23809badf5b: 8 public targets x 5 rounds per curve and arm "
                 "(40 verified runs each)", fontsize=9)
    fig.tight_layout()
    fig.savefig(sys.argv[2])


if __name__ == "__main__":
    main()

"""Plot the recorded ladder (raw/summary.json, written by the Rust
`analyze` subcommand).  Plots recorded values and the fitted lines that file
holds; computes nothing new.  python3 -I plot.py"""
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
S = json.loads((HERE / "raw" / "summary.json").read_text())
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)

ARMS = [("plain", "plain k = 4 (x-line, S₄ C₄)", "#2a6fdb"),
        ("d4", "D₄ Y-line (D₄ C₄)", "#d9480f"),
        ("dlp", "double large primes (x-line, S₄ C₄)", "#2b8a3e")]

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})


def pts(arm, key):
    xs, ys = [], []
    for r in S["rows"]:
        a = r["arms"].get(arm)
        if a:
            xs.append(r["log2_E"])
            ys.append(a[key])
    return xs, ys


def line(fit, x0, x1):
    b, a = fit["slope"], fit["intercept_log2"]
    return [x0, x1], [2 ** (a + b * x0), 2 ** (a + b * x1)]


# Figure 1: S_IC / S_rho-neg against log2 #E (charged), with the LA-only floor.
fig, ax = plt.subplots(figsize=(7.2, 4.4))
for arm, label, c in ARMS:
    xs, ys = pts(arm, "ratio_vs_rho_neg")
    f = S["fits"][arm]["ratio_vs_rho_neg"]
    ax.scatter(xs, ys, s=18, color=c, label=f"{label}: slope {f['slope']:+.3f}")
    lx, ly = line(f, min(xs), max(xs))
    ax.plot(lx, ly, color=c, lw=1)
    xs2, la = pts(arm, "la_ratio_vs_rho_neg")
    ax.scatter(xs2, la, s=14, marker="x", color=c, alpha=0.7)
ax.axhline(1.0, color="k", lw=0.8, ls="--")
ax.text(27, 1.25, "rho-negation (= 1)", fontsize=8)
ax.set_yscale("log", base=10)
ax.set_xlabel("log₂ #E(F_{p⁴})   (p = 101 … 1033, h = 4)")
ax.set_ylabel("S_IC / S_rho-neg  (GAE / √r)")
ax.set_title("k = 4 ladder: charged S against rho-negation (dots), linear algebra alone (×)")
ax.legend(fontsize=7.5, loc="center right", frameon=False)
fig.tight_layout()
for ext, kw in (("svg", {}), ("png", {"dpi": 150})):
    fig.savefig(FIG / f"ratio_vs_rho.{ext}", **kw)

# Figure 2: phase costs (GAE) against log2 #E, with fitted exponents.
fig, axs = plt.subplots(1, 2, figsize=(10, 4.2), sharex=True)
for arm, label, c in ARMS:
    for ax, key, fk in ((axs[0], "rel_gae", "relation_phase"), (axs[1], "la_gae", "linear_algebra")):
        xs, ys = pts(arm, key)
        f = S["fits"][arm][fk]
        ax.scatter(xs, ys, s=16, color=c, label=f"{arm}: #E^{f['slope']/4*4:.3f}")
        lx, ly = line(f, min(xs), max(xs))
        ax.plot(lx, ly, color=c, lw=1)
rx = [r["log2_E"] for r in S["rows"]]
rc = [r["rho_neg_S"] * 2 ** (r["log2_r"] / 2) for r in S["rows"]]
for ax in axs:
    ax.scatter(rx, rc, s=10, color="#555", marker="s", label="rho-negation, whole solve")
    ax.set_yscale("log", base=2)
    ax.set_xlabel("log₂ #E")
axs[0].set_ylabel("GAE (group-addition equivalents)")
axs[0].set_title("relation phase (trials × C₄/97 + scalar mults)")
axs[1].set_title("linear algebra (mod-r mults × 16/97)")
for ax in axs:
    ax.legend(fontsize=7.5, frameon=False)
fig.tight_layout()
for ext, kw in (("svg", {}), ("png", {"dpi": 150})):
    fig.savefig(FIG / f"phases.{ext}", **kw)
print("wrote", sorted(p.name for p in FIG.iterdir()))

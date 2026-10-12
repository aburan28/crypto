"""Plot recorded W10 values (layout only; every number comes from the run JSON).

usage: python -I plot_w10.py <run.json> <out_dir>
writes delta_vs_arity.svg, delta_vs_multiplicity.svg, finite_lattice.svg
"""
import json, sys, pathlib, collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

run = json.loads(pathlib.Path(sys.argv[1]).read_text())
out = pathlib.Path(sys.argv[2])
src = sys.argv[3] if len(sys.argv) > 3 else sys.argv[1]
recs = run["records"]
R = lambda part: [r for r in recs if r.get("part") == part]

fam = R("families")
clo = [r for r in R("closure") if r.get("encoding") == "x"]
eclo = [r for r in R("closure") if r.get("encoding") == "e" and not r.get("skipped")]
known = R("known")
jmbox = {r["m"]: r for r in known if r.get("case") == "jm_basic_box"}
names = {2: "S3", 3: "S4", 4: "S5"}

# ---------- graph 1: delta vs arity ----------
fig, ax = plt.subplots(figsize=(8.6, 5.4))
ms = [2, 3, 4]
lim = [jmbox[m]["limit"] for m in ms]
ub = [next(r["upper_bound"] for r in fam if r["m"] == m) for m in ms]
best_cl = []
for m in ms:
    rs = [r for r in clo if r["m"] == m]
    best_cl.append(max(rs, key=lambda r: (r["ell"], r["delta"])))
ax.plot(ms, ub, "v--", color="#7a4cc2", label="upper bound, whole family, every ℓ: 1/(m·d)  (charging argument)")
ax.plot(ms, lim, "o-", color="#2a6fdb", label="δ*: JM box limit ℓ→∞, 2/(m(m+1)d)  (= repo's 1/6, 1/24, 1/80)")
ax.plot(ms, [r["delta"] for r in best_cl], "s:", color="#3aa76d",
        label="exact closure optimum at largest ℓ run ("
        + ", ".join(f'{names[r["m"]]}: ℓ={r["ell"]}' for r in best_cl) + ")")
need = {2: 0.5, 3: 0.25, 4: 0.2}
ax.plot(ms, [need[m] for m in ms], "D-", color="#d0453b",
        label="δ IC needs (repo §2.4: 1/2 never wins, 1/4 ties, 1/5 first win)")
ax.errorbar([2.06], [0.12], yerr=[[0.02], [0.04]], fmt="*", color="#444", capsize=4,
            label="repo S3 probes B/C (40–56-bit): HG ≤ 0.10, heuristic fades 0.12–0.16")
for m, v in zip(ms, lim):
    ax.annotate(f"{v:.4f}", (m, v), xytext=(8, -14), textcoords="offset points", fontsize=7, color="#2a6fdb")
for m, v in zip(ms, ub):
    ax.annotate(f"{v:.4f}", (m, v), xytext=(6, 4), textcoords="offset points", fontsize=7, color="#7a4cc2")
ax.set_yscale("log")
ax.set_xticks(ms)
ax.set_xticklabels(["m=2\nS3 (box [0,2]²)", "m=3\nS4 (box [0,4]³)", "m=4\nS5 (box [0,8]⁴)"])
ax.set_ylabel("δ = log_p X  (log scale)")
ax.set_title("Coppersmith reach for S_{m+1}(x₁..x_m, z) ≡ 0 mod p, |x_i| < p^δ\n"
             "full-rank triangular single-polynomial lattices (JM basic/extended, MN/cuso-style shift sets), p → ∞",
             fontsize=9.5)
ax.legend(fontsize=6.6, frameon=False, loc="lower left")
fig.text(0.01, 0.01, f"Source: {src} (lit2026 W10, examples/lit2026_w10_coppersmith_bound.rs)", fontsize=6.5, color="#555")
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(out / "delta_vs_arity.svg")

# ---------- graph 2: delta_l vs l ----------
fig, axs = plt.subplots(1, 3, figsize=(12, 4.3))
for ax, m in zip(axs, ms):
    fr = sorted([r for r in fam if r["m"] == m], key=lambda r: r["ell"])
    ax.plot([r["ell"] for r in fr], [r["jm_basic"] for r in fr], "o-", ms=3, label="JM basic M=[0,dℓ]^m")
    ax.plot([r["ell"] for r in fr], [r["jm_extended"] for r in fr], "s-", ms=3, label="JM extended, best t")
    best = {}
    for r in clo:
        if r["m"] == m and (r["ell"] not in best or r["delta"] > best[r["ell"]]["delta"]):
            best[r["ell"]] = r
    xs = sorted(best)
    ax.plot(xs, [best[x]["delta"] for x in xs], "x", color="k", ms=8, mew=1.6,
            label="exact closure optimum (all M in box, all vertices)")
    er = sorted([r for r in eclo if r["m"] == m], key=lambda r: r["ell"])
    if er:
        ax.plot([r["ell"] for r in er], [r["delta"] for r in er], "+", color="#d0453b", ms=10, mew=1.6,
                label="closure, symmetric e-encoding")
    ax.axhline(jmbox[m]["limit"], color="#2a6fdb", ls="--", lw=0.8, label="limit 2/(m(m+1)d)")
    ax.axhline(fr[0]["upper_bound"], color="#7a4cc2", ls=":", lw=0.8, label="upper bound 1/(md)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("multiplicity ℓ")
    ax.set_title(f"{names[m]} (m={m}, d={2 ** (m - 1)})", fontsize=9)
axs[0].set_ylabel("δ_ℓ (p → ∞ at fixed ℓ)")
axs[2].legend(fontsize=6.3, frameon=False, loc="lower right")
fig.suptitle("Convergence in multiplicity ℓ (closed-form box families to ℓ = 128; exact closure optimum where run)", fontsize=10)
fig.text(0.01, 0.01, f"Source: {src}", fontsize=6.5, color="#555")
fig.tight_layout(rect=(0, 0.04, 1, 0.95))
fig.savefig(out / "delta_vs_multiplicity.svg")

# ---------- graph 3: finite lattices ----------
lat = R("lattice")
if lat:
    fig, ax = plt.subplots(figsize=(8, 4.4))
    groups = collections.defaultdict(list)
    for r in lat:
        groups[(r["poly"], r["ell"])].append(r)
    for (poly, ell), rs in sorted(groups.items()):
        rs.sort(key=lambda r: r["delta"])
        dim = max(r["dim"] for r in rs)
        line, = ax.plot([r["delta"] for r in rs], [r["successes"] / max(r["trials"], 1) for r in rs], "o-",
                        label=f"{poly}, ℓ={ell}, dim {dim} (≤{max(r['trials'] for r in rs)} trials/cell): ≥m rows vanish")
        ax.plot([r["delta"] for r in rs], [r["hg_first_m"] / max(r["trials"], 1) for r in rs], ".:",
                color=line.get_color(), alpha=0.7, label=f"{poly}, ℓ={ell}: HG met on first m rows")
    ax.axvline(0.2, color="#d0453b", lw=0.8, ls="--", label="δ = 1/5 (needed at m=4)")
    ax.set_xlabel("δ (planted x(P_i) ≤ p^δ, p ≈ 2^32)")
    ax.set_ylabel("fraction of trials")
    ax.set_title(f"Finite JM lattices with the repo LLL on {lat[0]['curve']} (p = {lat[0]['p']})", fontsize=9.5)
    ax.legend(fontsize=6.3, frameon=False)
    fig.text(0.01, 0.01, f"Source: {src}", fontsize=6.5, color="#555")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(out / "finite_lattice.svg")
print("ok")

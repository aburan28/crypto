"""Plot the recorded W11 values (layout only: every number is read from the
run's JSONL/TSV files written by the Rust examples; nothing is computed here
beyond medians/min/max of recorded samples for display).

usage: python -I plot.py <run_dir>
"""
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402

run = Path(sys.argv[1])
raw = run / "raw"
plt.rcParams.update({"font.size": 8, "svg.fonttype": "none"})
RID = run.name


def jl(name):
    p = raw / name
    if not p.exists():
        return []
    return [json.loads(l) for l in p.read_text().splitlines() if l.startswith("{")]


def save(fig, name):
    fig.tight_layout()
    fig.savefig(run / name)
    plt.close(fig)


# 1. BSGS counted ops per visited point
rows = [r for r in jl("bsgs_opcount.jsonl") if r.get("label", "").startswith("secp48")]
if rows:
    fig, axs = plt.subplots(1, 2, figsize=(8, 3.2))
    for ax, phase in zip(axs, ["baby", "giant"]):
        for kern, mk in [("paired", "o-"), ("coz", "s-")]:
            rs = sorted([r for r in rows if r["phase"] == phase and r["kernel"] == kern], key=lambda r: r["pair_k"])
            ax.plot([r["pair_k"] for r in rs], [r["mul_per_point"] + r["sqr_per_point"] for r in rs], mk, label=f"{kern}: M+S per point")
        ch = [r for r in rows if r["phase"] == phase and r["kernel"] == "chain"][0]
        ax.axhline(ch["mul_per_point"] + ch["sqr_per_point"], color="k", ls="--", label="chain (original): M+S per point")
        ax.set_xscale("log", base=2)
        ax.set_xlabel("pair window K (points per round = 2K+1)")
        ax.set_ylabel("field M+S per visited point (counted)")
        ax.set_title(f"{phase} phase, secp48k1 ladder curve, W=8 chains")
        ax.legend(fontsize=6)
    fig.suptitle(f"W11 (a) BSGS step kernels: counted field multiplications+squarings per point — {RID}", fontsize=8)
    save(fig, "bsgs_opcount.svg")

# 2. BSGS wall time per arm
rows = jl("bsgs_time.jsonl")
if rows:
    by = defaultdict(list)
    lab = {}
    for r in rows:
        by[(r["bits"], r["arm"])].append(r["total_s"])
        lab[r["bits"]] = r["label"]
    bits = sorted({b for b, _ in by})
    arms = ["original", "previous", "paired", "coz"]
    fig, ax = plt.subplots(figsize=(7, 3.4))
    w = 0.2
    for i, a in enumerate(arms):
        med = [statistics.median(by[(b, a)]) if by[(b, a)] else 0 for b in bits]
        lo = [m - min(by[(b, a)]) if by[(b, a)] else 0 for b, m in zip(bits, med)]
        hi = [max(by[(b, a)]) - m if by[(b, a)] else 0 for b, m in zip(bits, med)]
        ax.bar([j + (i - 1.5) * w for j in range(len(bits))], med, w, yerr=[lo, hi], capsize=2, label=a)
    ax.set_yscale("log")
    ax.set_xticks(range(len(bits)))
    ax.set_xticklabels([f"{lab[b]}\n(n < 2^{b})" for b in bits])
    ax.set_ylabel("wall seconds per solve batch (median, min–max bars)")
    n = len(rows) // (len(bits) * len(arms)) if bits else 0
    ax.set_title(f"W11 (a) BSGS full-group solve, L0 interleaved, {n} rounds — {RID}", fontsize=8)
    ax.legend(fontsize=6)
    save(fig, "bsgs_wall.svg")

# 3. binary x-only evaluators
rows = jl("binary_sect163.jsonl")
if rows:
    names = ["batch_inverse", "kernel_poly_given_h", "fraction", "fraction_multi", "product_normalised", "projective"]
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.4))
    for li, ell in enumerate([3, 5, 7]):
        rs = [r for r in rows if r["ell"] == ell]
        if not rs:
            continue
        c = rs[0]["counts_two_points"]
        axs[0].bar([j + (li - 1) * 0.27 for j in range(len(names))], [c[n]["mul"] for n in names], 0.27, label=f"ℓ={ell}")
        tm = [statistics.median([r["time_two_points"][n]["median_us"] for r in rs]) for n in names]
        axs[1].bar([j + (li - 1) * 0.27 for j in range(len(names))], tm, 0.27, label=f"ℓ={ell}")
    for ax, yl in zip(axs, ["F_2^163 multiplications per step (2 points, counted)", "µs per step (2 points), median over steps"]):
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels([n.replace("_", "\n") for n in names], fontsize=6)
        ax.set_ylabel(yl)
        ax.legend(fontsize=6)
    fig.suptitle(f"W11 (b) x-only binary isogeny evaluation over F_2^163 (inversions: 2,2,2,1,2,0 per step) — {RID}", fontsize=8)
    save(fig, "binary_xonly.svg")

# 4. biextension vs Miller
rows = jl("biext_k2.jsonl")
rows = [r for r in rows if "r_bits" in r]
if rows:
    by = defaultdict(list)
    for r in rows:
        by[r["r_bits"]].append(r)
    xs = sorted(by)
    fig, ax = plt.subplots(figsize=(6, 3.3))
    for key, lab in [("ops_miller_affine", "Miller, affine (repo port; + inversions)"), ("ops_miller_jacobian", "Miller, Jacobian, denominator-free"), ("ops_cubical", "cubical ladder (loop)")]:
        ax.plot(xs, [statistics.median([r[key]["m_plus_s"] for r in by[x]]) for x in xs], "o-", label=lab)
    ax.plot(xs, [statistics.median([r["ops_cubical"]["m_plus_s"] + r["ops_unsquare"]["m_plus_s"] for r in by[x]]) for x in xs], "s--", label="cubical + unsquare (t² → t)")
    ax.set_xlabel("log2 r (bits of the pairing order)")
    ax.set_ylabel("F_p M+S per pairing (loop, counted)")
    ax.set_title(f"W11 (c) k=2 Tate pairing on y²=x³+x: identical reduced values — {RID}", fontsize=8)
    ax.legend(fontsize=6)
    save(fig, "biext_k2.svg")

# 5. DP statistics
rows = jl("dp_stats.jsonl")
if rows:
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.4))
    for r in rows:
        lab = f'{r["curve"]} dp={r["dp_bits"]}'
        x = [0.25 * (i + 0.5) for i in range(len(r["hist"]))]
        tot = sum(r["hist"])
        axs[0].plot(x, [h / tot for h in r["hist"]], ".-", label=lab)
    r0 = rows[0]
    tot = sum(r0["hist_model_geo"])
    axs[0].plot([0.25 * (i + 0.5) for i in range(len(r0["hist_model_geo"]))], [h / tot for h in r0["hist_model_geo"]], "k--", label="geometric model")
    axs[0].set_yscale("log")
    axs[0].set_xlabel("DP trail length × θ (bins of 0.25/θ; last bin ≥ 8/θ)")
    axs[0].set_ylabel("fraction of DP-terminated walks")
    axs[0].legend(fontsize=5)
    labs = [f'{r["curve"]}\ndp={r["dp_bits"]}' for r in rows]
    axs[1].plot(range(len(rows)), [r["hit_rate_over_theta"] for r in rows], "o", label="DP hits per step / θ")
    axs[1].plot(range(len(rows)), [r["trail_mean_over_geo"] for r in rows], "s", label="mean trail / geometric mean")
    axs[1].plot(range(len(rows)), [r["trail_mean_over_cr"] for r in rows], "^", label="mean trail / competing-risk mean")
    axs[1].axhline(1.0, color="k", lw=0.6)
    axs[1].set_xticks(range(len(rows)))
    axs[1].set_xticklabels(labs, fontsize=5)
    axs[1].legend(fontsize=6)
    fig.suptitle(f"W11 (d) signed-Frobenius strong rho: distinguished-point statistics vs the geometric model — {RID}", fontsize=8)
    save(fig, "dp_stats.svg")

# 6. challenge ladder
rows = [r for r in jl("ladder_solve.jsonl") if "ratio" in r]
if rows:
    by = defaultdict(list)
    for r in rows:
        by[r["bits"]].append(r["ratio"])
    xs = sorted(by)
    fig, ax = plt.subplots(figsize=(6.5, 3.3))
    for x in xs:
        ax.plot([x] * len(by[x]), by[x], ".", color="#888")
    ax.plot(xs, [statistics.mean(by[x]) for x in xs], "o-", label="mean over seeds")
    ax.axhline(1.0, color="k", lw=0.6, label="√(πn/4)")
    ax.set_xlabel("challenge level k (secpkk1, Algorithm 3 descending)")
    ax.set_ylabel("rho iterations / √(πn/4)")
    ax.set_title(f"W11 (e) prime_fast rho (negation map) on the arXiv 2508.14011 ladder — {RID}", fontsize=8)
    ax.legend(fontsize=6)
    save(fig, "challenge_ladder.svg")
print("plotted", RID)

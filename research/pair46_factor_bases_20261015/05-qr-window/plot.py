"""Layout only: plots values recorded in results.json (no result is computed here)."""
import json, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = sys.argv[1] if len(sys.argv) > 1 else "."
r = json.load(open(f"{d}/results.json"))
lab = ["A0 plain", "A1 qr_window\nchi(f(x))=1", "A2 spec\nchi(x)!=-1", "A3 random\nhalf", "A4 twist\nchi(f(x))=-1"]
col = ["#7f8c8d", "#c0392b", "#2e86c1", "#8e44ad", "#27ae60"]

# Figure 1: point-density gain G = M1(arm)/M1(A0)
fig, ax = plt.subplots(figsize=(9, 4.2))
series = []
for c in r["pair46"]["curves"]:
    can = c["canonical"]
    ci = [None, can["G_A1_ci95"], can["G_A2_ci95"], can["G_A3_ci95"], None]
    series.append((f"pair46 {c['name']} [0,B)", can["G"], ci))
    series.append((f"pair46 {c['name']} 32x2^16 windows", [g["point"] for g in c["windows"]["G"]], [g["ci95"] for g in c["windows"]["G"]]))
tg = r["toy"]["groups"]["all"]["G_M1"]
series.append(("toy 200 curves (20-30 bit)", [g["point"] for g in tg], [g["ci95"] for g in tg]))
mk = ["o", "s", "^", "D", "v"]
for si, (name, pts, cis) in enumerate(series):
    for a in range(5):
        x = a + (si - 2) * 0.12
        y = pts[a]
        if cis[a] is not None:
            lo, hi = cis[a]
            ax.errorbar(x, y, yerr=[[y - lo], [hi - y]], fmt=mk[si], color=col[a], ms=5, capsize=2, label=name if a == 1 else None, mec="k", mew=0.4)
        else:
            ax.plot(x, y, mk[si], color=col[a], ms=5, mec="k", mew=0.4, label=name if a == 1 else None)
ax.axhline(2.0, ls="--", color="#c0392b", lw=0.8)
ax.axhline(1.0, ls=":", color="k", lw=0.8)
ax.text(4.45, 2.02, "predicted 2", fontsize=7, color="#c0392b", ha="right")
ax.set_xticks(range(5)); ax.set_xticklabels(lab, fontsize=8)
ax.set_ylabel("point-density gain G = M1(arm)/M1(A0)")
ax.set_title("05-qr-window: only the twist-aware filter doubles point density (95% intervals)", fontsize=10)
ax.legend(fontsize=7, loc="center right", title="marker = data set", title_fontsize=7)
ax.set_ylim(-0.1, 2.3)
fig.tight_layout()
for e in ("svg", "png"):
    fig.savefig(f"{d}/fig_density_gain.{e}", dpi=160)

# Figure 2: toy relations relative to A0
fig, axs = plt.subplots(1, 2, figsize=(9, 3.8))
for ax, key, title in [(axs[0], "ratio_M2", "relation rate M2 (>=1 decomposition) / A0"), (axs[1], "ratio_D", "decompositions per target D / A0")]:
    for gi, g in enumerate(["F", "N2"]):
        v = r["toy"]["groups"][g][key]
        for a in range(5):
            x = a + (gi - 0.5) * 0.25; y = v[a]["point"]; lo, hi = v[a]["ci95"]
            ax.errorbar(x, y, yerr=[[y - lo], [hi - y]], fmt="o" if g == "F" else "s", color=col[a], mec="k", mew=0.4, capsize=2,
                        label=("F family (j=-96^3)" if g == "F" else "N2 random prime order") if a == 0 else None)
    ax.axhline(1.0, ls=":", color="k", lw=0.8); ax.axhline(0.25, ls="--", color="#2e86c1", lw=0.6)
    ax.set_xticks(range(5)); ax.set_xticklabels([l.split("\n")[0] for l in lab], fontsize=7, rotation=20)
    ax.set_title(title, fontsize=9); ax.set_ylim(-0.05, 1.15)
axs[0].legend(fontsize=7)
fig.suptitle(f"toy Semaev m=2 at B=ceil(sqrt(2N)), T={r['toy']['T']} targets/curve; A1 identical to A0 on every target: {r['toy']['E2_all']}", fontsize=9)
fig.tight_layout()
for e in ("svg", "png"):
    fig.savefig(f"{d}/fig_toy_relations.{e}", dpi=160)

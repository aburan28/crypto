"""Summarise raw.json (written by the Rust instrument) into results.json and draw figures.
Summary statistics (means, t-intervals) of recorded values only; no new computation of counts."""
import json, math, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d0 = sys.argv[1] if len(sys.argv) > 1 else "."
raw = json.load(open(f"{d0}/raw.json"))
rows = raw["rows"]
T975 = {2: 12.706, 3: 4.303, 4: 3.182, 5: 2.776, 6: 2.571, 7: 2.447, 8: 2.365, 9: 2.306}

def ci(v):
    n = len(v); m = sum(v) / n
    if n < 2: return m, m, m
    s = math.sqrt(sum((x - m) ** 2 for x in v) / (n - 1))
    h = T975.get(n, 1.96) * s / math.sqrt(n)
    return m, m - h, m + h

summ = {}
for fam in ["RP", "RR"]:
    for d in range(2, 9):
        r = [x for x in rows if x["family"] == fam and x["d"] == d]
        rho = [x["rho"] for x in r]; delta = [x["rho"] - x["c_d"] for x in r]
        m, lo, hi = ci(rho); dm, dlo, dhi = ci(delta)
        summ[f"{fam}_d{d}"] = dict(n=len(r), rho_mean=m, rho_ci95=[lo, hi], c_d=r[0]["c_d"],
                                  delta_mean=dm, delta_ci95=[dlo, dhi],
                                  genus=sorted({x["genus"] for x in r}),
                                  min_exponent=min(x["exponent"] for x in r))
alg = [x for x in rows if x["family"] in ("RP", "RR", "PW", "DK")]
ext = [x for x in rows if x["family"] == "EXT"]
checks = dict(
    LB_all_ok=all(x["LB_ok"] in (True, None) for x in alg),
    LB_rows_checked=sum(1 for x in alg if x["LB_ok"] is not None),
    ZD_UB_all_ok=all(x["UB_ok"] for x in rows if x["family"] == "ZD"),
    ZD_max_F=max(x["F"] for x in rows if x["family"] == "ZD"),
    controls_ok=all(x.get("control_ok", True) for x in rows),
    ext_valid=all(x["control_valid"] for x in ext),
    fields_ok=all(f["prime_mr12"] and f["hasse_ok"] and f["order_check_4pts"] and f["p_mod4"] == 3 for f in raw["fields"]),
    squarefree_rows=sum(1 for x in alg if x["h_squarefree"]), alg_rows=len(alg),
)
min_alg_e = min(x["exponent"] for x in alg)
max_ext_e = max(x["exponent"] for x in ext)
pw_min = min((x for x in rows if x["family"] == "PW"), key=lambda x: x["rho"])
rf = [x["rho"] for x in rows if x["family"] == "RF"]
p46 = [x for x in rows if x["family"] == "P46"]
out = dict(protocol_sha256=raw["protocol_sha256"], instrument_wall_s=raw["wall_s"], fields=raw["fields"],
           checks=checks, summary_random_families=summ,
           headline=dict(min_exponent_algebraic=min_alg_e, max_exponent_ext_control=max_ext_e,
                         exponent_gap=min_alg_e - max_ext_e,
                         min_rho_algebraic=min(x["rho"] for x in alg),
                         min_rho_row=dict(family=pw_min["family"], k=pw_min["k"], d=pw_min["d"], rho=pw_min["rho"]),
                         RF_rho_mean_ci=ci(rf)),
           ext_control=ext, pair46=p46, rows=rows)
json.dump(out, open(f"{d0}/results.json", "w"), indent=1)

# Figure 1: rho vs d per family, with c_d and 1/d, and ext control
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ds = list(range(2, 9))
cols = {"RP": "#1f77b4", "RR": "#ff7f0e", "PW": "#2ca02c", "DK": "#9467bd"}
off = {"RP": -0.15, "RR": -0.05, "PW": 0.05, "DK": 0.15}
for fam, c in cols.items():
    for k, mk in [(16, "o"), (20, "s"), (24, "^")]:
        r = [x for x in rows if x["family"] == fam and x["k"] == k]
        ax[0].scatter([x["d"] + off[fam] for x in r], [x["rho"] for x in r], color=c, marker=mk, s=14,
                      label=f"{fam}" if k == 16 else None, alpha=0.8)
ax[0].plot(ds, [summ[f"RP_d{d}"]["c_d"] for d in ds], "k--", lw=1, label="c_d (generic S_d)")
ax[0].plot(ds, [1 / d for d in ds], "k:", lw=1, label="1/d (X1 fiber bound scale)")
p46d = [x["d"] for x in p46]; p46m = [x["rho_hat"] for x in p46]
ax[0].errorbar([d + 0.25 for d in p46d], p46m, yerr=[[m - x["ci95"][0] for m, x in zip(p46m, p46)], [x["ci95"][1] - m for m, x in zip(p46m, p46)]],
               fmt="D", color="#d62728", ms=4, label="P46 53-bit (MC, 95% CI)")
ax[0].axhline(sum(rf) / len(rf), color="grey", lw=0.8, label="RF random function")
ax[0].set_xlabel("map degree d"); ax[0].set_ylabel("rho = |F| / p"); ax[0].set_ylim(0, 1.05)
ax[0].set_title("Factor-base density, exact enumeration\n(markers o/s/^ = p~2^16/2^20/2^24; 3 seeds RP/RR)", fontsize=9)
ax[0].legend(fontsize=7, ncol=2)
# right: exponent vs log2 field size
for fam, c in cols.items():
    r = [x for x in rows if x["family"] == fam]
    ax[1].scatter([math.log2(x["p"]) + off[fam] for x in r], [x["exponent"] for x in r], color=c, s=10, label=fam)
zd = [x for x in rows if x["family"] == "ZD" and x["F"] > 0]
ax[1].scatter([math.log2(x["p"]) for x in zd], [x["exponent"] for x in zd], color="#8c564b", marker="x", s=14, label="ZD zero-dim (|F|<=2d)")
ax[1].scatter([math.log2(x["Q"]) for x in ext], [x["exponent"] for x in ext], color="#d62728", marker="*", s=90, label="EXT x in F_q over F_q^2")
ax[1].axhline(0.5, color="k", ls=":", lw=1); ax[1].text(16, 0.52, "|F| = sqrt(field)", fontsize=7)
ax[1].set_xlabel("log2(field size)"); ax[1].set_ylabel("log|F| / log(field size)")
ax[1].set_title("Exponent of |F|: prime-field families vs extension control", fontsize=9)
ax[1].legend(fontsize=7, loc="lower center", bbox_to_anchor=(0.5, 0.13), ncol=3); ax[1].set_ylim(0, 1.05)
fig.suptitle("10-bounded-degree: bounded-degree algebraic factor bases (protocol sha256 %s...)" % raw["protocol_sha256"][:12], fontsize=9)
fig.tight_layout()
fig.savefig(f"{d0}/fig_density.svg"); fig.savefig(f"{d0}/fig_density.png", dpi=150)

# Figure 2: P46 membership cost
fig, ax = plt.subplots(figsize=(5.5, 3.6))
ax.plot(p46d, [x["muls_per_membership"] for x in p46], "o-", color="#1f77b4")
ax.set_xlabel("degree d"); ax.set_ylabel("F_p mults per membership test (counted)")
ax.set_title("P46 (53-bit) membership cost: gcd(t^p - t, f - x)\nn = %d tests per d" % p46[0]["n"], fontsize=9)
fig.tight_layout(); fig.savefig(f"{d0}/fig_membership.svg"); fig.savefig(f"{d0}/fig_membership.png", dpi=150)
print(json.dumps(dict(checks=checks, headline=out["headline"]), indent=1))

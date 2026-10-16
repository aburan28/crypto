#!/usr/bin/env python3
"""Tabulate, fit and plot the recorded run (item 7). Usage: analyze.py <run_dir>"""
import glob, json, math, os, sys
import numpy as np

run = sys.argv[1]
recs = []
for f in sorted(glob.glob(os.path.join(run, "*.jsonl"))):
    for line in open(f):
        recs.append(json.loads(line))

rows = []  # one row per solve
inst_rows = []
for r in recs:
    inst_rows.append(dict(bits=r["bits"], q=r["q"], curve=r["curve"], N=r["N"], ell=r["ell"], c=r["c"],
                          **{k: v for k, v in r["checks"].items()}))
    for x in r["results"]:
        tc = x.get("total_counters", x.get("counters"))
        rows.append(dict(bits=r["bits"], q=r["q"], curve=r["curve"], ell=r["ell"], arm=x["arm"],
                         correct=x["correct"], cost=x["cost"], cost_cl9=x["cost_cl9"],
                         M=tc[0], I=tc[1], Ml=tc[2], Il=tc[3], gops=tc[4], steps=x["steps"],
                         wall=x["wall_total"],
                         nfb=x.get("nfb"), rows=x.get("rows"), restarts=x.get("restarts"),
                         full=x.get("full"), lp1=x.get("lp1"), lp2=x.get("lp2"), cycles=x.get("cycles"),
                         la_cost=(tc[2] if x["arm"] != "rho" else 0)))

# dedupe instance rows (ic and rho files both carry the instance)
seen = set(); inst_u = []
for d in inst_rows:
    k = (d["bits"], d["curve"])
    if k not in seen:
        seen.add(k); inst_u.append(d)
inst_u.sort(key=lambda d: (d["bits"], d["curve"]))

with open(os.path.join(run, "solves.csv"), "w") as f:
    keys = list(rows[0].keys())
    f.write(",".join(keys) + "\n")
    for d in sorted(rows, key=lambda d: (d["bits"], d["arm"], d["curve"])):
        f.write(",".join(str(d[k]) for k in keys) + "\n")
with open(os.path.join(run, "instances.csv"), "w") as f:
    keys = sorted({k for d in inst_u for k in d}, key=lambda k: list(inst_u[0].keys()).index(k) if k in inst_u[0] else 99)
    f.write(",".join(keys) + "\n")
    for d in inst_u:
        f.write(",".join(str(d.get(k, "")) for k in keys) + "\n")


def curve_points(arm, field="cost", transform=None):
    """per (bits, curve): mean of field over targets."""
    pts = {}
    for d in rows:
        if d["arm"] != arm:
            continue
        v = d[field] if transform is None else transform(d)
        pts.setdefault((d["bits"], d["curve"], d["q"]), []).append(v)
    out = [(q, b, c, float(np.mean(v))) for (b, c, q), v in pts.items()]
    return sorted(out)


def fit(pts):
    x = np.log([p[0] for p in pts]); y = np.log([p[3] for p in pts])
    A = np.vstack([x, np.ones_like(x)]).T
    (e, a), *_ = np.linalg.lstsq(A, y, rcond=None)
    return e, a


rng = np.random.default_rng(20261014)


def boot(pts, B=2000, fn=fit):
    bybits = {}
    for p in pts:
        bybits.setdefault(p[1], []).append(p)
    out = []
    for _ in range(B):
        s = []
        for b, ps in bybits.items():
            idx = rng.integers(0, len(ps), len(ps))
            s += [ps[i] for i in idx]
        out.append(fn(s))
    return out


summary = {}
arms = ["ic2lp", "basic", "rho"]
variants = {
    "primary": dict(field="cost"),
    "cl9": dict(field="cost_cl9"),
    "gops": dict(field="gops"),
}
fits = {}
for vname, v in variants.items():
    for arm in arms:
        pts = curve_points(arm, v["field"])
        if len(pts) < 3:
            continue
        e, a = fit(pts)
        bs = boot(pts)
        es = sorted(x[0] for x in bs)
        fits[(vname, arm)] = (e, a, bs)
        summary[f"{vname}:{arm}"] = dict(exponent=e, ci95=[es[int(0.025 * len(es))], es[int(0.975 * len(es))]],
                                         const=math.exp(a), n_points=len(pts),
                                         bits=sorted({p[1] for p in pts}))

# rho exponent in ell (instrument control): steps vs ell
pts_ell = {}
for d in rows:
    if d["arm"] == "rho":
        pts_ell.setdefault((d["bits"], d["curve"]), []).append((d["ell"], d["steps"]))
xs = np.log([v[0][0] for v in pts_ell.values()]); ys = np.log([np.mean([s for _, s in v]) for v in pts_ell.values()])
e_ell = np.polyfit(xs, ys, 1)[0]
summary["rho_steps_exponent_in_ell"] = e_ell
allrho = [d for d in rows if d["arm"] == "rho"]
summary["rho_steps_over_sqrt_ell_mean"] = float(np.mean([d["steps"] / math.sqrt(d["ell"]) for d in allrho]))
summary["rho_steps_over_sqrt_ell_sd"] = float(np.std([d["steps"] / math.sqrt(d["ell"]) for d in allrho], ddof=1))
summary["rho_cost_per_step_mean"] = float(np.mean([d["cost"] / d["steps"] for d in allrho if d["bits"] >= 12]))

# per-size table
table = []
for b in sorted({d["bits"] for d in rows}):
    t = dict(bits=b, q=[d["q"] for d in rows if d["bits"] == b][0])
    for arm in arms:
        ds = [d for d in rows if d["bits"] == b and d["arm"] == arm]
        if not ds:
            continue
        t[f"{arm}_n"] = len(ds)
        t[f"{arm}_correct"] = sum(d["correct"] for d in ds)
        t[f"{arm}_cost_mean"] = float(np.mean([d["cost"] for d in ds]))
        t[f"{arm}_cost_cl9_mean"] = float(np.mean([d["cost_cl9"] for d in ds]))
        t[f"{arm}_gops_mean"] = float(np.mean([d["gops"] for d in ds]))
        t[f"{arm}_wall_mean"] = float(np.mean([d["wall"] for d in ds]))
        if arm != "rho":
            t[f"{arm}_la_share"] = float(np.mean([d["Ml"] / d["cost"] for d in ds]))
            t[f"{arm}_nfb"] = float(np.mean([d["nfb"] for d in ds]))
            t[f"{arm}_steps"] = float(np.mean([d["steps"] for d in ds]))
    if "ic2lp_cost_mean" in t and "rho_cost_mean" in t:
        # ratio of per-curve means, bootstrap over curves
        ic = {d["curve"]: d["cost"] for d in rows if d["bits"] == b and d["arm"] == "ic2lp"}
        rh = {}
        for d in rows:
            if d["bits"] == b and d["arm"] == "rho":
                rh.setdefault(d["curve"], []).append(d["cost"])
        cs = sorted(set(ic) & set(rh))
        rat = []
        for _ in range(2000):
            idx = rng.choice(cs, len(cs))
            rat.append(np.mean([ic[c] for c in idx]) / np.mean([np.mean(rh[c]) for c in idx]))
        rat.sort()
        t["ratio_ic_rho"] = t["ic2lp_cost_mean"] / t["rho_cost_mean"]
        t["ratio_ci95"] = [rat[50], rat[1949]]
        t["ratio_ic_rho_neg"] = t["ic2lp_cost_mean"] / (t["rho_cost_mean"] / math.sqrt(2))
        t["ratio_ic_rho_cl9"] = t["ic2lp_cost_cl9_mean"] / t["rho_cost_cl9_mean"]
    table.append(t)
# rho model: measured per-step cost x ideal random-mapping step count (strongest rho reference)
for t in table:
    ds = [d for d in rows if d["bits"] == t["bits"] and d["arm"] == "rho"]
    ells = {d["curve"]: d["ell"] for d in ds}
    if ds and "ic2lp_cost_mean" in t:
        pstep = np.mean([d["cost"] / d["steps"] for d in ds])
        mean_sqrt = np.mean([math.sqrt(e) for e in ells.values()])
        t["rho_model_cost"] = pstep * math.sqrt(math.pi / 2) * mean_sqrt
        t["rho_model_neg_cost"] = pstep * math.sqrt(math.pi / 4) * mean_sqrt
        t["ratio_ic_rhomodel"] = t["ic2lp_cost_mean"] / t["rho_model_cost"]
        t["ratio_ic_rhomodel_neg"] = t["ic2lp_cost_mean"] / t["rho_model_neg_cost"]
summary["table"] = table
mt = [t for t in table if "rho_model_cost" in t]
em, am = np.polyfit(np.log([t["q"] for t in mt]), np.log([t["rho_model_cost"] for t in mt]), 1)
summary["rho_model_exponent"] = em
ei, ai, _ = fits[("primary", "ic2lp")]
summary["crossover_log2q_vs_rho_model"] = math.log2(math.exp((ai - am) / (em - ei)))
summary["crossover_log2q_vs_rho_model_neg"] = math.log2(math.exp((ai - am + 0.5 * math.log(2)) / (em - ei)))
# same with IC fit on top-5 sizes

# crossover from fits
def cross(va, vb, shift_b=0.0):
    ea, aa, _ = fits[va]; eb, ab, _ = fits[vb]
    if abs(ea - eb) < 1e-9:
        return None
    lq = (ab + shift_b - aa) / (ea - eb)
    return math.log2(math.exp(lq))

summary["crossover_log2q_primary"] = cross(("primary", "ic2lp"), ("primary", "rho"))
summary["crossover_log2q_rho_negation"] = cross(("primary", "ic2lp"), ("primary", "rho"), -0.5 * math.log(2))
summary["crossover_log2q_cl9"] = cross(("cl9", "ic2lp"), ("cl9", "rho"))
# bootstrap crossover (paired over curves is not possible: arms independent -> combine independent bootstrap draws)
bi = fits[("primary", "ic2lp")][2]; br = fits[("primary", "rho")][2]
cx = []
for (e1, a1), (e2, a2) in zip(bi, br):
    if e2 - e1 > 1e-6:
        cx.append(math.log2(math.exp((a1 - a2) / (e2 - e1))))
    else:
        cx.append(float("inf"))
cx.sort()
summary["crossover_log2q_primary_ci95"] = [cx[int(0.025 * len(cx))], cx[int(0.975 * len(cx))]]
summary["crossover_boot_frac_finite"] = sum(1 for c in cx if math.isfinite(c)) / len(cx)
# measured crossover: first size whose IC/rho ratio < 1
meas = [t for t in table if t.get("ratio_ic_rho", 9) < 1]
summary["measured_first_bits_ic_below_rho"] = meas[0]["bits"] if meas else None
# local exponents (consecutive sizes) for IC
loc = []
for t1, t2 in zip(table, table[1:]):
    if "ic2lp_cost_mean" in t1 and "ic2lp_cost_mean" in t2:
        dq = math.log(t2["q"] / t1["q"])
        loc.append(dict(bits=f"{t1['bits']}-{t2['bits']}",
                        ic2lp=math.log(t2["ic2lp_cost_mean"] / t1["ic2lp_cost_mean"]) / dq,
                        rho=math.log(t2["rho_cost_mean"] / t1["rho_cost_mean"]) / dq if "rho_cost_mean" in t2 else None))
summary["local_exponents"] = loc
# IC fit restricted to upper half sizes
pts = curve_points("ic2lp")
bmax = max(p[1] for p in pts)
up = [p for p in pts if p[1] >= bmax - 4]
e_up = fit(up)[0]
es = sorted(x[0] for x in boot(up))
summary["ic2lp_exponent_top5_sizes"] = dict(exponent=e_up, ci95=[es[50], es[1949]], bits=sorted({p[1] for p in up}))
ptsr = curve_points("rho"); upr = [p for p in ptsr if p[1] >= bmax - 4]
summary["rho_exponent_top5_sizes"] = fit(upr)[0]
e5, a5 = fit(up)
summary["crossover_log2q_top5_vs_rho_model"] = math.log2(math.exp((a5 - am) / (em - e5))) if em > e5 else None
summary["crossover_log2q_top5_vs_rho_model_neg"] = math.log2(math.exp((a5 - am + 0.5 * math.log(2)) / (em - e5))) if em > e5 else None
summary["crossover_log2q_top5_vs_rho_measured"] = math.log2(math.exp((a5 - fits[("primary", "rho")][1]) / (fits[("primary", "rho")][0] - e5)))
# random-multiplier rho control (deviation D2)
ctrl = os.path.join(os.path.dirname(run.rstrip("/")), "rho-randmult-20261014")
if os.path.isdir(ctrl):
    cv = []
    for f in glob.glob(os.path.join(ctrl, "rho_*.jsonl")):
        for line in open(f):
            r = json.loads(line)
            cv += [(r["bits"], x["steps"] / math.sqrt(r["ell"]), x["correct"]) for x in r["results"]]
    fib = [(d["bits"], d["steps"] / math.sqrt(d["ell"])) for d in rows if d["arm"] == "rho"]
    bset = sorted({b for b, _, _ in cv})
    fibm = [v for b, v in fib if b in bset]
    summary["rho_control_random_mult"] = dict(
        bits=bset, n=len(cv), correct=sum(c for _, _, c in cv),
        steps_over_sqrt_ell_mean=float(np.mean([v for _, v, _ in cv])),
        se=float(np.std([v for _, v, _ in cv], ddof=1) / math.sqrt(len(cv))),
        fibonacci_same_bits_mean=float(np.mean(fibm)),
        fibonacci_same_bits_se=float(np.std(fibm, ddof=1) / math.sqrt(len(fibm))), n_fib=len(fibm))
json.dump(summary, open(os.path.join(run, "summary.json"), "w"), indent=1, default=float)
print(json.dumps({k: v for k, v in summary.items() if k != "table"}, indent=1, default=float))
for t in table:
    print(t)

# ------------------------------------------------------------------ plot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10, "svg.fonttype": "none"})
colors = {"ic2lp": "#1f6fb4", "basic": "#8c8c8c", "rho": "#c2410c"}
labels = {"ic2lp": "genus-3 IC, double large primes (after transfer)",
          "basic": "genus-3 IC, full factor base (positive control)",
          "rho": "rho on E0(F_{q^3}), r-adding, no negation"}
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))
for arm in arms:
    pts = curve_points(arm)
    if not pts:
        continue
    q = np.array([p[0] for p in pts]); c = np.array([p[3] for p in pts])
    ax1.scatter(q, c, s=14, color=colors[arm], alpha=0.75, label=None)
    e, a, _ = fits[("primary", arm)]
    s = summary[f"primary:{arm}"]
    qq = np.geomspace(q.min(), q.max(), 50)
    ax1.plot(qq, np.exp(a) * qq ** e, color=colors[arm],
             label=f"{labels[arm]}: q^{e:.3f} [{s['ci95'][0]:.3f}, {s['ci95'][1]:.3f}]")
ax1.set_xscale("log", base=2); ax1.set_yscale("log", base=2)
ax1.set_xlabel("q (prime, base field of the genus-3 curve)")
ax1.set_ylabel("cost, F_q-multiplication equivalents (M + log2(q) I + M_ell)")
ax1.set_title("End-to-end DLP cost per solve (points = curves; rho = mean of targets)")
ax1.legend(fontsize=7.5, loc="upper left")
ax1.grid(True, which="major", alpha=0.3)
bits = [t["bits"] for t in table if "ratio_ic_rho" in t]
rat = [t["ratio_ic_rho"] for t in table if "ratio_ic_rho" in t]
lo = [t["ratio_ci95"][0] for t in table if "ratio_ic_rho" in t]
hi = [t["ratio_ci95"][1] for t in table if "ratio_ic_rho" in t]
rneg = [t["ratio_ic_rho_neg"] for t in table if "ratio_ic_rho" in t]
r9 = [t["ratio_ic_rho_cl9"] for t in table if "ratio_ic_rho" in t]
ax2.errorbar(bits, rat, yerr=[np.array(rat) - lo, np.array(hi) - rat], fmt="o-", color="#1f6fb4",
             capsize=3, label="IC-2LP / rho (primary unit; 95% bootstrap over curves)")
ax2.plot(bits, rneg, "s--", color="#c2410c", ms=4, label="IC-2LP / ideal negation-map rho (rho / sqrt 2)")
ax2.plot(bits, r9, "^:", color="#555555", ms=4, label="IC-2LP / rho with M_ell priced at 9 M")
rm = [t["ratio_ic_rhomodel"] for t in table if "ratio_ic_rho" in t]
rmn = [t["ratio_ic_rhomodel_neg"] for t in table if "ratio_ic_rho" in t]
ax2.plot(bits, rm, "D-", color="#6b21a8", ms=4, label="IC-2LP / rho model: measured cost per step x sqrt(pi ell/2)")
ax2.plot(bits, rmn, "D--", color="#a855f7", ms=4, label="IC-2LP / rho model with negation: x sqrt(pi ell/4)")
ax2.axhline(1, color="black", lw=0.8)
ax2.set_yscale("log", base=2)
ax2.set_xlabel("log2 q")
ax2.set_ylabel("cost ratio IC / rho")
ax2.set_title("Measured ratio; below 1 = transfer + IC cheaper than rho")
ax2.legend(fontsize=7.5)
ax2.grid(True, alpha=0.3)
fig.suptitle(f"Cubic GHS (m=3, r=5) -> genus 3 hyperelliptic: IC vs rho, run {os.path.basename(run)}", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(run, "..", "..", "cost_vs_q.svg"))
fig.savefig(os.path.join(run, "..", "..", "cost_vs_q.png"), dpi=150)

#!/usr/bin/env sage
# Analysis of the factor-base -> D* campaign: correlations, yield model fits,
# and work-per-relation optima. Reads only frozen JSON under raw/.
# Usage: sage analyze.sage > analysis.txt
import json, glob, math, os, sys

def _mean(xs):
    return float(sum(xs)) / len(xs)

def _sd(xs):
    m = _mean(xs)
    return math.sqrt(sum((float(x) - m) ** 2 for x in xs) / len(xs))

def pearson(xs, ys):
    n = len(xs)
    if n < 3: return None
    mx, my = _mean(xs), _mean(ys)
    sx, sy = _sd(xs), _sd(ys)
    if sx == 0 or sy == 0: return None
    return sum((float(x) - mx) * (float(y) - my) for x, y in zip(xs, ys)) / (n * sx * sy)

def load(path):
    with open(path) as f:
        return json.load(f)

# ── m=2 grid ─────────────────────────────────────────────────────────────
grid = {}   # (n, structure) -> list of rows
for p in sorted(glob.glob("raw/m2_n*_d*_s*.json")):
    d = load(p)
    for r in d.get("results", []):
        if "skipped" in r or "d_star_mean" not in r:
            row = {"n": d["field_degree"], "dim": r.get("dim"), "structure": r["structure"], "skipped": r.get("skipped")}
            grid.setdefault((d["field_degree"], r["structure"]), []).append(row)
            continue
        n = d["field_degree"]; dim = r["dim"]; st = r["structure"]
        works = []
        for t in r["targets"]:
            if t["d_star"] is None: continue
            sh = t.get("macaulay_shape", {}).get(str(t["d_star"]))
            if sh: works.append(sh[0] * sh[1])
        frac = r["decomposable_fraction"]
        row = {
            "n": n, "dim": dim, "structure": st, "nvars": 2 * dim,
            "d_mean": r["d_star_mean"], "d_min": r["d_star_min"], "d_max": r["d_star_max"],
            "censored": r["censored_above_dmax"], "unsat": r["unsat_targets"],
            "frac": frac, "work_mean": _mean(works) if works else None,
            "wpr": (_mean(works) / frac) if (works and frac > 0) else None,
            "base_size": 2 ** dim,
            "yield_model": (2 ** (2 * dim)) / float(2 ** n),   # |V|^2 / 2^n counting floor
        }
        grid.setdefault((n, st), []).append(row)

print("== m=2 grid: D*, yield, work-per-relation (mean over seeds) ==")
print("n   st         dim nvars D*mean[min..max] cens frac     W(D*)     W/rel     floor|V|^2/2^n  ratio")
coll = []
for (n, st), rows in sorted(grid.items()):
    dims = sorted(set(r["dim"] for r in rows))
    for dim in dims:
        rs = [r for r in rows if r["dim"] == dim and "d_mean" in r]
        if not rs: continue
        dmm = _mean([r["d_mean"] for r in rs]); dmin = min(r["d_min"] for r in rs); dmax = max(r["d_max"] for r in rs)
        cens = sum(r["censored"] for r in rs); frac = _mean([r["frac"] for r in rs])
        ws = [r["work_mean"] for r in rs if r["work_mean"]]
        wprs = [r["wpr"] for r in rs if r["wpr"]]
        W = _mean(ws) if ws else None; WPR = _mean(wprs) if wprs else None
        fl = 2 ** (2 * dim) / float(2 ** n)
        print("%-3d %-10s %3d %5d %5.2f [%d..%d] %4d %.4f %-9s %-9s %.3f" %
              (n, st, dim, 2 * dim, dmm, dmin, dmax, cens, frac,
               ("%.3g" % W) if W else "-", ("%.3g" % WPR) if WPR else "-", fl))
        if WPR: coll.append((n, st, dim, dmm, frac, WPR))

print("\n== correlations (m=2, non-subfield structures pooled) ==")
flat = [r for rows in grid.values() for r in rows if "d_mean" in r and r["structure"] != "subfield"]
print("pearson(D*mean, dim)      =", pearson([r["dim"] for r in flat], [r["d_mean"] for r in flat]))
print("pearson(D*mean, nvars)    =", pearson([r["nvars"] for r in flat], [r["d_mean"] for r in flat]))
print("pearson(D*mean, n)        =", pearson([r["n"] for r in flat], [r["d_mean"] for r in flat]))
wfl = [r for r in flat if r["wpr"]]
print("pearson(log W/rel, dim)   =", pearson([r["dim"] for r in wfl], [math.log2(r["wpr"]) for r in wfl]))
print("pearson(log frac, dim)    =", pearson([r["dim"] for r in flat], [math.log2(r["frac"] + 1e-12) for r in flat]))
# yield fit: frac vs |V|^2/2^n
print("pearson(log frac, log floor)=", pearson([math.log2(r["frac"] + 1e-12) for r in flat],
                                               [math.log2(r["yield_model"] + 1e-12) for r in flat]))

print("\n== structure contrasts at matched (n, dim) ==")
for (n, st), rows in sorted(grid.items()):
    if st != "subfield": continue
    for r in rows:
        if "skipped" in r: continue
        dim = r["dim"]
        for st2 in ("standard", "random"):
            r2 = next((x for x in grid.get((n, st2), []) if x["dim"] == dim and "d_mean" in x), None)
            if r2 and "d_mean" in r:
                print("n=%d dim=%d: subfield D*=%.2f frac=%.4f  vs %s D*=%.2f frac=%.4f" %
                      (n, dim, r["d_mean"], r["frac"], st2, r2["d_mean"], r2["frac"]))

print("\n== work-per-relation optima (best dim per (n, structure)) ==")
best = {}
for n, st, dim, dmm, frac, WPR in coll:
    k = (n, st)
    if k not in best or WPR < best[k][3]:
        best[k] = (dim, dmm, frac, WPR)
for k in sorted(best):
    dim, dmm, frac, WPR = best[k]
    print("n=%-3d %-10s best dim=%d (D*=%.2f, frac=%.4f, W/rel=%.3g)" % (k[0], k[1], dim, dmm, frac, WPR))

# ── delta sweeps ─────────────────────────────────────────────────────────
for p in sorted(glob.glob("raw/delta_*.json")):
    d = load(p)
    print("\n== %s: %s ==" % (os.path.basename(p), d.get("description", "")))
    base = None
    for r in d["rows"]:
        if r.get("delta") == 0 and "work_per_relation_proxy" in r: base = r["work_per_relation_proxy"]
    for r in d["rows"]:
        if "skipped" in r:
            print("delta=%+d skipped (%s)" % (r["delta"], r["skipped"])); continue
        wpr = r.get("work_per_relation_proxy")
        print("delta=%+d dim=%d nvars=%d D*mean=%s frac=%s W=%s W/rel=%s%s" % (
            r["delta"], r.get("dim", -1), r.get("nvars", -1),
            ("%.2f" % r["d_star_mean"]) if "d_star_mean" in r else "-",
            ("%.4f" % r["decomposable_fraction"]) if "decomposable_fraction" in r else "-",
            ("%.3g" % r["work_proxy"]) if "work_proxy" in r else float("nan"),
            ("%.3g" % wpr) if wpr else "-",
            ("  %.1fx vs delta=0" % (base / wpr)) if (wpr and base) else ""))

# ── m=3 chained and symmetric S4 ────────────────────────────────────────
for p in sorted(glob.glob("raw/m3_*.json")):
    d = load(p)
    print("\n== %s (m=3 chained) ==" % os.path.basename(p))
    for r in d.get("results", []):
        if "skipped" in r or "d_star_mean" not in r:
            print(r.get("structure"), r.get("skipped", "no D*")); continue
        print("n=%d %s dim=%d nvars=%s D*=%.2f [%d..%d] frac=%.4f cens=%d" % (
            d["field_degree"], r["structure"], r["dim"], r.get("nvars"),
            r["d_star_mean"], r["d_star_min"], r["d_star_max"],
            r["decomposable_fraction"], r["censored_above_dmax"]))
for p in sorted(glob.glob("raw/sym4_*.json")):
    d = load(p)
    print("\n== %s (symmetric S4) ==" % os.path.basename(p))
    for r in d.get("results", []):
        if "skipped" in r: print(r.get("structure"), r.get("skipped")); continue
        print("n=%d %s dim=%d nvars=%s D*=%s [%s..%s] point-frac=%.4f bool-unsat=%.4f spurious=%s eqdeg=%s" % (
            d["field_degree"], r["structure"], r["dim"], r.get("nvars"),
            ("%.2f" % r["d_star_mean"]) if "d_star_mean" in r else "-",
            r.get("d_star_min", "-"), r.get("d_star_max", "-"),
            r.get("point_decomposable_fraction", float("nan")),
            r.get("boolean_unsat_fraction", float("nan")),
            r.get("spurious_solutions"), r.get("equation_degree")))

# ── n=131 applicability ──────────────────────────────────────────────────
print("\n== frobwin_n131 (engine writes to stderr; log is the record) ==")
with open("raw/frobwin_n131.log") as f:
    print(f.read().strip())
for p in sorted(glob.glob("raw/stack_n*.json")):
    print("\n== %s ==" % os.path.basename(p))
    print(json.dumps(load(p), indent=1)[:2400])

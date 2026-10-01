#!/usr/bin/env python3
"""Summarise results/<tag>/*.jsonl into one table (Markdown + JSON) and fit
the calls-per-target exponent in |F| over the prime sizes 11..19.

    python3 analyse.py --tag main
"""
import argparse, glob, json, math, os, statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))

def load(tag):
    cells = {}
    manifests = {}
    for f in sorted(glob.glob(os.path.join(HERE, "results", tag, "n*_*.jsonl"))):
        with open(f) as fh:
            for line in fh:
                r = json.loads(line)
                if r["kind"] == "manifest":
                    manifests[r["n"]] = r
                    continue
                key = (r["n"], r["encoding"])
                cells.setdefault(key, []).append(r)
    return cells, manifests

def mean(xs):
    return st.fmean(xs) if xs else float("nan")

def summarise(cells, manifests):
    rows = []
    for (n, enc), rs in sorted(cells.items()):
        m = manifests[n]
        base = rs[0]["base"]
        fsize = m["sub_points"] if enc == "SUB" else m["wt_points"]
        agree = sum(r["agree"] == "ok" for r in rs)
        timeouts = sum(r["agree"] == "timeout" for r in rs)
        mism = sum(r["agree"] == "MISMATCH" for r in rs)
        yes = [r for r in rs if r["label"]]
        no = [r for r in rs if not r["label"]]
        depths = [d for r in rs for d in r["tame_depths"]]
        rows.append({
            "n": n, "encoding": enc, "base": base, "F_points": fsize, "F_x": m["sub_x"] if enc == "SUB" else m["wt_x"],
            "vars": rs[0]["vars"], "aux_vars": rs[0]["aux_vars"], "gens": rs[0]["gens"], "max_gen_degree": rs[0]["max_gen_degree"],
            "targets": len(rs), "decomposable": len(yes), "agree": agree, "timeouts": timeouts, "mismatches": mism,
            "encoding_violations": sum(r["encoding_violations"] for r in rs),
            "calls_mean": mean([r["calls"] for r in rs]), "calls_yes": mean([r["calls"] for r in yes]), "calls_no": mean([r["calls"] for r in no]),
            "tame_mean": mean([r["tame"] for r in rs]), "wild_mean": mean([r["wild"] for r in rs]),
            "budget_mean": mean([r["budget_calls"] for r in rs]), "cap_mean": mean([r["matrix_cap_hits"] for r in rs]),
            "tame_depth_mean": mean(depths), "coords_per_summand": len(rs[0]["tame_depths"]) and None,
            "xor_mean": mean([r["xor_words"] for r in rs]), "xor_per_call": mean([r["xor_words"] / r["calls"] for r in rs]),
            "rows_max": max(r["rows_max"] for r in rs), "cols_max": max(r["cols_max"] for r in rs), "max_f4_degree": max(r["max_f4_degree"] for r in rs),
            "solve_ms_mean": mean([r["solve_us"] / 1000 for r in rs]),
            "exhaustive_scanned": rs[0]["exhaustive_scanned"], "exhaustive_first_hit_mean": mean([r["exhaustive_first_hit"] for r in yes if r["exhaustive_first_hit"]]),
            "exhaustive_us_mean": mean([r["exhaustive_us"] for r in rs]),
            "calls_over_F": mean([r["calls"] for r in rs]) / fsize,
        })
    return rows

def fit(rows, enc, key="calls_mean"):
    pts = [(math.log2(r["F_points"]), math.log2(r[key])) for r in rows if r["encoding"] == enc and r["n"] >= 11 and r["timeouts"] == 0 and r[key] > 0]
    if len(pts) < 2:
        return None
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    xm, ym = mean(xs), mean(ys)
    sxx = sum((x - xm) ** 2 for x in xs)
    slope = sum((x - xm) * (y - ym) for x, y in zip(xs, ys)) / sxx
    return {"encoding": enc, "sizes": [r["n"] for r in rows if r["encoding"] == enc and r["n"] >= 11 and r["timeouts"] == 0], "exponent": slope, "points": pts}

def md_table(rows):
    cols = ["n", "encoding", "base", "F_points", "vars", "gens", "max_gen_degree", "targets", "decomposable", "agree", "timeouts",
            "calls_mean", "calls_no", "tame_mean", "wild_mean", "budget_mean", "tame_depth_mean", "xor_per_call", "solve_ms_mean", "calls_over_F"]
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join("--:" for _ in cols) + "|"]
    for r in rows:
        vals = []
        for c in cols:
            v = r[c]
            if isinstance(v, float):
                vals.append("nan" if math.isnan(v) else (f"{v:.3g}" if abs(v) >= 1000 or abs(v) < 0.01 else f"{v:.2f}"))
            else:
                vals.append(str(v))
        out.append("| " + " | ".join(vals) + " |")
    return "\n".join(out)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", required=True); a = ap.parse_args()
    cells, manifests = load(a.tag)
    rows = summarise(cells, manifests)
    fits = [f for f in (fit(rows, e) for e in ["C", "FC", "QFC", "MONO", "SUB"]) if f]
    xfits = [f for f in (fit(rows, e, "xor_mean") for e in ["C", "FC", "QFC", "MONO", "SUB"]) if f]
    out = {"rows": rows, "fits_calls": fits, "fits_xor": xfits, "manifests": manifests}
    # Round floats so that the file re-derives byte-identically across
    # platforms (the last digit of a fitted slope varies with the summation).
    def rounded(v):
        if isinstance(v, float):
            return round(v, 9)
        if isinstance(v, list):
            return [rounded(x) for x in v]
        if isinstance(v, tuple):
            return [rounded(x) for x in v]
        if isinstance(v, dict):
            return {k: rounded(x) for k, x in v.items()}
        return v
    with open(os.path.join(HERE, "results", a.tag, "summary.json"), "w") as f:
        json.dump(rounded(out), f, indent=1)
    print(md_table(rows))
    print()
    for f in fits:
        print(f"calls exponent {f['encoding']}: {f['exponent']:.3f} over n={f['sizes']}")
    for f in xfits:
        print(f"xor exponent {f['encoding']}: {f['exponent']:.3f} over n={f['sizes']}")

if __name__ == "__main__":
    main()

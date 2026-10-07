#!/usr/bin/env python3
"""Summarise results/*.jsonl into markdown tables (no third-party deps).
usage: scripts/report.py results/baseline.jsonl > results/baseline.md"""
import json, sys, statistics as st
from collections import defaultdict

recs = [json.loads(l) for l in open(sys.argv[1]) if l.strip()]
def fmt_ns(x):
    if x < 1e3: return f"{x:.0f} ns"
    if x < 1e6: return f"{x/1e3:.1f} us"
    if x < 1e9: return f"{x/1e6:.2f} ms"
    return f"{x/1e9:.2f} s"

out = []
def table(title, header, rows):
    out.append(f"\n### {title}\n")
    out.append("| " + " | ".join(header) + " |")
    out.append("|" + "|".join("---" for _ in header) + "|")
    for r in rows: out.append("| " + " | ".join(str(c) for c in r) + " |")

# kernel -> isogeny
k = [r for r in recs if r["group"] == "kernel"]
if k:
    for bits in sorted({r["p_bits"] for r in k}, key=int):
        for task in ["codomain_from_point", "codomain_given_h", "eval_point"]:
            sel = [r for r in k if r["p_bits"] == bits and (r["task"] == task or (task == "codomain_from_point" and r["task"].startswith("codomain_from_point")))]
            algos = []
            for r in sel:
                key = r["algo"] + ("(h given)" if r["task"] == "codomain_given_h" else "(poly+formulas)" if "poly+formulas" in r["task"] else "")
                if key not in algos: algos.append(key)
            ells = sorted({int(r["ell"]) for r in sel})
            rows = []
            for l in ells:
                row = [l]
                for a in algos:
                    m = [r for r in sel if int(r["ell"]) == l and (r["algo"] + ("(h given)" if r["task"] == "codomain_given_h" else "(poly+formulas)" if "poly+formulas" in r["task"] else "")) == a]
                    row.append(fmt_ns(m[0]["median_ns"]) if m else "-")
                rows.append(row)
            if rows: table(f"Kernel->isogeny, {bits}-bit p, task `{task}` (median)", ["l"] + algos, rows)

# find
f = [r for r in recs if r["group"] == "find"]
if f:
    for bits in sorted({r["p_bits"] for r in f}, key=int):
        sel = [r for r in f if r["p_bits"] == bits]
        rows = []
        for l in sorted({int(r["ell"]) for r in sel}):
            d = {r["algo"]: r for r in sel if int(r["ell"]) == l}
            g = lambda a: fmt_ns(d[a]["median_ns"]) if a in d else "-"
            rows.append([l, d.get("divpoly_factor", {}).get("kernels_found", "-"), g("divpoly_factor"), g("elkies_phi_bmss"), g("phi_roots_only"), g("phi_setup")])
        table(f"Kernel finding, {bits}-bit p (median; Phi setup is a single run)", ["l", "#kernels", "divpoly factor", "Phi roots+Elkies+BMSS", "Phi roots only", "Phi_l setup (once)"], rows)

# path
p = [r for r in recs if r["group"] == "path"]
if p:
    by = defaultdict(list)
    for r in p: by[(r["algo"], r["p_bits"])].append(r)
    rows = []
    for (algo, bits), rs in sorted(by.items(), key=lambda kv: (kv[0][0], int(kv[0][1]))):
        ok = sum(1 for r in rs if r["verified"])
        okr = [r for r in rs if r["verified"]]
        med = fmt_ns(st.median(r["ns"] for r in okr)) if okr else "-"
        def mm(key):
            v = [r[key] for r in okr if key in r]
            return f"{st.median(v):.0f}" if v else "-"
        work = " ".join(f"{k}={mm(k)}" for k in ["nodes_expanded", "steps", "nodes", "walk_steps", "bfs_nodes"] if any(k in r for r in okr))
        rows.append([algo, bits, f"{ok}/{len(rs)}", med, mm("path_len"), work])
    table("Path finding (median over verified instances)", ["algorithm", "p bits", "verified", "median time", "path len", "work counters (median)"], rows)
    bad = [r for r in p if not r["verified"]]
    if bad:
        out.append("\nUnverified / failed instances (kept in the raw data):\n")
        for r in bad:
            out.append(f"- {r['algo']} p_bits={r['p_bits']} instance={r.get('instance')} ns={r.get('ns')} note: {r['note']}")

# ---------------------------------------------------------------- V2 groups
def pivot(title, recs, rowkey, colkey, valkey="median_ns", fmt=fmt_ns, rowname=None, sortrow=int):
    rows = sorted({r[rowkey] for r in recs}, key=lambda v: sortrow(v))
    cols = []
    for r in recs:
        if r[colkey] not in cols: cols.append(r[colkey])
    body = []
    for rw in rows:
        line = [rw]
        for c in cols:
            m = [r for r in recs if r[rowkey] == rw and r[colkey] == c]
            line.append(fmt(m[0][valkey]) if m else "-")
        body.append(line)
    table(title, [rowname or rowkey] + cols, body)

k2 = [r for r in recs if r["group"] == "kernel2"]
if k2: pivot("V2 kernel -> isogeny, odd degree, 32-bit p (median): general Velu from P, x-only Velu from x(P), Kohel with h given", k2, "ell", "algo", rowname="l")
k2e = [r for r in recs if r["group"] == "kernel2_even"]
if k2e: pivot("V2 kernel -> isogeny, even / composite cyclic degree n, 32-bit p (median)", k2e, "degree", "algo", rowname="n")
k2m = [r for r in recs if r["group"] == "kernel2_mont"]
if k2m: pivot("V2 Montgomery x-only Velu vs Weierstrass Velu for the same CSIDH kernel (median)", k2m, "ell", "algo", rowname="l")
ch = [r for r in recs if r["group"] == "chain"]
if ch:
    rows = []
    for r in sorted(ch, key=lambda r: (int(r["ell"]), int(r["e"]), r["p_bits"], r["strategy"])):
        rows.append([r["ell"], r["e"], r["p_bits"], r["strategy"], fmt_ns(r["median_ns"]), int(r["l_mults"]), int(r["evals"]), int(r["builds"])])
    table("V2 l^e-kernel chains over F_p^2 (median); identical codomain for all strategies", ["l", "e", "p bits", "strategy", "time", "l-mults", "point evals", "isogeny builds"], rows)
bm = [r for r in recs if r["group"] == "bmss" and r["algo"] not in ("end_to_end_via_phi", "sigma_from_phi", "dual_isogeny")]
if bm: pivot("V2 (E, E~) -> isogeny, 40-bit p (median; sigma supplied where required)", bm, "ell", "algo", rowname="l")
ee = [r for r in recs if r["group"] == "bmss" and r["algo"] in ("sigma_from_phi", "dual_isogeny")]
if ee: pivot("V2 sigma from Phi and dual isogeny (median)", ee, "ell", "algo", rowname="l")
e2e = [r for r in recs if r["group"] == "bmss" and r["algo"] == "end_to_end_via_phi"]
if e2e:
    for r in e2e: r["algo2"] = "e2e_" + r["method"]
    pivot("V2 end-to-end from (E, l) with Phi_l given (median)", e2e, "ell", "algo2", rowname="l")
cs = [r for r in recs if r["group"] == "csidh"]
if cs:
    rows = []
    for r in cs:
        extra = {k: v for k, v in r.items() if k not in ("group", "algo", "verified", "note", "median_ns", "min_ns", "reps")}
        t = r.get("median_ns", r.get("ns"))
        rows.append([r["algo"], json.dumps(extra), fmt_ns(t) if t else "-", "yes" if r["verified"] else "NO"])
    table("V2 CSIDH-style action", ["what", "parameters / counters", "time", "verified"], rows)
en = [r for r in recs if r["group"] == "endo"]
if en:
    rows = [[r["p_bits"], r["instance"], r["height_3"], r["level_3"], fmt_ns(r["median_ns"]), "yes" if r["verified"] else "NO"] for r in en]
    table("V2 Kohel End(E) conductor (median)", ["p bits", "instance", "height of 3-volcano", "level of E", "time", "verified"], rows)
wk = [r for r in recs if r["group"] == "walks"]
if wk:
    by = defaultdict(list)
    for r in wk: by[(r["algo"], r["p_bits"])].append(r)
    rows = []
    for (algo, bits), rs in sorted(by.items(), key=lambda kv: (int(kv[0][1]), kv[0][0])):
        ok = [r for r in rs if r["verified"]]
        med = fmt_ns(st.median(r["ns"] for r in ok)) if ok else "-"
        steps = [r.get("steps", r.get("nodes_expanded")) for r in ok]
        rows.append([algo, bits, f"{len(ok)}/{len(rs)}", med, f"{st.median(steps):.0f}" if steps else "-"])
    table("V2 walks on identical instances (median over verified instances)", ["algorithm", "p bits", "verified", "median time", "median steps / nodes"], rows)
print("\n".join(out))

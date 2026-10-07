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
print("\n".join(out))

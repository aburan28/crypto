#!/usr/bin/env python3
"""Markdown tables for README.md from runs/*.jsonl.  python3 -I summarize.py"""
import json, os, statistics
R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "runs")


def load(name):
    p = os.path.join(R, name)
    return [json.loads(l) for l in open(p)] if os.path.exists(p) else []


XL = {(4, 4, 1): "5 5 5 5", (4, 5, 1): "5 5 5 5", (4, 6, 1): "5 5 5 5", (4, 5, 2): "6 6 6 6",
      (4, 6, 2): "6 6 6 6", (4, 7, 2): "6 6 6 6", (4, 8, 2): "6 6 6 6",
      (4, 9, 3): "≥7 ≥7 ≥7 ≥7", (4, 10, 3): "≥7 ≥7 ≥7 ≥7",
      (3, 9, 3): "6 6 6 6 (other draws)", (3, 13, 4): "6 6 6 6", (3, 10, 5): "7 7 7 7",
      (3, 15, 5): "≥7 ≥7 ≥7 ≥7"}
rows = load("c2_calibration.jsonl") + load("m_main.jsonl") + load("m_main_b.jsonl") + load("m_main_c.jsonl") + load("e_sat.jsonl")
cells = {}
sat = []
for r in rows:
    if "export" not in r:
        continue
    e = r["export"]
    if e.get("sat_ok"):
        sat.append(r)
        continue
    cells.setdefault((e["m"], e["n"], e["ell"]), []).append(r)
print("| m | (n, ℓ) | N | S | XL refutation degree (committed) | F4 degree per draw | GB = [1] | wall s | peak RSS MB | largest F4 matrix |")
print("|--:|---|--:|--:|---|---|---|---|---|---|")
for k in sorted(cells, key=lambda k: (-k[0], k[2], k[1])):
    rs = cells[k]
    m, n, ell = k
    e = rs[0]["export"]
    ok = all(r["exit"] == 0 for r in rs)
    def one_draw(r):
        if r["exit"] == 0:
            return str(r["d_f4"])
        why = "3 GB cap hit" if r["exit"] == -11 else str(r["exit"])
        return f"unresolved: {why} in a degree-{r['d_f4']} round"
    ds = "; ".join(one_draw(r) for r in rs)
    one = "yes" if all(r.get("gb_is_one") for r in rs) else "NO"
    secs = [r["secs"] for r in rs]
    mem = [r["peak_rss_kb"] / 1024 for r in rs]
    mm = max(rs, key=lambda r: r["max_matrix"][0] * r["max_matrix"][1])["max_matrix"]
    print(f"| {m} | ({n}, {ell}) | {e['n_vars']} | {n - m * ell:+d} | {XL.get(k, 'not measured')} | {ds} | {one} | "
          f"{min(secs):.1f}–{max(secs):.1f} | {min(mem):.0f}–{max(mem):.0f} | {mm[0]:,} × {mm[1]:,} |")
print()
print("Instrument correctness (arm E), satisfiable draws:")
for r in sat:
    e = r["export"]
    print(f"- m={e['m']} ({e['n']},{e['ell']}) draw {e['draw']}: exact count {e['solutions']}, "
          f"GB [1]? {r.get('gb_is_one')}, GB length {r.get('gb_len')}, F4 degree {r['d_f4']}, {r['secs']:.0f} s")
print()
print("Sensitivity control (arm D), random systems of the template's shape:")
for r in load("d_random.jsonl"):
    print(f"- {r['label']}: max F4 degree reached {r['d_f4']} ({'finished' if r['exit'] == 0 else r['exit']} at "
          f"{r['secs']:.0f} s), largest matrix {r['max_matrix'][0]:,} × {r['max_matrix'][1]:,}, "
          f"peak {r['peak_rss_kb']/1024:.0f} MB")
print()
print("Secondary (deviation): draws with no exact count, the GB deciding (solving degree of the full GB):")
for r in load("n_nocount.jsonl"):
    e = r["export"]
    st = "finished" if r["exit"] == 0 else ("3 GB cap hit" if r["exit"] == -11 else str(r["exit"]))
    print(f"- m={e['m']} (n={e['n']}, ℓ={e['ell']}) N={e['n_vars']} draw {e['draw']}: {st}, GB [1]? {r.get('gb_is_one')}, "
          f"GB length {r.get('gb_len')}, max F4 degree {r['d_f4']}, {r['secs']:.0f} s wall, {r['peak_rss_kb']/1024:.0f} MB, "
          f"largest matrix {r['max_matrix'][0]:,} × {r['max_matrix'][1]:,}")

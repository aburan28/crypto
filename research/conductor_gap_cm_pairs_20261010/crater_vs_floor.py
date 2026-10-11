#!/usr/bin/env python3
"""Crater vs floor, per arm, from a session's records.jsonl.

The spec lists the crater curve first (make_spec.py). For each arm:
  S_crater = mean S over the crater's workloads, S_floor = mean over all
  floor workloads; ratio = S_floor / S_crater; 95% interval by bootstrap
  over workloads (resample workloads within each level, 4000 draws).
Operation counts only (cost.s); wall time is never read.
Usage: python3 crater_vs_floor.py sessions/<dir> [...]
"""
import json, random, sys, statistics as st
from collections import defaultdict

def load(d):
    spec = json.load(open(f"{d}/spec.json"))
    crater = spec["workloads"]["curves"][0]
    key = lambda c: (c["p"], c["a"], c["b"])
    per = defaultdict(lambda: defaultdict(list))   # arm -> workload_id -> [S]
    level = {}
    for line in open(f"{d}/records.jsonl"):
        r = json.loads(line)
        if r["warmup"] or r["outcome"]["status"] != "verified":
            continue
        w = r["workload"]; wid = w["workload_id"]
        level[wid] = "crater" if key(w["curve_spec"]) == key(crater) else "floor"
        per[r["arm"]][wid].append(r["cost"]["s"])
    return spec, per, level

def boot(xs_c, xs_f, n=4000, seed=1):
    rng = random.Random(seed); out = []
    for _ in range(n):
        c = [rng.choice(xs_c) for _ in xs_c]; f = [rng.choice(xs_f) for _ in xs_f]
        out.append(st.mean(f) / st.mean(c))
    out.sort(); return out[int(0.025 * n)], out[int(0.975 * n)]

for d in sys.argv[1:]:
    spec, per, level = load(d)
    print(f"\n### {d}\n{spec['label']}\n")
    print("| arm | crater workloads | floor workloads | mean S crater | mean S floor | floor/crater | 95% bootstrap | sees gap? |")
    print("|---|---:|---:|---:|---:|---:|---|---|")
    for arm, wl in per.items():
        xc = [st.mean(v) for w, v in wl.items() if level[w] == "crater"]
        xf = [st.mean(v) for w, v in wl.items() if level[w] == "floor"]
        lo, hi = boot(xc, xf)
        ratio = st.mean(xf) / st.mean(xc)
        print(f"| {arm} | {len(xc)} | {len(xf)} | {st.mean(xc):.3f} | {st.mean(xf):.3f} | {ratio:.3f} | [{lo:.3f}, {hi:.3f}] | {'no' if lo <= 1 <= hi else 'interval excludes 1'} |")

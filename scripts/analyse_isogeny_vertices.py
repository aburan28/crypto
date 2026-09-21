#!/usr/bin/env python3
"""Recompute ANALYSIS.md's per-vertex statistics from the pooled CSV.

    python3 scripts/analyse_isogeny_vertices.py experiments/pooled_vertex_data.csv

Reproduces the three tables ANALYSIS.md reports about rho cost across the
isogeny class: the per-bit summary, the within-class coefficient of
variation, and the permutation test for class-structure leakage.  These
were previously computed ad hoc, which is why they could not be checked
when the underlying data changed.
"""
import csv, math, random, statistics, sys
from collections import defaultdict

TRIALS = 2000
SEED = 20260921


def load(path):
    with open(path) as fh:
        return [{k: int(v) for k, v in row.items()} for row in csv.DictReader(fh)]


def summary(rows):
    print("| bits | vertices | succ | median iters | sqrt(pi*n/2) | ratio |")
    print("|------|----------|------|--------------|--------------|-------|")
    for bits in sorted({r["bits"] for r in rows}):
        grp = [r for r in rows if r["bits"] == bits]
        ok = [r for r in grp if r["rho_success"]]
        med = statistics.median([r["rho_iters"] for r in ok]) if ok else float("nan")
        ref = statistics.median([math.sqrt(math.pi * r["n"] / 2) for r in ok]) if ok else float("nan")
        print(f"| {bits} | {len(grp)} | {len(ok)} | {med:.0f} | {ref:.0f} | {med/ref:.2f} |")


def classes(rows, bits):
    """Successful odd-n vertices grouped by their starting curve."""
    by = defaultdict(list)
    for r in rows:
        if r["bits"] == bits and r["rho_success"] and r["n_parity"] == 1:
            by[r["curve_idx"]].append(r["rho_iters"])
    return {k: v for k, v in by.items() if len(v) >= 2}


def cov_table(rows):
    print("\n| bits | classes (>=2 samples) | within-class CoV range |")
    print("|------|-----------------------|------------------------|")
    for bits in sorted({r["bits"] for r in rows}):
        cl = classes(rows, bits)
        covs = [statistics.stdev(v) / statistics.mean(v)
                for v in cl.values() if statistics.mean(v) > 0]
        rng = f"{min(covs):.2f} to {max(covs):.2f}" if covs else "n/a"
        print(f"| {bits} | {len(cl)} | {rng} |")


def permutation(rows):
    rng = random.Random(SEED)
    print(f"\n| bits | classes (>=2) | obs mean within-class var | null 5-95% | p(obs <= null) |")
    print("|------|---------------|---------------------------|------------|----------------|")
    for bits in sorted({r["bits"] for r in rows}):
        cl = classes(rows, bits)
        if len(cl) < 2:
            print(f"| {bits} | {len(cl)} | n/a | n/a | n/a |")
            continue
        sizes = [len(v) for v in cl.values()]
        pool = [x for v in cl.values() for x in v]
        obs = statistics.mean([statistics.pvariance(v) for v in cl.values()])
        null = []
        for _ in range(TRIALS):
            rng.shuffle(pool)
            i, parts = 0, []
            for s in sizes:
                parts.append(pool[i:i + s]); i += s
            null.append(statistics.mean([statistics.pvariance(p) for p in parts]))
        null.sort()
        lo, hi = null[int(0.05 * TRIALS)], null[int(0.95 * TRIALS)]
        p = sum(1 for x in null if obs <= x) / TRIALS
        print(f"| {bits} | {len(cl)} | {obs:,.0f} | [{lo:,.0f}, {hi:,.0f}] | {p:.3f} |")


if __name__ == "__main__":
    data = load(sys.argv[1] if len(sys.argv) > 1 else "experiments/pooled_vertex_data.csv")
    summary(data); cov_table(data); permutation(data)

#!/usr/bin/env python3
"""Recompute ANALYSIS.md's per-vertex statistics from the pooled CSV.

    python3 scripts/analyse_isogeny_vertices.py experiments/pooled_vertex_data.csv

Reproduces every per-vertex table ANALYSIS.md reports about rho cost across
the isogeny class: the per-bit summary (3.2), the odd-`n` conditional
summary (3.2.1), the trace-parity split (1), the conductor stratification
(3.3), the within-class coefficient of variation (3.4) and the permutation
test for class-structure leakage (3.4.1).  These were previously computed ad
hoc, which is why they could not be checked when the underlying data changed.

Round 9 covered only 3.2, 3.4 and 3.4.1, and 3.2.1, 3.3 and the trace-parity
table were left on pre-Velu-fix numbers -- self-contradicting ones, since
3.2.1 reported 24 successes at 14 bits where 3.2 reported 16 vertices in
total.  They are computed here now so the same gap cannot reopen.
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


def odd_summary(rows):
    """3.2.1 -- pooled statistics conditional on n odd."""
    print("\n| bits | n_succ | total odd-n | median iters | sqrt(pi*n/2) "
          "| median ratio | mean ratio | CoV |")
    print("|------|--------|-------------|--------------|--------------"
          "|--------------|------------|-----|")
    for bits in sorted({r["bits"] for r in rows}):
        grp = [r for r in rows if r["bits"] == bits and r["n_parity"] == 1]
        ok = [r for r in grp if r["rho_success"]]
        if not grp:
            print(f"| {bits} | 0 | 0 | n/a | n/a | n/a | n/a | n/a |")
            continue
        if not ok:
            print(f"| {bits} | 0 | {len(grp)} (0%) | n/a | n/a | n/a | n/a | n/a |")
            continue
        med = statistics.median([r["rho_iters"] for r in ok])
        ref = statistics.median([math.sqrt(math.pi * r["n"] / 2) for r in ok])
        ratios = [r["rho_iters"] / math.sqrt(math.pi * r["n"] / 2) for r in ok]
        mean = statistics.mean(ratios)
        cov = statistics.stdev(ratios) / mean if len(ratios) > 1 and mean else float("nan")
        print(f"| {bits} | {len(ok)} | {len(grp)} ({100 * len(ok) / len(grp):.0f}%) "
              f"| {med:.0f} | {ref:.0f} | {med / ref:.3f} | {mean:.2f} | {cov:.2f} |")


def parity_split(rows):
    """Section 1 -- success rate by order parity.

    For p odd, n = p + 1 - t with p + 1 even, so n parity EQUALS t parity:
    trace odd <=> n odd.  ANALYSIS.md's original header paired them the other
    way round ("trace even (n odd)"), contradicting its own derivation two
    lines below it; the columns are labelled by the measured n parity here.
    """
    print("\n| bit-width | n odd (trace odd) | n even (trace even) |")
    print("|-----------|-------------------|---------------------|")
    for bits in sorted({r["bits"] for r in rows}):
        cells = []
        for parity in (1, 0):
            grp = [r for r in rows if r["bits"] == bits and r["n_parity"] == parity]
            ok = [r for r in grp if r["rho_success"]]
            cells.append(f"{len(ok)} / {len(grp)} ({100 * len(ok) / len(grp):.0f}%)"
                         if grp else "0 / 0")
        print(f"| {bits}-bit | {cells[0]} | {cells[1]} |")


def conductor_table(rows):
    """3.3 -- conductor stratification over successful odd-n vertices."""
    for bits in sorted({r["bits"] for r in rows}):
        by = defaultdict(list)
        for r in rows:
            if r["bits"] == bits and r["rho_success"] and r["n_parity"] == 1:
                by[r["conductor"]].append(r["rho_iters"])
        print(f"\n{bits}-bit:")
        if not by:
            print("  (no successful odd-n vertices)")
            continue
        print("| conductor | n | median | mean | CoV |")
        print("|-----------|---|--------|------|-----|")
        for cond in sorted(by):
            v = by[cond]
            mean = statistics.mean(v)
            cov = statistics.stdev(v) / mean if len(v) > 1 and mean else 0.0
            label = "1 (max-order)" if cond == 1 else str(cond)
            print(f"| {label} | {len(v)} | {statistics.median(v):.0f} | {mean:.0f} | {cov:.2f} |")


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
    print("## 3.2 pooled statistics")
    summary(data)
    print("\n## 3.2.1 pooled statistics conditional on n odd")
    odd_summary(data)
    print("\n## 1 success rate by trace parity")
    parity_split(data)
    print("\n## 3.3 conductor stratification (odd n, successes)")
    conductor_table(data)
    print("\n## 3.4 within-class variation")
    cov_table(data)
    print("\n## 3.4.1 permutation test")
    permutation(data)

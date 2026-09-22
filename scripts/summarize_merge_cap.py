#!/usr/bin/env python3
"""Section 11.10's statistic for the merge-level cap: the *paired* difference
in fitted exponent.

`research/notes/index-calculus/RESEARCH_RESIDUAL_WALKS.md` §11.7 leaves one
loose end and names it — the large-prime linear algebra measures above its own
`4/9` because `LargePrimeEliminator` chains merges without bound, and "sieve
implementations cap the merge level for exactly this reason, and that is the
piece this module does not have."  §11.10 registers the experiment; this
summarises it.

**Why paired.**  The fitted exponent varies by about `0.04` between seeds,
which is the same size as the effect a cap produces.  A difference of means
across two arms therefore says very little.  Running the same seed capped and
uncapped makes each seed its own control, so the seed-to-seed spread cancels
and the statistic is the mean of the per-seed differences.

**Why per phase.**  The cap makes the matrix sparser and the eliminator busier
at the same time, and it pays for both in discarded relations.  A single
end-to-end number would hide a phase improving while the method gets worse,
which is the §5 mistake this whole experiment exists to avoid repeating.

Usage:
    python3 scripts/summarize_merge_cap.py BASELINE.json CAPPED.json
"""

from __future__ import annotations

import argparse
import json
import math
import statistics

#: `n` for each prime in the panel, as §11.7's table records it.
LOG2_N = {271: 24.2, 523: 27.1, 1039: 30.1, 2083: 33.1}
SIZES = (271, 523, 1039, 2083)


def large_prime_arm(path: str) -> dict:
    """`--protocol-la` emits three arms per (size, seed) — dense, Wiedemann,
    and Wiedemann with the double-large-prime rule.  Only the third has an
    eliminator, so only it can carry a cap."""
    rows = json.load(open(path))
    arm = [r["gaudry"] for i, r in enumerate(rows) if i % 3 == 2]
    by_size: dict = {}
    for g in arm:
        g["solve_ops"] = g["la_ops"] - g.get("lp_merge_ops", 0)
        by_size.setdefault(g["p"], []).append(g)
    return by_size


def slope(xs, ys) -> float:
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx


def exponents(by_size: dict, key: str, seeds: int) -> list:
    """One fitted exponent per seed, over the four sizes."""
    out = []
    for s in range(seeds):
        xs = [LOG2_N[p] for p in SIZES]
        ys = [math.log2(by_size[p][s][key]) for p in SIZES]
        out.append(slope(xs, ys))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("baseline", help="--merge-cap 0 panel JSON")
    ap.add_argument("capped", help="--merge-cap K panel JSON")
    args = ap.parse_args()

    base, cap = large_prime_arm(args.baseline), large_prime_arm(args.capped)
    seeds = min(len(base[SIZES[0]]), len(cap[SIZES[0]]))
    runs = [x for p in SIZES for x in base[p][:seeds] + cap[p][:seeds]]
    correct = all(x["correct"] for x in runs)
    k = cap[SIZES[0]][0].get("lp_max_merge_level")

    print(f"paired over {seeds} seeds x {len(SIZES)} sizes, cap k = {k}")
    print(f"every run recovered its logarithm: {correct}")
    if not correct:
        print("  -- a row without a verified answer is not a result (AGENTS.md section 2)")
    print()
    header = f"{'quantity':<24}{'uncapped':>10}{'capped':>10}{'paired d':>10}{'+/-':>8}{'worse/better':>14}"
    print(header)
    print("-" * len(header))
    rows = (
        ("end-to-end total", "total_ops"),
        ("linear algebra", "la_ops"),
        ("   of which solve", "solve_ops"),
        ("   of which merge", "lp_merge_ops"),
        ("residuals", "residuals"),
    )
    for label, key in rows:
        a = exponents(base, key, seeds)
        c = exponents(cap, key, seeds)
        d = [y - x for x, y in zip(a, c)]
        sem = statistics.stdev(d) / math.sqrt(seeds) if seeds > 1 else float("nan")
        up = sum(1 for x in d if x > 0)
        print(f"{label:<24}{statistics.fmean(a):10.3f}{statistics.fmean(c):10.3f}"
              f"{statistics.fmean(d):+10.3f}{sem:8.3f}{f'{up}/{seeds - up}':>14}")
    print()
    print("A POSITIVE paired difference means the cap made that quantity grow")
    print("faster in n, which is worse.  4/9 = 0.4444 is the exponent at stake.")
    print()
    print("Discarded relations, the cap's cost, by size:")
    for p in SIZES:
        ab = [x["lp_abandoned"] for x in cap[p][:seeds]]
        rs = [x["residuals"] for x in cap[p][:seeds]]
        share = statistics.fmean(a / r for a, r in zip(ab, rs) if r)
        print(f"  p = {p:<5} abandoned {statistics.fmean(ab):9.1f} of "
              f"{statistics.fmean(rs):9.1f} residuals  ({share * 100:.1f} %)")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Addendum model: the registered statistic under the count, and the sample it needs.

    python3 research/ic_triple_counted_20260923/model_gm.py

Written after the 4-fixture probe and before the confirmation run.  It uses
no probe data: only `model.py`'s count, the triple arm's committed
instructions per unit, and the fixed phases.  Trial counts are simulated as
what the count says they are -- each of `hits` successes a geometric wait with
the count's `p` -- and the cost of a fixture is fixed + table + |F|·trials,
priced at 1,132 instructions a unit.  The two arms draw different bases, so a
fixture's two costs are independent.
"""
import math
import random
import statistics

import model

N, R = model.CELLS["n43a1"]
FIXED, PER_UNIT = 2.03e6, 1132.0


def cost(rng, K, t, hits):
    s = 2 * N * K
    p = model.p_target(N, R, K, t)
    trials = sum(int(math.log(1 - rng.random()) / math.log(1 - p)) + 1 for _ in range(hits))
    return FIXED + PER_UNIT * (t * s * (s + 1) / 2 + s * trials)


def main():
    rng = random.Random(0)
    print("n43a1, counted/triple under the count (200,000 simulated pairs each):")
    for label, hits in (("K+1 hits", 4), ("one extra hit", 5)):
        logs = [math.log(cost(rng, 3, 1, hits) / cost(rng, 2, 1, 3)) for _ in range(200_000)]
        mean_ratio = statistics.fmean(cost(rng, 3, 1, hits) for _ in range(50_000)) / \
            statistics.fmean(cost(rng, 2, 1, 3) for _ in range(50_000))
        sd = statistics.pstdev(logs)
        gm = math.exp(statistics.fmean(logs))
        print(f"  {label:14s}: ratio of means {mean_ratio:.3f}; geometric mean of ratios {gm:.3f}; "
              f"sd of log ratio {sd:.3f}")
        for n in (32, 64, 128):
            half = 1.96 * sd / math.sqrt(n)
            print(f"      {n:>3d} fixtures: expected interval [{gm * math.exp(-half):.3f}, {gm * math.exp(half):.3f}]")


if __name__ == "__main__":
    main()

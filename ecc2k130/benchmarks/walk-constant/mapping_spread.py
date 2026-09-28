"""One mapping's constant against the average over mappings, for truly random
mappings: WALK-CONSTANT.md section 4.

The harnesses' --fixed-mapping rows stray from the mapping average by several
percent at W = 16.  This runs the same first-collision statistic on uniformly
random functions of N points, with no walk structure at all, and shows that
the spread belongs to the statistic: a few walks, each a sizeable fraction of
sqrt(N) long, sample a few neighbourhoods of one mapping's trees.  It shrinks
as the walks get more numerous and shorter.  A campaign's trails (2^33 of
them, 2^28 steps each against sqrt(N) ~ 2^61) are in the limit where it
vanishes, so the mapping average is the constant a campaign pays.
Run: python3 mapping_spread.py N W TRIALS MAPPINGS
"""
import random, sys
from math import sqrt
def rmr(N):
    t, s, k = 1.0, 1.0, 1
    while t > 1e-15 and k < N:
        t *= 1 - k / N; s += t; k += 1
    return s
def stat(f, N, W, T, rng):
    tot = tot2 = 0.0
    for _ in range(T):
        seen = set(); walks = []
        done = False
        for w in range(W):
            x = rng.randrange(N)
            if x in seen: done = True; break
            seen.add(x); walks.append(x)
        while not done:
            for i in range(W):
                y = f[walks[i]]; walks[i] = y
                if y in seen: done = True; break
                seen.add(y)
        v = len(seen); tot += v; tot2 += v * v
    m = tot / T
    return m, sqrt((tot2 / T - m * m) / T)
N, W, T, K = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
e = rmr(N)
rng = random.Random(1)
cs = []
for k in range(K):
    f = [rng.randrange(N) for _ in range(N)]
    m, se = stat(f, N, W, T, rng)
    cs.append(m / e)
    print(f"random mapping {k}: c = {m/e:.4f} +- {se/e:.4f}", flush=True)
mu = sum(cs) / K
print(f"N = {N}, W = {W}: mean {mu:.4f}, sd over mappings {sqrt(sum((c-mu)**2 for c in cs)/(K-1)):.4f}")

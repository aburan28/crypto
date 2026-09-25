"""The sigma walk on classes as pure arithmetic, WALK-CONSTANT.md section 4.

On the classes of <sigma, -1> the sigma walk is x -> (1 + s^j) x in
Z_l^* / <s, -1>, with j = 3 + b and b a branch of the class.  This runs the
same first-collision statistic as the harnesses with b a salted hash mapped
onto the n = 131 branch probabilities, with no curve arithmetic at all: a
third implementation, and the one that shows why n = 19 is degenerate (two of
its eight multipliers fall in one class, since 1 + s^10 = s^-9 (1 + s^9) when
n = 19, and short walks there return to their own classes).
Run: python3 sigma_classes.py N L S W1,W2,... TRIALS
"""
import random, sys
from math import sqrt
def setup(n, l, s):
    orbit = sorted(set([pow(s, k, l) for k in range(n)] + [l - pow(s, k, l) for k in range(n)]))
    rep = {}
    reps = []
    for x in range(1, l):
        if x in rep: continue
        cls = [(x * o) % l for o in orbit]
        r = min(cls)
        for y in cls: rep[y] = r
        reps.append(r)
    return rep, reps
def rmr(N):
    t, s, k = 1.0, 1.0, 1
    while t > 1e-15 and k < N:
        t *= 1 - k / N; s += t; k += 1
    return s
def trials(n, l, s, W, T, hashed=True):
    rep, reps = setup(n, l, s)
    N = len(reps)
    p = [0.1414, 0.1443, 0.1359, 0.1212, 0.1086, 0.1057, 0.1141, 0.1288]
    cdf = []; a = 0
    for x in p: a += x; cdf.append(a)
    mult = [(1 + pow(s, j, l)) % l for j in range(3, 11)]
    tot = 0.0; tot2 = 0.0
    for _ in range(T):
        salt = random.getrandbits(64)
        def br(r):
            u = random.Random(r ^ salt).random()
            for b, c in enumerate(cdf):
                if u < c: return b
            return 7
        seen = set(); walks = []
        done = False
        for w in range(W):
            x = rep[random.randrange(1, l)]
            if x in seen: done = True; break
            seen.add(x); walks.append(x)
        while not done:
            for i in range(W):
                x = walks[i]
                y = rep[(x * mult[br(x)]) % l]
                walks[i] = y
                if y in seen: done = True; break
                seen.add(y)
        v = len(seen); tot += v; tot2 += v * v
    m = tot / T; se = sqrt((tot2 / T - m * m) / T)
    e = rmr(N)
    print(f"n={n} W={W} T={T}: mean {m:.2f} (random {e:.2f}) c = {m/e:.4f} +- {se/e:.4f}", flush=True)
n, l, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
for W in [int(w) for w in sys.argv[4].split(',')]:
    trials(n, l, s, W, int(sys.argv[5]))

import math
PROBES = 3.0
def W2(K, n, r):
    size = 2 * n * K
    lam = size * (size + 1) / 2 / r
    best = None
    for t in range(1, K + 1):
        cov = 1 - ((K - t) / K) ** 2
        w = t * size + (K + PROBES) / (lam * cov)
        if best is None or w < best[0]: best = (w, t)
    return best
def W3(K, n, r):
    size = 2 * n * K
    lam = size * (size + 1) * (size + 2) / 6 / r
    best = None
    for t in range(1, K + 1):
        cov = 1 - ((K - t) / K) ** 3
        w = t * size * (size + 1) / 2 + (K + PROBES) / (lam * cov)
        if best is None or w < best[0]: best = (w, t)
    return best
# merged round 0022: measured IC/rho at each cell's measured best base (orbits)
cells = [("n23a1", 23, 4_196_903, 8, 0.831),
         ("n37a0", 37, 230_603_167, 16, 1.396),
         ("n43a1", 43, 4_644_189_029, 24, 2.207)]
print(f"{'cell':6s} {'K2':>3s} {'W2':>10s} | {'K3*':>3s} {'t3':>3s} {'W3*':>10s} {'W3/W2':>6s} | {'IC/rho v2':>9s} {'pred v3 (f=.65..81)':>20s}")
for name, n, r, K2, ratio in cells:
    w2, t2 = W2(K2, n, r)
    k3, (w3, t3) = min(((K, W3(K, n, r)) for K in range(1, 41)), key=lambda x: x[1][0])
    q = w3 / w2
    lo, hi = ratio * (0.81 * q + 0.19), ratio * (0.65 * q + 0.35)
    print(f"{name:6s} {K2:>3d} {w2:>10.0f} | {k3:>3d} {t3:>3d} {w3:>10.0f} {q:>6.3f} | {ratio:>9.3f} {lo:>9.3f} .. {hi:<8.3f}")
# the exponent, model-only: W* against r on a fixed n, both collectors at their own optimum
print("\nmodel-optimal W against r at n = 43 (rate check):")
for lr in (26, 30, 34, 38, 42):
    r = 2 ** lr
    w2 = min(W2(K, 43, r)[0] for K in range(1, 400))
    w3 = min(W3(K, 43, r)[0] for K in range(1, 400))
    print(f"  r=2^{lr}: W2*=2^{math.log2(w2):.2f}  W3*=2^{math.log2(w3):.2f}  W3/W2={w3/w2:.3f}")

print("\ncalibration: does W2 put its optimum where the tournament measured the best base?")
for name, n, r, K2, ratio in cells:
    k_opt = min(range(1, 200), key=lambda K: W2(K, n, r)[0])
    w_at_meas = W2(K2, n, r)[0]; w_opt = W2(k_opt, n, r)[0]
    print(f"  {name}: W2 argmin K={k_opt:>3d} (measured best {K2:>2d});  W2(measured)/W2(argmin) = {w_at_meas/w_opt:.3f}")
    for lab, fn, K in (("v2", W2, K2), ("v3", W3, min(range(1, 41), key=lambda K: W3(K, n, r)[0]))):
        w, t = fn(K, n, r); size = 2 * n * K
        table = t * size if lab == "v2" else t * size * (size + 1) / 2
        print(f"      {lab} K={K:>2d} t={t:>2d} size={size:>4d}: table {table:>9.0f}  scan {w - table:>9.0f}")

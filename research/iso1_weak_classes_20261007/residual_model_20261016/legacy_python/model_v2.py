"""Model V2 (PROTOCOL_V2.md): exact N_8 and thinning evaluation.

usage:
  python3 model_v2.py exact OMEGA8.csv LEVELS.csv            # test A
  python3 model_v2.py holdout P CENSUS.csv LEVELS.csv [...]  # test B
"""
import csv, math, sys
from collections import defaultdict
from fractions import Fraction

def v2(n):
    n = abs(n); k = 0
    while n % 2 == 0:
        n //= 2; k += 1
    return k

def n8(p, row):
    """Registered N_8(t) for a class_levels.gp row; None if k < 2."""
    t = int(row['t']); Q = p ** 6; k = v2(int(row['f_pi']))
    if k < 2:
        return None
    dk = int(row['D_K']); s = 0 if dk % 2 == 0 else (1 if dk % 8 == 1 else -1)
    deep = (t - (Q + 1)) % 32 == 0 or (t + Q + 1) % 32 == 0
    L = 5 * (2 - s) * 2 ** (k - 2) - 6 * (1 - s) if deep else (2 - s) * 2 ** (k - 2)
    return Fraction(row['H_odd']).limit_denominator(6) * L

def exact(om, lev):
    rows = list(csv.DictReader(open(om))); p = int(rows[0]['p'])
    assert rows[0]['e'] == '3' and all(r['in_omega'] == '1' for r in rows)
    N = defaultdict(int)
    for r in rows:
        N[abs(int(r['trace']))] += int(r['count'])
    bad = n = 0
    for r in csv.DictReader(open(lev)):
        t = int(r['t'])
        if t % p == 0:
            continue
        pred = n8(p, r) or 0
        n += 1
        if pred != N.get(t, 0):
            bad += 1
            if bad <= 5:
                print('  mismatch', t, pred, N.get(t, 0))
    listed = {int(r['t']) for r in csv.DictReader(open(lev))}
    nonord = sorted(t for t in N if t % p == 0)
    extra = set(N) - listed - set(nonord)
    print(f"p={p}: {n} ordinary rows, {bad} mismatches, {len(extra)} unexpected traces, "
          f"nonordinary traces with 8th-power lambda {[(t, N[t]) for t in nonord]} -> {'PASS' if bad == 0 and not extra else 'FAIL'}")

def auc(scores, labels):
    pairs = sorted(zip(scores, labels)); pos = sum(labels); neg = len(labels) - pos
    rs, i = 0.0, 0
    while i < len(pairs):
        j = i
        while j < len(pairs) and pairs[j][0] == pairs[i][0]:
            j += 1
        rs += (i + j + 1) / 2 * sum(l for _, l in pairs[i:j]); i = j
    return (rs - pos * (pos + 1) / 2) / (pos * neg)

def holdout(p, census, lev):
    q = p * p; r = 8 / (q - 1)
    W = {int(x['trace']): int(x['weak_representatives'])
         for x in csv.DictReader(open(census)) if x['trace_status'] == 'ordinary'}
    pz, z, chi, tot_pred = [], [], [], 0.0
    for x in csv.DictReader(open(lev)):
        t = int(x['t']); N = n8(p, x)
        if N is None:
            continue
        n = float(N) / 12; tot_pred += r * float(N)
        pz.append((1 - r) ** n); z.append(int(W[t] == 0))
        if n * r >= 5:
            chi.append((W[t] / 12 - n * r) ** 2 / (n * r * (1 - r)))
    exp0 = sum(pz); sd = math.sqrt(sum(a * (1 - a) for a in pz)); obs = sum(z)
    lo, hi = exp0 - 1.96 * sd, exp0 + 1.96 * sd
    tot_obs = sum(w for t, w in W.items() if t > 0)
    print(f"p={p} rows={len(z)} zeros obs={obs} exp={exp0:.1f} [{lo:.1f},{hi:.1f}] "
          f"{'PASS' if lo <= obs <= hi else 'FAIL'} dispersion={sum(chi)/len(chi):.3f} (n={len(chi)}) "
          f"auc={auc(pz, z):.4f} weak_total obs={tot_obs} pred={tot_pred:.0f}")
    pts = sorted((a, b) for a, b in zip(pz, z) if a > 1e-6); m = len(pts)
    for k in range(5):
        ch = pts[k * m // 5:(k + 1) * m // 5]
        print(f"   P(zero) in [{ch[0][0]:.2g},{ch[-1][0]:.2g}]: exp={sum(a for a,_ in ch):.1f} obs={sum(b for _,b in ch)}")
    print(f"   P(zero) <= 1e-6: obs={sum(b for a,b in zip(pz,z) if a <= 1e-6)}")

if __name__ == '__main__':
    if sys.argv[1] == 'exact':
        exact(sys.argv[2], sys.argv[3])
    else:
        a = sys.argv[2:]
        for i in range(0, len(a), 3):
            holdout(int(a[i]), a[i + 1], a[i + 2])

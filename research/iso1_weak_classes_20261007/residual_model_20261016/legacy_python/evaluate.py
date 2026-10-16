"""Evaluate the registered class-size model (PROTOCOL.md) on census CSVs.

usage: python3 evaluate.py P CENSUS.csv PREDICTORS.csv [P CENSUS PRED ...]
"""
import csv, math, sys

def auc(scores, labels):
    pairs = sorted(zip(scores, labels))
    pos = sum(labels); neg = len(labels) - pos
    rank_sum, i = 0.0, 0
    while i < len(pairs):
        j = i
        while j < len(pairs) and pairs[j][0] == pairs[i][0]:
            j += 1
        avg = (i + j + 1) / 2
        rank_sum += avg * sum(l for _, l in pairs[i:j])
        i = j
    return (rank_sum - pos * (pos + 1) / 2) / (pos * neg)

def run(p, census, pred):
    q = p * p
    W = {int(r['trace']): int(r['weak_representatives'])
         for r in csv.DictReader(open(census)) if r['trace_status'] == 'ordinary'}
    rows = [r for r in csv.DictReader(open(pred)) if float(r['H_sixteenth']) >= 0]
    R = q * q + q
    out = {}
    for name in ('H_sixteenth', 'H_quarter', 'h_K'):
        X = [float(r[name]) for r in rows]
        S = sum(X)
        mu = [R / 12 * x / S for x in X]
        e = [math.exp(-m) for m in mu]
        z = [W[int(r['t'])] == 0 for r in rows]
        exp0 = sum(e); sd = math.sqrt(sum(a * (1 - a) for a in e))
        obs0 = sum(z)
        big = [(W[int(r['t'])] / 12, m) for r, m in zip(rows, mu) if m >= 5]
        disp = sum((n - m) ** 2 / m for n, m in big) / len(big)
        out[name] = dict(obs=obs0, exp=exp0, lo=exp0 - 1.96 * sd, hi=exp0 + 1.96 * sd,
                         disp=disp, auc=auc(e, [int(v) for v in z]), e=e, z=z)
    o = out['H_sixteenth']
    print(f"p={p} rows={len(rows)} zero_obs={o['obs']} | " + " | ".join(
        f"{k}: exp={v['exp']:.1f} [{v['lo']:.1f},{v['hi']:.1f}] {'PASS' if v['lo'] <= v['obs'] <= v['hi'] else 'FAIL'} disp={v['disp']:.3f} auc={v['auc']:.4f}"
        for k, v in out.items()))
    # calibration deciles for the primary predictor
    pts = sorted((a, b) for a, b in zip(o['e'], o['z']) if a > 1e-6)
    n = len(pts); bins = []
    for k in range(10):
        chunk = pts[k * n // 10:(k + 1) * n // 10]
        bins.append(f"[{chunk[0][0]:.2g},{chunk[-1][0]:.2g}] exp={sum(a for a,_ in chunk):.1f} obs={sum(b for _,b in chunk)}")
    tail = sum(b for a, b in zip(o['e'], o['z']) if a <= 1e-6)
    print(f"  calibration (rows with p0>1e-6: {n}; zeros among the rest: {tail})")
    for b in bins:
        print("   ", b)

args = sys.argv[1:]
for i in range(0, len(args), 3):
    run(int(args[i]), args[i + 1], args[i + 2])

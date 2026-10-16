"""Forward predictions of model V2 from class data alone (no census input).

usage: python3 predict_v2.py P LEVELS.csv [P LEVELS.csv ...]
Prints the predicted number of ordinary positive high-depth rows with no weak
representative, its 95% interval, and the predicted total weak mass.
"""
import csv, math, sys
from model_v2 import n8

a = sys.argv[1:]
for i in range(0, len(a), 2):
    p = int(a[i]); q = p * p; r = 8 / (q - 1)
    pz, tot, rows = [], 0.0, 0
    for x in csv.DictReader(open(a[i + 1])):
        N = n8(p, x)
        if N is None:
            continue
        rows += 1; tot += r * float(N); pz.append((1 - r) ** (float(N) / 12))
    e = sum(pz); sd = math.sqrt(sum(z * (1 - z) for z in pz))
    print(f"p={p} high_depth_rows={rows} predicted_zero_rows={e:.1f} interval95=[{e-1.96*sd:.1f},{e+1.96*sd:.1f}] predicted_weak_total={tot:.0f} exact_weak_total={q*q+q}")

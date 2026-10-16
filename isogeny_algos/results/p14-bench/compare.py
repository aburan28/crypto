"""Per-record comparison of two benchmark arms: min over each arm's repetitions, geometric
mean of base/new per group, and the largest changes among records over 5 ms.
Usage (from isogeny_algos/): python3 results/p14-bench/compare.py results/p14-bench \
    main round5-pre-bmss sea phi bmss find csidh big kernel"""
import json, sys, math, glob
D, BASE, NEW = sys.argv[1:4]; groups = sys.argv[4:]
def tkey(r):
    for k in ("median_ns", "ns", "min_ns"):
        if isinstance(r.get(k), (int, float)): return k
    for k, v in r.items():
        if k.endswith("_ns") and isinstance(v, (int, float)): return k
def load(g, v, rep): return [json.loads(l) for l in open(f"{D}/{g}.{v}.{rep}.jsonl")]
allrows = []
for g in groups:
    B1, B2, N1, N2 = load(g, BASE, 1), load(g, BASE, 2), load(g, NEW, 1), load(g, NEW, 2)
    assert len(B1)==len(B2)==len(N1)==len(N2), g
    unver = sum(1 for r in N1+N2 if r.get("verified") is not True)
    rows = []
    for b1, b2, n1, n2 in zip(B1, B2, N1, N2):
        k = tkey(b1)
        if not k: continue
        b = min(b1[k], b2[k]); n = min(n1[k], n2[k])
        if not b or not n: continue
        # within-arm spread: how much the two base reps disagree
        spread = max(b1[k], b2[k]) / min(b1[k], b2[k])
        key = " ".join(f"{kk}={vv}" for kk, vv in b1.items() if isinstance(vv, str) and kk not in ("note","group"))
        rows.append((b / n, spread, b, n, key))
    gm = math.exp(sum(math.log(r[0]) for r in rows) / len(rows))
    big = [r for r in rows if r[2] > 5e6]  # records over 5 ms (less noise)
    gmb = math.exp(sum(math.log(r[0]) for r in big) / len(big)) if big else float("nan")
    print(f"== {g}: {len(rows)} records, new unverified={unver}, geomean {gm:.2f}x; records >5ms: {len(big)}, geomean {gmb:.2f}x")
    for r in sorted(big, key=lambda r: -r[0])[:6]:
        print(f"   {r[0]:5.2f}x (base spread {r[1]:.2f})  {r[2]/1e6:9.2f} -> {r[3]/1e6:9.2f} ms  {r[4]}")
    for r in sorted(big, key=lambda r: r[0])[:3]:
        print(f"   low {r[0]:5.2f}x (base spread {r[1]:.2f})  {r[2]/1e6:9.2f} -> {r[3]/1e6:9.2f} ms  {r[4]}")

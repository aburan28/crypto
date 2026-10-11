import json
LFD_LIBRARY = True
load("lfd_fast.sage")
rows = []
for n in [6, 7, 8, 9, 10]:
    for s in range(40):
        B, F, meta = descend_S3(n, n//2, 31000*n + s, "random", None)
        K, Vb, x3, b = LAST["K"], LAST["Vb"], LAST["x3"], LAST["b"]
        X = B.gens(); npr = n//2
        trV = [int(v.trace()) for v in Vb]
        if all(t == 0 for t in trV): continue
        Lh = sum(B(trV[j])*(X[j]+X[npr+j]) for j in range(npr))
        span = rref_polys(F, B)
        oks = [c for c in (0, 1) if len(rref_polys(span + [Lh + c], B)) == len(span)]
        ok = oks[0] if len(oks) == 1 else (None if not oks else "both")
        sb = b.sqrt()
        rows.append(dict(n=n, c=ok, trx3=int(x3.trace()), trb=int(b.trace()),
                         trsb_x3=int((sb/x3).trace()) if x3 != 0 else None,
                         trb_x3sq=int((b/x3**2).trace()) if x3 != 0 else None, tr1=int(K(1).trace())))
for r in rows:
    if r["c"] != r["trb_x3sq"]: print("MISMATCH", r, "one_in_span", len(oks))
from collections import Counter
print("found:", Counter(r["c"] is not None for r in rows))
for feat in ["trx3", "trb", "trsb_x3", "trb_x3sq", "tr1"]:
    agree = sum(1 for r in rows if r["c"] is not None and r[feat] is not None and r["c"] == r[feat])
    print(feat, agree, "/", sum(1 for r in rows if r["c"] is not None))
for combo in [("trx3","trb_x3sq"), ("trx3","tr1"), ("trx3","trb_x3sq","tr1"), ("trx3","trsb_x3")]:
    agree = sum(1 for r in rows if r["c"] is not None and None not in [r[f] for f in combo] and r["c"] == sum(r[f] for f in combo) % 2)
    print(combo, agree)

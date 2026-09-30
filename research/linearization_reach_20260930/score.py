"""Score results/results.json against PROTOCOL.md P1-P4 (thresholds as frozen)."""
import json
from collections import defaultdict

res = json.load(open("results/results.json"))
out = {}
p1 = []
p2 = []
for c in res["rank_defect"]:
    h = {int(k): v for k, v in c["hist"].items()}
    tot = sum(h.values())
    n, l, b, s = c["n"], c["l"], c["b"], c["lb_plus_l"]
    if s <= n - 4:
        p1.append((n, l, b, h.get(b, 0) / tot))
    if s >= n + 2:
        p2.append((n, l, b, sum(v for k, v in h.items() if k > b) / tot))
out["P1"] = dict(cells=len(p1), min_frac_defect_eq_b=min(f for *_, f in p1),
                 pass_=all(f >= 0.9 for *_, f in p1))
out["P2"] = dict(cells=len(p2), min_frac_defect_gt_b=min(f for *_, f in p2),
                 pass_=all(f >= 0.99 for *_, f in p2))
dev = [(r["n"], r["l"], r["b"], round(r["log2p"] - r["law"], 2)) for r in res["oracle"] if r["success"] >= 8]
groups = defaultdict(list)
for r in res["oracle"]:
    if r["success"] >= 8:
        groups[(r["n"], r["l"])].append((r["b"], r["log2p"]))
slopes = {}
for k, pts in groups.items():
    if len(pts) >= 2:
        mb = sum(b for b, _ in pts) / len(pts)
        my = sum(y for _, y in pts) / len(pts)
        slopes[f"{k}"] = round(sum((b - mb) * (y - my) for b, y in pts) / sum((b - mb) ** 2 for b, _ in pts), 3)
out["P3"] = dict(deviations=dev, slopes=slopes,
                 pass_=all(-2 <= d <= 2 for *_, d in dev) and all(0.6 <= s <= 1.4 for s in slopes.values()))
out["P4"] = dict(algebra_mismatch=sum(r["algebra_mismatch"] for r in res["oracle"]),
                 curve_reject=sum(r["curve_reject"] for r in res["oracle"]),
                 family_too_big=sum(r["family_too_big"] for r in res["oracle"]))
out["P4"]["pass_"] = out["P4"]["algebra_mismatch"] == 0
json.dump(out, open("results/score.json", "w"), indent=1)
print(json.dumps(out, indent=1))

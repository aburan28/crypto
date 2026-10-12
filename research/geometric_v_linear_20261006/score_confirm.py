"""Score results/confirm.json against PROTOCOL.md G1-G3 (thresholds as frozen)."""
import json

cells = json.load(open("results/confirm.json"))["cells"]


def mean_defect(c):
    tot = sum(c["rank_hist"].values())
    return sum((c["unknowns"] - int(k)) * v for k, v in c["rank_hist"].items()) / tot


geo_in = [c for c in cells if c["kind"] == "geometric" and 3 * c["l"] - 1 <= c["n"]]
g1_cells = [c for c in geo_in if c["success"] >= 8]
dev = [(c["n"], c["l"], round(c["log2p"] - c["law"], 2)) for c in g1_cells]
slopes = {}
for n in sorted({c["n"] for c in g1_cells}):
    pts = [(c["l"], c["log2p"]) for c in g1_cells if c["n"] == n]
    if len(pts) >= 3:
        ml = sum(x for x, _ in pts) / len(pts)
        my = sum(y for _, y in pts) / len(pts)
        slopes[n] = round(sum((x - ml) * (y - my) for x, y in pts) / sum((x - ml) ** 2 for x, _ in pts), 3)
g1 = bool(dev) and all(abs(d) <= 1.5 for *_, d in dev) and all(1.6 <= s <= 2.4 for s in slopes.values())

geo_def = [(c["n"], c["l"], round(mean_defect(c), 2)) for c in geo_in]
beyond = [c for c in cells if c["unknowns"] > c["n"]]
excess = [(c["kind"], c["n"], c["l"], c["unknowns"] - c["n"], round(mean_defect(c) - (c["unknowns"] - c["n"]), 2))
          for c in beyond]
g2 = all(d <= 2.0 for *_, d in geo_def) and bool(excess) and all(0 <= e <= 1.5 for *_, e in excess)

g3_cells = [(c["n"], c["l"], c["planted_recovered"], c["planted"]) for c in geo_in]
g3 = all(r >= p - 1 for *_, r, p in g3_cells)

score = dict(G1=dict(cells=len(g1_cells), deviations=dev, slopes=slopes, pass_=g1),
             G2=dict(geometric_mean_defect_in_reach=geo_def, beyond_reach_excess_over_N_minus_n=excess, pass_=g2),
             G3=dict(planted=g3_cells, pass_=g3))
json.dump(score, open("results/score_confirm.json", "w"), indent=1)
print(json.dumps(score, indent=1))

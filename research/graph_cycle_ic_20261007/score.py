"""Score results/n19.json and results/n23.json against PROTOCOL.md Y1-Y5."""
import json
import statistics as st

rows = json.load(open("results/n19.json")) + json.load(open("results/n23.json"))
cells = {}
for r in rows:
    for k, v in r.items():
        if k.startswith("l"):
            l, mode = k[1:].split("_")
            cells.setdefault((r["n"], int(l)), {"walk": [], "fresh": [], "rhoS": []})[mode].append(v)
    for (n, l), c in cells.items():
        if n == r["n"] and len(c["rhoS"]) < len(c["walk"]):
            c["rhoS"].append(r["S_rho"])
y1 = all(r["rho_ok"] for r in rows) and all(v["k_ok"] for c in cells.values() for m in ("walk", "fresh") for v in c[m] if v["solved_by"])
table, y2, y3, y4, y5 = [], True, True, True, True
for (n, l), c in sorted(cells.items()):
    rho_first = sum(v["solved_by"] == "rho" for v in c["walk"]) / len(c["walk"])
    fresh = [v for v in c["fresh"] if v["solved_by"] == "cycle"]
    epn = st.median(v["edges_per_node"] for v in fresh)
    import math
    dev = st.median(math.log2(v["tests"]) - (n - l) for v in fresh)
    s_fresh = st.median(v["S_tests"] for v in fresh)
    s_rho = st.median(c["rhoS"])
    y2 &= rho_first >= 0.9
    y3 &= 0.2 <= epn <= 2
    y4 &= -2 <= dev <= 3
    y5 &= s_fresh > s_rho
    table.append(dict(n=n, l=l, instances=len(c["walk"]), walk_rho_first=rho_first,
                      fresh_solved=f"{len(fresh)}/{len(c['fresh'])}", median_edges_per_node=round(epn, 3),
                      median_log2tests_minus_n_minus_l=round(dev, 2), median_S_tests_fresh=round(s_fresh, 2),
                      median_S_rho=round(s_rho, 3)))
score = dict(Y1=y1, Y2=y2, Y3=y3, Y4=y4, Y5=y5, cells=table)
json.dump(score, open("results/score.json", "w"), indent=1)
print(json.dumps(score, indent=1))

"""Score results/closures.json against PROTOCOL.md C1, C2 and the decisive outcome (as frozen)."""
import json

res = json.load(open("results/closures.json"))
ids = res["identities"]
c1 = dict(rows=[{k: r[k] for k in r} for r in ids],
          pass_=all(r[f"I{i}_ok"] == r[f"I{i}_total"] > 0 for r in ids for i in (1, 2, 3, 4))
          and sorted(r["n"] for r in ids) == [7, 11, 13])

cells = {}
for e in res["energy"]:
    c = cells.setdefault((e["n"], e["l"]), {})
    a = c.setdefault(e["kind"], dict(nontrivial=0, expected=0.0, seeds=0))
    a["nontrivial"] += e["nontrivial"]
    a["expected"] += e["expected_nontrivial"]
    a["seeds"] += 1
rows = []
for (n, l), c in sorted(cells.items()):
    g, rs, rv = c["geometric"], c["randomset"], c["randomV"]
    rows.append(dict(n=n, l=l, seeds=[g["seeds"], rs["seeds"], rv["seeds"]],
                     R=round(g["nontrivial"] / rs["nontrivial"], 4),
                     R_randomV=round(g["nontrivial"] / rv["nontrivial"], 4),
                     geo_vs_expect=round(g["nontrivial"] / g["expected"], 4),
                     randomset_vs_expect=round(rs["nontrivial"] / rs["expected"], 4)))
complete = len(rows) == 7 and all(r["seeds"] == [3, 3, 3] for r in rows)
c2 = dict(cells=rows, complete=complete,
          pass_=complete and all(0.67 <= r["R"] <= 1.5 for r in rows))
Rs = [r["R"] for r in rows]
decisive = complete and all(x >= 1.5 for x in Rs)
score = dict(C1=c1, C2=c2, decisive_idea_survives=decisive)
json.dump(score, open("results/score.json", "w"), indent=1)
print(json.dumps(score, indent=1))

"""Score results/confirm.json against PROTOCOL.md T1 and T2 (thresholds as frozen)."""
import json

res = json.load(open("results/confirm.json"))
viol_a = []
for r in res["part_a"]:
    # a violation: XL resolved at a LOWER t-degree than the count allows
    if r["measured_min_dT"] is not None and (r["count_min_dT"] is None or r["measured_min_dT"] < r["count_min_dT"]):
        viol_a.append(r)
viol_b = [r for r in res["part_b"] if r["measured_bmax"] > r["count_bmax"]]
n_viol = len(viol_a) + len(viol_b)
max_excess = max([r["measured_bmax"] - r["count_bmax"] for r in viol_b] +
                 [r["count_min_dT"] - r["measured_min_dT"] for r in viol_a if r["count_min_dT"] is not None] + [0])
agree_a = sum(1 for r in res["part_a"] if r["measured_min_dT"] == r["count_min_dT"])
over_a = sum(1 for r in res["part_a"] if r["measured_min_dT"] != r["count_min_dT"] and r not in viol_a)
incons = sum(r["inconsistent"] for r in res["raw"])
score = dict(
    T1=dict(cells_a=len(res["part_a"]), cells_b=len(res["part_b"]), violations=n_viol,
            max_excess=max_excess, part_a_exact=agree_a, part_a_count_optimistic_or_censored=over_a,
            part_b_gap=[(r["n"], r["l"], r["measured_bmax"], r["count_bmax"]) for r in res["part_b"]],
            pass_=n_viol <= 1 and max_excess <= 1),
    T2=dict(inconsistent_planted=incons, pass_=incons == 0),
)
json.dump(score, open("results/score_confirm.json", "w"), indent=1)
print(json.dumps(score, indent=1))

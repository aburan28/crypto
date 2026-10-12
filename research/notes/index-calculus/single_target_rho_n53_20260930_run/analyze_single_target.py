#!/usr/bin/env python3
"""Apply the registered rules of RESEARCH_SINGLE_TARGET_STRONG_RHO_N53_20260930.md."""
import json
import os
import random
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
J = lambda n: json.load(open(os.path.join(HERE, n)))
tim, seeds, cg = J("timing.json"), J("seeds.json"), J("callgrind.json")
B = tim["blocks"]
out, lines = {"gates": {}}, []
say = lambda *a: lines.append(" ".join(str(x) for x in a))
pct = lambda v, q: sorted(v)[min(len(v) - 1, int(q * len(v)))]

med = lambda key, arm, idx=None: st.median(B[i][arm][key] for i in (idx if idx is not None else range(len(B))))
W = {a: med("wall_s", a) for a in "DPS"}
C = {a: med("cpu_s", a) for a in "DPS"}
R2 = W["D"] / W["S"]
rng = random.Random(20260930)
boot = []
for _ in range(10000):
    idx = rng.choices(range(len(B)), k=len(B))
    boot.append(med("wall_s", "D", idx) / med("wall_s", "S", idx))
lo, hi = pct(boot, 0.025), pct(boot, 0.975)
out["timing"] = {"W_median_s": W, "C_median_s": C, "spread_maxmin_wall": {a: max(B[i][a]["wall_s"] for i in range(len(B))) / min(B[i][a]["wall_s"] for i in range(len(B))) - 1 for a in "DPS"},
                 "R2_W_D_over_W_S": R2, "R2_ci95": [lo, hi], "W_D_over_W_P": W["D"] / W["P"],
                 "C_D_over_C_S": C["D"] / C["S"], "C_D_over_C_P": C["D"] / C["P"],
                 "maxrss_kb_median": {a: med("maxrss_kb", a) for a in "DPS"}}
all_runs = [B[i][a] for i in range(len(B)) for a in "DPS"] + seeds["runs"] + cg["runs"]
out["gates"]["G1_all_complete_and_correct"] = all(r["gate"] and r["rc"] == 0 for r in all_runs)
if R2 and hi < 1:
    verdict = "SURVIVES (upper bound of W_D/W_S < 1)"
elif lo >= 1:
    verdict = "DIES (lower bound of W_D/W_S >= 1)"
else:
    verdict = "UNRESOLVED"

runs = {a: [r for r in seeds["runs"] if r["arm"] == a] for a in "PS"}
sw = {a: [r["wall_s"] for r in runs[a]] for a in "PS"}
Q = st.median(sw["P"]) / st.median(sw["S"])
bq = []
for _ in range(10000):
    p = rng.choices(sw["P"], k=len(sw["P"]))
    s = rng.choices(sw["S"], k=len(sw["S"]))
    bq.append(st.median(p) / st.median(s))
qlo, qhi = pct(bq, 0.025), pct(bq, 0.975)
ratio_steps = {a: st.median(r["extra"]["walk_steps"] / r["extra"]["ideal_steps"] for r in runs[a]) for a in "PS"}
out["gates"]["G2_walk_steps_band_0.7_3"] = all(0.7 <= v <= 3 for v in ratio_steps.values())
out["seeds"] = {"n": {a: len(runs[a]) for a in "PS"}, "median_wall_s": {a: st.median(sw[a]) for a in "PS"},
                "min_wall_s": {a: min(sw[a]) for a in "PS"}, "max_wall_s": {a: max(sw[a]) for a in "PS"},
                "Q_medianP_over_medianS": Q, "Q_ci95": [qlo, qhi], "median_steps_over_ideal": ratio_steps,
                "panel_seed_wall_s_P": runs["P"][0]["wall_s"], "panel_seed_rank_in_P_walls": sorted(sw["P"]).index(runs["P"][0]["wall_s"]) + 1,
                "median_walk_steps": {a: st.median(r["extra"]["walk_steps"] for r in runs[a]) for a in "PS"}}
if not out["gates"]["G2_walk_steps_band_0.7_3"]:
    strawman = "WITHHELD (G2 failed)"
else:
    strawman = "P IS A STRAWMAN (lower bound of Q >= 1.5)" if qlo >= 1.5 else "not established (lower bound of Q < 1.5)"
ir = {r["arm"]: r["Ir"] for r in cg["runs"]}
out["callgrind"] = {"Ir": ir, "Ir_D_over_Ir_S": ir["D"] / ir["S"], "Ir_D_over_Ir_P": ir["D"] / ir["P"], "Ir_P_over_Ir_S": ir["P"] / ir["S"]}
out["G3_panel_headline_reproduces_here"] = W["D"] / W["P"] < 1.0
out["verdict_wall_win_vs_strong_rho"] = verdict if out["gates"]["G1_all_complete_and_correct"] else "UNDETERMINED (G1 failed)"
out["strawman_test"] = strawman
for k, v in out.items():
    say(f"== {k}"); say(json.dumps(v, indent=1))
open(os.path.join(HERE, "analysis.txt"), "w").write("\n".join(lines) + "\n")
json.dump(out, open(os.path.join(HERE, "analysis.json"), "w"), indent=2)
print("\n".join(lines))

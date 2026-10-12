#!/usr/bin/env python3
"""Post-registration analysis for Amendments A1/A2 of RESEARCH_CACHEGRIND_N41_20260929.md.
Reads analysis.json (registered stages), explore_replicate.json, explore_interference.json,
explore_thp.json, explore_thp_a2.json. Writes analysis_explore.txt/json. Exploratory: no output
here can change the registered verdict.
"""
import json
import os
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
J = lambda n: json.load(open(os.path.join(HERE, n)))
reg, rep, itf, thp, a2, cal = J("analysis.json"), J("explore_replicate.json"), J("explore_interference.json"), \
    J("explore_thp.json"), J("explore_thp_a2.json"), J("calibration_summary.json")
nat = J("native.json")
ir = {a: reg["arms"][a]["Ir"] for a in ("ic", "rho")}
out, lines = {}, []
say = lambda *a: lines.append(" ".join(str(x) for x in a))
q = lambda v: st.quantiles(v, n=4, method="inclusive")

# ---- E3 pooled timing and post-hoc F ----
cpu = {a: [x[a]["cpu_s"] for x in nat["reps"]] + [x[a]["cpu_s"] for x in rep["reps"]] for a in ("ic", "rho")}
pooled = {a: st.median(v) for a, v in cpu.items()}
out["E3"] = {"n": {a: len(v) for a, v in cpu.items()}, "median_cpu_s": pooled,
             "spread_maxmin": {a: max(v) / min(v) - 1 for a, v in cpu.items()},
             "iqr_over_median": {a: (q(v)[2] - q(v)[0]) / st.median(v) for a, v in cpu.items()},
             "min_cpu_s": {a: min(v) for a, v in cpu.items()}}
cpi = {a: pooled[a] * 1e9 / ir[a] for a in pooled}
cpi_min = {a: min(cpu[a]) * 1e9 / ir[a] for a in pooled}
s_hi = {a: reg["arms"][a]["stall_ns_per_Ir_hi"] for a in ir}
s_lo = {a: reg["arms"][a]["stall_ns_per_Ir_lo"] for a in ir}
def F(c):
    d = c["ic"] - c["rho"]
    return {"delta_ns_per_Ir": d, "F_high": (s_hi["ic"] - s_hi["rho"]) / d, "F_low": (s_lo["ic"] - s_lo["rho"]) / d}
out["E3"]["post_hoc_F_pooled_median"] = F(cpi)
out["E3"]["post_hoc_F_min_of_reps"] = F(cpi_min)
out["E3"]["ns_per_Ir_pooled"] = cpi
out["E3"]["IPC_pooled"] = {a: 1 / (cpi[a] * 1e-9 * cal["f_hz_median"]) for a in cpi}

# ---- E1 ----
med = lambda k, kind: st.median(r[kind][k]["ns_per_step_median"] for r in itf["reps"])
e1 = {k: med(k, "loaded") / med(k, "alone") - 1 for k in ("pos", "neg")}
out["E1"] = {"slowdown_pos": e1["pos"], "slowdown_neg": e1["neg"], "informative": e1["pos"] >= 0.10 and e1["neg"] <= 0.03}

# ---- E2 as run (cpu = user+sys) ----
e2 = {a: {"cpu_4k": st.median(r["k4"][a]["cpu_s"] for r in thp["reps"]), "cpu_thp": st.median(r["thp"][a]["cpu_s"] for r in thp["reps"])} for a in ("ic", "rho")}
ch = [c["k4"]["ns_per_step_median"] / c["thp"]["ns_per_step_median"] for c in thp["chase_control"]]
out["E2_cpu_as_run"] = {"medians": e2, "chase_speedup_runs": ch,
                        "anon_huge_kb_ic_max": max(r["thp"]["ic"]["max_anon_huge_kb"] for r in thp["reps"]),
                        "ratio_4k_over_thp": {a: e2[a]["cpu_4k"] / e2[a]["cpu_thp"] for a in e2}}

# ---- A2: user-time share ----
m = lambda cond, a, key: st.median(r[cond][a][key] for r in a2["reps"])
u4 = {a: m("k4", a, "user_s") for a in ("ic", "rho")}
ut = {a: m("thp", a, "user_s") for a in ("ic", "rho")}
s4 = {a: m("k4", a, "sys_s") for a in ("ic", "rho")}
stt = {a: m("thp", a, "sys_s") for a in ("ic", "rho")}
d = {a: (u4[a] - ut[a]) * 1e9 / ir[a] for a in u4}
delta_user = u4["ic"] * 1e9 / ir["ic"] - u4["rho"] * 1e9 / ir["rho"]
share = (d["ic"] - d["rho"]) / delta_user
peak = st.median(r["k4"]["ic"]["peak_rss_bytes"] or 0 for r in a2["reps"]) or 338_350_080
huge = max(r["thp"]["ic"]["max_anon_huge_kb"] for r in a2["reps"]) * 1024
valid = huge >= 0.5 * peak and min(ch) >= 1.3
sys_inc = (stt["ic"] - s4["ic"]) / u4["ic"]
if not valid:
    reading = "not informative (validity failed)"
elif share >= 0.25:
    reading = "translation share >= 0.25"
elif share < 0.05 and sys_inc < 0.05:
    reading = "translation is not the mechanism"
else:
    reading = "unresolved"
spread = {(c, a): max(r[c][a]["user_s"] for r in a2["reps"]) / min(r[c][a]["user_s"] for r in a2["reps"]) - 1
          for c in ("k4", "thp") for a in ("ic", "rho")}
out["A2"] = {"user_s_4k": u4, "user_s_thp": ut, "sys_s_4k": s4, "sys_s_thp": stt,
             "user_ratio_4k_over_thp": {a: u4[a] / ut[a] for a in u4},
             "d_ns_per_Ir": d, "delta_user_ns_per_Ir": delta_user, "translation_share": share,
             "anon_huge_bytes_max": huge, "ic_peak_rss": peak, "valid": valid,
             "ic_sys_increase_over_user4k": sys_inc, "reading": reading,
             "user_spread_maxmin": {f"{c}/{a}": v for (c, a), v in spread.items()},
             "wall_s_median": {c: {a: m(c, a, "wall_s") for a in ("ic", "rho")} for c in ("k4", "thp")}}
# what the user-time gap looks like relative to the registered cpu-based delta
out["A2"]["user_share_of_total_cpu_4k"] = {a: u4[a] / (u4[a] + s4[a]) for a in u4}

for k, v in out.items():
    say(f"== {k}"); say(json.dumps(v, indent=1))
open(os.path.join(HERE, "analysis_explore.txt"), "w").write("\n".join(lines) + "\n")
json.dump(out, open(os.path.join(HERE, "analysis_explore.json"), "w"), indent=2)
print("\n".join(lines))

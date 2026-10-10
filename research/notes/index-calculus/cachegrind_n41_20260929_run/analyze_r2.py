#!/usr/bin/env python3
"""Apply Registration R2 of RESEARCH_CACHEGRIND_N41_20260929.md to r2.json (registered rules only)."""
import json
import os
import random
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
r2 = json.load(open(os.path.join(HERE, "r2.json")))
reg = json.load(open(os.path.join(HERE, "analysis.json")))
ir = {a: reg["arms"][a]["Ir"] for a in ("ic", "rho")}
B = r2["blocks"]
keys = [f"{a}/{m}" for a in ("rho", "ic") for m in ("k4", "thp")]
q = lambda v: st.quantiles(v, n=4, method="inclusive")
out = {"gates": {}}
lines = []
say = lambda *a: lines.append(" ".join(str(x) for x in a))


def share(idx):
    med = {k: st.median(B[i][k]["user_s"] for i in idx) for k in keys}
    d = {a: (med[f"{a}/k4"] - med[f"{a}/thp"]) * 1e9 / ir[a] for a in ("ic", "rho")}
    delta = med["ic/k4"] * 1e9 / ir["ic"] - med["rho/k4"] * 1e9 / ir["rho"]
    return (d["ic"] - d["rho"]) / delta, med, d, delta


# gates
out["gates"]["R-G1"] = all(B[i][k]["gates_G1"] and B[i][k]["rc"] == 0 for i in range(len(B)) for k in keys)
peak = st.median(B[i]["ic/k4"]["peak_rss_bytes"] for i in range(len(B)))
huge = st.median(B[i]["ic/thp"]["max_anon_huge_kb"] for i in range(len(B))) * 1024
ch = [c["k4"]["ns_per_step_median"] / c["thp"]["ns_per_step_median"] for c in r2["chase_pre"] + r2["chase_post"]]
out["gates"]["R-V"] = huge >= 0.5 * peak and st.median(ch) >= 1.3
out["gates"]["R-V_detail"] = {"median_anon_huge_bytes_ic_thp": huge, "ic_peak_rss": peak, "chase_speedups": ch, "median_chase_speedup": st.median(ch)}
iqr = {a: (q([b[f"{a}/k4"]["user_s"] for b in B])[2] - q([b[f"{a}/k4"]["user_s"] for b in B])[0]) / st.median(b[f"{a}/k4"]["user_s"] for b in B)
       for a in ("ic", "rho")}
out["gates"]["R-N"] = all(v <= 0.10 for v in iqr.values())
out["gates"]["R-N_detail"] = iqr
drift = {k: st.median(B[i][k]["user_s"] for i in range(8, 15)) / st.median(B[i][k]["user_s"] for i in range(0, 7)) - 1 for k in keys}
out["gates"]["R-D"] = all(abs(v) <= 0.10 for v in drift.values())
out["gates"]["R-D_detail"] = drift

# estimate + bootstrap
pt, med, d, delta = share(range(len(B)))
rng = random.Random(20260930)
boots = sorted(share(rng.choices(range(len(B)), k=len(B)))[0] for _ in range(10000))
lo, hi = boots[int(0.025 * len(boots))], boots[int(0.975 * len(boots))]
out["estimate"] = {"translation_share": pt, "ci95": [lo, hi], "median_user_s": med, "d_ns_per_Ir": d, "delta_user_ns_per_Ir": delta,
                   "median_sys_s": {k: st.median(B[i][k]["sys_s"] for i in range(len(B))) for k in keys},
                   "median_wall_s": {k: st.median(B[i][k]["wall_s"] for i in range(len(B))) for k in keys},
                   "ic_user_speedup_thp": med["ic/k4"] / med["ic/thp"], "rho_user_speedup_thp": med["rho/k4"] / med["rho/thp"]}
gates_ok = all(out["gates"][g] for g in ("R-G1", "R-V", "R-N", "R-D"))
if not gates_ok:
    verdict = "UNRESOLVED (gate)"
else:
    tags = []
    if lo >= 0.25:
        tags.append("STRONG (lower bound >= 0.25)")
    elif lo >= 0.10:
        tags.append("SUBSTANTIAL (lower bound >= 0.10)")
    if hi < 0.25:
        tags.append("BELOW A QUARTER (upper bound < 0.25)")
    verdict = " + ".join(tags) if tags else "UNRESOLVED"
out["verdict"] = verdict
out["H_mem"] = ("SUPPORTED (R2 STRONG)" if gates_ok and lo >= 0.25 else
                "partly supported by the measured translation share; remainder unresolved within [F_low, F_high] = [%.3f, %.3f]"
                % (reg["model"]["F_low"], reg["model"]["F_high"]))
for k, v in out.items():
    say(f"== {k}"); say(json.dumps(v, indent=1))
open(os.path.join(HERE, "analysis_r2.txt"), "w").write("\n".join(lines) + "\n")
json.dump(out, open(os.path.join(HERE, "analysis_r2.json"), "w"), indent=2)
print("\n".join(lines))

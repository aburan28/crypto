#!/usr/bin/env python3
"""Apply the registered rules of RESEARCH_CACHEGRIND_N41_20260929.md to the measured stages.

Inputs (same directory): calibration_summary.json, native.json, cachegrind.json, controls.json,
interference.json.  Output: analysis.txt / analysis.json.  Nothing here is tuned after the fact;
every threshold below is the protocol's.
"""
import json
import os
import statistics

HERE = os.path.dirname(os.path.abspath(__file__))
J = lambda name: json.load(open(os.path.join(HERE, name)))
MIB = 1024 * 1024

# registered thresholds
F_SUPPORT_LOW = 0.5      # F_low >= this  -> supported outright
F_SUPPORT_HIGH = 0.5     # F_high >= this and S == sensitive -> supported
F_REFUTE_HIGH = 0.25     # F_high < this -> refuted
F_REFUTE_HIGH_S = 0.5    # S == insensitive and F_high < this -> refuted
NATIVE_SPREAD_MAX = 0.10
NO_GAP_FRACTION = 0.05   # delta-cpi <= this x cpi_rho -> no throughput gap to explain
IR_SWEEP = {"ic": 12_441_253_531, "rho": 8_184_739_666}   # callgrind Ir in the merged sweep, cell n=41
IR_TOL = 1e-3
CTL_LO, CTL_HI = 0.7, 1.3
SENS_MIN, SENS_DELTA, INSENS_MAX = 0.10, 0.08, 0.05
INT_SPREAD_MAX = 0.05

cal = J("calibration_summary.json")
nat, cg, ctl, itf = J("native.json"), J("cachegrind.json"), J("controls.json"), J("interference.json")
out = {"gates": {}, "arms": {}, "notes": []}
lines = []


def say(*a):
    lines.append(" ".join(str(x) for x in a))


def counts(arm, ll):
    for r in cg["runs"]:
        if r["arm"] == arm and r["ll_bytes"] == ll:
            return r["counts"]
    raise KeyError((arm, ll))


def levels(c):
    d1 = c["D1mr"] + c["D1mw"]
    ll = c["DLmr"] + c["DLmw"]
    return d1 - ll, ll   # (D1 miss that hits LL, LL miss)


# ---- gates -------------------------------------------------------------------------------------
g1 = (all(r["gates_G1"] for r in cg["runs"]) and all(x[a]["gates_G1"] for x in nat["reps"] for a in ("rho", "ic"))
      and all(x[k][a]["gates_G1"] for x in itf["reps"] for k in ("alone", "loaded") for a in ("rho", "ic")))
out["gates"]["G1_correct_all_runs"] = g1
ir = {arm: counts(arm, 2 * MIB)["Ir"] for arm in ("ic", "rho")}
irs = {arm: {ll: counts(arm, ll)["Ir"] for ll in (2 * MIB, 8 * MIB, 32 * MIB)} for arm in ("ic", "rho")}
out["gates"]["G2_binary_identity"] = all(
    abs(irs[a][ll] / IR_SWEEP[a] - 1) <= IR_TOL for a in ("ic", "rho") for ll in irs[a])
out["gates"]["G2_detail"] = {a: {str(ll // MIB): irs[a][ll] for ll in irs[a]} for a in irs}
mono = True
for arm in ("ic", "rho"):
    m = [levels(counts(arm, ll))[1] for ll in (2 * MIB, 8 * MIB, 32 * MIB)]
    mono &= m[0] * 1.01 >= m[1] and m[1] * 1.01 >= m[2]
    out["arms"].setdefault(arm, {})["LL_misses_2_8_32MiB"] = m
out["gates"]["G4_monotone_LL_misses"] = mono

# G3 simulator sanity on controls
P = cal
def cg_ctl(label, s, ll):
    for r in ctl["cachegrind"]:
        if r["label"] == label and r["steps"] == s and r["ll_bytes"] == ll:
            return r["counts"]
    raise KeyError((label, s, ll))

nat_ctl = {r["label"]: r for r in ctl["native"]}
g3_rows, g3 = [], True
for lab, mode, w, k in [("chase_512K_k1", "chase", 0, 1), ("chase_8M_k1", "chase", 0, 1),
                        ("chase_256M_k1", "chase", 0, 1), ("chase_256M_k10", "chase", 0, 10)]:
    s1, s0 = ctl["steps"], ctl["s0"]
    diff = {}
    for ll in (2 * MIB, 8 * MIB):
        a, b = cg_ctl(lab, s1, ll), cg_ctl(lab, s0, ll)
        diff[ll] = {e: a[e] - b[e] for e in ("D1mr", "D1mw", "DLmr", "DLmw")}
    mid_h, far_h = levels(diff[2 * MIB])
    mid_l, far_l = levels(diff[8 * MIB])
    upper = (mid_h * P["P_mid_high_ns"] + far_h * P["P_far_high_ns"])
    lower = (mid_l * P["P_mid_low_ns"] + far_l * P["P_far_low_ns"]) / P["MLP_max"]
    base = nat_ctl["base_k%d" % (10 if k == 10 else 1)]["ns_per_step_median"]
    extra = (nat_ctl[lab]["ns_per_step_median"] - base) * (s1 - s0)
    ok = CTL_LO * lower <= extra <= CTL_HI * upper
    g3 &= ok
    g3_rows.append({"control": lab, "native_extra_ns": extra, "model_lower_ns": lower, "model_upper_ns": upper, "in_band": ok,
                    "mid_hi": mid_h, "far_hi": far_h, "far_lo": far_l})
out["gates"]["G3_simulator_sanity"] = g3
out["gates"]["G3_rows"] = g3_rows
sq = nat_ctl["seq_256M"]
seq_diff = {ll: {e: cg_ctl("seq_256M", 2, ll)[e] - cg_ctl("seq_256M", 1, ll)[e] for e in ("D1mr", "D1mw", "DLmr", "DLmw")}
            for ll in (2 * MIB, 8 * MIB)}
lines_per_pass = 256 * MIB // 64
out["seq_blind_spot_info"] = {
    "native_ns_per_line": sq["ns_per_line_median"],
    "model_upper_ns_per_line": levels(seq_diff[2 * MIB])[1] * P["P_far_high_ns"] / lines_per_pass,
    "model_lower_ns_per_line": levels(seq_diff[8 * MIB])[1] * P["P_far_low_ns"] / P["MLP_max"] / lines_per_pass,
}

# ---- native throughput --------------------------------------------------------------------------
cpu = {a: [x[a]["cpu_s"] for x in nat["reps"]] for a in ("ic", "rho")}
spread = {a: max(v) / min(v) - 1 for a, v in cpu.items()}
cpu_med = {a: statistics.median(v) for a, v in cpu.items()}
cpi_ns = {a: cpu_med[a] * 1e9 / ir[a] for a in cpu}          # ns per retired instruction
out["native"] = {"cpu_s_runs": cpu, "spread": spread, "cpu_s_median": cpu_med, "ns_per_Ir": cpi_ns,
                 "IPC": {a: 1.0 / (cpi_ns[a] * 1e-9 * P["f_hz_median"]) for a in cpi_ns},
                 "wall_s_median": {a: statistics.median(x[a]["wall_s"] for x in nat["reps"]) for a in cpu}}
out["gates"]["G5_native_spread"] = all(v <= NATIVE_SPREAD_MAX for v in spread.values())
d_cpi = cpi_ns["ic"] - cpi_ns["rho"]
out["native"]["delta_ns_per_Ir"] = d_cpi
no_gap = d_cpi <= NO_GAP_FRACTION * cpi_ns["rho"]

# ---- stall model --------------------------------------------------------------------------------
def stall_ns_per_ir(arm, hi):
    if hi:   # generous to H_mem: small LL, high prices, no overlap
        mid, far = levels(counts(arm, 2 * MIB))
        return (mid * P["P_mid_high_ns"] + far * P["P_far_high_ns"]) / ir[arm]
    mid, far = levels(counts(arm, 8 * MIB))   # stingy: bigger LL, low prices, full MLP
    return (mid * P["P_mid_low_ns"] + far * P["P_far_low_ns"]) / P["MLP_max"] / ir[arm]

s_hi = {a: stall_ns_per_ir(a, True) for a in ir}
s_lo = {a: stall_ns_per_ir(a, False) for a in ir}
F_high = (s_hi["ic"] - s_hi["rho"]) / d_cpi if d_cpi > 0 else float("nan")
F_low = (s_lo["ic"] - s_lo["rho"]) / d_cpi if d_cpi > 0 else float("nan")
c2 = {a: counts(a, 2 * MIB) for a in ir}
br = {a: c2[a]["Bcm"] / ir[a] * P["branch_ns"] for a in ir}
F_branch = (br["ic"] - br["rho"]) / d_cpi if d_cpi > 0 else float("nan")
for arm in ir:
    c8, c32 = counts(arm, 8 * MIB), counts(arm, 32 * MIB)
    out["arms"][arm].update({
        "Ir": ir[arm], "cpu_s_median": cpu_med[arm], "ns_per_Ir": cpi_ns[arm],
        "D1_misses": c2[arm]["D1mr"] + c2[arm]["D1mw"], "data_refs": c2[arm]["Dr"] + c2[arm]["Dw"],
        "D1_MPKI": 1000 * (c2[arm]["D1mr"] + c2[arm]["D1mw"]) / ir[arm],
        "LLd_MPKI_2MiB": 1000 * levels(c2[arm])[1] / ir[arm],
        "LLd_MPKI_8MiB": 1000 * levels(c8)[1] / ir[arm],
        "LLd_MPKI_32MiB": 1000 * levels(c32)[1] / ir[arm],
        "I1_misses": c2[arm]["I1mr"], "I1_MPKI": 1000 * c2[arm]["I1mr"] / ir[arm],
        "branches": c2[arm]["Bc"], "mispredicts": c2[arm]["Bcm"], "BMPKI": 1000 * c2[arm]["Bcm"] / ir[arm],
        "stall_ns_per_Ir_hi": s_hi[arm], "stall_ns_per_Ir_lo": s_lo[arm], "branch_ns_per_Ir": br[arm]})
out["model"] = {"F_high": F_high, "F_low": F_low, "F_branch_exploratory": F_branch, "no_throughput_gap": no_gap,
                "prices": {k: P[k] for k in ("P_mid_low_ns", "P_mid_high_ns", "P_far_low_ns", "P_far_high_ns", "MLP_max", "branch_ns")}}

# ---- interference -------------------------------------------------------------------------------
def med(kind, arm):
    return statistics.median(r[kind][arm]["cpu_s"] for r in itf["reps"])
def spr(kind, arm):
    v = [r[kind][arm]["cpu_s"] for r in itf["reps"]]
    return max(v) / min(v) - 1
sens = {a: med("loaded", a) / med("alone", a) - 1 for a in ("ic", "rho")}
noisy = any(spr("alone", a) > INT_SPREAD_MAX for a in ("ic", "rho"))
if noisy:
    S = "ambiguous"
elif sens["ic"] >= SENS_MIN and sens["ic"] - sens["rho"] >= SENS_DELTA:
    S = "sensitive"
elif sens["ic"] < INSENS_MAX:
    S = "insensitive"
else:
    S = "ambiguous"
out["interference"] = {"slowdown_minus_1": sens, "alone_spread": {a: spr("alone", a) for a in ("ic", "rho")},
                       "loaded_spread": {a: spr("loaded", a) for a in ("ic", "rho")}, "S": S,
                       "cpu_s_alone": {a: med("alone", a) for a in ("ic", "rho")},
                       "cpu_s_loaded": {a: med("loaded", a) for a in ("ic", "rho")}}

# ---- registered verdict --------------------------------------------------------------------------
gates_ok = all(out["gates"][g] for g in ("G1_correct_all_runs", "G3_simulator_sanity", "G4_monotone_LL_misses", "G5_native_spread"))
if not gates_ok:
    verdict = "UNDETERMINED (a registered gate failed)"
elif no_gap:
    verdict = "NO THROUGHPUT GAP TO EXPLAIN (H_mem moot at this cell)"
elif F_low >= F_SUPPORT_LOW or (F_high >= F_SUPPORT_HIGH and S == "sensitive"):
    verdict = "SUPPORTED"
elif F_high < F_REFUTE_HIGH or (S == "insensitive" and F_high < F_REFUTE_HIGH_S):
    verdict = "REFUTED as the main explanation"
else:
    verdict = "UNDETERMINED"
out["verdict"] = verdict
out["gates_all_blocking_ok"] = gates_ok

say("== gates ==")
for g, v in out["gates"].items():
    if g not in ("G3_rows", "G2_detail"):
        say(f"  {g}: {v}")
say("  G2 (non-blocking; binary identity vs sweep callgrind Ir):", json.dumps(out["gates"]["G2_detail"]))
say("  G3 rows:")
for r in g3_rows:
    say(f"    {r['control']:<16} native_extra={r['native_extra_ns']:.3e} ns  model=[{r['model_lower_ns']:.3e}, {r['model_upper_ns']:.3e}]  in_band={r['in_band']}")
say("  seq blind-spot (informational):", json.dumps(out["seq_blind_spot_info"]))
say("== native ==")
say("  ", json.dumps(out["native"]))
say("== per-arm simulated ==")
for a in ("ic", "rho"):
    say(f"  {a}:", json.dumps(out["arms"][a]))
say("== model ==")
say("  ", json.dumps(out["model"]))
say("== interference ==")
say("  ", json.dumps(out["interference"]))
say("== VERDICT (registered rule):", verdict)
open(os.path.join(HERE, "analysis.txt"), "w").write("\n".join(lines) + "\n")
json.dump(out, open(os.path.join(HERE, "analysis.json"), "w"), indent=2)
print("\n".join(lines))

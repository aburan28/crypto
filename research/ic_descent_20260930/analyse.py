#!/usr/bin/env python3
"""Ledger §22's analysis: every number the note and the page quote.

Reads runs/ (both arms' pricer reports, the controls, the probes) and
§20's analysis.json (its batch-rho price per size, which this round does
not re-measure), grades the five declared targets, and writes
analysis.json.  Computes nothing the reports do not carry except ratios,
means and intervals.

The primary speedup converts both arms at one unit: it is the paired
ratio of total time, since §21.4 found the unit itself moving between
binaries.  The ratio in each process's own unit is reported beside it.

    python3 analyse.py > analysis.json
"""
from __future__ import annotations

import json
import math
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
S20 = json.loads((HERE.parent / "ic_exponent_20260926" / "analysis.json").read_text())
SIZES = [(1, 19), (1, 23), (1, 45), (0, 37), (1, 43), (1, 47), (0, 41), (0, 53), (0, 61)]
PROBED = [(1, 19), (1, 23), (1, 45), (0, 37), (1, 43), (0, 41)]
SETS = (1, 2, 3, 4)
ROUNDS = 5
T95 = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365, 8: 2.306, 9: 2.262,
       10: 2.228, 11: 2.201, 12: 2.179, 13: 2.160, 14: 2.145, 15: 2.131, 16: 2.120, 17: 2.110,
       18: 2.101, 19: 2.093, 20: 2.086}
# Target 4's sizes: the probe put the fixed part at 25% or more of the descent.
TARGET4 = {(1, 19), (1, 23), (1, 45), (0, 37)}
# Target 3's sizes.
TARGET3 = {(1, 19), (1, 23), (1, 45)}


def load(path: Path) -> dict | None:
    return json.loads(path.read_text()) if path.exists() else None


def geo_ci(ratios: list[float]) -> dict:
    """Geometric mean and its 95% interval, `t` on the logarithms."""
    logs = [math.log(x) for x in ratios]
    m = statistics.fmean(logs)
    out = {"geomean": math.exp(m), "n": len(logs)}
    if len(logs) > 1:
        h = T95[len(logs) - 1] * statistics.stdev(logs) / math.sqrt(len(logs))
        out.update({"lo": math.exp(m - h), "hi": math.exp(m + h)})
    return out


def mean_ci(values: list[float]) -> dict:
    """§20's interval: arithmetic mean, `t` over the sets."""
    m = statistics.fmean(values)
    h = T95[len(values) - 1] * statistics.stdev(values) / math.sqrt(len(values))
    return {"mean": m, "lo": m - h, "hi": m + h, "n": len(values)}


def fit(xs: list[float], ys: list[float]) -> dict:
    """§20's declared fit, least squares with a `t` interval on the slope."""
    n = len(xs)
    xm, ym = statistics.fmean(xs), statistics.fmean(ys)
    sxx = sum((x - xm) ** 2 for x in xs)
    beta = sum((x - xm) * (y - ym) for x, y in zip(xs, ys)) / sxx
    alpha = ym - beta * xm
    dof = n - 2
    s2 = sum((y - alpha - beta * x) ** 2 for x, y in zip(xs, ys)) / dof
    h = T95[dof] * math.sqrt(s2 / sxx)
    return {"beta": beta, "lo": beta - h, "hi": beta + h, "points": n, "dof": dof}


def pinned(a: dict, b: dict) -> bool:
    return (a.get("status") == b.get("status") == "complete"
            and a.get("counts") == b.get("counts") and a.get("recovered") == b.get("recovered")
            and a.get("all_verified") is True and b.get("all_verified") is True)


def pairs_of(d: Path, tag: str = "", count: int = ROUNDS) -> list[tuple[dict, dict]]:
    out = []
    for i in range(1, count + 1):
        a = load(d / f"{tag}r{i}-baseline.price.json")
        b = load(d / f"{tag}r{i}-candidate.price.json")
        if a is not None and b is not None:
            out.append((a, b))
    return out


def med(rep: dict, key: str) -> float:
    return statistics.median(r[key] for r in rep["repetitions"])


def phase_ns(rep: dict, phases: tuple[str, ...]) -> float:
    return sum(statistics.median(r["phases_ns"][ph] for r in rep["repetitions"]) for ph in phases)


def s20_row(a: int, n: int) -> dict:
    return next(r for r in S20["sizes"] if (r["a"], r["n"]) == (a, n))


def size_row(a: int, n: int) -> dict:
    d0 = RUNS / "main" / f"k{a}n{n}"
    sets, primary, own, descent, unit = [], [], [], [], []
    all_pairs, pins_ok = 0, True
    s_before_sets, ratio_after_sets, s_own_after_sets = [], [], []
    prior = s20_row(a, n)
    for j in SETS:
        d = d0 / f"M{j}"
        first = pairs_of(d)
        double = pairs_of(d, "double-")
        figure = double if double else first
        all_pairs += len(first) + len(double)
        pins_ok &= all(pinned(x, y) for x, y in first + double)
        base_totals = [x["median"]["total_units"] for x, _ in first]
        wall = [med(x, "total_ns") / med(y, "total_ns") for x, y in figure]
        mine = [x["median"]["total_units"] / y["median"]["total_units"] for x, y in figure]
        stage = [phase_ns(x, ("descent",)) / phase_ns(y, ("descent",)) for x, y in figure]
        primary += wall
        own += mine
        descent += stage
        unit += [med(x, "unit_ns") / med(y, "unit_ns") for x, y in figure]
        s_before = statistics.median(x["median"]["s_per_target"] for x, _ in figure)
        s_after = s_before / geo_ci(wall)["geomean"]
        s_own_after = statistics.median(y["median"]["s_per_target"] for _, y in figure)
        rho = next(x["s_rho"] for x in prior["sets"] if x["set"] == f"M{j}")
        s_before_sets.append(s_before / rho)
        ratio_after_sets.append(s_after / rho)
        s_own_after_sets.append(s_own_after)
        sets.append({
            "set": f"M{j}", "rounds": len(first), "doubled": bool(double),
            "aa_spread": max(base_totals) / min(base_totals) if base_totals else None,
            "speedup": geo_ci(wall), "speedup_own_unit": geo_ci(mine), "descent_speedup": geo_ci(stage),
            "pins": all(pinned(x, y) for x, y in first + double),
            "s_before": s_before, "s_after_baseline_unit": s_after, "s_after_own_unit": s_own_after,
            "descent_share_before": statistics.median(
                x["median"]["phases_units"]["descent"] / x["median"]["total_units"] for x, _ in figure),
        })
    s_rho = prior["s_rho"]
    speed = geo_ci(primary)
    s_before = statistics.fmean(s["s_before"] for s in sets)
    s_after = statistics.fmean(s["s_after_baseline_unit"] for s in sets)
    control = load(RUNS / "control1" / f"k{a}n{n}-M1-control.json")
    return {
        "curve": prior["curve"], "a": a, "n": n, "log2_r": prior["log2_r"],
        "pairs": all_pairs, "pins": pins_ok,
        "control1": control.get("pass") if control else None,
        "speedup": speed,
        "speedup_own_unit": geo_ci(own),
        "descent_speedup_stage_diagnostic": geo_ci(descent),
        "unit_ns_baseline_over_candidate": geo_ci(unit),
        "s_before": s_before, "s_after": s_after,
        "s_after_own_unit": statistics.fmean(s_own_after_sets),
        "s_rho_s20": s_rho, "s20_ratio": prior["ratio"],
        "ratio_before": s_before / s_rho, "ratio_after": s_after / s_rho,
        "ratio_before_ci": mean_ci(s_before_sets), "ratio_after_ci": mean_ci(ratio_after_sets),
        "descent_share_before": statistics.fmean(s["descent_share_before"] for s in sets),
        "sets": sets,
    }


def probe_row(a: int, n: int) -> dict | None:
    base = load(RUNS / "probe" / f"k{a}n{n}-M1-baseline.json")
    cand = load(RUNS / "probe" / f"k{a}n{n}-M1-candidate.json")
    if not base or not cand:
        return None

    def parts(j: dict) -> dict:
        ph = {x["phase"]: x["mean_units_per_target"] for x in j["phases"]}
        return {"solve": j["solve_mean_units_per_target"], "target_query": ph["target_query"],
                "target_pdp": ph["target_pdp"], "target_descent": ph["target_descent"],
                "recovery_check": ph["recovery_check"],
                "relation_and_check": ph["target_descent"] + ph["recovery_check"]}

    b, c = parts(base), parts(cand)
    return {
        "curve": f"k{a}n{n}", "trials_equal": base["trials_per_target"] == cand["trials_per_target"],
        "trials_per_target": base["trials_total"] / len(base["trials_per_target"]),
        "unit_ns": {"baseline": base["unit_ns"], "candidate": cand["unit_ns"]},
        "baseline": b, "candidate": c,
        "solve_fall": b["solve"] / c["solve"],
        "relation_and_check_fall": b["relation_and_check"] / c["relation_and_check"],
        "target_query_fall": b["target_query"] / c["target_query"],
    }


def main() -> None:
    sizes = [size_row(a, n) for a, n in SIZES if (RUNS / "main" / f"k{a}n{n}").exists()]
    thread_pairs = pairs_of(RUNS / "threads" / "k0n41-M1", "", 3)
    threads = {"curve": "k0n41", "set": "M1", "pairs": len(thread_pairs),
               "pins": all(pinned(x, y) for x, y in thread_pairs),
               "speedup": geo_ci([med(x, "total_ns") / med(y, "total_ns") for x, y in thread_pairs])
               if thread_pairs else None}
    probes = [p for p in (probe_row(a, n) for a, n in PROBED) if p]
    rho = load(RUNS / "rho" / "k0n41-M1-verdict.json")
    by = {(s["a"], s["n"]): s for s in sizes}
    t1 = {"pairs": sum(s["pairs"] for s in sizes), "all_pinned": all(s["pins"] for s in sizes),
          "control1_all_pass": all(s["control1"] is True for s in sizes),
          "probe_trials_equal": all(p["trials_equal"] for p in probes)}
    t1["pass"] = t1["all_pinned"] and t1["control1_all_pass"] and t1["pairs"] >= 180 and t1["probe_trials_equal"]
    t2 = {"per_size": {p["curve"]: p["relation_and_check_fall"] for p in probes},
          "pass": len(probes) == 6 and all(p["relation_and_check_fall"] >= 3 for p in probes)}
    t3_rows = [p for p in probes if tuple(int(x) for x in p["curve"][1:].split("n")) in TARGET3]
    t3 = {"per_size": {p["curve"]: p["solve_fall"] for p in t3_rows},
          "pass": len(t3_rows) == 3 and all(p["solve_fall"] >= 1.8 for p in t3_rows)}
    t4_rows = [by[k] for k in sorted(TARGET4) if k in by]
    t4 = {"lower_bounds": {s["curve"]: s["speedup"]["lo"] for s in t4_rows},
          "pass": len(t4_rows) == 4 and all(s["speedup"]["lo"] > 1 for s in t4_rows)}
    regress = [s["curve"] for s in sizes if s["speedup"]["hi"] < 1]
    if threads["speedup"] and threads["speedup"].get("hi", 2) < 1:
        regress.append("K_0/GF(2^41) at four threads")
    t5 = {"regressions": regress, "pass": not regress}
    top = sorted(sizes, key=lambda s: s["log2_r"])[-4:]
    exponent = None
    if len(top) == 4:
        xs = [s["log2_r"] * math.log(2) for s in top]
        exponent = {
            "sizes": [s["curve"] for s in top],
            "before": fit(xs, [math.log(s["ratio_before"] * math.sqrt(s["n"])) for s in top]),
            "after": fit(xs, [math.log(s["ratio_after"] * math.sqrt(s["n"])) for s in top]),
        }
    least = min(sizes, key=lambda s: s["ratio_after"]) if sizes else None
    doc = {
        "what_this_is": "Ledger §22: baseline (main 57e7ce3a) against the candidate, 36 of §20's parameter "
                        "files, five ABAB rounds each; the primary speedup is the paired ratio of total time "
                        "(both arms at one unit); the ratio in each process's own unit is beside it.",
        "host": load(HERE / "host.json"),
        "sizes": sizes, "threads": threads, "rho_control": rho, "probes": probes,
        "exponent_refit": exponent,
        "least_ratio_after": {"curve": least["curve"], "log2_r": least["log2_r"], "ratio": least["ratio_after"],
                              "ci": least["ratio_after_ci"]} if least else None,
        "targets": {"1_identical_outputs": t1, "2_relation_and_check_3x": t2, "3_descent_1_8x": t3,
                    "4_speedup_excludes_1": t4, "5_no_regression": t5},
    }
    print(json.dumps(doc, indent=1))


if __name__ == "__main__":
    main()

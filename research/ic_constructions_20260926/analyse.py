#!/usr/bin/env python3
"""Ledger §21's analysis: every number the note and the page quote.

Reads runs/ (the pricer's reports of both arms, the controls, the
per-constructor prices) and §20's analysis.json (its batch-rho price per
size, which this round does not re-measure), grades the four declared
targets, and writes analysis.json.  Computes nothing the reports do not
carry except ratios, means and intervals.

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
SETS = (1, 2, 3, 4)
ROUNDS = 5
T95 = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365, 8: 2.306, 9: 2.262,
       10: 2.228, 11: 2.201, 12: 2.179, 13: 2.160, 14: 2.145, 15: 2.131, 16: 2.120, 17: 2.110,
       18: 2.101, 19: 2.093, 20: 2.086}
# §20 found the constructions at 15% of S or more at every size but these.
TARGET3_EXEMPT = {(0, 53), (0, 61)}


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


def pinned(a: dict, b: dict) -> bool:
    return (a.get("status") == b.get("status") == "complete"
            and a.get("counts") == b.get("counts") and a.get("recovered") == b.get("recovered")
            and a.get("all_verified") is True and b.get("all_verified") is True)


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


def pairs_of(d: Path, tag: str = "", count: int = ROUNDS) -> list[tuple[dict, dict]]:
    out = []
    for i in range(1, count + 1):
        a = load(d / f"{tag}r{i}-baseline.price.json")
        b = load(d / f"{tag}r{i}-candidate.price.json")
        if a is not None and b is not None:
            out.append((a, b))
    return out


# The constructions group PROTOCOL.md declares: the four phases that
# build the consumers around the work.  Not the report's own
# "constructions" group, which also holds the curve's setup — out of
# scope here, and unchanged.
CONSTRUCTIONS = ("select_projection", "collect_setup", "logs_setup", "descent_setup")


def constructions(rep: dict) -> float:
    return sum(rep["median"]["phases_units"][ph] for ph in CONSTRUCTIONS)


# Phases the change does not touch: the base, table, collection, logs and
# descent work, and the curve's setup.
UNCHANGED = ("setup", "select", "build", "collect", "verify", "la", "descent", "verify_final")


def med(rep: dict, key: str) -> float:
    return statistics.median(r[key] for r in rep["repetitions"])


def phase_ns(rep: dict, phases: tuple[str, ...]) -> float:
    return sum(statistics.median(r["phases_ns"][ph] for r in rep["repetitions"]) for ph in phases)


def s20_row(a: int, n: int) -> dict:
    return next(r for r in S20["sizes"] if (r["a"], r["n"]) == (a, n))


def size_row(a: int, n: int) -> dict:
    d0 = RUNS / "main" / f"k{a}n{n}"
    sets, speed, cons, all_pairs = [], [], [], 0
    wall, unit, still = [], [], []
    pins_ok = True
    s_arm = {"baseline": [], "candidate": []}
    share = {"baseline": [], "candidate": []}
    for j in SETS:
        d = d0 / f"M{j}"
        first = pairs_of(d)
        double = pairs_of(d, "double-")
        figure = double if double else first
        all_pairs += len(first) + len(double)
        pins_ok &= all(pinned(x, y) for x, y in first + double)
        base_totals = [x["median"]["total_units"] for x, _ in first]
        sp = [x["median"]["total_units"] / y["median"]["total_units"] for x, y in figure]
        cs = [constructions(x) / constructions(y) for x, y in figure]
        speed += sp
        cons += cs
        # Secondary, not declared: the same pairs in nanoseconds, the
        # unit's own shift between the binaries, and the untouched phases.
        wall += [med(x, "total_ns") / med(y, "total_ns") for x, y in figure]
        unit += [med(x, "unit_ns") / med(y, "unit_ns") for x, y in figure]
        still += [phase_ns(x, UNCHANGED) / phase_ns(y, UNCHANGED) for x, y in figure]
        for arm, k in (("baseline", 0), ("candidate", 1)):
            s = statistics.median(p[k]["median"]["s_per_target"] for p in figure)
            s_arm[arm].append(s)
            share[arm].append(statistics.median(constructions(p[k]) / p[k]["median"]["total_units"]
                                                for p in figure))
        sets.append({
            "set": f"M{j}", "rounds": len(first), "doubled": bool(double),
            "aa_spread": max(base_totals) / min(base_totals) if base_totals else None,
            "speedup": geo_ci(sp) if sp else None,
            "constructions_speedup": geo_ci(cs) if cs else None,
            "pins": all(pinned(x, y) for x, y in first + double),
            "s_baseline": s_arm["baseline"][-1], "s_candidate": s_arm["candidate"][-1],
            "units_per_process": {
                "baseline": [x["median"]["total_units"] for x, _ in figure],
                "candidate": [y["median"]["total_units"] for _, y in figure]},
            "unit_ns": {
                "baseline": statistics.median(r["unit_ns"] for x, _ in figure for r in x["repetitions"]),
                "candidate": statistics.median(r["unit_ns"] for _, y in figure for r in y["repetitions"])},
        })
    prior = s20_row(a, n)
    s_rho = prior["s_rho"]
    s_before = statistics.fmean(s_arm["baseline"])
    s_after = statistics.fmean(s_arm["candidate"])
    # Each set against §20's batch rho on that set's own 32 targets.
    rho_sets = [next(x["s_rho"] for x in prior["sets"] if x["set"] == f"M{j}") for j in SETS]
    control = load(RUNS / "control1" / f"k{a}n{n}-M1-control.json")
    return {
        "curve": prior["curve"], "a": a, "n": n, "log2_r": prior["log2_r"],
        "pairs": all_pairs, "pins": pins_ok,
        "control1": control.get("pass") if control else None,
        "speedup": geo_ci(speed), "constructions_speedup_stage_diagnostic": geo_ci(cons),
        "secondary_not_declared": {
            "wall_speedup": geo_ci(wall),
            "unit_ns_baseline_over_candidate": geo_ci(unit),
            "unchanged_phases_ns_baseline_over_candidate": geo_ci(still),
        },
        "s_before": s_before, "s_after": s_after, "s20_s_ic": prior["s_ic"], "s_rho_s20": s_rho,
        "ratio_before": s_before / s_rho, "ratio_after": s_after / s_rho, "s20_ratio": prior["ratio"],
        "ratio_before_ci": mean_ci([x / y for x, y in zip(s_arm["baseline"], rho_sets)]),
        "ratio_after_ci": mean_ci([x / y for x, y in zip(s_arm["candidate"], rho_sets)]),
        "constructions_share": {arm: statistics.fmean(v) for arm, v in share.items()},
        "s20_constructions_share": prior["shares"]["constructions"],
        "sets": sets,
    }


def threads_row(a: int, n: int) -> dict:
    d = RUNS / "threads" / f"k{a}n{n}-M1"
    pr = pairs_of(d, "", 3)
    sp = [x["median"]["total_units"] / y["median"]["total_units"] for x, y in pr]
    return {"curve": f"k{a}n{n}", "set": "M1", "pairs": len(pr), "pins": all(pinned(x, y) for x, y in pr),
            "threads": pr[0][0].get("threads") if pr else None,
            "speedup": geo_ci(sp) if sp else None}


def constructions_rows() -> list[dict]:
    out = []
    for a, n in [(0, 41), (0, 61)]:
        arms = {arm: load(RUNS / "constructions" / f"k{a}n{n}-M1-{arm}.json") for arm in ("baseline", "candidate")}
        if not all(arms.values()):
            continue
        names = list(arms["baseline"]["constructors"])
        rows = {}
        for name in names:
            b, c = arms["baseline"]["constructors"][name], arms["candidate"]["constructors"][name]
            rows[name] = {
                "baseline_first_per_point": b["first_units"] / arms["baseline"]["points"],
                "baseline_median_per_point": b["per_point"],
                "candidate_first_per_point": c["first_units"] / arms["candidate"]["points"],
                "candidate_median_per_point": c["per_point"],
            }
        out.append({"curve": f"k{a}n{n}", "points": arms["baseline"]["points"],
                    "signed_orbits": arms["baseline"]["signed_orbits"],
                    "cofactor_bits": arms["baseline"]["cofactor_bits"], "constructors": rows})
    return out


def main() -> None:
    sizes = [size_row(a, n) for a, n in SIZES if (RUNS / "main" / f"k{a}n{n}").exists()]
    threads = [threads_row(a, n) for a, n in [(0, 53), (0, 61)] if (RUNS / "threads" / f"k{a}n{n}-M1").exists()]
    rho = load(RUNS / "rho" / "k0n41-M1-verdict.json")
    t1 = {"pairs": sum(s["pairs"] for s in sizes), "all_pinned": all(s["pins"] for s in sizes),
          "control1_all_pass": all(s["control1"] is True for s in sizes),
          "pass": all(s["pins"] for s in sizes) and all(s["control1"] is True for s in sizes)
          and sum(s["pairs"] for s in sizes) >= 180}
    t2 = {"per_size": {s["curve"]: s["constructions_speedup_stage_diagnostic"]["geomean"] for s in sizes},
          "pass": all(s["constructions_speedup_stage_diagnostic"]["geomean"] >= 3 for s in sizes)}
    need = [s for s in sizes if (s["a"], s["n"]) not in TARGET3_EXEMPT]
    t3 = {"sizes": [s["curve"] for s in need],
          "lower_bounds": {s["curve"]: s["speedup"]["lo"] for s in need},
          "pass": len(need) == 7 and all(s["speedup"]["lo"] > 1 for s in need)}
    regress = [s["curve"] for s in sizes if s["speedup"]["hi"] < 1]
    regress += [t["curve"] + " (4 threads)" for t in threads if t["speedup"] and t["speedup"].get("hi", 2) < 1]
    t4 = {"regressions": regress, "pass": not regress}
    # §20's declared exponent fit, re-derived on each arm: log(ratio·√n)
    # against log r over the four largest sizes.  A derived figure, not a
    # target of this round.
    top = sorted(sizes, key=lambda s: s["log2_r"])[-4:]
    exponent = None
    if len(top) == 4:
        xs = [s["log2_r"] * math.log(2) for s in top]
        exponent = {
            "sizes": [s["curve"] for s in top],
            "before": fit(xs, [math.log(s["ratio_before"] * math.sqrt(s["n"])) for s in top]),
            "after": fit(xs, [math.log(s["ratio_after"] * math.sqrt(s["n"])) for s in top]),
            # Secondary, not declared: the candidate's cost in the baseline
            # binary's unit (S before over the paired wall-clock speedup),
            # which takes the unit's own shift between binaries out.
            "after_at_baseline_unit_secondary": fit(xs, [
                math.log(s["s_before"] / s["secondary_not_declared"]["wall_speedup"]["geomean"] / s["s_rho_s20"]
                         * math.sqrt(s["n"])) for s in top]),
        }
    least = min(sizes, key=lambda s: s["ratio_after"]) if sizes else None
    doc = {
        "what_this_is": "Ledger §21: baseline e7022b75 against the candidate, 36 of §20's parameter files, "
                        "five ABAB rounds each; speedup = baseline total units / candidate total units per pair.",
        "host": load(HERE / "host.json"),
        "sizes": sizes, "threads": threads, "rho_control": rho, "constructions": constructions_rows(),
        "exponent_refit": exponent,
        "least_ratio_after": {"curve": least["curve"], "log2_r": least["log2_r"], "ratio": least["ratio_after"],
                              "ci": least["ratio_after_ci"]} if least else None,
        "targets": {"1_identical_outputs": t1, "2_constructions_3x": t2, "3_speedup_excludes_1": t3,
                    "4_no_regression": t4},
    }
    print(json.dumps(doc, indent=1))


if __name__ == "__main__":
    main()

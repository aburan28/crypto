#!/usr/bin/env python3
"""Every figure ledger §23 and the scoreboard quote (PROTOCOL.md, "Reading").

Reads the checked rows (claims/, from claims.py) and the reports they
point to; computes nothing claims.py has not already checked.

Per size, over the R1 rows (T01-T64):

- the online speedup: mean rho online over mean index-calculus online,
  in each row's own unit, with a 95% bootstrap interval (targets
  resampled in pairs, 10,000 resamples, fixed seed);
- the same with rho at the canonical step (the model, secondary);
- the cold ratio, index calculus over rho, each arm's reusable set-up
  plus its mean online cost, with its interval;
- the break-even count of targets, set-up over the online saving;
- the precomputation boundary: Bernstein and Lange's generic walk with
  the index calculus's own set-up as its budget, a model;
- the A/A: T01-T04's second process against its first.

Across sizes: least-squares slopes of ln S on ln r with t intervals.

    IC_RUNS=runs python3 analyse.py > analysis.json
"""
from __future__ import annotations

import json
import math
import os
import random
import re
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_params  # noqa: E402

RUNS = (HERE / os.environ.get("IC_RUNS", "runs")).resolve()
BASE = (HERE / os.environ.get("IC_OUT", ".")).resolve()
CLAIMS = BASE / "claims"
BOOT = 10_000
BL_PRODUCT = 1.93 * 1.21
# t quantiles, 97.5%, by degrees of freedom.
T975 = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365, 8: 2.306}


def load_rows(a: int, n: int) -> tuple[list[dict], list[dict], list[dict]]:
    """(R1 claims, R2 claims, rows that did not pass)."""
    r1, r2, bad = [], [], []
    d = CLAIMS / f"k{a}n{n}"
    for path in sorted(d.glob("T*-R*.claim.json")):
        doc = json.loads(path.read_text())
        claim, ok = doc["claim"], doc["validation"]["status"] == "PASS"
        m = re.fullmatch(r"T(\d+)-R(\d+)\.claim\.json", path.name)
        claim["_target"], claim["_run"] = int(m[1]), int(m[2])
        if not (ok and claim["independent_validation"]):
            bad.append({"row": path.name, "status": doc["validation"]["status"],
                        "errors": doc["validation"]["validation_errors"]})
            continue
        (r1 if claim["_run"] == 1 else r2).append(claim)
    return r1, r2, bad


def report_of(claim: dict) -> dict:
    path = Path(claim["report"])
    return json.loads((path if path.is_absolute() else HERE.parents[1] / path).read_text())


def boot_ci(rows: list[dict], stat, seed: int) -> list[float]:
    rng = random.Random(seed)
    vals = []
    for _ in range(BOOT):
        sample = [rows[rng.randrange(len(rows))] for _ in rows]
        v = stat(sample)
        if v is not None and math.isfinite(v):
            vals.append(v)
    vals.sort()
    return [vals[int(0.025 * len(vals))], vals[int(0.975 * len(vals)) - 1]]


def mean(rows: list[dict], key: str) -> float:
    return statistics.fmean(r["units"][key] for r in rows)


def fit(points: list[tuple[float, float]]) -> dict:
    """Slope of y on x with its 95% t interval; None below two sizes."""
    k = len(points)
    if k < 2:
        return None
    xs, ys = [p[0] for p in points], [p[1] for p in points]
    mx, my = statistics.fmean(xs), statistics.fmean(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    slope = sum((x - mx) * (y - my) for x, y in points) / sxx
    icept = my - slope * mx
    if k < 3:
        return {"slope": slope, "points": k}
    resid = sum((y - icept - slope * x) ** 2 for x, y in points)
    se = math.sqrt(resid / (k - 2) / sxx)
    t = T975.get(k - 2, 1.96)
    return {"slope": slope, "ci95": [slope - t * se, slope + t * se], "points": k}


def isolation(a: int, n: int) -> dict:
    d = RUNS / f"k{a}n{n}"
    runs = contended = failed = 0
    for rec in d.glob("*.isolation.jsonl"):
        for line in rec.read_text().splitlines():
            run = json.loads(line)["run"]
            runs += 1
            contended += bool(run["contended"])
            failed += run["exit_status"] != 0
    return {"processes": runs, "contended": contended, "failed": failed,
            "retries": len(list(d.glob("*-retry*.price.json")))}


def size(a: int, n: int, seed: int) -> dict | None:
    r1, r2, bad = load_rows(a, n)
    if not r1:
        return None
    reports = [report_of(c) for c in r1]
    r = reports[0]["r"]
    sqrt_r = math.sqrt(r)
    canonical = statistics.median(rep["step_prices"]["canonical_step_units"] for rep in reports)
    steps = [rep["repetitions"][0]["rho_online"]["steps"] for rep in reports]
    floor_steps = math.sqrt(math.pi * r / (4 * n))

    def speedup(rows: list[dict]) -> float:
        return mean(rows, "rho_online_units") / mean(rows, "ic_online_units")

    def speedup_model(rows: list[dict]) -> float:
        return mean(rows, "rho_model_units") / mean(rows, "ic_online_units")

    setup = statistics.median(c["units"]["setup_units"] for c in r1)
    rho_setup = statistics.median(c["units"]["rho_setup_units"] for c in r1)

    def cold(rows: list[dict]) -> float:
        return (setup + mean(rows, "ic_online_units")) / (rho_setup + mean(rows, "rho_online_units"))

    def break_even(rows: list[dict]) -> float | None:
        saving = mean(rows, "rho_online_units") - mean(rows, "ic_online_units")
        return setup / saving if saving > 0 else None

    bl_online = BL_PRODUCT * r * canonical ** 2 / (2 * n * setup)
    per_row = [c["online_speedup"] for c in r1]
    aa = []
    for c2 in r2:
        c1 = next((c for c in r1 if c["_target"] == c2["_target"]), None)
        if c1:
            aa.append({"target": c2["_target"],
                       "online_speedup_R2_over_R1": c2["online_speedup"] / c1["online_speedup"],
                       "ic_online_R2_over_R1": c2["units"]["ic_online_units"] / c1["units"]["ic_online_units"],
                       "rho_online_R2_over_R1": c2["units"]["rho_online_units"] / c1["units"]["rho_online_units"]})
    return {
        "curve": r1[0]["curve"], "a": a, "n": n, "r": r, "log2_r": round(math.log2(r), 2),
        "curve_id": r1[0]["curve_id"], "candidate_id": r1[0]["candidate_id"],
        "rho_reference_uid": r1[0]["rho_reference_uid"],
        "rows": {"R1": len(r1), "R2": len(r2), "not_passing": bad},
        "all_rows_pass_the_check": not bad,
        "online_speedup": {"mean_ratio": speedup(r1), "ci95": boot_ci(r1, speedup, seed),
                           "median_of_rows": statistics.median(per_row),
                           "rows_above_one": sum(1 for s in per_row if s > 1)},
        "online_speedup_rho_model": {"mean_ratio": speedup_model(r1), "ci95": boot_ci(r1, speedup_model, seed + 1)},
        "cold_ratio_ic_over_rho": {"value": cold(r1), "ci95": boot_ci(r1, cold, seed + 2)},
        "break_even_targets": {"value": break_even(r1), "ci95": boot_ci(r1, break_even, seed + 3)},
        "s_ic_online_mean": mean(r1, "ic_online_units") / sqrt_r,
        "s_rho_online_mean": mean(r1, "rho_online_units") / sqrt_r,
        "s_rho_model_mean": mean(r1, "rho_model_units") / sqrt_r,
        "s_setup": setup / sqrt_r,
        "s_rho_setup": rho_setup / sqrt_r,
        "s_floor_steps": math.sqrt(math.pi / (4 * n)),
        "rho_steps_over_floor_mean": statistics.fmean(steps) / floor_steps,
        "rho_units_per_step_median": statistics.median(c["units"]["rho_units_per_step"] for c in r1),
        "canonical_step_units": canonical,
        "rho_step_over_model_median": statistics.median(c["units"]["rho_step_over_model"] for c in r1),
        "precomputation_boundary_model": {
            "what": "Bernstein-Lange: a generic walk given the index calculus's set-up as its budget, "
                    "online ~ 1.93*1.21 r c^2 / (2n P) units; a model",
            "bl_online_units": bl_online, "s_bl_online": bl_online / sqrt_r,
            "ic_online_over_bl": mean(r1, "ic_online_units") / bl_online},
        "ic_replay_ns_median": statistics.median(c["units"]["ic_replay_ns"] for c in r1),
        "rho_replay_ns_median": statistics.median(c["units"]["rho_replay_ns"] for c in r1),
        "aa": aa,
        "isolation": isolation(a, n),
        "diagnostics": diagnostics(a, n, r1, reports, setup, rho_setup),
    }


def diagnostics(a: int, n: int, r1: list[dict], reports: list[dict], setup: float, rho_setup: float) -> dict:
    """Read-outs beside the declared figures, labelled as such: where the
    online time goes, the curve's construction (in both arms' set-up), the
    cold ratio without it and with rho at the canonical step, and the online
    interval against §22's per-target descent at k = 32."""
    phases = {}
    for rep in reports:
        at = rep["median"]["ic_online_repetition"]
        r = rep["repetitions"][at]
        for k, v in r["ic_online"]["phases_ns"].items():
            phases.setdefault(k, []).append((v or 0) / r["unit_ns"])
    online = mean(r1, "ic_online_units")
    curve = statistics.median(statistics.median(x["setup_phases_ns"]["setup"] / x["unit_ns"]
                                                for x in rep["repetitions"]) for rep in reports)
    rho_online, rho_model = mean(r1, "rho_online_units"), mean(r1, "rho_model_units")
    out = {
        "online_phase_shares": {k: statistics.fmean(v) / online for k, v in phases.items()},
        "curve_construction_units": curve,
        "s_curve_construction": curve / math.sqrt(reports[0]["r"]),
        "cold_ratio_without_curve_construction": (setup - curve + online) / rho_online,
        "cold_ratio_rho_at_canonical_step": (setup + online) / (rho_setup + rho_model),
        "online_trials_mean": statistics.fmean(rep["repetitions"][0]["ic_online"]["trials"] for rep in reports),
    }
    s22 = sorted((HERE.parents[1] / "research" / "ic_descent_20260930" / "runs-isolated" / "main"
                  / f"k{a}n{n}").glob("M*/r*-candidate.price.json"))
    trials, per_trial, per_target = [], [], []
    for path in s22:
        d = json.loads(path.read_text())
        if d.get("status") != "complete":
            continue
        t = d["counts"]["descent"]["trials_per_target"]
        descent = (d["median"]["phases_units"]["descent"] + d["median"]["phases_units"]["verify_final"]) / d["targets"]
        trials.append(statistics.fmean(t))
        per_trial.append(d["median"]["phases_units"]["descent"] / d["targets"] / statistics.fmean(t))
        per_target.append(descent)
    if trials:
        ours = sum(rep["median"]["ic_online_units"] for rep in reports) / sum(
            rep["repetitions"][0]["ic_online"]["trials"] for rep in reports)
        out["against_s22_batch_descent"] = {
            "s22_files": len(trials),
            "online_over_s22_per_target_descent": online / statistics.median(per_target),
            "trials_ratio": out["online_trials_mean"] / statistics.fmean(trials),
            "units_per_trial_ratio": ours / statistics.median(per_trial),
        }
    return out


def main() -> None:
    sizes = []
    for k, (a, n) in enumerate(make_params.SIZES):
        s = size(a, n, 2309300 + 10 * k)
        if s:
            sizes.append(s)
    fits = {}
    for key in ("s_ic_online_mean", "s_rho_online_mean", "s_setup", "s_rho_model_mean"):
        fits[key] = fit([(math.log(s["r"]), math.log(s[key])) for s in sizes])
    fits["online_speedup"] = fit([(math.log(s["r"]), math.log(s["online_speedup"]["mean_ratio"])) for s in sizes])
    fits["cold_ratio"] = fit([(math.log(s["r"]), math.log(s["cold_ratio_ic_over_rho"]["value"])) for s in sizes])
    prediction = json.loads((HERE / "prediction.json").read_text())
    against = []
    for s in sizes:
        p = next(row for row in prediction["rows"] if (row["a"], row["n"]) == (s["a"], s["n"]))
        against.append({"curve": s["curve"],
                        "online_speedup": [p["online_speedup_probe_step"], s["online_speedup"]["mean_ratio"]],
                        "online_speedup_model": [p["online_speedup_canonical_step"],
                                                 s["online_speedup_rho_model"]["mean_ratio"]],
                        "cold_ratio": [p["cold_ratio_probe_step"], s["cold_ratio_ic_over_rho"]["value"]],
                        "ic_online_over_bl": [p["ic_online_over_bl_model"],
                                              s["precomputation_boundary_model"]["ic_online_over_bl"]]})
    print(json.dumps({"what_this_is": "Ledger §23's figures, from the checked rows.",
                      "sizes": sizes, "fits_ln_s_on_ln_r": fits,
                      "prediction_then_measurement": against}, indent=1))


if __name__ == "__main__":
    main()

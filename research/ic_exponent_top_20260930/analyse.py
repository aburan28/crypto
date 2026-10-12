#!/usr/bin/env python3
"""Ledger §23's analysis: every number the note and the page quote.

Reads the run directory (the pricer's reports with their isolation
records, the sweep choices, the Control 1 verdicts) and prediction.json,
grades the declared targets, and prints the analysis.  §20's arithmetic,
unchanged: a size's ratio is the mean index-calculus S over the mean
priced batch-rho S, with a 95% interval from the per-set ratios; the fit
is ln(ratio·√n) on ln r, least squares with a t interval on the slope.

A set's figure is its first clean run (uncontended, exit 0) among the
run and its retries; when §20's spread rule fired, the doubled rerun's
first clean run.  Contended runs are counted, never pooled.

    IC_RUNS=runs python3 analyse.py > analysis.json
"""
from __future__ import annotations

import json
import math
import os
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = (HERE / os.environ.get("IC_RUNS", "runs")).resolve()
PRED = json.loads((HERE / "prediction.json").read_text())
S20 = HERE.parent / "ic_exponent_20260926"
SETS = range(1, 9)
RETRIES = 2
T95 = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365, 8: 2.306, 9: 2.262,
       10: 2.228, 12: 2.179, 14: 2.145, 16: 2.120, 20: 2.086, 30: 2.042, 40: 2.021, 46: 2.013, 60: 2.000}


def t95(dof: int) -> float:
    keys = sorted(k for k in T95 if k <= dof)
    return T95[keys[-1]] if keys else float("nan")


def load(path: Path) -> dict | None:
    return json.loads(path.read_text()) if path.exists() and path.stat().st_size else None


def record_path(out: Path) -> Path:
    name = out.name
    for suffix in (".price.json", ".json"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return out.with_name(name + ".isolation.jsonl")


def state(out: Path) -> str:
    rec = record_path(out)
    if not rec.exists():
        return "unrecorded"
    run = json.loads(rec.read_text().splitlines()[-1])["run"]
    if run["exit_status"] != 0:
        return "failed"
    return "contended" if run["contended"] else "clean"


def attempts(base: Path) -> list[Path]:
    stem = base.name[: -len(".json")]
    return [base] + [base.with_name(f"{stem}-retry{k}.json") for k in range(1, RETRIES + 1)]


def first_clean(base: Path) -> tuple[dict | None, dict]:
    tally = {"runs": 0, "contended": 0, "failed": 0}
    for p in attempts(base):
        rep = load(p)
        if rep is None:
            break
        tally["runs"] += 1
        st = state(p)
        if st in ("contended", "failed"):
            tally[st] += 1
            continue
        return rep, tally
    return None, tally


def figure(base: Path) -> tuple[dict | None, dict]:
    """The set's figure and its accounting (§20's spread rule, then the first clean run)."""
    rep, tally = first_clean(base)
    doubled = False
    if rep is not None and rep.get("spread_max_over_min", 0) > 1.25:
        again, t2 = first_clean(base.with_name(base.name[: -len(".json")] + "-double.json"))
        for k in tally:
            tally[k] += t2[k]
        if again is not None:
            rep, doubled = again, True
    return rep, {**tally, "doubled": doubled}


def mean_ci(values: list[float]) -> dict:
    m = statistics.fmean(values)
    if len(values) < 2:
        return {"mean": m, "lo": None, "hi": None, "n": len(values)}
    s = statistics.stdev(values)
    h = t95(len(values) - 1) * s / math.sqrt(len(values))
    return {"mean": m, "lo": m - h, "hi": m + h, "sd": s, "n": len(values)}


def fit(xs: list[float], ys: list[float]) -> dict:
    n = len(xs)
    xm, ym = statistics.fmean(xs), statistics.fmean(ys)
    sxx = sum((x - xm) ** 2 for x in xs)
    beta = sum((x - xm) * (y - ym) for x, y in zip(xs, ys)) / sxx
    alpha = ym - beta * xm
    resid = [y - alpha - beta * x for x, y in zip(xs, ys)]
    dof = n - 2
    s2 = sum(e * e for e in resid) / dof if dof > 0 else float("nan")
    se = math.sqrt(s2 / sxx) if dof > 0 else float("nan")
    h = t95(dof) * se if dof > 0 else float("nan")
    return {"beta": beta, "alpha": alpha, "se": se, "lo": beta - h, "hi": beta + h, "points": n, "dof": dof,
            "residuals": resid}


def read(f: dict, lo: float, hi: float) -> dict:
    return {"interval": [lo, hi], "excludes_zero": not (lo <= 0 <= hi),
            "law_1_6_consistent": lo <= 1 / 6 <= hi}


def s20_choice(a: int, n: int) -> dict | None:
    p = S20 / "runs" / f"k{a}n{n}" / "sweep" / "chosen.json"
    return json.loads(p.read_text())["chosen"] if p.exists() else None


def size_row(prow: dict) -> dict:
    a, n, r = prow["a"], prow["n"], prow["r"]
    d = RUNS / f"k{a}n{n}"
    chosen = json.loads((d / "sweep" / "chosen.json").read_text())
    sets, acct = [], {"runs": 0, "contended": 0, "failed": 0, "doubled_sets": 0, "sets_without_a_clean_run": 0}
    for j in SETS:
        rep, tally = figure(d / "measure" / f"M{j}.price.json")
        for k in ("runs", "contended", "failed"):
            acct[k] += tally[k]
        acct["doubled_sets"] += int(tally["doubled"])
        if rep is None:
            acct["sets_without_a_clean_run"] += int(tally["runs"] > 0)
            continue
        ok = rep.get("status") == "complete"
        rho = rep.get("rho_batch", {})
        s_ic, s_rho = (rep["median"]["s_per_target"] if ok else None), rho.get("s_priced_per_target")
        sets.append({
            "set": f"M{j}", "status": rep.get("status"), "doubled": tally["doubled"],
            "spread": rep.get("spread_max_over_min"), "repetitions": len(rep.get("repetitions", [])),
            "all_verified": rep.get("all_verified"),
            "counts_identical": rep.get("counts_identical_across_repetitions"),
            "rho_verified": rho.get("all_verified"), "rho_agrees": rho.get("agrees_with_expected"),
            "rejected": rep.get("counts", {}).get("logs", {}).get("rejected") if ok else None,
            "points": rep.get("counts", {}).get("select", {}).get("points") if ok else None,
            "stored_pairs": rep.get("counts", {}).get("build", {}).get("stored_pairs") if ok else None,
            "tier": rep.get("counts", {}).get("build", {}).get("tier") if ok else None,
            "s_ic": s_ic, "s_rho": s_rho,
            "ratio": (s_ic / s_rho) if (s_ic and s_rho) else None,
            "phases_units": rep["median"]["phases_units"] if ok else None,
            "unit_ns": statistics.median([x["unit_ns"] for x in rep.get("repetitions", [])]) if ok else None,
            "s_cold": rep["median"].get("s_cold") if ok else None,
            "rho_cold_s_priced": rep.get("rho_cold", {}).get("s_priced"),
            "floor": rep.get("floor", {}).get("per_target_k"),
        })
    control = load(d / "measure" / "M1.params-control.json")
    good = [s for s in sets if s["ratio"] is not None]
    row = {
        "curve": prow["curve"], "a": a, "n": n, "r": r, "log2_r": prow["log2_r"], "new_in_23": prow["new_in_23"],
        "chosen": chosen["chosen"], "sweep": chosen, "s20_chosen": s20_choice(a, n),
        "control_1_M1": control.get("pass") if control else None,
        "accounting": acct, "sets": sets, "model": prow["model_ratio"], "law": prow["law_ratio"],
    }
    if good:
        s_ic = statistics.fmean(s["s_ic"] for s in good)
        s_rho = statistics.fmean(s["s_rho"] for s in good)
        per_phase = {ph: statistics.fmean(s["phases_units"][ph] for s in good) for ph in good[0]["phases_units"]}
        total = sum(per_phase.values())
        row.update({
            "s_ic": s_ic, "s_rho": s_rho, "ratio": s_ic / s_rho,
            "ratio_ci": mean_ci([s["ratio"] for s in good]),
            "ic_over_floor": s_ic / good[0]["floor"],
            "phase_shares": {ph: v / total for ph, v in per_phase.items()},
            "unit_ns": statistics.fmean(s["unit_ns"] for s in good),
        })
        cold = [s for s in good if s.get("rho_cold_s_priced")]
        if cold:
            row["cold_ratio_M1"] = cold[0]["s_cold"] / cold[0]["rho_cold_s_priced"]
    return row


def sweep_accounting() -> dict:
    tally = {"runs": 0, "contended": 0, "failed": 0}
    for size in PRED["rows"]:
        d = RUNS / f"k{size['a']}n{size['n']}" / "sweep"
        for p in sorted(d.glob("*.price*.json")):
            tally["runs"] += 1
            st = state(p)
            if st in ("contended", "failed"):
                tally[st] += 1
    return tally


def main() -> None:
    rows = [size_row(p) for p in PRED["rows"]
            if (RUNS / f"k{p['a']}n{p['n']}" / "sweep" / "chosen.json").exists()]
    done = sorted([r for r in rows if "ratio" in r], key=lambda r: r["r"])
    sets_all = [s for r in rows for s in r["sets"]]
    t1 = {
        "every_target_verified": all(s["all_verified"] and s["rho_verified"] and s["rho_agrees"] for s in sets_all),
        "no_rejected_relation": all(s["rejected"] == 0 for s in sets_all if s["rejected"] is not None),
        "counts_identical_across_repetitions": all(s["counts_identical"] for s in sets_all),
        "every_run_complete": all(s["status"] == "complete" for s in sets_all),
        "eight_sets_every_size": all(len(r["sets"]) == 8 for r in rows) and len(rows) == 6,
    }
    t1["met"] = all(t1.values())
    t2 = {"control_1_M1_every_size": all(r["control_1_M1"] is True for r in rows) and len(rows) == 6}
    exponent = {}
    if len(done) == 6:
        xs = [math.log(r["r"]) for r in done]
        ys = [math.log(r["ratio"] * math.sqrt(r["n"])) for r in done]
        f6 = fit(xs, ys)
        model6 = PRED["model_local_exponent_six_sizes_n_corrected"]
        exponent["primary_six_sizes"] = {
            "sizes": [r["curve"] for r in done], "fit": f6, **read(f6, f6["lo"], f6["hi"]),
            "model": model6, "model_consistent": f6["lo"] <= model6 <= f6["hi"],
            "aim_met": (not (f6["lo"] <= 0 <= f6["hi"])) or (not (f6["lo"] <= 1 / 6 <= f6["hi"])),
        }
        top = done[-4:]
        f4 = fit([math.log(r["r"]) for r in top], [math.log(r["ratio"] * math.sqrt(r["n"])) for r in top])
        model4 = PRED["model_local_exponent_top_four_n_corrected"]
        exponent["secondary_four_largest"] = {
            "sizes": [r["curve"] for r in top], "fit": f4, **read(f4, f4["lo"], f4["hi"]),
            "model": model4, "model_consistent": f4["lo"] <= model4 <= f4["hi"],
        }
        pts = [(math.log(r["r"]), math.log(s["ratio"] * math.sqrt(r["n"]))) for r in done for s in r["sets"]
               if s["ratio"]]
        exponent["diagnostic_per_set_points"] = {
            "note": "per-set points are not independent sizes; a diagnostic only, never the declared fit",
            "fit": fit([p[0] for p in pts], [p[1] for p in pts]),
        }
    least = min(done, key=lambda r: r["ratio"]) if done else None
    crossings = [r["curve"] for r in done if r["ratio_ci"]["hi"] is not None and r["ratio_ci"]["hi"] < 1]
    refusals = RUNS / "refusals.log"
    waits = RUNS / "psi-waits.log"
    wait_s = [float(l.split()[-1]) for l in waits.read_text().splitlines()] if waits.exists() else []
    report = {
        "what_this_is": "Ledger §23's analysis: six sizes, eight sets each, every timed process isolated; "
                        "§20's arithmetic unchanged.",
        "runs": RUNS.name,
        "host": load(RUNS / "host.json"),
        "sizes": rows,
        "isolation": {
            "sets": {k: sum(r["accounting"][k] for r in rows)
                     for k in ("runs", "contended", "failed", "doubled_sets", "sets_without_a_clean_run")},
            "sweep": sweep_accounting(),
            "refused_starts": len(refusals.read_text().splitlines()) if refusals.exists() else 0,
            "psi_wait_seconds": {"processes": len(wait_s), "median": statistics.median(wait_s) if wait_s else None,
                                 "max": max(wait_s) if wait_s else None},
        },
        "least_ratio": {"curve": least["curve"], "log2_r": least["log2_r"], "ratio": least["ratio"],
                        "ci": least["ratio_ci"]} if least else None,
        "targets": {"1_correct": t1, "2_control_1": t2, "3_4_exponent": exponent,
                    "5_crossing": {"sizes_with_interval_below_one": crossings}},
    }
    print(json.dumps(report, indent=1, default=str))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Ledger §20's analysis: every number the note and the page quote.

Reads runs/ (the pricer's reports, the sweep choices, the control
verdicts) and prediction.json, grades the five declared targets, and
writes analysis.json.  Computes nothing the reports do not carry except
means, intervals and the declared fit.

    python3 analyse.py > analysis.json
"""
from __future__ import annotations

import json
import math
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
PRED = json.loads((HERE / "prediction.json").read_text())
T95 = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365, 8: 2.306,
       10: 2.228, 12: 2.179, 14: 2.145, 16: 2.120, 20: 2.086, 30: 2.042}


def t95(dof: int) -> float:
    keys = sorted(k for k in T95 if k <= dof)
    return T95[keys[-1]] if keys else float("nan")


def figure(path: Path) -> dict | None:
    """The report a set's figure comes from: the doubled rerun when the
    spread rule fired, else the run itself."""
    if not path.exists():
        return None
    double = path.with_name(path.stem + "-double.json")
    return json.loads((double if double.exists() else path).read_text())


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
    return {"beta": beta, "alpha": alpha, "se": se, "lo": beta - h, "hi": beta + h, "points": n, "dof": dof}


def phase_row(rep: dict) -> dict:
    return rep["median"]["phases_units"]


def size_row(prow: dict) -> dict:
    a, n, r = prow["a"], prow["n"], prow["r"]
    d = RUNS / f"k{a}n{n}"
    chosen = json.loads((d / "sweep" / "chosen.json").read_text())
    sets = []
    for j in (1, 2, 3, 4):
        rep = figure(d / "measure" / f"M{j}.price.json")
        if rep is None:
            continue
        ctl_path = d / "measure" / f"M{j}.params-control.json"
        control = json.loads(ctl_path.read_text()) if ctl_path.exists() else None
        ok = rep.get("status") == "complete"
        s_ic = rep["median"]["s_per_target"] if ok else None
        rho = rep.get("rho_batch", {})
        s_rho = rho.get("s_priced_per_target")
        groups = rep["median"].get("groups", {}) if ok else {}
        k = rep.get("targets", 32)
        sqrt_r = math.sqrt(r)
        sets.append({
            "set": f"M{j}", "status": rep.get("status"), "message": rep.get("message"),
            "doubled": (d / "measure" / f"M{j}.price-double.json").exists(),
            "spread": rep.get("spread_max_over_min"),
            "repetitions": len(rep.get("repetitions", [])),
            "all_verified": rep.get("all_verified"),
            "counts_identical": rep.get("counts_identical_across_repetitions"),
            "rho_verified": rho.get("all_verified"), "rho_agrees": rho.get("agrees_with_expected"),
            "control_pass": control.get("pass") if control else None,
            "rejected": rep.get("counts", {}).get("logs", {}).get("rejected") if ok else None,
            "points": rep.get("counts", {}).get("select", {}).get("points") if ok else None,
            "columns": rep.get("counts", {}).get("select", {}).get("columns") if ok else None,
            "stored_pairs": rep.get("counts", {}).get("build", {}).get("stored_pairs") if ok else None,
            "tier": rep.get("counts", {}).get("build", {}).get("tier") if ok else None,
            "relations": rep.get("counts", {}).get("logs", {}).get("relations_accepted") if ok else None,
            "summands_scanned": rep.get("counts", {}).get("collect", {}).get("summands_scanned") if ok else None,
            "descent_trials": rep.get("counts", {}).get("descent", {}).get("trials_total") if ok else None,
            "s_ic": s_ic,
            "s_rho_counted": rho.get("s_counted_per_target"),
            "s_rho": s_rho,
            "s_rho_bailey_model": rho.get("s_bailey_model_per_target"),
            "canonical_step_units": rep.get("step_prices", {}).get("canonical_step_units"),
            "bailey_step_units": rep.get("step_prices", {}).get("bailey_step_units"),
            "ratio": (s_ic / s_rho) if (s_ic and s_rho) else None,
            "s_work": groups.get("work", 0) / (k * sqrt_r) if ok else None,
            "s_constructions": groups.get("constructions", 0) / (k * sqrt_r) if ok else None,
            "s_verification": groups.get("verification", 0) / (k * sqrt_r) if ok else None,
            "phases_units": phase_row(rep) if ok else None,
            "unit_ns": statistics.median([x["unit_ns"] for x in rep.get("repetitions", [])]) if ok else None,
            "first_over_median": (rep["first_repetition_total_units"] / rep["median"]["total_units"]) if ok else None,
            "s_cold": rep["median"].get("s_cold") if ok else None,
            "rho_cold_s_priced": rep.get("rho_cold", {}).get("s_priced"),
            "rho_cold_verified": rep.get("rho_cold", {}).get("all_verified"),
            "floor": rep.get("floor", {}).get("per_target_k"),
        })
    good = [s for s in sets if s["ratio"] is not None]
    row = {
        "curve": prow["curve"], "a": a, "n": n, "r": r, "log2_r": prow["log2_r"],
        "chosen": chosen["chosen"], "sweep": chosen,
        "sets": sets,
        "law": prow["law_ratio"], "model": prow["model_ratio"],
    }
    if good:
        s_ic = statistics.fmean(s["s_ic"] for s in good)
        s_rho = statistics.fmean(s["s_rho"] for s in good)
        s_work = statistics.fmean(s["s_work"] for s in good)
        row.update({
            "s_ic": s_ic, "s_rho": s_rho,
            "s_rho_counted": statistics.fmean(s["s_rho_counted"] for s in good),
            "ratio": s_ic / s_rho,
            "ratio_ci": mean_ci([s["ratio"] for s in good]),
            "ratio_bailey_model": s_ic / statistics.fmean(s["s_rho_bailey_model"] for s in good),
            "work_ratio_stage_diagnostic": s_work / s_rho,
            "shares": {
                "work": statistics.fmean(s["s_work"] / s["s_ic"] for s in good),
                "constructions": statistics.fmean(s["s_constructions"] / s["s_ic"] for s in good),
                "verification": statistics.fmean(s["s_verification"] / s["s_ic"] for s in good),
            },
            "phase_shares": {
                ph: statistics.fmean(s["phases_units"][ph] / sum(s["phases_units"].values()) for s in good)
                for ph in good[0]["phases_units"]
            },
            "ic_over_floor": s_ic / good[0]["floor"],
            "rho_counted_over_floor": statistics.fmean(s["s_rho_counted"] for s in good) / good[0]["floor"],
            "canonical_step_units": statistics.fmean(s["canonical_step_units"] for s in good),
            "bailey_step_units": statistics.fmean(s["bailey_step_units"] for s in good),
            "unit_ns": statistics.fmean(s["unit_ns"] for s in good),
        })
        cold = [s for s in good if s.get("rho_cold_s_priced")]
        if cold:
            row["cold_ratio_M1"] = cold[0]["s_cold"] / cold[0]["rho_cold_s_priced"]
            row["cold_s_ic_M1"] = cold[0]["s_cold"]
            row["cold_s_rho_M1"] = cold[0]["rho_cold_s_priced"]
    return row


def controls() -> dict:
    d = RUNS / "controls"
    out = {}
    frozen = figure(d / "frozen-headline.price.json")
    if frozen:
        c = frozen.get("counts", {})
        expected = {"relations": 197, "trials": 3450, "summands_scanned": 883200, "stored_pairs": 1519296,
                    "descent_trials": 52}
        got = {
            "relations": c.get("logs", {}).get("relations_accepted"),
            "trials": c.get("collect", {}).get("trials"),
            "summands_scanned": c.get("collect", {}).get("summands_scanned"),
            "stored_pairs": c.get("build", {}).get("stored_pairs"),
            "descent_trials": c.get("descent", {}).get("trials_total"),
        }
        ctl = d / "frozen-headline.params-control.json"
        out["frozen_headline"] = {
            "status": frozen.get("status"),
            "counts_expected": expected, "counts_now": got,
            "counts_reproduced": expected == got,
            "differences": {k: {"frozen": expected[k], "now": got[k]} for k in expected if expected[k] != got[k]},
            "control_1": json.loads(ctl.read_text()).get("pass") if ctl.exists() else None,
            "s_ic": frozen.get("median", {}).get("s_per_target"),
            "s_rho": frozen.get("rho_batch", {}).get("s_priced_per_target"),
            "ratio": frozen.get("ratios", {}).get("ic_over_batch_rho"),
            "groups": frozen.get("median", {}).get("groups"),
            "phases_units": frozen.get("median", {}).get("phases_units"),
            "was_section_19_5": 6.38,
        }
    for label in ("thread-n41", "thread-n53"):
        sets = []
        for j in (1, 2, 3, 4):
            rep = figure(d / f"{label}-M{j}.price.json")
            if rep is None:
                continue
            ctl = d / f"{label}-M{j}.params-control.json"
            sets.append({
                "set": f"M{j}", "status": rep.get("status"),
                "s_ic": rep.get("median", {}).get("s_per_target"),
                "s_rho": rep.get("rho_batch", {}).get("s_priced_per_target"),
                "ratio": rep.get("ratios", {}).get("ic_over_batch_rho"),
                "control_1": json.loads(ctl.read_text()).get("pass") if ctl.exists() else None,
                "all_verified": rep.get("all_verified"),
            })
        good = [s for s in sets if s["ratio"]]
        out[label] = {"sets": sets}
        if good:
            s_ic = statistics.fmean(s["s_ic"] for s in good)
            s_rho = statistics.fmean(s["s_rho"] for s in good)
            out[label].update({"s_ic": s_ic, "s_rho": s_rho, "ratio": s_ic / s_rho,
                               "ratio_ci": mean_ci([s["ratio"] for s in good])})
    return out


def main() -> None:
    rows = [size_row(p) for p in PRED["rows"] if (RUNS / f"k{p['a']}n{p['n']}" / "sweep" / "chosen.json").exists()]
    done = [r for r in rows if "ratio" in r]
    sets_all = [s for r in rows for s in r["sets"]]
    target1 = {
        "every_target_verified": all(s["all_verified"] and s["rho_verified"] and s["rho_agrees"] for s in sets_all),
        "no_rejected_relation": all(s["rejected"] == 0 for s in sets_all if s["rejected"] is not None),
        "counts_identical_across_repetitions": all(s["counts_identical"] for s in sets_all),
        "every_run_complete": all(s["status"] == "complete" for s in sets_all),
    }
    target1["met"] = all(target1.values())
    ctl = controls()
    target2 = {
        "control_1_every_measurement_run": all(s["control_pass"] for s in sets_all),
        "control_2_frozen_counts_reproduced": ctl.get("frozen_headline", {}).get("counts_reproduced"),
    }
    top = sorted(done, key=lambda r: r["r"])[-4:]
    exponent = None
    if len(top) == 4:
        f = fit([math.log(r["r"]) for r in top], [math.log(r["ratio"] * math.sqrt(r["n"])) for r in top])
        pts = [(math.log(r["r"]), math.log(s["ratio"] * math.sqrt(r["n"]))) for r in top for s in r["sets"] if s["ratio"]]
        f16 = fit([p[0] for p in pts], [p[1] for p in pts])
        law, model = 1 / 6, PRED["model_local_exponent_top_four_n_corrected"]
        exponent = {
            "sizes": [r["curve"] for r in top], "fit": f, "fit_per_set_points": f16,
            "law_1_6_consistent": f["lo"] <= law <= f["hi"],
            "model_consistent": f["lo"] <= model <= f["hi"],
            "model_exponent": model,
        }
    small = sorted(done, key=lambda r: r["r"])[:3]
    law_small = None
    if len(small) == 3:
        law_small = {
            "sizes": [r["curve"] for r in small],
            "ratio_ci_lo": [r["ratio_ci"]["lo"] for r in small],
            "twice_law": [2 * r["law"] for r in small],
            "law_falsified": all(r["ratio_ci"]["lo"] is not None and r["ratio_ci"]["lo"] > 2 * r["law"] for r in small),
        }
    minimum = None
    if done:
        least = min(done, key=lambda r: r["ratio"])
        minimum = {
            "curve": least["curve"], "log2_r": least["log2_r"], "ratio": least["ratio"], "ci": least["ratio_ci"],
            "model_minimum_confirmed": least["curve"] in ("K_0/GF(2^37)", "K_1/GF(2^45)", "K_1/GF(2^43)")
            and 2.85 / 2 <= least["ratio"] <= 2.85 * 2,
        }
    crossings = [r["curve"] for r in done if r["ratio_ci"]["hi"] is not None and r["ratio_ci"]["hi"] < 1]
    report = {
        "what_this_is": "Ledger §20's analysis, from runs/ and prediction.json.",
        "host": json.loads((HERE / "host.json").read_text()),
        "sizes": rows,
        "controls": ctl,
        "targets": {
            "1_correct": target1,
            "2_controls": target2,
            "3_exponent": exponent,
            "4_small_end": {"law": law_small, "model_minimum": minimum},
            "5_crossing": {"sizes_with_interval_below_one": crossings},
        },
    }
    print(json.dumps(report, indent=1, default=str))


if __name__ == "__main__":
    main()

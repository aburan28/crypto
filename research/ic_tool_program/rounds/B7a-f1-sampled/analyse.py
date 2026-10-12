#!/usr/bin/env python3
"""B7a's figures and decision (PROTOCOL.md, "Measurements" and
"Acceptance", with amendment 1), from the run tree only.

    IC_RUNS=<run tree> python3 run.py --steps <accepted>,B7a analyse

Measurements 2-4 come from `harness/bround.py`'s analysis. B7a adds:

- **5, F1 against F0.** Per size, F0's cold cost is the median over
  measurement 4's B7a-arm processes (two rows, five rounds), each the
  median repetition's set-up plus online interval, in that process's own
  units. F1's figures are the medians over its six processes (two rows,
  three rounds) of each of `extrapolated.ic.cold.units`'s median and
  bounds. The ratio is F1's median over F0's; F0 is inside when it lies
  in F1's predictive interval; the wall ratio is F1's process wall time
  over F0's, median against median.
- **6, the count's check and the carried constants.** Per size, F1's
  `κ_sample` with its interval; F0's `ρ` (relations accepted over
  columns) and `κ` (summands scanned per relation collected over
  `probes(μ, w)`, with F1's count `μ` and window `w`), each over F0's
  processes; and a least-squares fit of F0's `κ` against `n`.
- **7, partial tables.** Per size and fraction, against the full table's
  F1 runs on the same rows: the build per stored pair, the scan per
  summand and the descent per probe, each a ratio; and the partial
  table's distinct keys over the fraction squared times the full one's.

**The target** (the design's amendment 1): F0 inside F1's predictive
interval at 9 or more of the 11 sizes, and F1's ratio within
`[0.85, 1.18]` at `n = 53`, `59` and `61`.
"""
from __future__ import annotations

import math
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "harness"))
sys.path.insert(0, str(HERE))
import bench  # noqa: E402
import bround  # noqa: E402
import run as b7a  # noqa: E402
import stats  # noqa: E402

TIGHT_SIZES = (53, 59, 61)
TIGHT_BAND = (0.85, 1.18)
INSIDE_NEEDED = 9


def probes(mu: float, granule: float) -> float:
    return granule / -math.expm1(-granule / mu)


def f0_cold_units(rep: dict) -> float:
    return statistics.median((r["setup_ns"] + r["ic_online"]["wall_ns"]) / r["unit_ns"]
                             for r in rep["repetitions"])


def wall(out: Path) -> float | None:
    run = bench.run_record(out)
    return run["wall_seconds"] if run else None


def f0_reports(row: dict) -> list[tuple[Path, dict]]:
    out = []
    for k in range(1, bround.ROUNDS + 1):
        p = bench.figure_path(bround.runs() / "timing" / "cand" / row["id"] / f"r{k}.price.json")
        rep = bench.load(p)
        if bench.clean(p) and rep.get("status") == "complete":
            out.append((p, rep))
    return out


def f1_reports(row: dict) -> list[tuple[Path, dict]]:
    out = []
    for k in range(1, b7a.F1_ROUNDS + 1):
        p = b7a.figure(bround.runs() / "f1" / f"r{k}" / f"{row['id']}.price.json")
        rep = bench.load(p)
        if bench.clean(p) and rep.get("status") == "extrapolated":
            out.append((p, rep))
    return out


def med(values: list[float]) -> float | None:
    values = [v for v in values if v is not None and math.isfinite(v)]
    return statistics.median(values) if values else None


def measurement5(sizes: dict) -> list[dict]:
    out = []
    for (a, n), rows in sizes:
        f0 = [x for row in rows for x in f0_reports(row)]
        f1 = [x for row in rows for x in f1_reports(row)]
        f0_cold = med([f0_cold_units(rep) for _, rep in f0])
        bounds = {key: med([rep["extrapolated"]["ic"]["cold"]["units"][key] for _, rep in f1])
                  for key in ("median", "expectation_lo", "expectation_hi", "predictive_lo", "predictive_hi")}
        ratio = bounds["median"] / f0_cold if f0_cold and bounds["median"] else None
        inside = (f0_cold is not None and bounds["predictive_lo"] is not None
                  and bounds["predictive_lo"] <= f0_cold <= bounds["predictive_hi"])
        wall_f0, wall_f1 = med([wall(p) for p, _ in f0]), med([wall(p) for p, _ in f1])
        out.append({
            "slug": bench.curve_slug(a, n), "a": a, "n": n, "log2_r": round(math.log2(rows[0]["r"]), 3),
            "f0_processes": len(f0), "f1_processes": len(f1),
            "f0_cold_units": f0_cold, "f1_cold_units": bounds,
            "ratio_f1_over_f0": ratio, "f0_inside_predictive": inside,
            "wall_f1_over_f0": wall_f1 / wall_f0 if wall_f0 and wall_f1 else None,
            "in_tight_band": (TIGHT_BAND[0] <= ratio <= TIGHT_BAND[1]) if ratio and n in TIGHT_SIZES else None,
        })
    return out


def measurement6(sizes: dict) -> dict:
    rows_out, ns, kappas = [], [], []
    for (a, n), rows in sizes:
        f1 = [x for row in rows for x in f1_reports(row)]
        f0 = [x for row in rows for x in f0_reports(row)]
        if not f1:
            rows_out.append({"slug": bench.curve_slug(a, n), "n": n, "missing": "no F1 report"})
            continue
        count = f1[0][1]["phases"]["ic"]["collect"]["count"]
        sample = f1[0][1]["phases"]["ic"]["collect"]["kappa_sample"]
        rhos, ks = [], []
        for _, rep in f0:
            pc = rep["counts"]["pass"]
            rhos.append(pc["logs"]["relations_accepted"] / pc["select"]["columns"])
            found = pc["collect"]["relations_collected"]
            if found:
                ks.append(pc["collect"]["summands_scanned"] / found / probes(count["mu"], count["window"]))
        kappa_f0 = med(ks)
        if kappa_f0 is not None:
            ns.append(n)
            kappas.append(kappa_f0)
        rows_out.append({
            "slug": bench.curve_slug(a, n), "n": n,
            "kappa_sample": sample["kappa"], "kappa_sample_interval": sample["interval"],
            "relations_found": sample["relations_found"],
            "f0_rho": med(rhos), "f0_kappa": kappa_f0,
            "carried": f1[0][1]["phases"]["ic"]["collect"]["carried"],
        })
    return {"sizes": rows_out, "fit_f0_kappa_on_n": stats.fit([float(x) for x in ns], kappas)}


def measurement7(sizes: dict) -> list[dict]:
    out = []
    largest = {(r["a"], r["n"]) for r in b7a.largest_rows()}
    for (a, n), rows in sizes:
        if (a, n) not in largest:
            continue
        for g in b7a.PARTIAL_FRACTIONS:
            ratios: dict[str, list[float]] = {"build_per_pair": [], "scan_per_summand": [],
                                              "descent_per_probe": [], "keys_over_g2_full": []}
            for row in rows:
                full = f1_reports(row)
                p = b7a.figure(bround.runs() / "partial" / f"g{g}" / f"{row['id']}.price.json")
                part = bench.load(p)
                if not full or not (bench.clean(p) and part.get("status") == "extrapolated"):
                    continue
                fic, pic = full[0][1]["phases"]["ic"], part["phases"]["ic"]

                def per_pair(ic: dict) -> float:
                    return ic["build"]["ns"]["build"] / ic["build"]["stored_pairs"]

                def descent(ic: dict) -> float:
                    s = ic["descent"]["sample"]
                    return statistics.median(s.get("ns_per_probe") or s.get("ns_per_summand"))

                ratios["build_per_pair"].append(per_pair(pic) / per_pair(fic))
                ratios["scan_per_summand"].append(pic["collect"]["sample"]["ns_per_summand"]
                                                  / fic["collect"]["sample"]["ns_per_summand"])
                ratios["descent_per_probe"].append(descent(pic) / descent(fic))
                fraction = pic["build"]["partial_table"]["fraction"]
                ratios["keys_over_g2_full"].append(pic["build"]["distinct_keys"]
                                                   / (fraction ** 2 * fic["build"]["distinct_keys"]))
            out.append({"slug": bench.curve_slug(a, n), "n": n, "fraction_asked": g,
                        **{k: med(v) for k, v in ratios.items()}, "rows": len(ratios["build_per_pair"])})
    return out


def analyse(steps: str) -> dict:
    out = bround.analyse(steps)
    sizes = sorted(bround.by_size(bround.m1()).items(), key=lambda kv: kv[1][0]["r"])
    m5 = measurement5(sizes)
    out["f1_against_f0"] = {
        "what": "per size: F0's cold units (measurement 4's B7a arm, median of its processes) against F1's extrapolated cold units (median of its processes' figures)",
        "accounting": bround.accounting(bround.runs() / "f1") if (bround.runs() / "f1").exists() else None,
        "sizes": m5,
    }
    out["count_check"] = measurement6(sizes)
    out["partial_tables"] = measurement7(sizes) if (bround.runs() / "partial").exists() else None
    inside = sum(1 for s in m5 if s["f0_inside_predictive"])
    tight = [s for s in m5 if s["n"] in TIGHT_SIZES]
    target = {
        "inside_predictive": inside, "inside_needed": INSIDE_NEEDED, "sizes": len(m5),
        "tight_band": TIGHT_BAND,
        "tight": {s["slug"]: s["ratio_f1_over_f0"] for s in tight},
        "met": inside >= INSIDE_NEEDED and len(tight) == len(TIGHT_SIZES) and all(s["in_tight_band"] for s in tight),
    }
    out["falsification_target"] = target
    reasons = []
    incomplete = [s["slug"] for s in m5 if not s["f0_processes"] or not s["f1_processes"]]
    if incomplete:
        out["decision"] = {"accepted": False, "complete": False,
                           "reasons": [f"measurement 5 has no F0 or no F1 figure at {incomplete}"]}
        return out
    conf = out.get("conformance", {}).get("candidate", {})
    if conf.get("failed"):
        reasons.append(f"conformance cases failed: {conf['failed']}")
    for name in ("pin", "translate"):
        if name in out and not out[name]["held"]:
            reasons.append(f"the {name} check failed on {out[name]['failing_rows']}")
    if out.get("timing", {}).get("any_regression"):
        reasons.append("a size regresses beyond its A/A band at F0")
    if not target["met"]:
        reasons.append("the falsification target is missed")
    out["decision"] = {"accepted": not reasons, "complete": True, "reasons": reasons,
                       "note": "the tests (measurement 1) are run with cargo and recorded beside this analysis"}
    return out

#!/usr/bin/env python3
"""R05's figures and decision (PROTOCOL.md, "Success and stop"), from
runs/ only.

    tar -xJf runs.tar.xz && python3 analyse.py > analysis.json

Every ratio is the base over the candidate, so a ratio above 1 means the
candidate is faster. The tests run with cargo before the timed steps and
are recorded beside this analysis.
"""
from __future__ import annotations

import json
import math
import os
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "harness"))
sys.path.insert(0, str(HERE))
import bench  # noqa: E402
import run as r05  # noqa: E402
import stats  # noqa: E402

RUNS = (HERE / os.environ.get("IC_RUNS", "runs")).resolve()
R01 = HERE.parent / "R01-baseline-v0"
ACCEPT_LO = 1.10


def load(path: Path) -> dict | None:
    return json.loads(path.read_text()) if path.exists() else None


def figure(path: Path) -> dict | None:
    """The row's figure: its first clean, complete attempt, or None."""
    p = bench.figure_path(path)
    rep = load(p)
    return rep if rep is not None and rep.get("status") == "complete" and bench.clean(p) else None


def rounds_of(d: Path, arm: str, row_id: str) -> list[int]:
    return sorted(int(p.name[1:].split(".")[0]) for p in (d / arm / row_id).glob("r*.price.json")
                  if "-retry" not in p.name)


def paired(d: Path, rows: list[dict], measure) -> dict:
    ratios, missing = [], 0
    for r in rows:
        for k in rounds_of(d, "base", r["id"]):
            a = figure(d / "base" / r["id"] / f"r{k}.price.json")
            b = figure(d / "cand" / r["id"] / f"r{k}.price.json")
            if a is None or b is None:
                missing += 1
                continue
            ratios.append(measure(a) / measure(b))
    out = stats.geo_ci(ratios)
    out["missing_pairs"] = missing
    if "hi" in out:
        out["half_width"] = out["hi"] / out["geomean"] - 1
    return out


def accounting(d: Path) -> dict:
    counts = {"processes": 0, "contended": 0, "failed": 0}
    for rec in d.rglob("*.isolation.jsonl"):
        for line in rec.read_text().splitlines():
            run = json.loads(line)["run"]
            counts["processes"] += 1
            counts["contended"] += bool(run["contended"])
            counts["failed"] += run["exit_status"] != 0
    return counts


def by_size(rows: list[dict]) -> dict[tuple[int, int], list[dict]]:
    out: dict[tuple[int, int], list[dict]] = {}
    for r in rows:
        out.setdefault((r["a"], r["n"]), []).append(r)
    return out


def collect(rep: dict) -> float:
    return stats.setup_phase_ns(rep)["collect"]


def per_summand(rep: dict) -> float:
    return collect(rep) / rep["counts"]["pass"]["collect"]["summands_scanned"]


def r01_aa() -> dict[str, dict]:
    doc = json.loads((R01 / "analysis.json").read_text())
    return {bench.curve_slug(int(row["size"][1]), int(row["size"][3:])): row["cold"] for row in doc["aa"]}


def unit_v0() -> dict[str, float]:
    doc = json.loads((HERE.parents[1] / "baselines.json").read_text())
    v0 = next(b for b in doc["baselines"] if b["baseline"] == "v0")
    return {s["curve_alias"]: s["unit_ns"] for s in v0["sizes"]}


def s_cold(d: Path, arm: str, rows: list[dict], unit_ns: float) -> float | None:
    """S cold in v0's unit: per row the median over its rounds, then the
    median over the size's rows."""
    per_row = []
    for r in rows:
        cold = [stats.ic_cold_ns(rep) for k in rounds_of(d, arm, r["id"])
                if (rep := figure(d / arm / r["id"] / f"r{k}.price.json")) is not None]
        if cold:
            per_row.append(statistics.median(cold))
    return statistics.median(per_row) / unit_ns / math.sqrt(rows[0]["r"]) if per_row else None


def comparison(d: Path, rows: list[dict], aa: dict, units: dict) -> list[dict]:
    out = []
    for (a, n), rs in sorted(by_size(rows).items(), key=lambda kv: kv[1][0]["r"]):
        slug = bench.curve_slug(a, n)
        row = {"slug": slug, "a": a, "n": n, "log2_r": round(math.log2(rs[0]["r"]), 3), "rows": len(rs),
               "rounds": max((len(rounds_of(d, "base", r["id"])) for r in rs), default=0),
               "cold": paired(d, rs, stats.ic_cold_ns),
               "online": paired(d, rs, stats.ic_online_ns),
               "collect_stage": paired(d, rs, collect),
               "rho_online": paired(d, rs, stats.rho_online_ns)}
        if slug in units:
            unit = units[slug]
            row["s_cold_base"] = s_cold(d, "base", rs, unit)
            row["s_cold_cand"] = s_cold(d, "cand", rs, unit)
            reps = {arm: [rep for r in rs for k in rounds_of(d, arm, r["id"])
                          if (rep := figure(d / arm / r["id"] / f"r{k}.price.json")) is not None]
                    for arm in ("base", "cand")}
            if reps["base"] and reps["cand"]:
                row["collect_units_per_summand"] = {
                    arm: statistics.median(per_summand(x) for x in reps[arm]) / unit for arm in reps}
        band = aa.get(slug)
        if band and "lo" in band and "hi" in row["cold"]:
            row["aa_cold"] = {k: band[k] for k in ("geomean", "lo", "hi")}
            row["regresses_beyond_aa"] = row["cold"]["hi"] < band["lo"]
        out.append(row)
    return out


def decision(pin: dict | None, suite: list[dict], holdouts: list[dict]) -> dict:
    reasons = []
    if not (pin or {}).get("held"):
        reasons.append("an output differs from v0's")
    if not (pin or {}).get("names_agree"):
        reasons.append("a name disagrees with the registry")
    for a, n in r05.TARGETS:
        slug = bench.curve_slug(a, n)
        for name, rows in (("suite", suite), ("holdouts", holdouts)):
            row = next((r for r in rows if r["slug"] == slug), None)
            lo = (row or {}).get("cold", {}).get("lo")
            if lo is None or lo <= ACCEPT_LO:
                reasons.append(f"{slug} on the {name}: interval's lower end {lo} is not above {ACCEPT_LO}")
    for r in suite:
        if r.get("regresses_beyond_aa"):
            reasons.append(f"{r['slug']} regresses beyond its A/A band")
    return {"accepted": not reasons, "reasons": reasons,
            "note": "the tests run with cargo before the timed steps and are recorded beside this analysis"}


def main() -> None:
    aa, units = r01_aa(), unit_v0()
    pin = load(RUNS / "pin" / "pin.json")
    suite = comparison(RUNS / "compare", r05.suite_rows(), aa, units)
    holdouts = comparison(RUNS / "holdout", r05.holdout_rows(), {}, units)
    doc = {
        "what_this_is": "R05, the sharper presence filter: every figure the README, the ledger and the "
                        "scoreboard quote",
        "host": load(RUNS / "host.json"),
        "aa_source": load(RUNS / "aa-source.json"),
        "accounting": {step: accounting(RUNS / step) for step in ("compare", "holdout") if (RUNS / step).exists()},
        "pin": {k: pin[k] for k in ("held", "names_agree")} | {"rows": len(pin["rows"])} if pin else None,
        "extended": load(RUNS / "extended.json"),
        "suite": suite,
        "holdouts": holdouts,
    }
    doc["decision"] = decision(pin, suite, holdouts)
    print(json.dumps(doc, indent=1))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""R03's figures and decision (PROTOCOL.md, "Success and stop"), from
runs/ only.

    tar -xJf runs.tar.xz && python3 analyse.py > analysis.json

Every ratio is the base over the candidate, so a ratio above 1 means the
candidate is faster.
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
import bench  # noqa: E402
import stats  # noqa: E402

RUNS = (HERE / os.environ.get("IC_RUNS", "runs")).resolve()
R01 = HERE.parent / "R01-baseline-v0"
TARGET = (0, 57)
ACCEPT_LO = 1.3
COMPOSITE = [(1, 45), (0, 57)]


def load(path: Path) -> dict | None:
    return json.loads(path.read_text()) if path.exists() else None


def figure(path: Path) -> dict | None:
    p = bench.figure_path(path)
    rep = load(p)
    return rep if rep is not None and rep.get("status") == "complete" and bench.clean(p) else None


def rounds_of(d: Path, arm: str, row_id: str) -> list[int]:
    return sorted(int(p.name[1:].split(".")[0]) for p in (d / arm / row_id).glob("r*.price.json")
                  if "-retry" not in p.name)


def pairs(d: Path, rows: list[dict]):
    for r in rows:
        for k in rounds_of(d, "base", r["id"]):
            yield (figure(d / "base" / r["id"] / f"r{k}.price.json"),
                   figure(d / "cand" / r["id"] / f"r{k}.price.json"))


def paired(d: Path, rows: list[dict], measure) -> dict:
    ratios, missing = [], 0
    for a, b in pairs(d, rows):
        if a is None or b is None:
            missing += 1
            continue
        ratios.append(measure(a) / measure(b))
    out = stats.geo_ci(ratios)
    out["missing_pairs"] = missing
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


def construction(rep: dict) -> float:
    """The set-up phase that builds the curve (`setup` in the pricer's split)."""
    return stats.setup_phase_ns(rep)["setup"]


def r01_aa() -> dict[str, dict]:
    doc = json.loads((R01 / "analysis.json").read_text())
    return {bench.curve_slug(int(row["size"][1]), int(row["size"][3:])): row["cold"] for row in doc["aa"]}


def unit_v0() -> dict[str, float]:
    doc = json.loads((HERE.parents[1] / "baselines.json").read_text())
    v0 = next(b for b in doc["baselines"] if b["baseline"] == "v0")
    return {s["curve_alias"]: s["unit_ns"] for s in v0["sizes"]}


def s_cold(d: Path, arm: str, rows: list[dict], unit_ns: float) -> float | None:
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
               "composite": (a, n) in COMPOSITE,
               "cold": paired(d, rs, stats.ic_cold_ns),
               "setup": paired(d, rs, lambda rep: stats.median_rep(rep, "setup_ns")),
               "construction_stage": paired(d, rs, construction)}
        base_c = [construction(a) for a, _ in pairs(d, rs) if a is not None]
        cand_c = [construction(b) for _, b in pairs(d, rs) if b is not None]
        if base_c and cand_c:
            row["construction_ms_median"] = {"base": statistics.median(base_c) / 1e6,
                                             "cand": statistics.median(cand_c) / 1e6}
        if slug in units:
            row["s_cold_base"] = s_cold(d, "base", rs, units[slug])
            row["s_cold_cand"] = s_cold(d, "cand", rs, units[slug])
        band = aa.get(slug)
        if band and "lo" in band and "hi" in row["cold"]:
            row["aa_cold"] = {k: band[k] for k in ("geomean", "lo", "hi")}
            row["regresses_beyond_aa"] = row["cold"]["hi"] < band["lo"]
        out.append(row)
    return out


def decision(pin, suite, holdout) -> dict:
    reasons = []
    if not (pin or {}).get("held"):
        reasons.append("an output differs from v0's")
    slug = bench.curve_slug(*TARGET)
    for name, rows in (("suite", suite), ("holdouts", holdout)):
        row = next((r for r in rows if r["slug"] == slug), None)
        lo = (row or {}).get("cold", {}).get("lo")
        if lo is None or lo <= ACCEPT_LO:
            reasons.append(f"{slug} on the {name}: interval's lower end {lo} is not above {ACCEPT_LO}")
    for r in suite:
        if r.get("regresses_beyond_aa"):
            reasons.append(f"{r['slug']} regresses beyond its A/A band")
    return {"accepted": not reasons, "reasons": reasons}


def main() -> None:
    aa, units = r01_aa(), unit_v0()
    rows = [r for r in bench.slug_rows(bench.suite_rows("S"))
            if (r["a"], r["n"]) in COMPOSITE or r["recipe_seed"] == 201]
    holdout_rows = []
    for a, n in COMPOSITE:
        slug = bench.curve_slug(a, n)
        r0 = next(r for r in rows if (r["a"], r["n"]) == (a, n))
        holdout_rows += [{"id": f"{slug}/M5-T{i}", "a": a, "n": n, "r": r0["r"]} for i in (101, 102)]
    pin = load(RUNS / "pin" / "pin.json")
    suite = comparison(RUNS / "compare", rows, aa, units)
    held = comparison(RUNS / "holdout", holdout_rows, {}, units)
    doc = {
        "what_this_is": "R03, the curve's construction at composite degrees: every figure the README, the "
                        "ledger and the scoreboard quote",
        "host": load(RUNS / "host.json"),
        "aa_source": load(RUNS / "aa-source.json"),
        "accounting": {step: accounting(RUNS / step) for step in ("aa", "compare", "holdout")
                       if (RUNS / step).exists()},
        "pin": {k: pin[k] for k in ("held", "names_agree")} | {"rows": len(pin["rows"])} if pin else None,
        "suite": suite,
        "holdouts": held,
    }
    doc["decision"] = decision(pin, suite, held)
    print(json.dumps(doc, indent=1))


if __name__ == "__main__":
    main()

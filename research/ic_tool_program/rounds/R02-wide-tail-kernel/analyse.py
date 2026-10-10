#!/usr/bin/env python3
"""R02's figures and decision (PROTOCOL.md, "Success and stop", with
amendments 1–3), from runs/ only.

    tar -xJf runs.tar.xz && python3 analyse.py > analysis.json

Every ratio is v0′ over the candidate, or v0 over v0′ for the
accounting check, so a ratio above 1 means the second arm is faster.
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
TARGETS = [(1, 59), (0, 61)]
ACCEPT_LO = 1.10
CONTROL_THRESHOLD = 1.3


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


def paired(d: Path, arms: tuple[str, str], rows: list[dict], measure) -> dict:
    """The geometric mean of `measure(a) / measure(b)` over every row and
    every round both arms completed, with its t interval."""
    ratios, missing = [], 0
    for r in rows:
        for k in rounds_of(d, arms[0], r["id"]):
            a = figure(d / arms[0] / r["id"] / f"r{k}.price.json")
            b = figure(d / arms[1] / r["id"] / f"r{k}.price.json")
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


def build(rep: dict) -> float:
    return stats.setup_phase_ns(rep)["build"]


def per_summand(rep: dict) -> float:
    return collect(rep) / rep["counts"]["pass"]["collect"]["summands_scanned"]


def r01_aa() -> dict[str, dict]:
    """R01's A/A cold interval per size, keyed by slug."""
    doc = json.loads((R01 / "analysis.json").read_text())
    out = {}
    for row in doc["aa"]:
        a, n = int(row["size"][1]), int(row["size"][3:])
        out[bench.curve_slug(a, n)] = row["cold"]
    return out


def unit_v0() -> dict[str, float]:
    """v0's unit per size (baselines.json), keyed by slug."""
    doc = json.loads((HERE.parents[1] / "baselines.json").read_text())
    v0 = next(b for b in doc["baselines"] if b["baseline"] == "v0")
    return {s["curve_alias"]: s["unit_ns"] for s in v0["sizes"]}


def s_cold(d: Path, arm: str, rows: list[dict], unit_ns: float) -> float | None:
    """S cold in v0's unit: per row, the median over its rounds; then the
    median over the size's rows."""
    per_row = []
    for r in rows:
        cold = [stats.ic_cold_ns(rep) for k in rounds_of(d, arm, r["id"])
                if (rep := figure(d / arm / r["id"] / f"r{k}.price.json")) is not None]
        if cold:
            per_row.append(statistics.median(cold))
    if not per_row:
        return None
    return statistics.median(per_row) / unit_ns / math.sqrt(rows[0]["r"])


def comparison(d: Path, rows: list[dict], aa: dict, units: dict) -> list[dict]:
    out = []
    for (a, n), rs in sorted(by_size(rows).items(), key=lambda kv: kv[1][0]["r"]):
        slug = bench.curve_slug(a, n)
        row = {"slug": slug, "a": a, "n": n, "log2_r": round(math.log2(rs[0]["r"]), 3),
               "rows": len(rs), "rounds": max((len(rounds_of(d, "base", r["id"])) for r in rs), default=0),
               "cold": paired(d, ("base", "cand"), rs, stats.ic_cold_ns),
               "setup": paired(d, ("base", "cand"), rs, lambda rep: stats.median_rep(rep, "setup_ns")),
               "online": paired(d, ("base", "cand"), rs, stats.ic_online_ns),
               "collect_stage": paired(d, ("base", "cand"), rs, collect),
               "build_stage": paired(d, ("base", "cand"), rs, build),
               "rho_online": paired(d, ("base", "cand"), rs, stats.rho_online_ns)}
        if slug in units:
            unit = units[slug]
            row["s_cold_base"] = s_cold(d, "base", rs, unit)
            row["s_cold_cand"] = s_cold(d, "cand", rs, unit)
            reps = [rep for r in rs for k in rounds_of(d, "cand", r["id"])
                    if (rep := figure(d / "cand" / r["id"] / f"r{k}.price.json")) is not None]
            base_reps = [rep for r in rs for k in rounds_of(d, "base", r["id"])
                         if (rep := figure(d / "base" / r["id"] / f"r{k}.price.json")) is not None]
            if reps and base_reps:
                row["collect_units_per_summand"] = {
                    "base": statistics.median(per_summand(x) for x in base_reps) / unit,
                    "cand": statistics.median(per_summand(x) for x in reps) / unit}
        band = aa.get(slug)
        if band and "lo" in band and "hi" in row["cold"]:
            row["aa_cold"] = {k: band[k] for k in ("geomean", "lo", "hi")}
            row["regresses_beyond_aa"] = row["cold"]["hi"] < band["lo"]
        out.append(row)
    return out


def control() -> dict | None:
    d = RUNS / "control"
    if not d.exists():
        return None
    verdict = load(d / "verdict.json")
    return {"verdict": verdict,
            "stops_r02": verdict is not None
            and verdict["collect_ns_per_summand_scalar_over_simd"].get("geomean", 0) < CONTROL_THRESHOLD,
            "reading": "both AVX-512 scan kernels off together (amendment 3): a ratio below the threshold "
                       "stops R02; one above it does not confirm the addition's share"}


def callgrind() -> dict:
    d = RUNS / "callgrind"
    out = {}
    for f in sorted(d.glob("*.phases.json")) if d.exists() else []:
        tag = f.name[: -len(".phases.json")]
        doc = json.loads(f.read_text() or "{}")
        wf = load(d / f"{tag}.workflow.json") or {}
        out[tag] = {"total_ir": doc.get("total", {}).get("Ir"), "status": wf.get("status"),
                    "recovered": [s.get("recovered") for s in wf.get("solutions", {}).get("items", [])]}
    pairs = {}
    for tag, v in out.items():
        if tag.startswith("base-"):
            other = out.get("cand-" + tag[len("base-"):])
            if other and v["total_ir"] and other["total_ir"]:
                pairs[tag[len("base-"):]] = {"ir_base_over_cand": v["total_ir"] / other["total_ir"],
                                              "same_log": bool(v["recovered"]) and v["recovered"] == other["recovered"]
                                              and v["status"] == other["status"] == "complete"}
    return {"runs": out, "pairs": pairs,
            "control_held": bool(pairs) and all(abs(p["ir_base_over_cand"] - 1) < 1e-3 and p["same_log"]
                                                for p in pairs.values())}


def decision(ctrl, pin, suite, holdout) -> dict:
    reasons = []
    if ctrl is None or ctrl["stops_r02"]:
        reasons.append("the control did not clear 1.3")
    if not (pin or {}).get("held"):
        reasons.append("an output differs from v0's")
    for a, n in TARGETS:
        slug = bench.curve_slug(a, n)
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
    rows = bench.slug_rows(bench.suite_rows("S"))
    holdout_rows = []
    for (a, n), rs in by_size(rows).items():
        slug = bench.curve_slug(a, n)
        holdout_rows += [{"id": f"{slug}/M5-T{i}", "a": a, "n": n, "r": rs[0]["r"]} for i in (101, 102)]
    ctrl = control()
    pin = load(RUNS / "pin" / "pin.json")
    suite = comparison(RUNS / "compare", rows, aa, units)
    held = comparison(RUNS / "holdout", holdout_rows, {}, units)
    v0check = [{"slug": bench.curve_slug(a, n), "log2_r": round(math.log2(rs[0]["r"]), 3),
                "cold_v0_over_base": paired(RUNS / "v0check", ("v0", "base"), rs, stats.ic_cold_ns)}
               for (a, n), rs in sorted(by_size([r for r in rows if r["recipe_seed"] == 201]).items(),
                                        key=lambda kv: kv[1][0]["r"])]
    doc = {
        "what_this_is": "R02, the AVX-512 batched addition for wide-tail fields: every figure the README, "
                        "the ledger and the scoreboard quote",
        "host": load(RUNS / "host.json"),
        "aa_source": load(RUNS / "aa-source.json"),
        "accounting": {step: accounting(RUNS / step)
                       for step in ("control", "aa", "compare", "holdout", "v0check") if (RUNS / step).exists()},
        "control": ctrl,
        "pin": {k: pin[k] for k in ("held", "names_agree")} | {"rows": len(pin["rows"])} if pin else None,
        "extended_sizes": (load(RUNS / "compare" / "extended.json") or {}).get("extended_sizes"),
        "suite": suite,
        "holdouts": held,
        "v0_over_base": v0check,
        "callgrind": callgrind(),
    }
    doc["decision"] = decision(ctrl, pin, suite, held)
    print(json.dumps(doc, indent=1))


if __name__ == "__main__":
    main()

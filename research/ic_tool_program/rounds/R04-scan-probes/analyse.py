#!/usr/bin/env python3
"""R04's figures (PROTOCOL.md, "Measurements" and "What it decides"),
from runs/ only.

    tar -xJf runs.tar.xz && python3 analyse.py > analysis.json

The shares come from the probe arm's `scan_probes` totals. Each process
sums its in-process repetitions, so a share is a ratio of totals, and the
cycles per summand divide by the summands those repetitions scanned. The
overhead is the default arm's cold time over the probe arm's, so a ratio
below 1 means the probes cost time.
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
STAGES = ("subtract", "key", "filter", "admitted", "trial")
TOP = 4  # the top four sizes decide


def load(path: Path) -> dict | None:
    return json.loads(path.read_text()) if path.exists() else None


def figure(path: Path) -> dict | None:
    p = bench.figure_path(path)
    rep = load(p)
    return rep if rep is not None and rep.get("status") == "complete" and bench.clean(p) else None


def rounds_of(d: Path, arm: str, row_id: str) -> list[int]:
    return sorted(int(p.name[1:].split(".")[0]) for p in (d / arm / row_id).glob("r*.price.json")
                  if "-retry" not in p.name)


def accounting(d: Path) -> dict:
    counts = {"processes": 0, "contended": 0, "failed": 0}
    for rec in d.rglob("*.isolation.jsonl"):
        for line in rec.read_text().splitlines():
            run = json.loads(line)["run"]
            counts["processes"] += 1
            counts["contended"] += bool(run["contended"])
            counts["failed"] += run["exit_status"] != 0
    return counts


def r01_aa() -> dict[str, dict]:
    doc = json.loads((R01 / "analysis.json").read_text())
    return {bench.curve_slug(int(row["size"][1]), int(row["size"][3:])): row["cold"] for row in doc["aa"]}


def by_size(rows: list[dict]) -> dict[tuple[int, int], list[dict]]:
    out: dict[tuple[int, int], list[dict]] = {}
    for r in rows:
        out.setdefault((r["a"], r["n"]), []).append(r)
    return out


def shares(probes: list[dict]) -> dict:
    """Per stage: the median over processes of the share of the scan and of
    the nanoseconds per scanned summand."""
    per = {s: {"share": [], "ns_per_summand": []} for s in STAGES}
    scanned, admitted = [], []
    for rep in probes:
        p = rep["scan_probes"]
        cycles = p["stages"]
        scan = sum(cycles[s] for s in STAGES if s != "trial")
        summands = p["counts"]["summands"]
        if not scan or not summands:
            continue
        ns_per_cycle = 1e9 / p["tsc_hz"]
        for s in STAGES:
            per[s]["share"].append(cycles[s] / scan)
            per[s]["ns_per_summand"].append(cycles[s] * ns_per_cycle / summands)
        scanned.append(scan * ns_per_cycle / summands)
        admitted.append(p["counts"]["admitted"] / summands)
    out = {s: {k: statistics.median(v) if v else None for k, v in per[s].items()} for s in STAGES}
    out["scan_ns_per_summand"] = statistics.median(scanned) if scanned else None
    out["admitted_per_summand"] = statistics.median(admitted) if admitted else None
    out["processes"] = len(scanned)
    return out


def analyse_size(d: Path, rs: list[dict], band: dict | None) -> dict:
    ratios, probes, missing = [], [], 0
    for r in rs:
        for k in rounds_of(d, "default", r["id"]):
            a = figure(d / "default" / r["id"] / f"r{k}.price.json")
            b = figure(d / "probes" / r["id"] / f"r{k}.price.json")
            if a is None or b is None:
                missing += 1
                continue
            ratios.append(stats.ic_cold_ns(a) / stats.ic_cold_ns(b))
            if "scan_probes" in b:
                probes.append(b)
    overhead = stats.geo_ci(ratios)
    overhead["missing_pairs"] = missing
    out = {"overhead_default_over_probes": overhead, "stages": shares(probes)}
    if band and "lo" in band and "lo" in overhead:
        out["aa_cold"] = {k: band[k] for k in ("geomean", "lo", "hi")}
        out["overhead_inside_aa"] = band["lo"] <= overhead["lo"] and overhead["hi"] <= band["hi"]
    return out


def main() -> None:
    rows = [r for r in bench.slug_rows(bench.suite_rows("S")) if r["recipe_seed"] == 201]
    aa = r01_aa()
    sizes = []
    for (a, n), rs in sorted(by_size(rows).items(), key=lambda kv: kv[1][0]["r"]):
        slug = bench.curve_slug(a, n)
        row = {"slug": slug, "a": a, "n": n, "log2_r": round(math.log2(rs[0]["r"]), 3)}
        row.update(analyse_size(RUNS / "compare", rs, aa.get(slug)))
        sizes.append(row)
    top = sizes[-TOP:]
    leading = {}
    for row in top:
        st = row["stages"]
        ranked = sorted((s for s in STAGES if s != "trial" and st[s]["share"] is not None),
                        key=lambda s: st[s]["share"], reverse=True)
        leading[row["slug"]] = ranked[0] if ranked else None
    votes = [s for s in leading.values() if s]
    pin = load(RUNS / "pin" / "pin.json")
    doc = {
        "what_this_is": "R04, the scan's stages priced inside the scan: a stage diagnostic, never a speedup",
        "host": load(RUNS / "host.json"),
        "aa_source": load(RUNS / "aa-source.json"),
        "accounting": {step: accounting(RUNS / step) for step in ("aa", "compare") if (RUNS / step).exists()},
        "pin": {k: pin[k] for k in ("held", "probes_reported", "default_silent")} | {"rows": len(pin["rows"])}
        if pin else None,
        "sizes": sizes,
        "decision": {
            "leading_stage_by_top_size": leading,
            "leading_stage": max(set(votes), key=votes.count) if votes else None,
            "trusted_at_top_sizes": {row["slug"]: row.get("overhead_inside_aa") for row in top},
        },
    }
    print(json.dumps(doc, indent=1))


if __name__ == "__main__":
    main()

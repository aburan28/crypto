#!/usr/bin/env python3
"""Instruction-count S for the end-to-end index-calculus rungs.

The end-to-end gate (scripts/ic_e2e_benchmark.py) pins counters and gates a
same-host wall ratio, but leaves the index-calculus S null because its
counters are in mixed units.  This measures both arms in one hardware-
independent unit the repository already uses for the beats-rho thread:
**valgrind callgrind instructions (Ir)**, S = Ir / (targets * sqrt(r)).

Per rung, ONE callgrind pass of the full workflow (pipeline + in-process rho
on the same targets), single-threaded (RAYON_NUM_THREADS=1) so that rayon
workers do not spin under callgrind's serialised threading and the count is
work, not scheduling.  Ir_rho is the inclusive cost of the rho routine
``koblitz_signed_frobenius_rho_with_progress`` (all calls, setup included)
read from ``callgrind_annotate --inclusive=yes``; Ir_ic = total - Ir_rho, so
process start-up, curve construction, target resolution and JSON output are
charged to the index-calculus side, the conservative direction for the claim
under test.  Single-threaded Ir is reproducible to ~2e-5 relative.

The wall time of a callgrind run is meaningless and is not recorded as a
cost; the report's own timing fields are ignored here.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import math
import os
import platform
import re
import subprocess
import sys
from pathlib import Path

SCHEMA_VERSION = 1


def parse_callgrind_total(path: Path) -> int:
    """Total Ir from a callgrind output file (its 'summary:' or 'totals:' line)."""
    events = None
    total = None
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("events:"):
            events = line.split()[1:]
        elif line.startswith("summary:") or line.startswith("totals:"):
            nums = [int(x) for x in line.split()[1:]]
            if events and "Ir" in events:
                total = nums[events.index("Ir")]
            else:
                total = nums[0]
    if total is None:
        raise ValueError(f"no summary/totals line in {path}")
    return total


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def valgrind_version() -> str:
    return subprocess.run(["valgrind", "--version"], capture_output=True, text=True).stdout.strip()


RHO_SYMBOL = "koblitz_signed_frobenius_rho_with_progress"
# Inclusive cost of the largest matching symbol, across threads (rayon keeps
# one worker even at RAYON_NUM_THREADS=1 and hands it the table build).  A
# diagnostic split of Ir_ic, not a second accounting: the two rows may
# overlap or leave a remainder.
BREAKDOWN_SYMBOLS = {
    "pair_table_build": "PairSumTable>::build",
    "decomposition_scan": "PairSumTable>::decompose_fast",
    "collector_with_pair": "RelationCollector>::with_pair",
}


def inclusive_costs(callgrind_out: Path) -> list[tuple[int, str]]:
    """(inclusive Ir, symbol) rows from callgrind_annotate --inclusive=yes."""
    text = subprocess.run(["callgrind_annotate", "--inclusive=yes", "--threshold=100", str(callgrind_out)],
                          capture_output=True, text=True).stdout
    rows = []
    for line in text.splitlines():
        m = re.match(r"\s*([\d,]+)\s+\([\d. ]+%\)\s+\S+?:(.*)$", line)
        if m and "PROGRAM TOTALS" not in line:
            rows.append((int(m.group(1).replace(",", "")), m.group(2).strip()))
    return rows


def inclusive_of(rows: list[tuple[int, str]], needle: str) -> int | None:
    hits = [ir for ir, sym in rows if needle in sym]
    return max(hits) if hits else None


def run_callgrind(ic: Path, params: Path, run_dir: Path, report: Path, out: Path) -> dict:
    """One single-threaded callgrind pass of the full workflow (pipeline + rho)."""
    cmd = ["valgrind", "--tool=callgrind", f"--callgrind-out-file={out}", str(ic), "workflow",
           "--params", str(params), "--dir", str(run_dir), "--out", str(report), "--json"]
    env = dict(os.environ, RAYON_NUM_THREADS="1")
    started = _dt.datetime.now(_dt.timezone.utc)
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if proc.returncode != 0 or not report.is_file():
        raise RuntimeError(f"{' '.join(cmd)} failed ({proc.returncode}): {proc.stderr[-1500:]}")
    rows = inclusive_costs(out)
    rho = inclusive_of(rows, RHO_SYMBOL)
    if rho is None:
        raise RuntimeError(f"{out}: rho symbol {RHO_SYMBOL} not found; is the binary stripped?")
    return {"command": cmd, "env": {"RAYON_NUM_THREADS": "1"}, "ir_total": parse_callgrind_total(out), "ir_rho_inclusive": rho,
            "breakdown_inclusive": {k: inclusive_of(rows, v) for k, v in BREAKDOWN_SYMBOLS.items()},
            "callgrind_out": str(out),
            "wall_seconds_under_callgrind_not_a_cost": (_dt.datetime.now(_dt.timezone.utc) - started).total_seconds()}


def measure_rung(ic: Path, params: Path, work: Path) -> dict:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import ic_e2e_benchmark as bench  # noqa: E402
    name = params.stem
    d = work / name
    d.mkdir(parents=True)
    run = run_callgrind(ic, params, d / "run", d / "report.json", d / "callgrind.out")
    rep = json.loads((d / "report.json").read_text())
    m = bench.measure(rep)
    r = int(m["subgroup_order"])
    sqrt_r = math.sqrt(r)
    ir_rho = run["ir_rho_inclusive"]
    ir_ic = run["ir_total"] - ir_rho
    return {
        "params": str(params), "name": rep["name"], "degree": m["degree"], "subgroup_order": m["subgroup_order"],
        "subgroup_bits": m["subgroup_bits"], "targets": m["targets"], "ic_verified": m["ic_verified"], "rho_verified": m["rho_verified"],
        "counters": m["counters"],
        "ir_total": run["ir_total"], "ir_rho": ir_rho, "ir_ic": ir_ic,
        "S_ic": ir_ic / (m["targets"] * sqrt_r), "S_rho": ir_rho / (m["targets"] * sqrt_r),
        "rho_over_ic_ir": ir_rho / ir_ic,
        "rho_ir_per_group_addition": ir_rho / m["counters"]["rho_group_additions"],
        "rho_S_group_additions": m["rho_S"]["measured"],
        "breakdown_inclusive": run["breakdown_inclusive"],
        "run": run,
    }


def cmd_measure(args: argparse.Namespace) -> int:
    work = Path(args.output)
    if work.exists() and any(work.iterdir()):
        print(f"refusing to write into non-empty {work}", file=sys.stderr)
        return 2
    work.mkdir(parents=True, exist_ok=True)
    ic = Path(args.ic)
    result = {
        "schema_version": SCHEMA_VERSION, "unit": f"{valgrind_version()}-{platform.machine()}-Ir",
        "S_definition": "Ir / (targets * sqrt(r)) from ONE single-threaded (RAYON_NUM_THREADS=1) callgrind pass of the full workflow; Ir_rho = inclusive Ir of koblitz_signed_frobenius_rho_with_progress (all 32 calls, setup included), Ir_ic = total - Ir_rho, so process start-up, curve construction, target resolution and JSON output are charged to index calculus",
        "determinism": "single-threaded totals of two passes of the n=31 rung differed by 1.6e-5 relative and the rho inclusive by 1.2e-5; one pass is reported to that precision",
        "threads": "RAYON_NUM_THREADS=1 still leaves one rayon worker thread, to which the pair-table build is handed; PROGRAM TOTALS sums every thread, so Ir_total and Ir_ic are whole-process work; the rho routine runs on the main thread",
        "started_at": _dt.datetime.now(_dt.timezone.utc).isoformat(), "host": {"platform": platform.platform(), "machine": platform.machine()},
        "ic_binary": {"path": str(ic), "sha256": sha256(ic)},
        "git": subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(),
        "rungs": [],
    }
    for p in args.params:
        print(f"[ir] {Path(p).stem}: one single-threaded callgrind pass …", flush=True)
        rung = measure_rung(ic, Path(p), work)
        result["rungs"].append(rung)
        print(json.dumps({k: rung[k] for k in ("name", "ir_ic", "ir_rho", "S_ic", "S_rho", "rho_over_ic_ir", "breakdown_inclusive")}), flush=True)
        (work / "instructions.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n")
    result["finished_at"] = _dt.datetime.now(_dt.timezone.utc).isoformat()
    (work / "instructions.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n")
    return 0


def render_markdown(result: dict) -> str:
    lines = [f"Unit: `{result['unit']}`. `S = Ir / (targets · √r)` from one single-threaded callgrind pass; `Ir_rho` is the inclusive cost of the rho routine, `Ir_ic = total − Ir_rho` (start-up and output charged to IC).", "",
             "| rung | log₂ r | targets | Ir, index calculus | Ir, rho | S, index calculus | S, rho | rho / IC (Ir) | Ir per rho group addition |",
             "|:--|--:|--:|--:|--:|--:|--:|--:|--:|"]
    for r in result["rungs"]:
        lines.append("| `%s` | %.1f | %d | %s | %s | %.1f | %.1f | **%.3f** | %.0f |" % (
            r["name"], r["subgroup_bits"], r["targets"], f"{r['ir_ic']:,}", f"{r['ir_rho']:,}", r["S_ic"], r["S_rho"], r["rho_over_ic_ir"], r["rho_ir_per_group_addition"]))
    return "\n".join(lines) + "\n"


def cmd_report(args: argparse.Namespace) -> int:
    print(render_markdown(json.loads(Path(args.instructions).read_text())))
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)
    m = sub.add_parser("measure")
    m.add_argument("--ic", required=True)
    m.add_argument("--params", nargs="+", required=True)
    m.add_argument("--output", required=True, help="fresh directory; instructions.json is written here")
    m.set_defaults(func=cmd_measure)
    r = sub.add_parser("report")
    r.add_argument("--instructions", required=True)
    r.set_defaults(func=cmd_report)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())

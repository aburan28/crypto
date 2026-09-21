#!/usr/bin/env python3
"""E3, item 1: does a real solver's cost per call depend on the target?

Section 3.2 of `research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md`
localises a whole-base detector's witness by swapping: instead of querying
sub-bases, it asks the detector about `R - P + Q` for class-matched `Q`.  That
reduction is only worth anything if a detector charges the same for `R - P + Q`
as it charges for `R`.  The swap-localisation run measured the *query count*;
it never measured the *price of a query*, and the design listed the timing as
its first item.  This is that timing.

The lever is `--known-log`, which is the only thing that moves the descent
target while the curve, the factor base, the summand count and the seed stay
fixed.  Every solver in the repository's `ic run` is swept across the same
seven targets, and what is reported is solver work per unit of solver output,
because the totals also move with how many relations a run decided to collect
and that has nothing to do with the target.

This is a **stage diagnostic** in the sense of `AGENTS.md` section 8: it prices
one oracle call on one rung and infers nothing about a full discrete logarithm.

Usage:  python3 scripts/ecc2k130_e3_solver_panel.py [--degree 13] [--summands 3]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import statistics
import subprocess
import sys

REPO = Path(__file__).resolve().parents[1]
IC = REPO / "target/release/ic"
OUT = REPO / "experiments/ecc2k130_e3_solver_panel.json"

# Every solver `ic run` exposes.  `wdsat` needs an external binary that this
# repository does not vendor; it is swept anyway so the artefact records that
# it was attempted rather than skipped.
SOLVERS = ("groebner", "sat", "enumerate", "pair-table", "wdsat")

# Seven targets spread across the subgroup.  At degree 13 the subgroup has
# order 2003, so every one of these is a legal `--known-log`.
TARGETS = (53, 211, 499, 887, 1289, 1613, 1987)

# What counts as "one unit of solver output" for each solver, so that the
# ratio below is work per call rather than work per run.
UNIT = {
    "groebner": ("f4_word_ops", "f4_reductions", "F4 word operations per F4 reduction"),
    "sat": ("f4_word_ops", "independent_relations", "solver operations per independent relation"),
    "enumerate": (None, None, "cpu seconds for the whole fixed sweep"),
    "pair-table": (None, None, "cpu seconds for the whole fixed sweep"),
    "wdsat": (None, None, "not exercised: --solver wdsat requires --wdsat-binary"),
}


def one_run(solver: str, known_log: int, degree: int, summands: int) -> dict:
    cmd = [str(IC), "run", "--degree", str(degree), "--summands", str(summands),
           "--solver", solver, "--known-log", str(known_log), "--json"]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    try:
        d = json.loads(proc.stdout)
    except json.JSONDecodeError:
        return {"known_log": known_log, "status": "no_json",
                "stderr": (proc.stderr or "")[:200]}
    if d.get("status") != "complete":
        return {"known_log": known_log, "status": d.get("status", "?"),
                "message": (d.get("message") or "")[:200]}
    c, r = d.get("counts", {}), d.get("result", {})
    return {
        "known_log": known_log,
        "status": "complete",
        "verified": r.get("verified"),
        "f4_word_ops": c.get("f4_word_ops"),
        "f4_reductions": c.get("f4_reductions"),
        "independent_relations": c.get("independent_relations"),
        "cpu_seconds": d.get("resources", {}).get("cpu_seconds"),
    }


def spread(values) -> float:
    """Max over min, as a percentage above one.  The question is whether the
    per-call price moves with the target, so a peak-to-trough spread is the
    honest statistic: an average would hide exactly the swing being looked for."""
    lo, hi = min(values), max(values)
    return 100.0 * (hi / lo - 1.0) if lo > 0 else float("nan")


def panel(degree: int, summands: int) -> dict:
    out = {}
    for solver in SOLVERS:
        runs = [one_run(solver, k, degree, summands) for k in TARGETS]
        done = [r for r in runs if r["status"] == "complete"]
        num, den, unit = UNIT[solver]
        row = {"unit": unit, "runs": runs,
               "completed": len(done), "attempted": len(runs),
               "all_verified": bool(done) and all(r["verified"] for r in done)}
        if not done:
            row["exercised"] = False
            out[solver] = row
            continue
        row["exercised"] = True
        totals = [r[num] for r in done] if num else [r["cpu_seconds"] for r in done]
        row["total_spread_percent"] = round(spread(totals), 1)
        if num and den:
            per = [r[num] / r[den] for r in done if r[den]]
            row["per_call"] = [round(x, 3) for x in per]
            row["per_call_min"] = round(min(per), 3)
            row["per_call_max"] = round(max(per), 3)
            row["per_call_mean"] = round(statistics.fmean(per), 3)
            row["per_call_spread_percent"] = round(spread(per), 1)
        else:
            row["per_call_spread_percent"] = row["total_spread_percent"]
        out[solver] = row
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--degree", type=int, default=13)
    ap.add_argument("--summands", type=int, default=3)
    args = ap.parse_args()

    if not IC.exists():
        sys.exit(f"{IC.relative_to(REPO)} is not built; run `cargo build --release --bin ic`")

    solvers = panel(args.degree, args.summands)
    for name, row in solvers.items():
        if not row["exercised"]:
            print(f"{name:11s} not exercised ({row['unit']})")
            continue
        print(f"{name:11s} {row['unit']}")
        print(f"            per call spread {row['per_call_spread_percent']:5.1f}%"
              f"   whole-run total spread {row['total_spread_percent']:6.1f}%")

    report = {
        "schema": "ecc2k130_e3_solver_panel/v1",
        "question": "E3 item 1 -- is a real solver's cost per call a function "
                    "of the descent target?",
        "premise_under_test": "section 3.2's swap localisation queries "
                              "detector(R - P + Q); it is sound only if a "
                              "detector charges the same for that as for R",
        "kind": "stage diagnostic (AGENTS.md section 8): one oracle call on "
                "one rung, priced; nothing is inferred about a full ECDLP and "
                "no speedup is claimed",
        "command": f"./target/release/ic run --degree {args.degree} "
                   f"--summands {args.summands} --solver S --known-log K --json",
        "degree": args.degree,
        "summands": args.summands,
        "targets": list(TARGETS),
        "caveat": "the ratio is work per unit of solver output -- an F4 "
                  "reduction, an independent relation -- not strictly per "
                  "call, because the counters do not separate calls that "
                  "failed to decompose from calls that succeeded.  The "
                  "whole-run totals swing because a run decides how many "
                  "relations to collect, which is not a property of the target.",
        "solvers": solvers,
    }
    OUT.write_text(json.dumps(report, indent=2) + "\n")
    print(f"wrote {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""The programme's conformance runner, from B3 on: every step's cases, in
order, with the design's `until` rule applied
(`research/ic_tool_program/design/schema-v2.md` §9).

    python3 run.py --ic target/release/ic --through B3 [--build-commit SHA] [--out report.json]

The case sets, each frozen by its own step's `SHA256SUMS`:
- `v1/cases.json`, conformance suite v1: step B0's cases;
- `v2/cases.json`: B1's cases;
- every `v2-*/cases.json` (`v2-b3/`, …): the later steps' cases, each in
  its own directory so that no earlier step's pinned files change.

In a case, `{here}` is its own set's `params/`; `{cases}` is still
`v2/params/`, B1's frozen documents, which later cases copy and edit.

A case whose `until` step is at or before `--through` is not run; its
successor carries the new expectation.  Everything else — materialising
files, running, the expectations — is `v2/run.py`'s, unchanged.  This
runner's rules change only by a step's declaration.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("conformance_v2_run", HERE / "v2" / "run.py")
v2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(v2)

STEPS = ("B0", "B1", "B2", "B3", "B3b", "B4", "B5", "B6", "B7")


def case_sets() -> list[tuple[Path, list[dict]]]:
    sets = [(HERE / "v1", [{**c, "step": "B0"} for c in json.loads((HERE / "v1" / "cases.json").read_text())["cases"]]),
            (HERE / "v2", json.loads((HERE / "v2" / "cases.json").read_text())["cases"])]
    for d in sorted(p for p in HERE.glob("v2-*") if (p / "cases.json").is_file()):
        sets.append((d, json.loads((d / "cases.json").read_text())["cases"]))
    return sets


def cases_through(step: str) -> list[dict]:
    last = STEPS.index(step)

    def live(c: dict) -> bool:
        return STEPS.index(c["step"]) <= last and ("until" not in c or STEPS.index(c["until"]) > last)

    out = []
    for d, cases in case_sets():
        for c in cases:
            if live(c):
                out.append(json.loads(json.dumps(c).replace("{here}", str(d / "params"))))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ic", required=True, type=Path)
    ap.add_argument("--through", required=True, choices=STEPS)
    ap.add_argument("--build-commit")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    results = [v2.run_case(c, args.ic.resolve(), args.build_commit) for c in cases_through(args.through)]
    report = {"binary": str(args.ic), "through": args.through, "cases": len(results),
              "passed": sum(r["pass"] for r in results), "results": results}
    text = json.dumps(report, indent=1)
    if args.out:
        args.out.write_text(text + "\n")
    print(text)
    sys.exit(0 if report["passed"] == report["cases"] else 1)


if __name__ == "__main__":
    main()

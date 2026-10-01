#!/usr/bin/env python3
"""The conformance runner from B3 on: v1's and v2's cases with the design's
`until` rule applied, then B3's (`cases.json` here).

    python3 run.py --ic target/release/ic --through B3 [--build-commit SHA] [--out report.json]

A case whose `until` step is at or before `--through` is not run; its
successor carries the new expectation.  Everything else — materialising
files, running, the expectations — is `../v2/run.py`'s, unchanged.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
V2 = HERE.parent / "v2"
_spec = importlib.util.spec_from_file_location("conformance_v2_run", V2 / "run.py")
v2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(v2)

STEPS = ("B0", "B1", "B2", "B3", "B3b", "B4", "B5", "B6", "B7")


def with_b3_paths(case: dict) -> dict:
    return json.loads(json.dumps(case).replace("{b3}", str(HERE / "params")))


def cases_through(step: str) -> list[dict]:
    last = STEPS.index(step)
    v1 = [{**c, "step": "B0"} for c in json.loads((v2.V1 / "cases.json").read_text())["cases"]]
    b1 = json.loads((V2 / "cases.json").read_text())["cases"]
    b3 = json.loads((HERE / "cases.json").read_text())["cases"]

    def live(c: dict) -> bool:
        return STEPS.index(c["step"]) <= last and ("until" not in c or STEPS.index(c["until"]) > last)

    return [with_b3_paths(c) for c in v1 + b1 + b3 if live(c)]


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

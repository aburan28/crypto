#!/usr/bin/env python3
"""Apply the frozen Linux gate to first-qualified F5 segment receipts."""

import argparse
import hashlib
import json
from pathlib import Path


SEEDS = ("frozen", "holdout_a", "holdout_b", "holdout_c")
PRIMARY = "f5_n24_m24_d4"


def analyze(root: Path, threads: int) -> dict:
    index = json.loads((root / "INDEX.json").read_text())
    errors = []
    rows = {}
    if index.get("status") != "qualified" or index.get("threads") != threads:
        errors.append("segment qualification incomplete or thread count mismatched")
    for seed in SEEDS:
        chosen = index.get("seeds", {}).get(seed, {}).get("selected_attempt")
        if not isinstance(chosen, int):
            errors.append(f"{seed}: no qualified attempt")
            continue
        directory = root / seed / f"attempt{chosen}"
        status = json.loads((directory / "STATUS.json").read_text())
        raw = (directory / "paired.json").read_bytes()
        receipt = json.loads(raw)
        if not status.get("qualified") or receipt.get("status") != "complete":
            errors.append(f"{seed}: selected attempt not complete and qualified")
            continue
        summary = receipt["summary"][seed]
        for case, fields in summary.items():
            wall = fields["wall_ms"]
            median = wall["prior_new"]["median"]
            if case != PRIMARY and median < wall["aa_min"]:
                errors.append(f"{seed}/{case}: smaller-case regression")
            if case == PRIMARY:
                if threads == 1 and median < 2.0:
                    errors.append(f"{seed}: primary median below 2.00x")
                if threads == 1 and wall["prior_new"]["bootstrap_95pct"][0] <= 2.0:
                    errors.append(f"{seed}: primary interval lower bound at or below 2.00x")
                if threads == 2 and median < wall["aa_min"]:
                    errors.append(f"{seed}: two-thread primary regression")
                rows[seed] = {
                    "selected_attempt": chosen,
                    "paired_sha256": hashlib.sha256(raw).hexdigest(),
                    "complete_call_ratio": median,
                    "bootstrap_95pct": wall["prior_new"]["bootstrap_95pct"],
                    "aa_range": [wall["aa_min"], wall["aa_max"]],
                    "rank": receipt["runs"][0]["cases"][PRIMARY]["rank"],
                    "row_space_fp": receipt["runs"][0]["cases"][PRIMARY]["row_space_fp"],
                }
    return {
        "status": "pass" if not errors and len(rows) == len(SEEDS) else "fail",
        "threads": threads,
        "errors": errors,
        "primary": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--threads", type=int, choices=(1, 2), required=True)
    args = parser.parse_args()
    result = analyze(args.root, args.threads)
    (args.root / "CI_SUMMARY.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))
    return 0 if result["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())

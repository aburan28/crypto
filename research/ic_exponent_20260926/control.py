#!/usr/bin/env python3
"""Control 1 of ledger §20: the pricer's counts are the workflow's.

Compares an `ic price` report with an `ic workflow` report (and its
relation files) on the same parameter file, field by field:

- base: points, signed orbits, columns, selection cost, pair-table tier
  and stored pairs;
- collection: units, and per unit its relations, probes and summands
  scanned;
- logs: relations accepted, rejected and duplicate;
- descent: trials and recovered logarithm, target by target.

Prints a JSON verdict; exits non-zero on any difference.

    python3 control.py <price.json> <workflow.json> <workflow run dir>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path


def main() -> int:
    price = json.loads(Path(sys.argv[1]).read_text())
    wf = json.loads(Path(sys.argv[2]).read_text())
    run_dir = Path(sys.argv[3])
    pc = price["counts"]
    stages = {s["stage"]: s for s in wf["stages"]}
    fb = wf["factor_base"]
    units = []
    for path in sorted((run_dir / "relations").glob("unit-*.json")):
        doc = json.loads(path.read_text())
        units.append([doc["unit"], len(doc["relations"]), doc["count"], doc["summands_scanned"]])
    sol = sorted(wf["solutions"]["items"], key=lambda s: s["index"])
    checks = {
        "points": (pc["select"]["points"], fb["points"]),
        "signed_orbits": (pc["select"]["signed_orbits"], fb["signed_orbits"]),
        "columns": (pc["select"]["columns"], fb["columns"]),
        "selection_cost": (pc["select"]["selection_cost"], fb["selection_cost"]),
        "tier": (pc["build"]["tier"], fb.get("pair_table_tier")),
        "stored_pairs": (pc["build"]["stored_pairs"], fb.get("pair_table_stored_pairs")),
        "per_unit": (pc["collect"]["per_unit_index_relations_trials_summands"], units),
        "relations_accepted": (pc["logs"]["relations_accepted"], stages["logs"]["relations"]),
        "rejected": (pc["logs"]["rejected"], stages["logs"]["rejected"]),
        "duplicates": (pc["logs"]["duplicates"], stages["logs"]["duplicates"]),
        "units_extended": (pc["collect"]["units_extended"], stages["logs"]["units_extended"]),
        "descent_trials": (pc["descent"]["trials_per_target"], [s["descent_trials"] for s in sol]),
        "recovered": (price["recovered"], [s["recovered"] for s in sol]),
        "workflow_verified": (True, all(s["verified"] for s in sol) and wf["status"] == "complete"),
    }
    diffs = {k: {"price": a, "workflow": b} for k, (a, b) in checks.items() if a != b}
    print(json.dumps({"control": "price counts == workflow counts", "pass": not diffs,
                      "checked": sorted(checks), "differences": diffs}, indent=1))
    return 0 if not diffs else 1


if __name__ == "__main__":
    sys.exit(main())

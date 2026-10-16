#!/usr/bin/env python3
"""Independently replay the revised supervisor's source-only receipt."""

import hashlib
import json
from pathlib import Path
import subprocess


DIRECTORY = Path(__file__).resolve().parent
ROOT = DIRECTORY.parents[3]


def verify() -> dict:
    receipt = json.loads((DIRECTORY / "validation.json").read_text())
    assert receipt["schema"] == "n83.mixed-orbit-search-source-validation/v1"
    assert receipt["status"] == "SOURCE_ONLY_PASS"
    assert receipt["prior_total_charged_seconds"] == "12600"
    assert receipt["new_search_authorized"] is False
    assert receipt["retained_n83_search_executed"] is False
    assert receipt["total_index_calculus_runtime_ms"] is None
    assert receipt["selected_best_total_runtime"] is None
    assert receipt["study_python_tests_passed"] == 46
    commit = receipt["source_commit"]
    for relative, expected in receipt["sha256"].items():
        data = subprocess.run(["git", "show", f"{commit}:{relative}"],
                              cwd=ROOT, check=True, capture_output=True).stdout
        assert hashlib.sha256(data).hexdigest() == expected, relative
    log = (DIRECTORY / "study-python.log").read_text()
    assert "Ran 46 tests" in log and log.rstrip().endswith("OK")
    prior = json.loads((DIRECTORY.parent / "mixed-orbit-import-900s-20261012/audit.json").read_text())
    assert prior["schema"] == "n83.mixed-orbit-import-900-audit/v1"
    assert prior["status"] == "PASS" and prior["conservative_charge_seconds"] == 900
    assert prior["solver_search_executed"] is False
    return {"schema": "n83.mixed-orbit-search-source-replay/v1",
            "status": "PASS", "source_commit": commit,
            "checked_files": len(receipt["sha256"]),
            "study_python_tests_passed": 46,
            "new_search_authorized": False}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))

#!/usr/bin/env python3
"""Replay the source-bound residual-domain audit receipt."""

import hashlib
import json
from pathlib import Path
import subprocess


DIRECTORY = Path(__file__).resolve().parent
ROOT = DIRECTORY.parents[3]


def verify() -> dict:
    receipt = json.loads((DIRECTORY / "validation.json").read_text())
    assert receipt["schema"] == "n83.mixed-orbit-residual-domain-source-validation/v1"
    assert receipt["status"] == "SOURCE_ONLY_PASS"
    assert receipt["study_python_tests_passed"] == 46
    assert receipt["group_valid_outside_suffix_rejected"] is True
    assert receipt["retained_n83_search_executed"] is False
    assert receipt["additional_experimental_charge_seconds"] == 0
    assert receipt["total_index_calculus_runtime_ms"] is None
    assert receipt["selected_best_total_runtime"] is None
    for relative, expected in receipt["sha256"].items():
        source = subprocess.run(
            ["git", "show", f"{receipt['source_commit']}:{relative}"],
            cwd=ROOT, check=True, capture_output=True,
        ).stdout
        assert hashlib.sha256(source).hexdigest() == expected, relative
    log = (DIRECTORY / "study-python.log").read_text()
    assert "Ran 46 tests" in log and log.rstrip().endswith("OK")
    return {"schema": "n83.mixed-orbit-residual-domain-source-replay/v1",
            "status": "PASS", "source_commit": receipt["source_commit"],
            "checked_files": len(receipt["sha256"]),
            "study_python_tests_passed": 46,
            "retained_n83_search_executed": False}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))

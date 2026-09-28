#!/usr/bin/env python3
"""Compare derived decisions exactly, allowing only 4-ULP Python float drift."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
LEDGER_SHA = "84ad0e494880e51360491276e130aea93ec16ab095ca22fe6ee53f491612f207"
NEGATIVE = "NEGATIVE_FULLY_CHARGED_SAME_HOST_WALL_OBSERVATION"


def compare(expected, actual, path: str, drift: list[tuple[str, float]]) -> None:
    assert type(expected) is type(actual), f"type differs at {path}"
    if isinstance(expected, dict):
        assert set(expected) == set(actual), f"keys differ at {path}"
        for key in sorted(expected):
            compare(expected[key], actual[key], f"{path}.{key}", drift)
    elif isinstance(expected, list):
        assert len(expected) == len(actual), f"length differs at {path}"
        for index, (left, right) in enumerate(zip(expected, actual)):
            compare(left, right, f"{path}[{index}]", drift)
    elif isinstance(expected, float):
        assert math.isfinite(expected) and math.isfinite(actual), f"nonfinite at {path}"
        delta = abs(expected-actual)
        tolerance = 4*max(math.ulp(expected), math.ulp(actual))
        assert delta <= tolerance, f"float differs at {path}: {expected} vs {actual}"
        if delta:
            drift.append((path, delta))
    else:
        assert expected == actual, f"value differs at {path}: {expected!r} vs {actual!r}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("expected", type=Path)
    parser.add_argument("actual", type=Path)
    args = parser.parse_args()
    expected = json.loads(args.expected.read_bytes())
    actual = json.loads(args.actual.read_bytes())
    drift = []
    compare(expected, actual, "$", drift)
    for analysis in (expected, actual):
        assert analysis["n131_transfer"] is None
        assert analysis["common_group_addition_equivalent_S"] is None
        for n in ("37", "41"):
            for length in ("8", "32"):
                level = analysis["arms"][n]["levels"][length]
                assert level["classification"] == NEGATIVE
                assert level["all_three_pairs_complete"]
                assert len(level["blocks"]) == 3
    ledger = (HERE / "evidence/SHA256SUMS").read_bytes()
    assert hashlib.sha256(ledger).hexdigest() == LEDGER_SHA
    print(json.dumps({"classification": "DERIVED_DECISION_COMPARE_PASS",
                      "float_drift_fields": len(drift),
                      "maximum_abs_float_drift": max((delta for _,delta in drift), default=0.0),
                      "ledger_sha256": LEDGER_SHA}, sort_keys=True))


if __name__ == "__main__":
    main()

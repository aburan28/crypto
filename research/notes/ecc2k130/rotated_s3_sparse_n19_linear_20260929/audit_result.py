#!/usr/bin/env python3
"""Read back a complete new semantic replay without treating timing as attack cost."""
from __future__ import annotations

import argparse
import json

import verify as v


CONTROL_NAMES = {"sign_incomplete", "zero_zero_to_O", "O_plus_zero",
                 "zero_double_to_O", "wrong_O_target"}


def audit(path):
    frozen = json.loads((v.HERE / "FROZEN.json").read_text())
    v.check_freeze(frozen)
    row = json.loads(path.read_text())
    if (row.get("decision") != "PASS" or row.get("primary_paths") != 58825 or
        row.get("signed_point_tuples") != 117649 or
        row.get("producer_sha256") != frozen["producer_result_sha256"] or
        set(row.get("negative_controls", {})) != CONTROL_NAMES):
        raise AssertionError("semantic replay outcome/count/control mismatch")
    targets = row.get("target_rows")
    if (not isinstance(targets, list) or len(targets) != 33 or
        [t.get("id") for t in targets] !=
        [f"Q{i}T{j}" for i in range(8) for j in range(4)] + ["O"] or
        any(not isinstance(t.get("exact_signed_tuple_count"), int) or
            not isinstance(t.get("cnf_primary_path_count"), int) or
            not isinstance(t.get("assumption_literal"), int) for t in targets)):
        raise AssertionError("target replay schema")
    if (len(row.get("transition_cases", [])) != 5 or
        not all(x.get("stage") == f"S{i}" for i, x in enumerate(
            row["transition_cases"], start=2))):
        raise AssertionError("stage replay schema")
    for key in ("wall_seconds", "cpu_seconds", "peak_rss_bytes"):
        if not isinstance(row.get(key), (int, float)) or row[key] < 0:
            raise AssertionError(f"invalid {key}")
    return {"decision": "PASS", "primary_paths": row["primary_paths"],
            "signed_point_tuples": row["signed_point_tuples"],
            "targets": len(targets), "negative_controls": len(CONTROL_NAMES),
            "wall_seconds": row["wall_seconds"],
            "peak_rss_bytes": row["peak_rss_bytes"],
            "attack_cost_admitted": False}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", required=True, type=v.Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.result), sort_keys=True))


if __name__ == "__main__":
    main()

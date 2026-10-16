#!/usr/bin/env python3
"""Audit stage transport and classify later bounded N83 probe traces fail-closed."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCE = ROOT / "examples/koblitz_n83_fixed_witness_stage_probe.rs"
SOURCE_SHA256 = "f13751653b2470ade38279945b999ba825e10c39cae905682992d4f49a01ff33"
WORKER_SHA256 = "5c87e73d6e786588178e4c04ea447b2d115b2933c2cd77b675d7ee4c97c8c92c"
MARKER_STAGES = (
    "process_entry", "input_validated", "core_built", "domains_installed",
    "pair_ordered", "units_pinned", "solver_entered", "solver_returned",
    "model_checked",
)
HEX64 = re.compile(r"[0-9a-f]{64}\Z")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(raw):
    return hashlib.sha256(raw).hexdigest()


def classify(stdout, stderr, exit_code, *, self_test=False):
    """Return stage evidence without treating a prefix as SAT, UNSAT, or yield."""
    require(exit_code in (0, 124, 137), "unclassified process exit")
    lines = stderr.splitlines()
    require(lines and all(line for line in lines), "missing or blank marker line")
    try:
        rows = [json.loads(line) for line in lines]
    except (TypeError, ValueError) as error:
        raise ValueError("malformed stage JSON") from error
    expected = ("process_entry", "self_test_complete") if self_test else MARKER_STAGES
    require(len(rows) <= len(expected), "too many markers")
    previous_elapsed = -1.0
    input_hash = None
    previous_size = None
    for index, row in enumerate(rows):
        require(isinstance(row, dict)
                and row.get("schema") == "n83.fixed-witness-stage-marker/v1"
                and row.get("sequence") == index
                and row.get("stage") == expected[index],
                "stage order or schema differs")
        elapsed = row.get("elapsed_process_wall_ms_L0")
        require(isinstance(elapsed, (float, int)) and math.isfinite(elapsed)
                and elapsed >= previous_elapsed and elapsed >= 0,
                "stage wall time differs")
        previous_elapsed = elapsed
        if self_test or index == 0:
            require(row.get("input_blake3") is None
                    and row.get("solver_size") is None,
                    "unbound stage has input or solver")
        else:
            value = row.get("input_blake3")
            require(isinstance(value, str) and HEX64.fullmatch(value),
                    "stage input digest missing")
            if input_hash is None:
                input_hash = value
            require(value == input_hash, "input digest changed within run")
            size = row.get("solver_size")
            if index == 1:
                require(size is None, "input stage has solver")
            else:
                require(isinstance(size, dict)
                        and set(size) == {"variables", "clauses", "xors"}
                        and all(isinstance(size[key], int) and size[key] >= 0
                                for key in size), "solver size invalid")
                if previous_size is not None:
                    require(all(size[key] >= previous_size[key] for key in size),
                            "solver size decreased")
                previous_size = size
    if self_test:
        require(exit_code == 0 and len(rows) == 2 and len(stdout.splitlines()) == 1,
                "marker self-test did not complete")
        result = json.loads(stdout)
        require(result == {"schema": "n83.fixed-witness-marker-self-test/v1",
                           "status": "PASS", "stage_markers_emitted": 2},
                "marker self-test result differs")
        status = "PASS_marker_transport_only"
    elif exit_code == 124:
        require(not stdout, "wall-capped probe emitted terminal stdout")
        status = ("CENSORED_after_solver_entry_marker" if len(rows) >= 7
                  else "CENSORED_before_solver_entry_marker")
    elif exit_code == 137:
        require(not stdout, "resource-capped probe emitted terminal stdout")
        status = "CENSORED_resource_exit_unverified"
    else:
        require(len(rows) == len(MARKER_STAGES) and len(stdout.splitlines()) == 1,
                "successful process lacks terminal stage or receipt")
        result = json.loads(stdout)
        require(result.get("schema") == "n83.fixed-witness-sat-control/v2"
                and result.get("stage_markers_emitted") == 9
                and result.get("input_blake3") == input_hash,
                "terminal result differs")
        worker_status = result.get("status")
        if worker_status == "SAT_verified_witness_x":
            require(result.get("model_verified") is True, "SAT model not checked")
            status = "PASS_fixed_witness_encoding_only"
        elif worker_status == "UNKNOWN_conflict_cap":
            status = "CENSORED_conflict_cap"
        elif worker_status == "UNSAT_fixed_witness_encoding":
            status = "ENCODING_CONTRADICTION_requires_diagnosis"
        else:
            raise ValueError("unclassified terminal worker status")
    return {"status": status, "last_completed_stage": rows[-1]["stage"],
            "marker_count": len(rows), "input_blake3": input_hash,
            "solver_entry_marker_seen": len(rows) >= 7 and not self_test,
            "natural_relations": None, "selected_best_total_runtime": None}


def selftest_receipt():
    source = SOURCE.read_bytes()
    worker = (HERE / "worker").read_bytes()
    require(sha256(source) == SOURCE_SHA256 and sha256(worker) == WORKER_SHA256,
            "source or compiled worker changed")
    stdout = (HERE / "selftest.stdout").read_bytes()
    stderr = (HERE / "selftest.stderr").read_bytes()
    row = classify(stdout, stderr, 0, self_test=True)
    return {"schema": "n83.fixed-witness-stage-preflight/v1", **row,
            "source_sha256": sha256(source), "worker_sha256": sha256(worker),
            "stdout_sha256": sha256(stdout), "stderr_sha256": sha256(stderr),
            "solver_search_executed": False,
            "total_index_calculus_runtime_ms": None}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    row = selftest_receipt()
    rendered = json.dumps(row, sort_keys=True, indent=2) + "\n"
    path = HERE / "selftest-audit.json"
    if args.write:
        with path.open("x") as stream:
            stream.write(rendered)
    else:
        require(path.read_text() == rendered, "self-test audit differs")
    print(json.dumps({"status": row["status"], "marker_count": row["marker_count"]}))


if __name__ == "__main__":
    main()

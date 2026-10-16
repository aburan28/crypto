#!/usr/bin/env python3
"""Recheck focused validation without promoting an interrupted full suite."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BASELINE = "bbbf03e94782aca600de5345ce65e7739c9cf4fd"
SOURCE_DELTA_SHA256 = "b9b3b636cfb5e93669073977cdf1f92bbb3831f2933399ef5eb42188dc0e6175"
FOCUSED = "cryptanalysis::koblitz_large_prime::tests::mixed_orbit_domains_lift_disjoint_partial_and_fail_closed_on_caps"
LOGS = ("study-python", "boundary-python", "example-release", "focused-library",
        "release-lib", "release-lib-interrupted", "release-lib-background-launch-failed")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def check():
    logs = {}
    for label in LOGS:
        raw = (HERE / f"{label}.log.gz").read_bytes()
        logs[label] = {"compressed_sha256": digest(raw), "plain": gzip.decompress(raw).decode()}
    require("Ran 63 tests" in logs["study-python"]["plain"]
            and logs["study-python"]["plain"].rstrip().endswith("OK"),
            "N83 Python suite did not pass")
    require("Ran 43 tests" in logs["boundary-python"]["plain"]
            and logs["boundary-python"]["plain"].rstrip().endswith("OK"),
            "boundary Python suite did not pass")
    require("test result: ok. 0 passed; 0 failed" in logs["example-release"]["plain"],
            "current example release test did not exit cleanly")
    require(f"test {FOCUSED} ... ok" in logs["focused-library"]["plain"]
            and "test result: ok. 1 passed; 0 failed" in logs["focused-library"]["plain"],
            "focused current library test did not pass")
    partial = logs["release-lib"]["plain"]
    completed = len(re.findall(r"^test .* \.\.\. (?:ok|ignored|FAILED)", partial, re.M))
    require(completed == 1836 and "signal: 15, SIGTERM" in partial
            and "test result: ok." not in partial, "full rerun disposition differs")
    exit_record = json.loads((HERE / "release-lib-exit.json").read_text())
    require(exit_record["exit_code"] == 101
            and exit_record["log_bytes"] == len(partial.encode()),
            "full rerun exit receipt differs")
    attempts = json.loads((HERE / "validation-attempts.json").read_text())
    require([row["disposition"] for row in attempts["attempts"]] == [
        "INTERRUPTED_NO_TERMINAL_RECEIPT", "LAUNCH_FAILURE",
        "INTERRUPTED_BY_OPERATOR_UNDER_HOST_CONTENTION"],
        "infrastructure attempt labels differ")
    for row in attempts["attempts"]:
        label = row["partial_log"].removesuffix(".log.gz")
        require(row["compressed_sha256"] == logs[label]["compressed_sha256"],
                "infrastructure log digest differs")
    delta = subprocess.run(
        ["git", "diff", BASELINE, "HEAD", "--", "src", "Cargo.toml", "Cargo.lock"],
        cwd=ROOT, check=True, capture_output=True).stdout
    require(digest(delta) == SOURCE_DELTA_SHA256,
            "production-library delta from prior complete pass changed")
    prior_result = ROOT / "research/koblitz_n83_factor_base_sweep_20261008/verification/cross-stream-partial-worker-source-20261013/RESULT.md"
    require("4,001 tests with 113 ignored" in prior_result.read_text(),
            "documented prior full pass missing")
    return {
        "schema": "n83.fixed-witness-validation/v1",
        "status": "PASS_focused_checks_full_rerun_interrupted",
        "documented_prior_complete_library_commit": BASELINE,
        "documented_prior_complete_library": {"passed": 4001, "ignored": 113,
                                               "raw_log_retained": False},
        "source_delta_sha256": SOURCE_DELTA_SHA256,
        "current_focused_library_passed": 1,
        "current_example_compiled_tests": 0,
        "current_study_python_passed": 63,
        "current_boundary_python_passed": 43,
        "current_full_library_completed_lines": completed,
        "current_full_library_status": "INTERRUPTED_SIGTERM_under_host_contention",
        "log_compressed_sha256": {label: logs[label]["compressed_sha256"] for label in LOGS},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    row = check()
    rendered = json.dumps(row, sort_keys=True, indent=2) + "\n"
    path = HERE / "validation.json"
    if args.write:
        with path.open("x") as output:
            output.write(rendered)
    else:
        require(path.read_text() == rendered, "validation receipt changed")
    print(json.dumps({"status": row["status"],
                      "full_library": row["current_full_library_status"]}))


if __name__ == "__main__":
    main()

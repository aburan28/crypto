#!/usr/bin/env python3
"""Independently classify the bounded N83 fixed-witness SAT control receipt."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
STUDY = HERE.parents[1]
SOURCE_COMMIT = "55d5ed6771ad53ff6dcfe6e537ac3ec4598984bb"
DOMAINS_SHA256 = "63accc5b9dc9de6327c5b07e93b3e62717f8376a96c3f91f7bb70b77095f23bb"
SOURCE_SHA256 = "06730771fda845339fdd14908ded0e259e0294ebc37d62af39e5444bc6079088"
BINARY_SHA256 = "64ce185bfb1a4db69dcc576fb0749108b56e804efd57f14da561f07d41f51752"
REPLAY_SHA256 = "1155f9fbad9ff6b0b6d17200fd15880e0a7612a9db2647fa50a3ef5d148f2579"
S3_AUDIT_SHA256 = "589d1c2d8577d596624e75df1ea15dfa09b52d61d4b1081eb2e01093331da843"
SOURCE_RELATIVE = "examples/koblitz_n83_fixed_witness_control.rs"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check():
    source = ROOT / SOURCE_RELATIVE
    require(digest(source) == SOURCE_SHA256, "control source changed")
    committed = subprocess.run(["git", "show", f"{SOURCE_COMMIT}:{SOURCE_RELATIVE}"],
                               cwd=ROOT, check=True, capture_output=True).stdout
    require(hashlib.sha256(committed).hexdigest() == SOURCE_SHA256,
            "control source Git blob differs")
    require(digest(HERE / "domains.json") == DOMAINS_SHA256, "prepared domains changed")
    require(digest(HERE / "worker") == BINARY_SHA256, "archived native worker changed")
    domains = json.loads((HERE / "domains.json").read_text())
    require(len(domains["small_reps"]) == 1182
            and len(domains["residual_reps"]) == 866
            and domains["ordered_indices"] == [166, 0, 332, 196212, 196378]
            and domains["identity_mask"] == 0, "fixed domain shape differs")
    fixture = STUDY / "verification/mixed-orbit-known-witness-source-20261012/fixture.json"
    require(digest(fixture) == domains["fixture_sha256"], "fixture binding differs")
    for label in ("small", "envelope"):
        require(domains["archive_sha256"][label]
                == json.loads(fixture.read_text())["archive"][label]["sha256"],
                f"{label} archive binding differs")
    replay_path = STUDY / "verification/known-witness-first-pattern-replay-20261014/replay.json"
    require(digest(replay_path) == REPLAY_SHA256, "independent algebra replay changed")
    replay = json.loads(replay_path.read_text())
    require(replay["status"] == "PASS_ordered_group_and_S3_witness_for_first_identity_pattern"
            and replay["sat_order_indices"] == domains["ordered_indices"]
            and [row[0] for row in replay["sat_order_partial_sums"][1:-1]]
                == domains["partial_xs"], "algebra replay binding differs")
    s3_path = STUDY / "verification/s3-stream-integrity-20261014/audit.json"
    require(digest(s3_path) == S3_AUDIT_SHA256, "S3 content audit changed")
    s3 = json.loads(s3_path.read_text())
    require(s3["status"] == "PASS" and s3["object_count"] == 57,
            "S3 current-content audit invalid")
    launch = json.loads((HERE / "launch.json").read_text())
    require(launch["schema"] == "n83.fixed-witness-control-launch/v1"
            and launch["prior_charged_seconds"] == 21000
            and launch["conservative_new_charge_seconds"] == 1200
            and launch["cumulative_charged_seconds"] == 22200
            and launch["wall_cap_seconds"] == 840
            and launch["conflict_cap"] == 5000
            and launch["source_commit"] == SOURCE_COMMIT
            and launch["binary_sha256"] == BINARY_SHA256
            and launch["domains_sha256"] == DOMAINS_SHA256,
            "launch contract differs")
    raw_out = (HERE / "solver.stdout").read_bytes()
    raw_err = (HERE / "solver.stderr").read_bytes()
    exit_record = json.loads((HERE / "exit.json").read_text())
    require(exit_record["schema"] == "n83.fixed-witness-control-exit/v1"
            and exit_record["exit_code"] in (0, 124, 137)
            and exit_record["wall_cap_seconds"] == 840
            and exit_record["stdout_bytes"] == len(raw_out)
            and exit_record["stderr_bytes"] == len(raw_err), "exit record differs")
    result = None
    if exit_record["exit_code"] == 0:
        require(len(raw_out.splitlines()) == 1, "solver stdout not one JSON record")
        result = json.loads(raw_out)
        require(result["schema"] == "n83.fixed-witness-sat-control/v1"
                and result["scope"] == "known_answer_encoding_control_only"
                and result["fixture_sha256"] == domains["fixture_sha256"]
                and result["identity_mask"] == 0
                and result["fixed_x_blocks"] == 8
                and result["before_fixed_units"]["domain_clauses"] == 465738
                and result["conflict_cap"] == 5000
                and 0 <= result["conflicts"] <= 5000
                and 0 <= result["elapsed_wall_ms_L0"] <= 840000,
                "solver receipt counters differ")
        if result["status"] == "SAT_verified_witness_x":
            require(result["model_verified"] is True, "SAT model not checked")
            disposition = "PASS_fixed_witness_encoding_sat"
        elif result["status"] == "UNKNOWN_conflict_cap":
            require(result["model_verified"] is False, "unknown has model")
            disposition = "VALIDATED_CENSORED"
        elif result["status"] == "UNSAT_fixed_witness_encoding":
            require(result["model_verified"] is False, "UNSAT has model")
            disposition = "ENCODING_CONTRADICTION_REQUIRES_DIAGNOSIS"
        else:
            disposition = "PRODUCER_FAILURE"
    else:
        require(not raw_out, "timed out solver produced unexpected stdout")
        disposition = "VALIDATED_CENSORED"
    return {
        "schema": "n83.fixed-witness-control-audit/v1",
        "status": disposition,
        "scope": "known_answer_encoding_control_only",
        "source_commit": SOURCE_COMMIT,
        "prior_charged_seconds": 21000,
        "conservative_new_charge_seconds": 1200,
        "cumulative_charged_seconds": 22200,
        "exit_code": exit_record["exit_code"],
        "solver_status": result["status"] if result else None,
        "conflicts": result["conflicts"] if result else None,
        "model_verified": result["model_verified"] if result else None,
        "raw_sha256": {"solver.stdout": hashlib.sha256(raw_out).hexdigest(),
                       "solver.stderr": hashlib.sha256(raw_err).hexdigest(),
                       "launch.json": digest(HERE / "launch.json"),
                       "exit.json": digest(HERE / "exit.json")},
        "natural_relations": None,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    row = check()
    rendered = json.dumps(row, sort_keys=True, indent=2) + "\n"
    path = HERE / "audit.json"
    if args.write:
        with path.open("x") as stream:
            stream.write(rendered)
    else:
        require(path.read_text() == rendered, "audit receipt changed")
    print(json.dumps({"status": row["status"], "solver_status": row["solver_status"]}))


if __name__ == "__main__":
    main()

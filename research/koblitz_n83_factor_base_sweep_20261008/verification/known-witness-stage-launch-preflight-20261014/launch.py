#!/usr/bin/env python3
"""Pin and watch one native N83 stage probe; preserve and recheck raw output."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import sys
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PREP = HERE.parent / "known-witness-stage-preflight-20261014"
CHECKER = PREP / "check_trace.py"
WORKER = PREP / "worker"
SOURCE = ROOT / "examples/koblitz_n83_fixed_witness_stage_probe.rs"
DOMAINS = HERE.parent / "known-witness-ordered-sat-20261014/domains.json"
PINS = {
    "checker": "44d5eefb615ce4002840624094a1b7945a889c0ef8a9d1d5c5c5cf26c8c8fd0f",
    "worker": "5c87e73d6e786588178e4c04ea447b2d115b2933c2cd77b675d7ee4c97c8c92c",
    "source": "f13751653b2470ade38279945b999ba825e10c39cae905682992d4f49a01ff33",
    "domains": "63accc5b9dc9de6327c5b07e93b3e62717f8376a96c3f91f7bb70b77095f23bb",
}
PRIOR_CHARGE = 22200
CONFLICT_CAP = 5000
FULL_WALL_CAP = 840
SELFTEST_WALL_CAP = 5
MIN_FREE_BYTES = 256 * 1024 * 1024
sys.path.insert(0, str(PREP))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def check_pins():
    for name, path in (("checker", CHECKER), ("worker", WORKER),
                       ("source", SOURCE), ("domains", DOMAINS)):
        require(path.is_file() and sha256(path) == PINS[name],
                f"pinned {name} changed")
    require(os.access(WORKER, os.X_OK), "pinned worker is not executable")


def authorization_error(row):
    if (not isinstance(row, dict) or
            row.get("schema") != "n83.fixed-witness-stage-authorization/v1" or
            row.get("scope") != "known_answer_encoding_control_only" or
            not isinstance(row.get("approval_reference"), str) or
            not row["approval_reference"].strip() or
            row.get("prior_charged_seconds") != PRIOR_CHARGE or
            type(row.get("new_authorized_seconds")) is not int or
            row["new_authorized_seconds"] < FULL_WALL_CAP + 30 or
            row.get("wall_cap_seconds") != FULL_WALL_CAP or
            row.get("conflict_cap") != CONFLICT_CAP or
            row.get("domains_sha256") != PINS["domains"]):
        return "full-domain authorization receipt differs"
    return None


def save_json(path, row):
    with path.open("x") as output:
        json.dump(row, output, sort_keys=True, indent=2)
        output.write("\n")
        output.flush()
        os.fsync(output.fileno())


def group_gone(pgid):
    try:
        os.killpg(pgid, 0)
    except ProcessLookupError:
        return True
    except PermissionError:
        return False
    return False


def launch_child(command, wall, output):
    started = time.monotonic()
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    timed_out = False
    launch_error = None
    cleanup_error = None
    child_returncode = None
    pgid = None
    with (output / "stdout.raw").open("xb") as stdout, (
        output / "stderr.raw"
    ).open("xb") as stderr:
        try:
            child = subprocess.Popen(command, stdout=stdout, stderr=stderr,
                                     start_new_session=True)
            pgid = child.pid
            try:
                child_returncode = child.wait(timeout=wall)
            except subprocess.TimeoutExpired:
                timed_out = True
                try:
                    os.killpg(pgid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
                try:
                    child_returncode = child.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(pgid, signal.SIGKILL)
                    child_returncode = child.wait(timeout=5)
                if not group_gone(pgid):
                    try:
                        os.killpg(pgid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
        except (OSError, subprocess.SubprocessError) as error:
            launch_error = f"{type(error).__name__}: {error}"
        if pgid is not None and not group_gone(pgid):
            try:
                os.killpg(pgid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            for _ in range(20):
                if group_gone(pgid):
                    break
                time.sleep(0.05)
            if not group_gone(pgid):
                cleanup_error = "child process group still present"
        stdout.flush()
        stderr.flush()
        os.fsync(stdout.fileno())
        os.fsync(stderr.fileno())
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    return {
        "schema": "n83.fixed-witness-stage-child-exit/v1",
        "timed_out": timed_out,
        "child_returncode": child_returncode,
        "watchdog_exit_code": 124 if timed_out else
            (137 if child_returncode == -signal.SIGKILL else child_returncode),
        "launch_error": launch_error,
        "cleanup_error": cleanup_error,
        "process_group_cleanup_verified": cleanup_error is None,
        "wall_elapsed_ms_L0": (time.monotonic() - started) * 1000,
        "child_user_cpu_seconds_L0": after.ru_utime - before.ru_utime,
        "child_system_cpu_seconds_L0": after.ru_stime - before.ru_stime,
        "stdout_sha256": sha256(output / "stdout.raw"),
        "stderr_sha256": sha256(output / "stderr.raw"),
    }


def audit(output):
    check_pins()
    launch = json.loads((output / "launch.json").read_text())
    exit_row = json.loads((output / "exit.json").read_text())
    require(launch.get("schema") == "n83.fixed-witness-stage-launch/v1" and
            exit_row.get("schema") == "n83.fixed-witness-stage-child-exit/v1" and
            launch.get("pins") == PINS and
            launch.get("launcher_sha256") == sha256(Path(__file__)) and
            launch.get("prior_charged_seconds") == PRIOR_CHARGE and
            launch.get("scope") == "known_answer_encoding_control_only",
            "launch or exit schema/pins differ")
    mode = launch.get("mode")
    require(mode in ("self_test", "full_domain"), "launch mode differs")
    require(launch.get("wall_cap_seconds") ==
            (SELFTEST_WALL_CAP if mode == "self_test" else FULL_WALL_CAP),
            "wall cap differs")
    require(launch.get("logical_command") ==
            (["pinned_worker", "self-test-markers"] if mode == "self_test"
             else ["pinned_worker", "pinned_domains", str(CONFLICT_CAP)]),
            "logical worker command differs")
    if mode == "full_domain":
        auth_path = output / "authorization.json"
        require(sha256(auth_path) == launch.get("authorization_sha256"),
                "authorization bytes differ")
        auth = json.loads(auth_path.read_text())
        require(authorization_error(auth) is None and
                launch.get("conservative_new_charge_seconds") ==
                auth["new_authorized_seconds"] and
                launch.get("conflict_cap") == CONFLICT_CAP,
                "authorization or charge differs")
    else:
        require(launch.get("authorization_sha256") is None and
                launch.get("conservative_new_charge_seconds") == 0 and
                launch.get("conflict_cap") is None,
                "self-test launch charged or used full-domain inputs")
    stdout = (output / "stdout.raw").read_bytes()
    stderr = (output / "stderr.raw").read_bytes()
    require(sha256(output / "stdout.raw") == exit_row.get("stdout_sha256") and
            sha256(output / "stderr.raw") == exit_row.get("stderr_sha256"),
            "raw output differs from exit receipt")
    require(exit_row.get("watchdog_exit_code") ==
            (124 if exit_row.get("timed_out") else
             (137 if exit_row.get("child_returncode") == -signal.SIGKILL else
              exit_row.get("child_returncode"))),
            "watchdog exit disagrees with child exit")
    trace_error = None
    trace = None
    try:
        import check_trace  # Imported only after its source hash is checked.
        trace = check_trace.classify(stdout, stderr,
                                     exit_row["watchdog_exit_code"],
                                     self_test=mode == "self_test")
    except (ValueError, KeyError, TypeError) as error:
        trace_error = f"{type(error).__name__}: {error}"
    if exit_row.get("cleanup_error") or not exit_row.get("process_group_cleanup_verified"):
        status = "UNVERIFIED_process_cleanup"
    elif exit_row.get("launch_error"):
        status = "PRODUCER_FAILURE_launch"
    elif trace_error:
        status = "UNVERIFIED_trace"
    else:
        status = trace["status"]
    return {
        "schema": "n83.fixed-witness-stage-launch-audit/v1",
        "status": status,
        "mode": mode,
        "launch_sha256": sha256(output / "launch.json"),
        "exit_sha256": sha256(output / "exit.json"),
        "stdout_sha256": sha256(output / "stdout.raw"),
        "stderr_sha256": sha256(output / "stderr.raw"),
        "source_sha256": PINS["source"],
        "worker_sha256": PINS["worker"],
        "checker_sha256": PINS["checker"],
        "launcher_sha256": launch["launcher_sha256"],
        "domains_sha256": PINS["domains"],
        "trace": trace,
        "trace_error": trace_error,
        "conservative_new_charge_seconds":
            launch["conservative_new_charge_seconds"],
        "cumulative_charged_seconds":
            PRIOR_CHARGE + launch["conservative_new_charge_seconds"],
        "natural_relations": None,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }


def run(mode, output, authorization=None):
    check_pins()
    require(not output.exists(), "output directory must be fresh")
    require(shutil.disk_usage(output.parent).free >= MIN_FREE_BYTES,
            "less than 256 MiB free before launch")
    auth_bytes = None
    if mode == "full_domain":
        require(authorization is not None, "full-domain run needs authorization")
        auth_bytes = authorization.read_bytes()
        auth = json.loads(auth_bytes)
        require(authorization_error(auth) is None,
                "full-domain authorization differs")
    else:
        require(mode == "self_test" and authorization is None,
                "self-test cannot carry authorization")
    output.mkdir()
    if auth_bytes is not None:
        with (output / "authorization.json").open("xb") as target:
            target.write(auth_bytes)
            target.flush()
            os.fsync(target.fileno())
    command = ([str(WORKER), "self-test-markers"] if mode == "self_test"
               else [str(WORKER), str(DOMAINS), str(CONFLICT_CAP)])
    launch = {
        "schema": "n83.fixed-witness-stage-launch/v1",
        "mode": mode,
        "scope": "known_answer_encoding_control_only",
        "logical_command": (["pinned_worker", "self-test-markers"]
                            if mode == "self_test" else
                            ["pinned_worker", "pinned_domains", str(CONFLICT_CAP)]),
        "pins": PINS,
        "launcher_sha256": sha256(Path(__file__)),
        "wall_cap_seconds": SELFTEST_WALL_CAP if mode == "self_test"
                            else FULL_WALL_CAP,
        "conflict_cap": None if mode == "self_test" else CONFLICT_CAP,
        "prior_charged_seconds": PRIOR_CHARGE,
        "conservative_new_charge_seconds":
            0 if mode == "self_test" else auth["new_authorized_seconds"],
        "authorization_sha256": None if auth_bytes is None else
            sha256(output / "authorization.json"),
    }
    save_json(output / "launch.json", launch)
    exit_row = launch_child(command, launch["wall_cap_seconds"], output)
    save_json(output / "exit.json", exit_row)
    row = audit(output)
    save_json(output / "audit.json", row)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--self-test", action="store_true")
    group.add_argument("--run", action="store_true")
    group.add_argument("--check", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--authorization", type=Path)
    args = parser.parse_args()
    if args.check is not None:
        require(args.output is None and args.authorization is None,
                "check mode takes only an evidence directory")
        row = audit(args.check)
        require((args.check / "audit.json").read_text() ==
                json.dumps(row, sort_keys=True, indent=2) + "\n",
                "saved audit differs")
    else:
        require(args.output is not None, "fresh output directory required")
        row = run("self_test" if args.self_test else "full_domain",
                  args.output, args.authorization)
    print(json.dumps({"status": row["status"], "mode": row["mode"],
                      "marker_count": (row["trace"] or {}).get("marker_count")},
                     sort_keys=True))


if __name__ == "__main__":
    main()

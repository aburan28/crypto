#!/usr/bin/env python3
"""One bounded graceful TERM/start check for an already unhealthy Docker backend."""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import signal
import socket
import subprocess
import time


PARENT_PID = 22317
PARENT_COMMAND = "/Applications/Docker.app/Contents/MacOS/com.docker.backend"
SOCKET = "/Users/adamburan/.docker/run/docker.sock"
OBSERVED_CHILDREN = (22323, 22324, 22357)
START_CAP = 180
POST_START_PROBES = 12


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError(reason)


def utc() -> str:
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def ping() -> str:
    client = socket.socket(socket.AF_UNIX)
    client.settimeout(5)
    try:
        client.connect(SOCKET)
        client.sendall(b"GET /_ping HTTP/1.1\r\nHost: docker\r\nConnection: close\r\n\r\n")
        return client.recv(2048).decode("utf-8", "replace")
    except OSError as error:
        return f"{type(error).__name__}: {error}"
    finally:
        client.close()


def process(pid: int) -> dict | None:
    result = subprocess.run(["ps", "-o", "pid=,ppid=,state=,comm=", "-p", str(pid)],
                            capture_output=True, text=True, check=False, timeout=5)
    if result.returncode != 0 or not result.stdout.strip():
        return None
    fields = result.stdout.strip().split(None, 3)
    require(len(fields) == 4 and int(fields[0]) == pid, "unexpected process record")
    return {"pid": pid, "ppid": int(fields[1]), "state": fields[2], "command": fields[3]}


def snapshot() -> dict:
    return {str(pid): process(pid) for pid in (PARENT_PID, *OBSERVED_CHILDREN)}


def preflight() -> dict:
    parent = process(PARENT_PID)
    response = ping()
    require(parent is not None and parent["command"] == PARENT_COMMAND
            and parent["ppid"] == 1,
            "pinned Docker backend parent changed; signal withheld")
    require("503 Service Unavailable" in response
            and "Docker Desktop is unable to start" in response,
            "Docker socket is not in the recorded failed state; signal withheld")
    return {"utc": utc(), "socket_ping": response, "processes": snapshot()}


def run(output: Path) -> dict:
    require(not output.exists(), "output directory must be fresh")
    before = preflight()
    output.mkdir(parents=True)
    result = {"schema": "n83.docker-backend-graceful-term/v1",
              "started_utc": utc(), "before": before,
              "signal": "SIGTERM", "signalled_pid": PARENT_PID,
              "sigkill_sent": False, "vm_disk_modified_directly": False,
              "start_command": None, "start_exit_code": None,
              "start_timed_out": None, "solver_search_executed": False,
              "selected_best_total_runtime": None}
    started = time.monotonic()
    try:
        os.kill(PARENT_PID, signal.SIGTERM)
        result["signal_sent_utc"] = utc()
        for _ in range(15):
            if process(PARENT_PID) is None:
                break
            time.sleep(2)
        result["after_term"] = {"utc": utc(), "socket_ping": ping(),
                                 "processes": snapshot()}
        if result["after_term"]["processes"][str(PARENT_PID)] is None:
            command = ["docker", "desktop", "start", "--timeout", str(START_CAP)]
            result["start_command"] = command
            try:
                child = subprocess.run(command, capture_output=True, timeout=START_CAP + 20)
                stdout, stderr = child.stdout, child.stderr
                result["start_exit_code"] = child.returncode
                result["start_timed_out"] = False
            except subprocess.TimeoutExpired as error:
                stdout, stderr = error.stdout or b"", error.stderr or b""
                result["start_exit_code"] = None
                result["start_timed_out"] = True
            (output / "start.stdout").write_bytes(stdout)
            (output / "start.stderr").write_bytes(stderr)
            result["start_stdout_sha256"] = hashlib.sha256(stdout).hexdigest()
            result["start_stderr_sha256"] = hashlib.sha256(stderr).hexdigest()
            result["after_start_immediate"] = ping()
            for _ in range(POST_START_PROBES):
                if "200 OK" in result.get("after_start_polled", result["after_start_immediate"]):
                    break
                time.sleep(5)
                result["after_start_polled"] = ping()
            result["status"] = ("PASS_DOCKER_DAEMON_RECOVERED"
                                if "200 OK" in result.get("after_start_polled",
                                                      result["after_start_immediate"])
                                else "PRODUCER_FAILURE_DOCKER_DAEMON_STILL_503")
        else:
            result["status"] = "PRODUCER_FAILURE_BACKEND_TERM_STILL_RUNNING"
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        result["status"] = "PRODUCER_FAILURE_GRACEFUL_RECOVERY"
        result["error"] = f"{type(error).__name__}: {error}"
    result["elapsed_wall_seconds_L0"] = time.monotonic() - started
    result["ended_utc"] = utc()
    (output / "result.json").write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    if args.check:
        print(json.dumps(preflight(), sort_keys=True))
    else:
        require(args.output_dir is not None, "output directory required")
        value = run(args.output_dir)
        print(json.dumps({"status": value["status"],
                          "elapsed_wall_seconds_L0": value["elapsed_wall_seconds_L0"]},
                         sort_keys=True))
        raise SystemExit(0 if value["status"] == "PASS_DOCKER_DAEMON_RECOVERED" else 1)

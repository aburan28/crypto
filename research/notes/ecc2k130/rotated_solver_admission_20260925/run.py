#!/usr/bin/env python3
"""Cold fail-preserving runner for the static solver-interface admission gate."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import resource
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
FILES = {name: HERE / name for name in ("INPUT.json", "PROTOCOL.md", "audit.py", "run.py", "ci_replay.py")}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def freeze() -> dict:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert {name: sha(path) for name, path in FILES.items()} == frozen["files"]
    data = json.loads((HERE / "INPUT.json").read_text())
    assert data["domain"] == "ECC2K130-ROTATED-M56-SOLVER-ADMISSION-20260925-v1"
    assert [arm["name"] for arm in data["arms"]] == ["n13-m5", "n19-m6"]
    assert data["caps"] == {"preflight_wall_seconds": 30, "preflight_rss_bytes": 256 * 1024 * 1024}
    return data


def child_rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    started = utc()
    receipt = {"started_utc": started, "python": sys.version,
               "platform": platform.platform(), "source_commit": None,
               "freeze_sha256": sha(HERE / "FROZEN.json"), "status": "started"}
    try:
        data = freeze()
        receipt["source_commit"] = data["base_commit"]
        receipt["source_sha256"] = {name: sha(path) for name, path in FILES.items()}
        command = [sys.executable, str(HERE / "audit.py"), "--out", str(args.out / "result.json")]
        receipt["argv"] = command
        receipt["child_started_utc"] = utc()
        with (args.out / "stdout.txt").open("w") as stdout, (args.out / "stderr.txt").open("w") as stderr:
            try:
                child = subprocess.run(command, stdout=stdout, stderr=stderr, check=False,
                                       timeout=data["caps"]["preflight_wall_seconds"] + 2)
                receipt["child_exit_code"] = child.returncode
            except subprocess.TimeoutExpired:
                receipt["child_exit_code"] = "TIMEOUT"
        receipt["child_finished_utc"] = utc()
        receipt["peak_child_rss_bytes_upper"] = child_rss_bytes()
        assert receipt["child_exit_code"] == 0
        assert receipt["peak_child_rss_bytes_upper"] <= data["caps"]["preflight_rss_bytes"]
        outcome = json.loads((args.out / "result.json").read_text())
        assert outcome["decision"] in {"BLOCKED_BEFORE_SOLVER_TIMING", "ELIGIBLE_FOR_SEPARATE_TIMED_PR"}
        assert bool(outcome["admitted_solver_arms"]) == (outcome["decision"] == "ELIGIBLE_FOR_SEPARATE_TIMED_PR")
        receipt["status"] = "complete_preflight"
    except BaseException as error:
        receipt["status"] = "failed"
        receipt["failure"] = repr(error)
    finally:
        receipt["finished_utc"] = utc()
        receipt["artifact_sha256"] = {path.name: sha(path) for path in args.out.iterdir()
                                      if path.is_file() and path.name != "receipt.json"}
        (args.out / "receipt.json").write_text(json.dumps(receipt, sort_keys=True, separators=(",", ":")) + "\n")
    return int(receipt["status"] != "complete_preflight")


if __name__ == "__main__":
    raise SystemExit(main())

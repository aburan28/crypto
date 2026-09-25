#!/usr/bin/env python3
"""Frozen staged runner; preserve every child status, stdout and stderr."""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

import psutil

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
SOURCE_NAMES = ("INPUT.json", "count.py", "replay.py", "run.py", "ci_replay.py", "PROTOCOL.md")
PARENT = HERE.parent / "rotated_subspace_support_20260925" / "gate.py"
VERIFIER_PARENT = HERE.parent / "rotated_pdp_corpus_20260925" / "verify.py"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def utc() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def check_frozen() -> dict:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    for name in SOURCE_NAMES:
        assert digest(HERE / name) == frozen["sha256"][name], name
    assert digest(PARENT) == frozen["sha256"]["parent_gate.py"]
    assert digest(VERIFIER_PARENT) == frozen["sha256"]["parent_pdp_verify.py"]
    data = json.loads((HERE / "INPUT.json").read_text())
    assert data["domain"] == frozen["domain"]
    assert data["pilot_masks"] == 1 << 15 and data["full_masks"] == 1 << 21
    assert len(data["sample_masks"]) == len(set(data["sample_masks"])) == 64
    return frozen


def child(argv: list[str], stem: str, out: Path, timeout: int) -> dict:
    start_utc = utc()
    start = time.perf_counter()
    stdout_path = out / f"{stem}.stdout.txt"
    stderr_path = out / f"{stem}.stderr.txt"
    peak = 0
    expired = False
    with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
        proc = subprocess.Popen(argv, cwd=REPO, stdout=stdout, stderr=stderr)
        monitor = psutil.Process(proc.pid)
        while proc.poll() is None:
            try:
                peak = max(peak, monitor.memory_info().rss)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
            if time.perf_counter() - start > timeout:
                expired = True
                proc.kill()
                break
            time.sleep(0.05)
        code = proc.wait()
    return {"argv": argv, "started_utc": start_utc,
            "finished_utc": utc(), "wall_seconds": time.perf_counter() - start,
            "exit_code": code, "runner_timeout": expired,
            "sampled_peak_rss_bytes": peak,
            "stdout_sha256": digest(stdout_path),
            "stderr_sha256": digest(stderr_path)}


def run(out: Path) -> None:
    out.mkdir(parents=True, exist_ok=False)
    try:
        frozen = check_frozen()
    except BaseException as error:
        write_json(out / "preflight_failure.json", {"error": repr(error), "utc": utc()})
        raise
    receipt = {"status": "started", "started_utc": utc(),
               "frozen_sha256": digest(HERE / "FROZEN.json"),
               "source_sha256": frozen["sha256"],
               "source_commit": subprocess.check_output(
                   ["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
               "python": sys.version, "platform": platform.platform(),
               "psutil": psutil.__version__, "commands": [], "raw_sha256": {}}

    def checkpoint() -> None:
        receipt["raw_sha256"] = {
            str(path.relative_to(out)): digest(path)
            for path in sorted(out.rglob("*")) if path.is_file() and path.name != "receipt.json"
        }
        write_json(out / "receipt.json", receipt)

    try:
        for stage in ("pilot", "full"):
            stage_dir = out / stage
            stage_cmd = [sys.executable, str(HERE / "count.py"), "--input",
                         str(HERE / "INPUT.json"), "--stage", stage,
                         "--out", str(stage_dir)]
            stage_result = child(stage_cmd, f"{stage}-producer", out,
                                 120 if stage == "pilot" else 1830)
            receipt["commands"].append(stage_result)
            checkpoint()
            if stage_result["exit_code"] != 0 or stage_result["runner_timeout"]:
                raise RuntimeError(f"{stage} producer failed or hit cap")
            summary = json.loads((stage_dir / "summary.json").read_text())
            if stage_result["sampled_peak_rss_bytes"] > 768 * 1024 * 1024:
                raise MemoryError(f"{stage} sampled RSS exceeded cap")
            replay_cmd = [sys.executable, str(HERE / "replay.py"),
                          "--input", str(HERE / "INPUT.json"),
                          "--summary", str(stage_dir / "summary.json"),
                          "--output", str(stage_dir / "replay.json")]
            replay_result = child(replay_cmd, f"{stage}-replay", out,
                                  180 if stage == "pilot" else 3630)
            receipt["commands"].append(replay_result)
            checkpoint()
            if replay_result["exit_code"] != 0 or replay_result["runner_timeout"]:
                raise RuntimeError(f"{stage} independent replay failed or hit cap")
            if replay_result["sampled_peak_rss_bytes"] > 768 * 1024 * 1024:
                raise MemoryError(f"{stage} replay sampled RSS exceeded cap")
            if stage == "pilot":
                receipt["pilot_projected_full_wall_seconds"] = summary["total_wall_seconds"] * 64
                checkpoint()
                if summary["total_wall_seconds"] * 64 > 1200:
                    receipt["status"] = "full_censored_by_preregistered_pilot_projection"
                    receipt["finished_utc"] = utc()
                    checkpoint()
                    return
                if summary["peak_rss_bytes"] > 256 * 1024 * 1024:
                    receipt["status"] = "full_censored_by_preregistered_pilot_rss"
                    receipt["finished_utc"] = utc()
                    checkpoint()
                    return
        receipt["status"] = "complete"
    except BaseException as error:
        receipt["status"] = "failed"
        receipt["error"] = repr(error)
        raise
    finally:
        receipt["finished_utc"] = utc()
        checkpoint()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    run(args.out)


if __name__ == "__main__":
    main()

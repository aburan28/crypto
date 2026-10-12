#!/usr/bin/env python3
"""Run all frozen arms under child caps and retain an auditable receipt."""
from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import tarfile
import time
from datetime import datetime, timezone
from pathlib import Path

from ci_replay import HERE, check_freeze, hashes, sha


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def child(name: str, command: list[str], directory: Path, timeout: int) -> dict:
    started = now()
    clock = time.perf_counter()
    stdout = directory / f"{name}.stdout.txt"
    stderr = directory / f"{name}.stderr.txt"
    with stdout.open("wb") as out, stderr.open("wb") as err:
        try:
            completed = subprocess.run(command, stdout=out, stderr=err, timeout=timeout)
            code, timed_out = completed.returncode, False
        except subprocess.TimeoutExpired:
            code, timed_out = None, True
    return {"name": name, "command": command, "started_utc": started,
            "ended_utc": now(), "external_timeout": timed_out,
            "external_timeout_seconds": timeout, "exit_code": code,
            "wall_seconds": time.perf_counter() - clock,
            "stdout_sha256": sha(stdout), "stderr_sha256": sha(stderr)}


def check_parent_merged() -> None:
    """Require the exact pinned #793 head in freshly fetched origin/main."""
    repo = HERE.parents[3]
    pinned = "6d8d1542b5513986fa0a94c0e0e7e9b1e46012c5"
    subprocess.run(["git", "fetch", "origin", "main"], cwd=repo, check=True,
                   timeout=90, capture_output=True)
    result = subprocess.run(["git", "merge-base", "--is-ancestor", pinned,
                             "refs/remotes/origin/main"], cwd=repo,
                            timeout=30, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError("exact frozen #793 head is not in freshly fetched origin/main")


def main(out: Path) -> None:
    check_freeze()
    check_parent_merged()
    assert not out.exists(), f"refusing to overwrite {out}"
    raw = out / "raw"
    raw.mkdir(parents=True)
    receipt = {"status": "running", "started_utc": now(),
               "python": sys.version, "platform": platform.platform(),
               "freeze_sha256": sha(HERE / "FROZEN.json"), "children": []}
    try:
        for position in range(6):
            folder = raw / f"p{position}"
            folder.mkdir()
            for stage, cap, args in [
                ("produce", 330, ["--out", str(folder)]),
                ("verify", 630, ["--data", str(folder), "--out", str(folder / "verify.json")]),
            ]:
                command = [sys.executable, str(HERE / (stage + ".py")),
                           "--position", str(position), *args]
                result = child(f"p{position}-{stage}", command, raw, cap)
                receipt["children"].append(result)
                if result["exit_code"] != 0 or result["external_timeout"]:
                    raise RuntimeError(f"{result['name']} failed: {result['exit_code']}")
        analysis = child("analysis", [sys.executable, str(HERE / "analyze.py"),
                                      "--raw", str(raw), "--out", str(raw / "analysis.json")], raw, 30)
        receipt["children"].append(analysis)
        if analysis["exit_code"] != 0 or analysis["external_timeout"]:
            raise RuntimeError(f"analysis failed: {analysis['exit_code']}")
        receipt["status"] = "success"
    except Exception as exc:
        receipt["status"] = "failure"
        receipt["error"] = repr(exc)
    finally:
        receipt["ended_utc"] = now()
        receipt["raw_sha256"] = hashes(raw)
        (out / "receipt.json").write_text(json.dumps(receipt, sort_keys=True, separators=(",", ":")) + "\n")
    if receipt["status"] != "success":
        raise RuntimeError(receipt["error"])
    archive = out / "raw.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(raw, arcname="raw")
    (out / "raw.tar.gz.sha256").write_text(sha(archive) + "\n")
    print(json.dumps({"status": "success", "archive_sha256": sha(archive),
                      "decision": json.loads((raw / "analysis.json").read_text())["position_gate"]},
                     sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    main(args.out)

#!/usr/bin/env python3
"""Run frozen exact toy/density panels with durable success or failure receipt."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE_FILES = ["gate.py", "verify.py", "run.py", "ci_replay.py"]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def tree_hashes(root: Path) -> dict[str, str]:
    return {str(p.relative_to(root)): sha(p) for p in sorted(root.rglob("*")) if p.is_file()}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    raw = args.out / "raw"
    raw.mkdir()
    inputs = HERE / "inputs"
    receipt = {"status": "started", "started_utc": utc(),
               "python": sys.version, "platform": platform.platform(),
               "source_sha256": None, "input_sha256": None, "commands": []}
    try:
        frozen = json.loads((HERE / "FROZEN.json").read_text())
        current_sources = {name: sha(HERE / name) for name in SOURCE_FILES}
        current_inputs = tree_hashes(inputs)
        receipt["source_sha256"] = current_sources
        receipt["input_sha256"] = current_inputs
        assert current_sources == frozen["source_sha256"], "source freeze mismatch"
        assert current_inputs == frozen["input_sha256"], "input freeze mismatch"
        commands = [
            ("toy", [sys.executable, str(HERE / "gate.py"), "toy", "--out", str(raw / "toy")], 1800),
            ("density", [sys.executable, str(HERE / "gate.py"), "density", "--out", str(raw / "density"), "--inputs", str(inputs)], 900),
            ("verify", [sys.executable, str(HERE / "verify.py"), "--archive", str(raw),
                        "--inputs", str(inputs), "--output", str(raw / "verify_report.json")], 1800),
        ]
        for name, command, timeout in commands:
            item = {"name": name, "argv": command, "started_utc": utc(), "timeout_seconds": timeout}
            receipt["commands"].append(item)
            with (args.out / f"{name}.stdout.txt").open("w") as stdout, (args.out / f"{name}.stderr.txt").open("w") as stderr:
                try:
                    result = subprocess.run(command, stdout=stdout, stderr=stderr, timeout=timeout, check=False)
                    item["exit_code"] = result.returncode
                except subprocess.TimeoutExpired:
                    item["exit_code"] = "TIMEOUT"
            item["finished_utc"] = utc()
            if item["exit_code"] != 0:
                raise RuntimeError(f"{name} failed: {item['exit_code']}")
        receipt["status"] = "success"
    except Exception as error:
        receipt["status"] = "failed"
        receipt["failure"] = repr(error)
    finally:
        receipt["finished_utc"] = utc()
        receipt["raw_sha256"] = tree_hashes(raw)
        (args.out / "receipt.json").write_text(json.dumps(receipt, sort_keys=True, separators=(",", ":")) + "\n")
    return 0 if receipt["status"] == "success" else 1


if __name__ == "__main__":
    raise SystemExit(main())

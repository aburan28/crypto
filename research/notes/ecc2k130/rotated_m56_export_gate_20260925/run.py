#!/usr/bin/env python3
"""Run the frozen semantic gate sequentially and preserve every child attempt."""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ARMS = ("n13-m5", "n19-m6")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, data) -> None:
    path.write_text(json.dumps(data, sort_keys=True, separators=(",", ":")) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert sha(HERE / "INPUTS.json") == frozen["inputs_sha256"]
    assert sha(HERE / "gate.py") == frozen["gate_sha256"]
    assert sha(Path(__file__)) == frozen["runner_sha256"]
    receipt = {"protocol": frozen["domain"], "arms": {}}
    for arm in ARMS:
        result = args.out / f"{arm}.json"
        command = [sys.executable, str(HERE / "gate.py"), "--arm", arm, "--out", str(result)]
        start = dt.datetime.now(dt.timezone.utc).isoformat()
        proc = subprocess.run(command, capture_output=True, text=True, check=False)
        end = dt.datetime.now(dt.timezone.utc).isoformat()
        stdout = args.out / f"{arm}.stdout.txt"
        stderr = args.out / f"{arm}.stderr.txt"
        stdout.write_text(proc.stdout)
        stderr.write_text(proc.stderr)
        receipt["arms"][arm] = {
            "command": command,
            "start_utc": start, "end_utc": end,
            "exit_code": proc.returncode,
            "result_sha256": sha(result) if result.exists() else None,
            "stdout_sha256": sha(stdout), "stderr_sha256": sha(stderr),
        }
        save(args.out / "receipt.json", receipt)
        if proc.returncode != 0:
            raise SystemExit(f"{arm} failed with exit {proc.returncode}; raw output preserved")
    print(json.dumps({"decision": "PASS", "receipt_sha256": sha(args.out / "receipt.json")}, sort_keys=True))


if __name__ == "__main__":
    main()

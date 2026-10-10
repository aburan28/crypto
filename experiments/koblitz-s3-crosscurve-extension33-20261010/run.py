"""Run frozen T025-T033 cross-curve ABBA blocks with the original binaries."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OLD = ROOT / "experiments/koblitz-s3-crosscurve-20261009"
FREEZE = HERE / "FREEZE.json"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def preflight(binaries: Path) -> dict:
    raw = FREEZE.read_bytes()
    frozen = json.loads(raw)
    assert frozen["kind"] == "s3_crosscurve_extension33_freeze_v1"
    assert frozen["runner_sha256"] == sha(Path(__file__).read_bytes())
    assert frozen["old_runner_sha256"] == sha((OLD / "run_local_pilot.py").read_bytes())
    assert frozen["old_freeze_sha256"] == sha((OLD / "FREEZE.json").read_bytes())
    assert frozen["selected_indices"] == list(range(25, 34))
    assert frozen["planned_run_count"] == 72
    assert frozen["workers"] == 14 and frozen["per_run_timeout_seconds"] == 180
    for arm, name in (("reference", "s3_pair_candidate_frozen"),
                      ("prepared", "s3_partial_prepared_targetspan")):
        binary = binaries / name
        assert sha(binary.read_bytes()) == frozen["binary_sha256"][arm]
    return frozen


def run(binaries: Path, output: Path) -> None:
    frozen = preflight(binaries)
    if output.exists():
        raise FileExistsError(output)
    command = [sys.executable, str(OLD / "run_local_pilot.py"),
               "--first", "25", "--count", "9", "--binaries", str(binaries),
               "--output", str(output)]
    process = subprocess.run(command, check=False)
    ledger = output / "raw.jsonl"
    rows = [json.loads(line) for line in ledger.read_bytes().splitlines()] if ledger.exists() else []
    expected = [(row["n"], row["label"], row["arm"], row["run_id"])
                for row in frozen["run_order"]]
    actual = [(row["n"], row["label"], row["arm"], row["run_id"]) for row in rows]
    order_ok = actual == expected[:len(actual)]
    complete = (process.returncode == 0 and len(rows) == len(expected)
                and order_ok and all(row["verified_native_and_fixture"] for row in rows)
                and all(row["binary_sha256"] == frozen["binary_sha256"][row["arm"]]
                        for row in rows))
    receipt = {
        "kind": "s3_crosscurve_extension33_completion_v1",
        "freeze_sha256": sha(FREEZE.read_bytes()),
        "command": command,
        "runner_exit_code": process.returncode,
        "planned_run_count": len(expected),
        "completed_run_count": len(rows),
        "verified_run_count": sum(row["verified_native_and_fixture"] for row in rows),
        "order_prefix_verified": order_ok,
        "ledger_sha256": sha(ledger.read_bytes()) if ledger.exists() else None,
        "complete": complete,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "completion.json").write_text(json.dumps(receipt, indent=2,
                                                       sort_keys=True) + "\n")
    print(json.dumps(receipt, sort_keys=True), flush=True)
    if not complete:
        raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binaries", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(args.binaries.resolve(), args.output.resolve())

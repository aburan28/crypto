#!/usr/bin/env python3
"""Fail-preserving cold runner for the frozen phase-census children."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import resource
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
ROOT = HERE.parents[3]
FILES = {
    "input": HERE / "INPUT.json",
    "protocol": HERE / "PROTOCOL.md",
    "producer": HERE / "phase.py",
    "independent_replay": HERE / "verify.py",
    "runner": HERE / "run.py",
    "ci_replay": HERE / "ci_replay.py",
    "workflow": ROOT / ".github/workflows/ecc2k130-rotated-phase-selector.yml",
    "source_archive": NOTES / "rotated_pdp_corpus_20260925/evidence/raw.tar.gz",
    "producer_arithmetic": NOTES / "rotated_subspace_support_20260925/gate.py",
    "independent_arithmetic": NOTES / "rotated_pdp_corpus_20260925/verify.py",
    "independent_field": NOTES / "rotated_subspace_support_20260925/verify.py",
    "point_only": NOTES / "rotated_four_base_joint_rank_20260925/inputs/point_only.json",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def child_cpu() -> float:
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    return usage.ru_utime + usage.ru_stime


def child_peak_rss() -> int:
    value = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def freeze() -> dict:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert frozen["schema"] == "ecc2k130_rotated_phase_selector_freeze_v1"
    for label, path in FILES.items():
        assert path.is_file() and sha(path) == frozen["file_sha256"][label], label
    manifest = json.loads((HERE / "INPUT.json").read_text())
    assert manifest["schema"] == "ecc2k130_rotated_phase_selector_input_v1"
    assert manifest["q"] == 130873 and manifest["lambda"] == 41811
    assert manifest["field_polynomial"] == 0x80027
    assert manifest["phase_orders"] == {
        "spread": [(5 * j) % 19 for j in range(19)],
        "contiguous_control": list(range(19)),
    }
    assert manifest["caps"] == [1, 4, 8, 19]
    assert manifest["comparators"] == {
        "four_beta_cap4": {"support": 121851, "oracle_probes": 249827},
        "five_beta_cap19": {"support": 126265, "oracle_probes": 258849},
    }
    assert manifest["phase_producer_wall_seconds"] == 120
    assert manifest["independent_replay_wall_seconds"] == 300
    assert manifest["child_rss_cap_bytes"] == 512 * 1024 * 1024
    return frozen


def file_hashes(root: Path) -> dict[str, str]:
    return {str(path.relative_to(root)): sha(path)
            for path in sorted(root.rglob("*")) if path.is_file()
            and path.name != "receipt.json"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    receipt = {"status": "started", "started_utc": now(),
               "platform": platform.platform(), "python": sys.version,
               "commands": []}
    try:
        frozen = freeze()
        manifest = json.loads((HERE / "INPUT.json").read_text())
        receipt["freeze_sha256"] = sha(HERE / "FROZEN.json")
        receipt["file_sha256"] = frozen["file_sha256"]
        commands = [
            ("producer", [sys.executable, str(HERE / "phase.py"), "--out",
                           str(args.out / "producer")],
             manifest["phase_producer_wall_seconds"]),
            ("independent_replay", [sys.executable, str(HERE / "verify.py"),
                                    "--evidence", str(args.out / "producer"),
                                    "--out", str(args.out / "verification.json")],
             manifest["independent_replay_wall_seconds"]),
        ]
        for name, argv, timeout in commands:
            item = {"name": name, "argv": argv, "started_utc": now(),
                    "timeout_seconds": timeout}
            receipt["commands"].append(item)
            start = time.monotonic()
            cpu_start = child_cpu()
            with (args.out / f"{name}.stdout.txt").open("w") as stdout, \
                    (args.out / f"{name}.stderr.txt").open("w") as stderr:
                try:
                    process = subprocess.run(argv, stdout=stdout, stderr=stderr,
                                             timeout=timeout, check=False)
                    item["exit_code"] = process.returncode
                except subprocess.TimeoutExpired:
                    item["exit_code"] = "TIMEOUT"
            item["finished_utc"] = now()
            item["wall_seconds"] = time.monotonic() - start
            item["child_cpu_seconds"] = child_cpu() - cpu_start
            item["high_water_child_rss_bytes"] = child_peak_rss()
            item["stdout_sha256"] = sha(args.out / f"{name}.stdout.txt")
            item["stderr_sha256"] = sha(args.out / f"{name}.stderr.txt")
            if item["exit_code"] != 0:
                raise RuntimeError(f"{name} failed or timed out: {item['exit_code']}")
            if item["high_water_child_rss_bytes"] > manifest["child_rss_cap_bytes"]:
                raise RuntimeError(f"{name} exceeded frozen RSS cap")
        verified = json.loads((args.out / "verification.json").read_text())
        assert verified["status"] == "PASS" and verified["q_checked"] == 130873
        receipt["status"] = "success"
    except Exception as error:
        receipt["status"] = "failed"
        receipt["failure"] = repr(error)
    finally:
        receipt["finished_utc"] = now()
        receipt["result_sha256"] = file_hashes(args.out)
        (args.out / "receipt.json").write_text(
            json.dumps(receipt, sort_keys=True, separators=(",", ":")) + "\n")
    return int(receipt["status"] != "success")


if __name__ == "__main__":
    raise SystemExit(main())

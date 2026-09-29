#!/usr/bin/env python3
"""Fail-preserving runner for the preregistered five-base archive analysis."""
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
NOTES = HERE.parent
CORPUS = NOTES / "rotated_pdp_corpus_20260925/evidence/raw.tar.gz"
SWEEP = NOTES / "rotated_beta_sweep_20260925/evidence/raw.tar.gz"
INDEPENDENT = NOTES / "rotated_pdp_corpus_20260925/verify.py"
FILES = {"corpus_archive": CORPUS, "sweep_archive": SWEEP,
         "independent_curve": INDEPENDENT,
         **{name.replace(".py", ""): HERE / name
            for name in ("analyze.py", "verify.py", "run.py", "ci_replay.py")},
         "manifest": HERE / "input_manifest.json", "protocol": HERE / "PROTOCOL.md"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def peak_rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def freeze() -> dict:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    for label, path in FILES.items():
        assert sha(path) == frozen[label + "_sha256"], label
    manifest = json.loads((HERE / "input_manifest.json").read_text())
    assert manifest["domain"] == "ECC2K130-ROTATED-BETA-PORTFOLIO-20260925-v1"
    assert manifest["beta_order"] == [3, 338435, 303097, 464276, 42605]
    assert manifest["q"] == 130873 and manifest["m"] == 6
    assert manifest["caps"] == {"rss_bytes": 512 * 1024 * 1024, "wall_seconds": 120}
    return frozen


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    frozen = freeze()
    receipt = {"status": "started", "started_utc": utc(),
               "platform": platform.platform(), "python": sys.version,
               "freeze_sha256": sha(HERE / "FROZEN.json"),
               "input_sha256": {label: frozen[label + "_sha256"] for label in FILES},
               "commands": []}
    manifest = json.loads((HERE / "input_manifest.json").read_text())
    cap = manifest["caps"]
    commands = [("analysis", [sys.executable, str(HERE / "analyze.py"),
                              "--out", str(args.out / "outcome.json")]),
                ("independent_verify", [sys.executable, str(HERE / "verify.py"),
                                        "--expected", str(args.out / "outcome.json"),
                                        "--out", str(args.out / "verification.json")])]
    try:
        for name, argv in commands:
            item = {"name": name, "argv": argv, "started_utc": utc(),
                    "timeout_seconds": cap["wall_seconds"]}
            receipt["commands"].append(item)
            with (args.out / f"{name}.stdout.txt").open("w") as stdout, \
                    (args.out / f"{name}.stderr.txt").open("w") as stderr:
                try:
                    process = subprocess.run(argv, stdout=stdout, stderr=stderr,
                                             timeout=cap["wall_seconds"], check=False)
                    item["exit_code"] = process.returncode
                except subprocess.TimeoutExpired:
                    item["exit_code"] = "TIMEOUT"
            item["finished_utc"] = utc()
            item["peak_child_rss_bytes"] = peak_rss_bytes()
            if item["exit_code"] != 0 or item["peak_child_rss_bytes"] > cap["rss_bytes"]:
                raise RuntimeError(f"{name}: exit/RSS cap: {item['exit_code']}, {item['peak_child_rss_bytes']}")
        receipt["status"] = "success"
    except Exception as error:
        receipt["status"] = "failed"
        receipt["failure"] = repr(error)
    finally:
        receipt["finished_utc"] = utc()
        receipt["result_sha256"] = {path.name: sha(path) for path in args.out.iterdir()
                                    if path.is_file() and path.name != "receipt.json"}
        (args.out / "receipt.json").write_text(
            json.dumps(receipt, sort_keys=True, separators=(",", ":")) + "\n")
    return int(receipt["status"] != "success")


if __name__ == "__main__":
    raise SystemExit(main())

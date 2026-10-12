#!/usr/bin/env python3
"""Run the frozen two-arm corpus, preserving every child/failure receipt."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / "rotated_subspace_support_20260925"
SOURCES = ["producer.py", "verify.py", "run.py", "ci_replay.py"]
PARENT_SOURCES = ["gate.py", "verify.py"]
PARENT_ARCHIVE_SHA = "fad4b4cc48f980349d1ef472d18327368dd0f93f3d1284fc8bbdc20f0d3086d7"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def tree_hashes(path: Path) -> dict[str, str]:
    return {str(p.relative_to(path)): sha(p) for p in sorted(path.rglob("*")) if p.is_file()}


def check_freeze() -> dict:
    protocol = (HERE / "PROTOCOL.md").read_text()
    match = re.search(r'corpus `FROZEN\.json` SHA-256 `([0-9a-f]{64})`', protocol)
    assert match is not None and sha(HERE / "FROZEN.json") == match.group(1), "protocol freeze anchor"
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    manifest = json.loads((HERE / "input_manifest.json").read_text())
    sources = {name: sha(HERE / name) for name in SOURCES}
    parent_sources = {name: sha(PARENT / name) for name in PARENT_SOURCES}
    assert sources == frozen["source_sha256"], "corpus source drift"
    assert parent_sources == frozen["parent_source_sha256"], "parent source drift"
    assert sha(HERE / "input_manifest.json") == frozen["input_manifest_sha256"]
    assert manifest["domain"] == "ECC2K130-ROTATED-PDP-CORPUS-20260925-v1"
    assert [(a["name"], a["n"], a["poly"], a["q"], a["beta"], a["m"], a["d"])
            for a in manifest["arms"]] == [
                ("n13-m5", 13, 0x201b, 2003, 3, 5, 2),
                ("n19-m6", 19, 0x80027, 130873, 3, 6, 2)]
    assert manifest["parent_raw_archive_sha256"] == PARENT_ARCHIVE_SHA
    assert sha(PARENT / "evidence" / "raw.tar.gz") == PARENT_ARCHIVE_SHA
    return {"source_sha256": sources, "parent_source_sha256": parent_sources,
            "input_manifest_sha256": frozen["input_manifest_sha256"],
            "parent_archive_sha256": PARENT_ARCHIVE_SHA}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    raw = args.out / "raw"
    raw.mkdir()
    receipt = {"status": "started", "started_utc": utc(), "python": sys.version,
               "platform": platform.platform(), "commands": [], "freeze": None}
    try:
        receipt["freeze"] = check_freeze()
        for arm in ("n13-m5", "n19-m6"):
            producer_out = raw / arm
            verifier_out = raw / f"{arm}-verify.json"
            commands = [
                (f"{arm}-producer", [sys.executable, str(HERE / "producer.py"),
                                      "--arm", arm, "--out", str(producer_out)], 330),
                (f"{arm}-verify", [sys.executable, str(HERE / "verify.py"),
                                    "--arm", arm, "--archive", str(producer_out),
                                    "--out", str(verifier_out)], 630),
            ]
            for name, argv, timeout in commands:
                item = {"name": name, "argv": argv, "timeout_seconds": timeout,
                        "started_utc": utc()}
                receipt["commands"].append(item)
                with (args.out / f"{name}.stdout.txt").open("w") as stdout, \
                        (args.out / f"{name}.stderr.txt").open("w") as stderr:
                    try:
                        result = subprocess.run(argv, stdout=stdout, stderr=stderr,
                                                timeout=timeout, check=False)
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
        (args.out / "receipt.json").write_text(
            json.dumps(receipt, sort_keys=True, separators=(",", ":")) + "\n")
    return int(receipt["status"] != "success")


if __name__ == "__main__":
    raise SystemExit(main())

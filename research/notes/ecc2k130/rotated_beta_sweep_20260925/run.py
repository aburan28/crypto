#!/usr/bin/env python3
"""Frozen, fail-preserving four-arm n19 normal-beta support sweep runner."""
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
CORPUS = HERE.parent / "rotated_pdp_corpus_20260925"
PARENT = HERE.parent / "rotated_subspace_support_20260925"
SOURCES = ["preflight.py", "producer.py", "verify.py", "run.py", "ci_replay.py"]
REFERENCE_SOURCES = ["producer.py", "verify.py", "run.py"]
PARENT_SOURCES = ["gate.py", "verify.py"]
REFERENCE_ARCHIVE_SHA = "39f16990213e2c7525c6dae92c131127b12f914bbfd6182264249da72bc9546c"
REFERENCE_MERGE = "c77767c4a653734f428e110cca29985d721476b2"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def tree_hashes(path: Path) -> dict[str, str]:
    return {str(p.relative_to(path)): sha(p) for p in sorted(path.rglob("*")) if p.is_file()}


def check_freeze() -> dict:
    protocol = (HERE / "PROTOCOL.md").read_text()
    match = re.search(r'sweep `FROZEN\.json` SHA-256 `([0-9a-f]{64})`', protocol)
    assert match is not None and sha(HERE / "FROZEN.json") == match.group(1), "protocol freeze anchor"
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    manifest = json.loads((HERE / "input_manifest.json").read_text())
    selection = json.loads((HERE / "selection.json").read_text())
    sources = {name: sha(HERE / name) for name in SOURCES}
    reference_sources = {name: sha(CORPUS / name) for name in REFERENCE_SOURCES}
    parent_sources = {name: sha(PARENT / name) for name in PARENT_SOURCES}
    assert sources == frozen["source_sha256"], "sweep source drift"
    assert reference_sources == frozen["reference_source_sha256"], "#767 source drift"
    assert parent_sources == frozen["parent_source_sha256"], "#762 source drift"
    assert sha(HERE / "input_manifest.json") == frozen["input_manifest_sha256"], "manifest drift"
    assert sha(HERE / "selection.json") == frozen["selection_sha256"], "selection drift"
    assert sha(CORPUS / "evidence" / "raw.tar.gz") == REFERENCE_ARCHIVE_SHA, "#767 archive drift"
    assert manifest["domain"] == selection["domain"] == "ECC2K130-ROTATED-BETA-SWEEP-20260925-v1"
    assert manifest["reference_merge"] == REFERENCE_MERGE
    assert manifest["reference_archive_sha256"] == REFERENCE_ARCHIVE_SHA
    assert manifest["n"] == 19 and manifest["poly"] == 0x80027
    assert manifest["q"] == 130873 and manifest["m"] == 6 and manifest["d"] == 2
    assert manifest["reference_beta"] == 3 and manifest["reference_H"] == [385982, 301867]
    assert manifest["caps"] == {"preflight_wall_seconds": 60, "preflight_rss_bytes": 128 * 1024 * 1024,
                                "producer_wall_seconds": 300, "producer_rss_bytes": 512 * 1024 * 1024,
                                "verifier_wall_seconds": 600, "verifier_rss_bytes": 512 * 1024 * 1024}
    assert selection["status"] == "success"
    assert manifest["selected_betas"] == [row["beta"] for row in selection["selected"]]
    assert 3 <= len(selection["selected"]) <= 4
    assert selection["reference_f0_size"] == 7
    assert selection["reference_projected_signed"] == 7
    assert all(row["selection"] == "primary" for row in selection["selected"]) if not selection["fallback_used"] else True
    return {"source_sha256": sources, "reference_source_sha256": reference_sources,
            "parent_source_sha256": parent_sources,
            "input_manifest_sha256": frozen["input_manifest_sha256"],
            "selection_sha256": frozen["selection_sha256"],
            "reference_archive_sha256": REFERENCE_ARCHIVE_SHA,
            "reference_merge": REFERENCE_MERGE}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    raw = args.out / "raw"
    raw.mkdir()
    receipt = {"status": "started", "started_utc": utc(),
               "python": sys.version, "platform": platform.platform(),
               "commands": [], "freeze": None}
    try:
        receipt["freeze"] = check_freeze()
        selection = json.loads((HERE / "selection.json").read_text())
        commands = [("selection-verify", [sys.executable, str(HERE / "verify.py"),
                                          "--selection-only"], 90)]
        for row in selection["selected"]:
            beta = row["beta"]
            out = raw / f"beta-{beta}"
            commands.extend([
                (f"beta-{beta}-producer", [sys.executable, str(HERE / "producer.py"),
                                            "--beta", str(beta), "--out", str(out)], 330),
                (f"beta-{beta}-verify", [sys.executable, str(HERE / "verify.py"),
                                          "--beta", str(beta), "--archive", str(out),
                                          "--out", str(raw / f"beta-{beta}-verify.json")], 630),
            ])
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

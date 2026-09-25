#!/usr/bin/env python3
"""Check frozen bytes; with evidence, replay every archived semantic claim."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
DOMAIN = "k0-symbolic-dag-fullpoint-n2n3-v1"
SOURCES = ("PROTOCOL.md", "dag.py", "produce.py", "verify.py", "run.py", "ci_replay.py")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_freeze() -> dict:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert frozen["domain"] == DOMAIN
    assert frozen["curve"] == "y^2+xy=x^3+1"
    assert frozen["fields"] == [[2, 7], [3, 11]]
    assert frozen["model_width_formula"] == "7*n+3"
    assert frozen["child_wall_cap_seconds"] == 180
    assert frozen["external_child_cap_seconds"] == 195
    assert frozen["child_rss_cap_bytes"] == 512 * 1024 * 1024
    assert set(frozen["source_sha256"]) == set(SOURCES)
    assert sha(HERE.parents[3] / ".github/workflows/ecc2k130-symbolic-dag-fullpoint.yml") == frozen["workflow_sha256"]
    for name in SOURCES:
        assert sha(HERE / name) == frozen["source_sha256"][name], name
    return frozen


def check_evidence(path: Path, frozen: dict) -> dict:
    receipt = json.loads(path.read_text())
    assert receipt["domain"] == DOMAIN
    assert receipt["freeze_sha256"] == sha(HERE / "FROZEN.json")
    assert receipt["decision"] == "PASS"
    attempts = receipt["attempts"]
    assert [item["phase"] for item in attempts] == ["producer", "verifier"]
    archive = path.parent
    for item in attempts:
        phase = item["phase"]
        assert item["exit_code"] == 0 and item["external_timeout"] is False
        assert item["wall_seconds"] <= frozen["external_child_cap_seconds"]
        assert sha(archive / f"{phase}.stdout.txt") == item["stdout_sha256"]
        assert sha(archive / f"{phase}.stderr.txt") == item["stderr_sha256"]
        result_file = archive / ("producer/result.json" if phase == "producer" else "verify.json")
        assert sha(result_file) == item["result_sha256"]
        result = json.loads(result_file.read_text())
        assert result["decision"] == "PASS"
        assert result["wall_seconds"] <= frozen["child_wall_cap_seconds"]
        assert result["peak_rss_bytes"] <= frozen["child_rss_cap_bytes"]
        assert item["child_wall_seconds"] == result["wall_seconds"]
        assert item["child_cpu_seconds"] == result["cpu_seconds"]
        assert item["child_peak_rss_bytes"] == result["peak_rss_bytes"]
    producer = archive / "producer"
    produced = json.loads((producer / "result.json").read_text())
    verified = json.loads((archive / "verify.json").read_text())
    assert receipt["producer_sha256"] == sha(producer / "result.json")
    assert receipt["verifier_sha256"] == sha(archive / "verify.json")
    assert produced["rows_sha256"] == sha(producer / "rows.jsonl.gz")
    assert verified["producer_sha256"] == receipt["producer_sha256"]
    assert verified["rows_sha256"] == produced["rows_sha256"]
    with tempfile.TemporaryDirectory() as temporary:
        replay_file = Path(temporary) / "verify.json"
        subprocess.run([sys.executable, str(HERE / "verify.py"), "--producer",
                        str(producer), "--out", str(replay_file)],
                       cwd=HERE, check=True, capture_output=True, text=True,
                       timeout=frozen["external_child_cap_seconds"])
        replay = json.loads(replay_file.read_text())
    for key in ("domain", "decision", "fields", "producer_sha256", "rows_sha256"):
        assert replay[key] == verified[key], key
    return {"decision": "PASS", "fields": replay["fields"],
            "archive_verifier_sha256": receipt["verifier_sha256"]}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    frozen = check_freeze()
    evidence = check_evidence(args.evidence, frozen) if args.evidence else None
    print(json.dumps({"decision": "PASS", "freeze_sha256": sha(HERE / "FROZEN.json"),
                      "evidence": evidence}, sort_keys=True))


if __name__ == "__main__":
    main()

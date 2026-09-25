#!/usr/bin/env python3
"""Check committed evidence bytes, child receipts, and independently replayed labels."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
EVIDENCE = HERE / "evidence"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest = json.loads((EVIDENCE / "manifest.json").read_text())
    assert manifest["format"] == "ECC2K130-ROTATED-S3-EVIDENCE-v1"
    expected = manifest["files"]
    actual = {str(path.relative_to(EVIDENCE)) for path in EVIDENCE.rglob("*")
              if path.is_file() and path.name != "manifest.json"}
    assert actual == set(expected)
    for name, record in expected.items():
        path = EVIDENCE / name
        assert sha(path) == record["sha256"] and path.stat().st_size == record["bytes"]
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    receipt = json.loads((EVIDENCE / "receipt.json").read_text())
    assert receipt["protocol"] == frozen["domain"]
    assert receipt["freeze_sha256"] == sha(HERE / "FROZEN.json")
    assert receipt["decision"] == "PASS"
    assert [(row["arm"], row["phase"]) for row in receipt["attempts"]] == [
        ("n13-m5", "producer"), ("n13-m5", "verifier"),
        ("n19-m6", "producer"), ("n19-m6", "verifier")]
    for row in receipt["attempts"]:
        arm, phase = row["arm"], row["phase"]
        assert row["exit_code"] == 0 and not row["external_timeout"]
        result = EVIDENCE / arm / ("producer/result.json" if phase == "producer" else "verify.json")
        assert sha(result) == row["result_sha256"]
        assert sha(EVIDENCE / f"{arm}-{phase}.stdout.txt") == row["stdout_sha256"]
        assert sha(EVIDENCE / f"{arm}-{phase}.stderr.txt") == row["stderr_sha256"]
        if phase == "verifier":
            verifier = json.loads(result.read_text())
            assert verifier["decision"] == "PASS"
            assert verifier["producer_sha256"] == sha(EVIDENCE / arm / "producer/result.json")
    print(json.dumps({"decision": "PASS", "files": len(expected),
                      "manifest_sha256": sha(EVIDENCE / "manifest.json")}, sort_keys=True))


if __name__ == "__main__":
    main()

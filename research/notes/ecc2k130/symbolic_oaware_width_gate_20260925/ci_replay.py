#!/usr/bin/env python3
"""Frozen source checks and independent toy verifier replay from archives."""

import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

from audit import HERE, canonical, preflight, width_decision


def sha_file(path: Path) -> dict:
    raw = path.read_bytes()
    return {"sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}


def main() -> None:
    frozen, inp = preflight()
    evidence = HERE / "evidence"
    receipt_path = evidence / "receipt.json"
    if not receipt_path.exists():
        assert not evidence.exists() or not any(evidence.iterdir()), "partial outcome without receipt"
        print("FROZEN_HASHES_PASS; OUTCOME_HELD")
        return
    raw = receipt_path.read_bytes()
    assert len(raw) <= inp["caps"]["receipt_bytes"]
    r = json.loads(raw)
    assert raw == canonical(r)
    assert r["schema"] == "symbolic-oaware-width-gate-run-v1"
    assert r["status"] == "PASS", r["reason"]
    assert r["frozen_sha256"] == sha_file(HERE / "FROZEN.json")["sha256"]
    assert r["source_commit"] == frozen["parent_commit"]
    assert r["wall_seconds"] <= inp["caps"]["wall_seconds"]
    assert r["rss_cap_bytes"] == inp["caps"]["rss_bytes"]
    assert len(r["children"]) == 2
    for c in r["children"]:
        assert c["exit_code"] == 0 and not c["timed_out"]
        assert c["children_ru_maxrss_bytes_upper"] <= inp["caps"]["rss_bytes"]
        for channel in ("stdout", "stderr"):
            assert sha_file(evidence / c[channel]["path"]) == {k: c[channel][k] for k in ("sha256", "bytes")}
    assert set(r["files"]) == {p.name for p in evidence.iterdir() if p.is_file()} - {"receipt.json"}
    for name, meta in r["files"].items():
        assert sha_file(evidence / name) == meta
    producer_path, verify_path = evidence / "producer.json", evidence / "verify.json"
    producer = json.loads(producer_path.read_text())
    independent = json.loads(verify_path.read_text())
    assert producer["schema"] == "symbolic-oaware-width-gate-receipt-v1"
    assert producer["width"] == width_decision(inp)
    assert producer["toy_mu4_map"] == independent["toy_mu4_map"]
    assert producer["toy_mu4_map"]["status"] == "PASS"
    assert producer["cost"]["wall_seconds"] <= inp["caps"]["wall_seconds"]
    assert producer["cost"]["peak_rss_bytes"] <= inp["caps"]["rss_bytes"]
    with tempfile.TemporaryDirectory() as tmp:
        replay = Path(tmp) / "verify.json"
        cp = subprocess.run([sys.executable, str(HERE / "verify.py"),
                             "--producer", str(producer_path), "--output", str(replay)],
                            cwd=HERE, capture_output=True, timeout=inp["caps"]["wall_seconds"])
        assert cp.returncode == 0, cp.stderr.decode(errors="replace")
        assert replay.read_bytes() == verify_path.read_bytes()
    print("INDEPENDENT_ARCHIVE_REPLAY_PASS", hashlib.sha256(raw).hexdigest())


if __name__ == "__main__":
    main()

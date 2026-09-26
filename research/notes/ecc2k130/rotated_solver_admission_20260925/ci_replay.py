#!/usr/bin/env python3
"""Hash-only preregistration, then lightweight committed blocker replay."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from audit import audit
from run import HERE, FILES, freeze, sha


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    data = freeze()
    for name in ("audit.py", "run.py", "ci_replay.py"):
        compile((HERE / name).read_bytes(), name, "exec")
    for item in data["reference_files"].values():
        assert sha(HERE.parents[3] / item["path"]) == item["sha256"]
    if args.evidence is None:
        print("Frozen source, corpus, protocol and input hashes PASS; no solver or admission outcome read.")
        return
    receipt = json.loads((args.evidence / "receipt.json").read_text())
    assert receipt["status"] == "complete_preflight" and receipt["child_exit_code"] == 0
    assert receipt["freeze_sha256"] == sha(HERE / "FROZEN.json")
    assert receipt["source_sha256"] == {name: sha(path) for name, path in FILES.items()}
    assert receipt["peak_child_rss_bytes_upper"] <= data["caps"]["preflight_rss_bytes"]
    for name, digest in receipt["artifact_sha256"].items():
        assert sha(args.evidence / name) == digest
    assert {path.name for path in args.evidence.iterdir() if path.is_file()} == (
        set(receipt["artifact_sha256"]) | {"receipt.json"})
    recorded = json.loads((args.evidence / "result.json").read_text())
    fresh = audit(data, include_binaries=False)
    for key, value in fresh.items():
        assert recorded[key] == value, key
    assert recorded["wall_seconds"] <= data["caps"]["preflight_wall_seconds"]
    assert recorded["peak_rss_bytes"] <= data["caps"]["preflight_rss_bytes"]
    assert recorded["decision"] == "BLOCKED_BEFORE_SOLVER_TIMING"
    assert recorded["admitted_solver_arms"] == []
    print("Frozen corpus/interface blocker and raw receipt replay PASS; zero solver outcomes claimed.")


if __name__ == "__main__":
    main()

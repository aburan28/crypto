#!/usr/bin/env python3
"""Hash-only freeze preflight; archived outcome gets independent replay."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

from run import HERE, file_hashes, freeze, sha


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    frozen = freeze()
    for name in ("phase.py", "verify.py", "run.py", "ci_replay.py"):
        source = HERE / name
        compile(source.read_bytes(), str(source), "exec")
    if args.evidence is None:
        print("Frozen source/input/archive SHA-256 and syntax PASS; no phase outcome read.")
        return
    evidence = args.evidence.resolve()
    receipt = json.loads((evidence / "receipt.json").read_text())
    assert receipt["status"] == "success"
    assert receipt["freeze_sha256"] == sha(HERE / "FROZEN.json")
    assert receipt["file_sha256"] == frozen["file_sha256"]
    assert receipt["result_sha256"] == file_hashes(evidence)
    assert [item["name"] for item in receipt["commands"]] == [
        "producer", "independent_replay"]
    caps = json.loads((HERE / "INPUT.json").read_text())
    assert all(item["exit_code"] == 0
               and item["high_water_child_rss_bytes"] <= caps["child_rss_cap_bytes"]
               and item["stdout_sha256"] == sha(evidence / f"{item['name']}.stdout.txt")
               and item["stderr_sha256"] == sha(evidence / f"{item['name']}.stderr.txt")
               for item in receipt["commands"])
    archived = json.loads((evidence / "verification.json").read_text())
    assert archived["status"] == "PASS"
    with tempfile.TemporaryDirectory() as directory:
        fresh = Path(directory) / "verification.json"
        subprocess.run([sys.executable, str(HERE / "verify.py"), "--evidence",
                        str(evidence / "producer"), "--out", str(fresh)],
                       check=True, timeout=caps["independent_replay_wall_seconds"])
        assert json.loads(fresh.read_text()) == archived
    print("Independent all-q phase arrays and fixed-case row replay PASS.")


if __name__ == "__main__":
    main()

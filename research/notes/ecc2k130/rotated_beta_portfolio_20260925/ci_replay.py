#!/usr/bin/env python3
"""Hash-only PR preflight, then independent replay of committed evidence."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

from run import HERE, freeze, sha


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    frozen = freeze()
    for name in ("analyze.py", "verify.py", "run.py", "ci_replay.py"):
        source = HERE / name
        compile(source.read_bytes(), str(source), "exec")
    if args.evidence is None:
        print("Frozen protocol/code/reference SHA-256 and syntax PASS; no portfolio outcome read.")
        return
    receipt = json.loads((args.evidence / "receipt.json").read_text())
    assert receipt["status"] == "success"
    assert receipt["freeze_sha256"] == sha(HERE / "FROZEN.json")
    assert receipt["input_sha256"] == {key.removesuffix("_sha256"): value
                                       for key, value in frozen.items()}
    assert [row["name"] for row in receipt["commands"]] == ["analysis", "independent_verify"]
    assert all(row["exit_code"] == 0 and row["peak_child_rss_bytes"] <= 512 * 1024 * 1024
               for row in receipt["commands"])
    assert receipt["result_sha256"] == {path.name: sha(path) for path in args.evidence.iterdir()
                                        if path.is_file() and path.name != "receipt.json"}
    archived = json.loads((args.evidence / "verification.json").read_text())
    assert archived["status"] == "PASS"
    with tempfile.TemporaryDirectory() as directory:
        fresh = Path(directory) / "fresh.json"
        subprocess.run([sys.executable, str(HERE / "verify.py"), "--expected",
                        str(args.evidence / "outcome.json"), "--out", str(fresh)],
                       check=True, timeout=120)
        assert json.loads(fresh.read_text()) == archived
    print("Independent all-q group-index portfolio replay PASS.")


if __name__ == "__main__":
    main()

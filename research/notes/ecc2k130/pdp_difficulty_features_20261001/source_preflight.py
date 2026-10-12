#!/usr/bin/env python3
"""Check the pre-outcome source and protocol bytes without scoring Q."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    lock_bytes = (HERE / "SOURCE_LOCK.json").read_bytes()
    lock = json.loads(lock_bytes)
    assert lock["schema"] == "ecc2k130-pdp-difficulty-source-lock-v1"
    observed = {}
    for relative, expected in lock["file_sha256"].items():
        actual = sha((ROOT / relative).read_bytes())
        assert actual == expected, relative
        observed[relative] = actual
    receipt = {"schema": "ecc2k130-pdp-difficulty-source-preflight-v1",
               "status": "PASS", "source_lock_sha256": sha(lock_bytes),
               "files_checked": len(observed), "file_sha256": observed}
    args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "PASS", "files_checked": len(observed)}, sort_keys=True))


if __name__ == "__main__":
    main()

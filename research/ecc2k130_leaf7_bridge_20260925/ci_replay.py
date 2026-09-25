#!/usr/bin/env python3
"""Hash-only freeze check. Never imports or executes the bridge producer."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    spec = json.loads((HERE / "FROZEN.json").read_text())
    assert spec["schema"] == "ecc2k130-leaf7-bridge-frozen-v1"
    assert spec["status"] == "protocol_only_no_outcome"
    assert spec["parent_note_commit"] == "e5de58c8be2984139734c31acb3d8824532f1b14"
    assert spec["release_main_head"] is None
    assert spec["release_gate"] == "hold_parent_unmerged"
    assert spec["gate"]["no_pdp_or_dlp_claim"] is True
    for relative, expected in spec["input_sha256"].items():
        actual = sha(REPO / relative)
        assert actual == expected, f"{relative}: expected {expected}, got {actual}"
    for name, expected in spec["implementation_sha256"].items():
        actual = sha(HERE / name)
        assert actual == expected, f"{name}: expected {expected}, got {actual}"
    print("Frozen source and implementation hashes match; no bridge map or outcome was run.")


if __name__ == "__main__":
    main()

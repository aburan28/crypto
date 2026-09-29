#!/usr/bin/env python3
"""Hash-only freeze check. Never imports or executes the bridge producer."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    spec = json.loads((HERE / "FROZEN.json").read_text())
    assert spec["schema"] == "ecc2k130-leaf7-bridge-frozen-v1"
    assert spec["status"] == "protocol_only_no_outcome"
    assert spec["parent_note_commit"] == "e5de58c8be2984139734c31acb3d8824532f1b14"
    assert spec["release_main_head"] == "8d99c3287a1d01aba8b9827aea9bd15397eb97d3"
    assert spec["release_gate"] == "hold_host_preflight_and_review"
    assert spec["host_image"] == {
        "reference": "sagemath/sagemath:10.9",
        "manifest_sha256": "2401ffa8e9fc85c7ea17d3649bde5958b4dbf0858b3e504098c4102720151711",
        "platform": "linux/amd64",
    }
    assert subprocess.run(["git", "merge-base", "--is-ancestor",
                           spec["release_main_head"], "HEAD"], cwd=REPO).returncode == 0
    workflow = (REPO / ".github/workflows/ecc2k130-leaf7-bridge-freeze.yml").read_text()
    assert spec["host_image"]["manifest_sha256"] in workflow
    assert spec["gate"]["no_pdp_or_dlp_claim"] is True
    assert spec["gate"]["producer_status"] == "PRODUCER_PASS"
    assert spec["gate"]["fq_replay_status"] == "FQ_REPLAY_PASS"
    for relative, expected in (spec["input_sha256"] | spec["host_refusal_sha256"]).items():
        actual = sha(REPO / relative)
        assert actual == expected, f"{relative}: expected {expected}, got {actual}"
    for name, expected in spec["implementation_sha256"].items():
        actual = sha(HERE / name)
        assert actual == expected, f"{name}: expected {expected}, got {actual}"
    print("Post-merge source, implementation and hosted-image hashes match; bridge held without outcome.")


if __name__ == "__main__":
    main()

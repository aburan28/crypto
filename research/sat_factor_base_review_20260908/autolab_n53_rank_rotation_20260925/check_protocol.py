#!/usr/bin/env python3
"""Held hash/input preflight. It never runs an n53 IC producer."""
from __future__ import annotations

import ast
import gzip
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
ORBIT = HERE.parent / "autolab_orbit_extract_20260924"
SHARED = HERE.parent / "autolab_shared_log_n53_20260925"
PARENT = HERE.parent / "autolab_combined_l384_coldbase_n53_20260925"
BASE_GZ = ORBIT / "independent_replay_20260924_codex/base_header.jsonl.gz"
HEX40 = re.compile(r"^[0-9a-f]{40}$")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sources() -> dict[str, Path]:
    return {
        "rust_producer": REPO / "examples/koblitz_s5_sat_instance.rs",
        "cargo_toml": REPO / "Cargo.toml",
        "cargo_lock": PARENT / "Cargo.lock",
        "schedule_generator": HERE / "generate_schedule.py",
        "runner": HERE / "run_panel.py",
        "audit": HERE / "audit.py",
        "preflight": HERE / "check_protocol.py",
        "independent_group": ORBIT / "independent_replay_20260924_codex/replay.py",
        "independent_rank": ORBIT / "cold_batch_rank.py",
        "old_target_generator": SHARED / "generate_targets.py",
        "resource_measure": PARENT / "run_panel.py",
        "archive_sealer": PARENT / "archive.py",
        "archive_verifier": HERE / "verify_archive.py",
        "workflow": REPO / ".github/workflows/n53-target-cyclic-rank.yml",
    }


def inputs() -> dict[str, Path]:
    return {
        "target_scalars": HERE / "target_scalars.txt",
        "target_points": HERE / "target_points.jsonl",
        "old_public_q": PARENT / "points_L384.jsonl",
        "old_public_q_manifest": SHARED / "TARGET_MANIFEST.json",
        "certified_base_gzip": BASE_GZ,
    }


def preflight(*, require_release: bool = False) -> dict:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert frozen["schema"] == "n53_target_cyclic_rank_factorial_freeze_v1"
    assert frozen["status"] in ("held_no_outcome", "released_for_one_outcome")
    if require_release:
        assert frozen["status"] == "released_for_one_outcome"
        assert frozen["release_main_head"] is not None
        assert platform.system() == "Linux" and platform.machine() == "x86_64"
        assert sys.version_info[:2] == (3, 12)
        assert subprocess.check_output(["rustc", "--version"], text=True).startswith("rustc 1.93.1 (01f6ddf75 ")
    assert sha(HERE / "PROTOCOL.md") == frozen["protocol_sha256"]
    for name, path in sources().items():
        assert sha(path) == frozen["source_sha256"][name], f"source hash drift: {name}"
        if path.suffix == ".py":
            ast.parse(path.read_text(), filename=str(path))
    for name, path in inputs().items():
        assert sha(path) == frozen["input_sha256"][name], f"input hash drift: {name}"
    import generate_schedule
    scalar_bytes, point_bytes, report = generate_schedule.generate()
    assert scalar_bytes == inputs()["target_scalars"].read_bytes()
    assert point_bytes == inputs()["target_points"].read_bytes()
    assert report == frozen["schedule_proof"]
    header_bytes = gzip.decompress(BASE_GZ.read_bytes())
    assert header_bytes.count(b"\n") == 1
    header = json.loads(header_bytes)
    assert (header["base_hash"], header["orbit_columns"], header["factor_base_points"]) == (
        "d859319015ea405fd18aee41b51396ce4edcab64ef66265d8edcdeb5e040eb71",
        220, 23320,
    )
    assert HEX40.fullmatch(frozen["parent_pr_815_head"])
    subprocess.run(["git", "merge-base", "--is-ancestor", frozen["parent_pr_815_head"], "HEAD"],
                   cwd=REPO, check=True)
    base_main = frozen["frozen_main_head"]
    assert HEX40.fullmatch(base_main)
    subprocess.run(["git", "merge-base", "--is-ancestor", base_main, "HEAD"],
                   cwd=REPO, check=True)
    release_main = frozen["release_main_head"]
    if release_main is not None:
        assert HEX40.fullmatch(release_main)
        subprocess.run(["git", "merge-base", "--is-ancestor", base_main, release_main],
                       cwd=REPO, check=True)
    if require_release:
        expected = os.environ.get("KIC_ROTATION_EXPECTED_HEAD", "")
        actual = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
        assert HEX40.fullmatch(expected) and actual == expected, "checkout differs from reviewed head"
        assert subprocess.check_output(["git", "status", "--porcelain"], cwd=REPO, text=True).strip() == ""
        subprocess.run(["git", "fetch", "--no-tags", "origin", "main"], cwd=REPO, check=True,
                       stdout=subprocess.DEVNULL)
        current_main = subprocess.check_output(["git", "rev-parse", "origin/main"],
                                               cwd=REPO, text=True).strip()
        assert HEX40.fullmatch(current_main)
        # An unrelated main fast-forward does not change the exact reviewed
        # PR checkout or any frozen source/input bytes; reject rewrites.
        subprocess.run(["git", "merge-base", "--is-ancestor",
                        frozen["release_main_head"], current_main], cwd=REPO, check=True)
    return frozen


if __name__ == "__main__":
    value = preflight()
    print(json.dumps({"status": "HASH_INPUT_PREFLIGHT_PASS", "release_main_head": value["release_main_head"],
                      "targets": value["schedule_proof"]["target_count"], "arms": 4}, sort_keys=True))

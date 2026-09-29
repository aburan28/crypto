#!/usr/bin/env python3
"""Hash and lineage gate for the held native-cyclic L384 successor; no producer."""
from __future__ import annotations

import ast
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
ROOT = HERE.parent
ORBIT = ROOT / "autolab_orbit_extract_20260924"
ROTATION = ROOT / "autolab_n53_rank_rotation_20260925"
SHARED = ROOT / "autolab_shared_log_n53_20260925"
RHO_STUDY = ROOT / "autolab_matched_point_rho_n53_20260925"
OLD = ROOT / "autolab_combined_l384_certified_coldbase_n53_20260925"
SHA40 = re.compile(r"^[0-9a-f]{40}$")
BASE_ARRAYS_SHA256 = "2de8ec46916999dfe5e1ac95e68ea15a3f001eec27e013fa96d9525719b4f684"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sources() -> dict[str, Path]:
    return {
        "ic": REPO / "examples/koblitz_s5_sat_instance.rs",
        "rho": REPO / "examples/koblitz_rho_batch_ks.rs",
        "cargo_toml": REPO / "Cargo.toml",
        "cargo_lock": HERE / "Cargo.lock",
        "independent_group": ORBIT / "independent_replay_20260924_codex/replay.py",
        "independent_rank": ORBIT / "cold_batch_rank.py",
        "independent_rho_math": RHO_STUDY / "analyze.py",
        "independent_policy": ROTATION / "audit.py",
        "input_math": SHARED / "generate_targets.py",
        "generator": HERE / "generate_inputs.py",
        "materializer": HERE / "cold_base.py",
        "operational": HERE / "operational.py",
        "runner": HERE / "run_panel.py",
        "audit": HERE / "audit.py",
        "archive_sealer": HERE / "archive.py",
        "archive_verifier": HERE / "verify_archive.py",
        "dispatch_gate": HERE / "dispatch_gate.py",
        "build": HERE / "build.py",
        "preflight": HERE / "check_protocol.py",
        "tests": HERE / "test_gates.py",
        "workflow": REPO / ".github/workflows/n53-native-cyclic-l384.yml",
    }


def inputs() -> dict[str, Path]:
    return {
        "A_scalars": HERE / "training_A_scalars.txt",
        "B_scalars": HERE / "training_B_scalars.txt",
        "B_points": HERE / "training_B_points.jsonl",
        "public_points": HERE / "points_L384.jsonl",
        "validator_scalars": HERE / "validator_scalars_L384.txt",
        "parent_A_scalars": ROTATION / "target_scalars.txt",
        "parent_A_points": ROTATION / "target_points.jsonl",
        "parent_archive": ROTATION / "evidence/evidence.tar.gz",
        "parent_archive_manifest": ROTATION / "evidence/archive_manifest.json",
        "old_training_manifest": SHARED / "TARGET_MANIFEST.json",
        "old_public_points": OLD / "points_L384.jsonl",
    }


def parent_native_arrays_hash() -> tuple[str, str]:
    archive = inputs()["parent_archive"]
    with tarfile.open(archive, "r:gz") as tar:
        member = tar.extractfile("panel/native_cyclic/producer.stdout.jsonl")
        assert member is not None
        raw_bytes = member.read()
    raw = json.loads(raw_bytes)
    header = raw["compact_orbit_base_header"]
    from cold_base import POINT_KEYS, canonical
    arrays = {key: header[key] for key in POINT_KEYS}
    return hashlib.sha256(canonical(arrays)).hexdigest(), hashlib.sha256(raw_bytes).hexdigest()


def preflight(*, require_release: bool = False) -> dict:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert frozen["schema"] == "n53_native_cyclic_l384_freeze_v1"
    assert frozen["status"] in ("held_no_outcome", "released_for_one_outcome")
    if require_release and (frozen["status"] != "released_for_one_outcome"
                            or frozen["release_main_head"] is None):
        raise AssertionError("HELD: exact-head review and a separate release commit are required")
    assert sha(HERE / "PROTOCOL.md") == frozen["protocol_sha256"]
    assert set(frozen["source_sha256"]) == set(sources())
    for key, path in sources().items():
        assert sha(path) == frozen["source_sha256"][key], f"source drift: {key}"
        if path.suffix == ".py":
            ast.parse(path.read_text(), filename=str(path))
    assert set(frozen["input_sha256"]) == set(inputs())
    for key, path in inputs().items():
        assert sha(path) == frozen["input_sha256"][key], f"input drift: {key}"
    assert inputs()["A_scalars"].read_bytes() == inputs()["parent_A_scalars"].read_bytes()
    archive_manifest = json.loads(inputs()["parent_archive_manifest"].read_text())
    assert archive_manifest["archive_sha256"] == sha(inputs()["parent_archive"])
    assert archive_manifest["archive_sha256"] == frozen["parent_archive_sha256"]
    native_hash, native_raw_hash = parent_native_arrays_hash()
    assert native_hash == frozen["native_base_arrays_sha256"] == BASE_ARRAYS_SHA256
    assert native_raw_hash == frozen["parent_native_raw_sha256"]
    import generate_inputs
    files, proof = generate_inputs.generate()
    assert proof == frozen["input_proof"]
    for name, data in files.items():
        assert (HERE / name).read_bytes() == data
    base_head = frozen["frozen_main_head"]
    assert SHA40.fullmatch(base_head)
    subprocess.run(["git", "merge-base", "--is-ancestor", base_head, "HEAD"], cwd=REPO, check=True)
    release_head = frozen["release_main_head"]
    if release_head is not None:
        assert SHA40.fullmatch(release_head)
        subprocess.run(["git", "merge-base", "--is-ancestor", base_head, release_head], cwd=REPO, check=True)
    if require_release:
        assert platform.system() == "Linux" and platform.machine() == "x86_64"
        assert sys.version_info[:2] == (3, 12)
        assert subprocess.check_output(["rustc", "--version"], text=True).startswith("rustc 1.93.1 ")
        expected = os.environ.get("KIC_NATIVE_L384_EXPECTED_HEAD", "")
        actual = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
        assert SHA40.fullmatch(expected) and actual == expected, "checkout differs from reviewed release head"
        assert not subprocess.check_output(["git", "status", "--porcelain"], cwd=REPO, text=True).strip()
        subprocess.run(["git", "fetch", "--no-tags", "origin", "main"], cwd=REPO, check=True,
                       stdout=subprocess.DEVNULL)
        current_main = subprocess.check_output(["git", "rev-parse", "origin/main"], cwd=REPO, text=True).strip()
        subprocess.run(["git", "merge-base", "--is-ancestor", release_head, current_main],
                       cwd=REPO, check=True)
    return frozen


if __name__ == "__main__":
    value = preflight()
    print(json.dumps({"verdict": "HELD_HASH_INPUT_PREFLIGHT_PASS",
                      "status": value["status"],
                      "release_main_head": value["release_main_head"],
                      "A": 512, "B": 512, "Q": 384}, sort_keys=True))

#!/usr/bin/env python3
"""Hash-only L384 exact-certified-base preflight; never launches a measured child."""
from __future__ import annotations

import ast
import gzip
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
ORBIT = HERE.parent / "autolab_orbit_extract_20260924"
SHARED = HERE.parent / "autolab_shared_log_n53_20260925"
RHO = HERE.parent / "autolab_matched_point_rho_n53_20260925"
PRIOR = HERE.parent / "autolab_combined_l384_n53_20260925"
BASE_GZ = ORBIT / "independent_replay_20260924_codex/base_header.jsonl.gz"
HEX40 = re.compile(r"^[0-9a-f]{40}$")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sources() -> dict[str, Path]:
    return {
        "ic": REPO / "examples/koblitz_s5_sat_instance.rs",
        "rank_fixture_provenance": REPO / "examples/koblitz_rank_fixture.rs",
        "rho": REPO / "examples/koblitz_rho_batch_ks.rs",
        "cargo_toml": REPO / "Cargo.toml",
        "cargo_lock": HERE / "Cargo.lock",
        "workflow": REPO / ".github/workflows/n53-l384-certified-coldbase.yml",
        "runner": HERE / "run_panel.py",
        "cold_base": HERE / "cold_base.py",
        "materializer_gate_test": HERE / "test_materializer_gate.py",
        "operational": HERE / "operational.py",
        "audit": HERE / "audit.py",
        "archive_sealer": HERE / "archive.py",
        "archive_verifier": HERE / "verify_archive.py",
        "preflight": HERE / "check_protocol.py",
        "training_schedule_reference": ORBIT / "cold_batch_rank.py",
        "independent_group_reference": ORBIT / "independent_replay_20260924_codex/replay.py",
        "rho_audit_math": RHO / "analyze.py",
        "point_generation_math": SHARED / "generate_targets.py",
        "prior_runner": PRIOR / "run_panel.py",
    }


def inputs() -> dict[str, Path]:
    return {
        "target_points": HERE / "points_L384.jsonl",
        "validator_manifest": SHARED / "TARGET_MANIFEST.json",
        "source_points_b0": SHARED / "points_b0_L128.jsonl",
        "source_points_b1": SHARED / "points_b1_L128.jsonl",
        "source_points_b2": SHARED / "points_b2_L128.jsonl",
        "certified_base_gzip": BASE_GZ,
    }


def preflight(*, require_release: bool = False) -> dict:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    if frozen["schema"] != "n53_l384_certified_coldbase_freeze_v1":
        raise AssertionError("freeze schema drift")
    if frozen["status"] not in ("held_no_outcome", "released_for_one_outcome"):
        raise AssertionError("unknown release status")
    if require_release and (frozen["status"] != "released_for_one_outcome"
                            or frozen["release_main_head"] is None):
        raise AssertionError("HELD: exact-head CI and peer-reviewed release required")
    if sha(HERE / "PROTOCOL.md") != frozen["protocol_sha256"]:
        raise AssertionError("protocol hash drift")
    for key, path in sources().items():
        if sha(path) != frozen["source_sha256"][key]:
            raise AssertionError(f"source hash drift: {key}")
        if path.suffix == ".py":
            ast.parse(path.read_text(), filename=str(path))
    for key, path in inputs().items():
        if sha(path) != frozen["input_sha256"][key]:
            raise AssertionError(f"input hash drift: {key}")
    parts = [inputs()[f"source_points_b{i}"].read_bytes() for i in range(3)]
    target_bytes = inputs()["target_points"].read_bytes()
    if target_bytes != b"".join(parts):
        raise AssertionError("L384 is not the pinned three-block concatenation")
    points = [tuple(json.loads(line)) for line in target_bytes.splitlines()]
    if len(points) != 384 or len(set(points)) != 384:
        raise AssertionError("wrong target count or duplicates")
    base = json.loads(gzip.decompress(BASE_GZ.read_bytes()))
    if (base["base_hash"], base["orbit_columns"], base["factor_base_points"],
            base["field_x_values_scanned"]) != (
            "d859319015ea405fd18aee41b51396ce4edcab64ef66265d8edcdeb5e040eb71",
            220, 23320, 400):
        raise AssertionError("certified base fixture drift")
    if frozen["source_sha256"]["rho"] != "fedacb54e441979c8c32860e7b5639e43799234741d49677010164564536d2c8":
        raise AssertionError("rho source differs from #760")
    parent_head = frozen["parent_result_head"]
    if not HEX40.fullmatch(parent_head):
        raise AssertionError("invalid parent result SHA")
    subprocess.run(["git", "merge-base", "--is-ancestor", parent_head, "HEAD"],
                   cwd=REPO, check=True)
    base_head = frozen["frozen_main_head"]
    if not HEX40.fullmatch(base_head):
        raise AssertionError("invalid frozen main SHA")
    subprocess.run(["git", "merge-base", "--is-ancestor", base_head, "HEAD"],
                   cwd=REPO, check=True)
    release = frozen["release_main_head"]
    if release is not None:
        if not HEX40.fullmatch(release):
            raise AssertionError("invalid release main SHA")
        # Main may advance after this PR branches. The base must precede the
        # frozen current main, while the PR's own exact HEAD is checked below.
        subprocess.run(["git", "merge-base", "--is-ancestor", base_head, release],
                       cwd=REPO, check=True)
    if require_release:
        expected_head = os.environ.get("KIC_L384_EXPECTED_HEAD", "")
        actual_head = subprocess.check_output(["git", "rev-parse", "HEAD"],
                                              cwd=REPO, text=True).strip()
        if not HEX40.fullmatch(expected_head) or actual_head != expected_head:
            raise AssertionError("checkout differs from exact reviewed PR head")
        subprocess.run(["git", "fetch", "--no-tags", "origin", "main"],
                       cwd=REPO, check=True, stdout=subprocess.DEVNULL)
        current_main = subprocess.check_output(["git", "rev-parse", "origin/main"],
                                               cwd=REPO, text=True).strip()
        if current_main != release:
            raise AssertionError("origin/main moved after release freeze")
        subprocess.run(["git", "merge-base", "--is-ancestor", parent_head, current_main],
                       cwd=REPO, check=True)
    return frozen


if __name__ == "__main__":
    current = preflight()
    print(json.dumps({"status": "HASH_ONLY_PASS", "release_main_head": current["release_main_head"],
                      "targets": 384, "columns": 220}, sort_keys=True))

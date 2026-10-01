#!/usr/bin/env python3
"""Build the three binaries from the v2 source freeze plus the new generator."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(
    0, str(ROOT / "research/notes/ecc2k130/compact_frozen_source_replay_20260929")
)
from materialize import materialize  # noqa: E402

SOURCE_FREEZE = ROOT / "research/notes/ecc2k130/compact_s3_prefilter_20260930/FROZEN.json"
SOURCE_FREEZE_SHA256 = "3e9f67cc2cd6de5a8458badb3525983d561d9c3449c05118819fb424c5093b2b"
FROZEN_LOCK = ROOT / "research/notes/ecc2k130/compact_shared_log_20260925/Cargo.lock"
FROZEN_LOCK_SHA256 = "7d671f48c2da93f133d98802e80f858d1d9ea3b86996f7037f758990e1566627"
FAST_ARITH_SHA256 = "2a5c54bdedb6a0b1ffcf95dd7badd32c41d1482b53055e31fb9f28d514ac8c9e"
EXAMPLES = (
    "koblitz_base_window",
    "koblitz_orbit_dlp_s3_batch",
    "koblitz_rho_batch_ks_v3",
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def build(out: Path, profile: str, offline: bool) -> dict:
    assert profile in ("dev", "release")
    assert not out.exists(), "never overwrite a frozen build"
    assert not out.is_relative_to(ROOT), "build directory must be outside the checkout"
    assert sha(SOURCE_FREEZE) == SOURCE_FREEZE_SHA256
    assert sha(FROZEN_LOCK) == FROZEN_LOCK_SHA256
    cfg = json.loads((HERE / "CONFIG.json").read_text())
    source_freeze = json.loads(SOURCE_FREEZE.read_text())
    assert source_freeze["source_sha256"]["examples/koblitz_orbit_dlp_s3_batch.rs"] == (
        cfg["compact_source_sha256"]
    )
    assert source_freeze["source_sha256"]["examples/koblitz_rho_batch_ks_v3.rs"] == (
        cfg["rho_source_sha256"]
    )
    out.mkdir(parents=True)
    source = out / "source"
    snapshot = materialize(ROOT, [SOURCE_FREEZE], source)
    write_json(out / "materialization.json", snapshot)
    shutil.copy2(ROOT / "examples/koblitz_base_window.rs",
                 source / "examples/koblitz_base_window.rs")
    shutil.copy2(FROZEN_LOCK, source / "Cargo.lock")
    for name, expected in source_freeze["source_sha256"].items():
        assert sha(source / name) == expected, name
    assert sha(source / "src/cryptanalysis/koblitz_fast_arith.rs") == FAST_ARITH_SHA256
    assert sha(source / "Cargo.lock") == FROZEN_LOCK_SHA256
    assert sha(source / "examples/koblitz_base_window.rs") == sha(
        ROOT / "examples/koblitz_base_window.rs"
    )
    command = ["cargo", "build", "--locked"]
    if offline:
        command.append("--offline")
    if profile == "release":
        command.append("--release")
    for name in EXAMPLES:
        command.extend(("--example", name))
    environment = dict(os.environ)
    environment["CARGO_TARGET_DIR"] = str(out / "target")
    receipt = {
        "schema": "ecc2k130-base-window-frozen-build-v1",
        "status": "RUNNING",
        "profile": profile,
        "offline": offline,
        "command": command,
        "source_freeze_sha256": SOURCE_FREEZE_SHA256,
        "materialization_sha256": sha(out / "materialization.json"),
        "pinned_files": snapshot["pinned_files"],
        "generator_source_sha256": sha(source / "examples/koblitz_base_window.rs"),
        "compact_source_sha256": cfg["compact_source_sha256"],
        "rho_source_sha256": cfg["rho_source_sha256"],
        "fast_arith_sha256": FAST_ARITH_SHA256,
        "cargo_lock_sha256": FROZEN_LOCK_SHA256,
        "rustc_version_verbose": subprocess.check_output(
            ["rustc", "--version", "--verbose"], text=True
        ).strip(),
        "cargo_version": subprocess.check_output(["cargo", "--version"], text=True).strip(),
    }
    write_json(out / "BUILD_RECEIPT.json", receipt)
    stdout = out / "build.stdout.txt"
    stderr = out / "build.stderr.txt"
    try:
        with stdout.open("xb") as output, stderr.open("xb") as errors:
            code = subprocess.call(
                command, cwd=source, env=environment, stdout=output, stderr=errors
            )
        receipt["exit_code"] = code
        receipt["build_stdout_sha256"] = sha(stdout)
        receipt["build_stderr_sha256"] = sha(stderr)
        if code != 0:
            receipt["status"] = "FAIL"
        else:
            directory = out / "target" / ("release" if profile == "release" else "debug") / "examples"
            receipt["binaries"] = {
                name: {
                    "path": str(directory / name),
                    "sha256": sha(directory / name),
                    "bytes": (directory / name).stat().st_size,
                }
                for name in EXAMPLES
            }
            receipt["status"] = "PASS"
    except BaseException as error:
        receipt["status"] = "FAIL"
        receipt["error_type"] = type(error).__name__
        receipt["error"] = str(error)
        receipt["traceback"] = traceback.format_exc()
        write_json(out / "BUILD_RECEIPT.json", receipt)
        raise
    write_json(out / "BUILD_RECEIPT.json", receipt)
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--profile", choices=("dev", "release"), default="release")
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    result = build(args.out.resolve(), args.profile, args.offline)
    print(json.dumps({
        "status": result["status"],
        "profile": result["profile"],
        "binaries": result.get("binaries"),
    }, sort_keys=True))
    if result["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()

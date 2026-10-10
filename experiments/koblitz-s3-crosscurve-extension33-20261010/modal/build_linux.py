"""Build the frozen cross-curve S3 executables for Linux from one Git snapshot."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tarfile
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
COMMIT = "b562cdaed14254d4f6c5f2a4640b302592e8882b"
TARGET = "x86_64-unknown-linux-gnu"
EXAMPLES = ("s3_pair_candidate_frozen", "s3_partial_prepared_targetspan")
ARCHIVE_PATHS = (
    "Cargo.toml", "Cargo.lock", "src", "docs/ic", "docs/ecbench", "docs/curves",
    "examples/s3_partial_prepared_targetspan.rs",
    "experiments/koblitz-s3-pair-query-20261007-v4/baseline.rs",
    "experiments/koblitz-s3-pair-query-20261007-v4/candidate.rs",
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def command(*argv: str) -> str:
    return subprocess.check_output(argv, cwd=ROOT, text=True).strip()


def build(build_dir: Path, receipt_path: Path, log_path: Path) -> None:
    if build_dir.exists() or receipt_path.exists() or log_path.exists():
        raise FileExistsError("build directory or evidence already exists")
    plan = json.loads(PARENT.read_bytes())
    assert plan["kind"] == "s3_rank243_matched_curve_fresh_panel_v1"
    assert command("git", "rev-parse", COMMIT) == COMMIT
    build_dir.mkdir(parents=True)
    archive = build_dir / "source.tar"
    subprocess.run(["git", "archive", "--format=tar", f"--output={archive}",
                    COMMIT, *ARCHIVE_PATHS], cwd=ROOT, check=True)
    source = build_dir / "source"
    source.mkdir()
    with tarfile.open(archive, "r:") as tape:
        for item in tape.getmembers():
            destination = (source / item.name).resolve()
            if not destination.is_relative_to(source) or item.issym() or item.islnk():
                raise ValueError(f"unexpected source archive member: {item.name}")
        tape.extractall(source, filter="data")
    assert sha(source / "experiments/koblitz-s3-pair-query-20261007-v4/candidate.rs") == plan["reference_source_sha256"]
    assert sha(source / "examples/s3_partial_prepared_targetspan.rs") == plan["partial_source_sha256"]
    zig = shutil.which("zig")
    if zig is None:
        raise FileNotFoundError("zig cross-linker is required")
    sysroot = Path(command("rustup", "run", "stable", "rustc", "--print", "sysroot"))
    rustc, cargo = sysroot / "bin/rustc", sysroot / "bin/cargo"
    target_lib = sysroot / "lib/rustlib" / TARGET / "lib"
    if not (rustc.is_file() and cargo.is_file() and any(target_lib.glob("libcore-*.rlib"))):
        raise FileNotFoundError(f"installed Rust toolchain lacks {TARGET}")
    linker = build_dir / "zig-cc-x86_64-gnu"
    archiver = build_dir / "zig-ar"
    linker.write_text(
        "#!/usr/bin/env python3\nimport os, sys\n"
        f"zig = {zig!r}\n"
        "args = [arg for arg in sys.argv[1:] if arg != '--target=x86_64-unknown-linux-gnu']\n"
        "os.execv(zig, [zig, 'cc', '-target', 'x86_64-linux-gnu.2.28', *args])\n"
    )
    archiver.write_text(f'#!/bin/sh\nexec "{zig}" ar "$@"\n')
    linker.chmod(0o755)
    archiver.chmod(0o755)
    env = dict(os.environ)
    env.update({
        "CARGO_TARGET_DIR": str(build_dir / "target"),
        "RUSTC": str(rustc),
        "CARGO_BUILD_JOBS": "8",
        "ZIG_GLOBAL_CACHE_DIR": str(build_dir / "zig-global-cache"),
        "ZIG_LOCAL_CACHE_DIR": str(build_dir / "zig-local-cache"),
        "CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER": str(linker),
        "CC_x86_64_unknown_linux_gnu": str(linker),
        "AR_x86_64_unknown_linux_gnu": str(archiver),
    })
    argv = [str(cargo), "build", "--offline", "--locked", "--release",
            "--target", TARGET]
    for example in EXAMPLES:
        argv += ["--example", example]
    started = time.time_ns()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("wb") as log:
        result = subprocess.run(argv, cwd=source, env=env, stdout=log,
                                stderr=subprocess.STDOUT, check=False)
    finished = time.time_ns()
    binaries = {}
    if result.returncode == 0:
        for name in EXAMPLES:
            binary = build_dir / "target" / TARGET / "release/examples" / name
            binaries[name] = {"path": str(binary), "bytes": binary.stat().st_size,
                              "sha256": sha(binary), "file": command("file", str(binary))}
            assert "ELF 64-bit LSB" in binaries[name]["file"]
            assert "x86-64" in binaries[name]["file"]
    receipt = {
        "kind": "s3_crosscurve_modal_linux_build_v1",
        "source_commit": COMMIT,
        "source_tree_git_oid": command("git", "rev-parse", f"{COMMIT}:src"),
        "source_archive_sha256": sha(archive),
        "source_archive_paths": ARCHIVE_PATHS,
        "source_sha256": {"reference": plan["reference_source_sha256"],
                          "prepared": plan["partial_source_sha256"]},
        "cargo_lock_sha256": sha(source / "Cargo.lock"),
        "argv": argv, "exit_code": result.returncode,
        "started_unix_ns": started, "finished_unix_ns": finished,
        "log_sha256": sha(log_path), "host": platform.platform(),
        "rustc": command(str(rustc), "--version", "--verbose"),
        "rust_sysroot": str(sysroot), "zig": command(zig, "version"),
        "target": TARGET, "binaries": binaries,
    }
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"exit_code": result.returncode, "binaries": binaries}, sort_keys=True))
    if result.returncode:
        raise SystemExit(result.returncode)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--log", type=Path, required=True)
    args = parser.parse_args()
    build(args.build_dir.resolve(), args.receipt.resolve(), args.log.resolve())

"""Cross-build the C1200 and C1800 combined S3 examples for a Modal pilot."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TARGET = "x86_64-unknown-linux-gnu"
BUILD_DIR = Path("/private/tmp/s3-modal-combined-build")
EXAMPLES = (
    "s3_partial_prepared_fullrank_c2800mixed_targetspan",
    "s3_partial_prepared_fullrank_prefetch_projective_c1800prefix_targetspan",
)
SOURCE_HASHES = {
    EXAMPLES[0]: "a398b36071cd59f7d6775c8f4e2b03a6bf867319ba1636378d186f7fb738c679",
    EXAMPLES[1]: "79f77fef1857c9aed10b8d8913ea796e501f9b0afd615345746c453f8d4907b2",
}
LIBRARY_TREE_OID = "0a574800502e67882d8f8f88a846c1e80ee795fd"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def command(*argv: str) -> str:
    return subprocess.check_output(argv, cwd=ROOT, text=True).strip()


def build(receipt_path: Path, log_path: Path) -> None:
    if receipt_path.exists() or log_path.exists():
        raise FileExistsError("refusing to overwrite cross-build evidence")
    assert command("git", "rev-parse", "HEAD:src") == LIBRARY_TREE_OID
    for name, digest in SOURCE_HASHES.items():
        assert sha(ROOT / "examples" / f"{name}.rs") == digest
    zig = shutil.which("zig")
    if zig is None:
        raise FileNotFoundError("zig cross-linker is required")
    # The Homebrew rustc on PATH has no Linux standard library. Select the
    # rustup toolchain whose target library is actually installed.
    sysroot = Path(command("rustup", "run", "stable", "rustc", "--print", "sysroot"))
    rustc = sysroot / "bin" / "rustc"
    cargo = sysroot / "bin" / "cargo"
    target_lib = sysroot / "lib" / "rustlib" / TARGET / "lib"
    if not (rustc.is_file() and cargo.is_file() and
            any(target_lib.glob("libcore-*.rlib"))):
        raise FileNotFoundError(f"installed Rust toolchain lacks {TARGET}: {sysroot}")
    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    linker = BUILD_DIR / "zig-cc-x86_64-gnu"
    archiver = BUILD_DIR / "zig-ar"
    linker.write_text(
        "#!/usr/bin/env python3\n"
        "import os, sys\n"
        f"zig = {zig!r}\n"
        "args = [arg for arg in sys.argv[1:] "
        "if arg != '--target=x86_64-unknown-linux-gnu']\n"
        "os.execv(zig, [zig, 'cc', '-target', 'x86_64-linux-gnu.2.28', *args])\n"
    )
    archiver.write_text(f'#!/bin/sh\nexec "{zig}" ar "$@"\n')
    linker.chmod(0o755)
    archiver.chmod(0o755)
    env = dict(os.environ)
    env["CARGO_TARGET_DIR"] = str(BUILD_DIR / "target")
    env["RUSTC"] = str(rustc)
    env["CARGO_BUILD_JOBS"] = "8"
    env["ZIG_GLOBAL_CACHE_DIR"] = str(BUILD_DIR / "zig-global-cache")
    env["ZIG_LOCAL_CACHE_DIR"] = str(BUILD_DIR / "zig-local-cache")
    env["CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER"] = str(linker)
    env["CC_x86_64_unknown_linux_gnu"] = str(linker)
    env["AR_x86_64_unknown_linux_gnu"] = str(archiver)
    argv = [str(cargo), "build", "--offline", "--locked", "--release", "--target", TARGET]
    for name in EXAMPLES:
        argv += ["--example", name]
    started = time.time_ns()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("wb") as log:
        result = subprocess.run(argv, cwd=ROOT, env=env, stdout=log,
                                stderr=subprocess.STDOUT, check=False)
    finished = time.time_ns()
    binaries = {}
    if result.returncode == 0:
        for name in EXAMPLES:
            path = BUILD_DIR / "target" / TARGET / "release" / "examples" / name
            binaries[name] = {
                "path": str(path), "bytes": path.stat().st_size,
                "sha256": sha(path), "file": command("file", str(path)),
            }
    receipt = {
        "kind": "s3_modal_c1800_combined_cross_build_v1",
        "argv": argv, "exit_code": result.returncode,
        "started_unix_ns": started, "finished_unix_ns": finished,
        "log_sha256": sha(log_path),
        "host": platform.platform(),
        "rustc": command(str(rustc), "--version", "--verbose"),
        "rust_sysroot": str(sysroot),
        "cargo_build_jobs": 8,
        "zig": command(zig, "version"),
        "target": TARGET, "library_source_tree_git_oid": LIBRARY_TREE_OID,
        "source_sha256": SOURCE_HASHES,
        "cargo_lock_sha256": sha(ROOT / "Cargo.lock"),
        "binaries": binaries,
    }
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"exit_code": result.returncode, "binaries": binaries}, sort_keys=True))
    if result.returncode:
        raise SystemExit(result.returncode)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--log", type=Path, required=True)
    args = parser.parse_args()
    build(args.receipt, args.log)

#!/usr/bin/env python3
"""Build the two measured executables from one clean, pinned checkout."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import shutil
import time

from run import HERE, REPO, EXECUTABLES, SPEC, sha, write_json

LOCK_SHA = "7d671f48c2da93f133d98802e80f858d1d9ea3b86996f7037f758990e1566627"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    assert not out.is_relative_to(REPO.resolve())
    out.mkdir(parents=True, exist_ok=False)
    spec = json.loads(SPEC.read_bytes())
    snapshot = HERE / "Cargo.lock"
    assert sha(snapshot.read_bytes()) == LOCK_SHA
    root_lock = REPO / "Cargo.lock"
    if root_lock.exists():
        assert sha(root_lock.read_bytes()) == LOCK_SHA, "mutable root lock differs from frozen snapshot"
    else:
        shutil.copyfile(snapshot, root_lock)
    assert sha(root_lock.read_bytes()) == LOCK_SHA
    assert not subprocess.check_output(["git", "status", "--porcelain"],
                                       cwd=REPO, text=True).strip()
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"],
                                     cwd=REPO, text=True).strip()
    assert subprocess.run(["git", "merge-base", "--is-ancestor",
                           spec["base_commit"], commit], cwd=REPO).returncode == 0
    for path, expected in spec["source_sha256"].items():
        assert sha((REPO / path).read_bytes()) == expected
    command = ["cargo", "build", "--release", "--offline", "--locked",
               "--example", "koblitz_s5_sat_instance",
               "--example", "koblitz_rho_batch_ks"]
    started = time.monotonic_ns()
    with (out / "build.stdout.txt").open("wb") as stdout, (out / "build.stderr.txt").open("wb") as stderr:
        process = subprocess.run(command, cwd=REPO, stdout=stdout, stderr=stderr,
                                 env={**os.environ, "CARGO_NET_OFFLINE": "true"})
    receipt = {"schema_version": "1.0", "checkout_head": commit,
               "input_spec_sha256": sha(SPEC.read_bytes()),
               "cargo_lock_sha256": LOCK_SHA,
               "cargo_toml_sha256": sha((REPO / "Cargo.toml").read_bytes()),
               "cargo_net_offline": True,
               "source_sha256": spec["source_sha256"],
               "command": command, "returncode": process.returncode,
               "wall_ms": (time.monotonic_ns()-started)/1e6,
               "host": platform.node(), "platform": platform.platform(),
               "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
               "cargo": subprocess.check_output(["cargo", "--version"], text=True).strip(),
               "stdout_sha256": sha((out / "build.stdout.txt").read_bytes()),
               "stderr_sha256": sha((out / "build.stderr.txt").read_bytes())}
    if process.returncode == 0:
        receipt["executable_sha256"] = {
            key: sha((REPO / path).read_bytes()) for key, path in EXECUTABLES.items()
            if key != "train"}
    write_json(out / "build_receipt.json", receipt)
    assert process.returncode == 0, "pinned build failed; raw receipt retained"
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()

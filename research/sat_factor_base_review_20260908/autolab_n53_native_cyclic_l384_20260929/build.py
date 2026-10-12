#!/usr/bin/env python3
"""Preserve the excluded compilation phase, including its first failure."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time
import traceback

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
COMMAND = ["cargo", "build", "--locked", "--release", "--example",
           "koblitz_s5_sat_instance", "--example", "koblitz_rho_batch_ks"]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(panel: Path):
    assert panel.is_dir()
    assert (panel / "predispatch.json").is_file()
    receipt_path = panel / "build_receipt.json"
    assert not receipt_path.exists()
    stdout = panel / "build.stdout.txt"
    stderr = panel / "build.stderr.txt"
    receipt = {"schema": "n53_native_cyclic_l384_build_v1",
               "command": COMMAND,
               "status": "STARTED",
               "cargo_lock_sha256": sha(HERE / "Cargo.lock"),
               "cargo_toml_sha256": sha(REPO / "Cargo.toml"),
               "ic_source_sha256": sha(REPO / "examples/koblitz_s5_sat_instance.rs"),
               "rho_source_sha256": sha(REPO / "examples/koblitz_rho_batch_ks.rs"),
               "toolchain": {}}
    def save():
        receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    save()
    try:
        assert sha(REPO / "Cargo.lock") == receipt["cargo_lock_sha256"]
        for command in ("rustc", "cargo"):
            text = subprocess.check_output([command, "--version"], text=True).strip()
            receipt["toolchain"][command] = text
            assert text.startswith(command + " 1.93.1 "), text
        started = time.monotonic()
        with stdout.open("wb") as out_stream, stderr.open("wb") as err_stream:
            process = subprocess.run(COMMAND, cwd=REPO, stdout=out_stream, stderr=err_stream,
                                     timeout=900)
        receipt["wall_s"] = time.monotonic() - started
        receipt["returncode"] = process.returncode
        receipt["stdout_sha256"] = sha(stdout)
        receipt["stderr_sha256"] = sha(stderr)
        if process.returncode == 0:
            receipt["binary_sha256"] = {
                "ic": sha(REPO / "target/release/examples/koblitz_s5_sat_instance"),
                "rho": sha(REPO / "target/release/examples/koblitz_rho_batch_ks"),
            }
            receipt["status"] = "SUCCESS"
        else:
            receipt["status"] = "FAILED"
        save()
        assert process.returncode == 0, "frozen producer build failed"
    except Exception as exc:
        if receipt["status"] == "STARTED":
            receipt["status"] = "FAILED_OR_CENSORED"
        receipt["error"] = f"{type(exc).__name__}: {exc}"
        receipt["traceback"] = traceback.format_exc()
        save()
        raise
    print(json.dumps({"build": receipt["status"], "wall_s": receipt["wall_s"],
                      "binary_sha256": receipt["binary_sha256"]}, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, required=True)
    args = parser.parse_args()
    run(args.panel)


if __name__ == "__main__":
    main()

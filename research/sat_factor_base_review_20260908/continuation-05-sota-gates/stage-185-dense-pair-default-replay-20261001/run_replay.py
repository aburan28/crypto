#!/usr/bin/env python3
"""Run selected default F4 and same-binary direct MITM once."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


STAGE = Path(__file__).resolve().parent
REPO = next(parent for parent in STAGE.parents if (parent / "Cargo.toml").is_file())
SOURCE = Path("/Volumes/SSD990/crypto-kic-stage185-source")
MANIFEST = STAGE.parent / "stage-175-current-f4-single-target-20261001" / "input" / "manifest.json"
METER = REPO / "scripts" / "process_meter.py"
EXPECTED_BINARY = "3aefd5cd0daf73602bf8a38e7fdcb363e630ffd57bfebe47a37123b1d03429b1"
LOCK_SHA256 = "4f17b356fa7bac392b6d801d1c74fb9e36b6517f9465c8ebc19bb9a2792a84c5"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run_one(binary: Path, backend: str, output: Path) -> dict:
    root = output / backend
    root.mkdir()
    controls = [
        "VECLIB_MAXIMUM_THREADS=1",
        "OPENBLAS_NUM_THREADS=1",
        "OMP_NUM_THREADS=1",
        "MKL_NUM_THREADS=1",
        "BLIS_NUM_THREADS=1",
        "NUMEXPR_NUM_THREADS=1",
    ]
    if backend == "native-f4":
        controls[:0] = ["RAYON_NUM_THREADS=12", "PQ_F4_X1_BATCH=512"]
    command = [
        sys.executable,
        str(METER),
        "--cwd",
        str(SOURCE),
        "--timeout",
        "360" if backend == "native-f4" else "60",
        "--stdout",
        str(root / "stdout.json"),
        "--stderr",
        str(root / "stderr.txt"),
        "--metrics",
        str(root / "metrics.json"),
        "--exclusive-create",
        "--",
        "/usr/bin/env",
        *controls,
        str(binary),
        backend,
        str(MANIFEST),
    ]
    if backend == "native-f4":
        command.append("300")
    subprocess.run(command, check=False)
    metrics = json.loads((root / "metrics.json").read_text())
    report = json.loads((root / "stdout.json").read_text())
    extra = report.get("cost", {}).get("extra", {})
    summary = {
        "backend": backend,
        "returncode": metrics["returncode"],
        "status": report.get("status"),
        "wall_seconds": metrics["metrics"]["wall_seconds"],
        "core_seconds": metrics["metrics"]["total_core_seconds"],
        "single_core_seconds": metrics["metrics"]["single_core_seconds"],
        "peak_rss_bytes": metrics["metrics"]["peak_rss_bytes"],
        "equations": report.get("solver_equations_blake3"),
        "logical_xors": report.get("cost", {}).get("ops"),
        "performed_xors": extra.get("word_xors_performed"),
        "pair_dense_select_calls": extra.get("pair_dense_select_calls"),
        "pair_quadratic_select_calls": extra.get("pair_quadratic_select_calls"),
        "pair_dense_scratch_bytes_max": extra.get("pair_dense_scratch_bytes_max"),
        "full_m4ri_matrices": extra.get("full_m4ri_matrices"),
    }
    print(json.dumps(summary, sort_keys=True), flush=True)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=STAGE / "development" / "replay")
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    if sha256(binary) != EXPECTED_BINARY:
        raise SystemExit("binary SHA-256 mismatch")
    lock = SOURCE / "Cargo.lock"
    if sha256(lock) != LOCK_SHA256:
        raise SystemExit("supplied Cargo.lock SHA-256 mismatch")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    provenance = {
        "schema": "koblitz_stage185_replay_provenance.v1",
        "selection_commit": "8014149a2b55cd2cca202237a302643e39a50f6e",
        "runner_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
        ).strip(),
        "source_checkout": str(SOURCE),
        "binary": str(binary),
        "binary_sha256": sha256(binary),
        "supplied_lock_sha256": sha256(lock),
        "pair_selector_environment": "unset",
        "full_m4ri_environment": "unset",
    }
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    run_one(binary, "native-f4", output)
    run_one(binary, "direct-mitm", output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

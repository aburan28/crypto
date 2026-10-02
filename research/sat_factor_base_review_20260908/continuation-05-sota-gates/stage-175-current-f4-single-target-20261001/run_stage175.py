#!/usr/bin/env python3
"""Run the frozen Stage 175 single-target F4/direct panel."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


STAGE = Path(__file__).resolve().parent
REPO = next(parent for parent in STAGE.parents if (parent / "Cargo.toml").is_file())
MANIFEST = STAGE / "input" / "manifest.json"
METER = REPO / "scripts" / "process_meter.py"
EXPECTED_MANIFEST_SHA256 = "188dfdcb7e9398b0ab03193265f7b3035776d4a8a52de7dd4a51475dfc1a572c"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_one(binary: Path, backend: str, repeat: int, root: Path) -> None:
    out = root / f"{backend}-r{repeat}"
    out.mkdir(parents=True, exist_ok=False)
    environment = [
        "VECLIB_MAXIMUM_THREADS=1",
        "OPENBLAS_NUM_THREADS=1",
        "OMP_NUM_THREADS=1",
        "MKL_NUM_THREADS=1",
        "BLIS_NUM_THREADS=1",
        "NUMEXPR_NUM_THREADS=1",
    ]
    if backend == "native-f4":
        environment[:0] = ["RAYON_NUM_THREADS=12", "PQ_F4_X1_BATCH=512"]
    command = [
        sys.executable,
        str(METER),
        "--cwd",
        str(REPO),
        "--timeout",
        "360" if backend == "native-f4" else "60",
        "--stdout",
        str(out / "stdout.json"),
        "--stderr",
        str(out / "stderr.txt"),
        "--metrics",
        str(out / "metrics.json"),
        "--exclusive-create",
        "--",
        "/usr/bin/env",
        *environment,
        str(binary),
        backend,
        str(MANIFEST),
    ]
    if backend == "native-f4":
        command.append("300")
    subprocess.run(command, check=False)

    metrics = json.loads((out / "metrics.json").read_text())
    report = json.loads((out / "stdout.json").read_text())
    process = metrics["metrics"]
    extra = report.get("cost", {}).get("extra", {})
    print(
        json.dumps(
            {
                "backend": backend,
                "repeat": repeat,
                "returncode": metrics["returncode"],
                "status": report.get("status"),
                "wall_seconds": process["wall_seconds"],
                "core_seconds": process["total_core_seconds"],
                "peak_rss_bytes": process["peak_rss_bytes"],
                "ops": report.get("cost", {}).get("ops"),
                "word_xors_performed": extra.get("word_xors_performed"),
                "equation_fingerprint": report.get("solver_equations_blake3"),
            },
            sort_keys=True,
        ),
        flush=True,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=STAGE / "development" / "campaign")
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()

    if args.repeats <= 0:
        parser.error("--repeats must be positive")
    if sha256(MANIFEST) != EXPECTED_MANIFEST_SHA256:
        raise SystemExit("frozen manifest SHA-256 mismatch")
    binary = args.binary.resolve(strict=True)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    provenance = {
        "schema": "koblitz_stage175_runner_provenance.v1",
        "source_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
        ).strip(),
        "binary": str(binary),
        "binary_sha256": sha256(binary),
        "manifest": str(MANIFEST),
        "manifest_sha256": sha256(MANIFEST),
        "repeats": args.repeats,
        "host": subprocess.check_output(["uname", "-a"], text=True).strip(),
        "python": sys.version,
        "environment_controls": {
            "RAYON_NUM_THREADS": "12",
            "PQ_F4_X1_BATCH": "512",
            "BLAS_threads": "1",
        },
    }
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")

    for repeat in range(1, args.repeats + 1):
        run_one(binary, "native-f4", repeat, output)
        run_one(binary, "direct-mitm", repeat, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Run the Stage 182 single-core current/full screen."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


STAGE = Path(__file__).resolve().parent
REPO = next(parent for parent in STAGE.parents if (parent / "Cargo.toml").is_file())
MANIFEST = STAGE.parent / "stage-175-current-f4-single-target-20261001" / "input" / "manifest.json"
METER = REPO / "scripts" / "process_meter.py"
EXPECTED_BINARY = "fe17f002fedb2e4ef3037d100a4c541ffe0d2c7dcbfa6b2121a540804cc9050c"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=STAGE / "development" / "screen")
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    if sha256(binary) != EXPECTED_BINARY:
        raise SystemExit("binary SHA-256 mismatch")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output / "provenance.json").write_text(
        json.dumps(
            {
                "schema": "koblitz_stage182_screen_provenance.v1",
                "binary_source_commit": "850182efba8d8a9499e82aa6f11600941d3cb9d6",
                "runner_commit": subprocess.check_output(
                    ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
                ).strip(),
                "binary": str(binary),
                "binary_sha256": sha256(binary),
                "sequence": ["current", "full"],
                "single_core_controls": True,
            },
            indent=2,
        )
        + "\n"
    )
    controls = [
        "RAYON_NUM_THREADS=1",
        "PQ_F4_X1_BATCH=1",
        "VECLIB_MAXIMUM_THREADS=1",
        "OPENBLAS_NUM_THREADS=1",
        "OMP_NUM_THREADS=1",
        "MKL_NUM_THREADS=1",
        "BLIS_NUM_THREADS=1",
        "NUMEXPR_NUM_THREADS=1",
    ]
    for order, variant in enumerate(("current", "full"), 1):
        root = output / f"{order:02d}-{variant}"
        root.mkdir()
        mode = "0" if variant == "current" else "1"
        subprocess.run(
            [
                sys.executable,
                str(METER),
                "--cwd",
                str(REPO),
                "--timeout",
                "600",
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
                f"F4_F2_FULL_M4RI={mode}",
                str(binary),
                "native-f4",
                str(MANIFEST),
                "540",
            ],
            check=False,
        )
        metrics = json.loads((root / "metrics.json").read_text())
        report = json.loads((root / "stdout.json").read_text())
        extra = report.get("cost", {}).get("extra", {})
        print(
            json.dumps(
                {
                    "order": order,
                    "variant": variant,
                    "returncode": metrics["returncode"],
                    "status": report.get("status"),
                    "wall_seconds": metrics["metrics"]["wall_seconds"],
                    "core_seconds": metrics["metrics"]["total_core_seconds"],
                    "single_core_seconds": metrics["metrics"]["single_core_seconds"],
                    "peak_rss_bytes": metrics["metrics"]["peak_rss_bytes"],
                    "logical_xors": report.get("cost", {}).get("ops"),
                    "performed_xors": extra.get("word_xors_performed"),
                    "full_m4ri_matrices": extra.get("full_m4ri_matrices"),
                    "single_thread_requested": report.get("single_thread_requested"),
                    "equations": report.get("solver_equations_blake3"),
                },
                sort_keys=True,
            ),
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

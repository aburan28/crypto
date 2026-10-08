#!/usr/bin/env python3
"""Run three fresh interleaved current/full-M4RI confirmation pairs."""

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
SEQUENCE = ("current", "full", "full", "current", "current", "full")
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
    parser.add_argument("--output", type=Path, default=STAGE / "development" / "confirmation")
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    if sha256(binary) != EXPECTED_BINARY:
        raise SystemExit("binary SHA-256 mismatch")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output / "provenance.json").write_text(
        json.dumps(
            {
                "schema": "koblitz_stage181_confirmation_provenance.v1",
                "source_commit": "850182efba8d8a9499e82aa6f11600941d3cb9d6",
                "runner_commit": subprocess.check_output(
                    ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
                ).strip(),
                "binary": str(binary),
                "binary_sha256": sha256(binary),
                "sequence": SEQUENCE,
            },
            indent=2,
        )
        + "\n"
    )
    controls = [
        "RAYON_NUM_THREADS=12",
        "PQ_F4_X1_BATCH=512",
        "VECLIB_MAXIMUM_THREADS=1",
        "OPENBLAS_NUM_THREADS=1",
        "OMP_NUM_THREADS=1",
        "MKL_NUM_THREADS=1",
        "BLIS_NUM_THREADS=1",
        "NUMEXPR_NUM_THREADS=1",
    ]
    counts = {"current": 0, "full": 0}
    for order, variant in enumerate(SEQUENCE, 1):
        counts[variant] += 1
        repeat = counts[variant]
        root = output / f"{order:02d}-{variant}-r{repeat}"
        root.mkdir()
        mode = "0" if variant == "current" else "1"
        subprocess.run(
            [
                sys.executable,
                str(METER),
                "--cwd",
                str(REPO),
                "--timeout",
                "360",
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
                "300",
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
                    "repeat": repeat,
                    "returncode": metrics["returncode"],
                    "status": report.get("status"),
                    "wall_seconds": metrics["metrics"]["wall_seconds"],
                    "core_seconds": metrics["metrics"]["total_core_seconds"],
                    "peak_rss_bytes": metrics["metrics"]["peak_rss_bytes"],
                    "logical_xors": report.get("cost", {}).get("ops"),
                    "performed_xors": extra.get("word_xors_performed"),
                    "full_m4ri_matrices": extra.get("full_m4ri_matrices"),
                    "full_m4ri_blocks": extra.get("full_m4ri_blocks"),
                    "full_m4ri_table_word_xors": extra.get("full_m4ri_table_word_xors"),
                    "peak_table_bytes": extra.get("peak_table_bytes"),
                    "equations": report.get("solver_equations_blake3"),
                },
                sort_keys=True,
            ),
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

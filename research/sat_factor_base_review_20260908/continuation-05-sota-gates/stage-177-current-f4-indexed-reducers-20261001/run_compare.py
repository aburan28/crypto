#!/usr/bin/env python3
"""Run three interleaved exact-submask/linear-control pairs."""

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
SEQUENCE = ("linear", "indexed", "indexed", "linear", "linear", "indexed")
EXPECTED_BINARY = "07bc31f962de6329b3e1fe76ca5311ffbdb82c69df8650bab6eb99d318ee115b"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=STAGE / "development" / "paired")
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    if sha256(binary) != EXPECTED_BINARY:
        raise SystemExit("binary SHA-256 mismatch")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output / "provenance.json").write_text(
        json.dumps(
            {
                "schema": "koblitz_stage177_pair_provenance.v1",
                "binary": str(binary),
                "binary_sha256": sha256(binary),
                "binary_source_commit": "a337e2e901848ac440d70d88c838d37937467413",
                "runner_commit": subprocess.check_output(
                    ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
                ).strip(),
                "sequence": SEQUENCE,
                "manifest": str(MANIFEST),
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
    counts = {"linear": 0, "indexed": 0}
    for order, variant in enumerate(SEQUENCE, 1):
        counts[variant] += 1
        pair = counts[variant]
        root = output / f"{order:02d}-{variant}-r{pair}"
        root.mkdir()
        mode = "0" if variant == "linear" else "1"
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
                f"F4_F2_INDEXED_REDUCERS={mode}",
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
                    "repeat": pair,
                    "returncode": metrics["returncode"],
                    "status": report.get("status"),
                    "wall_seconds": metrics["metrics"]["wall_seconds"],
                    "core_seconds": metrics["metrics"]["total_core_seconds"],
                    "peak_rss_bytes": metrics["metrics"]["peak_rss_bytes"],
                    "divisor_tests": extra.get("divisor_tests"),
                    "submask_lookups": extra.get("divisor_submask_lookups"),
                    "linear_tests": extra.get("divisor_linear_tests"),
                    "ops": report.get("cost", {}).get("ops"),
                    "performed": extra.get("word_xors_performed"),
                    "equations": report.get("solver_equations_blake3"),
                },
                sort_keys=True,
            ),
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Run the frozen two-order four-/five-column screen."""

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
SEQUENCE = ("five", "four", "four", "five")
EXPECTED = {
    "five": "b156f3320564c1e9fed2b168aba7a0b4f3c875d677c1a9354751bba527ee5bd6",
    "four": "58870d396aecc33cd51ea51c4506301d14b8ea7d56044e87bad281e514af96af",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--five-binary", type=Path, required=True)
    parser.add_argument("--four-binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=STAGE / "development" / "screen")
    args = parser.parse_args()
    binaries = {
        "five": args.five_binary.resolve(strict=True),
        "four": args.four_binary.resolve(strict=True),
    }
    for variant, binary in binaries.items():
        if sha256(binary) != EXPECTED[variant]:
            raise SystemExit(f"{variant} binary SHA-256 mismatch")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output / "provenance.json").write_text(
        json.dumps(
            {
                "schema": "koblitz_stage180_screen_provenance.v1",
                "source_commit": "64e5272acc814b29b9c6c865d814c5bcdb4b439c",
                "runner_commit": subprocess.check_output(
                    ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
                ).strip(),
                "binary_sha256": EXPECTED,
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
    counts = {"five": 0, "four": 0}
    for order, variant in enumerate(SEQUENCE, 1):
        counts[variant] += 1
        repeat = counts[variant]
        root = output / f"{order:02d}-{variant}-r{repeat}"
        root.mkdir()
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
                str(binaries[variant]),
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
                    "ops": report.get("cost", {}).get("ops"),
                    "performed": extra.get("word_xors_performed"),
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

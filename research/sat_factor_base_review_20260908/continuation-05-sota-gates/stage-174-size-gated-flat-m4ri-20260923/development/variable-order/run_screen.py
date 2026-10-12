#!/usr/bin/env python3
import json
import os
from pathlib import Path
import subprocess
import sys

REPO = Path("/Volumes/SSD990/crypto-kic-stage174-native-f4")
ROOT = Path("/Volumes/SSD990/kic-stage174-f4-order-dev")
BINARY = Path("/Volumes/SSD990/kic-stage174-f4-target/release/examples/koblitz_pdp_backend")
MANIFEST = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-170-parallel-fixed-x1-construction-20260923/selected-run/tasks/000000-b-421e22a9c1c3b9d56396c8bbd0e46185bbebc32de0306f2bee6e9585703c2be4/instance/manifest.json"
METER = REPO / "scripts/process_meter.py"

orders = sys.argv[1:] or [
    "identity",
    "reverse",
    "basis-reverse",
    "interleave",
    "interleave-basis-reverse",
]

label = os.environ.get("RUN_LABEL", "screen")
for index, order in enumerate(orders):
    name = f"{label}-{index:02d}-{order.replace(',', '-')[:80]}"
    out = ROOT / name
    out.mkdir(parents=True, exist_ok=False)
    command = [
        sys.executable,
        str(METER),
        "--cwd",
        str(REPO),
        "--timeout",
        "360",
        "--stdout",
        str(out / "stdout.json"),
        "--stderr",
        str(out / "stderr.txt"),
        "--metrics",
        str(out / "metrics.json"),
        "--exclusive-create",
        "--",
        "/usr/bin/env",
        "RAYON_NUM_THREADS=12",
        "PQ_F4_X1_BATCH=512",
        f"PQ_F4_VARIABLE_ORDER={order}",
        "VECLIB_MAXIMUM_THREADS=1",
        "OPENBLAS_NUM_THREADS=1",
        "OMP_NUM_THREADS=1",
        "MKL_NUM_THREADS=1",
        "BLIS_NUM_THREADS=1",
        "NUMEXPR_NUM_THREADS=1",
        str(BINARY),
        "native-f4",
        str(MANIFEST),
        "300",
    ]
    subprocess.run(command, check=False)
    metrics = json.loads((out / "metrics.json").read_text())
    report = json.loads((out / "stdout.json").read_text())
    process = metrics["metrics"]
    extra = report.get("cost", {}).get("extra", {})
    print(
        json.dumps(
            {
                "name": name,
                "order": order,
                "returncode": metrics["returncode"],
                "status": report.get("status"),
                "wall_seconds": process["wall_seconds"],
                "core_seconds": process["total_core_seconds"],
                "peak_rss_bytes": process["peak_rss_bytes"],
                "ops": report.get("cost", {}).get("ops"),
                "build_ns": extra.get("build_ns"),
                "eliminate_ns": extra.get("eliminate_ns"),
                "pair_update_ns": extra.get("pair_update_ns"),
                "equation_fingerprint": report.get("solver_equations_blake3"),
            },
            sort_keys=True,
        ),
        flush=True,
    )

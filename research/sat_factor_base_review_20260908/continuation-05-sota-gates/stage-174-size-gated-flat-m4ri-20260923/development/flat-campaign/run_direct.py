#!/usr/bin/env python3
import json
from pathlib import Path
import subprocess
import sys

REPO = Path("/Volumes/SSD990/crypto-kic-stage174-native-f4")
ROOT = Path("/Volumes/SSD990/kic-stage174-flat-dev")
BINARY = Path("/Volumes/SSD990/koblitz-native-f4-build23-e51efb21/bin/koblitz_pdp_backend")
MANIFEST = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-170-parallel-fixed-x1-construction-20260923/selected-run/tasks/000000-b-421e22a9c1c3b9d56396c8bbd0e46185bbebc32de0306f2bee6e9585703c2be4/instance/manifest.json"
METER = REPO / "scripts/process_meter.py"

for index in range(3):
    out = ROOT / f"clean-direct-r{index + 1}"
    out.mkdir(parents=True, exist_ok=False)
    subprocess.run([
        sys.executable, str(METER),
        "--cwd", str(REPO), "--timeout", "60",
        "--stdout", str(out / "stdout.json"),
        "--stderr", str(out / "stderr.txt"),
        "--metrics", str(out / "metrics.json"),
        "--exclusive-create", "--",
        "/usr/bin/env",
        "VECLIB_MAXIMUM_THREADS=1", "OPENBLAS_NUM_THREADS=1", "OMP_NUM_THREADS=1",
        "MKL_NUM_THREADS=1", "BLIS_NUM_THREADS=1", "NUMEXPR_NUM_THREADS=1",
        str(BINARY), "direct-mitm", str(MANIFEST),
    ], check=False)
    metrics = json.loads((out / "metrics.json").read_text())
    report = json.loads((out / "stdout.json").read_text())
    print(json.dumps({
        "run": index + 1,
        "returncode": metrics["returncode"],
        "wall_seconds": metrics["metrics"]["wall_seconds"],
        "core_seconds": metrics["metrics"]["total_core_seconds"],
        "peak_rss_bytes": metrics["metrics"]["peak_rss_bytes"],
        "status": report.get("status"),
        "source_instance_id": report.get("source_instance_id"),
    }, sort_keys=True), flush=True)

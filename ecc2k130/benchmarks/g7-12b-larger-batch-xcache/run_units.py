#!/usr/bin/env python3
from pathlib import Path
import hashlib
import json
import subprocess
import time

from build import ROOT, OUT

rows = []
for label in ("batch24-cache6-min2", "batch32-cache8-min2"):
    binary = ROOT / "build" / f"g7-{label}-sharedx-test"
    log = OUT / f"{label}-sharedx.log"
    assert not log.exists()
    start = time.monotonic()
    with log.open("x") as stream:
        process = subprocess.run([str(binary)], stdout=stream,
                                 stderr=subprocess.STDOUT, timeout=600)
    text = log.read_text()
    row = {"label": label, "returncode": process.returncode,
           "elapsed_seconds": time.monotonic() - start,
           "passed": process.returncode == 0 and text.count("PASS:") >= 3
                     and "FAIL:" not in text,
           "log_sha256": hashlib.sha256(log.read_bytes()).hexdigest()}
    rows.append(row)
    print(json.dumps(row), flush=True)
    assert row["passed"], text[-4000:]
(OUT / "unit-runs.json").write_text(json.dumps(rows, indent=2) + "\n")

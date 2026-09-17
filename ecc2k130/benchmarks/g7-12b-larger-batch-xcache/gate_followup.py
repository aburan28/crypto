#!/usr/bin/env python3
from pathlib import Path
import hashlib
import json
import os
import subprocess
import time

from build import ROOT, OUT

rows = []
for label, batch in (("batch24-cache6-min3", 24),
                     ("batch32-cache8-min3", 32)):
    unit = ROOT / "build" / f"g7-{label}-sharedx-test"
    log = OUT / f"{label}-sharedx.log"
    start = time.monotonic()
    with log.open("x") as stream:
        process = subprocess.run([str(unit)], stdout=stream,
                                 stderr=subprocess.STDOUT, timeout=600)
    text = log.read_text()
    row = {"label": label, "gate": "actual-walk-shared-x",
           "returncode": process.returncode,
           "elapsed_seconds": time.monotonic() - start,
           "passed": process.returncode == 0 and text.count("PASS:") >= 3
                     and "FAIL:" not in text,
           "log_sha256": hashlib.sha256(log.read_bytes()).hexdigest()}
    rows.append(row); print(json.dumps(row), flush=True); assert row["passed"]
    binary = ROOT / "build" / f"g7-{label}"
    for tool in ("memcheck", "initcheck", "synccheck"):
        log = OUT / f"{label}-{tool}.log"
        command = ["compute-sanitizer", "--tool", tool, str(binary),
                   "--packed", "--threads", "129", "--steps", "2",
                   "--launches", "1", "--verify", "0", "--run-id",
                   str(61301 + batch), "--bench"]
        start = time.monotonic()
        with log.open("x") as stream:
            process = subprocess.run(command, stdout=stream,
                                     stderr=subprocess.STDOUT, timeout=300,
                                     env=dict(os.environ,
                                              CUDA_DISABLE_PTX_JIT="1",
                                              OMP_NUM_THREADS="8"))
        text = log.read_text()
        row = {"label": label, "gate": tool, "command": command,
               "returncode": process.returncode,
               "elapsed_seconds": time.monotonic() - start,
               "passed": process.returncode == 0
                         and "ERROR SUMMARY: 0 errors" in text
                         and "0 dropped" in text,
               "log_sha256": hashlib.sha256(log.read_bytes()).hexdigest()}
        rows.append(row); print(json.dumps(row), flush=True); assert row["passed"]
(OUT / "followup-gates.json").write_text(json.dumps(rows, indent=2) + "\n")

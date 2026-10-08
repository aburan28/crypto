#!/usr/bin/env python3
from pathlib import Path
import hashlib
import json
import os
import subprocess
import time

from build import ROOT, OUT

rows = []
for label, batch in (("batch24-cache6-min2", 24),
                     ("batch32-cache8-min2", 32)):
    binary = ROOT / "build" / f"g7-{label}"
    for tool in ("memcheck", "initcheck", "synccheck"):
        log = OUT / f"{label}-{tool}.log"
        assert not log.exists()
        command = ["compute-sanitizer", "--tool", tool, str(binary),
                   "--packed", "--threads", "129", "--steps", "2",
                   "--launches", "1", "--verify", "0", "--run-id",
                   str(61101 + batch), "--bench"]
        start = time.monotonic()
        with log.open("x") as stream:
            process = subprocess.run(command, stdout=stream,
                                     stderr=subprocess.STDOUT, timeout=300,
                                     env=dict(os.environ,
                                              CUDA_DISABLE_PTX_JIT="1",
                                              OMP_NUM_THREADS="8"))
        text = log.read_text()
        row = {"label": label, "tool": tool, "command": command,
               "returncode": process.returncode,
               "elapsed_seconds": time.monotonic() - start,
               "passed": process.returncode == 0
                         and "ERROR SUMMARY: 0 errors" in text
                         and "0 dropped" in text,
               "log_sha256": hashlib.sha256(log.read_bytes()).hexdigest()}
        rows.append(row)
        print(json.dumps(row), flush=True)
        assert row["passed"], text[-5000:]
(OUT / "sanitizers.json").write_text(json.dumps(rows, indent=2) + "\n")

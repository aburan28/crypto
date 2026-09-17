#!/usr/bin/env python3
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import json
import subprocess
import time

from build import ROOT, OUT, SRC

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def build(label):
    receipt = json.loads((OUT / f"{label}-build.json").read_text())
    command = list(receipt["command"])
    command[command.index(str(SRC / "src/main.cu"))] = str(SRC / "src/testsharedxcuda.cu")
    old_binary = str(ROOT / receipt["binary"])
    binary = ROOT / "build" / f"g7-{label}-sharedx-test"
    command[command.index(old_binary)] = str(binary)
    log = OUT / f"{label}-sharedx-build.log"
    assert not binary.exists() and not log.exists()
    start = time.monotonic()
    with log.open("x") as stream:
        process = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                                 timeout=900)
    row = {"label": label, "command": command,
           "returncode": process.returncode,
           "elapsed_seconds": time.monotonic() - start,
           "log_sha256": sha(log), "binary": str(binary.relative_to(ROOT)),
           "binary_sha256": sha(binary) if binary.exists() else None}
    (OUT / f"{label}-sharedx-build.json").write_text(json.dumps(row, indent=2) + "\n")
    print(json.dumps(row), flush=True)
    return row

labels = ("batch24-cache6-min3", "batch32-cache8-min3")
with ThreadPoolExecutor(max_workers=2) as pool:
    rows = list(pool.map(build, labels))
assert all(row["returncode"] == 0 for row in rows)

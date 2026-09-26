#!/usr/bin/env python3
"""Replay the tiny real-WDSat receipt without invoking a solver binary."""

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RECEIPT = json.loads((HERE / "receipt.json").read_text())
FILES = {
    "build_config_sha256": "config.h",
    "input_sha256": "contradiction.anf",
    "stdout_sha256": "stdout.txt",
    "stderr_sha256": "stderr.txt",
}
for field, filename in FILES.items():
    got = hashlib.sha256((HERE / filename).read_bytes()).hexdigest()
    assert got == RECEIPT[field], (field, got, RECEIPT[field])

assert (HERE / "contradiction.anf").read_text() == (
    "p cnf 1 2\nx 1 T 0\nx 1 0\n"
)
# WDSat ANF rows assert odd parity. The first equation forces x=0;
# the second forces x=1. Re-evaluate both assignments independently.
solutions = [x for x in (0, 1) if (x ^ 1) == 1 and x == 1]
assert solutions == RECEIPT["exhaustive_solution_bits"] == []
assert RECEIPT["returncode"] == 0
assert (HERE / "stdout.txt").read_text() == "UNSAT on XORGAUSS init\n"
assert (HERE / "stderr.txt").read_text() == ""
print("WDSat contradiction receipt: hashes, status, stdout and exhaustive UNSAT verified")

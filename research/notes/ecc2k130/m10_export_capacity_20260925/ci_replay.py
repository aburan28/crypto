#!/usr/bin/env python3
"""Hash-only release hold plus independent basis and n3 byte/branch preflight."""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_json(name: str) -> dict:
    output = subprocess.check_output([sys.executable, str(HERE / name)], cwd=ROOT, text=True)
    return json.loads(output)


def main() -> int:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    spec = json.loads((HERE / "INPUT.json").read_text())
    assert frozen["schema"] == "ecc2k130-m10-capacity-frozen-v1"
    assert frozen["source_base_head"] == "708e3c884c1707affbd186288639b67972047d2c"
    assert frozen["release_main_head"] is None or re.fullmatch(
        r"[0-9a-f]{40}", frozen["release_main_head"])
    assert spec["schema"] == "ecc2k130-m10-complete-chain-capacity-input-v1"
    assert spec["field_degree"] == 131 and spec["normal_beta"] == 3
    assert [row["dimensions"] for row in spec["arms"]] == (
        [[13] * 10, [14] + [13] * 9])
    assert spec["caps"] == frozen["caps"]
    for name, expected in frozen["source_sha256"].items():
        assert sha(HERE / name) == expected, name
    for relative, expected in frozen["input_sha256"].items():
        assert sha(ROOT / relative) == expected, relative
    for name, expected in spec["source_results_sha256"].items():
        assert sha(HERE / "inputs" / name) == expected, name
    basis = run_json("basis_verify.py")
    toy = run_json("selftest.py")
    assert basis["status"] == "PASS" and toy["status"] == "PASS"
    assert basis["arms"] == frozen["basis_expected"]
    assert {key: toy[key] for key in frozen["toy_expected"]} == frozen["toy_expected"]
    print(json.dumps({"status": "PASS", "release_main_head": frozen["release_main_head"],
                      "freeze_sha256": sha(HERE / "FROZEN.json"),
                      "basis_sha256": sha(HERE / "basis_verify.py"),
                      "toy_dag_nodes": toy["dag_nodes"]},
                     sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    sys.exit(main())

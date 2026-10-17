#!/usr/bin/env python3
"""Bind the archived original eight-arm N83 range to exact IC1 candidates."""

import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("n83_range_materializer", HERE / "materialize.py")
materialize = importlib.util.module_from_spec(spec)
spec.loader.exec_module(materialize)

materialize.RUN_ID = 38895743960
materialize.SOURCE_COMMIT = "e83a78be0b74d22fb34e315006d47e6bd82a1e4c"
materialize.ARTIFACT_ID = 12046277866
materialize.ARCHIVE_SHA256 = "41d4aa0d623967ab9c62cc3ded653e01e9ec049e78596d691d87b7f40f4cea45"
materialize.ANALYSIS_SHA256 = "5db680f0423f1e4b4c4138303135e826e543cb90f3b427147977ee9be192927c"
materialize.panel.SOURCE_COMMIT = materialize.SOURCE_COMMIT
materialize.runner.ARMS = (
    ("k2543root22_a", 2543, "root22_dynamic", 900),
    ("k1024root22", 1024, "root22_dynamic", 900),
    ("k2048root22", 2048, "root22_dynamic", 900),
    ("k2543root11", 2543, "root11_dynamic", 900),
    ("k3072root11", 3072, "root11_dynamic", 900),
    ("k3450root11fixed", 3450, "root11_fixed", 900),
    ("k3584root11fixed", 3584, "root11_fixed", 900),
    ("k2543root22_b", 2543, "root22_dynamic", 900),
)

if __name__ == "__main__":
    materialize.main()

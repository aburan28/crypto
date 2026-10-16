#!/usr/bin/env python3
"""Bind the archived current-source N83 range to exact IC1 candidates."""

import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("n83_range_materializer", HERE / "materialize.py")
materialize = importlib.util.module_from_spec(spec)
spec.loader.exec_module(materialize)

materialize.RUN_ID = 38897449665
materialize.SOURCE_COMMIT = "f0a5fce29a441566811f9afc260bc30a59f18144"
materialize.ARTIFACT_ID = 12050929059
materialize.ARCHIVE_SHA256 = "c6fbc8567699767a89095aae793ffe2394658cdfaaa33a76e10e8243ca85e594"
materialize.ANALYSIS_SHA256 = "49e71b6542b655b4097700fb6e8d98d595da23acac17061c9d1167f60b2f9965"
materialize.panel.SOURCE_COMMIT = materialize.SOURCE_COMMIT

if __name__ == "__main__":
    materialize.main()

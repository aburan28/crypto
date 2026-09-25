#!/usr/bin/env python3
"""Rehydrate the frozen base and replay the one-point IC witness from committed data."""

import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
IC = REPO / "research/sat_factor_base_review_20260908/autolab_orbit_extract_20260924"
SINGLE = HERE / "ic_point_1_raw"
BASE_GZ = IC / "independent_replay_20260924_codex/base_header.jsonl.gz"
TRAINING = IC / "cold_rank_20260924_codex"
VERIFIER = IC / "recover_point_batch.py"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    addendum = json.loads((HERE / "ic_single_addendum_protocol.json").read_text())
    manifest = json.loads((SINGLE / "manifest.json").read_text())
    assert digest(BASE_GZ) == manifest["base_gzip_sha256"]
    assert digest(SINGLE / "target_points.jsonl") == addendum["point_sha256"]
    assert manifest["source_sha256"] == addendum["point_source_sha256"]
    assert manifest["exe_sha256"] == addendum["point_executable_sha256"]
    assert digest(TRAINING / "manifest.json") == addendum["training_manifest_sha256"]
    assert digest(TRAINING / "validation.json") == addendum["training_validation_sha256"]
    with tempfile.TemporaryDirectory() as directory:
        temp = Path(directory)
        for path in SINGLE.iterdir():
            if path.is_file() and path.name != "base_header.jsonl":
                shutil.copy2(path, temp / path.name)
        with gzip.open(BASE_GZ, "rb") as packed:
            (temp / "base_header.jsonl").write_bytes(packed.read())
        process = subprocess.run(
            [sys.executable, str(VERIFIER), "--training", str(TRAINING),
             "--point-batch", str(temp)],
            cwd=REPO, capture_output=True, text=True,
        )
        assert process.returncode == 0, process.stderr
        assert (temp / "validation.json").read_bytes() == (SINGLE / "validation.json").read_bytes()
    report = json.loads((SINGLE / "validation.json").read_text())
    assert report["targets_recovered"] == 1
    assert report["rows"][0]["recovered_scalar"] == report["rows"][0]["validator_scalar"]
    print(json.dumps({
        "single_point_independent_replay": True,
        "validation_sha256": digest(SINGLE / "validation.json"),
        "recovered_scalar": report["rows"][0]["recovered_scalar"],
    }, sort_keys=True))


if __name__ == "__main__":
    run()

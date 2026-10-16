#!/usr/bin/env python3
"""Verify exact cloud ZIP and replay the N83 K-grid before sealing finalists."""

import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCREEN_RUN = ROOT / "research/n83_e2e_factor_base_20261016/run.py"
SCREEN_VERIFY = ROOT / "research/n83_e2e_factor_base_20261016/verify.py"

spec = importlib.util.spec_from_file_location("n83_e2e_screen_materialize", SCREEN_RUN)
screen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(screen)


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def write_identical_or_new(path, blob):
    if path.exists():
        if path.read_bytes() != blob:
            raise ValueError(f"existing frozen file differs: {path}")
    else:
        path.write_bytes(blob)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path, help="exact GitHub artifact ZIP")
    parser.add_argument("extract_dir", type=Path,
                        help="new durable directory for the raw extraction and replay")
    parser.add_argument("--github-sha256", required=True,
                        help="sha256 from the provider's artifact digest field")
    parser.add_argument("--run-id", required=True, type=int)
    parser.add_argument("--head-sha", required=True)
    args = parser.parse_args()
    archive = args.archive.resolve()
    extract = args.extract_dir.resolve()
    holdout = importlib.util.spec_from_file_location(
        "n83_e2e_holdout_materialize", HERE / "run.py")
    holdout_module = importlib.util.module_from_spec(holdout)
    holdout.loader.exec_module(holdout_module)
    if (args.run_id != holdout_module.SCREEN_RUN_ID
            or args.head_sha != holdout_module.SCREEN_HEAD_SHA):
        raise ValueError("provider run ID or executed commit differs from frozen screen")
    provider_sha = args.github_sha256.removeprefix("sha256:")
    if len(provider_sha) != 64 or any(c not in "0123456789abcdef" for c in provider_sha):
        raise ValueError("provider SHA-256 must be lowercase 64-digit hex")
    if digest(archive) != provider_sha:
        raise ValueError("downloaded artifact differs from GitHub digest")
    if extract.exists():
        raise ValueError("raw extraction directory already exists")
    extract.mkdir(parents=True)
    with zipfile.ZipFile(archive) as saved:
        if saved.testzip() is not None:
            raise ValueError("artifact ZIP integrity failed")
        names = saved.namelist()
        for name in names:
            parts = Path(name).parts
            if not parts or Path(name).is_absolute() or ".." in parts:
                raise ValueError("artifact contains an unsafe path")
        required = {"results/analysis.json", "results/host.json",
                    "results/source.json", "results/preflight.json"}
        if not required.issubset(names):
            raise ValueError("artifact lacks required result receipts")
        original_analysis = saved.read("results/analysis.json")
        saved.extractall(extract)
    source = json.loads((extract / "results/source.json").read_text())
    host = json.loads((extract / "results/host.json").read_text())
    if (source["source_sha256"] != screen.SOURCE_SHA256
            or source["runner_sha256"] != digest(SCREEN_RUN)
            or host["git_sha"] != args.head_sha):
        raise ValueError("artifact source or executed commit differs")
    command = [sys.executable, str(SCREEN_VERIFY), str(extract / "results")]
    replay = subprocess.run(command, capture_output=True, text=True,
                            check=False, cwd=ROOT)
    (extract / "local_replay.stdout.txt").write_text(replay.stdout)
    (extract / "local_replay.stderr.txt").write_text(replay.stderr)
    if replay.returncode != 0:
        raise ValueError(f"independent replay failed with exit {replay.returncode}")
    if (extract / "results/analysis.json").read_bytes() != original_analysis:
        raise ValueError("local analysis differs from exact cloud analysis")
    analysis = json.loads(original_analysis)
    if (analysis["schema"] != "n83.e2e-factor-base-analysis/v1"
            or analysis["status"] != "SCREEN_COMPLETE"
            or analysis["evidence_status"] != "PASS_INDEPENDENT_REPLAY"
            or analysis["aa_control"]["status"] != "PASS_10_PERCENT"
            or len(analysis["cold_process_ms_by_k"]) != len(screen.K_GRID)):
        raise ValueError("screen result did not admit finalist selection")
    challengers = [int(k) for k in analysis["descriptive_order"]
                   if int(k) != screen.CONTROL_K][:3]
    if len(challengers) != 3:
        raise ValueError("screen has fewer than three other completed Ks")
    selection = {
        "schema": "n83.e2e-holdout-selection/v1",
        "screen_run_id": args.run_id,
        "screen_source_commit": args.head_sha,
        "screen_source_sha256": screen.SOURCE_SHA256,
        "screen_archive_sha256": provider_sha,
        "screen_analysis_sha256": hashlib.sha256(original_analysis).hexdigest(),
        "challenger_ks": challengers,
    }
    selection_blob = (json.dumps(selection, indent=2, sort_keys=True) + "\n").encode()
    write_identical_or_new(HERE / "screen_analysis.json", original_analysis)
    write_identical_or_new(HERE / "selection.json", selection_blob)
    print(json.dumps(selection, sort_keys=True))


if __name__ == "__main__":
    main()

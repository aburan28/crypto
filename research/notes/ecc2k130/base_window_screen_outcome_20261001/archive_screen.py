#!/usr/bin/env python3
"""Seal the one-shot n41 screen, including every raw child and failure file."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_ir_ledger_20260930"))
from archive import equivalent, pack  # noqa: E402
from run_panel import sha  # noqa: E402
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/base_window_screen_20261001"))
from analyze_screen import analyze  # noqa: E402
from run_screen import FROZEN, INPUT_RECEIPT, config, schedule  # noqa: E402

RUN_ID = 36827211814
RUN_URL = f"https://github.com/aburan28/crypto/actions/runs/{RUN_ID}"
SOURCE_HEAD = "4978478120ae83a322c221e59d65a358d5510f29"
FROZEN_SHA256 = "6487e1c30ef1ae1db5cb99fd39297150b1a327baaff727e3966d3115bc88e2a7"
RAW_ROOT_NAME = "base-window-raw-36827211814"


def read(path: Path) -> dict:
    return json.loads(path.read_text())


def copy_record(source: Path, target: Path) -> dict:
    shutil.copyfile(source, target)
    return {"path": target.name, "bytes": target.stat().st_size, "sha256": sha(target)}


def archive(raw_root: Path, hosted_replay_path: Path, local_replay_path: Path,
            analysis_path: Path, out: Path) -> dict:
    assert raw_root.is_dir() and raw_root.name == RAW_ROOT_NAME
    assert all(path.is_file() for path in (hosted_replay_path, local_replay_path,
                                          analysis_path))
    assert not out.exists() and not out.is_relative_to(raw_root)
    assert sha(FROZEN) == FROZEN_SHA256
    cfg = config()
    report_path = raw_root / "base-window-screen/screen_run.json"
    report = read(report_path)
    assert report["schema"] == "ecc2k130-base-window-screen-run-v1"
    assert report["config"] == cfg and report["host"]["git_head"] == SOURCE_HEAD
    assert report["frozen_sha256"] == sha(FROZEN)
    assert report["input_receipt_sha256"] == sha(INPUT_RECEIPT)
    assert report["plan"] == [
        {"block": block, "arm": arm} for block, arm in schedule(cfg)
    ]
    build = read(raw_root / "base-window-screen/BUILD_RECEIPT.json")
    assert build["status"] == "PASS" and build["profile"] == "release"
    assert sha(raw_root / "build.stdout.txt") == build["build_stdout_sha256"]
    assert sha(raw_root / "build.stderr.txt") == build["build_stderr_sha256"]
    assert sha(raw_root / "isolation.jsonl") == sha(
        raw_root / "base-window-screen/isolation.jsonl"
    )
    hosted, local, analysis_row = (read(path) for path in (
        hosted_replay_path, local_replay_path, analysis_path
    ))
    assert hosted["schema"] == local["schema"] == "ecc2k130-base-window-screen-replay-v1"
    assert hosted["screen_run_sha256"] == local["screen_run_sha256"] == sha(report_path)
    assert hosted["isolation_sha256"] == local["isolation_sha256"] == sha(
        raw_root / "isolation.jsonl"
    )
    assert equivalent(hosted, local), "hosted and macOS independent replays differ"
    assert equivalent(analysis_row, analyze(report, hosted))

    (out / "raw").mkdir(parents=True)
    (out / "receipts").mkdir()
    archive_path = out / "raw/base-window-screen.tar.gz"
    member_hashes = pack(raw_root, archive_path)
    assert len(member_hashes) >= 400
    records = {
        "hosted_replay": copy_record(
            hosted_replay_path, out / "receipts/hosted_replay.json"
        ),
        "macos_replay": copy_record(
            local_replay_path, out / "receipts/macos_replay.json"
        ),
        "analysis": copy_record(
            analysis_path, out / "receipts/analysis.json"
        ),
    }
    for entry in records.values():
        entry["path"] = "receipts/" + entry["path"]
    manifest = {
        "schema": "ecc2k130-base-window-screen-hosted-archive-v1",
        "github_run_id": RUN_ID,
        "github_run_url": RUN_URL,
        "measured_main_head": SOURCE_HEAD,
        "source_lock_merge_commit": read(FROZEN)["source_lock"]["merge_commit"],
        "input_freeze_sha256": FROZEN_SHA256,
        "input_receipt_sha256": sha(INPUT_RECEIPT),
        "protocol_config_sha256": sha(HERE.parent / "base_window_screen_20261001/CONFIG.json"),
        "raw_root_name": RAW_ROOT_NAME,
        "raw_path": "raw/base-window-screen.tar.gz",
        "raw_bytes": archive_path.stat().st_size,
        "raw_sha256": sha(archive_path),
        "raw_member_sha256": member_hashes,
        "run_json_sha256": sha(report_path),
        "build_receipt_sha256": sha(raw_root / "base-window-screen/BUILD_RECEIPT.json"),
        "isolation_sha256": sha(raw_root / "isolation.jsonl"),
        "producer_status": report["status"],
        "hosted_replay_status": hosted["status"],
        "timing_eligible": hosted.get("timing_eligible"),
        "decision": analysis_row["decision"],
        "records": records,
        "extraction": "tar -xzf raw/base-window-screen.tar.gz -C DEST",
        "independent_replay": (
            "python3 research/notes/ecc2k130/base_window_screen_20261001/"
            "verify_screen.py --run-dir DEST/base-window-raw-36827211814/"
            "base-window-screen --out DEST/replay.json --relocated"
        ),
    }
    (out / "MANIFEST.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--hosted-replay", type=Path, required=True)
    parser.add_argument("--local-replay", type=Path, required=True)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = archive(*(path.resolve() for path in (
        args.raw_root, args.hosted_replay, args.local_replay, args.analysis, args.out
    )))
    print(json.dumps({
        "producer_status": result["producer_status"],
        "hosted_replay_status": result["hosted_replay_status"],
        "timing_eligible": result["timing_eligible"],
        "decision": result["decision"],
        "raw_bytes": result["raw_bytes"],
        "raw_sha256": result["raw_sha256"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

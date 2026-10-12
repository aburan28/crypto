#!/usr/bin/env python3
"""Rehash the sealed raw run and independently replay its complete verdict."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import tarfile
import tempfile
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_ir_ledger_20260930"))
from archive import equivalent  # noqa: E402
from run_panel import sha  # noqa: E402
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/base_window_screen_20261001"))
from analyze_screen import analyze  # noqa: E402
from run_screen import FROZEN, INPUT_RECEIPT, config  # noqa: E402
from verify_screen import verify  # noqa: E402
from archive_screen import FROZEN_SHA256, RAW_ROOT_NAME, RUN_ID, RUN_URL, SOURCE_HEAD  # noqa: E402

MAX_MEMBER_BYTES = 200 * 1024 * 1024
MAX_TOTAL_BYTES = 1024 * 1024 * 1024


def read(path: Path) -> dict:
    return json.loads(path.read_text())


def replay(archive_dir: Path) -> dict:
    manifest_path = archive_dir / "MANIFEST.json"
    manifest = read(manifest_path)
    assert manifest["schema"] == "ecc2k130-base-window-screen-hosted-archive-v1"
    assert manifest["github_run_id"] == RUN_ID and manifest["github_run_url"] == RUN_URL
    assert manifest["measured_main_head"] == SOURCE_HEAD
    assert manifest["raw_root_name"] == RAW_ROOT_NAME
    assert manifest["input_freeze_sha256"] == sha(FROZEN) == FROZEN_SHA256
    assert manifest["input_receipt_sha256"] == sha(INPUT_RECEIPT)
    assert manifest["protocol_config_sha256"] == sha(
        ROOT / "research/notes/ecc2k130/base_window_screen_20261001/CONFIG.json"
    )
    raw = archive_dir / manifest["raw_path"]
    assert raw.is_file() and raw.stat().st_size == manifest["raw_bytes"]
    assert sha(raw) == manifest["raw_sha256"]
    saved = {}
    for key, entry in manifest["records"].items():
        path = archive_dir / entry["path"]
        assert path.is_file() and path.stat().st_size == entry["bytes"]
        assert sha(path) == entry["sha256"]
        saved[key] = read(path)
    with tempfile.TemporaryDirectory(prefix="base-window-sealed-replay-") as tmp:
        dest = Path(tmp)
        expected = manifest["raw_member_sha256"]
        seen = set()
        total = 0
        with tarfile.open(raw, "r:gz") as packed:
            for member in packed:
                assert member.isfile() and not member.issym() and not member.islnk()
                parts = Path(member.name).parts
                assert len(parts) >= 2 and parts[0] == RAW_ROOT_NAME
                assert all(part not in ("", ".", "..") for part in parts)
                relative = str(Path(*parts[1:]))
                assert relative in expected and relative not in seen
                assert 0 <= member.size <= MAX_MEMBER_BYTES
                total += member.size
                assert total <= MAX_TOTAL_BYTES
                target = dest / Path(*parts)
                target.parent.mkdir(parents=True, exist_ok=True)
                with packed.extractfile(member) as source, target.open("xb") as output:
                    while chunk := source.read(1024 * 1024):
                        output.write(chunk)
                assert target.stat().st_size == member.size
                assert sha(target) == expected[relative], relative
                seen.add(relative)
        assert seen == set(expected) and len(seen) >= 400
        raw_root = dest / RAW_ROOT_NAME
        run_dir = raw_root / "base-window-screen"
        report_path = run_dir / "screen_run.json"
        report = read(report_path)
        assert report["schema"] == "ecc2k130-base-window-screen-run-v1"
        assert report["config"] == config() and report["host"]["git_head"] == SOURCE_HEAD
        assert sha(report_path) == manifest["run_json_sha256"]
        assert sha(run_dir / "BUILD_RECEIPT.json") == manifest["build_receipt_sha256"]
        assert sha(run_dir / "isolation.jsonl") == manifest["isolation_sha256"]
        independently_replayed = verify(run_dir, relocated=True)
        assert equivalent(independently_replayed, saved["hosted_replay"])
        assert equivalent(independently_replayed, saved["macos_replay"])
        independently_analyzed = analyze(report, independently_replayed)
        assert equivalent(independently_analyzed, saved["analysis"])
        assert manifest["producer_status"] == report["status"]
        assert manifest["hosted_replay_status"] == independently_replayed["status"]
        assert manifest["timing_eligible"] == independently_replayed.get("timing_eligible")
        assert manifest["decision"] == independently_analyzed["decision"]
    return {
        "schema": "ecc2k130-base-window-sealed-replay-v1",
        "status": "PASS",
        "github_run_id": RUN_ID,
        "manifest_sha256": sha(manifest_path),
        "raw_sha256": manifest["raw_sha256"],
        "members_rehashed": len(expected),
        "producer_status": report["status"],
        "replay_status": independently_replayed["status"],
        "timing_eligible": independently_replayed.get("timing_eligible"),
        "decision": independently_analyzed["decision"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists()
    result = replay(args.archive.resolve())
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

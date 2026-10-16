#!/usr/bin/env python3
"""Fail-closed supervisor for the remaining p61 selected-class vertex census."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import shutil
import subprocess
from typing import Any


TOTAL_COLUMNS = 19_894_128


@dataclass(frozen=True)
class SegmentPaths:
    segment: Path
    checkpoint: Path
    generation_receipt: Path
    verification_receipt: Path
    admission_manifest: Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise RuntimeError(f"{path}: expected a JSON object")
    return value


def segment_paths(root: Path, prefix: str, start: int, end: int) -> SegmentPaths:
    span = f"{start:012d}-{end:012d}"
    return SegmentPaths(
        segment=root / "segments" / f"{prefix}-segment-{span}.jsonl.gz",
        checkpoint=root / "checkpoints" / f"{prefix}-checkpoint-{end:012d}.json",
        generation_receipt=root / "receipts" / f"{prefix}-generate-{span}.json",
        verification_receipt=root / "verification" / f"{prefix}-verify-{span}.json",
        admission_manifest=root / "manifests" / f"{prefix}-admission-{span}.json",
    )


def require_exact_file_set(paths: list[Path], label: str) -> bool:
    present = [path.exists() for path in paths]
    if any(present) and not all(present):
        detail = ", ".join(f"{path}={'present' if exists else 'missing'}" for path, exists in zip(paths, present))
        raise RuntimeError(f"partial {label} evidence; refusing to continue: {detail}")
    return all(present)


def validate_generation(
    paths: SegmentPaths,
    start: int,
    end: int,
    previous_checkpoint: Path,
    source_commit: str,
    class_receipt_sha256: str,
) -> dict[str, Any]:
    receipt = load_json(paths.generation_receipt)
    checkpoint = load_json(paths.checkpoint)
    if receipt.get("status") != "pass":
        raise RuntimeError(f"{paths.generation_receipt}: generation status is not pass")
    inner = receipt.get("receipt", {})
    summary = inner.get("summary", {})
    if summary.get("x_start") != start or summary.get("x_end") != end:
        raise RuntimeError("generation receipt range mismatch")
    if summary.get("node_records") != (end - start) * 18:
        raise RuntimeError("generation receipt record-count mismatch")
    if inner.get("source_commit") != source_commit:
        raise RuntimeError("generation source-commit mismatch")
    if checkpoint.get("next_x") != end or checkpoint.get("x_start") != start:
        raise RuntimeError("checkpoint range mismatch")
    if checkpoint.get("source_commit") != source_commit:
        raise RuntimeError("checkpoint source-commit mismatch")
    if checkpoint.get("class_number_receipt_sha256") != class_receipt_sha256:
        raise RuntimeError("checkpoint class-number receipt mismatch")
    previous_digest = sha256_file(previous_checkpoint)
    if checkpoint.get("previous_checkpoint_sha256") != previous_digest:
        raise RuntimeError("checkpoint predecessor digest mismatch")
    segment_digest = sha256_file(paths.segment)
    if checkpoint.get("segment_sha256") != segment_digest:
        raise RuntimeError("checkpoint segment digest mismatch")
    if receipt.get("segment_digest", {}).get("sha256") != segment_digest:
        raise RuntimeError("generation receipt segment digest mismatch")
    checkpoint_digest = sha256_file(paths.checkpoint)
    if receipt.get("checkpoint_digest", {}).get("sha256") != checkpoint_digest:
        raise RuntimeError("generation receipt checkpoint digest mismatch")
    return receipt


def validate_verification(paths: SegmentPaths, start: int, end: int) -> dict[str, Any]:
    receipt = load_json(paths.verification_receipt)
    if receipt.get("status") != "pass":
        raise RuntimeError(f"{paths.verification_receipt}: verification status is not pass")
    inner = receipt.get("receipt", {})
    summary = inner.get("summary", {})
    if summary.get("x_start") != start or summary.get("x_end") != end:
        raise RuntimeError("verification receipt range mismatch")
    if summary.get("node_records") != (end - start) * 18:
        raise RuntimeError("verification receipt record-count mismatch")
    if inner.get("audit_points") != 2 or inner.get("audit_seed_x") != 7:
        raise RuntimeError("verification audit settings mismatch")
    if receipt.get("segment_digest", {}).get("sha256") != sha256_file(paths.segment):
        raise RuntimeError("verification segment digest mismatch")
    if receipt.get("checkpoint_digest", {}).get("sha256") != sha256_file(paths.checkpoint):
        raise RuntimeError("verification checkpoint digest mismatch")
    return receipt


def write_json_new(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(value, indent=2, sort_keys=True) + "\n"
    with path.open("x", encoding="utf-8") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())


def run(command: list[str]) -> None:
    print(json.dumps({"event": "exec", "argv": command}), flush=True)
    subprocess.run(command, check=True)


def require_mount(path: Path, stage: str) -> None:
    if not path.is_mount():
        raise RuntimeError(f"required mount is absent during {stage}: {path}")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--producer", type=Path, required=True)
    result.add_argument("--verifier", type=Path, required=True)
    result.add_argument("--root", type=Path, required=True)
    result.add_argument("--mount-root", type=Path, required=True)
    result.add_argument("--start-x", type=int, default=16_704)
    result.add_argument("--end-x", type=int, default=TOTAL_COLUMNS)
    result.add_argument("--chunk-columns", type=int, default=16_384)
    result.add_argument("--max-segments", type=int)
    result.add_argument("--threads", type=int, default=8)
    result.add_argument("--prefix", default="campaign-v1")
    result.add_argument("--previous-checkpoint", type=Path, required=True)
    result.add_argument("--generation-source-commit", required=True)
    result.add_argument("--verifier-source-commit", required=True)
    result.add_argument("--class-number-receipt-sha256", required=True)
    result.add_argument("--min-free-bytes", type=int, default=50_000_000_000)
    result.add_argument("--plan-only", action="store_true")
    return result


def validate_args(args: argparse.Namespace) -> None:
    if not (0 <= args.start_x < args.end_x <= TOTAL_COLUMNS):
        raise RuntimeError("require 0 <= start-x < end-x <= 19894128")
    if args.chunk_columns <= 0 or args.threads <= 0:
        raise RuntimeError("chunk-columns and threads must be positive")
    if args.max_segments is not None and args.max_segments <= 0:
        raise RuntimeError("max-segments must be positive")
    if len(args.generation_source_commit) != 40 or len(args.verifier_source_commit) != 40:
        raise RuntimeError("source commits must be full 40-hex identities")
    if len(args.class_number_receipt_sha256) != 64:
        raise RuntimeError("class-number receipt digest must be 64 hex characters")
    try:
        int(args.generation_source_commit, 16)
        int(args.verifier_source_commit, 16)
        int(args.class_number_receipt_sha256, 16)
    except ValueError as error:
        raise RuntimeError("source and receipt identities must be hexadecimal") from error


def main() -> int:
    args = parser().parse_args()
    validate_args(args)
    starts = range(args.start_x, args.end_x, args.chunk_columns)
    plan = []
    for index, start in enumerate(starts):
        if args.max_segments is not None and index >= args.max_segments:
            break
        end = min(start + args.chunk_columns, args.end_x)
        paths = segment_paths(args.root, args.prefix, start, end)
        plan.append({"start": start, "end": end, "records": (end - start) * 18, "paths": {key: str(value) for key, value in paths.__dict__.items()}})
    if args.plan_only:
        print(json.dumps({"status": "plan", "segments": plan}, indent=2, sort_keys=True))
        return 0

    require_mount(args.mount_root, "startup")
    if shutil.disk_usage(args.mount_root).free < args.min_free_bytes:
        raise RuntimeError("free-space floor is not met")
    for binary in (args.producer, args.verifier):
        if not binary.is_file() or not os.access(binary, os.X_OK):
            raise RuntimeError(f"missing executable: {binary}")
    if not args.previous_checkpoint.is_file():
        raise RuntimeError(f"missing predecessor checkpoint: {args.previous_checkpoint}")

    previous_checkpoint = args.previous_checkpoint
    producer_digest = sha256_file(args.producer)
    verifier_digest = sha256_file(args.verifier)
    for item in plan:
        require_mount(args.mount_root, f"segment x={item['start']} startup")
        if shutil.disk_usage(args.mount_root).free < args.min_free_bytes:
            raise RuntimeError(f"free-space floor crossed before x={item['start']}")
        start, end = item["start"], item["end"]
        paths = segment_paths(args.root, args.prefix, start, end)
        generation_exists = require_exact_file_set(
            [paths.segment, paths.checkpoint, paths.generation_receipt], "generation"
        )
        if not generation_exists:
            run([
                str(args.producer), "--threads", str(args.threads), "generate-segment",
                "--x-start", str(start), "--length", str(end - start),
                "--batch-steps", "256", "--audit-points", "1", "--audit-seed-x", "0",
                "--source-commit", args.generation_source_commit,
                "--class-number-receipt-sha256", args.class_number_receipt_sha256,
                "--previous-checkpoint", str(previous_checkpoint),
                "--output", str(paths.segment), "--checkpoint-output", str(paths.checkpoint),
                "--receipt-output", str(paths.generation_receipt),
            ])
        generation = validate_generation(
            paths, start, end, previous_checkpoint, args.generation_source_commit,
            args.class_number_receipt_sha256,
        )
        require_mount(args.mount_root, f"verification x={start}")
        if not paths.verification_receipt.exists():
            run([
                str(args.verifier), "--threads", str(args.threads), "verify-segment",
                "--input", str(paths.segment), "--checkpoint", str(paths.checkpoint),
                "--previous-checkpoint", str(previous_checkpoint),
                "--audit-points", "2", "--audit-seed-x", "7",
                "--receipt-output", str(paths.verification_receipt),
            ])
        verification = validate_verification(paths, start, end)
        require_mount(args.mount_root, f"admission x={start}")
        os.sync()
        if not paths.admission_manifest.exists():
            manifest = {
                "schema": "ec-hardness/p61-campaign-segment-admission/v1",
                "status": "pass",
                "completed_utc": datetime.now(timezone.utc).isoformat(),
                "range": {"x_start": start, "x_end": end, "records": (end - start) * 18},
                "producer": {"path": str(args.producer), "source_commit": args.generation_source_commit, "sha256": producer_digest},
                "verifier": {"path": str(args.verifier), "source_commit": args.verifier_source_commit, "sha256": verifier_digest},
                "segment": {"path": str(paths.segment), "sha256": sha256_file(paths.segment)},
                "checkpoint": {"path": str(paths.checkpoint), "sha256": sha256_file(paths.checkpoint)},
                "generation_receipt": {"path": str(paths.generation_receipt), "sha256": sha256_file(paths.generation_receipt), "elapsed_ms": generation["receipt"]["elapsed_ms"]},
                "verification_receipt": {"path": str(paths.verification_receipt), "sha256": sha256_file(paths.verification_receipt), "elapsed_ms": verification["receipt"]["elapsed_ms"]},
                "durability": {"per_file_sync": "receipt-recorded", "host_sync_before_manifest": "completed"},
            }
            write_json_new(paths.admission_manifest, manifest)
            os.sync()
        else:
            existing = load_json(paths.admission_manifest)
            if existing.get("status") != "pass" or existing.get("range", {}).get("x_end") != end:
                raise RuntimeError(f"invalid existing admission manifest: {paths.admission_manifest}")
        print(json.dumps({"event": "admitted", "x_start": start, "x_end": end, "manifest": str(paths.admission_manifest), "manifest_sha256": sha256_file(paths.admission_manifest)}), flush=True)
        previous_checkpoint = paths.checkpoint
    print(json.dumps({"status": "bounded_run_complete", "segments_admitted": len(plan), "next_x": plan[-1]["end"] if plan else args.start_x}), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

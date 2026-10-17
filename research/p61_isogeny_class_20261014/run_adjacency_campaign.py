#!/usr/bin/env python3
"""Fail-closed supervisor for p61 selected-degree adjacency shards."""

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
import time
from typing import Any, Iterable


DIRECTED_EDGES_PER_VERTEX = 7
KERNEL_COEFFICIENTS_PER_VERTEX = 115
EXPECTED_FINAL_X = 19_894_128


@dataclass(frozen=True)
class VertexAdmission:
    manifest: Path
    x_start: int
    x_end: int
    records: int
    segment: Path
    segment_sha256: str
    verification_receipt: Path
    verification_receipt_sha256: str


@dataclass(frozen=True)
class WindowPaths:
    shard: Path
    generation_receipt: Path
    verification_receipt: Path
    admission_manifest: Path


@dataclass(frozen=True)
class Window:
    source: VertexAdmission
    offset: int
    count: int
    paths: WindowPaths


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


def write_json_new(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(value, indent=2, sort_keys=True) + "\n"
    with path.open("x", encoding="utf-8") as handle:
        handle.write(text)
        handle.flush()
        try:
            os.fsync(handle.fileno())
        except OSError as error:
            if error.errno != 45:
                raise


def require_hex(value: Any, length: int, label: str) -> str:
    if not isinstance(value, str) or len(value) != length:
        raise RuntimeError(f"{label} must be {length} hexadecimal characters")
    try:
        int(value, 16)
    except ValueError as error:
        raise RuntimeError(f"{label} must be hexadecimal") from error
    if value != value.lower():
        raise RuntimeError(f"{label} must use lowercase hexadecimal")
    return value


def require_mount(path: Path, stage: str) -> None:
    if not path.is_mount():
        raise RuntimeError(f"required mount is absent during {stage}: {path}")


def require_free_space(path: Path, floor: int, stage: str) -> int:
    free = shutil.disk_usage(path).free
    if free < floor:
        raise RuntimeError(
            f"free-space floor crossed during {stage}: {free} < {floor} bytes"
        )
    return free


def run(command: list[str]) -> None:
    print(json.dumps({"event": "exec", "argv": command}), flush=True)
    subprocess.run(command, check=True)


def load_vertex_admission(path: Path, verify_files: bool) -> VertexAdmission:
    value = load_json(path)
    if value.get("schema") != "ec-hardness/p61-campaign-segment-admission/v1":
        raise RuntimeError(f"{path}: unsupported vertex admission schema")
    if value.get("status") != "pass":
        raise RuntimeError(f"{path}: vertex admission status is not pass")
    range_value = value.get("range", {})
    x_start = range_value.get("x_start")
    x_end = range_value.get("x_end")
    records = range_value.get("records")
    if not all(isinstance(item, int) for item in (x_start, x_end, records)):
        raise RuntimeError(f"{path}: malformed range")
    if not (0 <= x_start < x_end <= EXPECTED_FINAL_X):
        raise RuntimeError(f"{path}: range lies outside the frozen lattice")
    if records != (x_end - x_start) * 18:
        raise RuntimeError(f"{path}: vertex record count differs from 18 lanes")
    segment_value = value.get("segment", {})
    verification_value = value.get("verification_receipt", {})
    segment = Path(str(segment_value.get("path", "")))
    verification = Path(str(verification_value.get("path", "")))
    segment_sha = require_hex(segment_value.get("sha256"), 64, "segment digest")
    verification_sha = require_hex(
        verification_value.get("sha256"), 64, "verification receipt digest"
    )
    if verify_files:
        for evidence in (segment, verification):
            if not evidence.is_file():
                raise RuntimeError(f"{path}: missing admitted evidence {evidence}")
        if sha256_file(segment) != segment_sha:
            raise RuntimeError(f"{path}: source segment digest replay failed")
        if sha256_file(verification) != verification_sha:
            raise RuntimeError(f"{path}: source verification digest replay failed")
    return VertexAdmission(
        manifest=path,
        x_start=x_start,
        x_end=x_end,
        records=records,
        segment=segment,
        segment_sha256=segment_sha,
        verification_receipt=verification,
        verification_receipt_sha256=verification_sha,
    )


def discover_vertex_admissions(
    directory: Path,
    pattern: str,
    minimum_x: int,
    maximum_x: int,
    verify_files: bool,
) -> list[VertexAdmission]:
    found = []
    for path in sorted(directory.glob(pattern)):
        admission = load_vertex_admission(path, verify_files)
        if admission.x_end <= minimum_x or admission.x_start >= maximum_x:
            continue
        if admission.x_start < minimum_x or admission.x_end > maximum_x:
            raise RuntimeError(
                f"{path}: source segment straddles requested x boundary; use aligned bounds"
            )
        found.append(admission)
    for left, right in zip(found, found[1:]):
        if left.x_end > right.x_start:
            raise RuntimeError(
                f"overlapping admitted vertex ranges: {left.manifest} and {right.manifest}"
            )
    return found


def window_paths(root: Path, prefix: str, source: VertexAdmission, offset: int, count: int) -> WindowPaths:
    source_span = f"x{source.x_start:012d}-{source.x_end:012d}"
    record_span = f"r{offset:012d}-{offset + count:012d}"
    stem = f"{prefix}-{source_span}-{record_span}"
    return WindowPaths(
        shard=root / "shards" / f"{stem}.jsonl.gz",
        generation_receipt=root / "receipts" / f"{stem}-generate.json",
        verification_receipt=root / "verification" / f"{stem}-verify.json",
        admission_manifest=root / "manifests" / f"{stem}-admission.json",
    )


def windows_for_admissions(
    admissions: Iterable[VertexAdmission], root: Path, prefix: str, window_vertices: int
) -> list[Window]:
    windows = []
    for source in admissions:
        for offset in range(0, source.records, window_vertices):
            count = min(window_vertices, source.records - offset)
            windows.append(
                Window(
                    source=source,
                    offset=offset,
                    count=count,
                    paths=window_paths(root, prefix, source, offset, count),
                )
            )
    return windows


def exact_generation_state(paths: WindowPaths) -> bool:
    present = [paths.shard.exists(), paths.generation_receipt.exists()]
    if any(present) and not all(present):
        raise RuntimeError(
            "partial generation evidence; refusing to continue: "
            + ", ".join(
                f"{path}={'present' if exists else 'missing'}"
                for path, exists in zip(
                    (paths.shard, paths.generation_receipt), present
                )
            )
        )
    return all(present)


def validate_summary(summary: dict[str, Any], window: Window) -> None:
    if summary.get("source_record_offset") != window.offset:
        raise RuntimeError("adjacency source offset mismatch")
    if summary.get("source_record_count") != window.count:
        raise RuntimeError("adjacency source count mismatch")
    if summary.get("directed_edges") != window.count * DIRECTED_EDGES_PER_VERTEX:
        raise RuntimeError("adjacency directed-edge count mismatch")
    if summary.get("kernel_coefficients") != window.count * KERNEL_COEFFICIENTS_PER_VERTEX:
        raise RuntimeError("adjacency kernel-coefficient count mismatch")
    if summary.get("generation_failures") != 0:
        raise RuntimeError("adjacency receipt reports generation failures")
    expected = [(5, window.count), (13, 2 * window.count), (37, 2 * window.count), (59, 2 * window.count)]
    recorded = [
        (item.get("degree"), item.get("directed_edges"))
        for item in summary.get("degree_counts", [])
        if isinstance(item, dict)
    ]
    if recorded != expected:
        raise RuntimeError("adjacency per-degree counts mismatch")
    require_hex(summary.get("record_chain_sha256"), 64, "record-chain digest")


def validate_generation(window: Window, producer_commit: str) -> dict[str, Any]:
    receipt = load_json(window.paths.generation_receipt)
    if receipt.get("status") != "pass":
        raise RuntimeError("adjacency generation status is not pass")
    if receipt.get("source_segment_digest", {}).get("sha256") != window.source.segment_sha256:
        raise RuntimeError("generation receipt source segment mismatch")
    if (
        receipt.get("source_verification_receipt_digest", {}).get("sha256")
        != window.source.verification_receipt_sha256
    ):
        raise RuntimeError("generation receipt source verification mismatch")
    if receipt.get("adjacency_shard_digest", {}).get("sha256") != sha256_file(window.paths.shard):
        raise RuntimeError("generation receipt adjacency digest mismatch")
    inner = receipt.get("receipt", {})
    if inner.get("producer_source_commit") != producer_commit:
        raise RuntimeError("generation source commit mismatch")
    validate_summary(inner.get("summary", {}), window)
    return receipt


def validate_verification(
    window: Window, producer_commit: str, verifier_commit: str
) -> dict[str, Any]:
    receipt = load_json(window.paths.verification_receipt)
    if receipt.get("status") != "pass":
        raise RuntimeError("adjacency verification status is not pass")
    if receipt.get("source_segment_digest", {}).get("sha256") != window.source.segment_sha256:
        raise RuntimeError("verification receipt source segment mismatch")
    if (
        receipt.get("source_verification_receipt_digest", {}).get("sha256")
        != window.source.verification_receipt_sha256
    ):
        raise RuntimeError("verification receipt source verification mismatch")
    if receipt.get("adjacency_shard_digest", {}).get("sha256") != sha256_file(window.paths.shard):
        raise RuntimeError("verification receipt adjacency digest mismatch")
    inner = receipt.get("receipt", {})
    if inner.get("producer_source_commit") != producer_commit:
        raise RuntimeError("verification producer source commit mismatch")
    if inner.get("verifier_source_commit") != verifier_commit:
        raise RuntimeError("verification source commit mismatch")
    validate_summary(inner.get("summary", {}), window)
    return receipt


def validate_existing_admission(window: Window, producer_sha: str, verifier_sha: str) -> None:
    value = load_json(window.paths.admission_manifest)
    if value.get("schema") != "ec-hardness/p61-adjacency-window-admission/v1":
        raise RuntimeError("unsupported adjacency admission schema")
    if value.get("status") != "pass":
        raise RuntimeError("adjacency admission status is not pass")
    if value.get("producer", {}).get("sha256") != producer_sha:
        raise RuntimeError("adjacency admission producer binary mismatch")
    if value.get("verifier", {}).get("sha256") != verifier_sha:
        raise RuntimeError("adjacency admission verifier binary mismatch")
    if value.get("shard", {}).get("sha256") != sha256_file(window.paths.shard):
        raise RuntimeError("adjacency admission shard digest mismatch")
    if value.get("generation_receipt", {}).get("sha256") != sha256_file(
        window.paths.generation_receipt
    ):
        raise RuntimeError("adjacency admission generation receipt mismatch")
    if value.get("verification_receipt", {}).get("sha256") != sha256_file(
        window.paths.verification_receipt
    ):
        raise RuntimeError("adjacency admission verification receipt mismatch")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--producer", type=Path, required=True)
    result.add_argument("--verifier", type=Path, required=True)
    result.add_argument("--vertex-manifest-dir", type=Path, required=True)
    result.add_argument("--vertex-manifest-pattern", default="campaign-v1-admission-*.json")
    result.add_argument("--root", type=Path, required=True)
    result.add_argument("--mount-root", type=Path, required=True)
    result.add_argument("--producer-source-commit", required=True)
    result.add_argument("--verifier-source-commit", required=True)
    result.add_argument("--class-number-receipt-sha256", required=True)
    result.add_argument("--window-vertices", type=int, default=65_536)
    result.add_argument("--threads", type=int, default=4)
    result.add_argument("--prefix", default="adjacency-v1")
    result.add_argument("--minimum-x", type=int, default=0)
    result.add_argument("--maximum-x", type=int, default=EXPECTED_FINAL_X)
    result.add_argument("--max-windows", type=int)
    result.add_argument("--min-free-bytes", type=int, default=1_500_000_000_000)
    result.add_argument("--max-logical-bytes", type=int, default=6_000_000_000_000)
    result.add_argument("--follow", action="store_true")
    result.add_argument("--poll-seconds", type=int, default=60)
    result.add_argument("--plan-only", action="store_true")
    return result


def validate_args(args: argparse.Namespace) -> None:
    for value, label in (
        (args.producer_source_commit, "producer source commit"),
        (args.verifier_source_commit, "verifier source commit"),
    ):
        require_hex(value, 40, label)
    require_hex(args.class_number_receipt_sha256, 64, "class-number receipt digest")
    if args.producer_source_commit == args.verifier_source_commit:
        raise RuntimeError("producer and verifier source commits must differ")
    if args.window_vertices <= 0 or args.threads <= 0:
        raise RuntimeError("window-vertices and threads must be positive")
    if args.max_windows is not None and args.max_windows <= 0:
        raise RuntimeError("max-windows must be positive")
    if args.follow and args.max_windows is not None:
        raise RuntimeError("--follow cannot be combined with --max-windows")
    if not (0 <= args.minimum_x < args.maximum_x <= EXPECTED_FINAL_X):
        raise RuntimeError("require 0 <= minimum-x < maximum-x <= 19894128")
    if args.poll_seconds <= 0 or args.min_free_bytes < 0 or args.max_logical_bytes <= 0:
        raise RuntimeError("poll-seconds and max-logical-bytes must be positive")


def plan(args: argparse.Namespace, verify_files: bool) -> list[Window]:
    admissions = discover_vertex_admissions(
        args.vertex_manifest_dir,
        args.vertex_manifest_pattern,
        args.minimum_x,
        args.maximum_x,
        verify_files,
    )
    windows = windows_for_admissions(admissions, args.root, args.prefix, args.window_vertices)
    if args.max_windows is not None:
        windows = windows[: args.max_windows]
    return windows


def plan_record(args: argparse.Namespace, windows: list[Window]) -> dict[str, Any]:
    return {
        "status": "plan",
        "window_vertices": args.window_vertices,
        "source_admissions": len({window.source.manifest for window in windows}),
        "windows": [
            {
                "source_manifest": str(window.source.manifest),
                "source_x_start": window.source.x_start,
                "source_x_end": window.source.x_end,
                "source_record_offset": window.offset,
                "source_record_count": window.count,
                "directed_edges": window.count * DIRECTED_EDGES_PER_VERTEX,
                "paths": {
                    name: str(path) for name, path in window.paths.__dict__.items()
                },
            }
            for window in windows
        ],
    }


def execute_window(
    args: argparse.Namespace,
    window: Window,
    producer_sha: str,
    verifier_sha: str,
    logical_bytes_so_far: int,
) -> int:
    require_mount(args.mount_root, "adjacency generation")
    free_before = require_free_space(args.mount_root, args.min_free_bytes, "adjacency generation")
    if window.paths.admission_manifest.exists():
        validate_generation(window, args.producer_source_commit)
        verification = validate_verification(
            window, args.producer_source_commit, args.verifier_source_commit
        )
        validate_existing_admission(window, producer_sha, verifier_sha)
        return int(verification["receipt"]["uncompressed_bytes"])

    generated = exact_generation_state(window.paths)
    if not generated:
        run(
            [
                str(args.producer),
                "--threads",
                str(args.threads),
                "generate-shard",
                "--source-segment",
                str(window.source.segment),
                "--source-verification-receipt",
                str(window.source.verification_receipt),
                "--source-record-offset",
                str(window.offset),
                "--source-record-count",
                str(window.count),
                "--batch-vertices",
                str(min(window.count, 64)),
                "--producer-source-commit",
                args.producer_source_commit,
                "--class-number-receipt-sha256",
                args.class_number_receipt_sha256,
                "--output",
                str(window.paths.shard),
                "--receipt-output",
                str(window.paths.generation_receipt),
            ]
        )
    generation = validate_generation(window, args.producer_source_commit)
    logical_bytes = int(generation["receipt"]["uncompressed_bytes"])
    if logical_bytes_so_far + logical_bytes > args.max_logical_bytes:
        raise RuntimeError(
            "logical corpus cap would be crossed before verification; generated shard is retained but not admitted"
        )
    require_mount(args.mount_root, "adjacency verification")
    require_free_space(args.mount_root, args.min_free_bytes, "adjacency verification")
    if not window.paths.verification_receipt.exists():
        run(
            [
                str(args.verifier),
                "--threads",
                str(args.threads),
                "verify-shard",
                "--source-segment",
                str(window.source.segment),
                "--source-verification-receipt",
                str(window.source.verification_receipt),
                "--input",
                str(window.paths.shard),
                "--verifier-source-commit",
                args.verifier_source_commit,
                "--receipt-output",
                str(window.paths.verification_receipt),
            ]
        )
    verification = validate_verification(
        window, args.producer_source_commit, args.verifier_source_commit
    )
    if int(verification["receipt"]["uncompressed_bytes"]) != logical_bytes:
        raise RuntimeError("producer and verifier uncompressed byte counts differ")
    require_mount(args.mount_root, "adjacency admission")
    os.sync()
    free_after = require_free_space(args.mount_root, args.min_free_bytes, "adjacency admission")
    manifest = {
        "schema": "ec-hardness/p61-adjacency-window-admission/v1",
        "status": "pass",
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "source": {
            "vertex_admission_manifest": str(window.source.manifest),
            "vertex_admission_sha256": sha256_file(window.source.manifest),
            "segment": str(window.source.segment),
            "segment_sha256": window.source.segment_sha256,
            "verification_receipt": str(window.source.verification_receipt),
            "verification_receipt_sha256": window.source.verification_receipt_sha256,
            "x_start": window.source.x_start,
            "x_end": window.source.x_end,
            "record_offset": window.offset,
            "record_count": window.count,
        },
        "producer": {
            "path": str(args.producer),
            "source_commit": args.producer_source_commit,
            "sha256": producer_sha,
        },
        "verifier": {
            "path": str(args.verifier),
            "source_commit": args.verifier_source_commit,
            "sha256": verifier_sha,
        },
        "shard": {
            "path": str(window.paths.shard),
            "sha256": sha256_file(window.paths.shard),
            "stored_bytes": window.paths.shard.stat().st_size,
            "logical_uncompressed_bytes": logical_bytes,
        },
        "generation_receipt": {
            "path": str(window.paths.generation_receipt),
            "sha256": sha256_file(window.paths.generation_receipt),
            "elapsed_ms": generation["receipt"]["elapsed_ms"],
        },
        "verification_receipt": {
            "path": str(window.paths.verification_receipt),
            "sha256": sha256_file(window.paths.verification_receipt),
            "elapsed_ms": verification["receipt"]["elapsed_ms"],
        },
        "counts": {
            "vertices": window.count,
            "directed_edges": window.count * DIRECTED_EDGES_PER_VERTEX,
            "kernel_coefficients": window.count * KERNEL_COEFFICIENTS_PER_VERTEX,
        },
        "storage_guard": {
            "minimum_free_bytes": args.min_free_bytes,
            "free_bytes_before": free_before,
            "free_bytes_after": free_after,
            "maximum_logical_bytes": args.max_logical_bytes,
        },
        "durability": {
            "per_file_sync": "receipt-recorded",
            "host_sync_before_manifest": "completed",
        },
    }
    write_json_new(window.paths.admission_manifest, manifest)
    manifest_sha256 = sha256_file(window.paths.admission_manifest)
    print(
        json.dumps(
            {
                "event": "admitted",
                "source_x_start": window.source.x_start,
                "source_x_end": window.source.x_end,
                "source_record_offset": window.offset,
                "source_record_count": window.count,
                "manifest": str(window.paths.admission_manifest),
                "manifest_sha256": manifest_sha256,
            }
        ),
        flush=True,
    )
    return logical_bytes


def main() -> int:
    args = parser().parse_args()
    validate_args(args)
    if args.plan_only:
        print(json.dumps(plan_record(args, plan(args, False)), indent=2, sort_keys=True))
        return 0
    require_mount(args.mount_root, "startup")
    require_free_space(args.mount_root, args.min_free_bytes, "startup")
    for binary in (args.producer, args.verifier):
        if not binary.is_file() or not os.access(binary, os.X_OK):
            raise RuntimeError(f"missing executable: {binary}")
    producer_sha = sha256_file(args.producer)
    verifier_sha = sha256_file(args.verifier)
    if producer_sha == verifier_sha:
        raise RuntimeError("producer and verifier binaries must differ")

    completed: set[Path] = set()
    logical_bytes = 0
    while True:
        windows = plan(args, True)
        made_progress = False
        for window in windows:
            if window.paths.admission_manifest in completed:
                continue
            window_bytes = execute_window(
                args, window, producer_sha, verifier_sha, logical_bytes
            )
            logical_bytes += window_bytes
            completed.add(window.paths.admission_manifest)
            made_progress = True
        if not args.follow:
            break
        max_seen_x = max((window.source.x_end for window in windows), default=0)
        if max_seen_x >= args.maximum_x and windows and all(
            window.paths.admission_manifest.exists() for window in windows
        ):
            break
        if not made_progress:
            print(
                json.dumps(
                    {
                        "event": "waiting_for_vertex_admission",
                        "max_seen_x": max_seen_x,
                        "poll_seconds": args.poll_seconds,
                    }
                ),
                flush=True,
            )
        time.sleep(args.poll_seconds)
    print(
        json.dumps(
            {
                "status": "pass",
                "admitted_windows_this_process": len(completed),
                "logical_uncompressed_bytes": logical_bytes,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

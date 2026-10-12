#!/usr/bin/env python3
"""Publish a counts-only status snapshot for a separate DP-weight walk.

The production ECC2K-130 feed and the table-walk feed deliberately have
different namespaces.  This module reads only S3 object metadata, checkpoint
headers, and the public portions of slot manifests; it never copies a point,
seed, or walk coefficient into the published JSON.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
import re
import struct
import sys
from collections import defaultdict
from typing import Any


UTC = dt.timezone.utc
SLOT_RE = re.compile(r"(?:^|/)slot-(\d+)\.json$")
CKPT_RE = re.compile(r"(?:^|/)slot-(\d+)\.ck$")
DP_SUFFIX = ".bin"
CHECKPOINT_HEADER_BYTES = 40


def utc_stamp(value: dt.datetime | None) -> str | None:
    if value is None:
        return None
    return value.astimezone(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def as_datetime(value: Any) -> dt.datetime | None:
    if isinstance(value, dt.datetime):
        return value.astimezone(UTC)
    if isinstance(value, (int, float)) and math.isfinite(value):
        # Slot manifests use Unix seconds; accept milliseconds as well.
        seconds = float(value) / 1000 if value > 10_000_000_000 else float(value)
        return dt.datetime.fromtimestamp(seconds, UTC)
    if isinstance(value, str):
        text = value.strip()
        if text.endswith("Z"):
            text = text[:-1] + "+00:00"
        try:
            parsed = dt.datetime.fromisoformat(text)
        except ValueError:
            return None
        return parsed.replace(tzinfo=parsed.tzinfo or UTC).astimezone(UTC)
    return None


def list_objects(s3: Any, bucket: str, prefix: str) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    paginator = s3.get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        result.extend(page.get("Contents", []))
    return result


def read_json(s3: Any, bucket: str, key: str) -> dict[str, Any] | None:
    try:
        body = s3.get_object(Bucket=bucket, Key=key)["Body"].read()
        value = json.loads(body.decode("utf-8"))
    except (KeyError, OSError, ValueError, UnicodeDecodeError):
        return None
    return value if isinstance(value, dict) else None


def object_time(item: dict[str, Any]) -> dt.datetime | None:
    return as_datetime(item.get("LastModified"))


def object_records(item: dict[str, Any]) -> int:
    size = int(item.get("Size", 0) or 0)
    # A distinguished-point record is 32 bytes in the public uploader format.
    # Ignore a partial object rather than inventing a record.
    return size // 32 if size >= 32 else 0


def floor_hour(value: dt.datetime) -> dt.datetime:
    return value.astimezone(UTC).replace(minute=0, second=0, microsecond=0)


def checkpoint_iterations(s3: Any, bucket: str, key: str) -> int | None:
    try:
        body = s3.get_object(
            Bucket=bucket,
            Key=key,
            Range=f"bytes=0-{CHECKPOINT_HEADER_BYTES - 1}",
        )["Body"].read()
    except (KeyError, OSError):
        return None
    if len(body) < CHECKPOINT_HEADER_BYTES:
        return None
    # state.ck is: magic/version/threads/batch/lanes, then the uint64 base
    # iteration counter.  Only the counter is exposed in the result.
    try:
        _version, _curve, threads, batch, lanes, _run_id = struct.unpack_from("<6I", body, 8)
        base = struct.unpack_from("<Q", body, 32)[0]
    except struct.error:
        return None
    if not (1 <= threads <= 10_000_000 and 1 <= batch <= 1_000_000_000 and 1 <= lanes <= 4096):
        return None
    return int(base) + int(threads) * int(batch) * int(lanes)


def slot_number(key: str, pattern: re.Pattern[str]) -> str | None:
    match = pattern.search(key)
    return match.group(1) if match else None


def build_status(s3: Any, bucket: str, prefix: str, now: dt.datetime | None = None) -> dict[str, Any]:
    now = (now or dt.datetime.now(UTC)).astimezone(UTC)
    prefix = prefix.rstrip("/") + "/"
    campaign = read_json(s3, bucket, prefix + "campaign.json") or {}
    dp_objects = [item for item in list_objects(s3, bucket, prefix + "dp/") if str(item.get("Key", "")).endswith(DP_SUFFIX)]
    slot_objects = list_objects(s3, bucket, prefix + "slots/")
    checkpoint_objects = list_objects(s3, bucket, prefix + "ckpt/")

    records = sum(object_records(item) for item in dp_objects)
    times = [stamp for item in dp_objects if (stamp := object_time(item)) is not None]
    first_dp = min(times) if times else None
    last_dp = max(times) if times else None
    hour_start = floor_hour(now - dt.timedelta(hours=47))
    hourly: dict[dt.datetime, int] = defaultdict(int)
    last_hour = now - dt.timedelta(hours=1)
    last_day = now - dt.timedelta(days=1)
    dps_last_hour = 0
    dps_last_day = 0
    for item in dp_objects:
        stamp = object_time(item)
        if stamp is None:
            continue
        count = object_records(item)
        if stamp >= hour_start:
            hourly[floor_hour(stamp)] += count
        if stamp >= last_hour:
            dps_last_hour += count
        if stamp >= last_day:
            dps_last_day += count

    latest_checkpoint: dict[str, tuple[dt.datetime, str]] = {}
    for item in checkpoint_objects:
        key = str(item.get("Key", ""))
        number = slot_number(key, CKPT_RE)
        stamp = object_time(item)
        if number is None or stamp is None:
            continue
        previous = latest_checkpoint.get(number)
        if previous is None or stamp > previous[0]:
            latest_checkpoint[number] = (stamp, key)

    slot_rows: list[dict[str, Any]] = []
    checkpoint_total = 0
    active_slots = 0
    for item in slot_objects:
        key = str(item.get("Key", ""))
        number = slot_number(key, SLOT_RE)
        if number is None:
            continue
        manifest = read_json(s3, bucket, key) or {}
        updated = as_datetime(manifest.get("updatedAt", manifest.get("updated_at"))) or object_time(item)
        checkpoint = latest_checkpoint.get(number)
        iterations = checkpoint_iterations(s3, bucket, checkpoint[1]) if checkpoint else None
        if iterations is not None:
            checkpoint_total += iterations
        age = (now - updated).total_seconds() if updated else None
        active = age is not None and 0 <= age <= 1800
        if active:
            active_slots += 1
        dp_uploaded = int(manifest.get("dpUploaded", manifest.get("dp_uploaded", 0)) or 0)
        row: dict[str, Any] = {
            "worker_id": "slot-" + number,
            "dps": dp_uploaded,
            "first_dp_at": as_datetime(manifest.get("createdAt", manifest.get("created_at"))),
            "last_dp_at": updated,
            "state": "walking" if active else ("retired" if manifest.get("retired") else "stale"),
        }
        # GPU identity is useful to compare heterogeneous workers, but keep
        # the allowlist narrow: no pod IDs, owner names, or configuration.
        gpu = manifest.get("gpuName", manifest.get("gpu_name"))
        if isinstance(gpu, str) and gpu.strip():
            row["gpu"] = gpu.strip()[:120]
        slot_rows.append({key: utc_stamp(value) if isinstance(value, dt.datetime) else value for key, value in row.items()})

    slot_rows.sort(key=lambda row: row["worker_id"])
    generated = now
    recent = last_dp is not None and (now - last_dp).total_seconds() <= 3600
    campaign_id = str(campaign.get("campaignId") or campaign.get("campaign_id") or prefix.rstrip("/").split("/")[-1])
    dp_weight = int(campaign.get("dpWeight", campaign.get("dp_weight", 34)) or 34)
    curve_id = str(campaign.get("curveId") or campaign.get("curve_id") or "certicom-ecc2k-130")
    campaign_created = as_datetime(campaign.get("createdAt", campaign.get("created_at")))
    points = [{"hour": utc_stamp(hour), "dps": count} for hour, count in sorted(hourly.items())]
    status: dict[str, Any] = {
        "schema_version": 1,
        "generated_at": utc_stamp(generated),
        "source": "S3 object metadata and sanitized table-walk checkpoints",
        "campaign_id": campaign_id,
        "curve_id": curve_id,
        "dp_mask_bits": dp_weight,
        "campaign_created_at": utc_stamp(campaign_created),
        "first_dp_at": utc_stamp(first_dp),
        "last_dp_at": utc_stamp(last_dp),
        "dps": records,
        "dps_last_hour": dps_last_hour,
        "dps_last_day": dps_last_day,
        "collisions": None,
        "latest_collision_at": None,
        "hourly": points,
        "per_worker": slot_rows,
        "work": {
            "iterations": checkpoint_total or None,
            "walking_slots": active_slots,
            "slots": len(slot_rows),
            "campaign_dp_weight": dp_weight,
        },
        "walk_rate": None,
        "ingest": {
            "newest_object_at": utc_stamp(last_dp),
            "lag_seconds": max(0, int((now - last_dp).total_seconds())) if last_dp else None,
            "semantics": "uploaded DP records from S3 object sizes; not de-duplicated or collision-checked",
        },
        "count_semantics": "uploaded DP records from the separate DP34 table-walk namespace",
        "claim_boundary": (
            "This DP34 table-walk stream is a separate research campaign. Its uploaded records, "
            "iteration totals, and rate are not merged into the DP32 ECC2K-130 production campaign."
        ),
    }
    # JSON has no datetime values; the comprehensions above keep the public
    # document deterministic and make accidental secret/raw-object leakage
    # impossible through this result structure.
    status["state"] = "COLLECTING" if recent else ("IDLE_OR_STALE" if records else "EMPTY")
    return status


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bucket", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    try:
        import boto3
    except ImportError as exc:  # pragma: no cover - exercised in CI's install step
        raise SystemExit("boto3 is required to read the DP34 S3 telemetry") from exc
    session = boto3.session.Session()
    s3 = session.client("s3")
    status = build_status(s3, args.bucket, args.prefix)
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(status, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print("DP34 snapshot: %d records, %d slots" % (status["dps"], status["work"]["slots"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

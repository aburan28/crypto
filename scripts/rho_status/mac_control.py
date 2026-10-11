#!/usr/bin/env python3
"""Publish an allowlisted snapshot of the separate Mac control campaign.

The source objects are private. Never copy a raw worker manifest into Pages:
it includes storage keys and may grow private walk fields later.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import subprocess
import tempfile
from pathlib import Path


CAMPAIGN_ID = "ecc2k130-synthetic-metal-h128-v1"
SOURCE_SCHEMA = "ecc2k130-synthetic-metal-continuous-v1"
PUBLIC_SCHEMA = "ecc2k130-mac-control-public-v1"
BUCKET = "ecc2k130-status-590183823895"
PREFIX = f"campaigns/{CAMPAIGN_ID}/"
WORKER_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,63}\Z")
MAX_WORKERS = 64
PUBLIC_KEY = "mac-control.json"
PUBLIC_FIELDS = {"schema", "campaign_id", "generated_at", "source_state", "read_errors",
                 "worker_count", "dp_records", "walk_updates", "workers"}
PUBLIC_WORKER_FIELDS = {"worker_id", "updated_at", "sequence", "dp_records", "walk_updates"}


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def nonnegative_int(value: object) -> int:
    if type(value) is not int or value < 0:
        raise ValueError("invalid nonnegative counter")
    return value


def parse_utc(value: object) -> str:
    if not isinstance(value, str):
        raise ValueError("missing worker timestamp")
    try:
        stamp = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise ValueError("invalid worker timestamp") from exc
    if stamp.tzinfo is None:
        raise ValueError("worker timestamp has no timezone")
    return stamp.astimezone(dt.timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def read_worker(worker_id: str, manifest: object) -> tuple[dict, str, str]:
    """Return public fields and private identity hashes used only for grouping."""
    if not WORKER_ID.fullmatch(worker_id) or not isinstance(manifest, dict):
        raise ValueError("invalid worker record")
    if manifest.get("schema") != SOURCE_SCHEMA:
        raise ValueError("unexpected worker schema")
    cumulative = manifest.get("cumulative")
    if not isinstance(cumulative, dict):
        raise ValueError("missing cumulative counters")
    run_id = manifest.get("runIdentity")
    walk_id = manifest.get("walkIdentity")
    if not all(isinstance(x, str) and re.fullmatch(r"[0-9a-f]{64}", x) for x in (run_id, walk_id)):
        raise ValueError("invalid worker identity")
    row = {
        "worker_id": worker_id,
        "updated_at": parse_utc(manifest.get("updatedAt")),
        "sequence": nonnegative_int(manifest.get("sequence")),
        "dp_records": nonnegative_int(cumulative.get("dpRecords")),
        "walk_updates": nonnegative_int(cumulative.get("walkUpdates")),
    }
    return row, run_id, walk_id


def snapshot(records: list[tuple[str, object]], *, generated_at: str, errors: int = 0,
             source_available: bool = True) -> dict:
    """Build public data from verified worker manifests; omit all other fields."""
    workers = []
    run_ids = []
    walk_ids = []
    for worker_id, manifest in records:
        try:
            row, run_id, walk_id = read_worker(worker_id, manifest)
        except ValueError:
            errors += 1
            continue
        workers.append(row)
        run_ids.append(run_id)
        walk_ids.append(walk_id)
    workers.sort(key=lambda row: row["worker_id"])
    state = "available" if workers else "empty"
    if not source_available:
        state = "unavailable"
    elif len(set(walk_ids)) > 1:
        state = "mixed_walks"
    elif len(set(run_ids)) != len(run_ids):
        state = "ambiguous_runs"
    elif errors:
        state = "partial"
    can_sum = bool(workers) and state in ("available", "partial")
    return {
        "schema": PUBLIC_SCHEMA,
        "campaign_id": CAMPAIGN_ID,
        "generated_at": generated_at,
        "source_state": state,
        "read_errors": errors,
        "worker_count": len(workers),
        "dp_records": sum(row["dp_records"] for row in workers) if can_sum else None,
        "walk_updates": sum(row["walk_updates"] for row in workers) if can_sum else None,
        "workers": workers,
    }


def aws_bytes(*args: str) -> bytes:
    result = subprocess.run(("aws", *args), capture_output=True, timeout=25, check=False)
    if result.returncode:
        raise RuntimeError(f"AWS read failed (exit {result.returncode})")
    return result.stdout


def from_s3() -> dict:
    now = utc_now()
    try:
        listing = json.loads(aws_bytes("s3api", "list-objects-v2", "--bucket", BUCKET,
                                      "--prefix", PREFIX, "--delimiter", "/", "--output", "json"))
        prefixes = listing.get("CommonPrefixes", [])
        if not isinstance(prefixes, list) or len(prefixes) > MAX_WORKERS:
            raise ValueError("unexpected worker listing")
    except (OSError, RuntimeError, ValueError, json.JSONDecodeError, subprocess.TimeoutExpired):
        return snapshot([], generated_at=now, source_available=False)
    records = []
    errors = 0
    for item in prefixes:
        path = item.get("Prefix", "") if isinstance(item, dict) else ""
        worker_id = path[len(PREFIX):-1] if path.startswith(PREFIX) and path.endswith("/") else ""
        if not WORKER_ID.fullmatch(worker_id):
            errors += 1
            continue
        try:
            raw = aws_bytes("s3", "cp", f"s3://{BUCKET}/{PREFIX}{worker_id}/latest.json", "-",
                            "--only-show-errors")
            records.append((worker_id, json.loads(raw)))
        except (OSError, RuntimeError, ValueError, json.JSONDecodeError, subprocess.TimeoutExpired):
            errors += 1
    return snapshot(records, generated_at=now, errors=errors)


def publish_live(result: dict) -> None:
    """Publish only the validated public schema at one fixed, public S3 key."""
    if (set(result) != PUBLIC_FIELDS or result.get("schema") != PUBLIC_SCHEMA or
            result.get("campaign_id") != CAMPAIGN_ID or
            not isinstance(result.get("workers"), list) or
            any(not isinstance(row, dict) or set(row) != PUBLIC_WORKER_FIELDS
                for row in result["workers"])):
        raise ValueError("refusing to publish a non-public Mac control document")
    encoded = json.dumps(result, sort_keys=True, separators=(",", ":")).encode("utf-8") + b"\n"
    with tempfile.NamedTemporaryFile(prefix="mac-control-public-", suffix=".json") as output:
        output.write(encoded)
        output.flush()
        command = ("aws", "s3api", "put-object", "--bucket", BUCKET, "--key", PUBLIC_KEY,
                   "--body", output.name, "--content-type", "application/json",
                   "--cache-control", "public, max-age=60")
        completed = subprocess.run(command, capture_output=True, timeout=45, check=False)
        if completed.returncode:
            raise RuntimeError(f"public Mac status upload failed (exit {completed.returncode})")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--publish-live", action="store_true",
                        help="also upload the sanitized document to the public status bucket")
    args = parser.parse_args()
    result = from_s3()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n", encoding="utf-8")
    if args.publish_live:
        publish_live(result)
    print(f"Mac control snapshot: {result['source_state']}, {result['worker_count']} workers, "
          f"{result['dp_records']} DP records, {result['read_errors']} read errors")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

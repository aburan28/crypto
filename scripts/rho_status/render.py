#!/usr/bin/env python3
"""Merge hourly snapshots into history.json for the public status page.

The HTML/CSS in docs/ecc2k130-status/ is static. This script only maintains
the JSON the page fetches.
"""

from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone

HISTORY_LIMIT = 168


def load_json(path, default):
    if not path or not os.path.exists(path):
        return default
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def parse_time(value):
    if not value:
        return None
    text = str(value).replace("Z", "+00:00")
    return datetime.fromisoformat(text)


def merge_history(previous, snapshot, limit=HISTORY_LIMIT):
    points = []
    if isinstance(previous, dict):
        points = list(previous.get("points") or [])
    elif isinstance(previous, list):
        points = list(previous)
    points.append(
        {
            "generated_at": snapshot.get("generated_at"),
            "dps": int(snapshot.get("dps") or 0),
            "dps_last_hour": int(snapshot.get("dps_last_hour") or 0),
            "collisions": int(snapshot.get("collisions") or 0),
            "workers": int(snapshot.get("workers") or 0),
            "state": snapshot.get("state"),
        }
    )
    latest = {}
    order = []
    for point in points:
        key = point.get("generated_at")
        if key not in latest:
            order.append(key)
        latest[key] = point
    deduped = [latest[key] for key in order]

    def sort_key(item):
        stamp = parse_time(item.get("generated_at"))
        return stamp or datetime.fromtimestamp(0, timezone.utc)

    deduped.sort(key=sort_key)
    return {
        "schema_version": 1,
        "campaign_id": snapshot.get("campaign_id"),
        "points": deduped[-limit:],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--status", required=True)
    parser.add_argument("--history-in", default="")
    parser.add_argument("--history-out", required=True)
    args = parser.parse_args(argv)
    snapshot = load_json(args.status, None)
    if not isinstance(snapshot, dict):
        raise SystemExit("status JSON missing")
    history = merge_history(load_json(args.history_in, {}), snapshot)
    os.makedirs(os.path.dirname(os.path.abspath(args.history_out)) or ".", exist_ok=True)
    with open(args.history_out, "w", encoding="utf-8") as fh:
        json.dump(history, fh, indent=2, sort_keys=True)
        fh.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

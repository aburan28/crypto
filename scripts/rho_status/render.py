#!/usr/bin/env python3
"""Merge published snapshots into history.json for the public status page.

The HTML/CSS in docs/ecc2k130-status/ is static. This script only maintains
the JSON the page fetches.

It also measures the one number no single snapshot contains: the walk rate.
work_feed.py puts the campaign's checkpointed iteration total into each
snapshot, so the difference between two of them, over the time between them,
is iterations per second -- measured from the walk instead of derived from
the point count and a report interval the campaign no longer runs. The
history this script maintains is where those two snapshots come from, which
is why the measurement lives here and not in the pages.
"""

from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone

# Seven days of snapshots at the publish cadence in ecc2k130-status.yml
# (96 a day, every 15 minutes). This is a count, so it tracks the cron: at a
# slower cadence it covers more than a week, at a faster one less.
HISTORY_LIMIT = 672

# Smoothing window for the walk rate. A slot's iteration total only moves when
# it checkpoints, every 600 s, so a rate taken across a single 15-minute
# publish interval steps with which slots happened to land inside it. An hour
# covers every slot several times over and still reports the fleet that is
# running now rather than the campaign's average.
RATE_WINDOW_S = 3600
# Below this the checkpoint granularity dominates the difference.
MIN_RATE_SPAN_S = 600
# How far back to reach when the window holds no second sample -- after an
# outage in the publishing job, or on the first run after this field appears.
# The measured span is published next to the rate, so a long one is visible.
MAX_RATE_SPAN_S = 6 * 3600


def load_json(path, default):
    if not path or not os.path.exists(path):
        return default
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def parse_time(value):
    if not value:
        return None
    text = str(value).replace("Z", "+00:00")
    try:
        return datetime.fromisoformat(text)
    except ValueError:
        return None


def snapshot_iterations(snapshot):
    """The checkpointed iteration total in a snapshot, when it carries one.

    work_feed.py refuses to publish a stale or foreign total, so the field is
    absent rather than wrong when the campaign's feed is not answering. Points
    published before that field existed have none either.
    """
    work = snapshot.get("work")
    if not isinstance(work, dict):
        return None
    try:
        iterations = int(work.get("iterations"))
    except (TypeError, ValueError):
        return None
    return iterations if iterations > 0 else None


def rate_samples(points):
    """(time, iterations, timestamp as published) for points that carry both."""
    samples = []
    for point in points:
        stamp = point.get("generated_at")
        at = parse_time(stamp)
        try:
            iterations = int(point.get("iterations"))
        except (TypeError, ValueError):
            continue
        if at is not None and iterations > 0:
            samples.append((at, iterations, stamp))
    samples.sort(key=lambda sample: sample[0])
    return samples


def measure_rate(points, window_s=RATE_WINDOW_S, min_span_s=MIN_RATE_SPAN_S,
                 max_span_s=MAX_RATE_SPAN_S):
    """Iterations per second between two published snapshots, or None.

    The difference is taken against the oldest sample inside the window, so it
    spans about an hour of checkpoints rather than one publish interval; when
    an outage leaves the window with nothing to compare against it reaches
    back to the newest sample outside it, and publishes the span it used.

    None, not a number, whenever the difference would not mean anything: one
    sample, too short a span, or a total that went backwards. A slot whose
    checkpoint was reset contributes less than it did before, and no rate is
    more honest than a negative one.
    """
    samples = rate_samples(points)
    if len(samples) < 2:
        return None
    latest_at, latest, latest_stamp = samples[-1]
    inside = [s for s in samples[:-1] if (latest_at - s[0]).total_seconds() <= window_s]
    if inside:
        base_at, base, base_stamp = inside[0]
    else:
        outside = [s for s in samples[:-1] if (latest_at - s[0]).total_seconds() <= max_span_s]
        if not outside:
            return None
        base_at, base, base_stamp = outside[-1]
    span = (latest_at - base_at).total_seconds()
    if span < min_span_s or latest < base:
        return None
    return {
        "iterations_per_second": (latest - base) / span,
        "window_seconds": int(round(span)),
        "measured_from": base_stamp,
        "measured_to": latest_stamp,
        "iterations_from": base,
        "iterations_to": latest,
        "method": "difference of the campaign's checkpointed iteration totals "
                  "between two published snapshots",
    }


def stamp_published(snapshot, now=None):
    """Record when the publishing job wrote this document.

    `generated_at` is when the source counted -- the walker's query, or the
    ingest host's write, which may be up to its --status-every old when the
    hop is down. Without a second stamp a reader cannot tell a publisher that
    stopped from a source that did: Pages froze on an 11:25Z snapshot for 17
    hours on 2026-09-18 and the page had one sentence for both. The dashboard
    prints both and says which one is behind.
    """
    now = now or datetime.now(timezone.utc)
    snapshot["published_at"] = now.replace(microsecond=0).isoformat().replace("+00:00", "Z")
    return snapshot["published_at"]


def apply_rate(snapshot, rate):
    """Attach a measured rate to the snapshot, or clear a stale one. True when set."""
    if rate is None:
        snapshot.pop("walk_rate", None)
        return False
    snapshot["walk_rate"] = rate
    return True


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
            "iterations": snapshot_iterations(snapshot),
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
    parser.add_argument("--status-out", default="",
                        help="rewrite the snapshot with the measured walk rate")
    args = parser.parse_args(argv)
    snapshot = load_json(args.status, None)
    if not isinstance(snapshot, dict):
        raise SystemExit("status JSON missing")
    history = merge_history(load_json(args.history_in, {}), snapshot)
    os.makedirs(os.path.dirname(os.path.abspath(args.history_out)) or ".", exist_ok=True)
    with open(args.history_out, "w", encoding="utf-8") as fh:
        json.dump(history, fh, indent=2, sort_keys=True)
        fh.write("\n")
    # A historical rate is only valid when this snapshot carries the total it
    # was measured against; a missing feed must not publish one anyway.
    rate = (
        measure_rate(history["points"])
        if snapshot_iterations(snapshot) is not None else None
    )
    if args.status_out:
        apply_rate(snapshot, rate)
        stamp_published(snapshot)
        with open(args.status_out, "w", encoding="utf-8") as fh:
            json.dump(snapshot, fh, indent=2, sort_keys=True)
            fh.write("\n")
    if rate is None:
        print("no walk rate yet: needs two published snapshots with an iteration total")
    else:
        print("walk rate %.3f B it/s over %d s" % (
            rate["iterations_per_second"] / 1e9, rate["window_seconds"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

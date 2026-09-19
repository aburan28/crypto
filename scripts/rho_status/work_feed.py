#!/usr/bin/env python3
"""Copy the campaign's own iteration count into the public status snapshot.

snapshot.py reports what the distinguished-point store holds: points, and
when they arrived. The walk behind those points was derived on the published
pages from the point count times one report interval, and that interval is
not the one the campaign runs: the pages carry the weight-34 figure of
2^25.27 iterations per point (ecc2k130/aws/README.md) while campaign.json
runs `dpWeight` 32, a strictly sparser cutoff. Measured from each client's
own iteration and point counters the running interval is about 2^28.4
(ecc2k130/benchmarks/dp-interval), so the derivation understates the walk
by a factor near nine, and any rate taken from it is wrong by the same
factor.

Nothing has to be asked of the workers to fix that. Every slot checkpoints
its iteration base every 600 s, the campaign's ingest host sums those
checkpoints, and it publishes the sum in a feed whose only reader until now
was the operator. This script copies that total into the snapshot the pages
read, and render.py then measures the rate from the difference between two
published snapshots -- so the published walk rate comes from the walk.

The feed URL is not in the repository: the bucket is named for the AWS
account, which this public tree deliberately does not carry (the walker host
is a secret for the same reason). It is derived at run time from the calling
identity, the same `<stack>-status-<account>` convention infra.sh uses, or
given outright in RHO_WORK_FEED_URL.

Fail-soft on purpose when merging into a walker snapshot. A feed that is
missing, unreachable, stale, unparseable or from another campaign leaves
that snapshot exactly as snapshot.py wrote it and exits 0: the pages
already say they have no iteration count for that case, and losing the
walk total must not also lose the 15-minute point publication.

`--as-snapshot` is the other mode: the walker hop is down, so this file
*is* the publication. The ingest host already writes the same counts to
the status bucket every `--status-every` seconds. Republish that document
(without `per_slot`) rather than failing the job and leaving Pages on a
snapshot whose `generated_at` is a day old. Fail hard if the feed cannot
be read: there is nothing else to publish.

Usage:
    python3 scripts/rho_status/work_feed.py --status docs/ecc2k130-status/status.json
    python3 scripts/rho_status/work_feed.py --as-snapshot --status docs/ecc2k130-status/status.json
    RHO_WORK_FEED_URL=https://host/status.json python3 scripts/rho_status/work_feed.py --status s.json
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone

from snapshot import DEFAULT_CAMPAIGN, assert_public, iso, normalize, write_snapshot_file

# The feed is rewritten continuously by the ingest host, so an hour without a
# rewrite means that host stopped, not that the walk did. Publishing a frozen
# total would read as a campaign that halted, and would drag the measured rate
# towards zero as the snapshots caught up with it, so refuse it instead: no
# iteration count is a state the pages already render honestly.
MAX_FEED_AGE_S = 3600

# Slots checkpoint every 600 s (campaign.json `checkpointEvery`). Three missed
# checkpoints is a dead worker rather than a slow one.
FRESH_CHECKPOINT_S = 1800

DEFAULT_STACK = "ecc2k130"
DEFAULT_REGION = "us-west-2"
FETCH_TIMEOUT_S = 20


def parse_time(value):
    if not value:
        return None
    text = str(value).replace("Z", "+00:00")
    try:
        stamp = datetime.fromisoformat(text)
    except ValueError:
        return None
    return stamp if stamp.tzinfo else stamp.replace(tzinfo=timezone.utc)


def account_id():
    """The calling AWS account, or None when there is no usable credential."""
    try:
        proc = subprocess.run(
            ["aws", "sts", "get-caller-identity", "--query", "Account", "--output", "text"],
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    account = proc.stdout.strip()
    return account if proc.returncode == 0 and account.isdigit() else None


def feed_url(stack=DEFAULT_STACK, region=DEFAULT_REGION, account=None):
    account = account or account_id()
    if not account:
        return None
    return "https://%s-status-%s.s3.%s.amazonaws.com/status.json" % (stack, account, region)


def fetch(url, timeout=FETCH_TIMEOUT_S):
    with urllib.request.urlopen(url, timeout=timeout) as response:  # nosec B310 - https only
        return json.loads(response.read().decode("utf-8"))


def work_block(feed, campaign, now=None, max_age_s=MAX_FEED_AGE_S):
    """The publishable part of the feed, or None when it says nothing usable.

    Returns the iteration total the workers checkpointed, how many slots are
    behind it, and when the feed was written -- enough for a reader to see
    what the figure is and how current, and nothing per-worker beyond a count.

    `max_age_s` is how old a total may be when it is overlaid on a *fresh*
    walker snapshot. None disables that check: the feed is the snapshot, so
    a two-hour-old ingest document is still the current one, and the page
    reads freshness from `generated_at` rather than from us dropping the
    total.
    """
    if not isinstance(feed, dict) or feed.get("campaign_id") != campaign:
        return None
    work = feed.get("work")
    if not isinstance(work, dict):
        return None
    try:
        iterations = int(work.get("iterations"))
    except (TypeError, ValueError):
        return None
    if iterations <= 0:
        return None
    generated = parse_time(feed.get("generated_at"))
    if generated is None:
        return None
    now = now or datetime.now(timezone.utc)
    age = (now - generated).total_seconds()
    if max_age_s is not None and age > max_age_s:
        return None
    slots = work.get("per_slot") or []
    walking = [
        slot
        for slot in slots
        if isinstance(slot, dict)
        and not slot.get("retired")
        and isinstance(slot.get("checkpoint_age_s"), (int, float))
        and slot["checkpoint_age_s"] <= FRESH_CHECKPOINT_S
    ]
    block = {
        "iterations": iterations,
        "iterations_log2": round(math.log2(iterations), 6),
        "slots": len(slots),
        "walking_slots": len(walking),
        "feed_generated_at": feed.get("generated_at"),
        "feed_age_seconds": int(max(0.0, age)),
        "method": work.get("method") or "sum of the slots' checkpointed iteration bases",
        "source": feed.get("source") or "campaign work feed",
    }
    return block


def ingest_block(feed):
    """The publishable ingest-health fields, or None when the feed has none.

    `outstanding_objects` and `unrecognised_objects` are how the page tells a
    stopped fleet from a stopped ingest. On 2026-09-17 the store took nothing
    for five hours while the fleet walked, and the dashboard's only sentence
    for that was IDLE_OR_STALE — the same sentence as a dead walk. The feed
    carries the distinction; copy it, never invent it.
    """
    if not isinstance(feed, dict):
        return None
    ingest = feed.get("ingest")
    if not isinstance(ingest, dict):
        return None
    try:
        outstanding = int(ingest.get("outstanding_objects") or 0)
        unrecognised = int(ingest.get("unrecognised_objects") or 0)
    except (TypeError, ValueError):
        return None
    lag = ingest.get("lag_seconds")
    try:
        lag = int(lag) if lag is not None else None
    except (TypeError, ValueError):
        lag = None
    return {
        "outstanding_objects": outstanding,
        "unrecognised_objects": unrecognised,
        "newest_object_at": ingest.get("newest_object_at"),
        "lag_seconds": lag,
    }


def snapshot_from_feed(feed, campaign, now=None):
    """A public Pages snapshot built from the ingest host's status.json.

    The ingest document is already aggregates: campaign counts, a 48-hour
    hourly series, checkpointed work. It also carries `work.per_slot`, which
    names every walker's slot and run and must not be published. per_worker
    is empty: that table is DISTINCT worker_id from the DP heap, which this
    feed does not have. The GPU card still reads walking_slots from work.
    """
    if not isinstance(feed, dict) or feed.get("campaign_id") != campaign:
        return None
    generated = parse_time(feed.get("generated_at"))
    if generated is None:
        return None
    if feed.get("dps") is None:
        return None
    ingest = feed.get("ingest") if isinstance(feed.get("ingest"), dict) else {}
    snapshot = normalize(
        {
            "campaign_id": feed.get("campaign_id") or campaign,
            "curve_id": feed.get("curve_id"),
            "dp_mask_bits": feed.get("dp_mask_bits"),
            "campaign_created_at": feed.get("campaign_created_at"),
            "dps": feed.get("dps"),
            "workers": 0,
            "first_dp_at": feed.get("first_dp_at"),
            "last_dp_at": feed.get("last_dp_at"),
            "dps_last_hour": feed.get("dps_last_hour"),
            "dps_last_day": feed.get("dps_last_day"),
            "collisions": feed.get("collisions"),
            "latest_collision_at": feed.get("latest_collision_at"),
            "per_worker": [],
            "hourly": feed.get("hourly") or [],
            "ingest_outstanding": ingest.get("outstanding_objects"),
            "ingest_unrecognised": ingest.get("unrecognised_objects"),
        },
        campaign,
        "ingest-status-feed",
    )
    # normalize() stamps now(); the page's stale banner is about this field,
    # so it has to be when the ingest host counted, not when we copied it.
    snapshot["generated_at"] = iso(generated)
    snapshot["query_driver"] = "ingest-feed"
    merge_work(snapshot, feed, campaign, now=now, max_age_s=None)
    assert_public(snapshot)
    return snapshot


def merge_work(status, feed, campaign, now=None, max_age_s=MAX_FEED_AGE_S):
    """Attach the feed's iteration total to the snapshot. True when it did."""
    block = work_block(feed, campaign, now=now, max_age_s=max_age_s)
    if block is None:
        return False
    status["work"] = block
    ingest = ingest_block(feed)
    if ingest is not None:
        status["ingest"] = ingest
        if (
            status.get("state") == "IDLE_OR_STALE"
            and (ingest["outstanding_objects"] or ingest["unrecognised_objects"])
        ):
            # The store's last-hour count is about the ingest, not the walk.
            status["state"] = "INGEST_BEHIND"
    assert_public(status)
    return True


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--status", required=True, help="snapshot written by snapshot.py")
    parser.add_argument("--campaign", default=os.environ.get("RHO_CAMPAIGN", DEFAULT_CAMPAIGN))
    parser.add_argument("--feed-url", default=os.environ.get("RHO_WORK_FEED_URL", ""))
    parser.add_argument("--stack", default=os.environ.get("RHO_STACK", DEFAULT_STACK))
    parser.add_argument("--region", default=os.environ.get("AWS_DEFAULT_REGION", DEFAULT_REGION))
    parser.add_argument(
        "--as-snapshot",
        action="store_true",
        help="write a full public snapshot from the ingest feed; fail if it cannot",
    )
    args = parser.parse_args(argv)

    url = args.feed_url or feed_url(args.stack, args.region)
    if not url:
        print("no work feed URL and no AWS identity to derive one; leaving snapshot as it is")
        return 1 if args.as_snapshot else 0
    try:
        feed = fetch(url)
    except (urllib.error.URLError, OSError, ValueError, TimeoutError) as err:
        print("work feed unreadable (%s); leaving snapshot as it is" % err)
        return 1 if args.as_snapshot else 0

    if args.as_snapshot:
        snapshot = snapshot_from_feed(feed, args.campaign)
        if snapshot is None:
            print("ingest feed is not a usable %s snapshot" % args.campaign, file=sys.stderr)
            return 1
        write_snapshot_file(args.status, snapshot)
        generated = parse_time(snapshot.get("generated_at"))
        now = datetime.now(timezone.utc)
        age = int(max(0.0, (now - generated).total_seconds())) if generated else -1
        print(
            "wrote ingest-feed snapshot: dps=%d state=%s generated_at=%s age=%ds"
            % (snapshot["dps"], snapshot["state"], snapshot["generated_at"], age)
        )
        return 0

    with open(args.status, encoding="utf-8") as fh:
        status = json.load(fh)
    if not merge_work(status, feed, args.campaign):
        print("work feed carries no usable iteration count for %s; leaving snapshot as it is"
              % args.campaign)
        return 0
    with open(args.status, "w", encoding="utf-8") as fh:
        json.dump(status, fh, indent=2, sort_keys=True)
        fh.write("\n")
    work = status["work"]
    print("merged %d iterations (2^%.2f) from %d slots, feed %ds old"
          % (work["iterations"], work["iterations_log2"], work["walking_slots"],
             work["feed_age_seconds"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())

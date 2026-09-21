#!/usr/bin/env python3
"""Fail the run when the snapshot we just published is not current.

The publish job copies whatever the ingest feed last wrote. That copy
succeeding says nothing about whether the number moved: on 2026-09-20 the
ingest host stopped writing at 04:53Z and the workflow went on publishing
that same frozen snapshot every three minutes, green, for hours. 802
successful runs meant "I copied a file", not "the campaign is being
counted", and the two were indistinguishable from the outside.

This step is the end-to-end assertion that they are the same thing. It runs
*after* the Pages deploy on purpose. Everything else in this workflow fails
soft so that a view problem never freezes the dashboard -- the page renders
a stale banner honestly, and readers are better served by an old number
labelled old than by a page that stopped deploying. So publication happens
first and this only decides the colour of the run afterwards.

It is deliberately indifferent to *why* the feed stopped. A dead ingest
host, an expired credential, an RDS failover, a key shape nobody recognises
(which is how the store silently stopped for five hours on 2026-09-17): all
of them look like one thing here, a `generated_at` that is not advancing.
That is the property worth alarming on, because the list of causes is not
knowable in advance.

Threshold is `work_feed.MAX_FEED_AGE_S` (3600 s), the same constant
work_feed.py uses to refuse a frozen iteration total, so "too old to
publish as a total" and "too old to call the run green" cannot drift apart.
Note this is *tighter* than the dashboard's own STALE_AFTER_MS (90 min):
the run should go red slightly before the page starts apologising to
readers, not after.

Usage:
    python3 scripts/rho_status/check_feed_age.py --status docs/ecc2k130-status/status.json
    python3 scripts/rho_status/check_feed_age.py --status s.json --max-age-s 7200
    python3 scripts/rho_status/check_feed_age.py --status s.json --warn-only
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone

from work_feed import MAX_FEED_AGE_S, parse_time


def describe(status, now, max_age_s=MAX_FEED_AGE_S):
    """Return (ok, message) for a published snapshot.

    ok is True when the snapshot is current, or when there is legitimately
    nothing to be stale about yet. The message is written to be read in a
    workflow log by someone who has not seen this file.
    """
    if not isinstance(status, dict):
        return False, "snapshot is not a JSON object"

    generated = parse_time(status.get("generated_at"))

    # A campaign with no points has no feed to have stopped. The committed
    # placeholder is exactly this, so a fresh checkout does not read as an
    # outage. Guarded on dps too: a campaign that collected and then went
    # quiet must NOT be excused by this branch.
    dps = status.get("dps")
    if generated is None and not dps:
        return True, (
            "no snapshot yet (state=%s source=%s): nothing has been counted, "
            "so there is no feed age to check"
            % (status.get("state"), status.get("source"))
        )

    if generated is None:
        return False, (
            "snapshot carries %s points but no generated_at, so its age cannot "
            "be checked -- treating that as not current" % f"{dps:,}"
        )

    age = (now - generated).total_seconds()
    if age < 0:
        # A snapshot from the future is a clock problem, not freshness. Say so
        # rather than reporting a negative age as healthy.
        return False, (
            "snapshot generated_at is %.0f s in the future (%s); check the "
            "clock on whatever wrote it" % (-age, status.get("generated_at"))
        )

    published = parse_time(status.get("published_at"))
    detail = ""
    if published is not None:
        pub_age = (now - published).total_seconds()
        # Same split the dashboard draws: a fresh publish stamp on a stale
        # snapshot means the publisher ran and found nothing newer upstream,
        # so the ingest side is what stopped. An old stamp means the
        # publisher itself is not running.
        if pub_age <= max_age_s:
            detail = (
                " The publisher ran %.0f min ago and found nothing newer, so "
                "the ingest feed is what stopped, not this workflow." % (pub_age / 60)
            )
        else:
            detail = (
                " The publisher itself last wrote %.0f min ago, so this "
                "workflow is not running either." % (pub_age / 60)
            )

    if age > max_age_s:
        return False, (
            "feed is stale: generated_at %s is %.0f min old, over the %.0f min "
            "limit.%s" % (status.get("generated_at"), age / 60, max_age_s / 60, detail)
        )

    return True, (
        "feed is current: generated_at %s is %.0f min old, within the %.0f min limit"
        % (status.get("generated_at"), age / 60, max_age_s / 60)
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--status", required=True, help="the snapshot just published")
    parser.add_argument(
        "--max-age-s",
        type=float,
        default=float(os.environ.get("RHO_MAX_FEED_AGE_S", MAX_FEED_AGE_S)),
        help="age above which the run fails (default %d)" % MAX_FEED_AGE_S,
    )
    parser.add_argument(
        "--warn-only",
        action="store_true",
        help="report and exit 0; for rolling this out without failing runs first",
    )
    args = parser.parse_args(argv)

    try:
        with open(args.status) as handle:
            status = json.load(handle)
    except (OSError, ValueError) as exc:
        print("check_feed_age: cannot read %s: %s" % (args.status, exc), file=sys.stderr)
        return 1

    ok, message = describe(status, datetime.now(timezone.utc), args.max_age_s)
    stream = sys.stdout if ok else sys.stderr
    print("check_feed_age: %s" % message, file=stream)
    if ok or args.warn_only:
        if not ok:
            print("check_feed_age: --warn-only, not failing the run", file=sys.stderr)
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())

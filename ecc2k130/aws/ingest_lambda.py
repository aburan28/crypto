#!/usr/bin/env python3
"""Lambda entry point for dp_ingest: one scheduled pass, then exit.

The ingest loop never needed a host. `dp_ingest.py` already has the three
properties that make a scheduled invocation equivalent to a resident daemon,
and its own docstring is where they are argued:

  * idempotent -- every insert is ON CONFLICT DO NOTHING on
    (campaign_id, point_key), so a re-run or a kill mid-pass costs nothing;
  * no state outside the database it describes -- progress lives in
    `dp_ingest_progress`, committed in the same transaction as the points, so
    a fresh invocation resumes exactly where the last one stopped and there
    is no host-local file that can disagree;
  * bounded passes -- `--pass-objects` caps the work and `--once` publishes
    status.json before returning.

So this file holds no logic of its own. It builds an argv and calls
`dp_ingest.main()`, which is the same entry point `ingest.sh` runs. That is
deliberate and worth keeping: controlplane/README.md records that the
campaign's rules once lived in `worker.py`, in `dp_ingest.py` and in a shell
script at the same time, and that two of them could disagree without anybody
noticing until a GPU wrote into a run id somebody else was walking. A fourth
copy of "what a pass is" would be that mistake again.

Everything configurable is read from the environment by `dp_ingest`'s own
argparse defaults -- RHO_BUCKET, RHO_STATUS_BUCKET, RHO_INGEST_THREADS,
RHO_PASS_OBJECTS, RHO_METRIC_NAMESPACE -- so this passes only what is true
of *being a scheduled invocation* rather than restating that interface.

Two behaviours that are not incidental:

`--no-index` after the first run. dp_ingest creates the `found_at` index at
startup; that is right for a daemon that starts once and wrong for a function
that starts every two minutes. Set RHO_INGEST_INDEX=1 for the first
invocation against a new store, then leave it unset.

A failed pass raises. Returning quietly would make the function's `Errors`
metric useless, and a scheduled function nobody can alarm on is exactly as
silent as the dead EC2 instance this replaces -- which went unnoticed for
four hours on 2026-09-20 while a downstream workflow published its frozen
output, green, 802 times.
"""

from __future__ import annotations

import os

import dp_ingest


def build_argv(env=None):
    """The argv for one scheduled pass. Separate so a test can assert it."""
    env = os.environ if env is None else env
    argv = ["--once"]
    # dp_ingest's own env defaults cover bucket, status bucket, threads,
    # pass-objects and metric namespace. Do not restate them here.
    if not _truthy(env.get("RHO_INGEST_INDEX")):
        argv.append("--no-index")
    prefix = env.get("RHO_INGEST_PREFIX")
    if prefix:
        # Rehearsal lever: scope a run to one slot without a code change.
        argv += ["--prefix", prefix]
    return argv


def _truthy(value):
    return str(value or "").strip().lower() in {"1", "true", "yes", "on"}


def handler(event, context):
    argv = build_argv()
    rc = dp_ingest.main(argv)
    if rc != 0:
        # Raise rather than return: this is what puts a point on the Errors
        # metric, and Errors is what an alarm can see.
        raise RuntimeError("dp_ingest pass failed with rc=%d (argv=%r)" % (rc, argv))
    return {"rc": rc, "argv": argv}

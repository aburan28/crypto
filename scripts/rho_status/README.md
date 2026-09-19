# ECC2K-130 public status snapshot

GitHub Action, run every 15 minutes, that reads the private Pollard-ρ distinguished-point
store and publishes **aggregates only** to GitHub Pages.

The page never includes point keys, walk coefficients `(a, b)`, seeds, or
the database URL.

## Why the Action hops through the walker

RDS `rho-dp` is **not** publicly accessible. Security group `rho-dp-rds`
allows TCP 5432 only from VPC worker groups, including
`rho-ecc2k-walker`. GitHub-hosted runners therefore SSH to the running
instance tagged `Name=rho-ecc2k-walker` and run
`scripts/rho_status/snapshot.py` there, using `/opt/rho-ecc2k/env.sh`
already on that host.

`snapshot.py` does not scan `distinguished_points` on every run. It
creates `rho_dp_hour` / `rho_dp_recent`, backfills them once from the
heap, and installs a statement-level insert trigger so the ingest's
`INSERT ... SELECT` keeps the buckets current. Later snapshots read those
tables.

The first merged backfill (#448) did not unstick Pages. It set
`statement_timeout` to 800 s, took SHARE on the heap, grouped by
`(worker_id, found_at)`, and PostgreSQL cancelled it; `rho_dp_meta.ready`
stayed false and the hop kept retrying, so
https://aburan28.github.io/crypto/status/ froze on the 2026-09-18T11:25Z
snapshot. The hop now writes an index-only fallback *before* the worker
backfill: it counts one hour at a time through `(campaign_id, found_at)`
(no `worker_id`, no full-table `GROUP BY`) so a 187 M-row hash aggregate
cannot OOM the walker. The 2026-09-18T18:56Z publish died ~200s into that
GROUP BY (`Connection reset by peer`); counting by hour logs each step
and keeps SSH alive. The file is copied even if the remote `timeout`
fires, from a per-run path so leftover `/tmp/rho_status.json` is never
republished. Worker buckets backfill one hour at a time, resuming from
`rho_dp_meta.backfill_through`. The remote timeout is 1500 s; leftover
hours wait for the next scheduled run.

That still left Pages on 11:25Z: from 19:15Z on 2026-09-18 every scheduled
publish died in one second with `no running instance tagged
Name=rho-ecc2k-walker`. The 18:56Z hop reset SSH and the instance did not
come back. RDS is private; without that host the Action cannot query the
store. The ingest host, though, already writes the same campaign counts to
the public status bucket (`dp_ingest.py` `publishStatus`) every
`--status-every` seconds, and `work_feed.py` already fetched that file
*after* a successful hop. When the hop finds no walker — or finishes
without a file — `fetch_via_walker.sh` now republishes that ingest
document (`work_feed.py --as-snapshot`): counts, last-hour, hourly (48 h),
and the checkpointed walk total, with `work.per_slot` stripped. The
per-worker table is empty on that path; the GPU card still reads
`walking_slots`. `generated_at` is the ingest host's stamp, so the stale
banner tracks the feed rather than the copy.

Do not open `0.0.0.0/0` on the RDS security group for this dashboard.

## One-time GitHub + IAM setup

The publish job authenticates to AWS with a dedicated IAM user
`ecc2k130-status-gha` (describe walker + TCP/22 on the tagged walker
SG only). Do not upload the `adam` user keys.

1. **Create the service user and GitHub secrets** (administrator profile;
   `adam` cannot create IAM users or access keys):

   ```bash
   AWS_PROFILE=admin ./scripts/rho_status/gha_iam_user.sh --push-github
   ```

   Preferred longer-term: GitHub OIDC via
   `AWS_PROFILE=admin ./scripts/rho_status/gha_iam_role.sh` and switch the
   workflow back to `role-to-assume`.

2. **Secrets** (repo or environment `github-pages`):
   - `AWS_ACCESS_KEY_ID` / `AWS_SECRET_ACCESS_KEY` — from
     `ecc2k130-status-gha`, not from `adam`.
   - `RHO_WALKER_SSH_KEY` — the GHA deploy key (ed25519). You do **not**
     need the instance launch key (`meow34`). The matching public key is
     `scripts/rho_status/gha_walker.pub` and is installed on
     `ubuntu@walker`. The job opens TCP/22 from the runner's public IP
     for the duration of the hop, then revokes that rule.

3. **Pages**: Settings → Pages → Source = **GitHub Actions**.

4. If the walker is replaced, install `gha_walker.pub` on the new box
   and keep instance tag `Name=rho-ecc2k-walker` plus SG tags
   `Name=rho-ecc2k-walker` and `Purpose=ecc2k-dp-walker`.

The workflow does **not** need `DATABASE_URL`. The walker already has it.

## Local commands

```bash
# Offline tests (no database)
python3 scripts/rho_status/test_rho_status.py

# From a host that can reach RDS (the walker, or an SSH tunnel)
python3 scripts/rho_status/snapshot.py --out /tmp/status.json

# From a laptop with AWS creds + the GHA deploy key
RHO_WALKER_SSH_KEY=$HOME/.ssh/gha_walker \
  scripts/rho_status/fetch_via_walker.sh
```

After a snapshot, attach the campaign's iteration total and measure the
rate from it:

```bash
# Reads the campaign work feed; the bucket is derived from the calling AWS
# identity, or given outright with RHO_WORK_FEED_URL.
python3 scripts/rho_status/work_feed.py --status /tmp/status.json

python3 scripts/rho_status/render.py \
  --status /tmp/status.json \
  --history-in docs/ecc2k130-status/history.json \
  --history-out docs/ecc2k130-status/history.json \
  --status-out /tmp/status.json
```

The order matters: the rate is a difference between two snapshots, so the
total has to be in the snapshot before `render.py` records its history
point. Run the other way round and every point publishes without a total,
and no rate is ever measurable.

Open `docs/ecc2k130-status/index.html` next to the JSON files, or assemble
the whole site the way the Action publishes it:

```bash
python3 scripts/site/build.py --out _site
python3 -m http.server --directory _site 8000   # dashboard at /status/
```

## Published fields

`status.json` includes campaign id, curve id, DP counts, collision count,
per-worker counts, 7-day hourly buckets, and — when the campaign's work
feed is answering — the `work` and `walk_rate` blocks described below.
`state` is one of `COLLECTING`, `IDLE_OR_STALE`, `EMPTY`,
`COLLISION_RECORDED`, `INGEST_BEHIND`. `INGEST_BEHIND` is copied from the
ingest host's work feed when the store's last-hour count is zero but the
feed still has outstanding or unreadable objects — the 2026-09-17 outage
looked like a dead walk on this page until that distinction existed.

A recorded collision is **not** treated as a solved discrete log on the
page. Independent verification of `[k]P = Q` is still required.

## Iterations per second, and why it is not derived from points

`work_feed.py` adds a `work` block to `status.json` — the iteration total
the walkers checkpointed, how many slots are behind it, and how old the
feed was — and `render.py` adds `walk_rate`: that total's increase since an
earlier snapshot, the span it was measured over, and the two totals it
differenced. Both pages render the rate with its span, because a rate
without one is not checkable: on this fleet a 15-minute window swings by a
third while an hour does not, so `render.py` prefers the oldest sample
inside `RATE_WINDOW_S` (an hour) and reaches further back only when an
outage left the window empty.

It is a conservative figure by construction. A slot's iterations only count
once it has checkpointed them, so the walked-but-not-yet-checkpointed tail
is missing: measured against the live campaign on 2026-09-15 it read
72.5 B it/s where `ecc2k130/aws/status.py`, summing what the walkers
reported about themselves, read 86–101 B it/s over the same period. The
published number is the one that survives a worker dying mid-interval.

Neither field is guaranteed. `work_feed.py` refuses a feed that is stale,
foreign or unparseable and exits 0, because a frozen total does not read as
a dead ingest host, it reads as a dead campaign. With no total the pages
fall back to the point count at one point per `2^25.27` iterations at
`HW(x) <= 34` (`ecc2k130/aws/README.md`) and label the figure as such. That
fallback **reads low**: this campaign distinguishes at the sparser
`HW(x) <= 32`, about `2^28.4` iterations per point measured from each
client's own iteration and point counters
(`ecc2k130/benchmarks/dp-interval/`; the `2^27.9` first quoted was a ratio
of fleet-wide sums and is retired), so the derivation understates the walk
by a factor near nine and understated it on this page until the total was
published. Both exponents are
parameters of one campaign's distinguishing rule, so the pages apply them
to `ecc2k-130` only, and `scripts/site/test_build.py` pins every copy of
them to the campaign document.

The dashboard's progress bar is filled from the ratio of the **work**,
`2^(n - 60.9)`, never from the ratio of the exponents: `n = 39` is two
thirds of the way along the exponent and about a millionth of a millionth
of the way through the work. A test pins that, because it is the kind of
bar that gets "fixed" into a lie by anyone trying to make it look fuller.

The dashboard's ETA to that same `2^60.9` expected cost is the remaining
work divided by the measured rate above, so it moves with the fleet. A
snapshot with **no** iteration total falls back to the last-hour
distinguished-point amount, `(dps_last_hour × 2^25.27) / 3600` operations
per second, and says which of the two it used. The two accountings are
never mixed: a counted total divided by a derived rate would publish an ETA
about nine times too long, so a snapshot that has a total but not yet a
second one to difference shows no ETA rather than that. The ETA is an
expectation, not a deadline; `scripts/site/test_build.py` pins both
formulas, their order, and the refusal to mix them.

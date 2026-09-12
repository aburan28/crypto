# ECC2K-130 public status snapshot

Hourly GitHub Action that reads the private Pollard-ρ distinguished-point
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

Do not open `0.0.0.0/0` on the RDS security group for this dashboard.

## One-time GitHub setup

1. **Secrets** (repo or environment `github-pages`):
   - `RHO_WALKER_SSH_KEY` — private key for `ubuntu` on the walker
     (the `meow34` key the instance was launched with).
   - `AWS_ACCESS_KEY_ID` / `AWS_SECRET_ACCESS_KEY` — IAM user that can
     `ec2:DescribeInstances` in `us-west-2`, so the job can find the
     current spot public IP after a replacement.
   - Optional: `RHO_WALKER_HOST` if you want to pin a host and skip the
     lookup.
2. **Pages**: Settings → Pages → Source = **GitHub Actions**.
3. Keep SSH ingress on the walker SG limited to operators. The Action
   uses the same key; if you later rotate the spot box, retag the new
   instance `Name=rho-ecc2k-walker` and allow the operator IP again.

The workflow does **not** need `DATABASE_URL` in GitHub. The walker
already has it.

## Local commands

```bash
# Offline tests (no database)
python3 scripts/rho_status/test_rho_status.py

# From a host that can reach RDS (the walker, or an SSH tunnel)
python3 scripts/rho_status/snapshot.py --out /tmp/status.json

# From a laptop with AWS + walker SSH key
RHO_WALKER_SSH_KEY=$HOME/.ssh/meow34.pem \
  scripts/rho_status/fetch_via_walker.sh
```

After a snapshot:

```bash
python3 scripts/rho_status/render.py \
  --status /tmp/status.json \
  --history-in docs/ecc2k130-status/history.json \
  --history-out docs/ecc2k130-status/history.json
```

Open `docs/ecc2k130-status/index.html` next to the JSON files.

## Published fields

`status.json` includes campaign id, curve id, DP counts, collision count,
per-worker counts, and 7-day hourly buckets. `state` is one of
`COLLECTING`, `IDLE_OR_STALE`, `EMPTY`, `COLLISION_RECORDED`.

A recorded collision is **not** treated as a solved discrete log on the
page. Independent verification of `[k]P = Q` is still required.

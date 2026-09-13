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

After a snapshot:

```bash
python3 scripts/rho_status/render.py \
  --status /tmp/status.json \
  --history-in docs/ecc2k130-status/history.json \
  --history-out docs/ecc2k130-status/history.json
```

Open `docs/ecc2k130-status/index.html` next to the JSON files, or assemble
the whole site the way the Action publishes it:

```bash
python3 scripts/site/build.py --out _site
python3 -m http.server --directory _site 8000   # dashboard at /status/
```

## Published fields

`status.json` includes campaign id, curve id, DP counts, collision count,
per-worker counts, and 7-day hourly buckets. `state` is one of
`COLLECTING`, `IDLE_OR_STALE`, `EMPTY`, `COLLISION_RECORDED`.

A recorded collision is **not** treated as a solved discrete log on the
page. Independent verification of `[k]P = Q` is still required.

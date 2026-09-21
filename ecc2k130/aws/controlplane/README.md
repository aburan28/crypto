# The campaign control plane: RDS, ElastiCache and S3

**This package makes no performance claim.** Under AGENTS.md §3 it is
plumbing: no row on `docs/index-calculus-scoreboard.html` moves because of it
and no measurement is offered. What it is allowed to claim is the failure
modes it removes, and those are named below with the mechanism that removes
each one.

The campaign already ran on three stores. What it did not have was one place
that said how they relate, so the rules lived in `worker.py`, in
`dp_ingest.py` and in a shell script, and two of them could disagree without
anybody finding out until a GPU wrote into a run id somebody else was
walking.

| store | holds | if it is lost |
|---|---|---|
| RDS Postgres | slot leases, fencing tokens, the object ledger, checkpoint pointers, the worker roster | the campaign stops claiming and reporting; the corpus is intact |
| ElastiCache Redis | scan snapshots, fleet heartbeats, the DP hot set, fleet-wide locks | nothing; every call falls back to Postgres |
| S3 | the corpus and the checkpoints — the bytes | the campaign is over |

## The three invariants

**Postgres is the authority and Redis is a hint.** Nothing in this package
reads a Redis miss as evidence. `Cache.fleet()` returns `None` for "cannot
say" and `[]` for "nobody is beating", and the caller has to tell them apart;
`seenPoint()` returns `None`, `True` or `False` with the same discipline, and
a `False` never means a point is new — the hot set remembers one hour and the
corpus is 2^32.5 records. A cache that can change an answer is not a cache.

**S3 holds the bytes; the database holds what is true about them.** The
object is uploaded first and the ledger row is committed after, never the
reverse. An object with no row is invisible until the next pass, which is a
delay; a row with no object is a record of points that do not exist, which
sends the merge looking for them.

**A lease is fenced.** Every successful claim increments the slot's fence.
Every write a stale owner could still make — a ledger row, a checkpoint
pointer — carries the fence it was issued under, in the `WHERE` clause of the
statement that makes it, and the database rejects it if the slot has moved on.

That last one is the reason this package exists. A 180-second lease cannot
stop a process that was frozen past it — a hypervisor pause, memory pressure
on a spot instance, a partition that healed — from waking up and writing as
though it still owned a run id. `test_controlplane.py` runs exactly that
scenario: worker A is frozen, B claims and walks the slot, A comes back, and
every write A attempts raises `FenceError` because the database noticed, not
because A did.

## What each component does

### `slots.py` — the registry

A fourth backend under the same four-method contract as `worker.py`'s
DynamoDB, S3 and local-file registries (`scan`, `claim`, `heartbeat`,
`release`), so it is a drop-in: `ECC_SLOT_BACKEND=rds` and nothing else
changes. A test asserts the signatures still match, because that drift is
found on a GPU at 3 a.m. otherwise.

Three things it can do that the others cannot:

* **A claim is one statement.** DynamoDB scans and then conditionally updates
  per candidate; the S3 registry does a read-modify-write with `If-Match` per
  candidate. Here it is one conditional `UPDATE ... RETURNING`, so the
  database adjudicates and a fleet-wide spot interruption costs one round trip
  per worker instead of one per candidate it loses. The race is tested with
  sixteen threads against one database, twice: once racing on creation, once
  racing on resumption.
* **Leases are fenced** (above).
* **One store fewer.** The registry, the ledger and the dashboard's view are
  in one database, so "which worker wrote this object, under which lease" is a
  join and not a correlation by eye across DynamoDB and S3.

The claimability rules — which slots a CPU worker may resume, which GPU family
may take which checkpoint — are **imported from `worker.py`**, never restated.
Two copies of the rule that keeps a packed checkpoint away from the host
client is how a run id gets retired and its seeds lost.

### `cache.py` — ElastiCache

Everything fails open, with a circuit breaker: 0.5-second timeouts, no
retries, and after three consecutive failures the cache is skipped entirely
for 30 seconds. The locks are advisory and grant on a cache failure — the
fleet's work does not stop because the *cache* is down — which is why nothing
that must happen exactly once depends on holding one. Fencing is the slot's
job, not the lock's.

### `store.py` — S3

Objects are immutable and content-addressed, per `protocol.py`. Re-writing
identical bytes to a key is success (that is what a retried upload looks
like); writing *different* bytes to an existing key raises `ObjectExists`
rather than going through. `LocalObjectStore` is the same interface over a
directory for rehearsals.

### `db.py` — RDS

One connection per thread, never shared. A Multi-AZ failover kills every open
connection, so each statement reconnects and retries on a connection-shaped
error and on nothing else — a syntax error or a constraint violation is a bug
and is raised at once. The retries are safe because every statement in this
package is conditional or `ON CONFLICT DO NOTHING`; a caller adding a
`SET x = x + 1` must use `transaction()`, which does not retry. The backoff is
jittered, because a failover wakes the whole fleet at once.

### `schema.py` — the migrations

Append-only, idempotent, guarded by a Postgres advisory lock so a rollout that
runs `migrate` on thirty hosts serialises instead of racing through
`CREATE INDEX`. A migration that dies halfway re-runs harmlessly: every
statement is `IF NOT EXISTS` and the marker row is written last.

The dashboard's existing tables (`distinguished_points`, `dp_ingest_*`,
`rho_collisions`) are **not** touched. They are a derived view with a writer
already; a second writer on them would buy nothing.

## Configuration

```
DATABASE_URL=postgresql://user:pw@rho-dp...:5432/rho?sslmode=require
  # or RHO_DB_HOST plus RHO_DB_SECRET (default rho/dp-rds) in Secrets Manager
RHO_REDIS_URL=rediss://<primary-endpoint>:6379   # rediss:// => TLS
ECC_BUCKET=<campaign bucket>
ECC_PREFIX=                                       # optional key prefix
RHO_CAMPAIGN=ecc2k-130
RHO_CACHE=0                                       # run with no cache (the control)
ECC_LOCAL_STORE=<dir>                             # rehearsal: SQLite + a directory
ECC_S3_ENDPOINT=http://127.0.0.1:9010             # an S3-compatible endpoint
ECC_S3_REGION=us-east-1                           # or AWS_REGION
```

`ECC_S3_ENDPOINT` is empty on every deployed host, which means real S3. Set it
and the client is built against that endpoint with path-style addressing,
which is what MinIO needs and what an on-premises store or GCS's S3
interoperability API would be pointed at. Only MinIO and moto are tested (see
`integration/`); nothing here claims GCS works until something runs against
it.

An absent setting is an absent component, not an error. With no
`RHO_REDIS_URL` the coordinator reads Postgres for everything it would have
cached, which is also the configuration a cache claim would have to be
measured against.

`DATABASE_URL` and `RHO_DB_SECRET` are resolved exactly as `dp_ingest.py`
resolves them, on purpose: two spellings is how one campaign ends up with two
databases.

### ElastiCache settings this schema assumes

* **`maxmemory-policy: volatile-lru`.** Every key here carries a TTL, and
  `allkeys-lru` would let a burst of DP hints evict the fleet view.
* **Cluster mode:** keys carry the hash tag `{campaign}`, so one campaign's
  keys land on one shard and the fleet view is one round trip.
* **Encryption in transit on**, via `rediss://`.
* Sizing follows the fleet size and the DP rate, not the corpus: the hot set
  is one hour of points, and at the measured weight-32 rate (~40 points/s per
  GPU) that is small.

### RDS settings

* **Multi-AZ**, because `db.py`'s retry is written for a failover and a
  single-AZ instance turns one into an outage.
* The control-plane tables are small and hot; the corpus tables are not. They
  share an instance today. If the ingest's write load ever interferes with a
  claim, the control plane is the part that moves, because it is the part with
  no history to carry.
* IAM auth or Secrets Manager, never a password in a unit file.

### IAM

The worker's instance role needs `s3:GetObject`/`PutObject`/`ListBucket` on
the campaign bucket and reach to the database and the cluster. It does **not**
need DynamoDB when `ECC_SLOT_BACKEND=rds`, which is one service and one policy
fewer per host.

## Using it

```python
from controlplane import Coordinator
from controlplane.coordinator import workerInfo

co = Coordinator.fromEnv(migrateSchema=True)
lease = co.acquire(workerInfo(gpu=0, gpuName=name, gpuFamily=family,
                              instance=instanceId, instanceType=instanceType))
co.setCheckpoint(lease, "ck/slot-00007/000123.ck", iters=1_000_000)
co.publishObject(lease, "dp/slot-00007/<stream>-<offset>-<sha>.bin", data)
while walking:
    if not co.heartbeat(lease, {"iters": iters, "points": points, "rate": rate}):
        break            # the lease is gone; another worker holds this run id
co.release(lease, state="idle")
```

Command line:

```bash
python3 -m controlplane migrate     # idempotent; safe on every boot
python3 -m controlplane doctor      # RDS, ElastiCache and S3, reported separately
python3 -m controlplane status      # slots, corpus, fleet, cache
python3 -m controlplane slots --events --limit 50
python3 -m controlplane retire 140 --reason "checkpoint refused (exit 6)"
```

`doctor` reports the three stores separately because "the campaign is down" is
almost always one of them being down and the other two being fine.

## Tests

Two suites, and the split between them is deliberate.

```bash
python3 -m unittest discover -s aws -p test_controlplane.py -v   # anywhere
make -C ecc2k130 test-integration                                # real services
```

The first has no network, no AWS SDK and no Redis: the SQL runs against SQLite
(which takes the same `ON CONFLICT`, `RETURNING` and `INSERT ... SELECT ...
WHERE EXISTS` forms), the cache runs against a fake with the failure modes
that matter, and one test asserts in a clean subprocess that importing the
package pulls in neither `redis` nor `boto3`. The lease races are run with
threads, not asserted about. It must keep passing on a laptop with nothing
installed, which is why it is the one wired into `make test-production`.

The second runs against real Postgres, real Redis and a real S3 endpoint
(`integration/`). It re-runs the first suite's own test bodies against those
services — same assertions, different fixtures — and adds what a substitute
cannot show. Two bugs it found on its first run, both of them invisible to
SQLite and a fake:

* `migrate()` created its own bookkeeping table *before* taking the advisory
  lock. `CREATE TABLE IF NOT EXISTS` is not atomic in Postgres, so a rollout
  running `migrate` on every host at once left one of them holding a unique
  violation on `pg_type`. The lock is now taken first.
* redis-py 6 and later retry a failed connection **ten times with backoff by
  default**. A cluster that was not there cost 4.5 seconds per call against
  the 0.5 the package sets and documents — on the claim path, three times over
  before the cooldown had seen enough failures to skip. `Cache.fromEnv` now
  builds the client with no retries, and a test asserts the bound.

## Not here

* Terraform for the instance, the cluster and the parameter groups.
* Moving the live campaign's DynamoDB registry into RDS. `ECC_SLOT_BACKEND`
  chooses a backend for a *new* fleet; importing an existing slot table is a
  migration with a cutover, and a cutover that two backends could both serve
  during is exactly the double-claim this package exists to prevent.
* Collision detection. `merge.py` over the S3 corpus and `dp_ingest.py`'s
  unique-key check remain the two detectors; the hot set here is a hint and
  deliberately not a third one.
* Any measurement. See the first paragraph.

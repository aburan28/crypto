# Integration services

The control plane's own suite (`aws/test_controlplane.py`) runs its SQL
against SQLite, its cache against a fake and its object store against a
directory. That is what lets it pass on a laptop with nothing installed, and
it is also the reason it cannot see:

* what two hosts running `migrate` at the same moment do to each other —
  SQLite has no `pg_advisory_lock`, and Postgres's own `CREATE TABLE IF NOT
  EXISTS` is not atomic against a concurrent one;
* a `bigint` that must hold 2⁶² — SQLite is dynamically typed, so an
  `integer` column would pass there and overflow in the fleet's first week;
* a TTL the *server* owns, the unlock script Redis runs itself, or what a
  current redis-py does with a connection it cannot make;
* a listing longer than one 1,000-key S3 page, which is where the hand-rolled
  continuation-token loop in `store.py` lives.

`aws/test_integration.py` covers those, and re-runs the unit suite's own test
bodies against the real services while it is there.

## Running it

```bash
make -C ecc2k130 test-integration-deps    # psycopg, redis, boto3 — once
make -C ecc2k130 test-integration         # up, run, down
```

`services.sh up` starts Postgres, Redis and an S3 endpoint, waits until each
one answers, and writes `ecc2k130/build/integration/env`. `down` stops them
and takes their data with it. Both are idempotent: a service already
answering on its port is left alone.

Ports default to 5433, 6380 and 9010 rather than 5432, 6379 and 9000, so a
run cannot reach a developer's own database. Override with `ECC_IT_PG_PORT`,
`ECC_IT_REDIS_PORT`, `ECC_IT_S3_PORT`.

## Two ways to get the services

**Docker** (`docker-compose.yml`) is the intended way and the one CI uses:
`postgres:16-alpine`, `redis:7-alpine`, `minio/minio`. Nothing is mounted and
`down` removes the volumes, so there is no state to corrupt.

**Without a Docker daemon** the same script starts the host's own
`postgres`, `redis-server` and `moto_server` into `ecc2k130/build/integration/`
instead. It is not a fallback for its own sake: a suite that only runs in one
place is a suite that rots in the other, and this is the path the work was
developed on. `ECC_IT_NATIVE=1` forces it.

Either way the suite sees three URLs and nothing else, which is the point —
`ECC_S3_ENDPOINT` exists so the same code path that talks to S3 in production
is the one under test here.

## What the suite will not do

It will not skip quietly. With `ECC_INTEGRATION` unset the module skips as a
whole and says why; with it set, a service it cannot reach is an error. A
green run that tested nothing is worse than a red one.

It does not test GCS. `ECC_S3_ENDPOINT` is what a GCS S3-interoperability
endpoint would need, but nothing here has run against one, so nothing here
claims it works.

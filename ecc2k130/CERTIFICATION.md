# ECC2K-130 production-readiness gates

**Scope: local software evidence, not live production certification.**
This change is correctness engineering, not an arithmetic speedup, a new
attack, or a claim that the ECC2K-130 challenge has been solved. Source tests
do not certify a deployed binary, a database, or an existing corpus.

## Reproduce the local gates

From `ecc2k130/`, with a C++17 compiler, OpenMP and Python/NumPy installed:

```sh
make cpu build/test-production
make certify-local
bash aws/rehearse_merge.sh
```

`build/certification-local.json` records each command, return code, elapsed
time, log hash, source file hashes, and exact client/fixture binary hashes.
`productionCertified` remains **false**, even when every local check passes.
`python3 aws/certify.py --output build/review.json --require-production`
exits 2 after passing local tests: a local-only report cannot satisfy a
production deployment gate. This is deliberate, not a pending test hidden
behind a green result. The CI job is named **local certification gates**.

The [2026-09-14 local evidence](certification/local-20260914.json) and its
adjacent hashed logs record 34 passing storage/recovery tests, the reference
suite, C++ persistence fault tests, and packed-host/timing tests. This report
identifies the modified working tree through individual source hashes; it
does not describe an unchanged release commit or a deployed GPU binary.

### Acceptance contract

| Gate | Evidence required |
|---|---|
| Arithmetic/reference | Existing eight-field reference suite, including GF(2^131), passes |
| Collision recovery | Fixed GF(2^41) pair, separate worker corpora, independently replayed endpoints, recovered `369250562913`, verified `[k]P == Q` |
| Persistence | ENOSPC/EIO simulation and report overflow stop with no advanced checkpoint |
| Compatibility | Foreign campaigns, altered cutoff/build hashes, unpinned replay binaries and implicit legacy adoption are rejected |
| Integrity | Altered chunks/checkpoints/buckets, torn records, missing committed state and source truncation are rejected |
| Coordination | Single writer, stale lease rejection, conditional checkpoint publication, upload-order tests |
| Recovery | Interrupted ingest, retry after failed upload/solve, safe file rotation and a real CPU worker resuming on a replacement local host |
| Packed implementation | Host differential tests for the packed arithmetic variants; **device execution is separate** |

One failing command fails local readiness. No timeout, inconclusive random
rehearsal, report drop, corrupt input, or missing completion counts as a pass.
The fixed collision fixture is emitted by `build/test-production --fixture`;
its keys are computed by replay, not made equal by fabricating endpoint bytes.
The two frozen seeds are `0002000048880000` and `000100004cf60000`, with the
normal-basis curve-41 default cutoff. The test proves each separately reaches
the recorded orbit before giving the pair to the actual merge/client solver.

## Safety changes

The binary payload remains **32 bytes**: a little-endian 64-bit starting seed
and three little-endian 64-bit words of the canonical normal-basis x key.
There are **no stored scalar coefficients**. The resolver reconstructs them
by replay, aligns Frobenius/negation, and verifies `[k]P == Q`.

- The client requires aligned regular DP files and an exclusive writer lock.
  Failed open/write/flush/fsync/checkpoint writes return nonzero. A report
  overflow stops **before reseeding or advancing a checkpoint**, instead of
  logging a warning and continuing. Flush and fsync precede checkpoint
  publication; checkpoint rename is followed by directory fsync.
- Cycle-guard restarts no longer become false DP reports on the bitsliced
  CPU/CUDA path. The packed path already distinguished these cases.
- The 16-bit per-walk restart counter cannot wrap into the neighboring
  walk's seed range. Exhaustion returns 9 and stops the worker for review;
  allocate a fresh run ID rather than resuming that stream indefinitely.
- Strict manifests bind curve, cutoff, guard, seed/walk/key semantics,
  worker geometry, source hash, worker binary hash and replay binary hash.
  SHA-256 is over the full immutable payload. These checks detect accidental
  corruption and configuration mismatch; they are **not authentication
  against an actor able to rewrite both data and manifests**. IAM and
  protected campaign configuration remain part of the live threat model.
- DP payloads are uploaded before their `.bin.json` commit manifests, and
  both precede checkpoint publication. Uncommitted chunks remain pending;
  they are not included in accepted counts and are retried next pass.
- Checkpoints use immutable content-addressed `.ck` objects plus manifests.
  The slot registry's **conditional, owner-checked heartbeat** publishes
  `checkpointKey`. A stale worker can leave orphaned blobs but cannot commit
  the pointer. Resume follows that pointer and verifies the checksum before
  the client validates checkpoint shape. Owner IDs include a unique process
  epoch. Expired owners cannot renew their leases.
- Failed object access is not interpreted as a nonexistent checkpoint.
  Local supervisors and merge processes are exclusively locked.
- Merge uses a mixed hash of all key limbs, rather than raw low bits of a
  sparse coordinate. Ingest is chunked; a dirty bucket's sort still needs
  RAM proportional to that bucket. Checksums protect derived bucket files.
  A crash with partially updated derived buckets can require a **fresh merge
  work directory**, rebuilt from preserved immutable sources. Do not delete
  source corpora. This is fail-closed recovery, not uninterrupted availability.
- Transient solve failures remain retryable. A verified solution additionally
  requires successful client exit; a bare `k = ...` log line is insufficient.

### Failure codes

| Exit | Meaning | Action |
|---|---|---|
| 3 | Reference replay mismatch | Stop, preserve data, investigate arithmetic |
| 6 | Incompatible/incomplete checkpoint | Preserve it; correct configuration or audit recovery |
| 7 | Report overflow | Preserve old checkpoint; increase buffer or reduce launch size; replay lost launch |
| 8 | Persistence failure | Repair storage; inspect/quarantine a partial tail; restore from a verified checkpoint |
| 9 | Seed counter exhausted | Preserve corpus; stop that run ID and allocate a fresh stream |

Strict workers do not automatically retry these integrity failures. Returning
nonzero does not itself repair media or recover incomplete filesystem data.

## Deployment and compatibility

New campaign configuration uses `storageProtocol: "ecc2k-seed-orbit-v1"`.
`build.sh` fills in all three required 64-character hashes. `worker.py` checks
its executable against `binarySha256`; merge checks its CPU resolver against
`hostBinarySha256`. `protocol.py` must be deployed alongside `worker.py` and
`merge.py`; bootstrap/infra scripts now include it.

```sh
python3 aws/merge.py --work /data/merge-v2 \
  --s3 s3://CAMPAIGN_BUCKET/dp/ --campaign campaign.json \
  --client ./ecc2k130-cpu
```

**Do not turn this on in an existing campaign without a migration review.**
Old data has no manifest proving its cutoff, walk version or producer build.
The strict path intentionally refuses to bless it by merely writing a new
manifest. Preserve the old bucket/checkpoints; use a new campaign namespace
and fresh work directories for the strict deployment, or develop and audit a
migration that establishes the old data's provenance and replay invariants.

`merge.py --legacy ...` explicitly processes old raw corpora. The old bucket
layout requires a fresh derived merge directory because the bucket hash
changed. `ECC_ALLOW_LEGACY_STORAGE=1` permits the old worker path only for
legacy recovery/rehearsals; it is **not certified**. No live data was migrated,
removed, uploaded, or changed during this implementation.

The unmerged `rds_gpu.py` adapter is **not included in this certification**.
Its `(a=seed,b=0)` fields are not Pollard coefficients. A future RDS integration
must use a separate, typed **seed-orbit** record table/candidate queue, preserve
both seeds and protocol identity, and invoke this pinned replay resolver.
Do not import these manifests into a generic coefficient-based solver.

## Required live gates before a production sign-off

1. Build the exact pinned CUDA preset. Run the device arithmetic, compact
   storage, shared-sigma fixtures, `codegen/testpackedclient.py`, and bidirectional
   checkpoint tests on the actual GPU model. Include forced report overflow,
   guard restart and seed exhaustion. Host tests cannot cover PTX/SASS, races,
   driver behavior, or device faults. Recheck performance: this change makes
   no claim to preserve the old measured throughput.
2. Collect **DP34** on the production binary with zero drops. Independently
   replay complete sampled seeds through the pinned CPU resolver; verify
   endpoint, canonical key and first-DP semantics. Record replay latency.
   The quick high-cutoff GPU test is not a substitute for full-length replay.
3. Use an isolated cloud test campaign for lease expiry/takeover, failed and
   delayed uploads, stale writers, SIGKILL/replacement, interrupted file
   rotation, missing/corrupt objects, denied reads and simulated network
   partitions. Verify all committed DPs remain mergeable and checkpoint
   publication cannot pass missing DP chunks. Do not fault-inject against
   the user's live campaign without a separately agreed test scope.
4. Measure merge capacity with the **actual sparse canonical-key distribution**,
   max bucket size, RSS, disk headroom, ingest lag, sort/checksum cost and
   recovery time. Local toy tests are not a terabyte-scale load test.
5. Inspect deployed binary/source hashes, config, active leases, manifests,
   backlog and stored corpus. Sign off against that immutable deployment
   identity, not against `main` generally. Certify RDS separately if used.

For sizing only, at cutoff 34 the uniform trace-zero-coordinate approximation
gives `p = sum(C(131,k), even k <= 34) / 2^130 ≈ 2^-25.2697`. Combining this
with the historical planning budget `2^60.9` iterations gives approximately
53.2 billion records or **1.70 TB of raw 32-byte payload**, excluding metadata,
duplicates and tail-work effects. An S3 cache and derived buckets on one disk
already need roughly **3.4 TB before sort/workspace headroom**. These are
extrapolations, not measured capacity. Cutoff 33 is effectively cutoff 32
because subgroup x weights are even; lowering the cutoff produces fewer DPs
and longer walks, not the reverse.

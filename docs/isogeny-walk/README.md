# Isogeny walk: user guide

`isogeny_walk` explores the isogeny class of a prime-field elliptic curve
(P-256, P-224, P-192, or any short Weierstrass curve you give it).

- **Discovery.** It finds every curve reachable through small-degree
  isogenies, certifies each isogeny with an explicit kernel polynomial,
  and records every curve in the repository's curve formats: ICV1 slug,
  EC1 identity, `curves.yaml`, and `IW1` routes.
- **Traits.** It runs trait detection on every curve it finds.
- **Storage.** It can keep everything in S3.

This page is the how-to.  The design and the evidence are in
[`docs/curves/ic/README.md`](../curves/ic/README.md#walking-an-isogeny-class)
and the research notes
[`isogeny_walk_20261004`](../../research/isogeny_walk_20261004/README.md),
[`isogeny_walk_traits_20261004`](../../research/isogeny_walk_traits_20261004/README.md) and
[`isogeny_walk_p192_20261004`](../../research/isogeny_walk_p192_20261004/README.md).

## Is it ready?

**Ready to use:**

- **Walks of any size your machine's memory allows.**  Every edge is
  kernel-certified, every curve's order is proved for prime-order roots,
  and the outputs are deterministic, so the same arguments give the same
  bytes and the same run id.
  - The walker has run 20,000-curve walks of P-256, P-224 and P-192 with 0
    failures across 325,000 certified edges.
- **`verify`**, an independent replay of a whole walk from its records.
- **Trait detection** on every curve, as local shards or queued jobs.
- **S3 storage**: write-once and hash-checked, with fetch and later
  publish.

**Not there yet**, so plan around it:

- **One machine per walk.**  A walk is a single breadth-first search held
  in memory, at about 87 KB of RAM per curve (measured).  Only trait
  detection is split across machines.
- **No resume.**  A walk that is killed starts over.  Use `--max-curves` to
  keep runs to a size you can afford to repeat.
- **Prime fields and odd `ℓ` only.**  No binary or Koblitz curves (for
  example ECC2K-130) and no 2-isogenies.  The prime-order NIST curves have
  no rational 2-torsion, so they have no 2-isogenies anyway.
- **`ℓ ≤ ~61` in practice.**  `Φ_ℓ mod p` is rebuilt at every start, which
  takes about 10 s for `ℓ = 59` and grows like `ℓ⁵`.
- **No per-curve index-calculus detector yet.**  The traits it records
  (singularity, generator validity, `a = −3` model, coefficient sizes, the
  quadratic-residue prefix statistic, volcano level, endomorphism ring
  where provable) do not measure ECDLP cost.  Finding a *weaker* curve
  needs that detector; see [What to do with it](#what-to-do-with-it).
- **The taskq path is untested live.**  The job specs pass taskq's own
  validator, but no worker has run them yet.

## Quick start

```bash
git clone https://github.com/aburan28/crypto && cd crypto
cargo build --release --bin isogeny_walk
B=./target/release/isogeny_walk

$B class  --curve p256 --max-ell 61          # what the class looks like (seconds)
$B walk   --curve p256 --max-ell 61 --max-curves 2000 --out runs/p256-2k
$B verify --curve p256 --dir runs/p256-2k     # independent replay
```

`walk` writes three files into `--out`:

| file | contents |
|:--|:--|
| `curves.yaml` | every curve found, in the [`docs/curves/ic/curves.yaml`](../curves/ic/curves.yaml) format: ICV1 slug, EC1 id and UID, model, generator, traits with statuses, proved volcano levels, `IW1` routes |
| `isogeny_routes.json` | every curve node, every verified edge with its kernel polynomial, every `IW1` route |
| `walk.json` | the run id, configuration, class invariants, class audits, counts and trait distributions |

## Curves and primes

| `--curve` | root |
|:--|:--|
| `p256` | P-256 (registered) |
| `p224` | P-224 (registered) |
| `p192` | P-192 (registered) |
| `custom` | your curve: `--p --a --b --order --gx --gy [--cofactor H] [--name NAME]` (decimal or `0x` hex) |

`--order` is the subgroup order `r` of the generator; the group order is
`r × cofactor`.  For a custom curve the order audit is a proof only when
`r` is prime and larger than `4√p`.

Choose the primes `ℓ` with `--max-ell N` (every odd prime up to N) or
`--primes 13,23,37`.

- **Atkin primes are skipped automatically** (`walk.json` lists them).
  For those `ℓ` the curve has no rational `ℓ`-isogeny.
- **Run `class` first** to see which `ℓ` split, which ramify and which are
  Atkin, plus each volcano's depth:

  | | splitting `ℓ` (2 neighbours each) | ramified `ℓ` (1 neighbour each) |
  |:--|:--|:--|
  | P-256 | 11, 13, 17, … | 3, 5 |
  | P-224 | 11, 47, 59, 61 | 29; and 3, whose volcano has depth 1 |
  | P-192 | 13, 23, 37, 43 | 5, 11, 31 |

## `walk` options

| option | meaning |
|:--|:--|
| `--max-curves N` | stop expanding once N curves are known (default 256) |
| `--max-depth D` | do not expand curves D steps from the root |
| `--max-ell N` / `--primes …` | which `ℓ` to walk |
| `--threads N` | cap worker threads (default: all cores) |
| `--class-audit-sample K` | re-run the class audits on K walked curves (default 8) |
| `--no-class-audits` | skip the class audits |
| `--store s3://bucket/prefix` | publish the outputs to S3 (below) |
| `--prune-local` | after a successful publish, delete the local `curves.yaml` and `isogeny_routes.json` |
| `--include-atkin` | also try Atkin primes (they find nothing) |

## Compact P-256 grids

The generic `isogeny_walk walk` output keeps the explored multigraph and a
complete root route for every curve.  For a million-curve screening prefix,
use the specialised streaming certificate instead:

```bash
cargo build --release --bin p256_isogeny_million
M=./target/release/p256_isogeny_million
C=$(git rev-parse HEAD)

$M --threads 4 generate --side 1000 --source-commit "$C" \
  --output runs/p256-grid-1m.jsonl.gz > runs/GENERATE.json
$M --threads 4 verify --input runs/p256-grid-1m.jsonl.gz \
  --audit-points 2 --audit-seed-x 7 > runs/VERIFY.json

# After the receipts agree, derive deterministic on-demand cache slices.
$M stage-cache --input runs/p256-grid-1m.jsonl.gz \
  --generation-receipt runs/GENERATE.json \
  --verification-receipt runs/VERIFY.json \
  --run p256-grid-1m-20261006 --output runs/cache

# Uses only the AWS CLI's ambient worker identity. The exact certificate,
# receipts, and 127 cache parts are content-addressed; complete.json is last.
$M publish --store s3://BUCKET/p256-isogeny-grid --cache runs/cache \
  --input runs/p256-grid-1m.jsonl.gz \
  --generation-receipt runs/GENERATE.json \
  --verification-receipt runs/VERIFY.json --scratch runs/s3-scratch

# Hydrate the full canonical evidence or one low-latency eight-row slice.
$M fetch --store s3://BUCKET/p256-isogeny-grid \
  --run p256-grid-1m-20261006 --output cache/p256-grid-1m.jsonl.gz \
  --scratch cache/s3-scratch
$M fetch --store s3://BUCKET/p256-isogeny-grid \
  --run p256-grid-1m-20261006 --part part-0002-rows-0008-0015 \
  --output cache/rows-0008-0015.jsonl.gz --scratch cache/s3-scratch
```

The fixed 1,000 by 1,000 grid has a degree-13 spine and degree-11 rows.  It
records exactly one explicit kernel-certified parent isogeny for every
non-root curve, streams deterministic gzip JSON Lines, and independently
replays root choices, models, generators, identities, detector values, order
proofs, kernels and target isomorphisms.  Any repeated j-invariant aborts the
run; it is not replaced by an adaptively selected curve.  This compact grid is
a registered P-256 screening experiment, not a replacement for the generic
multigraph walker and not evidence of an ECDLP speedup.  Its frozen protocol is
[`research/p256_isogeny_million_20261006/PROTOCOL.md`](../../research/p256_isogeny_million_20261006/PROTOCOL.md).

The cache keeps the canonical gzip as the archival byte sequence and derives a
preamble, 125 eight-row slices, and a summary member.  Concatenating the 127
members after decompression reproduces the canonical uncompressed bytes.  S3
objects live at `objects/sha256/<first-two>/<digest>` and are uploaded with a
declared SHA-256 plus `If-None-Match: *`; `runs/<run>/complete.json` is the only
authoritative completion record and is created last.  Fetches use a partial
file, check both stored and decompressed SHA-256, then rename atomically.

`publish` can also commit the compact completion receipt through Cairn by
adding `--cairn URL --objective ID --identity FILE --cairn-state FILE` (or the
explicit Stage-0 `--submitter NAME`).  Run it again after the configured Cairn
epoch to reveal the persisted commitment.  Cairn carries a
`coordination-receipt` labelled `s3-content-address-only`; S3 remains the byte
authority, and the receipt is not an ECDLP result or a claim that Cairn replayed
the million curve certificates.

### Continuing the compact grid by global rows

Do not increase `--side` merely to continue the frozen million prefix. That
would regenerate and recount the complete square. A strip records a disjoint
half-open interval of global degree-13 spine rows while keeping the same
degree-11 row width:

```bash
$M --threads 4 generate-strip --width 1000 --y-start 1000 --height 64 \
  --source-commit "$C" --output runs/p256-strip-y1000-h64.jsonl.gz \
  > runs/STRIP_GENERATE.json
$M --threads 4 verify-strip --input runs/p256-strip-y1000-h64.jsonl.gz \
  --audit-points 2 --audit-seed-x 7 > runs/STRIP_VERIFY.json

# This checks both record chains and every full-width j value. It rejects an
# overlap instead of silently subtracting duplicates.
$M audit-j-union runs/p256-grid-1m.jsonl.gz \
  runs/p256-strip-y1000-h64.jsonl.gz > runs/J_UNION.json
```

For a nonzero `--y-start`, the certificate binds one preceding spine curve as
context and does not count it. The first emitted spine curve still carries its
explicit degree-13 kernel certificate; every horizontal curve carries its
degree-11 certificate. `verify-strip` independently reconstructs the boundary
from P-256 before replaying the emitted interval. `audit-j-union` is an exact
identity-accounting gate over already replayed certificates, not a substitute
for mathematical replay.

Adjacent row coordinates do not by themselves prove distinct curves. Report
cumulative coverage only after the union audit returns a unique count equal to
the sum of its inputs. A coordinate collision is retained as a discrepancy and
is never replaced adaptively. The frozen continuation protocol is
[`research/p256_j_windows_20261007/PROTOCOL.md`](../../research/p256_j_windows_20261007/PROTOCOL.md).

The executed 64-row continuation passed generation, independent replay and an
exact union audit: the prior million plus 64,000 new curves yielded 1,064,000
distinct full-width j-invariants. The
[result report](../../research/p256_j_windows_20261007/RESULTS.md) preserves
the receipts, transfer boundary, visuals and the still-unset ECDLP speedup.

### Structural trait census over certified grids

Use `trait-census` only after each source certificate has passed its matching
independent replay. The command re-audits source hashes, record chains, exact j
uniqueness, and frozen P-256 headers, then emits one compact trait record per
curve. `verify-trait-census` reopens those same sources, reconstructs the
entire uncompressed artifact, and byte-compares it with the recorded census.

```bash
$M trait-census \
  --source runs/p256-grid-1m.jsonl.gz \
  --source runs/p256-strip-y1000-h64.jsonl.gz \
  --source-commit "$C" \
  --output runs/p256-traits-1064k.jsonl.gz \
  > runs/TRAIT_GENERATE.json

$M verify-trait-census \
  --input runs/p256-traits-1064k.jsonl.gz \
  --source runs/p256-grid-1m.jsonl.gz \
  --source runs/p256-strip-y1000-h64.jsonl.gz \
  > runs/TRAIT_VERIFY.json
```

Both receipts include SHA-256 digests and byte counts for the stored gzip and
its decompressed content. The census stores class-wide arithmetic once and
curve-varying model, identity, detector, path, and incoming-edge traits once
per curve. An exact path degree is represented losslessly as
`11^x * 13^y`, together with its bit length and decimal-digit count.

The degree-11 **row** and degree-13 **spine** labels describe construction
axes. They are not isogeny-volcano directions. The separate
`volcano_direction` field is populated only from the Frobenius-discriminant,
endomorphism-conductor, and splitting proof; otherwise it remains typed as
unknown with a reason. A structural census does not measure discrete-log cost,
so its `ecdlp_speedup` remains null until a separate preregistered end-to-end
solver experiment supplies that evidence. The frozen census protocol is
[`research/p256_isogeny_traits_20261008/PROTOCOL.md`](../../research/p256_isogeny_traits_20261008/PROTOCOL.md).

## Sizing

Measured on a 14-core Apple M4 Pro, unisolated.  Treat these as estimates
for your machine.

| quantity | measured |
|:--|:--|
| peak memory | ≈ 87 KB per curve: 433 MB for 5,000 curves |
| raw output | ≈ 17 KB per curve: about 340–410 MB per 20,000 |
| stored in S3 (gzip) | about ¼ of raw |
| time, 20,000 curves, all `ℓ ≤ 61` | P-256 ≈ 2.5 min, P-192 ≈ 3.5 min, P-224 ≈ 6.5 min |
| `Φ_ℓ` startup | under 1 s for `ℓ ≤ 31`, ≈ 10 s for `ℓ = 59` |

| `--max-curves` | memory | raw disk | fits on |
|:--|:--|:--|:--|
| 20,000 | ~2 GB | ~0.4 GB | any laptop |
| 100,000 | ~9 GB | ~2 GB | a 16 GB laptop |
| 300,000 | ~26 GB | ~6 GB | a 32–48 GB machine |

Beyond that, use more machines on different roots or `ℓ` sets for full
multigraph output.  The compact P-256 grid above is streaming, but deliberately
retains only one parent edge per curve and has a different coverage contract.

## Running offline

Only S3 needs a network.

```bash
cargo fetch                                   # once, online
cargo build --release --offline --bin isogeny_walk
$B walk --curve p192 --max-ell 61 --max-curves 20000 --out runs/p192-20k
# … later, online:
$B publish --dir runs/p192-20k --store s3://crypto-autoresearcher/isogeny-walk
```

`publish` uploads the walk under the same run id `--store` would have used
(`walk.json#run.run_id`).  It refuses files whose hashes no longer match
`walk.json`, and it uploads nothing if that run is already stored.

## Storing in S3

### Point it at a bucket

Any S3 bucket works.  The team store is
`s3://crypto-autoresearcher/isogeny-walk` (account `590183823895`,
`us-west-2`; SSE-S3 encrypted, public access blocked).

The walker calls the **AWS CLI v2** (`aws s3api`) with your normal
credentials, so set them the usual way:

```bash
aws configure                      # or: export AWS_PROFILE=myprofile
export AWS_REGION=us-west-2        # the bucket's region
aws sts get-caller-identity        # check you are who you think you are
```

Use CLI v2.27 or newer: conditional writes (`--if-none-match`) need a recent
CLI.  The IAM permissions needed on the prefix are:

```json
{
  "Version": "2012-10-17",
  "Statement": [
    {"Effect": "Allow", "Action": ["s3:GetObject", "s3:PutObject"],
     "Resource": "arn:aws:s3:::crypto-autoresearcher/isogeny-walk/*"},
    {"Effect": "Allow", "Action": "s3:ListBucket",
     "Resource": "arn:aws:s3:::crypto-autoresearcher",
     "Condition": {"StringLike": {"s3:prefix": "isogeny-walk/*"}}}
  ]
}
```

`s3:ListBucket` matters.  Without it, S3 answers "access denied" instead of
"not found" for a missing object, and the walker cannot tell whether a run
already exists.

For your own bucket, substitute `s3://your-bucket/your-prefix`.  Recommended
settings are encryption on and public access blocked.  Versioning is
optional, because every object is written once.

### Commands

```bash
S=s3://crypto-autoresearcher/isogeny-walk

# Walk and publish; keep only the small files locally.
$B walk --curve p256 --max-ell 61 --max-curves 20000 --store $S --prune-local --out runs/p256-20k
cat runs/p256-20k/STORE.json                    # run id and marker URI

# Fetch a stored walk anywhere (every hash re-checked), then use it.
$B fetch  --from $S --run p256-aed06c3fc90f7645 --out w
$B verify --curve p256 --dir w

# Trait shards from S3 and back to S3, then merge from S3.
for i in $(seq 0 7); do
  $B traits --curve p256 --dir w --shard $i --of 8 --out t$i --store $S --run p256-aed06c3fc90f7645
done
$B collect --from $S --run p256-aed06c3fc90f7645 --of 8 --store $S --out traits
$B fetch --from $S --run p256-aed06c3fc90f7645 --collected --of 8 --out traits-again
```

### What is in the bucket

```text
<prefix>/runs/<curve>-<walk key>/
  attempts/<attempt id>/curves.yaml.gz, isogeny_routes.json.gz, walk.json.gz
  complete.json                                  ← the authoritative record
  traits/<N>/shard-<i>/attempts/…, complete.json
  traits/<N>/collected/attempts/…, complete.json
```

- **The run id** is the curve plus a hash of the walk's arguments.  The
  same arguments always give the same run id.
- **`complete.json`** lists every object by key, stored SHA-256, size and
  uncompressed SHA-256.  It is created only after every upload was checked
  against S3's own SHA-256.
- **Write-once.** Every object is created with `If-None-Match: *`, so
  nothing is ever overwritten.  If two machines publish the same run, one
  marker wins and the other run prints
  `another attempt already completed`.  That is not an error; the run is
  stored.

Runs stored so far:

| curve | run id | curves |
|:--|:--|--:|
| P-256 | `p256-aed06c3fc90f7645` | 20,000 |
| P-224 | `p224-03abf7dcdecb7b0a` | 20,000 |
| P-192 | `p192-71c04205135fcf8f` | 20,000 |

## Queued trait jobs (taskq)

```bash
$B plan --curve p256 --max-ell 61 --max-curves 20000 \
    --commit "$(git rev-parse HEAD)" --shards 16 \
    --store s3://crypto-autoresearcher/isogeny-walk --walk-from-store --out specs/
for f in specs/*.json; do taskq submit --spec "$f"; done
$B collect --from s3://crypto-autoresearcher/isogeny-walk --run <run id> --of 16 --out traits
```

- **The commit** must be pushed and must contain the walker.
- **`--walk-from-store`** makes each job fetch the stored walk instead of
  rebuilding it.  Publish the walk first, with `walk --store`.
- **Credentials.** The workers need AWS credentials for the prefix.

## Portable bounded P-256 task

The separate `p256_isogeny_task` binary packages PR #1351's exact
degree-11 P-256 prefix as one content-addressed job.  It is useful for testing
an untrusted worker boundary: the result is accepted only after every edge is
replayed and every digest is recomputed.

Run the local, no-network path with:

```bash
tools/run_p256_isogeny_task.sh ./runs/p256-prefix 4096 "$(git rev-parse HEAD)"
```

For a small integration check, replace `4096` with `2`.  The directory gets
the input task, compressed certificate, verified result and an offline Cairn
coordination receipt.  Nothing is submitted.  To emit a TaskQ spec from the
same task:

```bash
target/release/p256_isogeny_task plan-taskq \
  --task runs/p256-prefix/input-task.json --output p256-prefix-taskq.json
```

The binary also has `publish-cairn`, but it requires
`--ack-coordination-only`: the current receipt says only that its publisher
performed a local full replay.  It is not a self-contained Cairn checker and
must not be used on a paid objective.

This adapter is intentionally one semantic job.  The degree-11 rule produces
a sequential chain, so changing `run_id` or `task_id` leaves `work_sha256`
unchanged.  Racing several workers may reduce latency or add redundancy; it
does not enumerate more curves.  The frozen contract and the frontier-anchor
requirement for real fan-out are in
[`research/p256_isogeny_shared_task_20261005/PROTOCOL.md`](../../research/p256_isogeny_shared_task_20261005/PROTOCOL.md).

## Reading the results

- **Traits per curve** are in `traits.jsonl` (one JSON object per curve)
  or in each curve's `trait_status` in `curves.yaml`.
- **Distributions over the run** are in `collect.json#distributions`.
- **Class-wide audits** are in `walk.json#class_audits`: `ecc_safety`,
  the structural report and the PKM signals.  They are identical for every
  curve in a class; the `invariance_check` shows that on sampled curves.
- **Volcano positions** are in each curve's `position_alias`
  (`EC1…V3L1V11L0`), and directions are in the `IW1` ids (`E3d1` means
  down one level).

## What to do with it

The tool finds and certifies curves.  It does not yet say which one is
weaker.  The work in order:

1. **Walk the roots you care about** at 20k–100k curves, all `ℓ ≤ 61`, and
   keep everything in S3.  Three roots are already stored.
2. **Write a detector for the question.**  For index calculus that means a
   per-curve factor-base or decomposition statistic.  Add a `Detector` in
   `src/cryptanalysis/isogeny_walk/traits.rs`; every walk, shard and queued
   job then records it.  PR #1330's rule applies: freeze the statistic and
   its threshold *before* reading walk output.
3. **Screen, then measure.**  Run the detector over stored walks with
   `traits` shards.  Take only pre-registered finalists to the full
   one-target solver comparison against rho (`AGENTS.md` §8, §12).  A
   screen score is never a speedup.

## Troubleshooting

| message | cause |
|:--|:--|
| `cannot run the AWS CLI` | `aws` is not on `PATH` |
| `put … AccessDenied` / `get … AccessDenied` | credentials lack `s3:PutObject`/`s3:GetObject`, or `s3:ListBucket` for a missing key |
| `another attempt already completed` | that run is already stored; nothing to do |
| `… already exists; attempt ids must be fresh` | two publishes in the same nanosecond on one host; rerun |
| `walk.json has no run id` | the walk was made by a build older than `publish`; rerun the walk |
| `… does not match the hash walk.json records` | a file was edited after the walk; `publish` refuses it |
| `missing shards [...]` from `collect` | not every shard was run or published |
| slow start | building `Φ_ℓ` for large `ℓ`; use `--max-ell 31` |

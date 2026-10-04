# P-256 isogeny search: native S3 ledger and Cairn receipts

Status: **control-plane contract implemented; production search blocked**

Date frozen: 2026-10-04

Result class: **plan / stage diagnostic, not an ECDLP speedup**

## Question and falsification target

The research question is whether any curve in the `F_p`-isogeny class of
registered P-256 admits a lower **complete, verified one-target ECDLP cost**
than the matched P-256 baseline and Pollard-rho reference. The curve identity
is:

```text
icv1-fp256-t89188191154553853111372247798585809583-f188c491
```

The proposed search emits up to `2^40` walk records, applies a cheap structural
screen, then tests only independently held-out finalists with the whole solver.
A prefix-residue surplus, relation-yield improvement, or solver-stage
improvement is a stage diagnostic. It is not evidence that index calculus is
faster.

The hypothesis passes only if a finalist:

1. has a checked isogeny path from registered P-256;
2. solves the same planted one-target DLPs as the baseline and matched rho;
3. includes setup, relation collection, filtering, matrix work, individual
   logarithm, recovery and independent `[k]P = Q` verification;
4. has `speedup = baseline_total_operations / candidate_total_operations > 1`;
5. preserves that result on preregistered holdouts; and
6. for any wall-time claim, has a paired 95% confidence interval excluding no
   improvement on matched hardware.

Otherwise the result is negative or inconclusive. The screen is never itself
the success condition.

## Frozen scale and identity

| quantity | value | meaning |
|---|---:|---|
| shard bits | 20 | 1,048,576 independent shard ids |
| records per shard bits | 20 | 1,048,576 emitted records per committed shard |
| total | 1,099,511,627,776 | `2^40` emitted records |
| Batch array ceiling | 10,000 | at least 105 worker-array submissions |
| result class | `stage-diagnostic` | no speedup may be derived from this run |

`2^40` means **emitted records**, not `2^40` globally distinct curves. Every
committed shard must have no repeated `j`-invariant within that shard. Global
deduplication is a separate reduction. A production report must state all
three counts: emitted records, distinct `j`-invariants, and distinct registered
curve models.

The native contract is
`src/cryptanalysis/p256_isogeny_campaign.rs`. It rejects any production
manifest that is not exactly `2^20 x 2^20`, does not name registered P-256, or
does not classify itself as a stage diagnostic.

## Architecture and authority

```text
AWS Batch static array index
        |
        v
native P-256 isogeny worker                 [not implemented]
        |
        | immutable attempt objects, SHA-256 metadata
        v
S3 attempts/<attempt>/{raw,stage1}.jsonl.zst
        |
        | conditional PUT, If-None-Match: *
        v
S3 shards/<shard>/complete.json             [execution authority]
        |
        +--> reducers read only winning markers
        |
        +--> native Cairn transport commit/reveal
             coordination receipt            [audit reference only]
```

S3 is the execution ledger. Cairn is the immutable coordination and audit
ledger. A Cairn receipt names content by hash; it does not prove that S3 still
serves the bytes, that the bytes contain the claimed records, or that those
records are in the P-256 isogeny class.

No scheduler races on an S3 queue. Task `i` of `N` owns shards
`i, i+N, i+2N, ...`. A retry may start again, but only one attempt can create
the shard marker. A `409` or `412` from the conditional marker write means the
other attempt won; it is not a failed scientific run.

## S3 transaction

For run `p256-2p40-v1`, shard 7 and attempt `batch-42-attempt-1`:

```text
p256-isogeny-search/runs/p256-2p40-v1/
  RUN.json
  shards/00000007/
    attempts/batch-42-attempt-1/raw.jsonl.zst
    attempts/batch-42-attempt-1/stage1.jsonl.zst
    complete.json
```

A writer must:

1. read and hash the immutable `RUN.json` and native-worker configuration;
2. compute its deterministic shard;
3. upload the two attempt objects under a unique attempt id;
4. read their metadata back and compare SHA-256 values;
5. build a Merkle root over canonical raw records;
6. validate the completion marker with `ShardComplete::validate`; and
7. create `complete.json` with `If-None-Match: *`.

Reducers trust only a valid marker and the exact objects it names. Garbage
collection may remove unreferenced attempts after a retention interval; it may
never remove an object referenced by a winning marker. Bucket versioning,
public-access blocking and encryption at rest are deployment requirements.

At the frozen layout the primary layer alone creates at least 3,145,728 S3
objects: two attempt objects and one marker for every shard, before retries and
reductions. If a canonical raw record occupies `b` bytes after compression,
raw storage is `2^40 * b`: 64 TiB at 64 bytes, 128 TiB at 128 bytes and 256 TiB
at 256 bytes. A measured native pilot must replace `b` before any budget is
approved.

## Cairn commit-reveal receipt

After winning the S3 marker race, the publisher creates
`p256.isogeny.cairn.shard-receipt/v1` with:

- registered curve and run identity;
- shard id and record count;
- configuration, marker, raw-object, stage-one and Merkle-root SHA-256 values;
- the exact S3 marker key;
- `claim_class = coordination-receipt`; and
- `verification = s3-content-address-only`.

The stable local idempotency key is `run_id:shard_id`, not an attempt id or an
object hash. An AWS retry therefore cannot become a second logical Cairn unit.
`CairnTransport::publish_artifact` signs and persists the commitment using the
existing Cairn transport; `reveal_pending` reveals it only after the commitment
epoch closes.

This artifact **must not be attached to a paid piecework objective yet**. A
shape checker would prove only that the receipt is well formed. Before payment
is sound, a pinned native checker or independently bonded audit must establish
at least:

1. exact marker bytes and object hashes;
2. all record Merkle leaves or a proof system with an explicit soundness
   bound;
3. every claimed isogeny edge and the fixed P-256 starting identity;
4. deterministic shard ownership and record count; and
5. stage-one recomputation without trusting worker-supplied scores.

Random spot checks alone do not certify a `2^20`-record shard. If sampling is
used as an economic audit rather than a proof, the sample rate, miss
probability, bond and slash conditions must be frozen in the objective and the
claim must remain labelled probabilistic.

## Screening preregistration

The earlier prototype considered the number of quadratic residues among the
first 64 `x` values. Under a null `Binomial(64, 1/2)`,
`P(score >= 50) = 3.5347438e-6`, or about 3.89 million rows over `2^40` emitted
records before deduplication. This is a sizing calculation, not a frozen
production threshold.

The native pilot must freeze before looking at candidate identities:

- discovery statistic and threshold;
- translated or pseudorandom control probes that are never used for selection;
- per-shard top-K and global retention policy;
- global deduplication key;
- holdout statistic and multiplicity correction; and
- the exact end-to-end solver protocol for finalists.

Changing a threshold after inspecting the search population starts a new run
id and invalidates selection claims from the old run.

## Gates before scale-up

### Gate 0 — control plane (this change)

- Native manifest, assignment, S3 marker and Cairn receipt validation.
- Cairn generic artifact path uses its existing signed, persisted
  commit-reveal outbox.
- No AWS resource is created and no objective is posted.

### Gate 1 — native isogeny engine

- Implement P-256 field arithmetic and the selected small-degree isogeny
  action in Rust; no Python or Sage execution path.
- Reproduce frozen independent small-field and P-256 fixtures.
- Verify every emitted edge independently from the worker's choice logic.
- Establish that the walk preserves the P-256 trace and registered identity.

Until Gate 1 passes, no production worker image exists.

### Gate 2 — bounded pilot

- Run the same image and transaction protocol at `2^24` emitted records.
- Exercise simultaneous retries and prove one winning marker per shard.
- Measure curves/second, bytes/record, S3 request counts and reducer memory.
- Set explicit compute, storage and request-cost caps from those measurements.
- Independently replay a preregistered sample on a second implementation.

### Gate 3 — solver validation

- Run complete matched baseline/candidate/rho solves on toy and scaling
  analogues, then on feasible registered-curve workloads.
- Price every phase in one operation unit and publish the boundary/ratio table.
- Preserve failures, timeouts and out-of-memory outcomes.

### Gate 4 — expansion

Expansion is incremental and stops automatically if a cost cap is reached, a
correctness audit fails, global uniqueness is materially below the budget
model, the discovery effect disappears on controls/holdouts, or no retained
candidate improves complete solver cost. Passing the storage pilot alone does
not authorize `2^40` compute.

## Current conclusion

The S3/Cairn control contract is implementable and tested. The cryptographic
worker and a sound payable-work witness are not. Therefore the rigorous next
step is the native bounded pilot, not a `2^40` launch, and no index-calculus
speedup is established by this protocol.

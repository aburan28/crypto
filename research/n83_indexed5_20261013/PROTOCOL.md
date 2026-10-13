# N83 five-byte exact root index: frozen public-target gate

## Contribution and selection question

`KIC_INDEXED_ROOT5=1` stores a ten-bit key fingerprint and a thirty-bit
state/root reference in each five-byte slot. Every candidate match checks
the complete canonical 83-bit root in the retained 28-byte state before
returning a witness. Insertion retains the first witness for a repeated
root. The regular two-key request reaches `2^31` slots, or 10 GiB, for
the K2202–K2543 source-count maxima.

The first correctness gate is K512 on the one frozen public target. Its
exact base hash, state count, four-point relation, x codes, target probes
and recovered scalar must match the archived K512 packed-state arm. A
verified K1798 comparison is the next local gate if the host resource
preflight passes. K2500 is the larger-base primary candidate, intended for
an independently receipted Tier 1 cloud instance with enough memory.
These arms are run separately, each with a write-once directory. A
resource stop or timeout remains a censored row.

## Frozen instance and inputs

- Koblitz `a=b=1` over `F_2^83`, polynomial low terms `[0,2,4,7]`,
  subgroup order `8569786107849059`, cofactor `1128547018`, generator
  `(6614244651398184860292015,8868260410029249391539529)`.
- One previously unseen public point in `primary_target.jsonl`, SHA-256
  `21c8dfc7d488909e5c4a42c839fcc01bdefa9821c11ba001626a1d15697240e8`.
  The known scalar is in the separate `answer.jsonl` checker sidecar and
  is never passed to the producer.
- Deterministic prefix base `construct:83:1:K`, cofactor projection,
  signed-Frobenius closure, rank seed 29, twelve rank threads, one target
  thread. Bloom and sparse-canonical overrides are disabled. The primary
  online clock begins after reusable base/index/log preparation and ends
  after scalar replay. Setup, rank and final matrix algebra are separate.
- Same-point rho reference:
  `../n83_dense_root_20261013/runs/rho_primary/receipt.json`.
  Local CPU ratios are exploratory. A controlled Tier 1 claim requires
  the calibrated cloud-instance receipt; Tier 0 is physical-host isolation.

## Frozen code and limits

- Indexed-slot source commit: `3da76633c4661f9bc43fa23d662c7a8adb57e0a1`.
- Producer source SHA-256:
  `42f21942e5e85874dfe9017f720697ab07f4cb0546a5c7d27e448a5a61fc448b`.
- Producer release binary SHA-256:
  `7a0946cd6bfb419a4d3fa60c6d8dfd4a72197021f6c0d60e8339390d5f54136d`.
- Direct-child wrapper SHA-256:
  `1b1b175e9e06e1b27943a213b419b751ebaed109a553961529a644aea4907ec0`.
- Process caps: K512 900 s; K1798 1800 s; K2500 2400 s. Preserve
  wall time, child CPU, peak resident bytes, faults, context switches,
  status, exact summary/target rows and full point-defined base dump.
- K1798 local start requires at least 10 GiB free on the macOS data
  volume and no active memory-pressure guard. K2500 requires a separate
  Tier 1/0 host because the prior shared-host K2202 run exhausted swap.

`run.sh` verifies the source, binary, wrapper and public-input hashes
before it creates any run directory. The [independent checker](analyze.py)
replays the public target through generic binary-curve arithmetic and
checks every archived base orbit and the four-point relation.

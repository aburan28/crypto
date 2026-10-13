# N83 compact-state factor-base gate

## Contribution and question

The opt-in `KIC_STATE28=1` layout stores each N83 four-point index state
in 28 bytes, including both 83-bit roots, instead of the 48-byte wide
layout. The frozen test asks whether this makes the K2202 prefix base
complete with the 75% packed root-table request and whether its verified
single-target probe count improves on the K1798 dense control.

## Frozen instance and candidate arms

- Koblitz `a=b=1` over `F_2^83`, polynomial low terms `[0,2,4,7]`,
  subgroup order `8569786107849059`, cofactor `1128547018`, generator
  `(6614244651398184860292015,8868260410029249391539529)`.
- The public target is copied byte-for-byte from the prior frozen 64-point
  holdout; its SHA-256 is `21c8dfc7d488909e5c4a42c839fcc01bdefa9821c11ba001626a1d15697240e8`.
  The answer sidecar is excluded from producer input and timed intervals.
- Correctness control: K512 with wide states, then K512 with 28-byte states,
  each using the dense packed table. Compare exact base hash, state/key
  counts, rank, first verified target scalar, and first-witness relation.
- Primary arms after the control: K2202 with 28-byte states and the 75%
  root-table request, then K1798 with 28-byte states and the dense request.
  Each gets a 1,800-second process cap. Each run uses rank seed 29,
  12 rank threads, one target thread, no bloom, and no sparse-canonical
  override. The root-table slot limit is 8 GiB. Preserve timeout and
  resource failures as run rows.
- Pair the successful one-target online interval to the same-point strong
  rho receipt in `../n83_dense_root_20261013/runs/rho_primary`. CPU timing
  ratios on this shared host are exploratory. Tier 0/1 promotion needs its
  separate host receipt.
- After a complete K2202 run, archive the exact point-defined factor base
  via `KIC_DUMP_BASE`, check point/subgroup/orbit labels independently,
  and form full candidate IDs from canonical manifests. Do not infer an
  actual base count or candidate ID from K alone.

## Frozen code and resource envelope

- State-layout commit: `0609be6650fa66249502d76f70caaf760fbcb5cd`.
- Producer source SHA-256:
  `5f8e8ec971e159322f54639c9d8415024e0197360e36b3272e6f9e3ad4628596`.
- Release producer SHA-256:
  `13c223fdff4b5341b71277bed68f2bcfaf5d0f21013220bb9c75866f1fff2e15`.
- Host: Apple M4 Pro, 48 GiB RAM, macOS 26.6, shared CPU. A direct-child
  usage wrapper records wall time, child CPU, peak RSS, faults, context
  switches, status and cap. The wrapper is the frozen one in the load75
  study, hash `1b1b175e9e06e1b27943a213b419b751ebaed109a553961529a644aea4907ec0`.
- `cargo test --release --locked --offline --example koblitz_orbit_dlp_fast`
  passed all 13 focused tests. Build and test logs are archived here.

`run.sh` checks all source, binary, wrapper and input hashes before each
write-once arm. The run's primary online clock starts with target-dependent
computation after base, index, rank and log preparation and ends after
scalar replay. Setup, relation rank and final linear algebra are separate
reported phases.

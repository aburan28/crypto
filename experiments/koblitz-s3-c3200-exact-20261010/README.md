# C3200 exact quotient-root index with adaptive capacity

The C3200 implementation allocates its initial quotient-root table from one
expected distinct canonical key per regular state, while retaining exact
per-shard growth and first-witness reconstruction. This changes the initial
table from 16 GiB to **8 GiB** at C3200. The packed slot assigns 30 bits to
the state ordinal and 11 bits to linear-probe displacement; the remainder and
shard bits still reconstruct the complete 53-bit mixed key for exact checks.
The 16 release tests include a high-bit state witness in this slot layout and
the existing growth, collision, reconstruction, group-lift, and target-order
checks.

The earlier C2800 run inserted two raw roots per state but retained
415,077,662 distinct keys from 415,520,000 states. The C3200 initial
allocation uses 542,720,000 states / 64 shards, giving 16,777,216 slots
per shard and 1,073,741,824 total. If a shard has more distinct keys than
estimated, the table grows exactly. The frozen feasibility gate checks the
observed slot count and 32 GiB peak RSS, so a growth step remains a recorded
resource result.

| C3200 fixed payload | Bytes | GiB |
| --- | ---: | ---: |
| 542,720,000 states × 40 bytes | 21,708,800,000 | 20.22 |
| 1,073,741,824 root slots × 8 bytes | 8,589,934,592 | 8.00 |
| Combined | 30,298,734,592 | 28.22 |

`FEASIBILITY_FREEZE.json` binds the exact source and binary, N53 T025 public
point, relation seed, 14 preparation workers, one target worker, and 1200 s
limit before execution. `run_feasibility.py` verifies all hashes and records
stdout, report, stderr, child resource use, rank, layout, 32 GiB memory gate,
all five online phases, scalar replay, and group verification. The supplied
point is the target input; fixture construction remains outside the timing.
Local CPU wall time is exploratory until a host-level isolation receipt is
available.

## Verified T025 feasibility result

The frozen run completed in 177.7 s including reusable setup. It constructed
542,720,000 regular states and retained 541,971,126 distinct canonical keys
in exactly 1,073,741,824 slots. The factor base has 339,200 distinct usable
points and 3200 effective columns. Preparation reached full rank from 6300
verified relation attempts. The recovered scalar `429728658711` passed both
the native group check and a separate replay in the repository's checked
Sage runtime.

| Gate | Recorded result |
| --- | ---: |
| Solver exit and scalar replay | pass |
| Full rank before target input | 3200 |
| 30-bit state ordinal, 11-bit displacement, exact table size | pass |
| Peak RSS / 32 GiB cap | 23.45 GiB / pass |
| Charged T025 online interval | 5.551 ms, exploratory |
| Target S3 calls | 498 |

The T025 control produced the same first relation witness and 498 S3 calls
as C2800 and C3182. `PROMOTION.json` binds the raw receipt to candidate
`IC1N53Ckb1fb339200PDP4rootRCsampleLAgaussTDdirectISO0hdf57773efaf3`
and the T025 workload. The canonical candidate manifest, compressed report,
stdout, stderr, native receipt, Sage runtime receipt, and independent replay
are preserved here. `independent_replay.py` reconstructs the N53 polynomial
field and binary curve from the manifest and verifies the public point using
Sage elliptic-curve arithmetic.

The next measurement pairs fresh holdout points against C1200 and C2800 on
the same source and resource envelope, preserving all one-target failures and
phase costs. The local host has no exclusive CPU partition receipt, so its
timing ratios remain exploratory.

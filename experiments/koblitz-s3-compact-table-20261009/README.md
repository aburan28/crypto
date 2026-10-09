# Compact S3 root-table values at 396 factor-base columns

Packing each indexed pair state into 32 bits and storing keys and values in
separate arrays reduced the 396-column N53 root-table payload from
536,870,912 to 402,653,184 bytes. Exact-key lookup, first-inserted witness
choice, rank preparation, and scalar recovery were preserved. Six release
tests passed; one-worker N41 and N53 controls produced identical first
relations, S3-call counts, and prepared witnesses. On eight frozen N53
points, all 32 paired solves were verified by native code and independently
replayed by checked Sage, including 588 distinct relation witnesses.

## Exact candidates and method

Both candidates use curve `EC1N53Ckb1h738a572b91f2` over `GF(2^53)` with
modulus `x^53+x^6+x^2+x+1`, equation `y²+xy=x³+1`, and subgroup order
21044858204113. The same natural low-`x` factor-base construction yields
41,976 distinct subgroup-usable points and 396 folded columns. Both arms
request 395 target-independent relation rows; the fixed 14-row preparation
batch reaches rank 396. Each target starts 14 rotating first-winner workers
from the same index-state sequence. The [freeze](FREEZE.json) and
[candidate manifests](identities/candidates/) identify the sources, binary
builds, table layout, factor-base digest, target workloads, and flags.

| Arm | Candidate ID | Root-table slot layout | Payload for 33,554,432 slots |
| --- | --- | --- | ---: |
| Existing | `IC1N53Ckb1fb41976PDP4rootRCsampleLAgaussTDdirectISO0h1e2914fbd8c2` | `(u64 key, u64 value)` | 536.9 MB |
| Compact | `IC1N53Ckb1fb41976PDP4rootRCsampleLAgaussTDdirectISO0hc2162a1f0c52` | `u64` keys, separate `u32` values | 402.7 MB |

The 32-bit value assigns ten bits each to the two representative indices
and six bits each to their relative and absolute Frobenius shifts. The
constructor checks those bounds. A vacant key is still represented by the
same `u64::MAX` sentinel; duplicate keys retain the first inserted value.
This keeps the compact table's decoded witnesses equal to the existing
table's for a fixed query order. An exact synthetic collision and boundary
packing test covers this storage contract.

## Frozen one-target panel

The [parent plan](../koblitz-s3-crosscurve-20261009/FREEZE.json) fixes
T025–T032. Every point was solved twice per arm, alternating ABBA and BAAB
order. The online timer begins after reusable factor-base, index, and rank
preparation, then charges query generation, all target S3 work, relation
checking, descent, and scalar replay. The [raw ledger](local/raw_paired32.jsonl.gz)
retains each run, including phase costs, operation counts, witness traces,
and status. The [analysis](local/analysis.json) checks exact identities and
phase sums. Checked [Sage replay](local/sage_replay.json) with its
[runtime receipt](local/sage_runtime_info.json) verifies all 41,976 base
point labels, 588 distinct group equations, the prepared matrix rank, and
all recovered scalars.

| Measure | Existing | Compact |
| --- | ---: | ---: |
| Native and Sage verified runs | 16/16 | 16/16 |
| Median online ms | 3.750 | 5.704 |
| Median charged target PDP ms | 3.487 | 4.654 |
| Median target S3 calls | 9,346 | 4,784 |
| Median reusable setup ms, supplementary | 4,493 | 4,548 |
| A/A p95 absolute log online ratio | 0.751 | 2.427 |

The exploratory existing/compact online geometric ratio was **0.855×**
(target bootstrap 95% interval `[0.527, 1.437]`), with two of eight points
at least 2×. The existing/compact S3-call ratio was 1.791×; first-winner
thread scheduling changes the amount of completed work between invocations.
The A/A values exceed the 10% host-noise gate, and the host lacks an
exclusive-core isolation receipt, so controlled wall-time speedup remains
unknown. This table layout saves 25% of table payload but does not advance
the complete online result. A packed 12-byte inline slot is the next
layout to test because it retains a single slot read while reducing payload.

The compact source SHA-256 is
`e4d186a0e4d8d17986a4d76363a49a897a89718e0185d7b09fd5faa9dd51739b`;
its release binary SHA-256 is
`82687dbda39cb62e437610acb9f6bbb1e979aaa1b972080da66849ab9007f15f`.
The baseline release binary SHA-256 is
`8217647e609942bb1b2876cd1cd05dde4d6a7dd2b6f4d3713421607d2d3a3920`.
The plan SHA-256 is
`e4aadd161bdc0a1eb16fe6803d132c64fcb31f44615a49b4739da696e1af897e`,
and the uncompressed ledger SHA-256 is
`d833ca5ef4ac6ff3c588c9a7ecfd794470fa073dbb7d621f7f2882a288727cec`.

To reproduce from the exact source tree and rebuild both release examples,
use `CARGO_TARGET_DIR=/private/tmp/s3-prepared-target-20261009` with
`cargo build --release --offline --example s3_partial_prepared_parallel_targetspan
--example s3_partial_prepared_compact_table_targetspan`. Then run
`run_paired.py --first 25 --count 8 --output <new-directory>
--binaries /private/tmp/s3-prepared-target-20261009/release/examples`.
The runner checks both binary hashes and the code tree before every panel.
Use `analyze_paired.py --ledger <new-directory>/raw.jsonl --output
<new-directory>/analysis.json` for the paired result. For independent
replay, save `/Volumes/SSD990/cryptanalysis/sage --runtime-info` before
launching `/Volumes/SSD990/cryptanalysis/sage -python replay_sage.py` with
the new ledger, runtime receipt, and a new output receipt path.

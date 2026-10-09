# Eight-state batch inversion for S3 target extraction

The target extractor now shares one normal-basis field inversion across up
to eight ordinary indexed S3 states. Each state still produces its two
ordered S3 answers and undergoes exact root-table and group checks before
claiming a relation. Reusable relation preparation retains the existing
extractor. Six release tests passed, including batch-versus-scalar inverse
equality on N13, N41, and N53. One-worker N41 and N53 controls matched the
baseline's prepared witnesses, target relation, state probes, S3-call count,
and scalar exactly. The frozen N53 panel completed 32/32 verified native
single-target solves. Checked [Sage replay](local/sage_replay.json), with its
[runtime receipt](local/sage_runtime_info.json), independently verified all
scalars, 41,976 factor-base point labels, and 588 distinct relation
witnesses.

## Exact candidates and clock

The [freeze](FREEZE.json) fixes curve `EC1N53Ckb1h738a572b91f2` over
`GF(2^53)` with modulus `x^53+x^6+x^2+x+1`, curve `y²+xy=x³+1`, subgroup
order 21044858204113, and public targets T025–T032. The low-`x` factor
base has 41,976 subgroup-usable points and 396 signed-Frobenius columns.
The target-independent batch requests rank 395 and reaches rank 396.
Both arms use 14 rotating target workers, first verified relation,
cooperative cancellation, and worker drain. The online clock charges query
generation, all target S3 attempts (including batch lookahead), relation
checking, descent, and scalar replay. Factor-base, root-index, and rank
preparation are reusable setup outside this interval.

| Arm | Candidate ID | Target-state inversion | Root-table payload |
| --- | --- | --- | ---: |
| Existing | `IC1N53Ckb1fb41976PDP4rootRCsampleLAgaussTDdirectISO0h1e2914fbd8c2` | one per indexed pair | 536,870,912 bytes |
| Batch eight | `IC1N53Ckb1fb41976PDP4rootRCsampleLAgaussTDdirectISO0hff55c3b35e32` | one per eight regular pair products | 536,870,912 bytes |

The [candidate manifests](identities/candidates/) fix the source and binary
hashes, operation policy, and curve/base identities. The batch uses
Montgomery simultaneous inversion with the existing Itoh-Tsujii kernel;
zero denominators retain the scalar exceptional path. A cancellation may
stop evaluation during a prepared batch, and its preparation time remains
inside the target online interval.

## Paired result

Each point was solved twice per arm in alternating ABBA and BAAB order.
The [raw ledger](local/raw_paired32.jsonl.gz) retains every phase, operation
count, witness, and status. The [analysis](local/analysis.json) checks exact
identities, target pairing, scalar replay, and exclusive phase sums.

| Measure | Existing | Batch eight |
| --- | ---: | ---: |
| Native verified runs | 16/16 | 16/16 |
| Median target S3 calls | 7,413 | 5,413 |
| Median charged target PDP | 1.910 ms | 2.611 ms |
| Median online time | 2.250 ms | 2.863 ms |
| Median reusable setup, supplementary | 3,201 ms | 3,412 ms |
| A/A p95 absolute log online ratio | 2.620 | 1.158 |

The exploratory existing/batch-eight online geometric ratio was **1.161×**
with target bootstrap 95% interval `[0.796, 1.712]`. Four of eight points
favored batching and two exceeded 2×. First-winner worker scheduling changes
the amount of completed S3 work between invocations. A/A repeat noise
exceeds the 10% host-noise gate, and the host has no exclusive-core
isolation receipt, so controlled speedup remains unknown. This panel does
not establish the requested 2× online gain. The next candidate moves
Frobenius orientation to a precomputed target conjugate, eliminating four
normal-to-polynomial conversions per target state.

The batch-eight source SHA-256 is
`2a3ecdd3309bda625d8528994b3e3b58ccbbdfeeb7c04c94087212529a7ae100`;
its release binary SHA-256 is
`b81249301e11f817ef93c628bba77e5ab48a96913411146dddd20216eb25b7f2`.
The baseline binary SHA-256 is
`8217647e609942bb1b2876cd1cd05dde4d6a7dd2b6f4d3713421607d2d3a3920`.
The plan SHA-256 is
`353ace203fe94e2ab446148748e840f35acb0bb156b079f5a9dd05ca403bb95c`,
and the uncompressed ledger SHA-256 is
`b71f6a937273efec1a602e6ec7c1896c7cfe2448a0d5bb9e59fce36327faf4a8`.

To reproduce, build both release examples using
`CARGO_TARGET_DIR=/private/tmp/s3-prepared-target-20261009 cargo build
--release --offline --example s3_partial_prepared_parallel_targetspan
--example s3_partial_prepared_batch8_targetspan`. Run
`run_paired.py --output <new-directory> --binaries
/private/tmp/s3-prepared-target-20261009/release/examples` and then
`analyze_paired.py --ledger <new-directory>/raw.jsonl --output
<new-directory>/analysis.json`. For independent replay, save
`/Volumes/SSD990/cryptanalysis/sage --runtime-info` before launching
`/Volumes/SSD990/cryptanalysis/sage -python replay_sage.py` with a new
receipt path.

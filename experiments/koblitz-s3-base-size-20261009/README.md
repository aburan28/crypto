# S3 factor-base size near a root-table capacity step

Increasing the N53 Frobenius-closed factor base from 244 to 396 folded
columns reduced target S3 calls by a **3.02× geometric mean** on eight
frozen public points while retaining a 33,554,432-slot root table. The
complete one-target online interval improved by an exploratory **1.56×**.
All 48 paired runs recovered and verified their scalars. Checked Sage
independently verified every factor-base orbit and label, 1,414 distinct
relation witnesses, all three prepared relation traces, and all 48 scalar
recoveries.

## Frozen candidates and method

The curve is `EC1N53Ckb1h738a572b91f2`, with `GF(2^53)` represented modulo
`x^53+x^6+x^2+x+1`, equation `y²+xy=x³+1`, and prime subgroup order
21044858204113. The exact 14-worker S3 binary, relation seed, target
query policy, and public point are paired. The variable is the number of
signed Frobenius orbit representatives obtained by the same natural low-`x`
scan. Each representative contributes 106 distinct subgroup-usable points.
The [freeze](FREEZE.json) and [candidate manifests](identities/candidates/)
fix the construction, actual base digests, rank thresholds, root-table
capacities, code identity, and workload IDs. The 244-column candidate's
manifest remains in the [parallel experiment](../koblitz-s3-parallel-target-20261009/identities/candidates/).

| Folded columns | Actual usable points `B` | Root-table slots | Candidate ID | Requested rank / actual rank |
| ---: | ---: | ---: | --- | ---: |
| 244 | 25,864 | 16,777,216 | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0he884cf865574` | 243 / 243 |
| 320 | 33,920 | 33,554,432 | `IC1N53Ckb1fb33920PDP4rootRCsampleLAgaussTDdirectISO0he065e8af00b1` | 319 / 319 |
| 396 | 41,976 | 33,554,432 | `IC1N53Ckb1fb41976PDP4rootRCsampleLAgaussTDdirectISO0h5be6b5627281` | 395 / 396 |

The 396-column preparation batch reached full rank before the target was
loaded, one row beyond its requested threshold. A separate construction
check at 400 columns allocated 67,108,864 slots. Selecting 396 columns
therefore kept nearly the same pair-sum coverage with half as many table
slots as that next allocation step. Setup and all target-independent rows
remain outside the primary online interval, with their costs recorded
separately.

## Paired one-target panel

The [parent plan](../koblitz-s3-crosscurve-20261009/FREEZE.json) fixes the
T025–T032 public points. Each point was solved twice with each candidate;
the three-arm order rotated by target and reversed on repetition. The
online clock starts after factor-base, root-index, and relation preparation,
then charges target query, all S3 decomposition work, relation checking,
target descent, and scalar replay. The [compressed raw ledger](local/raw_paired48.jsonl.gz)
preserves every invocation, status, rank witness, phase cost, root-call
count, and run ID. The [analysis](local/analysis.json) verifies point and
candidate identity, actual base size and digest, binary and source hashes,
and exclusive phase sums. The checked [Sage replay](local/sage_replay.json)
uses the saved [runtime receipt](local/sage_runtime_info.json) and both
new [factor-base dumps](local/).

| Measure | 244 columns | 320 columns | 396 columns |
| --- | ---: | ---: | ---: |
| Native and Sage verified runs | 16/16 | 16/16 | 16/16 |
| Median online ms | 6.692 | 4.430 | 3.302 |
| Median charged target PDP ms | 6.560 | 4.267 | 3.009 |
| Median target S3 calls | 29,661 | 16,409 | 7,208 |
| Median reusable setup ms, supplementary | 3,045 | 2,240 | 3,624 |
| A/A p95 absolute log online ratio | 1.473 | 2.196 | 1.707 |

| Paired target-level ratio | 244 / 320 | 244 / 396 |
| --- | ---: | ---: |
| Exploratory online geometric mean | 1.172× | 1.563× |
| Target bootstrap 95% interval | [0.742, 1.853]× | [1.036, 2.197]× |
| Points at least 2× | 2/8 | 3/8 |
| Target S3-call geometric mean | 1.972× | 3.018× |
| Reusable setup geometric ratio, supplementary | 1.377× | 0.884× |

The A/A values exceed the 10% noise gate (`log(1.10) = 0.0953`), and this
host has no exclusive-core isolation receipt. The wall ratios are
exploratory and the controlled speedup remains unknown. The 396-column
method clears 2× in target S3-call reduction; its complete online result
does not yet clear the frozen 2× timing threshold. The next experiment is
to reduce random root-table lookup cost at 396 columns or test whether a
smaller target-worker count improves that larger index under paired inputs.

The source SHA-256 is
`71f527cf7d4f401ad252a6f5ad3fa504c3f8becac9567e254e10d0b48cd75cfa`;
the release binary SHA-256 is
`ff090babc4477419dac6bac051670fd75e1e3832b34da1886839042a5d7c87c1`.
The plan SHA-256 is
`af593053ec19055ff6e72732e99eafc8f69e6be2da893487127dd6643bad1d2b`,
and the uncompressed ledger SHA-256 is
`bc3d9063af05d15d6c3a20d849618ec4fcaac864023299bcdac6aaa51bf3b590`.

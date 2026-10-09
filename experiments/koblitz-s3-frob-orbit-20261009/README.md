# Frobenius-orbit S3 target search

The target-query kernel can evaluate every Frobenius orientation of an indexed
point pair with one field inversion. It transforms the public target into its
inverse Frobenius conjugates, shares the indexed pair's polynomial coordinates,
and maps each returned root back to its absolute canonical shift. Seven
release example tests checked the ordered root and shift identities on both
curves. On the frozen eight-point panel, all 64 one-target runs recovered and
verified their scalars; checked Sage independently replayed all 64 and 854
distinct relation witnesses.

## Exact pipeline

Both curves use `y²+xy=x³+1`. N41 is
`EC1N41Ckb1hcbf4a37f6e66`, over `GF(2^41)` with modulus
`x^41+x^3+1`, subgroup order 549756390943, 20,008 usable factor-base
points, and 244 folded columns. N53 is `EC1N53Ckb1h738a572b91f2`,
over `GF(2^53)` with modulus `x^53+x^6+x^2+x+1`, subgroup order
21044858204113, 25,864 usable points, and 244 folded columns. The
[parent freeze](../koblitz-s3-crosscurve-20261009/FREEZE.json) fixes the
public targets and exact workload IDs. The [experiment freeze](FREEZE.json)
fixes the 243 prepared relation rows, source and binary hashes, workers,
target order, and candidate identities. `GF2_ARM_PMULL_FOLD` was unset for
both arms.

| Curve | Baseline candidate | All-orientation candidate | Target workers |
| --- | --- | --- | ---: |
| N41 | `IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0h0f4d1332b718` | `IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0h11dd73f3c165` | 1 |
| N53 | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0hf313888ef69e` | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0h0bf60c673463` | 14 |

For each index state the kernel computes `F^-k(Q)` for every `k`, then
evaluates the two ordered group-root queries against the original indexed
pair. Montgomery inversion shares one inverse across all nonzero denominator
pairs; zero denominators follow the scalar fallback. Every candidate root
is mapped through its canonical Frobenius shift before lookup, and a full
four-point group equation is checked before a relation is accepted. The
first verified relation cooperatively stops other target workers. Rank
preparation retains the baseline's 16-state batch kernel.

## Frozen single-target panel

Each of eight previously supplied public points per curve was solved twice
per arm in ABBA or BAAB order. One invocation contains one target. The online
interval begins after factor-base, root-index, and rank preparation, then
includes target query, all charged S3 attempts, relation check, target
descent, and scalar replay. Process launch, input load, and known-scalar
fixture construction are outside it. The [raw ledger](local/raw_paired64.jsonl.gz)
retains every invocation, phase, memory peak, and status. The
[analysis](local/analysis.json) checks candidate and run IDs, hashes,
target pairings, phase sums, and A/A stability. The [Sage replay](local/sage_replay.json)
and checked [runtime receipt](local/sage_runtime_info.json) independently
check the group equations, prepared rank, and recovered scalars.

| Measure | N41 | N53 |
| --- | ---: | ---: |
| Native and Sage verified runs | 32/32 | 32/32 |
| Baseline/orbit online geometric ratio | 1.222× | 1.500× |
| Target bootstrap 95% interval | [0.571, 3.185]× | [0.556, 4.044]× |
| Points at least 2× | 2/8 | 3/8 |
| Median online ms, baseline / orbit | 0.419 / 0.404 | 3.526 / 2.800 |
| Median charged target PDP ms, baseline / orbit | 0.328 / 0.302 | 3.278 / 2.688 |
| Median target S3 calls, baseline / orbit | 1,152 / 984 | 24,896 / 37,206 |
| Baseline/orbit target-call geometric ratio | 0.695× | 0.799× |
| A/A p95 absolute log wall ratio, baseline / orbit | 0.755 / 1.743 | 4.870 / 1.028 |

The A/A values exceed the 10% noise gate (`log(1.10) = 0.0953`), and the
host has no auditable exclusive-core isolation receipt. These wall ratios
are exploratory, with controlled speedup unknown. The all-orientation
method increases target S3 calls on this panel, and its complete one-target
online result remains below the frozen 2× advancement threshold on both
curves. Future target-query changes should reduce redundant orientation
work and be compared on these same points; the first-winner parallel path
remains the stronger N53 candidate for an isolated measurement.

The orbit [source](../../examples/s3_partial_prepared_frob_orbit_targetspan.rs)
SHA-256 is `aedf1e534c5332edf45f28ea48312ed9afd79a0e78af3ba4bc737c3e24d33489`;
the release binary SHA-256 is
`8a12ff6fa603745beebfc5229f47e928ddf9add8d2f6eabdb8cf26578eba7693`.
The plan SHA-256 is
`157ceb51830039add3a558dd19bda34b852e1f9ec2adceacc58359a384cb228d`,
and the uncompressed 64-run ledger SHA-256 is
`b02c005f346854e3a7d6f7f65c56f9303b0b45e24115942389e093eff0f5c735`.

# Batched field inversions for S3 target decomposition

Montgomery batch inversion solves the two ordered S3 group-root queries for
up to 16 consecutive index states with one field inversion. It preserves the
original query-seeded orientation and serial order within each target worker.
On the frozen N41/N53 T025–T032 points, all 64 one-target runs recovered and
verified their scalars; checked Sage independently replayed all 64 and 756
distinct relation witnesses. The exploratory paired online ratios were
1.27× on N41 and 1.15× on N53, below the 2× advancement threshold.

## Exact pipeline

Both curves use `y²+xy=x³+1` and polynomial-basis binary fields. N41 is
`EC1N41Ckb1hcbf4a37f6e66`, with modulus `x^41+x^3+1`, prime subgroup
order 549756390943, 20,008 usable factor-base points, and 244 folded
columns. N53 is `EC1N53Ckb1h738a572b91f2`, with modulus
`x^53+x^6+x^2+x+1`, subgroup order 21044858204113, 25,864 usable points,
and 244 folded columns. Both candidates start the target after 243 prepared
relation rows. The [parent freeze](../koblitz-s3-crosscurve-20261009/FREEZE.json)
fixes every public point, subgroup, input law, and workload ID. The
[experiment freeze](FREEZE.json) fixes the binary, source, target order,
workers, batch size, and candidate identities.

| Curve | Immediate parallel candidate | Batch inversion candidate | Target workers |
| --- | --- | --- | ---: |
| N41 | `IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0h8a0e7bb5a56a` | `IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0hd011289346c8` | 1 |
| N53 | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0he884cf865574` | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0h884c6acc438b` | 14 |

The batch method multiplies the nonzero denominator pairs for up to 16
states, inverts their total product once, then recovers each inverse in
reverse order. A zero denominator follows the existing scalar fallback;
states without verified prelifted group roots follow the existing S3 root
solver. A test checks ordered outputs on both curves, including zero
denominators and skipped states. N41's one-worker arm produced the same
first relation on every paired point; batching rounded its charged root-call
count up to a batch boundary. N53's first-winner race changes work completed
with thread scheduling, while every reported target relation is verified.

## Frozen single-target panel

Each of eight public points on each curve was solved twice per arm in ABBA
or BAAB order. One invocation contains one previously supplied target point.
The online timer begins after factor-base, root-index, and rank preparation;
it includes target query, all charged S3 decomposition attempts, relation
check, target descent, and scalar replay. Process launch, input load, and
known-scalar fixture construction are outside the interval. The
[raw ledger](local/raw_paired64.jsonl.gz) keeps every invocation, phase
measurement, memory peak, and status. The [analysis](local/analysis.json)
checks run IDs, code hashes, targets, phase sums, and point-block ratios.
The [Sage replay](local/sage_replay.json) with its checked
[runtime receipt](local/sage_runtime_info.json) independently checks group
equations, prepared rank, and recovered scalars.

| Measure | N41 | N53 |
| --- | ---: | ---: |
| Native and Sage verified runs | 32/32 | 32/32 |
| Immediate/batched online geometric ratio | 1.269× | 1.152× |
| Target bootstrap 95% interval | [0.904, 1.886]× | [0.761, 1.841]× |
| Points at least 2× | 1/8 | 1/8 |
| Median online ms, immediate / batched | 1.130 / 0.594 | 5.973 / 4.496 |
| Median charged target PDP ms, immediate / batched | 0.991 / 0.410 | 5.795 / 4.332 |
| Median target S3 calls, immediate / batched | 1,139 / 1,152 | 20,941 / 29,248 |
| Immediate/batched target-call geometric ratio | 0.952× | 0.908× |
| A/A p95 absolute log wall ratio, immediate / batched | 1.236 / 1.358 | 1.159 / 1.173 |

The A/A values exceed the 10% noise gate (`log(1.10) = 0.0953`). These
local wall ratios are exploratory; the controlled speedup field stays
unknown pending host-wide CPU isolation evidence. The method reduces the
field-inversion stage, but its complete one-target online interval does not
meet the 2× threshold on this panel. The next optimization should target
the remaining per-state canonical lookup and the 14-worker race's variable
work, with identical target and timing boundaries.

The batched source SHA-256 is
`d948d542448fcd2456da0a30905182e0e05fbfe347c620d9fbfcfb4b7491cd6e`,
the local release binary SHA-256 is
`83528bd170cc17db25b49c983d1e2e1e960b86cace6189176be2470ed71de531`,
the plan SHA-256 is
`f8dcd9aa09dc6aa3d46a7f67139ddc4358baa966dfce2f74114649a2901d7d3d`,
and the uncompressed raw ledger SHA-256 is
`9ce4893b441b06856d34f5e86d451a4e8e3004db32b967b488f7b1222c0443f3`.

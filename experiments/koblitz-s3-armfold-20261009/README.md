# Guarded ARM PMULL sparse-tail folding

The binary-field kernel now reduces sparse-modulus products with two PMULL
folds on ARM64 when `GF2_ARM_PMULL_FOLD=1`. The default ARM table kernel and
the portable CPU path remain available. The opt-in path passed 492,000
independent product and square checks across widths 2–63; all six S3 example
tests passed with folding enabled. On frozen N41/N53 T025–T032 public points,
all 64 single-target executions verified, and checked Sage independently
recovered every scalar and checked 756 distinct relation witnesses.

## Candidate and workload identity

The [freeze](FREEZE.json) fixes the shared source tree at commit
`3e46e47526b19b9b5cef576050ef78efffcc58c4`, the same release binary
for both modes, public points, target order, two repetitions per mode and
point, and 243 prepared relation rows. The exact curve and factor-base
records are inherited from the [parent freeze](../koblitz-s3-crosscurve-20261009/FREEZE.json).
Both curves use `y²+xy=x³+1`. N41 uses `GF(2^41)` modulo
`x^41+x^3+1`, subgroup order 549756390943 and
20,008 usable points; N53 uses `GF(2^53)` modulo `x^53+x^6+x^2+x+1`,
subgroup order 21044858204113 and 25,864 usable points. Both have 244
folded columns and use the same 16-state target group-inversion batch.

| Curve | Table candidate ID | PMULL-fold candidate ID | Target workers |
| --- | --- | --- | ---: |
| N41 | `IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0h0f4d1332b718` | `IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0h85a4418304d9` | 1 |
| N53 | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0hf313888ef69e` | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0hbd234768f9b2` | 14 |

Both arms invoke the same binary on the same point and resource envelope.
The table arm unsets `GF2_ARM_PMULL_FOLD`; the fold arm sets it to `1`.
The target online clock starts after reusable base, index, and rank
preparation and charges query, all target-dependent S3 attempts, relation
check, descent, and scalar replay. Process launch, input load, and target
fixture construction are outside it. The [raw ledger](local/raw_paired64.jsonl.gz)
preserves one row per candidate, workload, and run ID with phase costs,
memory, and statuses. The [analysis](local/analysis.json) checks identities,
hashes, phases, and point-block ratios; the [Sage replay](local/sage_replay.json)
and checked [runtime receipt](local/sage_runtime_info.json) establish the
independent equations and scalar recovery.

## Paired full-operation result

| Measure | N41 | N53 |
| --- | ---: | ---: |
| Native and Sage verified runs | 32/32 | 32/32 |
| Exploratory table/fold online geometric ratio | 0.870× | 1.061× |
| Target bootstrap 95% interval | [0.431, 1.625]× | [0.616, 1.759]× |
| Points at least 2× | 2/8 | 1/8 |
| Median online ms, table / fold | 0.573 / 0.942 | 10.770 / 7.485 |
| Median charged target PDP ms, table / fold | 0.355 / 0.837 | 10.078 / 7.123 |
| Table/fold target S3-call geometric ratio | 1.000× | 0.751× |
| A/A p95 absolute log wall ratio, table / fold | 1.288 / 1.362 | 3.147 / 1.559 |

N41's one-worker runs have identical first relations and S3-call counts
across modes. N53's first-winner race changes completed work with scheduling.
The A/A wall ratios exceed the 10% gate (`log(1.10) = 0.0953`), and the
[host record](hardware.json) lacks an exclusive CPU partition, affinity,
fixed frequency, IRQ and throttling evidence. The controlled speedup field
remains unknown. Neither curve meets the 2× advancement threshold in the
full single-target panel.

The [fixed-operand kernel diagnostic](local/kernel_analysis.json) processes
1,048,576 multiplications per repetition, nine repetitions per process, in
table/fold/fold/table order. Its exploratory table/fold geometric ratios are
1.164× at N41 and 2.831× at N53. The N53 table A/A absolute log ratio is
1.289, so the larger kernel figure is not a stable measured gain. All modes
produced the same nonzero checksum. The raw process records and
[validator outputs](local/validate_fold.txt) are retained under `local/`.

The shared arithmetic source SHA-256 is
`3a2bad8583d5c36bc17ba28e625943da85af090eb2e46f8cf11d75108f62a7ae`,
the release binary SHA-256 is
`c63f0bf9bd637b29620233e7a2d1f34648de0dbb9ab9a49ee69e9c738dffc34c`,
the plan SHA-256 is
`ebb0175e784ee68c6654ddb43401c2803784622848e48afe94b60e8d85407340`,
and the uncompressed raw ledger SHA-256 is
`c37d327aff3d00663f61b4c67249cca1b8b02e4ef6807c12a0cc346cf45da1e2`.

The repository's complete library test target could not compile from this
sparse checkout because unrelated `include_str!` fixtures are absent. The
standalone polynomial reference check, six release S3 tests, 64 native
single-target solves, and Sage replay exercise the changed arithmetic path
on this physical Apple M4 Pro. Physical x86-64 and non-Apple ARM devices
remain to be validated before any platform-specific performance claim.

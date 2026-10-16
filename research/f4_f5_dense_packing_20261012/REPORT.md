# Boolean F5 packing with bounded direct indexing

The enhancement reuses F4's bounded direct column-index construction in the
shared Boolean row packer used by F5. On the fixed eleven-system suite,
complete F5 calls use **3.9471779705% fewer instruction references** under the
preregistered paired estimator, with 95% bootstrap interval
**[3.9239302372%, 3.9573654824%]**, and request **2.0332272879% fewer bytes**.
F4 is neutral. Exact reduced outputs, inputs, ranks, reports, criterion
membership and word-XOR accounting match; no measured memory metric regresses.
These are instruction and allocator results on bounded Boolean systems.
Elapsed-time, process RSS and whole-solver gains remain unmeasured.

The [protocol](PROTOCOL.md) fixes the F5-specific gate before production edits
and measurements. This extends [PR #1723](https://github.com/aburan28/crypto/pull/1723)
and preserves the previous [index/readback/fusion study](../f4_f5_fused_readback_20261012/REPORT.md),
including its rejected attempts and gate-scope clarification.

## What changed

[`pack_rows`](../../src/cryptanalysis/koblitz_groebner.rs) previously built a
mixed-mask hash map and looked up every term. It now tries the existing bounded
`u16` direct table, using the same construction as F4. It selects the packing
loop once per matrix. Generic specialization keeps the dense and hashed loops
separate, avoiding representation dispatch for each fallback term.

The table is admitted only when its bytes fit within the minimum bytes of the
hash entries. Empty, sparse/high-bit, overflow and greater-than-65,535-column
layouts retain hashing; the missing sentinel is checked before writing a bit.
`KIC_F4_HASH_COLUMN_INDEX=1` forces the shared hash path for compatibility
checks. Row order, column order, duplicates' OR packing, caps, pruning and
elimination are unchanged. The same sorted columns and cancellation semantics
remain required by the matrix builders.

## Observed comparison

[Native comparison](records/comparison.json), from 18 fresh processes and 792
audited calls: three primary repetitions per variant, three A/A repetitions,
five interleaved A/B/A′ rounds, all eleven systems and three measured arms.
The separate hash-only candidate control adds 33 calls; its
[native verification](records/hash-control-verification.json) passes.

| interval | baseline median suite Ir | candidate median suite Ir | paired percent fewer Ir | native verdict |
|:--|--:|--:|--:|:--|
| complete F4 | 93,241,594 | 93,258,282 | −0.0198628515% | neutral |
| complete F5 | 82,213,799 | 78,971,828 | 3.9471779705% | fewer instructions |
| standalone criterion | 3,383,582 | 3,361,188 | 0.6618429818% | neutral |

The percentages come from median paired log-ratios, not arithmetic on the two
separately reported endpoint medians. The 2% verdict floor, A/A MAD control,
2,000-resample percentile bootstrap and fixed seed are unchanged. Native
`totals` are reconciled with function self-costs; wrapped Callgrind `summary`
fields are retained rather than used as costs. Per-case values and intervals
are fully recorded and plotted.

Nine F5 cases have fewer-instruction verdicts. The degree-three twelve-variable
case (−0.4554944249%) and shifted high-bit case (−0.2581906405%) are neutral;
their small increases remain in the data. All eleven F4 cases are neutral.
There is no more-instruction verdict among the 22 complete-call cases.

One standalone criterion case, `n8_d4` at seed `0x9e3779b97f4a7c15`, records
**15.4494279343% more instructions**. The two twelve-variable criterion cases
have fewer-instruction verdicts; the other eight are neutral. The criterion
source is unchanged, but its generated instruction cost in the complete
driver's execution context is not invariant; the cause of the changed phase
costs is not established. This regression is preserved and visible in the
per-case graph. The frozen gate treats criterion Ir as a separate diagnostic
and still guards every criterion allocation interval.

`f5_with_complete_call_guards_met=true`: F5 aggregate fewer, F4 aggregate not
more, no complete-call case more, and no memory regression in any primary or
paired interval. `acceptance_met=true`. The older, stricter
`two_arm_acceptance_met=false` and `complete_call_acceptance_met=false` remain
unchanged fields: both require F4 itself to have a fewer-instruction verdict.
They are not the preregistered acceptance rule for this F5-specific enhancement.

| primary totals, 33 calls per arm | baseline requests | candidate requests | baseline requested bytes | candidate requested bytes |
|:--|--:|--:|--:|--:|
| F4 | 51,783 | 51,783 | 9,129,570 | 9,129,570 |
| F5 | 51,174 | 51,174 | 9,660,012 | 9,463,602 |
| criterion | 6,831 | 6,831 | 784,656 | 784,656 |

No allocation-request, cumulative-byte, absolute-peak or peak-above-entry
regression occurs in primary or paired calls. F5 peak-above-entry medians fall
on nine cases; the two neutral fallback cases are unchanged. For example,
the twelve-variable degree-four F5 interval falls from 378,364 to 368,622 live
bytes above entry. Inputs and retained warm outputs remain in entry bytes.
These are allocator bytes, not RSS. Absolute peaks, entry values and every
interval remain in the raw JSONL and native comparison.

## Reproduction and provenance

- Production baseline: `0677b1eb0eb4275f5c3366a48689d5620c0bf0b8`.
- Protocol/common harness commit: `28608a4d6879448764bfafeff7a98cb7123fbe78`.
- Measured production candidate: `d56c96b926569c8cd4bfe027802ed9342c2ca36d`.
- Test-only parallel-fixture adjustment: `29d153ab97b92811b42abdab3126b18be8669896`.
- Frozen driver SHA-256: `2d221b46f04ab815ad74fade62ddcf9c89fc62d7f452e399648c3db942a3c0b8`.
- Baseline binary SHA-256: `7c91eeb582e389d57ad441aab233addeb9fcc3d0b50879d63512778e038f95fa`.
- Measured candidate binary SHA-256: `bfc93c705b347dae74e3992d7d7eceee44ce176fcdf2914b5a94af053cba2a17`.

The final normal driver was rebuilt after the test-only fixture adjustment.
Its [native ELF verification](records/runtime-elf-verification.json) confirms
identical code/data, all 30 allocated-section layouts and program headers,
apart from permitted GNU build-ID contents. The measured and rebuilt
[whole-file hashes](records/rebuilt-driver-sha256.txt) remain separate;
runtime equivalence does not assert whole-file identity. The checker source
is the frozen [native ELF helper](../f4_f5_fused_readback_20261012/records/verify_elf.rs).

[`reproduction.md`](records/reproduction.md) gives build, measurement,
analysis and rendering commands. [`run-profiles.sh`](run-profiles.sh) only
orchestrates compiled Rust execution, refuses existing run directories and
retains commands, timestamps, environment, JSONL, raw dumps and actual exit
codes. [`verify.sh`](verify.sh) retains independent checks and failures.
No crypto execution, correctness, replay or result computation uses Python.
The [renderer](render.py) lays out recorded fields only.

Host: AMD EPYC 9V74 Linux x86-64, five exposed logical CPUs, cgroup quota four;
Rust 1.99.0, baseline release target, Callgrind 3.24, one driver Rayon worker.
Runtime AVX2/BMI2 dispatch remains intact; no global target-CPU flag is added.
[CPU](records/cpu.txt), [OS](records/os.txt), [compiler](records/compiler.txt),
[quota](records/cpu-quota.txt), [memory](records/memory.txt) and hashes are
retained. CPU-pressure telemetry is unavailable, so controlled timing was not
attempted. Ir does not require uncontended wall-time eligibility.

All old workloads remain unchanged. The three previously frozen holdouts are
now known regression fixtures, not newly unseen inputs. The shifted input has
64 polynomial coordinates supported on bits 52–63 and six multiplier variables.
The default and hash-only variants return identical values; no workload, seed,
cap, control or gate was removed after observing the results.

## Verification and delivery status

| requested requirement | status | evidence / remaining check |
|:--|:--|:--|
| enhance existing implementation | verified complete | source change, guarded native comparison, release tests |
| compare with prior implementation | verified complete | 792 calls, native exact-output and allocation guards |
| retain fallback behavior | verified complete | 33-call hash control; independent library boundary test |
| repository checks | verified complete | [validation](records/validation.json), [successful full repeat](records/verification-network/library.log), focused/example/lint logs |
| plot groups, graphs, PDFs | verified complete at publication | [six-page report](reports/run-20261012T1503Z/report.pdf), editable sources, inspection and reproduction receipts |
| commit/push/update PR | delivered through PR #1723 | implementation and evidence commits on the existing feature branch |

The independent packing test reconstructs bits from column positions, covering
empty rows, constants, arbitrary order, holes, duplicate terms, word boundaries,
padding, sparse/high-bit fallback, the sentinel limit and the Rayon path.
Existing F4/F5 row-space and solver checks remain intact. Native comparator
tests and strict comparator/driver Clippy pass. The focused suites pass 47
Boolean Gröbner tests, 18 F5 tests and 21 PQ Gröbner tests. Cargo's two example
targets pass five and zero tests respectively.

The first full run passes 3,987 and fails nine mock-network tests because the
default sandbox denies their loopback socket creation; 112 are ignored.
The [environment correction](VERIFICATION-ENVIRONMENT.md) retains this failure
and repeats the same checks with network permission. The same test binary's
37 mock-network checks then pass without source or ignore changes. The same
mandatory full Cargo command then passes **3,996**, with **zero failures and
112 ignored**, exit zero, in **238.05 seconds** under its 1,200-second budget.
All 4,108 catalog names are accounted for; no unfinished or unexpected tests
remain. Ignored tests stay ignored and named; no ignore flags are broadened.
This run does not repeat the prior study's separately passed long isogeny check,
which is unrelated to row packing and remains linked historical evidence.

Both release driver builds retain nine pre-existing library warnings. No
unrelated warning or test is changed. All nineteen native profile processes
complete with exit zero; all unsuccessful verification evidence remains intact.

The six graph groups have separate one-page PDFs and a combined six-page PDF.
All twelve PDF pages are inspected; an independent render reproduces the
SVGs, 150-dpi PNGs, Typst sources and PDFs byte for byte. Artifact hashes and
inspection receipts are retained beside the records. Benchmark browsing copies
match the source artifacts.

## Canonical visual coverage

The new history uses three separately matched episodes, with gaps rather than
compounding ratios across different baselines. The complete-call history,
criterion history, per-case instruction intervals, request counts, requested
bytes and additional peak each have editable sources, SVG/PNG and inspected
group/combined PDFs. Every core/compatibility process is covered by this one
bounded campaign report, with raw identifiers and controls beside it.

`docs/perf/PERFORMANCE_INDEX.md` and `benchmarks/README.md` receive the new result.
Both this study and the prior study in the same PR receive benchmark browsing
copies. The prior frozen SVGs/PDFs are reproduced byte for byte; new 150-dpi PNGs
come from those same figures via [render-prior.py](render-prior.py).

Existing `docs/perf/f4-f5-allocation-history.*`, group registry and imported CI
plots remain frozen to their original eight-system/legacy-profiler inputs.
Their fingerprints, interval and worker coverage differ from this campaign.
The previous study-local figures remain snapshots, represented as separately
matched segments in the new history. The planning `NEXT_STEPS` diagram remains
its dated proposal. `docs/index-calculus-scoreboard.html`,
`docs/ic/progress-timeline.json`, curve diagrams and whole-method dashboards
remain unchanged because this is a matrix-packing result with no IC/rho or
full-solver ratio, curve discovery or promoted ledger claim.

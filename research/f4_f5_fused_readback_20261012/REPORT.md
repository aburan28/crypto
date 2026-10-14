# Boolean F4/F5 construction and indexing: native follow-up

The corrected implementation records **5.2805% fewer complete-F4 instructions**
and **11.0184% fewer complete-F5 instructions** on the frozen eleven-system
suite. All 22 complete-call cases pass the native instruction verdict, outputs
and accounting match, and no allocation or peak-live metric regresses.
The revised complete-call delivery gate passes. **The original all-33-case
gate fails:** three separately invoked criterion cases use more instructions.
The criterion aggregate is neutral. Runtime and process RSS were not measured.

This is Boolean matrix-engine engineering, including construction, packing,
reduction and readback, rather than a measured complete Gröbner solver or
cryptanalytic method improvement. Elimination, pruning rules, matrix order,
size caps and counted word-XOR work are unchanged. No degree, exponent,
crossover or coefficient-field generalization follows from this suite.

## Code and measured reference

- Production baseline: `c2963dcff359814bf8203b946b88e9bb62431795` (current
  `main` after [PR #1677](https://github.com/aburan28/crypto/pull/1677) merged).
- Corrected production candidate: `d0a1d2d88121d868351760b1285c57c9e835a1fc`.
- Common driver/protocol revision: `d4281243b9bff596837db20bfcb9160f363fae1d`;
  the measured driver is byte-identical source to the earlier hashing study.
- Final native analysis revision: `c8e1bcb6`; exact revision and hashes are
  [recorded separately](records/attempt-2-fixed/analysis-revision.txt).
- Host: AMD EPYC 9V74, cgroup CPU quota 4 and exposed CPUs 0–4; recorded
  Rust 1.99 and locally extracted Callgrind 3.24. One Rayon worker per driver.
  [Compiler](records/compiler.txt), [host](records/host-cgroup.txt),
  [profiler](records/profiler-version-2.txt) and binary/source hashes are retained.

F4 adopts unique descending-column readback terms directly, avoiding an
unnecessary sort and second term vector. Its column layout uses direct `u16`
lookup only when the full mask-range table fits within the minimum hash-entry
storage. Missing entries use `u16::MAX`; oversized, sparse or high-bit ranges
retain the mixed-mask hash map. Standard hashing remains an explicit control.
The nested fast representation preserves the original heap-resident metadata
size. `KIC_F4_HASH_COLUMN_INDEX=1` selects the mixed-mask control.

F5 enables the existing fused full-F4-row count and selected-row construction
by default. It generates each Boolean product once, still counting nonempty
pruned rows toward the original full row cap. `KIC_F5_FUSED_BUILD=0` restores
the two-pass control. Empty-call behavior, reports and direct-pack handling
are preserved.

The fused pass can retain selected rows before an oversized call reaches its
full row-count limit; the former count-only preflight could reject without
retaining them. Caps and `None` behavior remain unchanged, but oversized-call
allocation peaks are not profiled in this cohort. The explicit two-pass flag
retains that preflight. Memory nonregression claims here cover the frozen suite,
not all possible inputs or failure paths.

## Frozen comparison and its scope

[PROTOCOL.md](PROTOCOL.md) freezes eight prior seeded systems plus three
holdouts before production changes. They cover 8, 10 and 12 multiplier
variables, degrees 3–5 and a 64-bit-coordinate holdout with support at
`[52,54,57,60,62,63]` but six multiplier variables. Exact seeds and the native
driver are linked there. The source hashes show unchanged multiplication,
elimination and workload definitions.

Each cohort retains three primary repetitions per variant, a three-repetition
baseline A/A control, and five rotated A/B/A′ rounds of all eleven systems
and three arms. This is 792 audited core calls plus 33 compatibility calls.
The two original baseline processes are reused byte for byte in later
cohorts; they are **shared samples**, not extra independent repetitions.
All profile processes finished successfully under the 180-second budget.
Every command, UTC start/end, exit and dump remains beside its JSONL records.

The interval includes allocator instrumentation, construction, packing,
elimination and readback; input creation, warmup, digesting, verification and
serialization are outside it. Merged dumps include workers. The native Rust
reader requires positive `totals`, sums function self-costs independently,
excludes inclusive call edges, and retains Callgrind's wrapped `summary`
fields unchanged. All input/output and criterion-membership digests, reports,
ranks and XOR accounting match, including F4/F5 cross-output equality.

The [native final comparison](records/attempt-2-fixed/comparison.json) computes
median paired log-ratios and 95% bootstrap intervals from five rounds, with
2,000 resamples and frozen seed `0x243f6a8885a308d3`. Fewer-instruction verdicts
require an interval excluding zero and an effect above both twice the A/A
noise estimate and 2%. Medians of each arm and paired effects are separate
statistics; do not derive a different verdict from the ratio of two medians.

| Arm | Baseline median Ir / 11-system suite | Candidate median Ir | Paired instruction reduction | 95% bootstrap interval | Native verdict |
| --- | ---: | ---: | ---: | ---: | --- |
| Complete F4 | 98,440,694 | 93,248,010 | 5.2805% | 5.2743–5.2843% | Fewer instructions |
| Complete F5 | 92,393,939 | 82,211,306 | 11.0184% | 11.0157–11.0315% | Fewer instructions |
| Standalone criterion | 3,357,250 | 3,384,189 | −0.7968% | −0.8081–−0.7910% | Neutral |

Per-input values below are the native paired instruction reductions, rounded
for display. `s1=0x2545f4914f6cdd1d`, `s2=0x9e3779b97f4a7c15`; all holdouts
use `0xd1b54a32d192ed03`. All F4/F5 entries have the fewer-instruction verdict;
the [complete native record](records/attempt-2-fixed/comparison.json) retains
every case's interval and A/A control, including neutral phase cases.

| Input | Complete F4 | Complete F5 | Standalone criterion | Criterion verdict |
| --- | ---: | ---: | ---: | --- |
| n8/d3/s1 | 16.11% | 13.66% | −2.31% | More instructions |
| n8/d3/s2 | 15.61% | 15.03% | −5.16% | More instructions |
| n8/d4/s1 | 6.31% | 11.00% | −0.13% | Neutral |
| n8/d4/s2 | 6.32% | 11.97% | −0.25% | Neutral |
| n10/d4/s1 | 6.07% | 11.21% | −0.88% | Neutral |
| n10/d4/s2 | 5.49% | 11.84% | −0.73% | Neutral |
| n10/d5/s1 | 4.12% | 9.88% | −0.35% | Neutral |
| n10/d5/s2 | 4.16% | 10.17% | −0.38% | Neutral |
| Holdout n12/d3 | 16.57% | 17.39% | −2.57% | Neutral |
| Holdout n12/d4 | 5.17% | 11.68% | −4.07% | More instructions |
| Holdout shifted/d4 | 5.69% | 6.41% | −1.45% | Neutral |

The criterion source/rules are unchanged in this iteration. Its standalone
invocations nevertheless record the differences above; their cause is not
established by this campaign. Complete F5 already includes criterion cost.
The separately invoked diagnostic is not added again to that complete-call
cost. It establishes no criterion improvement.

## Allocation and live-byte results

These totals cover 33 primary calls per arm, not elapsed time or RSS.
Requests, cumulative requested bytes and live-at-entry/absolute/additional
peak live bytes remain distinct. Every primary and paired workload is checked
for request, byte and both peak regressions; none occurs in the corrected cohort.

| Arm | Baseline requests | Candidate requests | Request reduction | Baseline cumulative bytes | Candidate cumulative bytes | Byte reduction |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Complete F4 | 63,498 | 51,783 | 18.4494% | 10,231,140 | 9,129,570 | 10.7668% |
| Complete F5 | 51,744 | 51,174 | 1.1016% | 9,802,092 | 9,660,012 | 1.4495% |
| Criterion | 6,831 | 6,831 | 0% | 784,656 | 784,656 | 0% |

The corrected high-bit F4 holdout has identical entry bytes (43,714), absolute
peak (97,818) and additional peak (54,104) on both sides. Earlier metadata
increased that last value by eight bytes; its rejection is retained below.
Live-byte intervals include retained inputs/warm outputs at entry. Subtracting
entry from that same interval's native live-byte peak does not measure RSS.
The [timing eligibility record](records/timing-eligibility.txt) records missing
CPU pressure telemetry; no native elapsed-time or fresh-process RSS campaign
was admitted. Those roadmap items remain pending.

## Rejected attempts and explicit protocol changes

1. [First candidate](records/comparison.json), `07151ff0`: F4 1.4363% fewer
   instructions (neutral), F5 11.1305% fewer, no memory regression. It failed
   the original gate. Two criterion cases also regressed. All raw profiles,
   controls and the 45-test focused pass remain unchanged.
2. [Second candidate](records/attempt-2/comparison.json), `07e568b1`: complete
   F4 5.5122% fewer instructions and F5 11.0349% fewer; high-bit F4 peak
   live bytes grew eight bytes consistently. Three criterion cases also
   regressed. The original gate failed. Its 46 Gröbner, 18 F5 and 21 polynomial
   focused tests passed; its superseded full-library run was deliberately
   interrupted (exit 143), **not passed or timed out**.
3. The [native metadata probe](records/layout-probe.json) measures original /
   flat-three-way / nested index sizes at 48 / 56 / 48 bytes. The
   [bounded correction](METADATA-FIX.md) explicitly expands the two-attempt
   stopping rule by one correction and complete rerun, with unchanged inputs,
   statistics and memory guards. No new optimization hypothesis follows it.
4. [ACCEPTANCE-CLARIFICATION.md](ACCEPTANCE-CLARIFICATION.md) explicitly narrows
   the delivery gate from the extra all-33-case phase guard to complete calls,
   before corrected measurement. `complete_call_acceptance_met=true`;
   the preserved original `two_arm_acceptance_met=false`. This is a gate-scope
   change after earlier phase regressions were observed, not a pass of the
   original protocol. Statistical thresholds and all-arm memory checks are
   unchanged. Every phase regression remains visible.

Tool failures are retained: the comparison compiler's initial type-inference
error and a missing `VALGRIND_LIB` version query, followed by successful
recorded fixes. No profiler workload, seed, failure or regression was dropped.
Conductor was unavailable; [the scope attempt](records/conductor.log) is retained.

## Verification and reproduction

The corrected-source focused tests pass: 46 Gröbner/F4, 18 F5 and 21 polynomial
tests. New checks cover word-boundary readback, padding, bit 63, zero/constant
and cancelled products, all row-cap boundaries, missing direct-index entries,
sentinel/overflow fallback, and original metadata-size preservation. Existing
random F5 row-space/order and exact-layout/solver tests remain active.

The first corrected-source full-library check failed: 3,994 passed, one failed,
112 ignored (303.77 seconds, exit 101). Its isolated reproducer also failed.
[TEST-FIX.md](TEST-FIX.md) records the inherited stale Montgomery fixture,
the empty baseline source diff and the test-only repair. The first repair exposed
a second stale detail: the extension fixture lacked its required modulus
(seven coordinate tests passed, one failed). Both failed runs are retained.
The complete fixture now proves the supported order, rejects the old wrong
order, preserves the unavailable-extension status and rejects its incomplete
serialization. Production curve arithmetic and registered records are unchanged.

Final verification **passed**:

- Locked release-library check: **3,995 passed, zero failed, 112 ignored**,
  314.69 seconds, exit 0, within the existing 20-minute limit. Named coverage
  accounts for all **4,107** catalogue tests, with no unfinished or unexpected
  names. [Final log](records/verification-final-2/library.log) and
  [coverage](records/verification-final-2/library-completed-names.txt).
- Unmodified long isogeny test, explicitly invoked despite its current-main
  ignore: **one passed**, zero failed, 822.65 seconds, exit 0 within its separate
  20-minute limit. [Exact command/log](records/verification-final/jv.log).
  The standard suite still reports 112 ignored; this separate pass is not
  silently folded into that count.
- Coordinate repair: all **24** native module tests and all **eight** focused
  Cargo coordinate tests passed. The temporary native module build initially
  lacked an `isogeny` re-export; its failed build and corrected pass are retained.
- Native comparison: all **five** tests passed, recomputation byte-identical;
  Cargo example targets passed **five + zero** tests. Both example Clippy checks
  passed. [Example receipts](records/verification-final-2/example-tests.log).

The fixture-only normal rebuild has a different ELF file hash. Its byte
comparison is retained as a mismatch, not relabeled identical. The
[native ELF verifier](records/verify_elf.rs) confirms identical loader headers
and all **30** allocated-section layouts, code and data, except the GNU
build-ID contents. Nonruntime symbol metadata differs. The
[verification receipt](records/verification-final-2/runtime-elf-verification.json)
therefore establishes runtime equivalence to the measured driver; both hashes
remain recorded. The measured production revision and comparison are unchanged.

The [two-flag compatibility control](records/attempt-2-fixed/legacy-control-verification.json)
matches all outputs and accounting in a fresh process; its environment file
records both legacy flags. The [cross-study native check](records/history-bridge-verification.json)
matches 99 prior hashing-candidate and fresh-baseline outputs/accounting.
Recomputing the final comparison gives byte-identical JSON. Rust executes,
verifies, parses and computes all results; shell only builds/runs; Python only
lays out the recorded values and renders PDFs. Reproduction commands and
named validation coverage are recorded with the final checks.

## Graph groups and PDFs

Campaign/report id: `run-20261012T0543Z`. Editable renderer:
[render.py](render.py); pinned dependencies:
[groebner-render-requirements.txt](../../scripts/perf/groebner-render-requirements.txt).
The final rendering keeps each SVG, Typst source and one-page group PDF beside
the [combined report PDF](reports/run-20261012T0543Z/report.pdf).

All six combined pages and all six one-page group PDFs were rasterized and
opened for final visual inspection. The instruction revision labels were
shortened after draft inspection to avoid overlap. An independent final render
reproduces every SVG, Typst source and PDF byte for byte. The
[inspection receipt](records/report-inspection.json) and hashes are retained.

- [Complete-call instruction history](reports/run-20261012T0543Z/complete-instructions.svg)
  and [PDF](reports/run-20261012T0543Z/complete-instructions.pdf).
- [Criterion instruction history](reports/run-20261012T0543Z/criterion-instructions.svg)
  and [PDF](reports/run-20261012T0543Z/criterion-instructions.pdf).
- [Allocation requests](reports/run-20261012T0543Z/allocation-requests.svg)
  and [PDF](reports/run-20261012T0543Z/allocation-requests.pdf).
- [Additional peak live bytes](reports/run-20261012T0543Z/peak-live.svg)
  and [PDF](reports/run-20261012T0543Z/peak-live.pdf).
- [First rejected attempt](reports/run-20261012T0543Z/rejected-attempt.svg)
  and [PDF](reports/run-20261012T0543Z/rejected-attempt.pdf).
- [Rejected metadata regression](reports/run-20261012T0543Z/rejected-memory.svg)
  and [PDF](reports/run-20261012T0543Z/rejected-memory.pdf).

The instruction histories show the earlier hashing pair and this corrected
pair as disconnected measured segments. Do not compound a ratio across
unrelated intervening `main` changes. All points, samples and intervals come
from their native comparisons. This is the latest native group linked from
[the performance index](../../docs/perf/PERFORMANCE_INDEX.md).

Graphs checked and left frozen: the original allocation study, the historical
CI n=17/n=23 instruction/runtime groups, the previous hashing study and original
planning-snapshot diagram/PDF. Their runs, suites or recording contracts have
not changed. The aggregate history, IC scoreboard and progress timeline gain
no new solver, curve, relation-yield or whole-method measurement in this change.
No old measurements or diagnostic timings are replaced by these instruction
results.

## Requirement status

| Requested or required deliverable | Current status / evidence |
| --- | --- |
| Optimize the existing Boolean F4/F5 implementation | Implemented; complete-call comparisons and all output/accounting controls verified |
| Retain original gates and all rejected evidence | Verified; two rejected comparisons and explicit bounded/gate-scope amendments retained |
| Validate memory and high-bit fallback | Verified; no final regression, exact high-bit peak, independent lookup and metadata checks |
| Required library and example checks | Verified: 3,995 library passes, 112 ignored; explicit long test passed; focused/example/lint checks passed |
| Graph groups, PDFs and editable sources | Verified: six SVG/group-PDF sets, six-page combined PDF, all 12 pages inspected, byte-identical independent render |
| Commit, push and ready follow-up PR | Complete repository deliverable; the feature branch's ready PR is the publication receipt |
| Controlled runtime, fresh-process RSS, whole-solver scaling/crossover | Not established; separate roadmap obligations remain pending |

# Boolean F4 readback and fused F5 row construction

Requested: continue optimizing the existing F4/F5 implementation and retain
the PR, graph groups and PDF reporting convention. The previous PR #1677 is
merged. This follow-up starts from production baseline
`c2963dcff359814bf8203b946b88e9bb62431795` on a fresh feature branch.

Hypotheses fixed before production edits: F4 can adopt the already ordered,
unique terms emitted by packed-column readback instead of sorting and copying
them through `from_monos`; F5 can use its existing fused full-row count and
selected-row builder by default instead of generating products twice. Keep
`KIC_F5_FUSED_BUILD=0` as the explicit legacy control. No changes to elimination,
pruning rules, row ordering, column ordering, size caps or word-XOR accounting.
The previous production profile is retained as a source-linked bottleneck
diagnostic, not a measurement of the new baseline.

Freeze the identical native driver `examples/f4_f5_row_build_bench.rs` from
the previous hashing study. Eleven workloads: eight planted quadratic Boolean
systems `(n,m,D)=(8,8,3),(8,8,4),(10,10,4),(10,10,5)`, each with seeds
`0x2545f4914f6cdd1d` and `0x9e3779b97f4a7c15`; holdouts `(12,12,3),(12,12,4)`
and the degree-four, six-equation system supported on bits
`[52,54,57,60,62,63]`, seed `0xd1b54a32d192ed03`. The shifted input has
64-bit polynomial coordinates and six low-bit multiplier variables.

Build the common driver before production edits and retain the binary and
source hashes. Rebuild that identical driver with candidate production source.
Fix host, Rust compiler, flags, Callgrind 3.24 and one Rayon worker. The
comparison tool gains an optional inspection arm; its execution/statistics
contract stays unchanged. Record production and harness revisions separately.

Measure complete F4 calls, complete F5 calls and standalone criterion construction
in instruction references (Ir), with common allocation instrumentation. Include
packing, elimination and readback. Exclude input construction, warmup, output
digesting, verification and serialization. Use client start/dump/stop requests
and `--separate-threads=no`, including workers. Native Rust performs execution,
verification, profile parsing, result arithmetic and statistics. Shell only
builds/runs; Python may lay out recorded values and render PDFs.

Retain three primary repetitions per variant, a three-repetition baseline A/A
control, and five interleaved A/B/A′ rounds (one repetition per system and arm,
rotating process order). Eighteen processes and 792 audited calls, each process
with a 180-second limit. Preserve every failed or incomplete run. Do not remove
seeds or workloads adaptively. Callgrind's `summary` can wrap after zero-stats;
use positive `totals` and independently sum self-cost records, excluding
inclusive call edges. Keep wrapped fields unchanged.

Additional compatibility control, declared before candidate measurement:
one fresh candidate process with `KIC_F5_FUSED_BUILD=0` and one repetition
of all eleven systems/three arms. The native checker compares its outputs
and accounting to the primary baseline. These 33 extra calls are outside
the fixed 792-call, five-pair acceptance campaign.

Correctness: identical input/output and criterion-membership digests, reports,
rank, complete and criterion word-XOR accounting; identical reduced F4/F5
outputs on each input. Independently test sorted readback over multiple packed
words, padding and bit 63; check fused row counting at all cap boundaries with
zero, constant, cancelled and pruned rows using polynomial multiplication as
the reference. Retain existing solver, blocked F4 and F5 row-space checks.

Acceptance: both complete-call arms must pass the native instruction gate;
no per-workload instruction regression above that gate; no allocation-request,
cumulative-byte, absolute-peak or additional-peak-live-byte regression on any
workload. Use median paired log-ratios, the A/A median absolute deviation
noise floor, a 2% minimum verdict threshold and 95% percentile bootstrap with
2,000 resamples and seed `0x243f6a8885a308d3`. Report neutral cases and all
regressions. The existing comparator's F5/nonregressing-F4 aggregate gate is
supplemented by the stricter two-arm requirement here. If either hypothesis
fails acceptance, retain its rejected evidence and revise or remove that change.

Native counters record requests, cumulative requested bytes, live-at-entry and
absolute peak live bytes; derive additional peak above that same interval's
entry. These are allocator bytes, not RSS. Controlled native timing is not
earned on this host (CPU pressure telemetry unavailable); no elapsed-time
speedup or process RSS claim. General field/prime-field F4/F5 and whole-solver
timing remain outside this bounded Boolean engineering comparison.

Verification: focused Gröbner/F4/F5 tests, the changed example's tests and
Clippy; then `cargo test --locked --release --lib -j2 -- --test-threads=3`
with the existing 20-minute budget and normal Rayon defaults. Preserve failures,
timeouts and unfinished-test coverage; do not weaken a test to assert success.

Delivery: versioned native records, source/binary hashes, comparison, Markdown
report, separate instruction and allocation graphs, inspected group and combined
PDFs and editable sources; link from the performance index. Existing frozen
historical/diagnostic plots retain their original inputs and units. Open a new
ready PR and stop without monitoring CI. Conductor is unavailable; the failed
scope reservation is retained in `records/conductor.log`.

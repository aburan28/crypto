# Generic exclusive phase controls

Classification: **accounting**. This change partitions existing work; it does
not improve a solver, qualify a performance comparison, or consume an
improvement round. Round one retained the incumbent and two rounds remain.

The opt-in worker separates reusable preparation, ordinary-query construction,
PDP, relation checking, matrix construction, final subgroup-field LA and the
target-dependent stages. Its single-target online interval excludes reusable
preparation and ends after scalar replay. Failed attempts remain charged.
Unentered phases remain null. An external process interval assigns launch,
input and reporting costs to setup exactly once. Boolean Macaulay reduction
belongs to PDP, not final relation LA.

The measured worker requires one Rayon thread and keeps phase-changing
collection work on its controlling thread. Its strict session rejects a
boundary executed on another thread. Scope guards restore containing phases
and disable tracing on errors. The ordinary library API retains parallel
collection when tracing is disabled. As a result, enabled/disabled comparisons
measure a whole-mode change, including collection scheduling; they cannot
isolate timer overhead.

The frozen local panel contains 49 jobs and three process repetitions per job:
seven PDP backends with dense/sparse final LA on both n9 curves, n9/n13
subgroup-orbit collection-window controls, seven deliberate preparation
failures, and both n9 rho references. The initial pass verified all 147 pairs:
126 complete pairs and 21 intended incomplete pairs. Query/counter histories,
factor-base census and independently replayed certificates agree between modes.
Every exclusive native phase partition closes exactly.

The final source, including the strict foreign-thread guard, passes the same
147-pair panel with the same 126/21 completion split. Its raw reports are in
`final/`; they supplement the initial pass. The final six worker controls also
pass, with the isolated strict tests explicitly enabled.

The initial descriptive exclusive/legacy online ratio had median 0.799942 and
range 0.020024–578.530638. These highly unstable clocks come from an unpinned
shared macOS host. They establish neither low overhead nor a speedup. Raw
clocks remain in the evidence, while comparative online/cold cost, instruction
cost, S, boundary ratios and speedup stay unset. The final pass remains unstable:
median 0.821059, range 0.011945–9.820866. Calibration is inconclusive in both
passes; choosing the more favorable pass would not repair it.

The independent phase adapter rejects missing, duplicate, mixed-schema and
misattributed intervals. It checks both the contiguous Callgrind part sequence
and the whole-process instruction checksum, including the termination tail.
CI pairs all 44 existing native controls with exclusive mode and profiles 13
exclusive cases covering every backend, sparse LA, rho, a collection walk and
an unstarted target. CI results belong to the PR's checks, not to the local
macOS evidence.

## Evidence and reproduction

`PROTOCOL.md` fixes the panel and accounting. The initial `worker-raw.jsonl`
contains all 294 reports; `summary.json`, `inputs.json`, fixture files and
`environment.json` retain the jobs, outcomes and identities. The four
`initial-*` source artifacts reproduce the source bytes before strict-session
hardening and later test extensions. Their uncompressed hashes are checked
against the initial environment manifest. The initial compile failure is
retained separately; it was a missing closure parameter type, fixed before
the successful controls.

`final/` contains another 294 raw worker reports, fresh fixtures and source,
binary and input hashes. Both sets of 147 exclusive reports are replayed by
the regression suite without trusting their saved phase receipts.

`final/initial-regressions/` retains a failed test expectation: the singular
dense matrix produces a provisional zero vector, then fails group certification;
sparse LA rejects uncovered columns before certification. The initial test
incorrectly expected no relation-check interval for either path. The corrected
test requires the observed distinction. Instrumentation was correct and the
fix changes only the test. `final/measured-library-source.rs.gz` preserves the
exact source used by the final worker panel before that test-only correction.

Final local validation passes 82 distinct release Rust tests (84 executions,
because the descent filter overlaps two accounting controls): seven exclusive
controls, four query controls, two LA controls, 64 descent regressions, the
all-target three-oracle equivalence control and six worker tests. Clippy
succeeds with existing warnings. The complete Python suite passes 143 tests
with two platform skips; the final source-preservation follow-up passes all
seven phase tests. All 39 site tests pass. Logs and file identities are in
`final/` and `controls.json`.

Use the pinned dependency lockfile from `ci/Cargo.lock`. From the repository
root:

```sh
cp research/ic_candidate_tournament_20260915/ci/Cargo.lock Cargo.lock
RAYON_NUM_THREADS=1 cargo test --locked --release --example ic_tournament_worker -- --include-ignored --test-threads=1
RAYON_NUM_THREADS=1 cargo test --locked --release --lib exclusive_ -- --include-ignored --test-threads=1
cargo build --locked --release --example ic_tournament_worker
python3.12 research/ic_candidate_tournament_20260915/goal_20260924/generic-exclusive-phases/run_controls.py \
  --worker target/release/examples/ic_tournament_worker --out /tmp/ic-exclusive-fresh
python3.12 -m unittest discover -s research/ic_candidate_tournament_20260915 -p 'test_*.py' -v
```

Choose a new output directory each time. The strict-session tests require an
isolated process; their default `ignore` attributes are explicitly overridden
by the commands above and CI. The local toolchain is Rust 1.93.1, Python 3.12.13,
macOS arm64; Linux CI pins Rust 1.94.1 and Valgrind 3.22. Exact source/dispatch,
base/matrix admission and qualification against the optimized incumbent remain
separate requirements before a tournament comparison.

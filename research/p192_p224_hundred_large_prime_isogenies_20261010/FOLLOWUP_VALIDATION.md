# Required library gate follow-up, 2026-10-10

The scientific campaign is unchanged: 100 independently verified maps across
two prime degrees, with 50 linked edges in each public P192/P224 class. The
frozen producer is source 49a280dc6; the standalone 449-model archive is source
3c39c5f23. A fresh native audit again matched all 891 files in MANIFEST.json.
The archived executable again passed a two-step chain replay with complete=true.

The previous full-library run stopped before a final summary. Its retained
root-complete-prefix.log records at least ten failures; it is not a passing
suite. The old temporary log and process were unavailable when this follow-up
started. The committed checkpoint and every original map receipt remain intact.

## Reproductions and repairs

The saved release test executable reproduced three degree-71 construction
failures and the isogenous-model factor-base failure in isolation. The
linearised-kernel elimination packed field images and preimages into one u64,
discarding coordinates above bit 63. Packing both words into u128 restores the
admitted wide-field range. A regression checks X^2+X, the identity and the zero
map at degree 71, including all high-coordinate basis vectors.

The packed orbit walk applied Frobenius even when the model's coefficients
were not fixed by it. Models without that endomorphism now use the existing
general orbit walk, whose step is the model's declared identity operation.
The existing fast point lifts remain available, and ordinary Koblitz models
retain their packed orbit path. Both existing isogenous-model tests pass.

The j=0 stage-accounting fixture failed on the second isolated repetition.
It demanded full relation rank with only two extra random rows. The fixture
now uses 32 redundant rows and checks the corresponding exact raw-row count;
the solver and its failure behavior are unchanged. This remains a stochastic
fixture, rather than a deterministic rank guarantee.

The three cover-solver failures, the Edwards four-point failure and the gate
benchmark failure passed in isolated runs. Their full-run failures may depend
on host contention, timing or shared process state; the isolated passes do not
prove which cause applied. Their arithmetic and acceptance tests were retained.
The full rerun uses two concurrent test threads with localhost access.

## Checks completed so far

- 57 wide-arithmetic tests pass, including the new high-coordinate regression.
- The degree-71 curve constructor and subgroup-order check pass.
- Both isogenous-model tests pass.
- The revised j=0 fixture passes 50 consecutive isolated repetitions; all
  individual test outputs are retained in one concatenated log.
- All 12 active cover-module tests pass in isolation; its existing ignored
  exhaustive test remains ignored.
- The Edwards four-point test and gate benchmark pass in isolation.
- The original artifact audit matches all 891 recorded files.
- The archived 449-model executable replays two saved P224/1117 steps.
- Formatting and git diff --check pass for the follow-up changes.

## Catalogue and all-target integration

The generated trait index contained 121 older records and an older registry
digest. The native curve_traits executable now derives all 449 current models
with the existing factor budget. Its fresh rebuild/check passes. The guide
identifies the older 121-model classification panel as a historical snapshot;
current records carry individual certificate statuses. A checked Curve25519
record explicitly reports not_evaluated because its coordinate/order adapter
is unavailable. The scientific artifact audit still matches all 891 original
manifest bindings after this additional generated index is updated.

The all-target check first encountered a sparse-checkout omission of the
committed Weil-composition contract. The exact tracked directory was restored.
The next check found malformed syntax in koblitz_rho_batch_ks. Comparing the
example to its last clean functional revision, 76869fd7e, showed exactly 36
stale added lines from a union merge. Removing those additions restores that
revision byte for byte. The restored example builds, is formatted, and passes
three known-answer fixtures in each of ordinary, negation-only and signed
Frobenius modes. These runs validate behavior; their timings are not a matched
performance result. The example has no unit tests, so the fixture runs provide
its runtime check.

The subsequent all-target check exposed 98 stale ic-binary API errors. Existing
upstream repair 3d8866c03 was cherry-picked unchanged as c36350601. It changes
only src/bin/ic.rs, src/bin/ic/experiment.rs and src/bin/hyperelliptic_cover.rs;
their contents match that upstream repair exactly. The final all-target and
binary-test checks are recorded separately from the preserved failed attempts.

Both repaired binaries pass their release tests: 6 hyperelliptic-cover tests and
8 ic tests. The example's cargo test target compiles successfully and contains
zero unit tests; the nine public known-answer smoke fixtures provide its
runtime validation. The all-target check also required materializing the exact
tracked prime-Fourier renderer modules omitted by sparse checkout. A remaining
SAT-example union residue is restored to the already merged 9e976658b revision.

## Complete library run and final calibration fixture correction

The complete two-thread library run finished with 3807 passed, 1 failed and
110 ignored tests in 2378.30 seconds. All earlier algorithm and integration
failures passed; the single remaining failure was the gate benchmark's
assertion that an elapsed-time enumeration/addition ratio exceeded 0.5.
The benchmark's operation counts use thread-local counters; this assertion
is a wall-clock performance floor, not a deterministic operation-count
invariant. Its isolated pass does not establish a scheduling-independent
lower bound for a contended full-suite run.

The test retains exact gate/oracle agreement, verdict totals and positive work
counts. It now additionally requires finite, positive calibration values,
including the enumeration/addition ratio. Production benchmark code and every
scientific performance artifact remain unchanged. Performance thresholds are
assessed in matched benchmark panels, separately from these correctness checks.
The focused correction and a subsequent complete library run are recorded
separately; the failed complete run remains preserved.

The final workspace compile passes after the remaining already merged repair
files are restored: 9e976658b's rank fixture, unknown-scalar panel and benchmark
file, and a2c157c22's compact-orbit example. The compact-orbit example also
receives formatting only. No new solver feature is introduced. Compile-time
CryptoPro and prepared-tournament fixture imports were restored byte for byte
from the tracked tree. Compiler warnings remain visible in the successful log.

The repaired compact-orbit and rank examples pass 3 and 12 release tests. The
SAT example passes 12 tests and exposes one missing tracked n53 fixture header;
the exact header directory is restored for its recheck. The restored example
source is unchanged by that fixture restoration.

The next complete library rerun reached another stochastic fixture failure in
ec_index_calculus::tests::staged_dlp_recovers_and_splits_stages, which passed
in the preceding complete run. It was stopped after the known failure; its
durable log and SIGTERM/exit receipt remain retained. It is not a passing run.
The ordinary accounting fixture now uses 64 extra random rows and retains its
exact row-count, known-answer verification and stage checks. Production solver
behavior is unchanged. Like the earlier j=0 fixture, it remains stochastic.
Both staged fixtures and the subsequent complete suite are checked separately.

With the certified n53 header restored, all four touched example targets pass:
3 compact-orbit tests, 12 rank-fixture tests, 13 SAT-instance tests and 2 public
point-input panel tests (30 total). Their successful terminal receipt is retained
separately from the earlier missing-fixture failure. The all-target workspace
check passes with the exact static imports materialized. The ordinary staged
fixture's focused check passes. Both staged fixtures then pass 100 consecutive
paired invocations (200 successful fixture executions). Their complete outputs
are retained. This repetition is a regression check, not a deterministic rank
guarantee. The complete corrected-library run is recorded separately.

## Passing committed-version closeout

The complete release library run on 6742ecd0b37d70e02706214ed8f51368305bc37e
finished successfully: **3808 passed, 0 failed, 110 existing tests ignored**,
in 1660.57 seconds. It used two concurrent test threads and localhost access.
All tests were included except the suite's existing ignored tests; no new
skip was introduced. The terminal log and exit code 0 are retained as
root-committed-final-20261010.log and .exit.

The same committed version passes cargo check --workspace --all-targets.
Its terminal output and exit code 0 are retained separately. The 30 touched
example tests, 14 repaired-binary tests, generated 449-model trait freshness
check, 200 paired staged-fixture executions, dedicated isogeny checks and
original 891-file artifact audit provide the scoped validation described above.
Hosted CI remains a separate gate for the PR and ReleaseMe publication.

This closeout supersedes the pending library/push status in the frozen campaign
report. The scientific maps, timings, source/binary archives and original
manifest retain their historical bytes.

## Main integration and frozen canonical context

After the passing branch was pushed as 12a846a56, GitHub reported conflicts
against main 8b6b1384e. The release conflict concerns equivalent cancellation
comments; both sides use the same cancellation predicate. The merged workflow
retains the dedicated isogeny platform jobs and pinned ReleaseMe publication.
The compact-orbit example retains main's current code with formatting applied.
Main's committed research, ledger and browser additions are preserved.

Main changes the current browser context. Before integration, all ten
canonical outputs bound by the original scientific manifest were archived
from clean 12a846a56 as evidence/canonical-at-run-12a846a56.tar. The original
native audit runs in that extracted canonical context with the same sealed
study files and passes all 891 bindings, retaining manifest hash
3b7ba9aeac62495be74d4e6f7e5089af12834e85eb5c968b6aa1bedaed5558b5.
The receipt is evidence/frozen-canonical-context-audit.json. The original
scientific manifest is a frozen-context manifest; it is not a claim that a
later global browser file still has its earlier hash.

The merged browser refresh uses the native canonical-roster helper, which
asserts that only its four source hashes change. Incoming browser candidate
records and every other value are preserved. Workflow Actionlint and example
formatting pass. The merged-version library, workspace and example checks
are recorded separately from the pre-merge passing checks.

The merged workspace all-target compile and three compact-orbit example tests
pass. The merged library build includes 3925 tests. Its first run exposed the
S4 staged-format fixture's six-extra-row rank variance; that test passes in
isolation, and the failed run was stopped with its SIGTERM/exit receipt retained.
The S4 fixture now uses 64 extra rows in both compared formats and requires
the exact raw-row count. Its mathematical verification assertions and production
solver remain unchanged. The focused repetition and subsequent full run are
separate checks; the interrupted run is not presented as passing.

The corrected S4 fixture passes its focused release test and 100 consecutive
isolated repetitions. Every output is retained in s4-stage-100-repetitions-20261010.log.
The full merged-tree run follows this committed correction. This is additional
fixture validation and does not change the frozen map or performance evidence.

## Passing merged-source closeout

Committed merged source bff6dcb274b6fc29c46d49a5603f5bb801cc2db5 passes the
complete release library suite: **3815 passed, 0 failed, 110 existing tests
ignored**, in 1926.77 seconds. The command includes every active library test
and uses two concurrent test threads with localhost access. Its stdout/stderr
and exit code 0 are retained separately from all earlier failed and interrupted
runs. The committed merged source also passes the workspace all-target compile.

The release workflow passes Actionlint, the resolved example passes three
release tests, and the original native artifact audit passes 891 bindings in
the preserved canonical-output context. Only documentation and terminal
receipts are appended after this validated source revision. Hosted CI and
ReleaseMe publication remain separate gates.

The full command is:

```sh
cargo test --release --lib --target-dir target/isogeny-root-gate-20261010 -- --test-threads=2
```

Its durable log is evidence/validation/root-repaired-full-20261010.log and its
terminal status is written to root-repaired-full-20261010.exit. Until that
command finishes successfully, the required library gate remains pending and
the new commits are not pushed. Hosted CI and ReleaseMe publication follow
their existing gates. These follow-up logs are additional integration evidence
bound by their Git commit; the frozen scientific manifest is not regenerated.

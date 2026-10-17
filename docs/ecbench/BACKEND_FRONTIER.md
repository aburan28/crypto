# Backend leaderboard and improvement loop

Open the [interactive backend leaderboard](../backend-leaderboard/index.html),
its [Markdown table](../backend-leaderboard/LEADERBOARD.md), or its
[machine-readable data](../backend-leaderboard/leaderboard.json). The
[IC recipe leaderboard](../ic/LEADERBOARD.md) carries the admitted pipeline
comparisons, including evidence outside this report's explicit session manifest.

This report puts F4 matrix, inherited F4, F5 matrix, F6-IC, custom root finding,
SAT and reference methods in one view. Each row retains the complete combination:
method and implementation, oracle, solver options, factor-base identity, linear
algebra, exact curve/EC1/UID, subgroup, public target, algorithm seeds and resource
conditions. Engine labels come from configuration, not an arm's display name.
The report is a standalone native Rust consumer of `ecbench` sessions; it runs
no solvers and builds independently of the cryptography library.

## What the frontier means

An observed count frontier contains the entries for which no other eligible arm
on the **same exact target in the same session, with the same rounds and seeds**
has strictly lower mean S at every weight 0, 0.1 and 1 GAE per native operation.
It is a finite-workload sensitivity result. Equal or crossing bands remain
unresolved. The count/memory frontier also keeps process-RSS tradeoffs visible.
The raw native counters keep their own units: a conflict, a learned clause,
a word XOR and a table lookup are not asserted to cost the same amount.

The **IC recipe frontier** applies the same comparison within IC arms only.
It identifies the incumbent an IC agent should improve even while a reference
method remains cheaper. The JSON retains both `ic_count_frontier` /
`ic_count_dominated_by` and the cross-method `count_frontier` /
`count_dominated_by`; selecting a backend does not silently recalculate either
boundary against a smaller set of competitors.

The tool checks file hashes, raw record seals, sequence, exact plan/record
identity, units and checked-answer status. It does not reconstruct a curve,
rerun an arithmetic verifier or execute a solver. Admission also requires a
hash-bound historical native `ecbench` audit receipt: that audit checked the
arithmetic and derived costs of the exact files. Missing counts, absent audit
receipts, incomplete plans, failures, nondeterministic charged work
and invalid costs block count ranking. A historical replay receipt must bind all
five session files by hash before its reproduced record IDs can supply a replay
badge. This badge says what the old receipt checked; it is not a fresh or
independent-host performance replay. Audit failures remain visible.

Older counted SAT sessions placed their solver observation under
`detail.decomposition.solver` while the top-level `solver` field was null.
For explicitly unpriced solver work the read-only projection exposes that
recorded observation to the shared band code and records its JSON pointer as
`solver_origin`. It does not replace a clock-based charge or invent a price.
The original decomposition/system/solver and calibration context is retained,
including original/learned clauses, decisions, propagations, budget outcomes and
legacy priced fields. Bulky auxiliary tables stay in the hash-bound source
records instead of being duplicated in the report.
Native solver work, SAT conflicts and learned clauses also appear directly in
the table; absent counters remain unknown.

Cold solve time covers the method's setup through recovery and internal
verification. Process time is fork-to-reap, including startup. The runner's
independent scalar check is outside these clocks. Online time excludes reusable
setup and appears only when every recorded window has a matching phase sum.
The JSON retains exclusive phase clocks, the named boundary events and all
attempts; absent clocks stay null. A budget-limited or failed solve is not a
completed end-to-end measurement.

A timing verdict requires every paired attempt to verify, every round to meet
at least L2, at least five matched rounds, and a same-target A/A control interval.
Its paired 95% bootstrap interval must exclude parity and the A/A noise interval.
Intervals resample rounds within **one target**, using 2,000 resamples and fixed
seed 20261017; they do not infer performance on other targets or growing sizes.
Raw wall observations remain readable when these gates fail.

## Field sizes and coverage

The manifest explicitly asks for extension degrees 3, 4, 5, 31, 51, 53 and 83,
and the report includes every other measured degree too. Thus fields of size
2³, 2⁴ and 2⁵ are shown separately from degrees 51, 53 and 83. Subgroup order r
is a different field. A missing cell means **no imported measurement**, not an
algorithmic impossibility or an unsupported-security claim. Planned arms without
measured records are retained separately; warmups do not count as measurements.

Degree 51 has proper intermediate subfields GF(2³) and GF(2¹⁷). Degree 31 has a
different Frobenius-module structure from degrees 53, 83 and 131. A result that
uses either distinction cannot silently transfer between these sizes. Degree 83
on the exact E₀ curve, field representation and subgroup remains the smaller
confidence gate for ECC2K-130. A degree-83 rho-only calibration is not coverage
of F4, F5, F6-IC, SAT or root finding there. The `ecbench` wide-curve method limits
are an implementation gap, not a mathematical obstruction; see the harness
[limits](README.md#12-limits) and the linked IC recipe evidence.

## Rebuild and check

From the repository root:

```sh
cargo build --locked --release --manifest-path tools/backend-leaderboard/Cargo.toml
tools/backend-leaderboard/target/release/backend-leaderboard \
  --manifest docs/backend-leaderboard/manifest.json \
  --out docs/backend-leaderboard
tools/backend-leaderboard/target/release/backend-leaderboard \
  --manifest docs/backend-leaderboard/manifest.json \
  --out docs/backend-leaderboard --check
```

In a sparse checkout, materialize the session directories and receipt files
named by the manifest with `scripts/sparse-checkout.sh add --path ...` first.
The generator does not silently skip an unreadable input. `--check` never repairs
an output. The JSON, HTML and Markdown come from one build and must be committed
together. The tool directly shares the harness's native pricing-band, bootstrap
and identity-codec source, and extracts its native counter vocabulary at build
time. Calibration and curve-alias inputs are compiled in and hashed in the
report. The HTML is standalone and works offline; the published route is
`/backend-leaderboard/`. Source links resolve to the repository, with the exact
content hashes in the data beside them.

## An agent's next iteration

1. Select an exact target and incumbent in `leaderboard.json`. Read
   `count_dominated_by`, `count_blockers`, native solver counters, cold phases,
   replay status and the source session. A cell with missing evidence first
   needs that evidence; it has no established winner to beat.
2. Freeze a protocol before running: baseline binary/source, workloads and
   hashes, output checks, seeds, resource limits, operation accounting, timing
   boundaries, A/A and positive controls, holdouts, success and stop conditions.
   Choose an existing session's spec as a starting point, writing a **new** spec;
   do not change frozen evidence. Keep the requested size matrix visible when
   a constructor, backend or available memory cannot complete a cell.
3. Add a challenger changing one lever: F4/F5/F6, root finding, SAT XOR encoding,
   conflict or learned-clause policy, solver budget, factor base, or a combination.
   For a direct solver comparison hold the factor base and query stream fixed;
   for a recipe comparison price the changed construction and setup as well.
4. Run the baseline and candidate on identical single-target workloads and fresh
   holdouts through `ecbench run`, under the repository isolation rules. Use at
   least five interleaved rounds and same-target A/A controls for timing. Keep
   every timeout, OOM, failed verification and incomplete collection in the run.
   Save the native `ecbench resources`, `band` and paired `compare` reports too.
5. Verify with `ecbench verify --replay-all`, independently replay on another
   host, and compare the fastest available external tool on the same inputs and
   host. For Gröbner systems export with `gb_cas_export`; identify Magma/F5,
   msolve, FGb and other unrun references explicitly. For roots use the matching
   FLINT/NTL/Magma computation. No external result is inferred from our own arms.
6. Accept only the preregistered full-operation improvement with preserved
   correctness, complete costs, matched references and holdouts. Fewer SAT
   conflicts or clauses is a stage result until it moves that total. Update the
   manifest, generated views and affected canonical IC panels in the same PR.
   Keep the prior baseline as evidence, then name the next measurable target.

Typical commands for a **new frozen spec and output directory**:

```sh
target/release/ecbench plan --spec NEW-SPEC.json --json
target/release/ecbench run --spec NEW-SPEC.json --out NEW-SESSION
target/release/ecbench verify --dir NEW-SESSION --replay-all --exit-code
target/release/ecbench resources --dir NEW-SESSION --out NEW-RESOURCES.json
target/release/ecbench band --dir NEW-SESSION --ranks --json NEW-BAND.json
target/release/ecbench compare --dir NEW-SESSION --a BASELINE --b CANDIDATE --json
```

An `ecbench bound` / `frontier` / `challenge` can additionally freeze a fitted
cross-size objective; see [bounds and challenges](../bounds/README.md). This
backend view does not replace those protocols or promote a ledger row.

## Requested scope and evidence

| Requirement | Status | Implementation / evidence |
|---|---|---|
| Backend and combination leaderboard | Verified complete | Native tool, generated views, full configurations and exact-target rows |
| F4, F5, custom F6 and root finding | Verified for imported sessions | Explicit manifest; SAT and reference arms included; historical source/binary hashes retained |
| Field-size ladder, including small powers and 51/53/83 | Grid complete; measured coverage partial | Absent and unattempted cells remain visible; 3/4/5/51/83 have no imported measurements |
| Establish an auditable frontier | Verified for the imported panel | IC and cross-method count-band frontiers, memory tradeoffs and timing gates; current/global fastest remains unestablished |
| Guide agents toward improvements | Implemented and validated | Frozen incumbent/challenger loop, counters, failure records, replay, external-tool requirements and CI freshness checks |
| End-to-end operation time | Historical coverage partial | Cold, online, process and exclusive phases retain their boundaries; no admitted timing winner or fresh campaign |

This change reindexes historical observations. It does not run a fresh performance
campaign, repair wide-field solver implementations, complete absent degree-51/83
backend comparisons, or establish an end-to-end timing winner where the evidence
does not meet the gates. These open cells remain in the report for subsequent
experiments. The existing IC leaderboard and progress chart are unchanged because
no new empirical measurement or ledger promotion is introduced.

Conductor was checked before editing but its CLI and `.conductor/WORKFLOW.md`
are absent in this workspace. Work is isolated on a feature branch; no
reservation is represented as obtained.

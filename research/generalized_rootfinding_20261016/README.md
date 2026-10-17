# Generalized rootfinding: frozen research checkpoint

This checks in the previously local source and evidence for the 2026-10-16
research cycle. **Historical Python development calibration; no qualified
native, end-to-end index-calculus, or F4/F5/F6 performance claim.** The requested
5x target and the unchanged 2x release gate remain open.

The 112 original payloads are preserved under `historical/`. Files larger than
1 MB use deterministic gzip transport; decoding restores the exact original
bytes. `TRANSPORT.json` records both transport and original SHA-256 digests,
lengths, and paths. No frozen algorithm, input, measurement, or failure was edited.
The main source modules are:

- `historical/rootfinder/`: finite-field arithmetic and exact bounded solving.
- `historical/research/fast_curve_fibers_20261015_v1/`: the prior cached reference.
- `historical/research/separator_batch_20261016_v1/`: shared complete left table.
- `historical/research/support_cost_20261016_v1/`: repeated-root support criteria.
- `historical/lab/ledger/attempts/`: all three current-cycle raw attempts.

These snapshots are evidence, not a new Python execution path for this repository.
Further solver, oracle, harness and performance work must use a native port,
with matched fresh measurements. The source snapshot's older README/skill and
runner commands describe historical tooling, not an exception to that requirement.

## Recorded observations

| Historical comparison | Workload and samples | Recorded ratio, prior/candidate | Scope |
| --- | --- | ---: | --- |
| Cached separator / shared-left separator | 20 fields, 38 curve contexts, 3 cold pairs; source domain at most 4 | 1.3830038111 | Whole-cold calibration panel; isolation unqualified |
| Nilpotence power / squarefree shortcut | GF(2^83), 16 fixed-parameter relaxed quintics, 3 internal rounds | 1.1408942780 | Support verification, not search or roots |
| Nilpotence power / squarefree shortcut | P192/P224/P256/P384/P521, mixed repeated-root samples, 3 internal rounds | 0.5948–0.7034 | Blanket shortcut regresses; no adoption |

The separator panel covers all five requested NIST prime fields, binary
extensions through degree 571 including 83, Koblitz and ordinary binary curve
models, and odd extensions. This is finite calibration coverage, not all fields,
all curve models, or all scales. Outputs retain the complete original eight
coordinates. Support experiments classify supplied quintics and leave original
root counts unset. None of these comparisons uses a genuine matched Gröbner
baseline. No IC/rho cost or ratio is available.

Historical independent oracle algorithms passed 15 complete small-system cases
and 4,448 support cases (characteristics 2, 3, 5, repeated/inseparable roots).
That is sole-model work with independent algorithms, not external review.
The written algebraic contracts are archived, not newly proved or Lean-certified
by this PR. See [the result decision](RESULT_DECISION.md),
[results PDF](RESULTS.pdf), and [run reports](reports/run-20261016T0125Z/REPORT.md).

## Native archival verification

From the repository root:

```sh
cargo test --manifest-path research/generalized_rootfinding_20261016/Cargo.toml --locked
cargo run --release --manifest-path research/generalized_rootfinding_20261016/Cargo.toml --locked
```

The Rust auditor verifies all 112 decoded payload hashes and lengths, the two
pinned dependency closures, and 190 recorded context-output comparisons across
three cold pairs. It never invokes Python. This is transport/recorded-output
verification, **not independent equation solving or performance qualification**.
Gzip decoding is bounded at 512 MiB per payload. The auditor includes path-escape
and tampered-transport negative controls and the existing native SHA-256 tests.

To inspect a compressed raw result without modifying the snapshot:

```sh
gzip -dc research/generalized_rootfinding_20261016/historical/lab/ledger/attempts/batch-separator-manyfield-exact-coldpairs-20261016-v1/0001/scientific_result.json.gz > /tmp/rootfinding-record.json
```

The full older multi-gigabyte journal is identified by its frozen digest, not
included. Older attempts and native runtime dependencies are not all bundled.
This checkpoint is therefore not a complete resumable checkout of the original
lab. Every current-cycle raw record and pinned source/input payload is included.

## Requirement status and next work

| Requirement | Status and remaining condition |
| --- | --- |
| Check in tested code and current evidence | Included, subject to PR review/CI |
| Native archival audit | Implemented; receipts in `validation/` |
| General prime/extension exact solving | Partial bounded prototypes; native port pending |
| Actual hard GF(2^83) all-roots search | Unfinished; the explicit left-table bound is 2^56 |
| Qualified >=2x and requested >=5x over genuine matched F4/F5/F6 | Unmet; preserve every condition of `historical/performance/GATE.md` |
| Novelty and formal proof | No priority claim; primary-source review and full formalization remain open |

Next: port the solver and complete-output oracle to Rust before new comparisons;
test an inversion-free squarefree criterion on a frozen stratified holdout; and
profile sparse contraction width before adding two F8 places to the hard support
search. Do not enumerate 2^73 parameters, or materialize a dense 2^44 transform.

This archival addition creates a separate descriptive dashboard panel and a
benchmark browsing entry. The IC/rho progress timeline, best-of leaderboard,
boundary facts and exponents are unchanged because no admissible IC/rho,
pipeline-cost, or exponent measurement was added or superseded.

Conductor was attempted before editing, but its CLI is absent and
`.conductor/WORKFLOW.md` is not present at the pinned base. The exact diagnostic
is retained under `validation/`; no reservation is claimed.

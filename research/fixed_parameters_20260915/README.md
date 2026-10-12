# Fixed-parameter persistent IC

## Restore the research evidence

Large text artifacts are stored as readable chunks so each publication request
fits the automatic review limit. Restore their exact bytes before training on
the saved panels or rerunning an audit:

```sh
python3 research/ic_solver_selection_20260915/restore_evidence.py
```

The script verifies every artifact's SHA-256 hash and refuses to replace a
changed file. Use `--verify` to validate chunks without writing files.

This extends the existing arbitrary-width ONB arithmetic with a fixed-parameter
campaign and exposes it as `ic fixed`. It accepts the actual ECC2K-130 polynomial
basis, generator and target, persists pair tables and verified relations, and
supports both target-independent precomputation and direct target equations.

Classification: **engineering capability**, with accounting diagnostics. The
mathematical counting boundary is unchanged. A completed cached target requires
no new pairs or probes, but all original cold work remains part of its cost.
No operation-cost, runtime-speedup, exponent or practical degree-131 completion
claim is made. The old Rust field limit is unchanged; the new command calls the
CPU Python engine.

The [contract](contract.json) was recorded before measurement. The equivalent
suite exception applies because this persistent ONB pair-table interface is not
the frozen WDSat ANF interface. It compares matched complete ephemeral, durable,
and interrupted/resumed workloads and includes fresh holdouts. Every witness is
verified and every completed scalar satisfies the group equation. Negative pair
queries are independently checked by exhaustive group enumeration on small fields.

See [usage and limitations](../../docs/ic/FIXED_PARAMETERS.md).

```sh
/path/to/python ecc2k130/codegen/test_indexcalc_fixed.py -v
cargo test --release --test ic_fixed --test ic_framework --test ic_progress
/path/to/python research/fixed_parameters_20260915/run.py --out research/fixed_parameters_20260915/run-01
```

Each evidence output directory must be new. SQLite campaign files are temporary
in this regression suite; raw relation transcripts, inputs, results and accounting
are preserved in its JSON artifacts. The actual usage commands keep their SQLite
files across runs.

## Saved validation

The [accepted comparison](run-01/summary.json) contains 72 matched cases:
24 parameter/seed combinations with three repetitions each. Each case has an
ephemeral reference, a durable cold run, and a staged resumed run: **216 complete
computations and 432 target verifications**, with identical relation transcripts
in every matched case. These counts include repetitions; there are 48 target
instances across the parameter/seed combinations. Fresh holdouts use seeds 1009,
1013 and 1019. All 72 completed-result replays built zero new pair candidates
and performed zero new probes.

One primary unit throughout: calibrated common operations, currently unmeasured.

| Variant | Common operations | S | Ratio to rho | Ratio to counting floor | Verified targets | Class |
|---|---|---|---|---|---|---|
| Ephemeral reference | null | null | null | null | 144/144 | engineering reference |
| Durable cold | null | null | null | null | 144/144 | engineering capability |
| Staged resume | null | null | null | null | 144/144 | engineering capability |

Supplementary complete API times across 72 runs are 0.486248362 seconds for the
ephemeral reference, 8.659221777 seconds for durable cold, and 13.023860283 seconds
for staged resume. Persistence and revalidation cost work on these tiny examples;
these figures are diagnostics on a shared host, without paired confidence
intervals or a speedup claim. Each API time includes report archival and database
close. All intermediate invocations and cold costs remain in `raw.jsonl`.

[Degree-131 evidence](run-01/degree131.json) validates the actual fixed parameters,
14 signed orbit columns and 3668 signed points. Two invocations commit 11 and 13
new pair candidates without repeating the first chunk. A bounded natural-target
SAT attempt saves its outcome and no relation; its status is correctly
`insufficient_relations`. No full-size scalar is reported. Its uniform-target
coverage ceiling is 8231743940 / 680564733841876926932320129493409985128.

Ten Python tests pass, including exhaustive small-group decomposition checks,
parameter/basis validation, direct target equations, fresh-target reuse,
transaction rollback, corruption rejection, SAT and writer exclusion. Twenty-two
Rust integration tests pass: two new CLI tests, 16 existing framework tests and
four existing progress tests. `cargo check --bin ic` also passes; existing
unrelated compiler warnings remain. See [validation](validation/summary.json).

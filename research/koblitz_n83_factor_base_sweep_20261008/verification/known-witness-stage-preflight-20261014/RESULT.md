# Source-only staged successor to the ordered N83 control

The prior [fixed-witness control](../known-witness-ordered-sat-20261014/RESULT.md)
hit its external wall cap with no native receipt, leaving model construction
versus SAT search unresolved. The new
[`koblitz_n83_fixed_witness_stage_probe` example](../../../../examples/koblitz_n83_fixed_witness_stage_probe.rs)
retains its public ordered witness and full-domain input path while emitting
flushed JSON-line markers at process entry, input validation, core build,
domain installation, first-pair ordering, fixed-unit installation, solver
entry, solver return, and model check. Each post-input marker binds the input
by BLAKE3 and reports monotone wall time; construction markers report SAT
variables, clauses, and XOR rows. A future externally capped run can therefore
identify the last completed stage even if it emits no terminal result.

This is **source preparation and marker transport only**. The release example
[compiled and passed its test harness](example-release.log.gz) (zero unit tests),
then `self-test-markers` emitted two ordered markers and a `PASS`
result without reading the N83 domains or building a solver. The
[self-test receipt](selftest-audit.json) pins the source and binary SHA-256s
and both raw [stdout](selftest.stdout) and [stderr](selftest.stderr).
The [trace checker](check_trace.py) replays that receipt. Six focused tests
exercise prefix-only wall-cap classification and fail-closed handling of
malformed order, changed input identity, and a claimed UNSAT encoding result.
No N83 solver search, natural relation collection, rank, or target logarithm
was performed, and the conservative local ledger remains **22,200 seconds**.

Reproduce from the repository root:

```sh
cargo build --release --example koblitz_n83_fixed_witness_stage_probe --offline
target/release/examples/koblitz_n83_fixed_witness_stage_probe self-test-markers
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/known-witness-stage-preflight-20261014/check_trace.py --check
python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008/verification/known-witness-stage-preflight-20261014 -p 'test_*.py'
```

A later full-domain run needs a separately frozen grant, a source-matched
Linux build, the pinned K1182/K2048 input, one-CPU/4-GiB/zero-swap isolation,
an external watchdog, and an audit that preserves every raw marker and process
exit. The marker `solver_entered` is emitted immediately before `solve()`;
seeing it is evidence of reaching that boundary, while a cap remains
censored. This preflight neither selects a factor base nor supplies a
complete one-target runtime comparison.

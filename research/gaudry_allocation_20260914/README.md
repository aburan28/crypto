> Frozen experiment archive: the source change is preserved in `optimization.patch`;
> the benchmark harnesses are saved under `harnesses/`. This UI PR does not apply
> that runtime patch.

# Gaudry elimination and normal-form allocation optimization

**Engineering:** the frozen comparison finds lower runtime with exactly the
same solver outputs and field-operation counters. No arithmetic, relation-yield,
normalized-cost or asymptotic improvement is claimed.

## Change and boundary

`src/cryptanalysis/gaudry_cubic.rs` contains three related overhead reductions:

* Forward elimination borrows the normalized pivot row using disjoint slices
  instead of allocating and cloning that row for every pivot.
* Normal-form recursion computes dependencies into the memo table, then borrows
  their slices. Cache hits and returns no longer clone quotient vectors.
* Dense column IDs index vectors of pivot/standard indices instead of performing
  repeated hash lookups in back-substitution.

The arithmetic loops, operation-counter increments, dependency order and missing
normal-form behavior are preserved. Failed border lookups remain uncached, so
retry/fallback work is counted exactly as before. This changes representation
costs; the measured field-operation ratio to the unmodified reference is **1.0**.
No generic-group boundary moves. The change is recorded in
[optimization.patch](optimization.patch).

[contract.json](contract.json) was frozen before the candidate was measured.
Its target requires identical solutions/counters and a paired 95% interval below
one before claiming a runtime improvement on the measured workload. This uses
the equivalent-suite exception in `AGENTS.md` §8: the Gaudry `F_(p³)` symmetrized
S4 solver does not implement the binary-ANF WDSat protocol. The parent common-cost
contract still applies. It does not alter or replace the frozen WDSat results.

## Correctness and timing

Two isolated source trees used identical frozen Cargo manifests, locks and build
commands. Only the Gaudry implementation differs. This prevents concurrent edits
to other solvers and dependencies in the shared checkout from entering the
comparison. Both variants use the same new residual harness and the unmodified
existing full-DLP harness. Timing subprocesses are pinned to CPU 0 on this host;
see the recorded affinity in the provenance file for the authoritative value.
Variants run in randomized pairs on the same inputs, with three repetitions.

The residual corpus has 60 targets in each of four prime/seed cases: (67,11),
(271,1), and fresh holdouts (67,401), (271,503). Every complete decomposition set
matches the other version and the independent meet-in-the-middle oracle. As in
the existing Gaudry test, MITM's two-term relations are excluded from the S4
triple-set comparison; every returned decomposition is independently checked by
curve addition. Complete baseline/candidate outputs, including exceptional
results and all arithmetic counters, are compared without that filtering.

The cold full-DLP corpus uses p=67 and seeds 11,401,503, including both holdouts.
Each process builds the instance and factor base, collects relations, solves,
recovers the planted scalar, verifies it, and runs matched rho. The cross-check
flag audits **every residual** with MITM, including empty answers. All 18 scalar
recoveries and 3,234 residual cross-checks passed. All non-timing full reports
match exactly, including relation statistics, solve counters and rho results.
The stage corpus adds 1,440 residual checks (240 distinct inputs, six replays).

One unit in the supplementary timing table: **seconds of whole child-process
wall time**. This includes cold instance/factor-base setup, allocation, oracle
construction, collection, elimination, result verification, output, process
startup, and the independent MITM/rho checks. The harness provides no warm reuse.
Ratios and intervals use paired log time ratios, clustered by input (all three
repetitions of an input stay together). They are not ratios of displayed totals.

| Variant | Workload | Sum seconds | Paired geometric time / baseline | Paired cluster-bootstrap 95% interval | Field counters / baseline | Correctness | Class |
|---|---|---:|---:|---|---:|---|---|
| Unmodified | Residual corpus, 12 runs | 3.178126 | 1.000000 | reference | 1.0 | 720/720 residuals | engineering reference |
| Borrowed rows/cache | Residual corpus, 12 runs | 3.013655 | 0.945600 | 0.934631–0.956698 | 1.0 | 720/720 residuals | engineering |
| Unmodified | Cold full DLP, 9 runs | 5.500866 | 1.000000 | reference | 1.0 | 9/9 scalars | engineering reference |
| Borrowed rows/cache | Cold full DLP, 9 runs | 5.113203 | 0.929583 | 0.928868–0.930822 | 1.0 | 9/9 scalars | engineering |

Evidence: [comparison including every pair](run-01/comparison.json),
[raw commands, status and timings](run-01/raw.jsonl), and
[source/binary hashes and dependency provenance](run-01/provenance.json).
There were 42 successful runs, zero missing results, and zero output/counter
mismatches. Each per-run JSON is retained beside these files.

These observations support approximately 5.4% less elapsed time for the residual
corpus and 7.0% less for the tested cold full-DLP workload. The full-DLP interval
uses only three input clusters and the stage interval only four; neither supports
extrapolation to larger fields, different hosts or other index-calculus engines.
The whole-process time includes audit overhead. Perf sampling was unavailable
because of host permissions; selection was based on the source's repeated
allocations and the baseline's existing phase counters, not a sampled CPU profile.

## Operation accounting limits

The existing report converts counted Fp multiplications using 63 multiplications
per affine group addition. Those reported solver S values and ratios to matched
rho are identical before/after: S=3317.53, 6054.55, 3837.07 and ratios=816.30,
1201.32, 827.41 for seeds 11,401,503. These are **existing partial solver metrics**;
they do not newly price cold instance/factor-base setup, allocation or bookkeeping.
The complete common-operation total, S, cost/rho and cost/floor therefore remain
null for this iteration. A runtime improvement is not a total-operation or
cryptanalytic advance. No new scaling exponent or rho crossover is claimed.

## Tests and reproduction

All 14 Gaudry unit tests pass in the isolated candidate, including the 600-input
border retry comparison, MITM equivalence, dense/sparse/large-prime scalar
recovery, and new quotient/cache tests. The new tests include dependent rows,
free leading columns, rectangular matrices and unavailable border columns.
[Full test log](gaudry-tests.log).

Build both isolated copies using `cargo build --release --locked --example
 gaudry_cubic_bench --example gaudry_allocation_stage`. The reference uses the
unmodified Gaudry source at repository commit
`01fa8426d963fefcb8ac048783863f15d9352d91`; copy the new stage harness into both.
Use the saved `run-01/cargo-manifest.toml` as Cargo.toml and
`run-01/cargo-lock.txt` as Cargo.lock to reproduce dependency versions. Apply
`optimization.patch` only to the candidate.

```sh
python3 research/gaudry_allocation_20260914/run.py \
  --baseline /path/to/baseline --candidate /path/to/candidate \
  --output /tmp/gaudry-fresh-run
cargo test --release --locked --manifest-path /path/to/candidate/Cargo.toml \
  --lib cryptanalysis::gaudry_cubic::tests -- --test-threads=1
```

The runner refuses an existing output directory and preserves failed/timed-out
processes. A missing completion or mismatched result/counter rejects the run.

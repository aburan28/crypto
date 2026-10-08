# Stage 23: public unknown-scalar IC versus signed-Frobenius rho

Stage 23 turns the degree-23 relation-yield work into a five-target end-to-end
toy panel. Targets are generated as public subgroup points by a
domain-separated BLAKE3 hash-to-affine stream followed by cofactor projection.
No target scalar is constructed, stored, or passed to either solver process.

The index-calculus arm uses the frozen target-independent `K_0` divisor
`[0,2]` factor base. Factor-base logarithms and the target scalar are learned
only from collected relation rows and modular linear algebra. The matched rho
arm receives the exact same public point and uses the signed-Frobenius/negation
quotient walk with a fresh deterministic 16-entry jump table per restart.

## Producer modes

The Rust producer has three process modes:

```text
koblitz_unknown_scalar_panel targets production
koblitz_unknown_scalar_panel ic production TARGET_ID X Y IC_SEED
koblitz_unknown_scalar_panel rho production TARGET_ID X Y RHO_SEED
```

`targets` emits five unique public points and independent IC/rho seeds. `ic`
retains one record for every generated relation candidate, including every
refutation, conflict-capped `Unknown`, invalid model, and skipped direct
relation. It also retains the learned relation matrix, rank history, modular
solve attempts, recovered scalar when present, and exact progress events.
`rho` retains setup, walk, restart, collision, canonicalization, Frobenius,
negation, partition, and point-verification counters.

The `smoke` profile uses two public targets on `K_1 / F_(2^7)` with divisor
`[0,2]`. Smoke output is operational only and is never scientific evidence.

## Production plan

Generate the exact write-once plan from a clean checkout:

```text
python3 scripts/run_koblitz_unknown_scalar_panel.py plan \
  --profile production \
  --output /Volumes/SSD990/koblitz-stage23-unknown-scalar-production-20260910-01
```

The plan prints an outer `process_meter.py` command. Run that exact command so
the whole driver encloses the fresh locked build, two public discoveries,
target generation, and five IC/rho pairs. The production task order is fixed:

1. one fresh `--locked`, single-job release build;
2. public factor-base discovery for `a=0` and `a=1`;
3. one hash-to-subgroup target-panel process;
4. five IC/rho process pairs in target order.

Every child has a fresh process receipt. The driver uses one worker and requests
one thread. Nonzero exits, watchdogs, missing logs, rank failures, and exhausted
rho walks are retained as incomplete and remain inconclusive. A process whose
meter terminated an orphan process group is also incomplete even if its direct
return code is zero. They are never converted to negative results, retried in
place, or omitted from accounting.

Verify a terminal run with a new output directory:

```text
python3 scripts/run_koblitz_unknown_scalar_panel.py verify \
  --run-root /Volumes/SSD990/koblitz-stage23-unknown-scalar-production-20260910-01 \
  --outer-metrics /Volumes/SSD990/koblitz-stage23-unknown-scalar-production-outer-20260910-01/driver.metrics.json \
  --output /Volumes/SSD990/koblitz-stage23-unknown-scalar-verification-20260910-01
```

## Completion and evidence boundary

One row is complete only when both algorithms return a scalar, independently
verify `[d]G = Q` for the exact public target, and agree on the recovered value.
The IC row must retain every attempted target, the admitted relation rows,
matrix dimensions and rank, modular linear-algebra history, zero invalid models,
and no direct-relation bypass. An IC process that reaches its trial cap still
computes, charges, and retains its terminal matrix rank before returning an
incomplete result. The rho setup ledger must balance against the number of fresh
jump tables. Rho walk timing stops before collision-candidate verification, so
setup, walk, and verification timing intervals do not overlap. Any other row is
incomplete and inconclusive.

The relation matrix has one column for each projected signed-Frobenius orbit
logarithm and a final column for the target scalar `d`. Completion does not
require every factor-base logarithm to be unique. The verifier independently
row-reduces the retained matrix modulo the prime subgroup order and requires
the affine solution space to have one invariant target coordinate: equivalently,
every homogeneous nullspace vector has zero in the final coordinate. It then
requires that coordinate to equal the reported scalar and retains the exact
`[d]G = Q` point check. A deficient total rank is therefore admissible only
when its free directions affect factor-base logs but cannot change `d`.

The project verifier reconstructs the exact permitted task prefix from the
archived protocol and target panel. For every task it reopens the intent, raw
meter output, receipt, stdout, stderr, inputs, and parsed result; checks the
command, watchdog, environment, and file locations; then recomputes rows,
resources, ratios, critical-failure state, and panel status. Updating hashes and
re-sealing a coherently altered command or summary is therefore rejected.

Panel completion requires all five rows. Even a complete local panel remains a
public synthetic degree-23 result. Project-authored verification confirms
custody and derivable relationships but does not independently replay the
finite-field arithmetic, relation matrix, linear algebra, or rho walk. The
control plane therefore retains `scientific_measurement_admitted: false` until
an independent payload replay exists.

Measured single-core elapsed time remains unavailable; the process meter's
legacy `single_core_seconds` field is aggregate user plus system CPU. Peak RSS
is a maximum process high-water mark, not summed or simultaneous aggregate
parallel memory. The inclusive outer receipt must not be added to child sums.

This stage is not a cryptographic-size attack, private-key recovery,
asymptotic crossover, external reproduction, novelty result, or Koblitz
index-calculus SOTA claim. The frozen machine-readable protocol is
`stage-23-unknown-scalar-protocol.json`.

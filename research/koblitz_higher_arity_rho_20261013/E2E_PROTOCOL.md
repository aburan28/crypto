# End-to-end timing addendum

**Status:** frozen before execution, 2026-10-13.

## Requested scope

| User requirement | Measurement evidence | Acceptance state before execution |
|---|---|---|
| Run the current ECDLP candidate end to end | One native `ecbench` session on the registered n=41 Koblitz public target, from cold factor-base setup through independently verified scalar recovery | pending |
| Measure elapsed time | Five measured rounds after one warmup, with B2, B3, matched strong rho, and a duplicate-rho A/A control interleaved in one session | pending |
| Explain the bottleneck | Exclusive phase clocks, online-window clocks, native counters, process wall time, and peak memory for every measured run | pending |

## Boundary, reference, and question

The boundary is the registered n=41 curve's generic signed-Frobenius floor
`sqrt(pi/(2*82))` in `S`. The reference arm is
`rho.signed_frobenius_strong` on the exact public point used by the IC arms.

The question is narrowly empirical: on this host and this one registered n=41
workload, how much cold end-to-end elapsed time is spent in factor-base setup,
oracle setup, relation collection, linear algebra, and verification for B2 and
B3, and how do their complete solve clocks compare with matched strong rho?

No n=41 result is an exponent fit or an ECC2K-130 transfer result. Operation
counts remain the primary cross-host result. Wall ratios are results only if
both paired runs earn L2 and the comparison lies outside the admitted rho A/A
noise interval.

## Frozen inputs

- Source revision at freeze: `5f83882f976499a55a1f6851281131a85f7611db`.
- Curve: Koblitz `a=0`, n=41, registered as
  `icv1-f2m41-tm2308219-7f48b14a`.
- Target: one public target, `target_seed=2026101301`.
- Algorithm seed: `2026101302`.
- Factor base: `compact-orbit-scan:columns=20,raw_x_cap=100000000`.
- Relation cap: 256 trials.
- Arms: matched strong rho, duplicate strong-rho A/A control, B2 quartic
  six-sum extraction, and B3 two-S3-chain six-sum extraction.
- Schedule: alternate ordering, one warmup, five measured rounds, one thread,
  900-second timeout per run, requested isolation L2.

The exact machine-readable input is [`E2E_SPEC.json`](E2E_SPEC.json). Workload,
seeds, algorithms, accounting intervals, and stop conditions will not be
changed after observing results.

## Accounting intervals

- **Cold solve wall:** the harness's complete solver interval, including all
  cold setup, relation work, linear algebra, and internal verification.
- **Fork-to-reap wall:** cold solve plus process launch/teardown and record
  transport overhead.
- **IC exclusive phases:** factor base, oracle setup, relations, linear
  algebra, and verification, whose clocks are taken directly from the sealed
  records.
- **Online window:** target-dependent work only, as defined by `ecbench`; this
  is reported separately and is not substituted for the cold interval.

## Success and stop conditions

The run succeeds if all 24 planned executions complete and verify, phase clocks
sum consistently, and exact replay reproduces every measured deterministic run.
It stops on a wrong answer, timeout, memory/resource refusal, or the harness's
fail-closed isolation refusal.

If the host refuses L2 because required telemetry or quietness is unavailable,
the refusal is retained. The unchanged spec may then be run with the explicit
`--allow-busy` downgrade to answer the user's diagnostic timing question; every
such elapsed-time number is descriptive L0/L1 evidence, not an admitted wall
claim.

## Pre-execution environment note

The required `conductor check` was attempted with the study, benchmark, and
blocker paths reserved. This runtime has no `conductor` executable (exit 127),
so no scope-conflict certificate is available.

# Protocol E-2: batch verification of column logs and relations

Frozen 2026-10-07, before any instrument is built.  **Engineering,
accounting-preserving**: every verification is still performed and the
soundness error is stated; only its cost changes.  Status: **PENDING**.

## Derivation (stated before measuring)

The producers certify every factor-base column by a scalar multiplication,
`[x_o]G = R_o`, and replay every relation under `assert_eq!`.  The autolab
remeasurement in
[`../../docs/ic/BOUNDARY_TARGETS.md`](../../docs/ic/BOUNDARY_TARGETS.md)
found that `solution_validation_ms` was 7.78 of the 14.78 charged
milliseconds per target on the `partition_walk` arm at `n = 37`, and noted
that dropping it on one side only would turn a 1.15× loss into a 1.65× win
by accounting alone.  Verification therefore stays, on both arms, and gets
cheaper.

**Batch verification.**  Draw independent 64-bit multipliers `s_o` and
check

    Σ_o s_o · x_o · G  =  Σ_o s_o · R_o

with one multi-scalar multiplication on each side (Straus or Pippenger).
If any `x_o` is wrong, the error vector `e` is nonzero modulo `r`, and
`Σ s_o e_o ≡ 0 (mod r)` holds with probability at most `2^{−64}` over the
`s_o`.  The relations are checked the same way: a random combination of the
relation rows, evaluated on the base points and the targets, must vanish.
Pippenger's cost for `K` terms with 64-bit scalars is about
`64·K / log₂ K` group additions against `K · 64 · 1.5` for `K` separate
scalar multiplications at this scalar size, so roughly `log₂ K` times
fewer additions: near 9× at `K = 600`.  A deterministic full check stays
available as an audit mode.

## Instrument (Rust, to build in the follow-on PR)

1. Add `verify_columns_batch` and `verify_relations_batch` to
   `src/cryptanalysis/koblitz_relation_solver.rs`, with a 64-bit CSPRNG
   seeded per run and the seed recorded, and a `--verify full|batch` flag
   in the producers.
2. Give the rho arm the same batch check for its own validation, so the
   comparison is matched.
3. Run `icv1-f2m41-tm2308219-7f48b14a`, `icv1-f2m53-tm56619371-dac20a85`
   and `icv1-f2m83-t6151469093347-cdcc5432`, three runs each, both modes,
   and record the charged verification interval and its instruction count.
4. **Soundness test.**  Plant one wrong column log in 1,024 trials, and one
   wrong relation coefficient in 1,024 trials, and record how many the
   batch check accepts.

## Predictions (pass/fail)

- **V1 (cost).**  The charged verification interval falls by at least 5×
  at `n = 41` and `n = 53`, and by at least 8× at `n = 83`, in
  instructions.
- **V2 (soundness).**  Zero of 2,048 planted faults are accepted.
- **V3 (share).**  Verification is below 10% of the charged rank-stage
  cost and below 10% of the charged online cost on every curve.
- **V4 (verdicts).**  Every verdict, scalar and relation is identical
  between the two modes on every run.

## Decision rule (registered)

V2 and V4 passing admits batch mode as the default charged verification,
with the full check kept for independent replay.  V1 and V3 size the gain.
Class **engineering**.  A single accepted fault in V2 fails the protocol.

## Stop condition and inadmissible moves

Bounded: three curves, two modes, 2,048 fault trials.

Inadmissible: removing verification; applying batch mode to one arm only;
reusing the multipliers across runs; reporting the saving in wall time
without instruction counts; counting the audit-mode check as part of the
charged online interval when it was not run there.

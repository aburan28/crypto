# Fixing the reference rho: a strong single-target backend for `koblitz_rho_fixture`

Follow-on to `RESEARCH_SINGLE_TARGET_STRONG_RHO_N53_20260930.md` (PR #1090), which found that
the rho behind the ledger's selected-panel n = 53 "wall crossover" was a strawman. This note
records the fix. It is engineering on the reference, not a measurement campaign; the numbers
it cites were measured in PR #1090.

## What was wrong

44 stage harnesses invoke `examples/koblitz_rho_fixture.rs`, most of the recent ones as
`koblitz_rho_fixture 53 0 signed_frobenius 1 packed <seed> <scalar>`. Both of its backends
(`reference`, on `BinaryPoint`, and `packed`, on packed `u64` words):

- store **every** step in a hash map (no distinguished points);
- canonicalize the signed-Frobenius orbit by an O(n) polynomial-basis scan;
- invert once per step.

At n = 53 on the panel's target, `packed` cost 43× the instructions and about 60× the wall
time (median over 16 walk seeds) of the strong single-target rho, and the direct index-calculus
arm that "beat" it was 44× slower than the strong rho in wall time (PR #1090).

## The fix

1. **Library module `crypto_lib::cryptanalysis::koblitz_strong_rho`** (`StrongRho`): the
   strong walk as a reusable, tested component — distinguished points, signed-Frobenius
   canonical form in a normal basis (least rotation of the x normal coordinates, sign by y,
   multiplier `±λ^k` from a table), library `Gf2` arithmetic with Itoh–Tsujii inversion by
   rotations, `lanes` lockstep walks with one batched inversion per step. Scalar products mod
   `r` use a float-quotient path below 2^50 and `u128` above, so it covers every rung the
   fixture accepts (including n = 59).
2. **`koblitz_rho_fixture … strong`**: a third backend, signed-Frobenius only, emitting the
   same `rho_public_fixture` JSON row the stage verifiers read (`verified`, `published_q`,
   `recovered_fixture_scalar`, `walk_steps`, `restarts`, `setup_ms`, …) plus
   `"reference_grade":"strong"`, `lanes`, `distinguished_bits`, `walks`, `fruitless_cycles`,
   `capped_walks`, `wasted_merges`. Environment: `KIC_RHO_LANES` (default 32),
   `KIC_RHO_DP_BITS` (default 4).
3. **`reference` and `packed` are unchanged**, so archived stage runs still reproduce, but now
   print a one-line note to stderr that they are not a valid `vs_rho` reference. No stage
   script is edited: the frozen protocols stay as run.
4. **Ledger rule** (`docs/ic/BOUNDARY_TARGETS.md`, `boundary_targets.json`
   `measurement_schema.vs_rho.rho_reference_minimum`): a `vs_rho` comparison must use the
   strong reference, or a rho shown to be at least as fast per step on the same target.

## Validation

- **Unit tests** (library, 6): group law against the library curve on every small rung;
  inversion; both modular-multiply paths against `u128` (moduli up to 2^57); the canonical
  form is the same for all `2n` signed Frobenius images and its multiplier satisfies
  `canon(P) = a'·G + b'·Q`; planted targets recovered at 1, 8 and 32 lanes; median steps
  within a generous factor of `√(πr / 2A)` at n = 37/41. Fixture test: seeded and explicit
  scalars recovered at six rungs with the verifier fields present.
- **Walk identity** (`fixed_reference_rho_20261001/check_walk_identity.py`,
  `walk_identity.json`): the fixture's `strong` backend and
  `koblitz_rho_batch_ks_strong` rung 3 (32 lanes, `dp_bits` 4) were written separately and
  walk **the same trajectory, step for step**, in 13 of 13 cases: n = 37, 41, 53 under batch
  seeds 531310, 1, 2, 3, plus the selected-panel target (scalar 476811900269, seed 531310,
  470,279 steps — the walk PR #1090 measured as its strong arm). So the PR #1090 strong-arm
  measurements (0.132 s, 1.09 G instructions on that target) are measurements of this
  backend's walk.
- **n = 59** (accepted by the fixture, not by the batch example): recovers the planted scalar
  for `a = 0` (r ≈ 2^33.2) and `a = 1` (r ≈ 2^44.5); the `a = 0` process spends about 6.5 s
  constructing the curve, the same for every backend, and 4.8 ms walking (`packed`: 131 ms).
- `rustfmt` clean; `cargo clippy --release --lib --example koblitz_rho_fixture -- -D
  warnings` clean apart from an unknown-lint notice from the repository's own `-A` flag on
  this toolchain.

## Not changed, not claimed

- No existing ratio in the ledger is re-measured here; PR #1090's erratum stands as written.
- The rows at n = 41 that compared against `koblitz_rho_fixture` remain unaudited; rerunning
  them against `strong` is the obvious next measurement.
- `koblitz_signed_frobenius_rho_with_progress` in `koblitz_index_calculus.rs` (the walk the
  `ic boundary` operation counts use) is not touched or audited here.
- The strong backend is single-threaded. Equal-core comparisons need a parallel
  distinguished-point rho (shared table), which this does not add.
- Strong is the best rho built here, not the best possible: after the normal-basis
  canonical form, the least-rotation scan is still about half of its instructions.

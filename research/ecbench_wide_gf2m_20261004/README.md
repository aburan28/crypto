# Wide GF(2^m) Koblitz curves in ecbench: the m = 83 gate becomes runnable

**What changed.** `ecbench` now builds the m = 83 confidence-gate curve of
AGENTS.md §8a, `icv1-f2m83-tm6151469093347-debefd74`
(`E_0: y² + xy = x³ + 1` over `GF(2^83)` modulo `z^83 + z^45 + z^2 + z + 1`,
`#E = 4 · 2417851639230796216685689`). It runs `rho.signed_frobenius_strong`
and `claw.pair_table` on it, on the same one-target workloads and into the same
sealed records as every other curve. This is infrastructure. No measurement at
m = 83 is reported here, and the cost of one is priced below.

| | |
|---|---|
| field | `src/cryptanalysis/wide_gf2m.rs`: `GF(2^n)`, `n ≤ 127`, on `u128` words |
| curve | `src/cryptanalysis/koblitz_wide.rs`: `K_a` for `64 ≤ n ≤ 127` |
| strong rho | `src/cryptanalysis/koblitz_strong_rho.rs`, generic over word and scalar width |
| claw | `src/cryptanalysis/ecbench/claw.rs`, generic over `ClawGroup` |
| harness | `ecbench/workload.rs` (`Instance::Wide`, `_v2` target laws), `ecbench/methods.rs` (`solve_wide`), `ecbench/canonical.rs` (`compat_u128`) |
| smoke session | [`sessions/smoke-m83`](sessions/smoke-m83) `ECBS1h487c13e3a766`, spec [`spec-smoke-m83.json`](spec-smoke-m83.json), audit [`audit-smoke-m83.json`](audit-smoke-m83.json) |

## How it was built so the existing evidence cannot move

Every committed session replays exactly in CI, so the work was done in
three steps, each gated on that replay.

1. **Wide scalars without new bytes.** Subgroup orders, the planted scalar
   and the recovered scalar became `u128` under `canonical::compat_u128`.
   It writes a JSON number while the value fits a `u64` and a decimal string
   above that. One-word curves keep `uniform_scalar_sha256_v1` and
   `hash_to_subgroup_v1` exactly. A wide order gets
   `uniform_scalar_sha256_v2` and `hash_to_subgroup_v2`, which hash the
   same inputs as decimal strings with a 128-bit draw.
   **Gate:** all 26 committed sessions re-audit `ok`, with record hashes,
   workload ids and run ids unchanged.
2. **One strong rho, two widths.** The walk is now generic over the field
   word (`u64`/`u128`) and the scalar (`u64`/`u128`). `StrongRho` and its
   companion names are the one-word instantiation, so the committed code is
   that instantiation. The hash mixers fold a wide word to 64 bits by a map
   that is the identity below `2^64`.
   **Gates:**
   - 6,449 committed replays are identical. Every session is replayed in
     full except the pair-claw and yield-sweep sessions, which are sampled
     with 60 replays each.
   - A new unit test runs the `u128`-word walk on the one-word curves from
     n = 17 to 41. It equals the `u64`-word walk step for step and counter
     for counter.
3. **One claw, two widths.** The pair claw is written once over
   `ClawGroup`. `NarrowClaw` is the committed one-word code. `WideClaw` uses
   the strong reference's `u128` arithmetic and normal basis.
   **Gates:**
   - 160 of 160 sampled replays of the committed pair-claw session are
     identical.
   - The wide claw recovers planted logs at both shapes on
     `icv1-f2m41-tm2308219-7f48b14a` and `icv1-f2m47-t22705043-f4e44623`.

Field and curve checks:

- `WideGf2` agrees with `F2mElement` at n = 67, 83 and 127, on both the
  hardware carry-less kernel and the portable one.
- `WideKoblitz` reproduces §8a's order and the registered slug at m = 83.
- Its group order agrees with the one-word constructor wherever both
  exist.
- The `u128` scalar ring agrees with `BigUint` up to a modulus of
  `2^127 − 1`.

## The smoke session

Two public targets on the m = 83 curve, both arms with budgets far below a
solve: strong rho capped at 10^6 steps (`step_cap_factor 0`), and the claw
at `cap_multiple 0`. Both arms exhausted, as designed. The run exercises
workload build, the `_v2` public-target law, wide dispatch, sealed records
and the audit. It is not a measurement. Its records carry `r` as the string
`"2417851639230796216685689"`.

## What an m = 83 measurement costs

The figures below are rates measured on this host (Apple silicon, one
core, isolation L0), then extrapolated.

- **Strong rho.** 0.49 µs per step at 32 lanes, from the smoke run's
  10^6-step cap. A solve needs `√(πr/4n) ≈ 2^37.1 ≈ 1.5·10^11` steps,
  about **20 single-core hours per target**. A panel of 8 targets is about
  160 core-hours, run in parallel across cores.
- **The claw.** Its cost-minimising table at m = 83 has `√(r/n) ≈ 2^37.3`
  classes, far beyond one host's memory. The PR's own shape
  (`table_scale ≈ 1/32`) still holds about `2^32` classes. At the measured
  ratio of 32.6× strong rho, that shape needs about `2^42` additions, or
  roughly a month of one core per target.
  [aburan28/cryptanalysis#175](https://github.com/aburan28/cryptanalysis/pull/175)
  reached its n = 83 hit with a sharded Bloom-filter table on a CI fleet,
  which this single-process port does not have.

These are extrapolations from the measured per-step rate and from the flat
claw-to-rho ratio of
[`ecbench_pair_claw_20261003`](../ecbench_pair_claw_20261003/README.md)
(slope −0.015 over r = 2^32 to 2^47).

## Limits

- Of the registered wide Koblitz curves, only m = 83 and m = 97
  (`icv1-f2m97-t378251973071909-f27192f3`) have a prime `#E/h`. m = 67, 71,
  73, 79 and 127 are refused: for m = 79, `#E/4` has the factor 149627. So
  the cheapest wide solve is the gate curve itself.
- The generator is the constructor's (the first `x = 1, 2, …` on the curve,
  cofactor cleared). The registry has no EC1 representation for it yet, so
  records carry `ec1: null`. A cross-repository join with the PR's
  `EC1N83Ckb1h876c2921cb64` needs that representation registered.
- No other method has a wide path. `ic.pipeline`, BSGS and the kangaroo
  refuse a wide curve.
- While validating, the committed `ecbench_yield_sweep_20261004/sessions/koblitz`
  session failed to replay identically on this host with **main's own
  binary**: 6, 13 and 16 of 60 replays differed in three runs, all in
  `descent-algebraic` IC arms. This is pre-existing, unrelated to this
  change, and was raised as a separate task.

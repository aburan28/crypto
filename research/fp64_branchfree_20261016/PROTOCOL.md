# Protocol: branch-free selections in production `fp_backend::Fp64` and `prime_fast::Fp64`

Frozen before any run below. The protocol commit is the last commit that
changes this file or the code it tests before the first run; it is recorded
in `runs/protocol_commit.txt` when the runs start, and both binaries are
built from it. A first version covering
`fp_backend::Fp64` alone was committed at `7332be2e1`; it was superseded
before any run, when `prime_fast::Fp64` was added to the scope, and no
measurement was taken under it.

This protocol follows `hwcost` round 2 (`research/hwcost_round2_20261016/`,
PR #2003). Round 2's addendum 1 rewrote `Fp64`'s data-dependent selections
with `std::hint::select_unpredictable` (`hwcost::kernel::Fp64Cm`). That cut
the rho kernel to 0.40–0.50× of `Fp64`'s time at every B ≥ 4096 on this
host's baseline-ISA build, with results bit for bit the same.

## The change

Two production Montgomery `u64` fields, each with `R = 2⁶⁴`:

- `fp_backend::Fp64` (odd `p < 2⁶³`), the `PrimeField` backend used by the F4
  engine, `two_sylow_search`, the `Fp52` IFMA backend's scalar operations
  and `hwcost`;
- `prime_fast::Fp64` (odd `p < 2⁶²`), the arithmetic of the production
  prime-field rho engine (`rho_parallel`) and of `prime_fast`'s
  index-calculus relation search.

Three selections change in each:

- the final conditional subtraction of the Montgomery reduction (`redc` /
  `reduce`);
- the final conditional subtraction in `add`;
- the wrap-around in `sub`.

Each becomes one `overflowing_sub` and one `select_unpredictable`, exactly as
in `Fp64Cm`. `neg` keeps its branch: `a == 0` is almost never true on data,
so the branch predicts. Every user of either type inherits the change.

**Correctness gate (before any timing).**

- **New tests.** Two unit tests, both named
  `fp64_selects_match_branchy_formulas` (one in `fp_backend`, one in
  `prime_fast`), compare reduction, `mul`, `sqr`, `add` and `sub` (and
  `fp_backend`'s `double`) with the previous branchy formulas, bit for bit.
  They cover:
  - every pair of an edge set (0, 1, 2, p − 2, p − 1, ⌊p/2⌋, ⌊p/2⌋ + 1);
  - 20,000 random pairs;
  - the reduction on 20,000 random inputs across its whole range
    `t < p·2⁶⁴`.

  The primes are the existing `fp_backend` lists, and 5, 7, 13, 10⁹ + 7,
  2³¹ − 1, 2⁶¹ − 1 and 2⁶² − 57 for `prime_fast`.
- **Existing tests.** These must pass unchanged: `fp64_matches_biguint`,
  the `f4_tower_generic` `Fp64` tests, the `Fp52` tests, `two_sylow_search`,
  and `prime_fast`'s known-answer rho and index-calculus tests.
- **Run gates.** Every `hwcost measure` run checks each field's walk against
  an unbatched reference, and the timed arithmetic against `fp_backend::Fp64`
  at the same prime (for `pf64`, the representations coincide). Every
  `pfrho` solve asserts its logarithm.

The results are bit for bit the previous ones, so no recorded output of any
other study changes.

## Arms and runs

- **Binaries.**
  - *Candidate (B):* `cargo build --release --bin hwcost` at the protocol
    commit.
  - *Baseline (A):* the same tree with `src/cryptanalysis/fp_backend.rs` and
    `src/cryptanalysis/prime_fast.rs` restored to their content at `main`
    `f25d5e520`, so the measuring code is identical and only the two field
    implementations differ.

  Both use the default target (x86-64 baseline), the repository's production
  build. Each run records its binary's SHA-256.
- **Host.** This session's cloud container: Intel Xeon @ 2.10 GHz, family 6
  model 207, 4 vCPUs, KVM. Isolation is L0. There is no PMU, so there are no
  branch-miss counts.
- **Runs.** Four runs back to back: A1 baseline, B1 candidate, A2 baseline,
  B2 candidate. Each is
  `hwcost measure --fields fp64,fp64cm,pf64,w1 --prime-fast-rho-bits 44`
  with round-1 parameters (7 repetitions, 100 ms, batches 2⁴ … 2²², all arm
  groups), writing `runs/base-1`, `runs/cand-1`, `runs/base-2` and
  `runs/cand-2`. Each run is fitted with `hwcost fit`.
- **Fields.**
  - `fp64`: production `fp_backend::Fp64` at 2⁶³ − 25, under test.
  - `pf64`: production `prime_fast::Fp64` at 2⁶² − 57, through
    `hwcost::kernel::PrimeFast`, under test.
  - `fp64cm`: `hwcost::kernel::Fp64Cm`, whose selections were already
    branch-free and whose code does not change. It is a control and the
    reference for P3.
  - `w1`: `FpWide<1>`, untouched, a control.
- **`pfrho.44`.** The production rho `prime_fast::rho_parallel` end to end,
  one known-answer solve per sample on `find_a3_curve(44, 3)`. It uses one
  thread, 1024 walkers, the negation map, a 2¹⁰ table and 6 DP bits, and
  reports wall ns per walk step with table set-up included.

`hwcost compare` computes, per arm, `r_k = B_k / A_k` and the drift
`d = |A2/A1 − 1|`. The claims are computed by
`hwcost compare --max-ratio/--max-dev` and `hwcost pair` (Rust), not by the
report.

| id | role | claim | rule |
| --- | --- | --- | --- |
| P1 | gate | the rho kernel on `fp_backend::Fp64` gets faster | `kernel.fp64` at every B ≥ 4096 (6 sizes): max(r1, r2) ≤ 0.85 |
| P1pf | gate | the rho kernel on `prime_fast::Fp64` gets faster | `kernel.pf64` at every B ≥ 4096: max(r1, r2) ≤ 0.85 |
| P4 | gate | the production rho gets faster | `pfrho.44`: max(r1, r2) ≤ 0.95 |
| P2a, P2b | guard | dependent multiply and squaring chains on `fp_backend::Fp64` do not get slower | `fp64.mul.latency`, `fp64.sqr.latency`: max(r1, r2) ≤ 1.10 each |
| P2c, P2d | guard | the same on `prime_fast::Fp64` | `pf64.mul.latency`, `pf64.sqr.latency`: max(r1, r2) ≤ 1.10 each |
| P3 | match | production `fp_backend::Fp64` matches the study's `Fp64Cm` | in cand-1 and cand-2, `fp64 / fp64cm` kernel ratio ≤ 1.10 at every B ≥ 4096 |
| C1a | control | unchanged `Fp64Cm` kernel does not move | `kernel.fp64cm`, B ≥ 4096: max \|r_k − 1\| ≤ 0.10 |
| C1b | control | unchanged `FpWide<1>` kernel does not move | `kernel.w1`, B ≥ 4096: max \|r_k − 1\| ≤ 0.10 |
| D1 | descriptive | every arm's class (faster / slower / no detected effect beyond max(2d, 2%)) | counts per group, no threshold |
| D2 | descriptive | the machine code changed as intended | instructions, conditional jumps and `cmov` of `step_fp64`, `step_fp64cm` and `step_pf64` in both binaries (`objdump -d`) |

P2 guards the latency-bound users (exponentiation chains in `pow`, `sqrt`,
Legendre and the batch inversion's running product). On such a chain a
`cmov` adds latency where a well-predicted branch does not.

P4's 5% bar is smaller than P1's because `rho_parallel` spends its time on
more than field selections: the table lookup, the negation map's
data-dependent sign choice (unchanged) and distinguished-point handling.

Each type's change lands only if its gates and guards pass:
`fp_backend::Fp64` on P1, P2a and P2b; `prime_fast::Fp64` on P1pf, P4, P2c
and P2d. P3 and C1 qualify the reading: a control outside ±10% means the
builds differ in more than the change, and the gates are reported as
confounded. Any failure is reported as measured.

## Outputs

- `runs/{base,cand}-{1,2}/`: `measurements.json`, `samples.csv`,
  `medians.csv`, `validation.{json,csv}`, `profile.json`, `fit.txt`,
  `stderr.txt`, `binary.txt` and timestamps.
- `runs/compare/`, `runs/cm-pair/`, `runs/disasm/`.
- A run report and a results report, with graphs and PDFs.
- `benchmarks/20261016-fp64-branchfree/`.

## Commands

```sh
S=<scratch dir>; D=research/fp64_branchfree_20261016/runs; F=src/cryptanalysis
cargo build --release --bin hwcost && cp target/release/hwcost $S/hwcost-cand
git show f25d5e520:$F/fp_backend.rs > $F/fp_backend.rs; git show f25d5e520:$F/prime_fast.rs > $F/prime_fast.rs
cargo build --release --bin hwcost && cp target/release/hwcost $S/hwcost-base
git checkout HEAD -- $F/fp_backend.rs $F/prime_fast.rs
for r in base-1 cand-1 base-2 cand-2; do
  mkdir -p $D/$r; $S/hwcost-${r%-*} measure --fields fp64,fp64cm,pf64,w1 --prime-fast-rho-bits 44 \
    --out $D/$r 2> $D/$r/stderr.txt
done
for r in base-1 cand-1 base-2 cand-2; do $S/hwcost-cand fit --run $D/$r > $D/$r/fit.txt; done
mkdir -p $D/compare $D/cm-pair
$S/hwcost-cand compare --base $D/base-1,$D/base-2 --cand $D/cand-1,$D/cand-2 --out $D/compare \
  --max-ratio P1:kernel.fp64.:4096:0.85 --max-ratio P1pf:kernel.pf64.:4096:0.85 \
  --max-ratio P4:pfrho.44:0:0.95 \
  --max-ratio P2a:fp64.mul.latency:0:1.10 --max-ratio P2b:fp64.sqr.latency:0:1.10 \
  --max-ratio P2c:pf64.mul.latency:0:1.10 --max-ratio P2d:pf64.sqr.latency:0:1.10 \
  --max-dev C1a:kernel.fp64cm.:4096:0.10 --max-dev C1b:kernel.w1.:4096:0.10 > $D/compare/summary.txt
$S/hwcost-cand pair --runs $D/cand-1,$D/cand-2 --a fp64cm --b fp64 --min-batch 4096 \
  --threshold 1.10 --id P3 --out $D/cm-pair > $D/cm-pair/summary.txt
for b in base cand; do research/fp64_branchfree_20261016/disasm.sh $S/hwcost-$b $b $D/disasm; done
```

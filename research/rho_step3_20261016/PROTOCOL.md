# Strong rho, round 3 (the walk step): protocol

Declared on 2026-10-16, before any change in this round was written or measured.
Parent: `ed282c5b9` (round 2, #2019).

## Where the time is

A `perf record -e cpu-clock` sample of a line-table build of the parent (n61,
8 fixtures; `runs/profile/`) shows:
- `Gf2::mul_fold` at 18.5% of the walk's time: every field multiplication is a
  call. `mul_fold` is a `pclmulqdq` target-feature function, so it does not
  inline into the generic walk code.
- `Gf2::inv_clmul`, the batch inversion's single inverse, at 8.4%.
- `StrongRhoG::inv`, the walk starts' inverses, at 4.2%.

## Candidate changes (no value changes)

1. **Run the walk where the field's arithmetic inlines.**
   - `RhoField` gains `with_fast_arith(f)`. For `Gf2`, when the field folds
     (which implies `pclmulqdq`), it runs `f` inside a
     `#[target_feature(enable = "pclmulqdq")]` function; otherwise, and for
     `WideGf2`, it calls `f` as before.
   - `solve`'s body moves to an `#[inline(always)]` `solve_impl` run through
     it, as `koblitz_fast` and `Gf2::batch_inv_clmul` already do.
   - The hot helpers (`add_with_slope`, the multiply and square wrappers) are
     `#[inline(always)]`, so the multiplications inline into the step.
2. **Walk-start batching**, measured separately and only if change 1 lands.
   The start points `cursor + j·stride` for the next B starts are formed
   together, with one batched inversion instead of one inversion each. The
   same points come out in the same order, and the same group-addition count
   is charged.

## Identity gate

As in round 2:
- every output line of `koblitz_rho_fixture n 0 signed_frobenius 4 strong`
  at n = 41, 53 and 61 is identical, timing aside;
- all strong-rho unit tests pass with debug assertions on.

## Measurements

As in round 2:
- **Primary: Callgrind instructions per walk step** at n53, 4 fixtures.
- **Native time:** A/A (baseline × 6), then 5 interleaved A/B rounds under
  `isolated_bench`, at n = 53 and 61.
- **Call it established only outside the A/A spread.**

Target-feature inlining mainly removes call overhead and lets the compiler
schedule independent products together. Its effect on time can therefore
differ from its effect on instruction counts, and both are reported.

## Addendum 1 (2026-10-16, after change 1 was built and counted, before any timing)

- **Change 1 as built differs in form, not in effect.** The declared closure
  trampoline (`with_fast_arith(f)`) did not inline: the closure compiled as
  its own function without the target feature, `mul_fold` was still called
  per product, and the count rose (2,474,427,299 instructions at n53, 4
  fixtures). It was replaced by `solve_clmul`, a `#[target_feature(enable =
  "pclmulqdq")]` method that `solve` calls when `RhoField::clmul_walk()` holds
  (a folding `Gf2` on x86-64), around an `#[inline(always)]` `solve_impl`.
  Outputs are identical; the count is 2,323,045,497 (−2.0%). Change 1 lands.
- **Change 2 is measured next, as declared.** `starts_ahead` forms
  `base + (j+1)·stride` for `j < lanes` with one batched inversion and falls
  back to the group law for infinity or equal x. The `lanes` multiples of the
  stride are formed once per solve, outside the charged counts, as the stride
  itself is. A new unit test checks the batch against the sequential chain,
  including infinity, doubling and an infinite base.
- **Timing compares the round-2 build (`ed282c5b9`) with changes 1+2**,
  under the declared A/A and A/B schedule.

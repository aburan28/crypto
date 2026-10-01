# R02: the AVX-512 batched addition for fields with a wide reduction tail

**Declared 2026-10-01, before any candidate code and before R01's
results were read.** This comes from ledger §23's reports and a reading
of the code; R01's profile checks it rather than sources it. The plan
is `research/notes/index-calculus/IC_TOOL_PROGRAM.md` (Track A, A2).
Nothing below changes after the first run except by a dated amendment
appended at the end.

## Hypothesis

**The scan is slower at the two largest sizes because the 8-lane kernel
is off there.** At `K_1/GF(2^59)` and `K_0/GF(2^61)` the collection
scan's batched addition falls back to the scalar path. The SIMD kernel
`Simd512` in `src/cryptanalysis/koblitz_fast.rs` refuses those fields:
- It reduces a product `L + z^n·H` through the tail `t` of the field's
  polynomial (`z^n ≡ t`).
- It needs `H·t` to fit in one word, which means `n + deg t ≤ 65`.
- Both fields have `n + deg t = 66`:
  - `GF(2^59)` is `z^59 + z^7 + z^4 + z^2 + 1`;
  - `GF(2^61)` is `z^61 + z^5 + z^2 + z + 1`.
- Every other suite size qualifies, at 44–61.

**The cost per summand follows that split exactly.** In §23's reports,
at `0bf67f16`, one thread, isolated, a scanned summand costs:

| curve | `n + deg t` | kernel | ns per summand | units per summand | collection's share of the set-up |
|:--|--:|:--|--:|--:|--:|
| `K_1/GF(2^47)` | 52 | on | 38.7 | 1.44 | 42% |
| `K_0/GF(2^57)` | 61 | on | 44.8 | 1.48 | 31% |
| `K_0/GF(2^41)` | 44 | on | 33.4 | 1.37 | 67% |
| `K_0/GF(2^53)` | 59 | on | 48.7 | 1.59 | 65% |
| `K_1/GF(2^59)` | 66 | **off** | 71.7 | 2.38 | 71% |
| `K_0/GF(2^61)` | 66 | **off** | 68.9 | 2.30 | 89% |

**The candidate.** The kernel carries the bits of `H·t` that pass `z^63`
into its second fold. Those bits are `h >> (64 − t_i)` for each tail
term `t_i ≥ 1`. The second fold's input becomes
`(f >> n) | (o << (64 − n))`. That input still has degree at most
`deg t − 2`, and the fold's output lands below `z^n` while
`2·deg t ≤ n + 1`, which is the kernel's existing second condition.

So the field elements are the same as the scalar path's, lane for lane.
The condition `n + deg t ≤ 65` is dropped.

- **Fields that already qualify keep their exact instruction stream.**
  The carry is a const generic, monomorphised away when off.
- **The kernel serves four callers:** the scan, the folded pair-table
  build, the full-scan variants and the descent's walk. All four take it
  at the two sizes.

## Class

Engineering, if it gains: the counts are unchanged by construction, so
the ratio to the floor is flat.

## A control first, on v0, before any candidate is timed

v0 at `K_0/GF(2^53)`, `M1`'s two rows, five rounds ABAB. One arm sets
`KIC_SCAN_SIMD=0`, the existing same-binary switch to the scalar path;
the other is the default. That is 20 isolated processes.

- **Confirmed** if the scalar arm's collection time per scanned summand
  is at least 1.3× the kernel arm's (geometric mean over the 10 pairs).
- **Otherwise the mechanism is wrong.** R02 stops before the candidate
  is timed and records the control.

## Pinned outputs

- **Every row, untimed.** The candidate runs on all 88 S rows and both
  smoke rows. Every row's outputs must equal v0's from R01's profile
  pass: counts, both arms' scalars, rho's counts, and the verification
  flags.
- **Tests:**
  - the kernel against the scalar field on every degree 2–62 the
    relaxed condition accepts;
  - the batched addition against the scalar batched addition on the
    `K_1/GF(2^59)` and `K_0/GF(2^61)` curves themselves;
  - the two suite fields marked wide.

## Timed comparison

- **Rows.** v0 against the candidate on all 88 S rows, five rounds,
  with the arms' order alternating (IC_TOOL_PROGRAM.md §5). That is 880
  isolated processes.
- **Holdouts.** Two holdout rows per size, on seed 205 with targets
  `T101` and `T102`: `public_hash_seed` 23101 and 23102, rho seeds
  `0x230000 + 101` and `+ 102`. Five rounds give another 220 processes.
- **The A/A.** R02 runs in R01's container session. It uses R01's A/A
  if the host manifest's CPU, kernel, memory and THP settings are
  unchanged, and runs its own on `M1`'s rows otherwise.
- **Callgrind, untimed.** Both arms run at `K_0/GF(2^61)` and
  `K_0/GF(2^41)` `M1-T01`, `--repeats 1`, as in R01. The figures are
  the instruction counts and the collection's share.
- **If R02 is accepted,** the new baseline gets the rule's comparison:
  §23's protocol at its six sizes with 64 targets, on the candidate.

## Prediction

- **`K_0/GF(2^61)`: 1.37×** [1.18, 1.53], the paired cold-time ratio,
  v0 over the candidate. Collection, 89% of the set-up, falls from 2.30
  units a summand to `K_0/GF(2^53)`'s 1.59. The range runs from 1.9 to
  1.4 units. The build's and the descent's gains are left out, so the
  prediction is conservative.
- **`K_1/GF(2^59)`: 1.31×** [1.15, 1.45], by the same arithmetic at
  71%.
- **The other nine sizes: 1.00×.** Their kernels are unchanged.

## Success and stop

**Accepted** if all of the following hold:
1. The control confirms the mechanism.
2. Every pinned output is identical.
3. At `K_1/GF(2^59)` and `K_0/GF(2^61)`, the paired cold-time ratio's
   95% interval lies above 1.10, on the suite rows and on the holdouts
   separately.
4. No size regresses beyond its A/A band. A regression is a candidate
   interval whose upper end lies below the lower end of that size's A/A
   interval.

**Rejected** otherwise. The code is kept on record and the numbers go
in the ledger.

**Stopped early** if the control fails, if any output differs, or if a
verification fails.

## Inadmissible

- Changing the field polynomial: that changes every coordinate, so it
  is a different instance.
- Changing the suite, the recipes or the targets.
- Pooling contended runs.
- Quoting the stage figures (collection, build) as the speedup: the
  speedup is cold time.
- Claiming a gain at a size whose kernel did not change.

## Cost

- The control: 20 processes.
- The pin: 90 untimed.
- The comparison: 1,100 processes, about 2.5 hours.
- Callgrind: four runs, about 40 minutes.
- If accepted, the rule's comparison: 384 processes, about 40 minutes.

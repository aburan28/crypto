# lit2026 W12: flagged ECDLP claims — zero-minor guess-and-determine, and Hu ePrint 2025/015

**Workstream.** W12 of the [lit2026 speedup program](../lit2026_speedup_program/README.md)
(kind `audit`).
**Date.** 2026-10-11 (run started 2026-10-12 00:09 UTC).
**Kind of result.** Diagnostic, with a negative outcome. No speedup is claimed. No
deployed-curve or asymptotic claim is made; the repo's claim-hygiene rules apply.
**PDF.** [`REPORT.pdf`](REPORT.pdf). Run report:
[`reports/run-20261012T0009Z/report.pdf`](reports/run-20261012T0009Z/report.pdf).

## Questions

1. **Zero-minor method.** Mahalanobis, arXiv 2607.09814v2, "A Guess and Determine
   Attack on the ECDLP", builds on arXiv 1703.07544 and arXiv 2310.04132 ("Minors solve
   the ECDLP").
   - Does a faithful implementation solve planted prime-field instances?
   - How does its counted cost scale with the group order `r`, against the
     repository's `rho.negation` measured in the same unit?
   - Falsifier for any speedup claim: a log-cost slope ≥ 0.5, or a constant far above
     rho.
2. **Hu, ePrint 2025/015.** Does it specify an algorithm that can be implemented and
   measured? See [`HU_EPRINT_2025_015.md`](HU_EPRINT_2025_015.md).

## Answer in one paragraph

The reduction and the guess-and-determine search are implemented in Rust:
[`examples/lit2026_w12_zero_minor.rs`](../../examples/lit2026_w12_zero_minor.rs).

- **Correctness.** They solved all 48 planted targets on the six registered
  prime-order curves from 16 to 26 bits, and every answer was verified outside the
  solver.
- **Scaling.** Measured in counted `F_p` operations:
  - the search phases scale with slope 0.905 (95% bootstrap 0.769–1.021) in `log₂ r`;
  - the whole method, including the fixed set-up, scales with slope 0.783
    (0.669–0.879);
  - the number of candidate minors examined before a zero minor scales with slope
    0.929 (0.793–1.048), and its median is 0.37–2.1 × `r` per curve;
  - `rho.negation` on the same targets scales with slope 0.377 (0.350–0.404) in field
    ops. That is below 0.5 because the fixed set-up of the walk weighs on 16–26-bit
    curves.
- **Paired cost ratio.** Zero-minor over rho grows from 260× (215–319) at 16 bits to
  6715× (2263–15336) at 26 bits.

Both falsifier conditions hold. The slope is above 0.5, and the constant is 2–4 orders
of magnitude above rho. The data matches the heuristic that each candidate maximal
minor is zero with probability about `1/r`, so the search is exhaustive in `r`. The
paper itself concludes the same: its Table 2 says the attack "is no better than
exhaustive search".

The paper also contains an underspecified "Pollard rho on penultimate intersections".
It was not implemented; see Ambiguities.

**Hu 2025/015** is blocked:

- The full text sits behind an IACR Cloudflare challenge, which was not bypassed.
- The abstract leaves "IICA" undefined.
- The abstract concerns `F_p^*` discrete logarithms, not the ECDLP.

## Papers (read on arXiv, full text)

| id | what it contributes | used here |
|---|---|---|
| arXiv 1703.07544 (Mahalanobis–Mallick, INDOCRYPT 2018) | Algorithm 1: the matrix `M` of degree-`n'` monomials at `3n' + l` points; its left kernel `K` has dimension `l` (Thm 2.3); a vector of `K` with `l` zeros means its support sums to `O` (Cor 2.4) | point set, `M`, `K`, decoding |
| arXiv 2310.04132 (Abdullah–Mahalanobis) | anti-diagonal form `K = [A \| J]`; a zero minor of `A` gives a maximal zero minor of `K` (Thm 3.1); search heuristics (all 2-minors, GESC, APM) | anti-diagonal form, zero-entry check, Theorem 3.1 correspondence (tested) |
| arXiv 2607.09814v2 (Mahalanobis, 2026-08-06) | guess `a` (ℓ' dense columns) and recursion to `K' = K[a]ᵀ` in anti-diagonal form (Thm 1, Alg 1 lines 3–8); determine step: `T` = left kernel of `K'[b]`, signature matrix, hashed penultimate intersections, collision gives a zero minor (Thms 3–5, Alg 1 lines 10–22, Alg 2); defect `d` from Eq. (4) | implemented as specified |

## Method (what the code does)

Per target, with the paper's parameter rules:

- **Kernel.**
  - Set `n' = ` the largest even integer `≤ log₂ r`, so that `ℓ = 3n'` is even
    (2607.09814 §1–2).
  - Draw `ℓ` points `n_i P` and `ℓ` points `−n_j Q`, all distinct and random.
  - Build `M` (`2ℓ × C(n'+2, 2)`) and its left kernel `K` (`ℓ × 2ℓ`; `dim K = ℓ` is
    asserted).
  - Reduce `K` to anti-diagonal form by an RREF that visits the sparse columns first,
    in reverse order.
  - Check every dense entry for zero. A zero entry is a 1×1 zero minor and solves the
    instance.
- **Guess.**
  - Draw a random `a ⊂ [ℓ]` with `|a| = ℓ' = ℓ/2`.
  - Form `K' = K[a]ᵀ` and reduce it to anti-diagonal form. A singular block is itself
    a zero minor.
- **Determine.**
  - Enumerate every `b ⊂ [ℓ']` with `|b| = ℓ' − d` in lexicographic order.
  - Compute `T` = left kernel of `K'[b]` (dimension `d`).
  - Build the signature matrix, `(ℓ' + d) × d`:
    - `T·c` for the `d` dense columns `c ∉ b`;
    - the columns of `T` for the `ℓ'` sparse columns, which are free.
  - For every `(d−1)`-subset of its rows, compute the 1-dimensional right kernel and
    hash it. A repeat whose union has `d` rows is a zero `d × d` minor (Thm 5,
    Alg 2). It gives a zero maximal minor of `K'`, hence of `K` (Thm 1).
- **Decode.**
  - Find the vector of the row space of `K` that vanishes on the zero minor's columns.
  - Its support is a set of points summing to `O`, so
    `m = Σ n_i / Σ n_j mod r`.
  - The runner checks `[m]G = Q` in a separate step.
- **Defect.** `d` is the least value with Eq. (4) success ≥ 0.9; this is the paper's
  rule (Tables 2–3). It gives `d = 3` at 16–24 bits and `d = 4` at 26 bits.
- **Counting.**
  - Every `F_p` multiplication and inversion goes through a counting `Fp`. Additions
    are free, the convention of `ic_boundary::FieldOps`.
  - The `2ℓ` scalar multiplications use `PrimeCurve`'s counting group law.
  - Hash probes are recorded as an uncharged counter.
- **Reference.** `ic_boundary::rho_walk_with(…, NegationClasses, RhoWalk::negation())`
  is the code path of ecbench's `rho.negation`. It ran on the same targets with 4
  walk seeds each, and is charged the field ops its group law counts.

### Why a standalone example, not an ecbench method

- ecbench's unit is the group operation, with field ops as a secondary tally of the
  group law. This method's work is linear algebra over `F_p` and has no group
  operations after set-up. As an ecbench method it would be charged only its `2ℓ`
  scalar multiplications, with everything else in `_uncharged`, and that would
  misstate its cost.
- The comparison is therefore made in the one unit both methods share: counted field
  operations.
- Adding a method would also change `methods.rs` and the IC candidate identity, and it
  would need a sealed session, for a method whose verdict is clear without one.
- This is a scope decision, recorded here. It is not a blocker.

## Ambiguities in arXiv 2607.09814, and how each was resolved

1. **Signature normalisation.** The paper takes the last row of an RREF of
   `[Hᵀ | Tᵀ]` (Thm 3). That is `T·c` up to a scalar, and only the line matters for the
   kernels. The code uses `T·c`. A unit test checks that `T` is the left kernel.
2. **Representative of a penultimate intersection.** The paper says "first non-zero
   entry 1". The code uses the RREF free-variable vector, whose last non-zero entry is
   1. Both are canonical representatives of the same line, so the set of collisions is
   the same.
3. **Collisions whose union exceeds `d` rows.** Algorithm 2 returns "No result" for
   these. The code follows that, counts them (`collision_union_gt_d`, 0 in every run),
   and keeps enumerating.
4. **`(d−1)`-row submatrices of rank `< d−1`** ("general position" is assumed and not
   checked, §4). The code treats them as zero minors, which they are, and counts them
   (`sig_rank_deficient`, 0 in every run). The same holds for `dim T ≠ d`
   (`t_dim_not_d`, 0).
5. **What happens after all `b` fail for one `a`.** The paper's Las Vegas loop is
   implicit. The code draws a new random `a`, and new points on a singular block.
   Every run in this session succeeded on its first `a` and first kernel.
6. **`n'`.** The paper says "approximately `log₂ p`", and `ℓ` must be even. The code
   uses the largest even integer `≤ log₂ r`; `--nprime` overrides it.
7. **Not implemented: the "Pollard rho on penultimate intersections" (§5, Table 1).**
   The paper defines no walk, iteration function or collision rule, and says "Further
   research is needed for a proper implementation". No faithful implementation is
   possible.
8. **Gray-code and rank-one kernel updates** (suggested in §2 and §5) are not
   implemented. The phase tallies bound what they could save:
   - at 26 bits, `determine_t` is 7.7% of the ops and `determine_kernels` 91%;
   - at 24 bits they are 68% and 24%.

   Removing the `T` computation entirely would still leave the `(d−1)`-kernel count,
   which scales with the covered minors (slope 0.93).

## Results (run-20261012T0009Z, L0, counted field ops; full table in the run report)

| bits | curve | r | n', ℓ', d | rho.negation median ops (32 runs) | zero-minor median ops (8 runs) | paired ratio zm/rho, geomean (95% bootstrap) |
|---|---|---|---|---|---|---|
| 16 | icv1-fp16-t295-8d3c3165 | 53117 | 14, 21, 3 | 1.92e3 | 5.53e5 | 260 (215–319) |
| 18 | icv1-fp18-tm337-d28d5e09 | 212437 | 16, 24, 3 | 2.98e3 | 1.48e6 | 455 (326–626) |
| 20 | icv1-fp20-t727-cd198a38 | 562427 | 18, 27, 3 | 4.62e3 | 2.54e6 | 563 (350–883) |
| 22 | icv1-fp22-tm1385-475dcb5f | 2435809 | 20, 30, 3 | 7.93e3 | 9.03e6 | 785 (349–1633) |
| 24 | icv1-fp24-t1059-6df599df | 16482251 | 22, 33, 3 | 1.67e4 | 2.16e7 | 1109 (706–1743) |
| 26 | icv1-fp26-tm1775-7e8fb6df | 42653927 | 24, 36, 4 | 2.49e4 | 2.71e8 | 6715 (2263–15336) |

Slopes of `log₂ cost` against `log₂ r` (OLS over runs; 95% bootstrap over runs within
curves; [`fits.json`](reports/run-20261012T0009Z/fits.json)):

| series | slope | 95% bootstrap |
|---|---|---|
| rho.negation, field ops (n = 192) | 0.377 | 0.350–0.404 |
| rho.negation, group ops (n = 192) | 0.400 | 0.371–0.430 |
| zero-minor, field ops, all phases (n = 48) | 0.783 | 0.669–0.879 |
| zero-minor, field ops, search phases only (n = 47) | 0.905 | 0.769–1.021 |
| zero-minor, candidate `d`-subsets covered (n = 47) | 0.929 | 0.793–1.048 |
| zero-minor, `d = 3` runs only, 16–24 bits (n = 40) | 0.631 | 0.542–0.717 |

**Identity checks.**

- All 48 zero-minor and 192 rho runs recovered the planted logarithm. Each answer was
  checked by `[m]G = Q` outside the solvers.
- A second full execution reproduced all 240 rows with identical counts, seeds and
  recovered values (`replay-compare.txt`).
- The defect sweep (`d = 2, 3, 4` on the 20-bit curve, 8 targets each, all verified)
  gave identical per-run counts under two binaries (`c70d57ec7`, `06b032f0c`). Its
  medians were 1.21e7, 2.54e6 and 3.32e6 ops. The paper's rule picks `d = 3`, and
  `d = 3` has the lowest median.

**Tests** (`cargo test --release --example lit2026_w12_zero_minor`, 7 tests):

- the ladder curves match the registry ICV1 slugs;
- `dim K = ℓ`, with the anti-diagonal form;
- Theorem 1 / Theorem 3.1: the maximal minor equals `± det A[R, C]`, checked by
  determinant;
- planted instances are solved end to end at 12 bits;
- the default parameters work at 14 bits;
- the left kernel of the signature step is correct;
- Eq. (4) is monotone in `d`;
- rho solves.

### Interpretation (labelled; not part of the measurement)

- Every zero maximal minor of `K` is an `ℓ`-subset of the `2ℓ` points that sums to
  `O`. Each candidate is such a subset with probability about `1/r`.
- The hyperplane and signature machinery tests candidates in bulk with linear
  algebra. It does not change that probability. So the expected number of candidates
  is `Θ(r)` and the cost is linear in `r`, up to the `poly(log r)` factor per
  candidate.
- The measured covered/`r` medians (0.37–2.1) and the slope of the search phases near
  1 are what this predicts.
- Fitted slopes below 1 for the whole method reflect the fixed set-up (the kernel of
  `M`, 51% of the ops at 16 bits). They also reflect the jump from `d = 3` to `d = 4`
  at 26 bits.
- Nothing measured here suggests a sub-rho method. This is consistent with the
  paper's own Table 2 and §1 ("no better than the exhaustive search"; "poses no
  security risk").

## Limits

- **Sizes.** The workstream title says 20–30 bits; the brief says 16–32.
  - The measured ladder is 16–26 bits. Its 26-bit curve has `r ≈ 2^25.35`. It is the
    largest prime-order curve registered under `find_prime_order_curve(·, 59297)`.
  - Larger sizes would need new curves. `find_prime_order_curve` counts points in
    `O(p)`, which is impractical at 28–32 bits on this host, and new curves would also
    need registering.
  - Extrapolating the measured cost to 30 bits gives about 32× the 26-bit cost per
    target, roughly 15–30 CPU-minutes. That is affordable, but it was not run.
  - **Partial.**
- **One seed per zero-minor target.** Eight targets per curve. The intervals are
  bootstrap intervals over those runs.
- **The rho reference's field-op slope (0.377) is below its asymptotic 0.5.** The
  walk's fixed set-up (jump table, start stride) weighs on 16–26-bit curves. The
  comparison is per target at each size, so it does not depend on the slope.
- **Wall times are L0** on a loaded shared host (load 79–106). They are reported only
  as medians in the run report and are not used for any conclusion.
- **The `original` arm** (program baseline `4b6e09b78`) cannot run this method; the
  example is new. Journey entries therefore carry `original: null`.

## Exact commands

```bash
export CARGO_TARGET_DIR=/Volumes/SSD990-2/crypto/targets/lit2026-w12
export TMPDIR=/Volumes/SSD990-2/crypto/targets/lit2026-w12/tmp   # system disk was full
cargo test --release --example lit2026_w12_zero_minor
cargo build --release --example lit2026_w12_zero_minor
cd research/lit2026_w12_flagged_claims_20261011/reports/run-20261012T0009Z
lit2026_w12_zero_minor run --bits 16,18,20,22,24,26 --targets 8 --rho-seeds 4 --out scaling.json
lit2026_w12_zero_minor dsweep --bits 20 --ds 2,3,4 --targets 8 --out dsweep.json
lit2026_w12_zero_minor refit --in scaling.json --out fits.json
lit2026_w12_zero_minor run ... --out <scratch>/replay.json
lit2026_w12_zero_minor compare --a scaling.json --b <scratch>/replay.json
python -I ../../plot_scaling.py .     # layout only
```

**Commits.**

- Branch base: `3c239129d`, the program commit stacked on #1681.
- Measured solver: `c70d57ec7`.
- Analysis modes and fit guard: `06b032f0c`.

**Host.** Shared Apple M4 Pro, macOS (Darwin 25.6.0), isolation L0.

## Canonical graphs checked and left unchanged

- `docs/ic/progress-timeline.json`
- `docs/index-calculus-scoreboard.html`
- `docs/performance-gains/`

None changes. This is an audit of an external claim on a method that is not index
calculus and not a repository pipeline, and it produced no new incumbent and no
speedup.

## Journey entries

`research/lit2026_speedup_program/journey/entries/20261011-W12-*.json`:

- one `negative` entry per measured curve (6);
- one `diagnostic` entry for the slope fit;
- one `negative` entry for Hu 2025/015 (blocked).

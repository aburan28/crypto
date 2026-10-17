# Audit of the generic-DLP code: findings and status

Audit of 2026-10-17, over the code that solves the elliptic-curve discrete
logarithm without index calculus and the arithmetic under it:

- BSGS: `bsgs_fast`, `ecdlp_variants`, and ecbench's `generic.rs`;
- Pollard rho: `prime_fast`, `pollard_rho`, `pollard_collab`,
  `aut_folded_rho`, `preprocessing_rho`, `koblitz_strong_rho`, and
  `ic_boundary`'s walk;
- the kangaroo and Gaudry–Schost;
- Pohlig–Hellman;
- the `ecc/` field, point and constant-time stacks.

Line numbers are as of `main` at `b1bf5b893`.

Each finding says how it is known:

- **measured**: by this study's runs;
- **verified**: by reading the code against a worked case, or by a test in
  the fixing change;
- **reported**: found by a read-only audit pass and not yet reproduced
  here. A reported cost figure is an estimate, not a measurement. Some came
  from scratch models outside the repository and are not evidence.

Status is one of: *fixed in #2048*, *fixed in the correctness PR* (pending,
see below), or *open*. Open findings are not claims. They are leads, each
with the change that would settle it.

## Fixed in #2048 (this study, engineering, outputs pinned)

| # | where | finding | how known | effect |
|:--|:--|:--|:--|:--|
| 1 | `bsgs_fast.rs` kernels and `BabyTable` | BSGS on tables past L2 waited on probe misses. At load ≈ ½ the probe's "empty?" branch depends on the missing load, so misses serialised. | measured | 1.30–1.98× lower one-target wall time at 44–52 bits ([README](README.md)) |
| 2 | `prime_fast.rs:540` `mul`, `:554` `combine` | Affine double-and-add paid one binary-GCD inversion per group operation. `rho_parallel` calls `combine` for its 1024-entry table and every walk start, all outside `total_steps`. | measured | start-up 6.8–16.1× cheaper; whole solve 1.8–5.8× faster at 32–40 bits (1 thread) |

## Fixed in the correctness PR (branch `claude/generic-dlp-correctness-fixes`, pending)

These change counted outputs by design: walk transcripts, step counts. They
are therefore re-pinned deliberately in `tests/ec_vartime_walks_pinned.rs`,
with before/after counts, in their own PR.

| # | where | finding | how known | effect |
|:--|:--|:--|:--|:--|
| 3 | `pollard_collab/walk.rs:115` | **The negation-map 2-cycle escape never fires.** `last == Some(&np)` is tested before `np` is normalised. In a fruitless 2-cycle `R0 → R1 = −(R0+M) → −R0`, the unflipped sum is `−R0`, which never equals the normalised `R0`, so walks loop until `step_cap` and end as dead trails. | verified (worked case; the fix adds `fruitless_two_cycles_are_escaped_deterministically`) | negation mode usable at real `dp_bits`; escape made deterministic (smaller-x member doubled) |
| 4 | `pohlig_hellman.rs:157` | Each digit is a linear search, `O(q)`, while the module claims `O(√q)`. | verified | per-prime BSGS with the table shared across the `e` digits: `≤ ⌈√q⌉(e+1)` points, same digits |
| 5 | `pohlig_hellman.rs:85` via `j0_twists::factorise_small` | The whole order is trial-divided in `BigUint` up to √(cofactor), before the smoothness bound is consulted. For an order with a large prime factor (a cryptographic `curve.n`) that does not finish. | verified | trial division bounded by `smoothness_bound` |
| 6 | `pohlig_hellman.rs:110` | The CRT result is never checked against `Q`, and order 1 returns `None`. | verified | `[d]G = Q` checked; order 1 gives `d = 0` |

## Open

### Correctness and claim fidelity

| # | where | finding | how known | what would settle it |
|:--|:--|:--|:--|:--|
| 7 | `ecc/ct.rs:39` `scalar_mul_secret`; `ecc/point.rs:347` `scalar_mul_ct` | **Coarse timing leak in protocol code.** Curves other than secp256k1 and P-256 fall back, unblinded, to the affine `BigUint` ladder. It starts at `r0 = O`, where `O + P` and `2·O` are shortcuts, so each iteration above the scalar's top set bit costs less and wall time reveals the scalar's bit length. The nonce bit length is what lattice attacks use. The doc says the coarse operation-count leak is closed. Callers: `sm2.rs`, `gost_3410_2012.rs`, `ec_kcdsa.rs`, `asymmetric/ec_elgamal.rs`. (`SECURITY.md` already disclaims limb-level timing.) | verified by reading; leak size not measured | Blind as `scalar_mul_secret_blinded` does, or fix the bit length (`k + n` or `k + 2n`), plus a test that times leading-zero scalars against full-length ones. Better: run the RCB general-`a` formulas over `ct_bignum`'s runtime Montgomery context. |
| 8 | `ecdlp_variants/gaudry_schost.rs:118` `branch` | The branch is `x`'s low byte mod `r`, and the DP test (`:131`) zeroes the low `dp_bits`. With `r = 32` and `dp_bits = 4`, branches 0 and 16 are never taken from a non-DP point: a 30-way walk. | verified by reading | Take the branch from bits above the DP bits, then re-pin the GS transcripts. Impact on the iteration count is expected to be small, but it is not measured. |
| 9 | `pollard_collab/walk.rs:79` `branch_index` | The same overlap: `x mod r` against a low-bits DP test. Branch 0 is never taken from a non-DP point once `dp_bits ≥ 5`. | reported | as #8 |
| 10 | `ecdlp_variants/bsgs.rs:326` (#8), `:390` (#9) | The negation-map interleaved and grumpy variants use stride `m`, not `2m − 1`, so the x-folded tables do not widen the giant stride. They pay for folding and do not get its coverage. | verified by reading | Stride `2m − 1` with tables of `m` entries, the counts compared against #3 and #4 over a target panel |
| 11 | `ecdlp_variants/bsgs.rs:214` `grumpy_giants` | Giants `Q − iS` and `Q + iS` (`m = ⌈√n⌉`) cannot meet each other, so this is not Bernstein–Lange's three-walk search (giants `Q + i·mG` and `2Q − k(m+1)G`). The doc's `≈ 1.25√n` is BL's figure, not this code's. | reported (the auditor's scratch model gave 1.59√n, which is not evidence) | A native panel of targets measuring the mean ops/√n of #4 as written, then BL's construction against the same panel |
| 12 | `ecdlp_variants/gaudry_schost.rs:262–307` | Tame and wild walks both start uniformly over `[0, n)`. In the negation variant, wild–wild collisions with opposite signs are discarded, though they solve the DLP. | reported | One table; recover whenever `dv ≠ 0`; a native panel |
| 13 | `aut_folded_rho.rs:473`, `:519` | The walk steps the raw point and uses the canonical form only for the bucket, so walks that meet in an orbit do not merge. The `√6` saving does not happen. The one test that would see it is `#[ignore]` and asserts only "under 4× unfolded". | reported | Measure folded against unfolded mean ops on a panel; or replace the solver with `rho_walk_with(CubicClasses, RhoWalk::negation())`, which already runs as ecbench's `rho.cubic` |
| 14 | `pollard_rho.rs:265` | The Floyd result is returned without checking `g^x = h` when the gcd is 1. `RhoSolution.iterations` counts iterations (3 operations each), not operations as documented. | reported | One exponentiation before returning; fix the doc or the unit |
| 15 | `preprocessing_rho.rs:246`, `:311` | The online walk never restarts, so it can cycle for all `max_online_steps`. `expected_online_cost` ignores `dp_bits`. | reported | A restart cap of about `20·2^d`; the online cost formula with DPs |

### Cost (counted or wall), all on reference arms

Every item in this table makes a **reference** cheaper: rho, BSGS or the
kangaroo. That widens the gap an index calculus has to close. It is not an
index-calculus improvement, and it needs a fresh preregistered calibration
before any matched-rho figure moves. `rho.frozen_reference` stays as it is.

| # | where | finding | how known | what would settle it |
|:--|:--|:--|:--|:--|
| 16 | `ic_boundary.rs:1254` `rho_jumps_for` | The jump count is capped at 16; only {4, 8, 16} were calibrated. The departure from a random walk costs about `1/√(1 − 1/r)` (Teske; literature). | reported | A calibration of 32 and 64 jumps at `2^28` to `2^40`, mean steps over a panel |
| 17 | `ic_boundary.rs:3435` `PrimeCurve::add_raw` | ecbench's prime rho is one affine addition with an extended-Euclid inversion per step, with `%` reduction. Index-calculus relation walks already batch through `prime_fast::BatchAddScratch`; rho does not. | reported | Lanes plus simultaneous inversion. `S` is unchanged (GAE); field operations and wall time fall. |
| 18 | `koblitz_strong_rho.rs:464`, `:1206`, `:1251` | `dp_bits = 4` makes walk starts about 6.7% of `S` at `n = 41`. Strides that are multiples of `G` alone waste collisions with probability about `1/(2n)`. Setup dominates at small `n`. | reported (from committed counters) | Size-dependent `dp_bits` and a stride `sG + tQ`, measured as an ecbench session |
| 19 | `ecbench/generic.rs:203` (kangaroo) | The power-of-two jump mean lands in `(target/2, target]`, so the expected cost is `√r(f + 1/f)`, 2.0–2.5√r. | reported | Rescale the jumps so the mean hits the target, measured over a panel |
| 20 | `prime_fast.rs:846` (rho step) | Each step updates `(a, b)` mod `n` with two `add_mod` and the negation's `neg_mod`. Re-walking from the seeds on a collision would drop that work. | derived (this study) | Engineering with pinned output: identical `total_steps` and log at one thread |
| 21 | `bsgs_fast.rs:386` `FastCurve::scalar_mul` | Affine double-and-add with a Fermat inversion per operation, used for chain seeding and verification. Small at 40–52 bits; dominant on tiny intervals and Pohlig–Hellman subproblems. | derived | The same Jacobian ladder as #2 |
| 22 | one-target BSGS | Non-interleaved negation-map BSGS costs `≈ 1.00√n` on average. Interleaving with stride `√(2n)` gives `≈ 0.943√n`, at about 2.8× the table memory. | derived | A new protocol: a counted-cost change, and multi-target reuse is lost |

### Protocol-code speed (not on any ECDLP hot path)

| # | where | finding | how known |
|:--|:--|:--|:--|
| 23 | `secp256k1_point.rs:265`, `p256_point.rs:241`, `point.rs:356` | Each ladder bit computes two doublings; a `cswap` ladder needs one (about 27–30% per bit). | reported |
| 24 | `p256_point.rs:148–201` | P-256 uses general-`a` formulas although `a = −3` (RCB Algorithms 4 and 6). | reported |
| 25 | `secp256k1_field.rs:150`, `p256_field.rs:109` | Inversion is a cmov ladder over the public exponent `p − 2` (768 operations). A fixed addition chain is equally constant-time and about 2.8× cheaper. | reported |
| 26 | `ct_bignum.rs:374` | `mont_sqr` is `mont_mul`. A dedicated squaring saves 6 of 16 limb products. | reported |

## What the audit checked and found correct (reported, not re-derived)

The audit found the following correct; this study did not re-derive them:

- the RCB complete formulas in `p256_point.rs` and `secp256k1_point.rs`;
- `ct_bignum::mont_mul`'s carries;
- the Jacobian `dbl-2007-bl` and `madd-2007-bl` formulas in `point.rs`;
- every `O`, `P + (−P)` and `y = 0` path;
- `fp_backend`'s REDC bounds, `FpP192::reduce7`, the CIOS carries and the
  IFMA wide reduction;
- `crt_combine`;
- the BSGS index arithmetic.

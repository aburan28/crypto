# ECDLP across every NIST curve

`src/cryptanalysis/ecdlp_nist/` and the `ecdlp-nist` binary: one front end
that takes any of the fifteen FIPS 186-4 curves — `P-192 … P-521` over
`F_p`, `K-163 … K-571` (Koblitz, the "anomalous binary curves") and
`B-163 … B-571` (pseudo-random) over `F_{2^m}` — audits it for the
structural conditions each known ECDLP attack needs, and runs the branch
the curve calls for.

## What runs where

| branch | precondition | module | NIST curves |
|---|---|---|---|
| Smart–Semaev–Satoh–Araki | `#E(F_p) = p` (trace 1) | `anomalous.rs` | none is anomalous; demonstrated on CM-constructed anomalous curves of any size (`anomalous-<bits>-<seed>`) |
| Pohlig–Hellman | generator order composite | `pohlig.rs` | every `n` is prime; demonstrated on toy full groups (`toy-k19a0-full`) |
| MOV / Frey–Rück | embedding degree `k ≤ 100` | audit only (the `k = 2` transfer is `mov_attack.rs`) | no curve has `k ≤ 100` |
| Weil descent / GHS | composite extension degree `m` | audit only (`ghs_*`, `binary_ecc::cover_transfer`) | every `m` is prime |
| BSGS, parallel kangaroo | scalar in a known interval | `interval.rs` | runs on all fifteen at full width (tests: `2^16`–`2^26` wide) |
| parallel rho | otherwise | `rho.rs` | the real cost; refused as infeasible unless forced |

The rho folds the negation map on every curve and the Frobenius
`τ(x, y) = (x², y²)` on Koblitz curves: the walk moves on classes
`{±τ^i(P) : i < m}` of size `2m`, carrying the coefficients `(a, b)` of
`P = [a]G + [b]Q` exactly through the eigenvalue `λ` (`τ = [λ]` on the
subgroup, `λ` the root of `λ² − μλ + 2 ≡ 0 mod n` matched against
`τ(G)`).  Expected additions are `√(πn / (2·class))`:

| curve | `n` bits | class | log₂ expected (unfolded) | log₂ expected (folded) |
|---|---:|---:|---:|---:|
| P-192 | 192 | 2 | 96.3 | 95.8 |
| P-224 | 224 | 2 | 112.3 | 111.8 |
| P-256 | 256 | 2 | 128.3 | 127.8 |
| P-384 | 384 | 2 | 192.3 | 191.8 |
| P-521 | 521 | 2 | 260.8 | 260.3 |
| K-163 | 163 | 326 | 81.3 | 77.2 |
| K-233 | 232 | 466 | 115.8 | 111.4 |
| K-283 | 281 | 566 | 140.8 | 136.3 |
| K-409 | 407 | 818 | 203.8 | 199.0 |
| K-571 | 570 | 1142 | 284.8 | 279.7 |
| B-163 | 163 | 2 | 81.3 | 80.8 |
| B-233 | 233 | 2 | 116.3 | 115.8 |
| B-283 | 282 | 2 | 141.3 | 140.8 |
| B-409 | 409 | 2 | 204.3 | 203.8 |
| B-571 | 570 | 2 | 285.3 | 284.8 |

(`ecdlp-nist list` prints this table from the live parameters; the bit
counts are those of the subgroup order `n`, so `K-233` shows 232 and
`B-283` shows 282.  Cofactors: `1` on `P-*`, `2` on `B-*` and `K-163`
(`a = 1`), `4` on the other four `K-*` (`a = 0`); the audit also checks
that every `n` is prime, no embedding degree is `≤ 100`, no trace is `1`,
every binary `m` is prime, and each Koblitz `h·n` equals the Lucas-recurrence
order.)

## The Smart attack, end to end

`canonical_lift.rs` carried the Hensel lift and the formal logarithm but
returned `Err` from `smart_attack_anomalous` because affine arithmetic
modulo `p²` cannot represent `[p]P̂`, which is the point at infinity modulo
`p`.  `ecdlp_nist/anomalous.rs` closes that: `[p]P̂` and `[p]Q̂` are
computed in Jacobian coordinates over `Z/p²Z` (ring operations only, so
precision is exact), their `Z` coordinates are divisible by `p`, the
formal parameter `t = −x/y = −XZ/Y` is `≡ ψ(t) mod p²`, and
`k = (t_Q/p)(t_P/p)⁻¹ mod p`.  The lift `b̃ = b + jp` is retried for the
next `j` when `t_P/p ≡ 0` (the canonical lift, probability `1/p`), and
every answer is verified by `[k]P = Q` on `E(F_p)`.  The old entry point now
delegates to this one.

Anomalous test curves are built without point counting: for the
class-number-one discriminants `D ∈ {−11, −19, −43, −67, −163}` (those `≡ 5 mod 8`, so `p` is odd),
`p = (1 + |D|v²)/4` prime gives Frobenius trace `±1`, so the CM curve
`y² = x³ + 3j(1728−j)x + 2j(1728−j)²` or its quadratic twist has exactly
`p` points; `[p]P = O` on a random point picks the right one.  The
`ecdlp-nist smart --bits 256` command constructs one and breaks it.

## Commands

```text
ecdlp-nist list                                   # audit table of the fifteen curves
ecdlp-nist audit K-233 [--json]                   # one curve, every check
ecdlp-nist solve P-256 --secret 0x1234567 --interval-bits 32
ecdlp-nist solve B-571 --random-bits 24 --method kangaroo --threads 4
ecdlp-nist solve K-283 --point <x>,<y> --interval-bits 40   # external point, bounded
ecdlp-nist solve toy-k23a1 --random-bits 21 [--no-fold]      # whole-group rho, Frobenius folded or not
ecdlp-nist solve toy-k19a0-full --secret 99999               # Pohlig–Hellman on a composite order
ecdlp-nist smart --bits 256 --seed 7                         # construct an anomalous curve, break it
ecdlp-nist bench P-521                                       # additions per second on this host
```

`solve` picks the branch from the audit unless `--method` names one.  A
whole-group rho whose expected cost exceeds `2^40` additions is refused;
`--force --max-iterations N` runs it to the budget and reports the budget
exhausted.  Exit status is non-zero when no scalar was found.

## Measurements (this container, release build)

Filled from the runs recorded in the PR.  Bounded instances are planted
known-answer scalars; the arithmetic is the crate's affine `BigUint` /
`F2mElement` group law, not a tuned field implementation, so the rates
below are what this code does, not what the curve admits.

Build: `cargo build --release`; one container vCPU class, `--threads 4`
where stated, a clippy run sharing the machine during the interval and rho
timings (so treat wall-clock as indicative; the operation counts are exact).

**Affine additions per second, one thread** (`ecdlp-nist bench`):

| curve | additions / s |
|---|---:|
| P-256 | 31 900 |
| P-521 | 11 300 |
| K-233 | 86 000 |
| K-571 | 20 700 |
| B-571 | 20 900 |

**Interval solves on real curves** (`solve <curve> --random-bits b --interval-bits b`, seed 11, 4 threads; `ops` are group additions, `exp` the textbook `2√(W/2)` for BSGS and `2√W + N·2^dp` for kangaroo):

| curve | width | method | ops | exp | wall |
|---|---:|---|---:|---:|---:|
| P-256 | 2^20 | bsgs | 1 441 | 1 448 | 43 ms |
| P-256 | 2^28 | bsgs | 22 139 | 23 170 | 689 ms |
| P-256 | 2^32 | kangaroo | 63 008 | 133 120 | 810 ms |
| K-283 | 2^20 | bsgs | 1 441 | 1 448 | 50 ms |
| K-283 | 2^28 | bsgs | 22 139 | 23 170 | 341 ms |
| K-283 | 2^32 | kangaroo | 26 436 | 133 120 | 181 ms |
| B-409 | 2^20 | bsgs | 1 441 | 1 448 | 48 ms |
| B-409 | 2^28 | bsgs | 22 139 | 23 170 | 637 ms |
| B-409 | 2^32 | kangaroo | 724 736 | 133 120 | 8.7 s |

Kangaroo cost has a heavy tail: six seeds on P-256 at width `2^32` took
111 592, 281 020, 126 308, 292 212, 125 692 and 389 120 additions (mean
≈ 1.7 × the `2√W + N·2^dp` figure; the `B-409` row above is a 5 × outlier).
The planted scalar is recovered and verified in every run.  All fifteen
curves are covered at width `2^16` by both methods in the unit tests and at
`2^26` by the CLI tests.

**Whole-group rho on toy Koblitz curves, Frobenius folding on and off**
(`solve toy-k… --random-bits 40 --threads 4`, three seeds each):

| curve | `n` | class | expected | measured additions (3 seeds) |
|---|---:|---:|---:|---|
| toy-k23a1 | 4 196 903 | 46 | 379 | 326, 359, 617 |
| toy-k23a1 (`--no-fold`) | 4 196 903 | 1 | 2 568 | 1 681, 4 089, 3 338 |
| toy-k41a0 | 549 756 390 943 | 82 | 102 621 | 144 994, 114 753, 117 060 |
| toy-k41a0 (`--no-fold`) | 549 756 390 943 | 1 | 929 277 | 1 297 945, 424 371, 1 222 044 |

Folding buys the predicted `√(2m)` (≈ 6.8 × at `m = 23`, ≈ 9 × at
`m = 41`) in additions; each folded step also pays `m` squarings for the
orbit minimum, which at these sizes costs about a quarter of an addition
(0.45 s vs 2 s wall at `m = 41`).

**Smart attack on CM-constructed anomalous curves** (`ecdlp-nist smart --bits b --seed 3`):

| `p` bits | construct | attack | lifts tried |
|---:|---:|---:|---:|
| 129 | 2 ms | 3 ms | 1 |
| 257 | 2 ms | 10 ms | 1 |
| 522 | 16 ms | 33 ms | 1 |

Polynomial time, as the theory says: the whole run is a few scalar
multiplications modulo `p²`.


## Scope and claims

* Nothing here recovers a full-width scalar on a NIST curve, and the
  reports say so.  Interval solves exercise the full-width arithmetic and
  the bookkeeping on real parameters; whole-group searches run only on toy
  or constructed curves.
* The audit certifies the *absence* of each structural precondition on the
  fifteen curves (trace ≠ 1, prime `n`, no embedding degree ≤ 100, prime
  `m`); it is not a security proof, and the generic cost it reports is the
  textbook expectation, not a measured bound (`docs/bounds/`, `ecbench`
  are the measured-cost tooling).
* "Anomalous" in the NIST naming (`K-xxx`, anomalous binary curves) means
  Koblitz; the `p`-adic Smart attack is for trace-1 prime-field curves and
  the audit reports both facts separately.
* Koblitz orders are cross-checked against the Lucas recurrence from the
  trace over `F_2`; `B-163` is `sect163r2`, per FIPS 186-4.

## Tests

`cargo test --release --lib cryptanalysis::ecdlp_nist` (unit: folding is a
class function with exact multipliers, rho with and without Frobenius,
BSGS/kangaroo on every NIST curve, Smart on brute-forced tiny curves and
CM curves from 40 to 256 bits, Pohlig–Hellman on composite toy groups,
audit invariants on all fifteen), `cargo test --release --lib
cryptanalysis::canonical_lift` (the delegating entry point recovers every
`d` on tiny anomalous curves), `cargo test --release --test ecdlp_nist`
(the binary across curves, with JSON).

# Lead 16.2: a fast 4E test as a subgroup check (defender side)

Round 4 of `RESEARCH_ISOGENY_FIBERS_ECDLP.md` (§14.2 Prop. 14.2, §16.2).
Date 2026-10-11. Runs started 2026-10-11T19:50:38Z, worktree commit
`440c0e163` (branch `research/isogeny-fibers-round4`, plus the untracked
round-4 files). Host: x86_64, 4 cores shared with about 6 other agents. Every
run was single-threaded.
Code: `examples/r4_16_2_fast_4e.rs`.

## Question, hypothesis, falsifier (from §16.2)

* **Question.** Take the a = 0 Koblitz curve y² + xy = x³ + 1 over F_{2^n},
  n odd, with #E = 4r and r prime. Does the two-trace closed form of
  Prop. 14.2 decide membership in the order-r subgroup G = 4E exactly? How
  much cheaper is it than the [r]P = O check?
* **Hypothesis.** The test is exact. It costs one multiplication, O(n)
  squarings (or table lookups) and three traces. [r]P costs hundreds of
  multiplications.
* **Falsifier.** Any point where the closed form and [r]P = O disagree.

## 1. The test, re-derived (derived; then verified exhaustively and by sampling)

The group law is the library's (`src/binary_ecc/curve.rs`). Doubling Q = (u, v), u ≠ 0, gives

    λ = u + v/u,   x = λ² + λ,   y = u² + (λ + 1)·x.

**Halving.** P = (x, y) is 2Q for a rational Q iff λ² + λ = x has a root in
F_{2^n}, iff Tr(x) = 0 (Hilbert 90 for Artin–Schreier). For odd n the
half-trace H(x) = Σ_{i=0}^{(n−1)/2} x^{2^{2i}} satisfies H(x)² + H(x) = x + Tr(x).
So when Tr(x) = 0 the two roots are λ = H(x) and H(x) + 1. The two
rational halves are Q and Q + T, where T = (0, 1) is the only rational
2-torsion point. Each half has u² = y + (λ + 1)x.

**Second level.** Q ∈ 2E iff Tr(u) = 0, applying the same criterion to Q. Now
Tr(u) = Tr(u²) = Tr(y) + Tr(λx) + Tr(x), and Tr(x) = 0, so
Tr(u) = Tr(y) + Tr(x·H(x)) for λ = H(x). For λ = H(x) + 1 the value changes
by Tr(x) = 0, so the two halves agree. That has to happen, because
T = 2·(1, 0) ∈ 2E. P ∈ 4E iff some half lies in 2E. Therefore:

> **P ∈ 4E ⟺ P = O, or [Tr(x) = 0 and Tr(y) + Tr(x·H(x)) = 0].**
> The traced elements are x, y, and the product x·H(x).

This matches §14 Prop. 14.2.

**Special points.**
* O is accepted.
* T = (0, 1) needs no special case: Tr(0) = 0, then Tr(1) + Tr(0) = n mod 2 = 1,
  so T is rejected. That is correct, since T has order 2.
* The order-4 points (1, 0) and (1, 1) have Tr(1) = 1 and are rejected at the
  first level.
* The halving formula needs u ≠ 0. A half with u = 0 is T, and T doubles to O,
  which is not affine, so the formula covers every affine P.

**Why 4E = G.** For odd n, (1, 0) has order 4 and Tr(1) = 1, so it has no
half. So the 2-Sylow subgroup of E(F_{2^n}) is exactly Z/4 for every odd n.
This is confirmed by `two_part_exactly_4` for n = 5…23 and by #E ≡ 4 (mod 8)
at n = 131 and n = 233. When #E = 4r with r odd, 4E is the order-r subgroup.
When r is composite (n = 11, 15, 17, 21), the test still decides 4E exactly,
but 4E is no longer of prime order.

**Preconditions.**
* n is odd. For even n there is no half-trace, and Tr(1) = 0 changes the
  2-part.
* a = 0. For a = 1 the criterion reads Tr(x) = Tr(a). That curve has
  cofactor 2, and the second level does not apply.
* P is already known to lie on the curve. The on-curve check (2M + 2S) is
  needed by both methods and is not counted.

## 2. What ran (exact commands)

```sh
export CARGO_TARGET_DIR=/tmp/claude-0/-home-user/717e4727-b5a7-5b25-924a-7b31edb43331/scratchpad/r4target
cargo build --release --example r4_16_2_fast_4e
B=$CARGO_TARGET_DIR/release/examples/r4_16_2_fast_4e
$B sample --n 131 --uniform 600000 --balanced 100000 --lib 2000 --seed 131 --bench-batch 2048 --bench-reps 5 --out sample_n131.json   # > sample_n131.txt
$B sample --n 233 --uniform 600000 --balanced 100000 --lib 2000 --seed 233 --bench-batch 2048 --bench-reps 5 --out sample_n233.json   # > sample_n233.txt
$B exhaustive --nmax 23 --out exhaustive.json                                                                                          # > exhaustive.txt
```

A first exhaustive run with `--nmax 25` finished n = 5…23. It was stopped
during n = 25 to stay inside the brief's limit of about 10 minutes per run.
Its log is `exhaustive_nmax25_stopped_during_n25.txt`, and n = 25 was not
checked. The run was repeated with `--nmax 23` to get the JSON.

**Implementation.**
* The example carries its own fixed four-word GF(2^n) field: carry-less
  `pclmulqdq` products and generic sparse-modulus reduction. It is
  cross-checked against `crypto_lib::binary_ecc::F2mElement` on 2000 random
  multiplications and squarings, and the half-trace identity is checked too
  (`field_vs_library_and_halftrace_identity: true`).
* Moduli: n = 131 uses `IrreduciblePoly::deg_131` (z^131 + z^8 + z^3 + z^2 + 1).
  n = 233 uses `deg_233` (z^233 + z^74 + 1) with the library's sect233k1
  generator and order. Small n uses the first irreducible trinomial or
  pentanomial found by a Ben-Or test.
* Methods compared, all on the same field code:
  1. **[r]P, binary López–Dahab (LD)**: left-to-right double-and-add, with LD
     doubling (3M + 5S) and mixed addition (8M + 5S).
  2. **[r]P, τ-NAF LD**: r is reduced modulo τ^n − 1 in Z[τ] (τ² = −τ − 2),
     which is valid for every point because τ^n − 1 kills E(F_{2^n}). The
     check N(τ^n − 1) = #E holds. Each Frobenius costs 3S.
  3. **Closed form, H by squarings**: n − 1 counted squarings.
  4. **Closed form, H by byte table**: H is F₂-linear, so it is ceil(n/8)
     lookups in precomputed 256-entry tables (17 KB for n = 131, 30 KB for
     n = 233), with no field operations.
* Arithmetic check: the LD and τ-NAF scalar multiplications are checked
  against the library's affine `point_add`, `point_double` and `scalar_mul`
  (`ld_tnaf_vs_library_affine: true`).
* **Ground truth.** The coset k ∈ Z/4 of P in E/4E is read from R = [r]P
  (binary LD). R is one of O, (1,0), T or (1,1), and k = j·r mod 4.
  Independently, the library's `BinaryCurve::is_valid_public_point` (on-curve
  check plus affine [r]P) was run on 2000 points per n.
* **Sample.** Per n:
  * 600 000 uniform points (random x, solve for y, random choice of the two
    y);
  * 100 000 "balanced" quadruples [4]P + {O, P4, T, −P4}, which hit each
    coset exactly 100 000 times;
  * 9 edge points: O, T, ±P4, ±G, G+T, G±P4.

  Total 1 000 009 points per n.

## 3. Measured results

### 3.1 Exactness (measured)

**Exhaustive, every point on the curve, n = 5, 7, …, 23** (`exhaustive.json`). The
enumerated #E equals 2^n + 1 − t_n in every case. |4E| = #E/4 by [#E/4]P = O.
The closed form (both half-trace variants) accepts exactly the members of 4E:
**0 false accepts and 0 false rejects** at every n. The largest case is
n = 23 with 8 383 412 points.

**Sampled, n = 131 (the ECC2K-130 curve) and n = 233 (sect233k1)**, 1 000 009
points each. The table gives the true coset (from [r]P) against the test
verdict:

| n | true coset | in G | test accept | test reject | bit pattern (Tr x, Tr y + Tr xH(x)) |
|--:|:--:|:--:|--:|--:|:--|
| 131 | 0 | yes | 249 897 | 0 | (0,0) ×249 897 |
| 131 | 1 | no | 0 | 249 504 | (1,0) ×249 504 |
| 131 | 2 | no | 0 | 250 049 | (0,1) ×250 049 |
| 131 | 3 | no | 0 | 250 559 | (1,1) ×250 559 |
| 233 | 0 | yes | 250 584 | 0 | (0,0) ×250 584 |
| 233 | 1 | no | 0 | 249 533 | (1,1) ×249 533 |
| 233 | 2 | no | 0 | 249 993 | (0,1) ×249 993 |
| 233 | 3 | no | 0 | 249 899 | (1,0) ×249 899 |

* **Errors: 0 of 1 000 009** at each n.
* The τ-NAF and binary [r]P agree on all points (0 disagreements). The
  squaring-chain and table variants agree on all points.
* The library oracle agrees on 2000 of 2000 points per n.
* r(131) = 680564733841876926932320129493409985129. It is a probable prime,
  and it equals the ECC2K-130 order quoted in the literature. That was
  cross-checked by web search only; the paper was not opened.
* r(233) equals the library's sect233k1 order. It is a probable prime, with
  cofactor 4.

**Observation (measured, not derived).** On these samples, the pair of bits
also separates the two cosets of order 4. When Tr(x) = 1, the second
expression is still computed, with H(x) a root of λ² + λ = x + 1. It came out
constant on each odd coset, but which of the two patterns goes to coset 1
swaps between n = 131 and n = 233. It may be a full Z/4-valued invariant;
this is not proved here.

### 3.2 Counted field operations per check (measured; members of G)

The counts are constant over all 249 896 (n = 131) and 250 583 (n = 233)
members of G. For binary [r]P the minimum equals the maximum.

| n | method | M | S | traces | H-table lookups | inversions |
|--:|:--|--:|--:|--:|--:|--:|
| 131 | [r]P binary LD | 661 | 816 | 0 | 0 | 0 |
| 131 | [r]P τ-NAF LD | 346 | 609 | 0 | 0 | 0 |
| 131 | closed form, H by squarings | 1 | 130 | 3 | 0 | 0 |
| 131 | closed form, H by table | 1 | 0 | 3 | 17 | 0 |
| 233 | [r]P binary LD | 1255 | 1506 | 0 | 0 | 0 |
| 233 | [r]P τ-NAF LD | 618 | 1085 | 0 | 0 | 0 |
| 233 | closed form, H by squarings | 1 | 232 | 3 | 0 | 0 |
| 233 | closed form, H by table | 1 | 0 | 3 | 30 | 0 |

**Analytic check (derived).**
* Binary: 3(L−1) + 8(w−2) + 2 multiplications and 5(L−1) + 5(w−2) + 1 squarings,
  with L = bitlen(r) and w = popcount(r). That is (L, w) = (130, 36) and
  (232, 72), giving 661/816 and 1255/1506, which matches the measurement. The final
  addition (r−1)P + P exits early with the result O after 2M + 1S.
* τ-NAF (lengths 133 and 235, with 45 and 79 nonzero digits): 8(z−2) + 2
  multiplications. There is one fewer Frobenius than len − 1, because the
  last digit is 0 and the Frobenius of O is free. This matches 346/609 and
  618/1085.
* The closed form costs exactly 1M + (n−1)S + 3 traces, matching.
* The JSON field `analytic_members_of_G` holds the simpler formula, which
  charges the final addition in full. It is 6 or 7 operations higher; the
  difference is explained by the early exit above.

**Ratios (derived from the table).**
* Multiplications: [r]P τ-NAF / closed form = **346×** (n = 131) and
  **618×** (n = 233). Binary [r]P gives 661× and 1255×.
* Squarings, τ-NAF against the squaring-chain variant: 609/130 = 4.7× and
  1085/232 = 4.7×. The table variant needs no squarings.
* Over all points, about half of them leave at the first trace, so the
  closed form averages 0.50 M per point.

### 3.3 Contended wall time per check (measured; label: contended)

Median of 5 repetitions over 2048 members of G. This was a single thread on a
4-core host shared with other agents. Minimum and maximum are in the JSON.

| n | [r]P binary LD | [r]P τ-NAF LD | closed form, squarings | closed form, table | library `is_valid_public_point` (2000 pts, 1 pass) |
|--:|--:|--:|--:|--:|--:|
| 131 | 232 µs | 160 µs | 22.8 µs | 0.124 µs | 2.06 ms |
| 233 | 165 µs | 98 µs | 11.8 µs | 0.122 µs | 3.98 ms |

**Wall-time ratios.**
* τ-NAF against the table variant: about 1300× (n = 131) and about 800×
  (n = 233).
* τ-NAF against the squaring chain: about 7× and about 8×.

**Caveats.**
* Contention is visible. The n = 131 rows came out slower than the n = 233
  rows, and the uniform-set runs (see JSON) differ from the members-of-G runs
  by up to 4× for the chain variant. The ratios are therefore order-of-magnitude
  only.
* The example's squaring is a generic carry-less product with a generic
  reduction, not a bit-interleave table with a hand-coded trinomial
  reduction. That makes the squaring-chain variant slower than an optimised
  implementation would be, and it also inflates the S terms of both [r]P
  methods.
* The counted operations are the primary metric.

## 4. Verdict

**Confirmed.** On y² + xy = x³ + 1 over F_{2^n}, n odd, P ∈ 4E iff P = O or
[Tr(x) = 0 and Tr(y) + Tr(x·H(x)) = 0].
* The test was exact on every point for n = 5…23, and on 10⁶ sampled points
  each at n = 131 and n = 233, with zero errors in every coset.
* For cofactor-4 instances (ECC2K-130's curve, sect233k1) this is the
  prime-order subgroup check.
* Counted cost on members of G: 1 multiplication, 3 traces, and either n − 1
  squarings or ceil(n/8) table lookups. A τ-NAF [r]P costs 346 M + 609 S
  (n = 131) and 618 M + 1085 S (n = 233).

**Literature status: not novel as far as checked, and no novelty is claimed.**
* The halvability criterion Tr(x) = Tr(a) is classical: Knudsen 1999,
  Schroeppel, and the point-halving literature of Fong–Hankerson–López–Menezes.
* A web search (2026-10-11) found prior uses of halving as a cofactor-2/4
  validation step:
  * US patent 7,257,709 / US 2004/0114760, "Method and apparatus for
    performing validation of elliptic curve public keys" (priority 2002).
    Search snippets say it validates points on binary curves with cofactor 2
    or 4 by testing whether the point has a half, checking Tr(x) = Tr(a) and
    avoiding a scalar multiplication.
  * B. King, "A point compression method for elliptic curves defined over
    GF(2^n)", PKC 2004, LNCS 2947.
  * Eagle–Galbraith–Ong, "Point compression for Koblitz elliptic curves", IACR
    ePrint 2009/086. The snippet: for #E = 4ℓ, halve P and check
    Tr(x_Q) = Tr(a²).

  That is a two-level halving test, the same criterion as here. The closed
  form only inlines the half's trace, as Tr(u) = Tr(u²) = Tr(y) + Tr(x·H(x)),
  so the square root and the second half-trace are not needed.
* Pornin, ePrint 2022/1164, does subgroup membership by successive halvings
  on twisted Edwards curves. ePrint 2026/749 ("Divide-and-Pair") is a recent
  survey of faster subgroup membership tests.
* Limits of this check: only search-result snippets were read. The patent and
  paper PDFs could not be fetched (HTTP 403 and Cloudflare blocks), so the
  exact formulations were not compared.

## 5. Limits

* Toy fields for the exhaustive check (n ≤ 23). At n = 131 and n = 233 the
  evidence is a 10⁶-point sample plus the derivation, not an enumeration.
* The test assumes an on-curve affine input and odd n. a = 1 curves
  (cofactor 2) and even n were not tested.
* The example is not constant-time. The early exit on Tr(x) and the data-
  dependent table lookups leak which coset a rejected point is in, which is
  public information. On accepted points, the cost and access pattern depend
  on x only through table indices; it was not audited for side channels.
* Wall times are contended, and the field code is not tuned. No `identity.certificate/v1`
  record was produced for the algebraic identity, which `research-visuals`
  requires for an identity claim. The identity is a derivation plus
  exhaustive and sampled checks here.
* `cargo clippy --release --example r4_16_2_fast_4e -- -D warnings` cannot
  pass in this worktree for a reason outside this lead. With `-D warnings`
  the *library* at `440c0e163` fails with 28 clippy errors in `src/` (for
  example `src/cryptanalysis/semaev_decomp.rs:1070`, unnecessary cast), so the
  example is never reached. Running clippy without `-D` reports zero lints
  located in `examples/r4_16_2_fast_4e.rs`. The only example-crate warning is
  the crate-wide `unknown lint: clippy::chunks_exact_to_as_chunks`.
  `rustfmt --edition 2021 --check` passes.

## 6. Next step

* Move the test into `src/` (see the coordinator note) and use it in
  `BinaryCurve::is_valid_public_point` when a = 0, n is odd and h = 4. Keep
  the [r]P path as a debug cross-check.
* Derive (or refute) the observed full Z/4 coset label from the two bits.
* Write the identity certificate.
* Compare against an optimised squaring kernel and against Fong et al.'s
  half-trace-with-odd-bits trick.

## Files

* Inputs and raw outputs: `sample_n131.json`, `sample_n131.txt`,
  `sample_n233.json`, `sample_n233.txt`, `exhaustive.json`, `exhaustive.txt`,
  `exhaustive_nmax25_stopped_during_n25.txt`, `run_started_utc.txt`.
* Graph: `costs.svg`, produced by `plot_costs.py` (matplotlib, recorded values
  only) from `plot_data.json`.
* Report: `report.typ` and `report.pdf`.

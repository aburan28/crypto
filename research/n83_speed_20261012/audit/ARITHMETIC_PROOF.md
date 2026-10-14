# Exact scalar arithmetic and independent audit

The candidate changes the scalar representation inside the configured worker
only. The production `StrongRhoG` implementation and the earlier frozen worker
remain unchanged. Both arms use the same `Gf2_128` field, exact legacy curve,
generator, signed-Frobenius quotient, jump-seed derivation and walk code.

The exact fixture is `icv1-f2m83-t6151469093347-9eba6cd5`: field
`z^83+z^7+z^4+z^2+1`, curve `y^2+xy=x^3+x^2+1`, subgroup order
`r=8569786107849059`, cofactor `1128547018`, and generator
`(6614244651398184860292015,8868260410029249391539529)`. The scalar changes
are constant-factor engineering; they introduce no new ECDLP exponent or
curve transfer.

## Proposition 1: high-product reconstruction

Put `W=2^64`, `B=W^2=2^128`, and decompose two unsigned 128-bit words as
`a=a0+W*a1`, `b=b0+W*b1`, with all four limbs in `[0,W)`.
Distributivity gives

`a*b=p00+W*(p01+p10)+B*p11`,

where `pij=ai*bj`. Writing each `pij=W*hiij+loij`, with its low limb in
`[0,W)`, the upper 128 bits are

`p11+hi01+hi10+floor((hi00+lo01+lo10)/W)`.

This is the worker's `mul_hi_128`. The middle accumulator is below `3W`
and therefore fits in 128 bits. Each partial sum of the final nonnegative
terms is at most the final high word, which is below `B`; no intermediate
addition overflows. The exact floor/carry reconstruction and Rust casts are
supported here by this written proof and independent native tests. They are
not yet a formal proof of the compiled Rust function.

The distributive polynomial expansion is checked in Lean as
`N83ScalarArithmetic.limb_product_expansion`. Certificate
[`IDC1hc23ae1d80fc5fefa`](limb_expansion.certificate.json) also replays that
polynomial over `F_(2^61-1)` at 64 statement-bound points and 64 additional
checker-seeded points. Its `proof` field is null. That evaluated certificate
does not certify integer overflow, the carry/floor formula, or the Barrett
algorithm; the readable integer argument above has the larger scope.

## Proposition 2: the reciprocal is exact and representable

For `2<=m<=2^64-1`, the worker stores `mu=floor(B/m)`. It computes this
without representing `B` in a 128-bit word:

`mu=floor((B-1)/m)+[ (B-1) mod m = m-1 ]`.

Indeed, if `B-1=u*m+v` with `0<=v<m`, then `B=u*m+(v+1)`. The quotient
increases exactly when `v=m-1`. Since `m>=2`, `mu<=B/2`, so the extra one
cannot overflow the 128-bit representation. This specific representability
and successor-division bridge is written and tested; its formalization
status is listed separately in the proof map.

## Proposition 3: one correction computes the exact remainder

Let `m>0`, `B>0`, `0<=x<B`, `mu=floor(B/m)` and
`q=floor(x*mu/B)`. Define remainders `s=B-m*mu` and `t=x*mu-B*q`.
Then `0<=s<m`, `0<=t<B`, and

`B*x=B*q*m+m*t+s*x`.

Since `mu*m<=B` and `q*B<=x*mu`, multiplication and cancellation by
positive `B` give `q*m<=x`. Put `u=x-q*m`. The identity above becomes

`B*u=m*t+s*x`.

The first term is below `m*B`. The second term is at most `m*x`, which is
below `m*B` because `x<B` and `m>0`. Thus `0<=u<2m`.
Subtracting `m` exactly when `u>=m` produces a value in `[0,m)` congruent
to `x`, hence exactly `x mod m`.

Lean checks this full general natural-number statement, with actual floor
division definitions rather than an assumed quotient-error conclusion, as
`barrett_residual_bounds` and `barrett_one_correction`. For worker inputs
`a,b<m<=2^64-1`, `x=a*b<B`. Also `q*m<=x<B`, the residual is below
`2m<2^65`, and the corrected result is below `m`. Consequently its
multiplication, subtraction and final 64-bit cast are exact under the
canonical-input invariant. The final bridge from this natural-number
theorem to the Rust word operations is manual plus native verification.

## Proposition 4: canonical addition and overflow-aware subtraction

For `a,b<m`, the integer sum `S=a+b` lies in `[0,2m)`. One subtraction
when `S>=m` therefore returns `S mod m`. Lean checks this as
`canonical_add_one_subtraction`.

For 64-bit words, put `W=2^64` and assume `m<W`. If `S<W`, the machine
sum equals `S` and the conditional test selects the same branch. If
`S>=W`, the overflowing sum is `S-W`; `S>=m` and `0<=S-m<m<W`.
Wrapping subtraction of `m` from `S-W` therefore produces exactly `S-m`.
This establishes the worker's `(carry || sum>=m)` implementation under
the invariant, including moduli above `2^63`. The word-model bridge is
manual plus native tests, not a BitVec proof of the Rust operation.

## Canonical-input and randomness audit

All `StrongRhoG` hot-path scalar inputs are canonical: samples are in
`[0,m)` or `[1,m)`, lambda powers and collision values use exact modular
arithmetic, and each canonicalization or state update maps canonical
residues to canonical residues. `mul_mod`, `sub_mod` and `inverse_mod`
delegate to the exact production `u64` code. The optimized `mul_in` and
`add_mod` are intended for those canonical inputs; no unrestricted input
equivalence is claimed.

`ExactUniform` delegates constructed sampling and both one-shot sampling
methods to `UniformInt<u64>`, preserving rejection draws and RNG positions.
The worker's tests compare exclusive/inclusive and constructed/one-shot
streams across four ranges, including ranges above `2^63`, and compare
the next RNG word after each cycle. The exact `12345*G` solve checks
answers, jumps, selected canonical states, exported counters, charges,
and the post-solve RNG position. No full per-step trace is exported.

## Independent native checks

[`cpp_check.json`](cpp_check.json) records 100,049 exact high-product
checks and 1,048,575 multiply/add checks across 1,023 moduli, including
2, subgroup-order neighbors, `2^50` neighbors, `2^63` neighbors and
`2^64-1`. The high-product and reciprocal oracles use GMP arbitrary
precision; remainders use independent exact unsigned-128 division.
41,936 tested products require the correction branch. This finite test
suite supports implementation correctness; it is not a universal theorem
or an independent replay of solver records.

The first C++ build failed because Boost headers were absent. That failed
build is retained, and the implementation uses the installed GMP backend.
All source, commands, receipts and formal verification status are retained
beside this report. The final benchmark report links the exact campaign
replay receipt separately.

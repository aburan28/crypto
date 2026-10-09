# Wide-order ISO-1 assessment (2026-10-09)

This note separates three inputs that have similar subgroup bit lengths but
different fields. The implementation changes in this branch widen sparse
linear algebra and parts of the `F_{p^6}` tower. They do **not** make the
complete ISO-1 instance generator or sieve accept a wide subgroup order.

## Scope and evidence

The existing ISO-1 path uses `p: u64`, `Spec.l: u64`, `Spec.d: u64`,
`SparseRel` coefficients and rhs as `u64`, and a `u64` Wiedemann solver.
`group_order_4_prime` uses a `u128` group order but rejects `l > u64::MAX`.
The end-to-end walk and sieve reports and descent coefficients also use
`u64`. Those are active interface limits, not just serialization limits.

This branch adds `WideSparseRel` and `wiedemann_big` with arbitrary-precision
modular multiplication and inversion, plus `Fq3::pow_big`, a wide field-order
path for square testing and Tonelli–Shanks, and `mul_big` for both elliptic
points and Jacobian divisors. The old `u64` Wiedemann path remains available.
The tower's noncube search now skips the `F_p` slice when `p ≡ 2 (mod 3)`:
every nonzero element of that slice is a cube, so scanning it costs `p`
fruitless exponentiations. Candidate order for existing instances is
preserved at the first viable element.
Tests exercise exact linear systems modulo a 113-bit prime and the P-192
subgroup order; a tower test uses base primes 524309 and 4820937851.
These tests check components, not relation collection or logarithm recovery.

## Extension-field curve with a roughly 113-bit subgroup

For this route, the elliptic curve lives over `F_{p^6}`, with `p` around
`2^(115/6)` because the weak curve has cofactor four. One representative
base prime is 524309, for which `p^6` has 115 bits and `p^6/4` has 113
bits. The tower fits its present `u64` coefficients, and a full
group order fits `u128`. A prime subgroup order above 64 bits needs wide
primality testing, wide instance generation, wide descent coefficients,
and a wide sparse-system call wired into the actual sieve. None of those
four steps is complete here.

The point-count search in `group_order_4_prime` is a baby-step giant-step
search in a Hasse interval of width about `4 p^3`. Its table has about
`2 p^(3/2)` points, roughly 7.6e8 entries at p=524309; that alone is
beyond a routine local test. An externally certified group order could
avoid that search, but subgroup and transfer checks would still be needed.
The factor base has order `p/2`, about 2.6e5 classes. A bounded component
test is feasible; a complete measured attack requires additional resources
and wide-order integration. In the present sequential Wiedemann code, the
Krylov sequence alone needs `2N` matrix-vector products. For a roughly
six-entry row and `N ≈ p/2`, that is on the order of `12N²`, roughly
8e11 wide modular products before the solution reconstruction and
verification. This is a code-derived workload estimate, not a timing.

## Extension-field curve with a roughly 192-bit subgroup

Here `p` is around `2^(194/6)` and `p^6` around `2^194`. The prime
4820937851 exercises the wide field-order path, but it is only a field
fixture, not a certified weak elliptic curve with a prime-order subgroup.
Full order and trace arithmetic exceed `u128`, and the current
`group_order_4_prime` and isogeny-walk trace code must be redesigned.
Its Hasse-interval baby-step giant-step table is on the order of
`2 p^(3/2)`, about 6.7e14 entries. The factor base alone is on the order
of 2.4 billion classes. Arbitrary-precision coefficients resolve an
arithmetic limit but do not make that workload practical.

## NIST P-192 itself

P-192 is a curve over a 192-bit **prime field** `F_P`, with prime order
and cofactor one, not a curve over `F_{p^6}` for a 32-bit `p`.
Its subgroup order can be used as a modulus in the wide linear algebra
test. That test gives no map from P-192 into the ISO-1 cover setting.
An isogeny defined over `F_P` preserves the number of `F_P`-rational
points; hence it cannot take a prime-order P-192 curve to one with full
`F_P`-rational 2-torsion. Extending scalars to `F_{P^6}` changes that
particular obstruction but makes the ISO-1 factor-base parameter `p=P`,
about `2^192`. The factor base and relation work would then be vastly
larger than the generic Pollard-rho scale `sqrt(n)`, about `2^96` group
operations. This is an assessment of this specific construction; it does
not rule out other ECDLP methods.

## Standardized 113-bit binary curves

The SEC 2 `sect113r1` and `sect113r2` curves live over `F_{2^113}`.
Their field characteristic, curve equation, and subgroup structure differ
from this odd-prime `F_{p^6}` implementation. A 113-bit scalar modulus
working in `wiedemann_big` does not provide a transfer or decomposition
algorithm for either binary curve. They would require a separate
characteristic-two construction and cost model.

Sources: [NIST SP 800-186](https://nvlpubs.nist.gov/nistpubs/SpecialPublications/NIST.SP.800-186.pdf)
(P-192 is a prime-field Weierstrass curve, legacy use, and the pseudorandom
prime-field curves have cofactor one),
[FIPS 186-4](https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-4.pdf)
(archived P-192 parameters), and
[Joux–Vitse](https://eprint.iacr.org/2011/020.pdf)
(cover and decomposition over composite-degree extension fields), and
[SEC 2](https://www.secg.org/SEC2-Ver-1.0.pdf)
(113-bit binary curves).

## Acceptance state

| Requirement | State |
| --- | --- |
| Wide modular sparse solve at 113 and 192 bits | Verified: `cargo test --lib wiedemann_wide_prime_orders` passed, 1 test, on this branch |
| `F_{p^6}` exponentiation and square root with `p^6 > u128` | Verified: `cargo test --lib tower_field_and_curve_scalars_beyond_u64` passed at both listed base primes |
| Big scalar multiplication on existing tower/Jacobian types | Verified: above test and `cargo test --lib jacobian_scalar_beyond_u128` passed |
| 113-bit subgroup instance generation and complete ISO-1 solve | Unresolved: point count, primality, `Spec` and sieve/descent widths |
| 192-bit extension-field instance and complete ISO-1 solve | Unresolved: above plus full-order and walk trace widths, scale |
| Actual P-192 ISO-1 attack | Unresolved: no applicable cover/transfer route established |
| `sect113r1`/`sect113r2` ISO-1 attack | Unresolved: odd-prime ISO-1 does not apply; the binary GHS route is separately screened and only its trace branch is executable |

The next implementation gate is a **certified wide-order input** for a
specific weak curve, with independent checks of `[l]G = O`, `Q = [d]G`,
the transferred points, and the subgroup order. Only then can the sieve
be widened and run on a real 113-bit instance without relying on the
infeasible order-search path.

Validation on this branch used the unoptimized Cargo test profile. The
targeted tests `wiedemann_wide_prime_orders`,
`tower_field_and_curve_scalars_beyond_u64`, and
`jacobian_scalar_beyond_u128` each passed (one test each). The existing
`tower_arithmetic` regression test and both tests selected by
`wiedemann_solves_random_sparse_systems` also passed. The initial sparse-checkout
test build failed because embedded fixtures were absent; restoring the
tracked fixtures resolved that setup failure. A first wide tower test was
interrupted after it exposed the noncube-search scan; the corrected search
was then tested successfully. The complete repository test suite was not
run for this change.

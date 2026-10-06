# P-256 exact-endomorphism and scalar-orbit factor-base screen, round 25: result

Date run: 2026-10-06

The screen found an exact one-log-class factor base, but it also proves why
that does not improve P-256 index calculus.  The closest attainable scalar
orbit has 139,592 negation-folded columns and every one of its 139,592
transport edges replays in the P-256 group.  With a known anchor it is a
generic rho/claw encoding; with an unknown anchor it needs at least two target
relations and is already 1.447 times rho at the free-oracle boundary.  Its
coordinate membership degree is 139,592, not a residual degree at most five.

No full-depth unplanted relation was attempted.

## Exact endomorphism-ring result

The previously unresolved 216-bit CM-discriminant cofactor splits completely.
The binary verifies the product and supplies recursive Lucas certificates; no
external factor database is trusted by the result:

```text
|Delta| = 455213823400003756884736869668539463648899917731097708475249543966132856781915
        = 3 * 5 * 456597257999
          * 1428624589419343516204097
          * 46523541035814968339936406074986559003387.
```

All five factors are distinct primes.  Since `Delta=-|Delta|` is congruent to
1 modulo 4 and squarefree, it is a fundamental discriminant.  Consequently

```text
Z[pi] = End(E) = O_K
```

and the Frobenius conductor is exactly one.  The minimum degree of a
non-integer endomorphism follows from

```text
N((a+b*sqrt(Delta))/2) = (a^2 + |Delta|*b^2)/4:

min degree = (1+|Delta|)/4
           = 113803455850000939221184217417134865912224979432774427118812385991533214195479
           = 2^255.9750076747.
```

More importantly for the factor-base question, every endomorphism is
`a+b*pi`, and `pi(P)=P` for every `P` in `E(F_p)`.  Its rational-point action
is therefore the known scalar `[a+b]`.  P-256 has no low-degree non-scalar
endomorphism and no non-generic rational-point transport hidden in its CM
ring.

## Exhaustive scalar-width screen

The complete factorisation of `n-1` was verified, its 91-bit factor and `n`
were proved prime by replayable Lucas certificates, and the least primitive
root modulo `n` is 7.  Exhaustive enumeration covered all 10,240 divisors.

| scalar order | folded columns | distance from 131,458 |
|--:|--:|--:|
| 279,184 | **139,592** | **8,134** |
| 146,589 | 146,589 | 15,131 |
| 293,178 | 146,589 | 15,131 |
| 114,567 | 114,567 | 16,891 |
| 229,134 | 114,567 | 16,891 |

The selected order is `279184 = 16*17449`; its halfway power is `-1`, so the
full orbit folds to 139,592 columns.  All 139,584 exact-order generators of
that subgroup were screened.  The minimum signed binary-cost representative
is

```text
m = 379483656589291393959605088584129393463128953089898693989044430844151895048
```

with 248 bits, Hamming weight 96, and 344 binary double/add operations per
application.

## Native FB1 materialisation and replay

The SHA-256 try-and-increment anchor lifted at counter 1.  Its discrete log
relative to `G` is not part of the construction.  The full orbit was emitted
using the repository's native `ecbench.factor_base/v1` identity and
`ecbench.factor_base_dump/v1-wide` row format:

```text
FB1he6b6b6e25de6
columns       139,592
signed points 279,184
points SHA-256 7b2df8a9acc38cd338847abb215c007bd03a28243244451acb36f93131e3fde3
dump bytes     70,236,610
dump SHA-256   1613190ec5cc3e0bd83c1c51984334cae2dababf116bb1480f23920f7080cb4d
```

All points are nonidentity, on curve, and unique up to sign.  Every adjacent
equation `P_(i+1)=[m]P_i` and the final equation `[m]P_(B-1)=-P_0` was replayed
directly in the P-256 group: 139,592 checks, zero failures.  The run charged
13,400,832 replay additions, 34,618,816 replay doublings, 4,449,455 fixed-base
output additions, 8,160 table additions, 256 table doublings, 697,942
normalisation multiplications, and 18 batch inversions.

The exact log quotient is one dimensional.  Ordinary zero relations do not
determine its anchor: after substituting `log(P_i)=m^i h`, every such row is
the already-known coefficient identity times `h`.

## End-to-end boundary

Using round 24's exact negation-folded collision law and the selected width:

| case | target collisions | oracle S | / rho | log2 oracle work | log2 direct work | direct / rho | log2 one-list bytes |
|:--|--:|--:|--:|--:|--:|--:|--:|
| known anchor `H=G` | 1 | 1.253637 | 0.964336 | 128.326 | 131.326 | 7.714692 | 133.326 |
| unknown hash anchor | 2 | 1.880456 | 1.446505 | 128.911 | 131.911 | 11.572038 | 133.911 |

The known-anchor oracle constant is below the repository's `1.3*sqrt(n)` rho
reference, but all factor-base logs are then already known and the target
representation is a generic claw search.  Direct S17 sample construction is
7.715 times rho and materialisation needs about `2^133.326` bytes.  With the
unknown anchor, two independent target equations are necessary to eliminate
`h`; even an unimplemented free oracle is above rho before sample generation.

The exact support polynomial has degree 139,592.  The rational multiplication
map `[m]` has degree `m^2`, approximately `2^495.493`.  A short scalar
evaluation chain is not a low-degree membership equation.  Structured
residual degree remains unknown, so the degree-at-most-five gate fails.

## Decision

This is a positive transport result and a negative attack result.  Exact
elliptic-log transport can reduce a P-256 factor base to one class, but only by
choosing a scalar orbit.  That construction either exposes all logs and is rho
under another name, or retains an anchor whose elimination costs more than
rho.  It does not lower the S17 summation-polynomial degree, direct work,
memory, or complete one-target cost.  Every promotion gate except correctness
therefore fails.

The remaining search target is narrower: a factor base with low-degree
coordinate membership and exact log transport whose point action is not a
known scalar.  The complete CM-ring result rules that out among P-256
endomorphisms over the base field.

## Reproduction

```bash
cargo test --bin p256_scalar_orbit_factor_base
cargo run --release --bin p256_scalar_orbit_factor_base -- \
  --round24 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_negation_cross_colour_round24_20261006/negation-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_scalar_orbit_round25_20261006/scalar-orbit-result.json \
  --dump target/p256-scalar-orbit-round25.factor-base.json
```

Canonical result artifact: `scalar-orbit-result.json`, 14,033 bytes, SHA-256
`dcb191f7d5d7a660e3ce13128b25ce1c017c7d89ec2d01a242ab536651647faf`.

# Frozen selection protocol: isogeny-class trio (2026-10-17)

Written before the full B- and P-shape searches ran. The Koblitz table is an
exhaustive finite enumeration (24 degrees x 2 values of `a`), so it has no
search freedom; it was computed before this file was written (disclosed
here), and the K rule below is stated against it. Two 400-seed calibration
runs (`binary --m 47 --count 400`, `prime --bits 47 --count 400`) measured
timing before this file was written; the full runs below include those seeds
and supersede them.

## Request (user, 2026-10-17)

Pick three curves, one binary, one Koblitz and one prime-field, between 41
and 66 bits (preferably above 45), with shapes that emulate standardized
curves, whose isogeny classes have Frobenius conductor `f > 1` so that the
volcanoes have more than one level, as roots for an experiment that iterates
over every curve on the isogeny graph.

## Definitions

- `t` is the Frobenius trace, `Delta = t^2 - 4q = f^2 d_K` with `d_K`
  fundamental. `f` is the conductor of `Z[pi]` in `O_K`.
- For each prime `l | f` the `l`-volcano has height `v_l(f)`; its crater is
  split, ramified or inert by the Kronecker symbol `(d_K / l)`.
- The strata of the class are the orders `O_c`, `c | f`; stratum `c` holds
  `h(c^2 d_K)` curves and the class holds `H = sum_{c | f} h(c^2 d_K)`
  F_q-isomorphism classes (Schoof 1987, Thm 4.6).
- "Walkable" means `f` is 61-smooth: every vertical step has degree at most
  61, the practical limit of the repository's modular-polynomial walker
  (`docs/isogeny-walk/README.md`).

## Shapes

| Family | Standard it emulates | Model | Field | Cofactor |
|:--|:--|:--|:--|:--|
| B (binary) | FIPS 186 B-163 ... B-571 | `y^2 + xy = x^3 + x^2 + b`, `b` from a public seed | `GF(2^m)`, `m` prime, lowest-weight trinomial else pentanomial (FIPS rule) | 2 |
| K (Koblitz) | FIPS 186 K-163 ... K-571 | `y^2 + xy = x^3 + a x^2 + 1`, `a in {0,1}` | same | `#E(F_2)`: 4 for `a = 0`, 2 for `a = 1` |
| P (prime) | FIPS 186 P-192 ... P-521 | `y^2 = x^3 - 3x + b`, `b` from a public seed | `GF(p)`, `p = 2^k - c` the largest prime below `2^k` | 1 |

## K rule (exhaustive table, `selector koblitz --min 41 --max 64`)

1. Only prime `m` (the standard family).
2. Primary pick: the standard cofactor (`#E = 4r` for `a = 0`, `2r` for
   `a = 1`, `r` prime) with `m >= 46` if one exists; otherwise the
   standard-cofactor curve in `[41, 64]` with the largest `r`.
3. Alternate pick for "everything above 45 bits": the prime-`m` curve with
   the largest prime subgroup `r` among those with at least two distinct
   primes in `f` (at least four strata).

## B and P rules

- Primary size: `m = 47` and `k = 47`, the smallest prime degree above 45,
  which keeps the full class near `10^7` curves.
- Alternate size: `m = 61` and `k = 61` (`p = 2^61 - 1`, the Mersenne prime,
  as in P-521 = `2^521 - 1`), for a trio that is entirely above 45 bits.
- Seeds: `b = SHA-256("isogeny-class-trio/2026-10-17/B-<m>/b/<i>")`, low `m`
  bits; for P, the low 64 bits mod `p` of
  `SHA-256("isogeny-class-trio/2026-10-17/P-<k>/b/<i>")`. Indices
  `i in [0, 131072)` at 47 bits and `i in [0, 32768)` at 61 bits.
- Order: BSGS over the Hasse window on seeded points (with `#E = 2 mod 4`
  for B). A curve is kept only if the counted order is the unique match and
  `r` is prime; `r > 4 sqrt(q)` then proves `#E` exactly.
- Qualifying: `f > 1`, `f` 61-smooth, `max_l v_l(f) >= 2`.
- Ranking: largest `max_l v_l(f)`, then most strata (divisors of `f`), then
  lowest index `i`.

## Recorded for each pick

Curve, field, generator from a public seed and its order check, `t`,
`Delta`, `d_K`, `f`, every volcano (height, crater type, degree of the field
holding the `l`-kernel x-coordinates), exact `h(d_K)` by counting reduced
forms, `h` of every stratum, `H`, and the root's level where it is provable:
K by `tau in End(E)` (crater), B for `l = 3` by `psi_3` root counts and
characteristic-2 Velu walks, P by the repository walker
`src/bin/isogeny_walk.rs class` as an independent cross-check.

## Not claimed

No ECDLP-cost, weakness or security statement about any curve. The selection
is exact arithmetic (orders, factorizations, class numbers); there is no
empirical or comparative claim, so no control population is needed.

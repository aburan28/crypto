# Curve testbed: recommended curves for testing

These curves are the recommended starting points for testing isogeny-class,
volcano and ECDLP tooling. Every one has a standardized shape (FIPS 186
K-, B- or P-curve form), a proved group order and generator, and a
Frobenius conductor `f > 1`, so its isogeny class has volcanoes with more
than one level. Each class has been counted exactly. The full curve catalog
is the registry, [`../registry.json`](../registry.json); this page only
highlights a few curves from it. Machine-readable list:
[`testbed.json`](testbed.json).

## ★ Featured: the 61-bit trio

| | Koblitz | Binary | Prime |
|:--|:--|:--|:--|
| Curve | `icv1-f2m61-t158598901-ab42b6c5` | `icv1-f2m61-tm2657393969-b9ff0fb0` | `icv1-fp61-tm813574577-a3b71255` |
| EC1 | `EC1N61Ce0hd8c15448081e` | `EC1N61Crbh854aa3d90d2d` | `EC1P61Cfph303b4b060ff9` |
| Model | `y^2 + xy = x^3 + 1` over `GF(2^61)`, `x^61 + x^5 + x^2 + x + 1` | `y^2 + xy = x^3 + x^2 + 0x9c60f836d06fc8`, same field | `y^2 = x^3 − 3x + 0x54aa66b458b77ec` over `p = 2^61 − 1` |
| Prime subgroup | 48 bits (cofactor 14156) | 61 bits (cofactor 2) | 62 bits (prime order) |
| Conductor `f` | `1951 · 587551` | `7^2` | `3^2 · 5^3` |
| Volcanoes | heights 1, 1 at large primes (inert) | height 2 (split crater) | heights 2 (inert) and 3 (split) |
| Strata / class size | 4 / 1,147,491,009 | 3 / 851,960,648 | 12 / 920,783,750 |
| Where it sits | crater; the only curve with `End = Z[τ]` | floor (`End = Z[π]`) | floor of both volcanoes |
| What it is good for | Frobenius (τ) endomorphism; large prime conductor gaps; the class is reachable only through degree-1951 and degree-587551 steps | char-2 7-isogeny volcanoes; a near-prime order like B-163 | the P-521 prime form; many strata reachable with 3- and 5-isogenies (the repository walker supports both) |

## ★ Featured: large prime conductor gap (61 bits)

These curves have Frobenius conductor `f = l`, a single large prime, so the
class is a crater (`End = O_K`) and a floor (`End = Z[π]`) at index `l`.
They were built by class-number-one CM, and each has prime order.

| Curve | EC1 | Form | Gap `l` | Class size | Good for |
|:--|:--|:--|:--|:--|:--|
| `icv1-fp61-t205-983ecb51` | `EC1P61Cfph54bb8886c801` | `y^2 = x^3 + 5` over `p = 0x1ffffff43f9a9a39` (secp256k1 form, `j = 0`, order-6 automorphisms on the crater) | 1,753,413,037 | 584,471,013 | crater-vs-floor gaps with the `sqrt(3)` automorphism mechanism of `research/drastic_prime_conductor_gap_round2_20261012/` |
| `icv1-fp61-t7175-7e38f752` | `EC1P61Cfphde4753195a5f` | `y^2 = x^3 − 3x + 0x14eb628f40fdf1f3` over `p = 2305842959800795691` (P-curve form, `d = −19`) | 696,735,691 | 696,735,693 | the largest prime gap with an `a = −3` model |
| `icv1-f2m61-t158598901-ab42b6c5` | `EC1N61Ce0hd8c15448081e` | the Koblitz curve above | 1951 and 587551 | 1,147,491,009 | binary: for prime `m`, class-number-one CM gives only Koblitz classes |

Also built: `icv1-fp61-t727-ddf032b6` (`d = −43`, `l = 463,138,201`),
`icv1-fp61-t2279-df6b5e41` (`d = −67`, `l = 371,028,661`),
`icv1-fp61-t95-4d1bb712` (`d = −163`, `l = 237,876,229`). No prime-order
curve exists for `d = −4, −7, −8, −11`; the report gives the reasons.

Caution: from the crater, the floor is reachable only by a degree-`l`
isogeny, which no tool here computes. The floor curves have to be
constructed directly.

## Also useful: the 41/47-bit trio (smaller classes)

| Curve | Shape | `f` | Strata / class size | Root |
|:--|:--|:--|:--|:--|
| `icv1-f2m41-tm2308219-7f48b14a` | Koblitz, cofactor 4 (the only standard-cofactor prime-degree Koblitz curve from 41 to 66 bits) | `409 · 1721` | 4 / 708,153 | crater |
| `icv1-f2m47-t4276051-826c5755` | binary, cofactor 2 | `31^2` | 3 / 11,754,700 | floor |
| `icv1-fp47-tm8785349-43e36ce5` | prime, `p = 2^47 − 115`, prime order | `3^5` | 6 / 13,461,228 | floor |

## Notes

- Binary curves over odd-degree fields only have volcanoes at primes
  `l ≡ ±1 (mod 8)` (7, 17, 23, 31, …), never at 3 or 5.
- For a prime-degree Koblitz curve every prime dividing `f` is `≡ ±1 (mod m)`.
- Selection, proofs, figures and the PDF are in
  [`research/isogeny_class_trio_20261017/`](../../../research/isogeny_class_trio_20261017/README.md).

Status: the parameters are exact and verified (native selector tests,
independent walker cross-check for the prime curves). Nothing on this page
claims an ECDLP cost, weakness or security property.

## Adding a curve

To add a curve, give it a registry entry and evidence that it has the
properties its row claims. Mark it `featured` in `testbed.json` only when
its owner recommends it.

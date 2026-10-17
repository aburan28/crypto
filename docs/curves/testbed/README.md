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

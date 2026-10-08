# Curve traits: grouping registered curves by structure

[`registry.json`](registry.json) says *which* curve a name denotes. It does
not say which curves are *alike*: that every Koblitz curve shares one CM
field and one descended Frobenius whatever its degree, that a curve over
`GF(2^18)` is the quadratic twist of one, or how far `Z[π]` sits below the
maximal order. [`traits.json`](traits.json) records those traits for every
registry curve, each value with how well it is established, and
`curve_traits` groups curves by them across sizes or ranks curves by
likeness to one.

```bash
cargo build --release --bin curve_traits
./target/release/curve_traits keys                          # every grouping key, explained
./target/release/curve_traits show ECC2K-130                # one record (slug, standard name or legacy spelling)
./target/release/curve_traits group --by cm,descent --min-sizes 2
./target/release/curve_traits similar ECC2K-130 --other-sizes --top 10
./target/release/curve_traits build                         # registry.json → traits.json
./target/release/curve_traits check                         # CI: traits.json is current
```

Every query reads the committed `traits.json`; nothing is recomputed.
`--json` gives machine-readable output.

## What each record holds

All of it is derived from the registry's model and recorded order
`#E = q + 1 − t`. The implementation is
[`src/cryptanalysis/curve_traits/`](../../src/cryptanalysis/curve_traits/mod.rs).

| Field | Meaning | How |
|:--|:--|:--|
| `order_check` | how well `#E` itself is established | exhaustive count for `q ≤ 2^22`; for a binary curve defined over `GF(2^k)`, `k ≤ 16`, a count over `GF(2^k)` lifted by `t_n = V_{n/k}(t_k, 2^k)`; otherwise a generator certificate (`G` on the curve, `[r]G = O`, `r` prime and `r > 4√q`, so one multiple of `r` lies in the Hasse interval) |
| `trace_ratio` | `t / 2√q ∈ [−1, 1]` | size-free |
| `frobenius.disc` | `Δ = t² − 4q`, its primes and any unfactored composites | trial division to 2¹⁶, perfect powers, Brent's rho with a fixed iteration budget |
| `frobenius.cm_disc` | `d_K`, the fundamental discriminant of `Q(π)` | `Δ = v²·d_K` |
| `frobenius.conductor` | `v = [O_K : Z[π]]`: the gap between `Z[π]` and the maximal order, with its primes | the square part of `Δ` |
| `frobenius.conductor_fraction` | `log v / log √|Δ|`: 0 when `Z[π]` is maximal, near 1 when `d_K` is tiny beside `Δ` | size-free |
| `frobenius.class_number` | `h(d_K)` for `|d_K| ≤ 10^7` | reduced forms |
| `small_primes` | for each `ℓ ≤ 31`: `v_ℓ(v)` (the depth of the `ℓ`-isogeny volcano), how `ℓ` splits in `Q(π)`, and the number of rational `ℓ`-isogenies when the depth is 0 | exact from `Δ` modulo powers of `ℓ`, without factoring `Δ` |
| `subfield` (binary) | `j_field_degree`: least `k` with `j ∈ GF(2^k)`; `definition_degree`: least `k` such that the curve is the base change of one over `GF(2^k)`; `base_trace`: that curve's trace `t_k` | see below |
| `subgroup` | prime subgroup order `r` and cofactor | the registry representation, else the largest prime of `#E` |
| `embedding` | `ord_r(q)` and the bit length of `(r − 1)/k` | exact from the factorisation of `r − 1`; else `q^k ≢ 1` for `k ≤ 1000` gives a bound |
| `twist` | the quadratic twist's order `q + 1 + t`, its largest prime and cofactor | factored with the same budget |
| `keys` | the grouping keys below | |

**Subfields.** `E_{a,b} : y² + xy = x³ + ax² + b` over `GF(2^n)` is the base
change of a curve over `GF(2^k)` exactly when `b ∈ GF(2^k)` and either `n/k`
is odd or `Tr(a) = 0`. So `j` can lie in a subfield the curve does not
descend to: the quadratic twist of a Koblitz curve over an even degree has
`j = 1` but descends only to `GF(4)`. When the curve descends, `t_k` is
found as the unique Hasse-bounded root of `V_{n/k}(t_k, 2^k) = t_n`, or by
counting points over `GF(2^k)` when the root is not unique. When `n/k` is
even, both twists over `GF(2^k)` descend and only `|t_k|` is an invariant
(`base_trace_sign_free`).

**What is not decided.** `End(E)` lies between `Z[π]` and `O_K`. Which
order it is, and so each curve's level in its volcanoes, needs a walk
(`src/bin/isogeny_walk.rs`); the record gives the depths, not the levels.
Except at depth 0, the rational `ℓ`-isogeny count depends on that level and
is left out. Two registry models of one curve under different moduli are
two records, as ICV1 requires; isomorphism across moduli is not detected.

## Status of a value

| Status | Meaning |
|:--|:--|
| `proved` | exact; every primality claim under it is deterministic (Miller–Rabin to the bases 2…41 is a proof below 3.3·10²⁴) |
| `probable` | exact if every probable prime under it is prime |
| `bounded` | only a bound is known: `conductor` is then a lower bound and `cm_disc` is absent |
| `unknown` | attempted and not determined within the factoring budget |
| `not_evaluated` | skipped by a declared limit (`h(d_K)` above `|d_K| = 10^7`, a subfield above degree 40) |
| `not_applicable` | does not apply (no proper subfield of definition, a prime field's subfields) |

An unfactored composite of even multiplicity still fixes `d_K`. Every
Koblitz curve has `Δ = −7·v²`, so its CM field and conductor are exact even
where `v` itself does not factor (the degree-233, -409 and -577 Koblitz
curves). The budget counts rho iterations, not seconds, so a rebuild is
byte-identical on any machine. A registry entry the model contradicts
(an order that the count refutes, a generator `[r]G ≠ O`, a certified
`End(E)` outside `Q(π)`) stops the build; it is an error to fix, not a
trait.

## Grouping keys

Each key is discrete and does not grow with the field, so equal keys can
join curves of any size. `curve_traits keys` prints this list.

| Key | Values |
|:--|:--|
| `char` | `2` or `p` |
| `family` | the registry's family label |
| `ordinary` | `yes` when `p ∤ t` |
| `cm` | `d_K` when `|d_K| ≤ 10^6`; `large` above; `unknown` when `Δ` did not factor |
| `descent` | `k<k>\|t\|=<|t_k|>` for the smallest subfield of definition `GF(2^k)` (Koblitz: `k1\|t\|=1`); `none` |
| `descent_signed` | `descent` with the sign of `t_k`, or `±` when both twists descend |
| `jfield` | `k<k>` when `j` lies in a proper subfield `GF(2^k)`; `full` |
| `cofactor` | `#E / r` |
| `twist_cofactor` | the twist's order over its largest prime |
| `embedding` | `k=<k>` for an embedding degree at most 20; `large`; `unknown` |
| `conductor` | `1` when `Z[π]` is maximal; `smooth` when every prime of `v` is below 2¹⁶; `rough`; `unknown` |
| `split` | how 3, 5, 7, 11, 13 split in `Q(π)`: `s`, `i` or `r` for split, inert, ramified |
| `depth` | `v_ℓ(v)` for `ℓ = 2, 3, 5, 7` |

A curve's **signature** is its values of `char, ordinary, cm, descent,
cofactor`, the default of `group --by`.

## Similarity

`similar <curve>` ranks every other curve by the summed weight of the keys
on which it differs from the target: `cm` and `descent` 4, `char` and
`ordinary` 3, `jfield` and `cofactor` 2, the rest 1. Ties are broken by
`|Δ trace_ratio| + |Δ conductor_fraction|`. A key that is `unknown` on
either side counts as a difference and is shown with a `?`: two curves
whose discriminants did not factor are not thereby alike. `--keys` compares
on a chosen set instead, each weighted 1; `--other-sizes` leaves out curves
over a field of the target's size.

## What the current registry shows

From `group --by cm,descent --min-sizes 2` on the 121 registry curves:

- **CM by −7, defined over `GF(2)`: 59 curves at 41 degrees**, 7 to 577:
  the `K_a / GF(2^n)` family, every member with `|t₁| = 1`;
  `descent_signed` separates `K_0` (`t₁ = −1`) from `K_1` (`t₁ = 1`).
  Adding `cofactor` picks out the `K_0` curves of order 4·prime (14 curves
  from `GF(2^7)` to `GF(2^571)`, ECC2K-130 among them) and the `K_1`
  curves of order 2·prime (8 curves up to sect163k1); 35 others have a
  composite `#E/4` or `#E/2`, and two (degrees 157 and 577) an order that
  did not factor within the budget.
- **CM by −15, defined over `GF(4)` with `t₂ = 1`: 3 curves**, over
  `GF(2^6)`, `GF(2^10)` and `GF(2^14)`.
- `icv1-f2m18-t999-40751283` has CM by −7 and `j = 1` but descends only
  to `GF(4)` (`t₂ = 3`): it is the quadratic twist of the Koblitz curve
  over `GF(2^18)`. It shares the CM key with the Koblitz family and differs
  in `descent`.
- secp256k1 has `d_K = −3` (`j = 0`) and a 128-bit conductor;
  P-256's `Δ` does not factor within the budget, so its `cm` is `unknown`
  and its conductor only bounded.
- Every order in the registry is certified: 60 by exhaustive count, 41 by
  a descent count, the rest by generator certificates (4 of those rest on
  a probable prime).

## Keeping it current

`traits.json` is derived from `registry.json`. When the registry changes,
rebuild it in the same pull request:

```bash
cargo run --release --bin curve_traits -- build
```

The `curve-traits` workflow runs the module's tests and `curve_traits
check`, which rebuilds with the committed file's budget and fails on any
difference.

Not yet covered: curves the isogeny walker writes (`isogeny_walk` output
directories) and the IC crosswalk's `ic/curves.yaml`; the DiSSECT trait
slots in crypto-autoresearcher's curve catalog. Each needs an adapter that
keeps its own identity rules.

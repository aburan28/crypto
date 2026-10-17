# Curve traits: grouping registered curves by structure

[`registry.json`](registry.json) says *which* curve a name denotes. It does
not say which curves are *alike*: that every Koblitz curve shares one CM
field and one descended Frobenius whatever its degree, that a curve over
`GF(2^18)` is the quadratic twist of one, or how far `Z[π]` sits below the
maximal order. [`traits.json`](traits.json) records those traits for every
registry curve and for the curves of the other sources in
[`trait-sources.json`](trait-sources.json), each value with how well it is
established, and
`curve_traits` groups curves by them across sizes or ranks curves by
likeness to one.

```bash
cargo build --release --bin curve_traits
./target/release/curve_traits keys                          # every grouping key, explained
./target/release/curve_traits show ECC2K-130                # one record (slug, standard name or legacy spelling)
./target/release/curve_traits group --by cm,descent --min-sizes 2
./target/release/curve_traits similar ECC2K-130 --other-sizes --top 10
./target/release/curve_traits build                         # registry + sources + hints → traits.json
./target/release/curve_traits check                         # CI: traits.json is current
./target/release/curve_traits unfactored                    # composites left unsplit (input of factor-hints.gp)
```

Every query reads the committed `traits.json`; nothing is recomputed.
`--json` gives machine-readable output.

## What each record holds

All of it is derived from the registry's model and recorded order
`#E = q + 1 − t`; reading the registry checks that identity and the Hasse
bound `|t| ≤ 2√q` before any trace is classified. The implementation is
[`src/cryptanalysis/curve_traits/`](../../src/cryptanalysis/curve_traits/mod.rs).

| Field | Meaning | How |
|:--|:--|:--|
| `order_check` | how well `#E` itself is established | exhaustive count for supported models with `q ≤ 2^22`; for a binary curve defined over `GF(2^k)`, `k ≤ 16`, a count over `GF(2^k)` lifted by `t_n = V_{n/k}(t_k, 2^k)`; otherwise a generator certificate when the model's group law is implemented (`G` on the curve, `[r]G = O`, `r` prime and `r > 4√q`, so one multiple of `r` lies in the Hasse interval); over `GF(p^k)` with no usable generator, a point-order certificate (see below) |
| `trace_ratio` | `t / 2√q ∈ [−1, 1]` | size-free |
| `frobenius.disc` | `Δ = t² − 4q`, its primes and any unfactored composites | trial division to 2¹⁶, perfect powers, Brent's rho with a fixed iteration budget |
| `frobenius.cm_disc` | `d_K`, the fundamental discriminant of `Q(π)` | `Δ = v²·d_K` |
| `frobenius.conductor` | `v = [O_K : Z[π]]`: the gap between `Z[π]` and the maximal order, with its primes | the square part of `Δ` |
| `frobenius.conductor_fraction` | `log v / log √|Δ|`: 0 when `Z[π]` is maximal, near 1 when `d_K` is tiny beside `Δ` | size-free |
| `frobenius.class_number` | `h(d_K)` for `|d_K| ≤ 10^7` | reduced forms |
| `source` | `registry`, `isogeny-walk` or `ic-crosswalk`; the file; the source's `curve_uid`; whether the slug is registered | see "Sources" |
| `endomorphism` | what is known of `End(E)`: `conductor_divides` (a bound `f₀` on its conductor), `conductor` (exact, when resolved), `maximal` (`End(E) = O_K` established), the volcano `levels` at each prime of `v`, and the orders it is known to `contain` | see below |
| `small_primes` | for each `ℓ ≤ 31`: `v_ℓ(v)` (the depth of the `ℓ`-isogeny volcano), how `ℓ` splits in `Q(π)`, and the number of rational `ℓ`-isogenies when the depth is 0 | exact from `Δ` modulo powers of `ℓ`, without factoring `Δ` |
| `subfield` (binary) | `j_field_degree`: least `k` with `j ∈ GF(2^k)`; `definition_degree`: least `k` such that the curve is the base change of one over `GF(2^k)`; `base_trace`: that curve's trace `t_k` | see below |
| `subgroup` | prime subgroup order `r` and cofactor | the registry representation, else the largest prime of `#E` |
| `embedding` | `ord_r(q)` and the bit length of `(r − 1)/k` | exact from the factorisation of `r − 1`; else `q^k ≢ 1` for `k ≤ 1000` gives a bound |
| `twist` | the quadratic twist's order `q + 1 + t`, its largest prime and cofactor | factored with the same budget |
| `dissect` | the DiSSECT trait families this module can define exactly (see "DiSSECT families") | |
| `keys` | the grouping keys below | |

Montgomery and Edwards registry equations over `GF(p)` are mapped to short
Weierstrass form for every check; other catalogue models retain their exact
model JSON, characteristic, degree, field size, trace, and order identity,
and their checks are `not_evaluated`.

**Curves over `GF(p^k)`, `k > 1`** (`y² = x³ + ax + b` in the registry's
polynomial basis, `3 < p < 2^63`) have Jacobian group arithmetic
([`torsion.rs`](../../src/cryptanalysis/curve_traits/torsion.rs), `FpK`).
A recorded generator certifies the order as above. Without one, the
**point-order certificate** takes points at fixed pseudo-random
x-coordinates, checks `[N]P = O` for the recorded order `N`, and computes
each point's exact order from the factorisation of `N`; once the lcm `L` of
those orders exceeds `4√q`, `N` is the only multiple of `L` in the Hasse
interval, so it is `#E`, as strong as the primes of `N`. A point `N` does
not kill stops the build; a lcm that stays below `4√q` gives `bounded`.
Exhaustive counting over `GF(p^k)` is not implemented, so a curve whose
order does not factor stays `not_evaluated`.

**Subfields.** `E_{a,b} : y² + xy = x³ + ax² + b` over `GF(2^n)` is the base
change of a curve over `GF(2^k)` exactly when `b ∈ GF(2^k)` and either `n/k`
is odd or `Tr(a) = 0`. So `j` can lie in a subfield the curve does not
descend to: the quadratic twist of a Koblitz curve over an even degree has
`j = 1` but descends only to `GF(4)`. When the curve descends, `t_k` is
found as the unique Hasse-bounded root of `V_{n/k}(t_k, 2^k) = t_n`, or by
counting points over `GF(2^k)` when the root is not unique. When `n/k` is
even, both twists over `GF(2^k)` descend and only `|t_k|` is an invariant
(`base_trace_sign_free`).

**Three orders, kept apart.** `Δ = t² − 4q` determines `Z[π]` and its
conductor `v`; `d_K` determines the maximal order `O_K`; `End(E)` lies
between them, and its conductor `f` divides `v` but is not fixed by `Δ`
alone. (For an ordinary curve every geometric endomorphism is defined over
the working field, so the two rings coincide.) The `endomorphism` record
bounds `f` by every order known to lie in `End(E)`: `Z[π]`; the Frobenius
`π_k` of a field of definition `GF(2^k)`, since `π = π_k^{n/k}`, giving
`f | v_k` with `t_k² − 4·2^k = v_k²·d_K`; `Z[ζ₃]` for `j = 0` and `Z[i]` for
`j = 1728` on a prime-field model; and a registry certificate, which fixes
it. When the bound reaches 1, `End(E) = O_K` is established: every Koblitz
curve (`Z[τ]` is already maximal), the registry curves defined over `GF(4)`,
secp256k1.

**Volcano levels.** Otherwise each prime `ℓ` of `f₀` is resolved on its own
(`endomorphism.levels`; the level is `v_ℓ(f)`, 0 on the crater and the depth
`v_ℓ(v)` on the floor):

- a level the source proves (the walker's `proved_volcano_levels`, the
  crosswalk's `volcano_levels`) is taken, and checked against the test
  below wherever the test runs; a disagreement stops the build;
- otherwise the **scalar test**
  ([`torsion.rs`](../../src/cryptanalysis/curve_traits/torsion.rs)): with
  `ℓ | v`, Frobenius has the double eigenvalue `λ ≡ t/2` on `E[ℓ]` and is
  either the scalar `λ` or a Jordan block, and it is the scalar exactly when
  `ℓ | [End(E) : Z[π]]`, i.e. above the floor. With `k` the order of `λ`
  mod `ℓ`, the scalar case is `x^{q^k} ≡ x (mod ψ_ℓ)`; a Jordan block
  fails it. Not scalar fixes the level at the depth; scalar at depth 1
  fixes it at 0; scalar at greater depth only bounds it by `depth − 1`.
  The test runs over `GF(p)`, `GF(2^n)` and `GF(p^k)` for `ℓ ≤ 1000`
  within a fixed cost budget (`SCALAR_TEST_COST`), and is compared with a
  brute-force count of rational `ℓ`-torsion over small prime, binary and
  `GF(p^k)` fields in the tests.

When every prime is resolved, `conductor` is exact and `end_order` is
`maximal` or `non-maximal`; otherwise `End(E)` stays `bounded` between
`Z[π]` and `O_K`, with `conductor_divides` the product of the per-prime
bounds. Primes above the test's reach (the large primes of `v` on most
large curves) need a walk (`src/bin/isogeny_walk.rs`).

**What kind of claim each value is.** Every value is an exact invariant or
a bound on one, in the sense of
[`docs/endomorphism-rules.md`](../endomorphism-rules.md): none is a proposed
route, an implemented map, a verified subgroup action or a measured
advantage. A small class number, embedding degree, cofactor or conductor
prime is a lead with specific hypotheses, not a verdict. Two curves with
equal keys are alike in those invariants only: a shared signature supplies
no isomorphism, isogeny or transfer between them, and two registry models
of one curve under different moduli stay two records, as ICV1 requires.

## Status of a value

| Status | Meaning |
|:--|:--|
| `proved` | exact; every primality claim under it is deterministic (Miller–Rabin to the bases 2…41 is a proof below 3.3·10²⁴) |
| `probable` | exact if every probable prime under it is prime |
| `bounded` | only a bound is known: `conductor` is then a lower bound and `cm_disc` is absent |
| `unknown` | attempted and not determined within the factoring budget |
| `not_evaluated` | skipped by a declared limit (`h(d_K)` above `|d_K| = 10^7`, a subfield above degree 40), by an unimplemented group law, or (over `GF(p^k)`) by an order with no generator that did not factor |
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
| `end_order` | `maximal` when `End(E) = O_K` is established; `non-maximal` when a smaller order is established; `unresolved`; `unknown` |
| `isogenies` | rational `ℓ`-isogenies over `GF(q)` for `ℓ = 2, 3, 5, 7` (`?` unknown, `-` the characteristic) |
| `source` | `registry`, `isogeny-walk` or `ic-crosswalk` |

A curve's **signature** is its values of `char, ordinary, cm, descent,
cofactor`, the default of `group --by`.

## Similarity

`similar <curve>` ranks every other curve by the summed weight of the keys
on which it differs from the target: `cm` and `descent` 4, `char` and
`ordinary` 3, `jfield` and `cofactor` 2, `twist_cofactor`, `embedding`,
`conductor`, `split` and `depth` 1 (`family` and `descent_signed` are not
compared by default). Ties are broken by
`|Δ trace_ratio| + |Δ conductor_fraction|`. A key that is `unknown` on
either side counts as a difference and is shown with a `?`: two curves
whose discriminants did not factor are not thereby alike. `--keys` compares
on a chosen set instead, each weighted 1; `--other-sizes` leaves out curves
over a field of the target's size.

## Current catalogue and earlier classification panel

The generated index covers every registered model. It includes the
100 independently verified large-prime isogeny codomains. Each record retains
the status of its order certificate and every derived trait. Models whose
coordinate/order adapter is unavailable explicitly report `not_evaluated`;
Curve25519 is one such record. A complete catalogue entry does not by itself
certify its coordinates or group order.

In the current file (2112 registry curves and 176 from the sources: 175
walked curves of the P-256 and P-224 isogeny classes and the ECC2K-130
degree-263 descendant), `End(E)` is established for 1014 curves (670
maximal, 344 non-maximal) and unresolved for 1239; 35 have a `Δ` that did
not factor. Of the 1324 curves over `GF(p^k)`, `k > 1`, 1319 have a
certified order (9 by generator, 1310 by point orders: 484 proved, 826
probable) and 301 a resolved `End(E)`. The unresolved are mostly levels the
scalar test only bounds (1184 levels "above the floor" at depth 2 or more,
almost all at `ℓ = 2`), then 108 levels beyond its cost limit and 72 at a
prime of `v` above its range. The volcano levels behind it: 1246 decided
by the scalar test, 13 walker certificates (each agreeing with the test),
111 from a certified `End(E)`, 228 where `f₀` is already prime to `ℓ`. P-256's `Δ` factors with
the committed hints as `−3·5·456597257999·p81·p136` (squarefree), so
`v = 1`, `d_K = Δ` and `End(E) = O_K`, probable through its 136-bit prime
(the 81-bit one is below the Miller–Rabin proof bound); every walked curve
of its class shares that. The ECC2K-130
descendant has `End(E)` of conductor 263: level 1 in its 263-volcano, level
0 at the other prime of `v`.

The following classification panel describes the earlier 121-model registry
snapshot, with registry digest
`a43b9dd3529815ae68a90449f1ce76007a529421a1b0c1b078f99d62d698beed`.
Its counts refer to that snapshot. The current generated file provides the
expanded records and their individual certification statuses.

From `group --by cm,descent --min-sizes 2` on that 121-model snapshot:

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
  P-256's `Δ` did not factor within the budget then (it does now, with the
  factor hints).
- `End(E) = O_K` is established for 114 of the 121 curves, 61 of them
  where `Z[π]` itself is not maximal: all 59 Koblitz curves through `Z[τ]`,
  the four curves defined over `GF(4)` through their `GF(4)` Frobenius, and
  secp256k1 through `Z[ζ₃]`. Six remain unresolved between `Z[π]` and
  `O_K`; P-256 is unknown with its `d_K`.
- Every order in that snapshot was certified: 60 by exhaustive count, 41 by
  a descent count, the rest by generator certificates (4 of those rest on
  a probable prime).

## Sources

Besides the registry, the build reads every file named in
[`trait-sources.json`](trait-sources.json):

- **`isogeny-walk`**: the `curve_nodes` of an isogeny walk's
  `isogeny_routes.json` (the P-256 `ℓ = 31` and P-224 `ℓ = 61` radius-2
  runs under `research/isogeny_walk_20261004/`). Each node's ICV1 slug is
  recomputed from its model and order and must equal the walker's; its
  generator certifies the order like a registry representation, and its
  proved volcano levels feed `endomorphism.levels`.
- **`ic-crosswalk`**: the entries of [`ic/curves.yaml`](ic/curves.yaml) that
  are not registry models (today the ECC2K-130 degree-263 descendant at
  level 1). The binary model `y² + xy = x³ + a₂x² + a₄x + a₆` is taken to
  the registry form by `y ↦ y + a₄` (`b = a₆ + a₄²`, `G ↦ (x, y + a₄)`);
  its proved `End(E)` discriminant and volcano level are kept. The file is
  read by a strict YAML-subset reader that refuses anything it does not
  understand, including a repeated key inside an entry it reads.

These slugs are computed identities, not registrations (`source.registered`
is `false`): cite them in prose only after registering them (`AGENTS.md`
§11). A source curve whose slug is already registered is skipped, so the
registry record wins.

## DiSSECT families

[`dissect.rs`](../../src/cryptanalysis/curve_traits/dissect.rs) computes the
families of DiSSECT (Sedlacek et al.) that follow from `t`, `q`, `r` and the
volcano levels, under the definitions written there; parity with DiSSECT's
own code is not claimed. Every value carries the order check's status
(`dissect.status`).

| Field | DiSSECT family | Value |
|:--|:--|:--|
| `torsion[ℓ]`, `ℓ ≤ 19` | `torsion_extension`, `isogeny_extension`, `isogeny_neighbors`, `division_polynomials` (`ℓ ≤ 5`), `volcano` | Frobenius mod `ℓ` (`split`, `inert`, `jordan`, `scalar`); least `k` with an `ℓ`-torsion point, with all of `E[ℓ]`, with one `ℓ`-isogeny, with all `ℓ + 1` of them over `GF(q^k)`; rational `ℓ`-isogenies; the degrees of `ψ_ℓ`'s irreducible factors |
| `conductor_ratio[d]`, `d = 2, 3, 4` | `conductor` | `|U_d(t, q)| = [Z[π] : Z[π^d]]` |
| `kn_factorization` | `kn_factorization` | `k·r ± 1`, `k ≤ 8` |
| `trace_factorization[d]` | `trace_factorization` | `t_d = V_d(t, q)`, `d = 1, 2` |
| `twist_order[d]` | `twist_order` | `q^d + 1 + t_d`, `d = 1, 2` |
| `small_prime_order_index[ℓ]` | `small_prime_order` | `(r − 1)/ord_r(ℓ)`, `ℓ ≤ 13`, when `r − 1` factors |
| `pow_distance` | `pow_distance` | `#E` and `r` against the nearest power of 2 |
| `x962_invariant` | `x962_invariant` | `a³/b² mod p` |
| `hamming_x[w]` | `hamming_x` | how many `x` of Hamming weight `w` are x-coordinates of points, where at most 600 values are tried |

Frobenius on `E[ℓ]` is fixed up to conjugacy by `X² − tX + q mod ℓ` except
at a double root, where the volcano level decides between the scalar and a
Jordan block; an unresolved level leaves that `ℓ`'s entry empty. The
factorisations use trial division to 2¹⁶, the factor hints, 2¹⁰ rho
iterations and a primality test (the trace `t` gets the full budget), and
are written `2^2·3·p250`: a large prime as `p<bits>`, an unsplit
composite as `c<bits>`. `not_evaluated_families` in the file lists the
families left out and why (`brainpool_overlap`, `multiples_x`, `q_torsion`,
`square_4p1`, class-number bounds, `weierstrass`).

## Factor hints

[`factor-hints.txt`](factor-hints.txt) lists primes, one per line, that
PARI/GP found in composites the rho budget leaves. The build re-checks
each (Miller–Rabin, above 2¹⁶) and uses it only as a trial divisor before
rho, so a hint can complete a factorisation but cannot change what is
proved. Regenerate with [`factor-hints.gp`](factor-hints.gp):

```bash
./target/release/curve_traits unfactored > /tmp/unfactored.txt
UNFACTORED=/tmp/unfactored.txt HINT_SECONDS=60 \
  gp -q -s 1000000000 docs/curves/factor-hints.gp < /dev/null > /tmp/hints.txt
mv /tmp/hints.txt docs/curves/factor-hints.txt
./target/release/curve_traits build
```

## Keeping it current

`traits.json` is derived from `registry.json`, the factor hints, the source
manifest and its files; the file records each input's SHA-256. When any of
them changes, rebuild it in the same pull request:

```bash
cargo run --release --bin curve_traits -- build
```

The `curve-traits` workflow runs the module's tests and `curve_traits
check`, which rebuilds with the committed file's budget and fails on any
difference.

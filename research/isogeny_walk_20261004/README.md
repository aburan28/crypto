# Isogeny walks from P-256 and P-224

Status: **tooling and enumeration; no ECDLP cost measured**

Date: 2026-10-04

Result class: **not an index-calculus measurement.** No IC/rho ratio, no
phase is priced, and the scoreboard and leaderboard carry no row for this
work (`AGENTS.md` §7a lists no panel for a curve enumeration).  The
registry gains P-224 (§11), so the leaderboard roster and the lab browser
are regenerated in the same change.

## What this is for

PR #1330 preregisters a `2^40`-record P-256 isogeny search.  Its Gate 1
blocks production until a **native Rust P-256-scale isogeny engine**
exists, reproduces fixtures, verifies every emitted edge independently of
the worker's choice logic, and preserves the P-256 trace and identity.
This directory lands that engine and its first runs.  It also gives any
thread a way to enumerate a prime-field curve's isogeny class in the
repository's curve formats:

- ICV1 slugs (`docs/curves/ICV1.md`);
- EC1 identities (`docs/curve-identities.md`);
- `docs/curves/ic/curves.yaml` records, the format mirrored from
  cryptanalysis in PR #1332 and integrated here;
- `IW1` routes (cryptanalysis `VOLCANO_NAMING.md`).

## Method

`src/cryptanalysis/isogeny_walk/`, driven by `src/bin/isogeny_walk.rs`.

1. **Field.** `GF(p)` for any odd `p < 2^256` in four-limb Montgomery
   form.  The same code serves P-256, P-224 and the small primes the tests
   brute-force against.
2. **`Φ_ℓ mod p`.** Built from the `q`-expansion of
   `j = 1728·E4³/(E4³ − E6²)`.
   - The symmetric unknown `c_ab` first appears at `q^{−(a+ℓb)}` with
     coefficient 1, and those pole orders are distinct, so the system is
     triangular and needs no division.
   - Every other exponent from `q^{−ℓ(ℓ+1)}` to `q^8` must vanish, and
     does.  The runs record the count; for `Φ_59` it is 1,720.
3. **Neighbours.** The `F_p`-roots of `Φ_ℓ(j, Y)`, by
   `gcd(f, Y^p − Y)` and Cantor–Zassenhaus.
4. **Kernel.** Elkies' construction turns each root into an explicit
   kernel polynomial.
   - Differentiating `Φ_ℓ(j(τ), j(ℓτ)) = 0` once gives the normalised
     codomain.
   - Differentiating twice gives `E2 − ℓE2'`.
   - The `q`-expansion of `℘` at the `ℓ`-torsion gives
     `Σ x(Q) = −ℓ(E2 − ℓE2')`.
   - Vélu's `℘_{E'} = ℘_E + Σ[℘_E(z+Q) − ℘_E(Q)]` gives the remaining
     power sums.
   - The derivation is in `kernel.rs`.
5. **Certificate.** `kernel::verify_kernel` never looks at `Φ_ℓ`.  It
   accepts `h` only if `h`:
   - is squarefree of degree `(ℓ−1)/2`;
   - divides the `ℓ`-division polynomial;
   - has roots closed under `[g]` for a generator `g` of `(Z/ℓ)^*/±1`.

   Together these make the roots exactly the `x`-coordinates of one cyclic
   subgroup of order `ℓ`, Galois-stable because `h ∈ F_p[x]`.  Vélu's
   codomain must then be `F_p`-isomorphic to the recorded target.
6. **Order.** The isogeny is `F_p`-rational, so the codomain has the
   root's order.  The audit `[n]P = O` additionally proves `#E = n`,
   because `n` is prime and exceeds the Hasse interval's width.
7. **Levels and directions.** For odd `ℓ`, the `ℓ`-volcano has depth
   `⌊v_ℓ(t² − 4p)/2⌋`.
   - Depth 0: every curve is at level 0.
   - Depth 1: `ℓ + 1` rational `ℓ`-isogenies means the surface, one means
     the floor.

   An `IW1` edge is `h`, `d` or `u` only when both endpoint levels are
   proved; otherwise it is `x`.
8. **Models and generators.** These follow fixed, recorded rules:
   - `icwalk-canon/v1`: the registered model for the root, else an `a = −3`
     model with the least `b`, else a `j`-model;
   - `icwalk-gen/v1`: the least-`x` point.

   The generator is **not** the image of the root's generator, so the walk
   establishes no discrete-log transport.

### Validation

`cargo test --release --lib isogeny_walk` runs 10 tests of the walker
(the filter also picks up 2 older `jv_isogeny_walk` tests; all 12 pass).
They cover:

- **`Φ_2`, `Φ_3`:** equal the literature tables reduced mod the P-256
  prime.
- **`j` series:** equals `q^{-1} + 744 + 196884q + …`.
- **Brute force over `F_1009`:** for `ℓ = 3, 5, 7`, kernels built from
  rational `ℓ`-torsion points with independent `u64` arithmetic pass the
  certificate and land on a root of `Φ_ℓ(j, Y)`.  The Elkies construction
  from that root reproduces the brute-force kernel exactly (at least 8
  curves per `ℓ`).
- **P-256 root counts:** the number of roots of `Φ_ℓ(j, Y)` equals what
  the splitting of `ℓ` in `Z[π]` predicts for `ℓ ≤ 13`.
- **Certificate rejections:** a perturbed kernel, and a kernel mixing two
  subgroups of a curve over `F_211` with full rational 5-torsion.
- **P-224 replay:** a short P-224 walk re-verifies from its own records.
  A tampered kernel, a relabelled direction and a forged route id each
  fail the replay.  P-224's surface position and its 1 horizontal plus 3
  descending 3-edges are pinned.
- **EC1:** the P-256 EC1 identity equals the registry's
  (`EC1P256Cp256h0523b774e066`).

Both committed runs' `curves.yaml` validate against
`docs/curves/ic/curves.schema.json`.  This was an ad hoc local check; no
validator is committed.  The Python registry builder independently gives
P-224 the slug and EC1 identity the Rust walker emits.

## Frozen inputs

| curve | ICV1 slug | EC1 |
|:--|:--|:--|
| P-256 | `icv1-fp256-t89188191154553853111372247798585809583-f188c491` | `EC1P256Cp256h0523b774e066` |
| P-224 | `icv1-fp224-t4733100108545601916421827343930821-fec01e99` | `EC1P224Cp224h8e110c585480` |

The primes are every odd `ℓ ≤ 61` (31 for the committed P-256 run).  Atkin
primes have no rational `ℓ`-isogeny and are skipped; each run's
`walk.json` lists them with the reason.

**Success condition.** Every root of `Φ_ℓ(j, Y)` yields a certified edge,
and every root count equals the splitting prediction (depth 0) or is
`ℓ + 1` or 1 (depth 1).  Every curve's order is proved, and the
independent replay passes.

**Stop condition.** Any failure.  None occurred.

## Class invariants

From `walk.json#class`; they hold for every curve in each class.

| | P-256 | P-224 |
|:--|:--|:--|
| trace `t` | 89188191154553853111372247798585809583 | 4733100108545601916421827343930821 |
| `t² − 4p`, trial factors `< 2^20` | `−3 · 5 · C₂₅₅` (255-bit cofactor) | `−3³ · 29 · 79 · 7523 · 40927 · C₁₈₂` (182-bit cofactor) |
| cofactor `C` | composite, unfactored | composite, unfactored |
| `End(E)` | `unk` (depends on `C`) | `unk` |
| twist order | `3·5·13·179·q`, `q` a 241-bit probable prime | `3²·11·47·c`, `c` a 212-bit composite, unfactored |
| embedding degree | `> 1000` | `> 1000` |
| Elkies `ℓ ≤ 61` | 11 13 17 23 29 37 41 43 47 59 | 11 47 59 61 |
| `ℓ \| t² − 4p` | 3, 5 (depth 0, one horizontal edge) | 3 (depth **1**), 29 (depth 0) |
| Atkin `ℓ ≤ 61` | 7 19 31 53 61 | 5 7 13 17 19 23 31 37 41 43 53 |

**P-224 sits on the surface of a depth-1 3-volcano.**  `v_3(t² − 4p) = 3`
forces `3 | [O_K : Z[π]]` exactly once.  P-224 has four rational
3-isogenies:

- one horizontal (3 is ramified in `O_K`);
- three descending to curves with `3 | [O_K : End(E)]`.

Those floor curves (`V3L1`) are in P-224's isogeny class.  No walked P-256
volcano has depth above 0.

## Runs

| run | curves | expanded | edges | failures | orders proved | replay |
|:--|--:|--:|--:|--:|--:|:--|
| [`p256-ell31-radius2`](runs/p256-ell31-radius2/) | 84 | 13 | 156 | 0 | 84 | pass |
| [`p224-ell61-radius2`](runs/p224-ell61-radius2/) | 93 | 14 | 173 | 0 | 93 | pass |
| [`p256-ell61-20k`](runs/p256-ell61-20k/) | 20,000 | 3,576 | 78,672 | 0 | 20,000 | pass |
| [`p224-ell61-20k`](runs/p224-ell61-20k/) | 20,000 | 10,602 | 117,840 | 0 | 20,000 | pass |

The radius-2 runs commit all three outputs.  The 20k runs commit
`walk.json` and the replay receipt.  Their 333 MB (P-256) and 407 MB
(P-224) of raw output are recorded in [`MANIFEST.json`](MANIFEST.json) by
SHA-256 and byte count.  The walk is deterministic, so the
command there regenerates the same bytes from this commit: the outputs
were produced twice, before and after a summary-only change, with
identical hashes.  They are not uploaded anywhere, so the
regeneration command is the only durable route to them.

Root counts on every expanded curve match the prediction:

- **P-256:** 1 root for `ℓ = 3, 5` and 2 for each Elkies prime, on all
  3,576 curves.
- **P-224, `ℓ = 3`:** 4 roots on 3,940 curves (surface) and 1 root on
  6,662 curves (floor); 9,398 unexpanded curves stay unresolved.

Traits across the 20k runs (`walk.json`):

| | P-256 | P-224 | expectation |
|:--|--:|--:|:--|
| curves with an `a = −3` model | 10,195 (51.0%) | 4,937 (24.7%) | ½ for `p ≡ 3 (4)`, ¼ for `p ≡ 1 (4)` |
| `qr_prefix_64` range | 16–48 | 16–50 | `Binomial(64, ½)`: mean 32, sd 4 |

The `qr_prefix_64` histograms in `walk.json` look like the binomial.  The
one P-224 curve at 50 meets PR #1330's sizing threshold.  At
`P(≥ 50) = 3.5·10⁻⁶` per curve, that is about 0.07 expected in 20,000
draws.  The statistic depends on the model, so the figure is a property
of the recorded model only.  It is not evidence about any DLP.

Practicality note, not a measurement: these timings are from
unisolated runs, not through `tools/isolated_bench.py`.  Host: Apple M4
Pro, 14 cores, 48 GiB, macOS 26.6, rustc 1.93.1.

- **`Φ_ℓ mod p` build:** 10.1 s for `Φ_59`, with the primes built in
  parallel.
- **Walks:** P-256 took 115 s for 20,000 curves with 12 primes; P-224
  took 231 s.
- **Replay:** 31 s and 71 s.

## What this does not establish

- No discrete-log transport between curves.  The recorded generators
  follow `icwalk-gen/v1`, not the isogenies.
- No change in DLP cost on any walked curve.  Any such claim needs the
  complete one-target protocol of `AGENTS.md` §8 and §12 against matched
  rho.
- **`End(E)` beyond the walked primes.**  Both discriminants keep a
  composite cofactor, so `End(E)`, the full volcano depths and the
  component structure stay `unk`.
- **Registered identities.**  Walked curves carry computed ICV1 slugs and
  EC1 identities but are not registered.  Register one before citing it in
  prose (§11).
- **Depth-1 level claims.**  They rest on the root counts the walk
  observed.  The replay checks them for consistency with edge directions
  and the class depth, but does not recompute them; that would need
  `Φ_ℓ`.

## Next

PR #1330's Gate 1 asked for:

- a native engine;
- independent per-edge verification;
- preservation of the trace and the registered identity;
- frozen small-field and P-256 fixtures.

This change supplies all four (the tests and the radius-2 runs are the
fixtures).  Gate 2, a bounded `2^24`-record pilot with storage and cost
measurements, is next.  It should freeze its screening statistic before
reading any walk output.

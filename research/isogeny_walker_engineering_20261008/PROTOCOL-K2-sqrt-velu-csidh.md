# Protocol K-2: √élu in the CSIDH toy for large prime degree

Frozen 2026-10-08, before any instrument is built.  **Engineering.**  No
security number changes; the toy's parameters stay toy.  Status: **PENDING**.

## Derivation (stated before measuring)

Vélu's formulas evaluate an `ℓ`-isogeny with kernel `⟨P⟩` in about `ℓ`
field operations by visiting every kernel point.  Bernstein, De Feo,
Leroux and Smith (ANTS 2020) evaluate the same isogeny in `Õ(√ℓ)`
operations: on a Montgomery curve the product `h_S(X) = ∏_{s∈S}(X − x([s]P))`
over `S = {1, 3, …, ℓ − 2}` is rewritten through an index-set
decomposition `S ⊇ I ± J` as a resultant of two polynomials of degree
about `√ℓ`, using the biquadratic relation between `x(P+Q)`, `x(P−Q)` and
`x(P), x(Q)`; the codomain coefficient and point images follow from
`h_S(1)`, `h_S(−1)` and `h_S(x(Q))`.  The crossover with plain Vélu is
near `ℓ ≈ 100` in the authors' measurements; CSIDH-512's largest prime is
587.

The toy at `p = 419` uses `ℓ ∈ {3, 5, 7}`, below any crossover, so the
protocol first adds a parameter set with large `ℓ` while keeping
everything else toy: `p = 4·∏ℓ_i − 1` with `ℓ_i` the first 20 odd primes
plus 587, chosen so that `p` is prime and under 128 bits, with the field
in the repository's existing multi-word arithmetic.

## Instrument (Rust, to build in the follow-on PR)

1. `sqrt_velu_isogeny(A, x(P), ℓ, points)` on Montgomery curves, with
   the index sets, the biquadratic `F₀, F₁, F₂`, a resultant over
   `F_p[X]` of degree about `√ℓ` by product trees, and the codomain and
   image formulas; plain Vélu kept as the reference.
2. Equality tests: for every `ℓ_i` of the large parameter set and 100
   random kernel points each, the √élu codomain and the images of 10
   random points agree with plain Vélu exactly.
3. A paired benchmark, instructions and wall, per `ℓ_i`, both
   implementations, under `tools/isolated_bench.py` where the host
   allows.

## Predictions (pass/fail)

- **K2-1 (equality).**  Zero disagreements with plain Vélu on every
  `(ℓ_i, kernel)` tested.
- **K2-2 (crossover).**  √élu is slower than Vélu at `ℓ ≤ 31` and faster
  at `ℓ ≥ 149`, in instructions.
- **K2-3 (gain).**  At `ℓ = 587` √élu costs at most `0.4×` Vélu's
  instructions.
- **K2-4 (action).**  One full class-group action on the large set,
  with √élu for `ℓ ≥ 149` and Vélu below, costs at most `0.7×` the
  all-Vélu action.

## Decision rule and inadmissible moves

K2-1 is a gate; K2-2 to K2-4 size the gain.  Class **engineering**.  The
toy remains non-constant-time and non-interoperable, and says nothing
about CSIDH's security or about any hybrid KEM.  Inadmissible: changing
the parameter set between arms; counting kernel-point sampling in one
arm only; reporting wall time without instruction counts.

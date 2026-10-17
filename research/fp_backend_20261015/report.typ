#set page(paper: "a4", margin: (x: 1.6cm, y: 1.6cm), numbering: "1")
#set text(size: 9pt)
#set par(justify: true)
#show heading.where(level: 1): set text(size: 14pt)
#show heading.where(level: 2): set text(size: 11pt)

#let meta = read("runs/final/meta.txt")
#let claims = csv("runs/final/claims.csv")
#let summary = csv("runs/final/summary.csv")

= Variable-time prime-field backends (`fp_backend`)

Study `research/fp_backend_20261015` · 2026-10-15 · module `src/cryptanalysis/fp_backend.rs` ·
harness `examples/fp_backend_bench.rs`, analysis `examples/fp_backend_report.rs` ·
protocol `PROTOCOL.md` (committed before the measured run).

== What was built

One trait, `PrimeField`, with element arithmetic, inverse, batch inverse, exponentiation,
Legendre symbol, square root and a *lazy accumulator* (`F::Wide`: unreduced sums of
products, one reduction at the end) with row kernels (`row_load`, `row_sub`, `row_axpy`,
`row_scatter_axpy`, `row_reduce`, `row_reduce_scaled`, `dot`). Backends: `Fp32` (canonical
`u32`, Barrett; the lazy kernel of `f4_fp_tower` below 2#super[31], Shoup products with a
carry fold above), `Fp64` (Montgomery, `p < 2^63`), `FpWide<N>` (Montgomery CIOS, any
`p < 2^(64N)`), `FpWide<N, true>` (carry-free CIOS for primes with a spare top bit, alias
`Fp224`), `FpP192` (Solinas reduction for NIST P-192). Inverse: Pornin's binary GCD
(wide) and Kaliski's almost inverse (one word); Legendre: binary Jacobi; square root:
`a^((p+1)/4)` or Tonelli–Shanks with 8/12-bit windowed 2-Sylow logarithms when
`v₂(p−1) ≥ 16`. Lazy-accumulator capacities are derived in the module documentation and
asserted in debug builds: unbounded for `Fp32`/`Fp64` (self-renormalising), `2^64 − 2`
products for `FpWide`/`FpP192`.

== Correctness

`cargo test --release --lib fp_backend`: 8 tests, every operation of every backend
against `BigUint` on edge and random values at 3–17, 31/32-bit, 61–63-bit (two with
`v₂(p−1) = 40`), P-192, a 190-bit prime, P-224, P-256, a 200-bit prime, 2#super[127]−1,
2#super[128]−159, one-limb primes to 2#super[64]−59 and a 512-bit prime; lazy rows,
scatter rows and dot products against eager reduction including all-(p−1) inputs; the
binary GCD never fell back. Same suite passes with debug assertions and overflow checks on.
The harness re-checks every timed kernel against `BigUint` or the f4_fp_tower copy before
timing.

== Run

#raw(meta)

Wall-clock ns per operation, 9 interleaved repetitions of every arm (≈ 40 ms each), one
core (`taskset -c 3`) of a *shared, noisy* 4-vCPU host with another lane running: no
isolation (L0). Ratios are paired by repetition (baseline/candidate, > 1 = candidate
faster); verdicts on the median paired ratio against the frozen thresholds.

#figure(image("fig_ns_per_op.svg", width: 100%),
  caption: [Median ns per operation (whiskers: min–max over 9 repetitions), run
  `runs/final`. Blue: this module; gray: existing code (`f4_fp_tower` copy,
  `prime_fast::Fp64`, `ct_bignum`, `BigUint`).])

#pagebreak()
== Pre-registered claims

#let rows = claims.slice(1)
#table(
  columns: (auto, 1fr, auto, auto, auto, auto, auto, auto),
  inset: 3pt,
  stroke: 0.4pt + gray,
  align: (left, left, right, right, right, right, right, center),
  table.header([*id*], [*claim*], [*base ns*], [*cand ns*], [*ratio*], [*min–max*], [*thr*], [*verdict*]),
  ..rows.map(r => (
    r.at(0), text(size: 7.5pt, r.at(1)),
    str(calc.round(float(r.at(4)), digits: 2)),
    str(calc.round(float(r.at(5)), digits: 2)),
    str(calc.round(float(r.at(6)), digits: 3)),
    str(calc.round(float(r.at(7)), digits: 2)) + "–" + str(calc.round(float(r.at(8)), digits: 2)),
    r.at(10),
    if r.at(11) == "PASS" [*PASS*] else [*#text(fill: red)[FAIL]*],
  )).flatten()
)

#figure(image("fig_claims.svg", width: 92%),
  caption: [Median paired ratio per claim with min–max over repetitions and the
  frozen threshold (tick). Run `runs/final`.])

#pagebreak()
== All arms (median ns per operation, min, max)

#table(
  columns: (1fr, auto, auto, auto),
  inset: 2.5pt,
  stroke: 0.3pt + gray,
  align: (left, right, right, right),
  table.header([*arm*], [*median*], [*min*], [*max*]),
  ..summary.slice(1).map(r => (
    text(size: 7pt, raw(r.at(0))),
    text(size: 7pt, str(calc.round(float(r.at(5)), digits: 2))),
    text(size: 7pt, str(calc.round(float(r.at(8)), digits: 2))),
    text(size: 7pt, str(calc.round(float(r.at(9)), digits: 2))),
  )).flatten()
)

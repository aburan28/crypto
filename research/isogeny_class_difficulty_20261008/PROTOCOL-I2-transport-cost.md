# Protocol I-2: the cost of moving a DLP across a prime-degree isogeny, against `√r`

Frozen 2026-10-08, before any instrument is built.  **Pricing; stage
diagnostic.**  `S`, end-to-end cost and speedup are **unset**.

## Derivation (stated before measuring)

Transporting `(P, Q)` from `E` to an `ℓ`-isogenous `E′` costs

    T(E → E′) = F(ℓ) + 2·V(ℓ),

where `F` is the one-time cost of *finding* the edge (the kernel, or a
route that certifies the codomain) and `V` is the cost of evaluating `φ`
on one point.  In `F_q` multiplications:

| representation | `F(ℓ)` | `V(ℓ)` | validity |
|:--|:--|:--|:--|
| Vélu from a rational kernel point | cofactor multiplication in `E(F_{q^k})`, `k` the eigenvalue order; `k² log q` per step | `O(ℓ)` | kernel rational over `F_{q^k}` with small `k` |
| √élu [BDLS20] from the same kernel | same `F` | `Õ(√ℓ)` | same |
| `Φ_ℓ` route (walker, companion thread C1–C3) | `ℓ³`–`ℓ⁵` once per `ℓ`, then `ℓ²` per curve for the kernel | `O(ℓ)` via the kernel polynomial | `p > 4ℓ`, Elkies `ℓ` |
| class-group route (companion C5) | class-group precomputation once per class, then a product of `≈ h` small-degree steps | `Σ_i O(q_i)` over the route's small degrees, independent of `ℓ` | horizontal `ℓ` only |

The DLP itself costs `≈ 1.3·√r` group operations by rho.  So the transport
is worth considering only while `T(E → E′) < √r`, and *finding* a weak
neighbour by exhaustive walking costs the class size times the screen,
which is the real budget.  This protocol produces the table of
`T(E → E′)/√r` by representation and `ℓ`, at toy sizes where every column
can be run to the end, and the `ℓ` at which each representation's
transport alone exceeds rho's whole solve.

The point of the exercise is not the number but the shape: horizontal
large-`ℓ` transport is bounded by the class-group route and is therefore
*flat in `ℓ`*, while vertical large-`ℓ` transport has no such route and
is the one place where a large prime degree is a genuine cost barrier.

## Instrument (Rust, follow-on PR)

- Extend `isogeny_walk` with a `transport` subcommand: given a certified
  edge (or C5 route) and a point, evaluate the isogeny and count field
  operations; verify the image lies on `E′` and has order `r`.
- The √élu evaluator exists in `src/cryptanalysis/` for the binary Vélu
  work (`oriented-binary-velu`); the prime-field one is a follow-on.
- Charge in `ecbench`'s unit so the ratio to `√r` is in one column.

## Frozen inputs

| item | value |
|:--|:--|
| curves | the I-1 prime-field classes at `p ≈ 2^{20}` and `2^{28}`; one `D = −3` class from I-5 |
| degrees | every Elkies `ℓ ≤ 61`; from the companion E1 window, every `ℓ ≤ 2^{20}` with eigenvalue order `≤ 6`; `ℓ ∈ {101, 211, 1009}` through the `Φ_ℓ` route where it lands; one vertical `ℓ > 61` per class that has one |
| points | 8 random order-`r` points per edge |
| seed | 20261012 |

## Predictions (pass/fail)

- **Q1 (shape).**  `V(ℓ)` fits `c·ℓ` for Vélu and `c·√ℓ·log ℓ` for √élu
  with `R² ≥ 0.98` over the degree sweep.
- **Q2 (flatness of the class-group route).**  For horizontal `ℓ` in the
  sweep, the route cost varies by at most `3×` across two orders of
  magnitude of `ℓ`.
- **Q3 (crossover).**  At `p ≈ 2^{28}` (`√r ≈ 2^{14}`), the `Φ_ℓ` route's
  `F(ℓ)` exceeds `√r` group operations for every `ℓ ≥ 61`, and the
  rational-torsion route does so only for `k ≥ 4`.
- **Q4 (correctness).**  Every transported pair satisfies `φ(Q) = [d]φ(P)`
  on `E′` with the hidden `d` checked after the run; zero failures.

## Decision rule (registered)

- Q1–Q4 pass: **boundary**.  The transport table is recorded and read by
  the weak-curve methodology as `T`.  Horizontal large `ℓ` is a cost
  question only; vertical large `ℓ` is a reach question and is handed to
  I-5.
- Q2 fails: the class-group route does not flatten at toy size; record
  the norms, and the methodology charges the `Φ_ℓ` route instead.
- Q4 fails: a bug; halt.

## Stop condition and inadmissible moves

Bounded: three classes, the degree sweep above, 8 points per edge.

Inadmissible: charging only `V` and not `F`; charging a route found by a
different method than the one priced; counting a transport whose image
was not verified; extrapolating `F(ℓ)` beyond the sweep.

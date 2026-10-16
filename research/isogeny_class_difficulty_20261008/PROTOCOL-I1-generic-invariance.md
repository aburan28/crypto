# Protocol I-1: generic-attack invariance across an isogeny class, large prime degrees included

Frozen 2026-10-08, before any instrument is built.  **Control experiment;
stage diagnostic of nothing.**  `S`, end-to-end cost and speedup are
**unset**.

## Derivation (stated before measuring)

A separable `ℓ`-isogeny `φ: E → E′` over `F_q` with `ℓ ≠ r` restricts to a
group isomorphism `E[r] → E′[r]`.  Pollard rho, BSGS and the kangaroo use
only the group law and the negation map (and `Aut` folding where
`j ∈ {0, 1728}`), so their expected cost in group operations is a function
of `r` and `|Aut|` alone.  Prediction: the `S` column (`operations / √r`)
is the same on `E` and on every `E′` in its class, with `S_rho ≈ 1.3` under
the repository's accounting, up to the `√|Aut|/√2` factor at the special
`j`.  Any residual difference is field-arithmetic cost (which `ecbench`
splits out) or an implementation effect, and must vanish under an
isomorphic model.

Why this is worth a run at all: it is the null every other protocol in
this directory is read against, and it is the test of the harness's
claim that it compares curves by `r` and not by accident of model.

## Instrument (Rust, follow-on PR)

- `isogeny_walk walk` to enumerate the class at each prime-field toy size
  up to the walker cap; the C4 rational-torsion window of
  [`../large_prime_isogeny_degree_20261008/PROTOCOL.md`](../large_prime_isogeny_degree_20261008/PROTOCOL.md)
  to add, per class, every large prime `ℓ` whose eigenvalue order is
  `≤ 6` at toy size; one vertical edge per class where `f_π` has a prime
  factor `> 61`.
- `ecbench` spec per class: arms `rho-strong` on `E` and on each selected
  `E′` (small-`ℓ` neighbour, large-`ℓ` neighbour, vertical neighbour, and
  the isomorphic-model control of `E′`), one hidden target per cell, same
  `r`, interleaved, A/A arm first, L2 isolation.
- The binary ladder uses the native binary rho path on the registered
  prime-degree curves and their `Φ_ℓ mod 2` neighbours.

## Frozen inputs

| item | value |
|:--|:--|
| prime fields | three random classes each at `p ≈ 2^{20}, 2^{24}, 2^{28}`; cofactor `≤ 4` |
| binary | `icv1-f2m13-t181-515ee569`, `icv1-f2m19-t797-b6cf2467`, `icv1-f2m23-t5197-69e76b73` and three `ℓ`-neighbours each, `ℓ ∈ {3, 5, 7}` |
| targets per arm | 40 at `2^{20}`, 24 at `2^{24}`, 12 at `2^{28}`; binary 24 |
| seed | 20261011 |
| unit | ecbench group operations; field multiplications and inversions split |

## Predictions (pass/fail)

- **Q1.**  For every class and every selected `E′`, the paired bootstrap
  95% interval of `S(E′)/S(E)` contains 1, and its width is within the
  A/A interval of the class.
- **Q2.**  Large-`ℓ` and vertical neighbours are indistinguishable from
  small-`ℓ` neighbours on Q1: no ordering by `ℓ` or by edge direction
  reaches `p < 0.05` under a Kruskal–Wallis test across the class.
- **Q3.**  The isomorphic-model control of each `E′` matches `E′` on
  group operations exactly and differs only in the field-operation split.
- **Q4.**  Where a class contains a `j ∈ {0, 1728}` member (I-5's CM
  classes), `S` on it is `S(E)·√(2/|Aut|)` within the A/A interval, and
  nowhere else.

## Decision rule (registered)

- Q1–Q4 pass: **boundary (control confirmed)**; the harness compares by
  `r` and generic difficulty is class-invariant at these sizes.  Every
  other protocol here reads against this null.
- Q1 or Q2 fails on some `E′`: not a result.  It is a bug report against
  the harness or the walker (a wrong `r`, a non-separable edge, a target
  not in the subgroup); the run halts and the cell is traced before
  anything is classed.
- Q3 fails: implementation effect; recorded, and the model with the
  cheaper arithmetic is noted as such, never as an easier curve.

## Stop condition and inadmissible moves

Bounded: 9 prime-field classes, 3 binary curves, at most 6 arms per
class, one run.

Inadmissible: comparing cells with different `r`; a target not of order
`r`; pooling contended and uncontended wall times; reading a wall-time
difference as a difficulty difference; omitting the A/A arm.

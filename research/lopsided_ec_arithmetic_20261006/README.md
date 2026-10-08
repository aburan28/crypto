# Protocol: EC-arithmetic speedups from lopsided thin products? (round 5)

Frozen 2026-10-06. **Screening round**: no performance claim, no `S`
accounting, no `ecbench` sessions, no dashboard edits. New question from
the operator: can the Alman–Vassilevska Williams thin-product technique
(`arXiv:2610.06783v1`) speed up *elliptic-curve arithmetic itself*
(point/field/scalar-mul/MSM/pairing), rather than the index-calculus
pipeline of rounds 1–4 (merged PR #1485; open PRs #1493, #1494, #1501)?

Conductor note: no `conductor` binary exists in this environment, so no
scope reservation could be registered (as in rounds 2–4).

## Questions

1. Which EC-arithmetic kernels in this repository have the technique's
   required shape (many wanted inner products over one shared thin
   middle), ring (small-integer products), and scale (large N)?
2. What is the nearest miss, quantified — and what exactly kills it?
3. If nothing qualifies, what would a qualifying kernel look like (exact
   hook for a future revisit)?

## Scope and method

- Inventory the repo's EC kernels by file: point add/double and scalar
  mul (`src/ecc/point.rs`), batch Schnorr verification (`src/ecc/schnorr.rs`),
  bulletproof IPA (`src/zk/bulletproofs.rs`), KZG/polynomial
  (`src/zk/kzg.rs`, `src/zk/polynomial.rs`), extension fields and pairing
  (`src/bls12_381/`), field arithmetic, batch transforms.
- Screen each against the three gates in order (shape → ring → size),
  using the round-1 interface definitions and the round-4 envelope. A
  kernel failing shape has no (N, D, W) mapping at all; ring and size are
  only reached with one.
- Quantify the nearest miss with closed-form arguments (one-shot reading
  bound; entry-size bound), not measurements.
- No `src/` changes: a negative screen writes no code. A passing kernel
  would earn a measurement protocol per round 2's promotion rules; none
  is expected (pre-registered).

## Reference and boundaries

EC-arithmetic speed is counted in field/group operations on fixed-size
(256-bit-class) objects; the index-calculus unit `S` does not apply to
these kernels, and no cross-method claim is made. Boundaries for any
future positive: matched baseline/candidate operation counts on identical
inputs with identical outputs. Pre-registered expectation: zero kernels
pass shape, so no measurement round opens.

## Success and stop conditions

- Success: every inventoried kernel screened with the blocking gate
  named; nearest miss quantified; qualifying-kernel hook stated exactly.
- Stop: protocol, dated report with the screening table, figure plus
  sources, and PDF committed; PR opened. No speedup claimed.

## Deliverables

1. This protocol (`README.md`), the dated screening (`REPORT.md`).
2. `figure.mmd` (editable source) + `figure.svg`: kernel funnel and the
   (N, entry-size) feasibility plane.
3. `report.html` + `report.pdf` (report with visual included).

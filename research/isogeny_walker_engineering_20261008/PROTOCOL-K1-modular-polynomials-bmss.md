# Protocol K-1: precomputed modular polynomials and BMSS kernel polynomials for the class walker

Frozen 2026-10-08, before any instrument is built.  **Engineering.**  No
ECDLP cost, security number or ledger row changes.  Status: **PENDING**.

## Derivation (stated before measuring)

The walker finds the `ℓ`-isogenous neighbours of a curve as roots of
`Φ_ℓ(X, j(E))` in `F_p` and certifies each edge with an explicit kernel
polynomial.  Two costs bound it today:

1. `Φ_ℓ mod p` is rebuilt at every start, at a cost growing like `ℓ⁵`,
   which is why `ℓ` stops near 61.  `Φ_ℓ` over `Z` is a fixed object: its
   coefficients can be computed once per `ℓ` (Bröker–Lauter–Sutherland
   by CRT, or read from a published table where the licence permits),
   stored, and reduced mod `p` at start-up in time linear in the table
   size.  Sutherland's tables reach `ℓ` in the thousands; the size grows
   like `ℓ³ log ℓ` bits.
2. Given a root `j′`, the kernel polynomial is currently obtained by
   methods of cost about `ℓ²` field operations.  Bostan–Morain–Salvy–Schost
   compute the kernel polynomial of a normalised `ℓ`-isogeny in
   `Õ(ℓ)` operations from the two curve equations, and Elkies' procedure
   supplies the normalised codomain from `Φ_ℓ` and its partial
   derivatives.

Both are standard and both are absent here.  Together they change the
walker's per-step cost from about `ℓ²` to `Õ(ℓ)` plus a table read, and
its reach from `ℓ ≤ 61` to whatever the stored tables cover.

## Instrument (Rust, to build in the follow-on PR)

1. A table format for `Φ_ℓ` coefficients (sparse by monomial, integers
   as little-endian limbs), a loader that reduces mod `p`, and a
   generator that computes `Φ_ℓ` for `ℓ ≤ 199` by the CRT method and
   checks each against the known `Φ_2`, `Φ_3`, `Φ_5` and against the
   identity `Φ_ℓ(j(E), j(E/G)) = 0` on random curves with rational
   `ℓ`-torsion over small prime fields.
2. Elkies' normalised codomain from `Φ_ℓ` and its derivatives, and the
   BMSS kernel polynomial, each checked against the existing kernel
   certificate on every edge of a P-256 walk of 2,000 curves at
   `ℓ ≤ 61`.
3. The walker gains `--phi-tables DIR`; without the flag it behaves as
   today.

## Predictions (pass/fail)

- **K1-1 (start-up).**  With tables, start-up at `ℓ_max = 59` is below
  0.5 s (today about 10 s), and at `ℓ_max = 199` below 5 s.
- **K1-2 (step).**  Per-step wall at `ℓ = 59` on P-256 falls by at least
  `2×` against today's walker on the same host, with instruction counts
  reported beside wall time; at `ℓ = 199` the step costs at most `4×` the
  `ℓ = 59` step.
- **K1-3 (certificates).**  Every edge of a 2,000-curve P-256 walk at
  `ℓ ≤ 61` reproduces today's kernel certificate exactly; every edge at
  `61 < ℓ ≤ 199` certifies.
- **K1-4 (class coverage).**  A 20,000-curve walk at `ℓ ≤ 199` reaches
  at least `3×` as many distinct `j` as the same budget at `ℓ ≤ 61`.

## Decision rule and inadmissible moves

K1-1 to K1-3 passing makes tables the default; K1-4 sizes the gain.
Class **engineering**; the frontier and ledger do not move.  Inadmissible:
reading a faster walker as any statement about P-256's security; shipping
a table without the identity check; counting table generation as walk
cost.

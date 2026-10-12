# Three index-calculus tricks from an abliterated model, and how they die

**Date:** 2026-10-02
**Curve:** `E_0: y^2 + xy = x^3 + 1` over `GF(2^m)`, the `a = 0` Koblitz curve.
Trace recurrence as in `koblitz_point_count`: `s_0 = 2`, `s_1 = -1`,
`s_k = -s_{k-1} - 2 s_{k-2}`, `#E(F_{2^m}) = 2^m + 1 - s_m`.
**Model:** `abliterated-model-large-v2` on Abliteration.ai, `reasoning_effort: high`.
Two calls, `max_tokens` 4000 and 12000. Both finished `length` with an empty
message: every completion token was reasoning. The three tricks below are
taken from that reasoning trace, not from a finished answer. The trace is
not a result. The checks below are.

The prompt listed the dead ends already on the ledger (Semaev Weil descent
end to end, generic factor-base yield, subfields at composite `m`, and a
faster PDP stage) and asked for tricks those do not cover. The model did
not produce a method whose exponent beats rho. It produced three classical
ideas and, inside the same trace, the reasons they fail.

No row below is an `S` measurement. The scoreboard is unchanged. There is
no speedup, and none is claimed.

## Boundaries, before the checks

- **Floor.** A generic algorithm on a group of prime order `r` with
  automorphism group of order `A` needs about `sqrt(pi r / (2A))` group
  operations. On this Koblitz curve `A = 2m` (negation and Frobenius), so
  `S_floor = sqrt(pi / (4m))`.
- **Reference.** Pollard rho on the same curve, same operation accounting.
  The ledger's Koblitz rho sits near `S ≈ 1.3` once the automorphisms are
  used. Anything that does not beat that, with every phase priced and the
  logarithm verified, is not faster than rho.

## The three tricks

| rank (cheapest kill first) | trick | phase it would have to move | what was checked | class |
| --- | --- | --- | --- | --- |
| 1 | MOV / Frey–Rück: transfer the DLP through the Weil pairing into `F_{2^{mk}}^*` | the whole method; needs embedding degree `k` small enough that the finite-field DLP beats rho's `1/2` | embedding degree of the largest prime factor of `#E`, `m ≤ 41` | killed |
| 2 | Fourier diagonalization of a Frobenius-circulant relation matrix | linear algebra only, by a factor `m` for sparse Wiedemann (`m` systems of size `N/m` cost `N^2 w / m`) | counting, not a pipeline run | not an advance |
| 3 | Silverman xedni: lift the instance to a curve over `Q` whose Mordell–Weil rank makes the logarithm readable | the whole method; the dream is polynomial time if a random lift has rank 1 | not run | untested here |

### 1. Embedding degree

`r` is the largest prime factor of `#E(F_{2^m})`. `k` is the order of `2^m`
modulo `r`, which is the extension degree a pairing transfer lands in.
Supersingular curves have `k ≤ 6`. This curve is ordinary (`s_1 = -1`).

Recompute with `python3 research/notes/index-calculus/abliterated_tricks_20261002/embedding_degree.py`.

| m | bits of r | k | log2(k) |
| ---: | ---: | ---: | ---: |
| 5 | 4 | 2 | 1 |
| 7 | 5 | 4 | 2 |
| 11 | 5 | 1 | 0 |
| 13 | 11 | 22 | 4 |
| 15 | 10 | 25 | 4 |
| 17 | 8 | 7 | 2 |
| 19 | 17 | 492 | 8 |
| 23 | 21 | 91124 | 16 |
| 29 | 14 | 554 | 9 |
| 31 | 21 | 11608 | 13 |
| 41 | 40 | 6704346231 | 32 |

The small-`k` rows are the rows whose subgroup is at most 8 bits. Rho there
is a few dozen operations; a pairing into `F_{2^{mk}}` does not win. From
`m = 19` upward, `k > 6` in every row, and at `m = 41` the pairing lands in
an extension of degree about `2^{32}`. That is the kill. A toy pass at
`m ≤ 17` (a few rows with `k ≤ 7` and a tiny `r`) does not survive to
`m = 31`, where `k = 11608` and `r` is only 21 bits, so rho is the smaller
computation by an enormous margin.

`m = 31` is the exploratory Koblitz degree in this repository. It is not
the `m = 83` gate, and this check does not need that gate: the transfer is
already worse than rho at `m = 31`.

### 2. Frobenius-circulant linear algebra

If the factor base is a union of Frobenius orbits, the relation matrix is
block-circulant of order `m`. Over a field whose characteristic does not
divide `m` (the relation modulus is the odd prime `r`), the group algebra
splits and one `N × N` sparse solve becomes `m` solves of size `N/m`.
Sparse Wiedemann then costs about `1/m` of the unsplit solve.

That factor is `m`, a constant in the security parameter relative to the
exponential. It does not change the exponent of relation collection, and
rho already spends the same Frobenius action on its `sqrt(1/(2m))` factor.
The trace itself notes that collection is untouched. No pipeline was run,
so there is no `S` and no speedup. A toy matrix that diagonalizes at
`m ≤ 15` would still leave the exponent alone at `m = 31`.

### 3. Xedni

The trace proposes lifting the two points to a curve over `Q` and reading
the logarithm off Mordell–Weil, and then says a random lift does not have
the rank that would make this polynomial. This note does not run an xedni
solver and does not treat that argument as a measurement. The idea stays
untested here. It is not an index-calculus variant and it has no `S`.

## What the model did not do

It did not name a decomposition, a factor base, or a descent whose
relation yield beats the counting ceiling, and it did not name an exponent
below `1/2` that survives the three checks above. The empty final message
is why this note quotes the trace's conclusions and then checks them,
instead of quoting a finished answer.

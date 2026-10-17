# Protocol (lead 3 of the prime-degree survey): `decompose_by_e3` against the pair-table oracle, counted

Registered 2026-10-15, before the run below, as an addendum to this study.

**Lead 2 needs no new run.** A weighted or block grading of the
symmetrised model was measured on 2026-09-12
(`research/notes/ecc2k130/RESEARCH_KOBLITZ_INDEX_CALCULUS.md`, "The block
structure, measured"). The systems are multilinear in their blocks, so a
per-block degree bound either reproduces the plain Macaulay matrix
(1.00×) or enlarges it (0.21× at `n = 15`). The working technique is
multihomogeneous F5 with Jacobian syzygies; it is a different algorithm and
is not built.

**Lead 3.** `koblitz_sym_s4::decompose_by_e3`, the exact `m = 3`
decomposition by enumerating `e₃ ∈ W₃`, now counts its `F_{2^n}`
multiplications (setup, columns, right-hand sides, splitting tests) and
its word XORs. `examples/e3_vs_table.rs` runs it on uniform target
abscissae and compares against the root-finding pair-table oracle's
`5|F| − 3` multiplications per target, `|F| = 2^{ℓ−1}`
(`research/pdp_root_finder_20261012` §0, the same unit).

Setup:
- Koblitz curve with `b = 1`;
- `V = ⟨1, z, …, z^{ℓ−1}⟩`;
- cells `(n, ℓ)` = (13, 4), (13, 5), (17, 4), (17, 5), (19, 5), (19, 6),
  (23, 5), (23, 6);
- 16 targets each;
- `cargo run --release --example e3_vs_table -- 16 > raw/e3_vs_table.jsonl`.

**Predictions.**
- E1: e₃ ÷ table is at least 100× at every cell.
- E2: at fixed `n`, the ratio grows by at least 4× per unit of `ℓ`
  (`2^{3ℓ−2}` against `2^{ℓ−1}`).

Decision: if both hold, lead 3 is a sealed negative at stage level, and no
ecbench engine is built for it.

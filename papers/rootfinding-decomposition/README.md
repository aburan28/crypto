# Point Decomposition by Root Finding: Replacing Gröbner-Basis Solvers in Elliptic-Curve Index Calculus

`rootfinding_decomposition.tex` writes up the replacement of the F4/F5
Gröbner-basis decomposition solve in this repository's elliptic-curve index
calculus by root finding: direct subtraction over prime fields (the roots of
Semaev's `S₃(x_R, x_i, X)` are `x(R ∓ F_i)`, so no polynomial system is
solved) and a meet-in-the-middle pair-sum table over binary fields. It is the
typeset form of `../../RESEARCH_SEMAEV_ROOTFIND_VS_GROEBNER.md`.

**Status**: research draft, anonymized, not peer-reviewed, not submitted, not
posted, **not cleared for either** until the items under "Before submission"
are closed. Restore authors with `make deanonymize` after creating an
`AUTHORS` file (see `../AUTHORS.example`).

**Change class** (per `AGENTS.md`): *engineering*. Every ratio in the paper is
a stage measurement of the per-target decomposition oracle, and the paper
says so in the abstract, §6 and §7. It is not reported as an end-to-end
speedup; the one end-to-end comparison at small size (§6.1) favours F4.

## What it says

- Proposition 1 (with proof): for `y² = x³ + ax + b` and affine `R`, `F`
  with `x_R ≠ x_F`, the roots of `X ↦ S₃(x_R, x_F, X)` are exactly
  `x(R + F)` and `x(R − F)`. Hence the `m = 2` decomposition system has a
  parametrized variety and needs no Gröbner basis.
- Algorithm 1 (prime, direct subtraction with batched inversion, `Θ(|F|)`
  field operations and `2|F|` lookups) and Algorithm 2 (binary, pair-sum
  table, `1 / |F| / |F|²` lookups for `m = 2 / 3 / 4`), both exact and
  complete for the stated `m`.
- A fail-closed verification protocol (group re-addition, `S₃` vanishing,
  Gröbner-variety membership, agreement on refutation, positive controls)
  under which every reported number was produced.
- Measurements against the codebase's own matrix-F4 and Buchberger solvers on
  identical targets: Koblitz `F_2^9` / `F_2^13` (383× / 1535× per
  decomposition), small prime curves (3–5 × 10⁴×), prime fields 40–83 bits,
  and NIST P-192/224/256/384/521 with standard parameters (3021–5802× per
  target against Buchberger on the same ideal in the same field; P-256 also on
  the dedicated field, 11–12× over generic arithmetic).

Every number cites a committed artefact (`crypto:research/semaev_*_20261012/`)
in its table caption.

## Lean formalization (Proposition 1)

`lean/S3Roots.lean` proves, over an arbitrary field with the two curve
equations as hypotheses, the factorization
`S₃(xR, xF, X) = (xR − xF)² (X − x₊)(X − x₋)` (`SemaevS3.S3_factor`) with
`x± = λ±² − xR − xF`, `λ± = (±yF − yR)/(xF − xR)` the affine chord formulas,
and derives `S₃(xR, xF, X) = 0 ↔ X = x₊ ∨ X = x₋` (`SemaevS3.S3_root_iff`).
No `sorry`; `#print axioms` gives only `propext`, `Classical.choice`,
`Quot.sound`. Toolchain Lean 4.34.1, Mathlib `v4.34.1` pinned in
`lean/lake-manifest.json`. Check with

```
cd lean && lake exe cache get && lake build
```

(the Mathlib cache is ≈7.6 GB; `.lake/` is git-ignored).

## What it does not say

- No ECDLP, key-recovery, sub-rho or deployed-curve-security claim. The NIST
  measurements are of a subroutine on a factor base of 6 or 64 points; the
  relation-collection and linear-algebra stages remain exponentially out of
  reach and the asymptotics of index calculus are unchanged (§7).
- No change to the Gaudry–Diem verdict that prime-field index calculus does
  not beat Pollard rho.
- No change to the exponential Frobenius-invariant factor base on Koblitz
  curves (`|V| = 2^{ord₂(n)}`); the binary measurements are on bounded
  factor bases and subspaces by design.
- No claim that the meet-in-the-middle ordering persists at large `m`: its
  cost is `Θ(|F|^{m−1})`; all measurements are at `m = 2` (binary table to
  `m ≤ 4`).
- No claim against state-of-the-art F4/F5 implementations: the baselines are
  this codebase's own, and where the Buchberger baseline capped its step
  budget the ratios understate its cost.

## Building

```
make           # PDF (pdflatex, no bibliography; thebibliography inline)
make check     # unresolved refs, LaTeX errors, overfull boxes, page count
make clean
```

The rules are the shared `../paper.mk`. The committed PDF was built with
`tectonic -X compile rootfinding_decomposition.tex` (same source, same
packages); it has no overfull boxes and all references resolve. ≈10 pages.

## Reproduction

- Code: `src/cryptanalysis/prime_fast.rs` (Algorithm 1),
  `src/cryptanalysis/koblitz_index_calculus.rs` (`PairSumTable`, Algorithm 2;
  `DecompositionStrategy::PairTable` is the default).
- Benchmarks: `examples/semaev_root_vs_groebner_bench.rs` (§6.1–6.2),
  `examples/semaev_rootfind_scaling_bench.rs` (§6.3),
  `examples/semaev_rootfind_nist_bench.rs` (§6.4; knobs `F4_FB`,
  `F4_CONSTRUCTED`, `F4_RANDOM`, `SCALE_FB`, `SCALE_RANDOM`).
- Artefacts: `research/semaev_root_vs_groebner_20261012/`,
  `research/semaev_rootfind_scaling_20261012/`,
  `research/semaev_rootfind_nist_20261012/` (console transcripts + JSON).
- `cargo test --release --lib`: 2247 passed at the measured commit.
- Lean: `lean/S3Roots.lean` (see above).

## Before submission

Owed, in order:

1. **Novelty audit.** The `m = 2` observation (Proposition 1) is folklore in
   spirit; the paper claims only the concrete statement, the implementation,
   the protocol and the measurements. Before any posting, a literature pass
   must confirm how Amadori–Pintore–Sala (FFA 2018) and any later
   root-finding-vs-Gröbner work relate, and the Related Work paragraph must be
   corrected if the characterization there is off.
2. **Lean formalization of Proposition 1** — done, see below; what remains
   owed is tying `slope`/`addX` to Mathlib's `WeierstrassCurve.Affine` group
   law by name rather than by the (identical) chord formulas, so the
   statement is literally about `R + F` in Mathlib's `Point` type.
3. **Stronger F4 baseline.** Re-run §6.2 and §6.4 against an external
   state-of-the-art F4/F5 (e.g. Magma or msolve) on the same ideal, so the
   "not faster than Θ(|F|)" argument in §7 is backed by a measurement rather
   than an argument.
4. **Per-curve fields for P-192/224/384/521.** §6.4-C is P-256 only; the
   other four run on generic arithmetic. Instantiate their Solinas fields and
   fill the column.
5. **Timing hygiene.** Re-measure Tables 1–6 on an idle machine with repeat
   counts and report spread; the P-521 control row shows run-to-run jitter
   (3.66 ms vs 1.13 ms) that a single run cannot resolve.
6. **Independent replay** of every table by someone other than the author of
   the benchmarks, per the repository's ledger rules.

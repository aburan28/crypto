# Index-calculus technique integration, 2026-10-12

**Scope.** The nine technique gaps named on 2026-10-12 for the ECC2K-130 /
Koblitz index-calculus flow, each mapped to the module that carries it on
`main`, what was missing, what this round built, and where the evidence is.
Every row below is a toy-ladder measurement on public synthetic or
known-answer instances (`n ≤ 41` Koblitz rungs, `K_0: y² + xy = x³ + 1`);
none is a key-recovery, sub-rho, or deployed-curve-security claim, and the
total-work verdict of `RESEARCH_ECC2K130_IC_FEASIBILITY.md` §3 is not
changed by anything here unless a row says so explicitly.

**Inventory correction.** An earlier inventory in this program was taken
against a stale checkout and reported the large-prime variation, the
filtering/sparse-LA stage, the `Z[τ]` closure and the symmetric-S₄ encoding as
absent. On `origin/main` they were already present (`koblitz_large_prime`,
`large_prime_filter`, `ecbench_large_prime`, `koblitz_sparse_la`,
`compact-orbit-scan:norm=N`, `koblitz_sym_s4`); the real gaps were narrower
and are the ones built here.

## Status table

| # | Technique | Carried by (on `main` before this round) | Gap on 2026-10-12 | Built this round | Evidence |
|---|---|---|---|---|---|
| 1 | Double-large-prime relation graph | `koblitz_large_prime` (1- and 2-LP exact elimination, wide), `large_prime_filter` (hypergraph peel + Markowitz two-core), `ecbench_large_prime`, `gaudry_cubic::LargePrimeEliminator`; lead-ladder §7.4 closed the *single*-LP table variant by analysis | `koblitz_large_prime` had no caller in the Koblitz pipeline; no measured double-LP graph on the Frobenius-folded pair table | see §1 | `research/ic_techniques_20261012/double_large_prime/` |
| 2 | Joux–Vitse F4 trace reuse (G1) | nothing (`LEAD_LADDER` §L7: "not implemented") | all of it | see §2 | `research/ic_techniques_20261012/g1_f4_trace_replay/` |
| 3 | Sparse Macaulay + Krylov refutation (G3) | `sparse_macaulay` (structured sparse elimination; argues a Krylov method gives no reduced basis) | the BooleanSolve consistency test itself (1 ∈ rowspace via Wiedemann, no basis needed) and the n=41 falsifier numbers | see §3 | `research/ic_techniques_20261012/g3_sparse_macaulay_wiedemann/` |
| 4 | Faugère–Svartz equivariant F4 (G2) | `examples/g2_equivariant_orbit_probe.rs`, `examples/g2_wide_probe.rs`: block decomposition exact, 31–128× rank-time ratio vs dense symmetric (`RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md` §3) | steps 4–5 of the scope: inverse transform, `reduce_system` engine, agreement vs the gauge-fixed baseline | see §4 | `research/ic_techniques_20261012/g2_equivariant_f4/` |
| 5 | Shantz–Teske δ method | dimension is a parameter everywhere; nothing sweeps δ and reports the net | the sweep | see §5 | `sage/factor_base_dor.sage --delta-sweep`, `research/ic_techniques_20261012/a2_cyclotomic_and_delta/` |
| 6 | Galbraith–Gebregiyorgis special subspace, direct symmetric S₄ at m=3 | `koblitz_sym_s4` (quadratic symmetric S₄ with geometric `W₁,W₂,W₃`; agrees with chain and brute force) | the D\* measurement of the symmetric encoding against the chain, and the special-vs-generic subspace comparison | see §6 | `sage/factor_base_dor.sage --symmetric-s4`, `sage/factor_base_dor_results.json` |
| 7 | Higher-arity extraction (6-sum) | Boolean 6-sum systems (`wide_sixsum`, `koblitz_rotated_chain`, `native_signed_mitm`), the 4-sum pair-table kernel; feasibility note §4.7 calls the extraction kernel "the concrete open design" | the extraction kernel and its probe accounting | see §7 | `research/ic_techniques_20261012/sixsum_extraction/` |
| 8a | `Z[τ]` small-norm multiplier closure (A1) | `compact-orbit-scan:norm=N`, measured: 9.2× one-target total at n=41, m=3 (`research/ic_leads_20261010/RESULTS_ROUND2.md` §A1) | none — landed and measured | — | `RESULTS_ROUND2.md` |
| 8b | Cyclotomic-sparse bases (A2) | nothing (`cyclotomic_cosets` is a different object) | all of it | see §8 | `research/ic_techniques_20261012/a2_cyclotomic_and_delta/` |
| 9 | Filtering + sparse LA for a 129-bit modulus | `koblitz_sparse_la` (dedup, singleton, clique cascade, merge, block Wiedemann, block Lanczos) with `MAX_MODULUS_BITS = 63` | the wide-modulus path | see §9 | `research/ic_techniques_20261012/wide_sparse_la/` |

Sections 1–9 below are filled from the per-item evidence directories.

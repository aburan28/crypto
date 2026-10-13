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
| 2 | Joux–Vitse F4 trace reuse (G1) | nothing (`LEAD_LADDER` §L7: "not implemented") | all of it | **WIP** (uncommitted worktree `ict-g1-trace-20261012`) | — |
| 3 | Sparse Macaulay + Krylov refutation (G3) | `sparse_macaulay` (structured sparse elimination; argues a Krylov method gives no reduced basis) | the BooleanSolve consistency test itself (1 ∈ rowspace via Wiedemann, no basis needed) and the n=41 falsifier numbers | **WIP** (uncommitted worktree `ict-g3-wiedemann-20261012`) | — |
| 4 | Faugère–Svartz equivariant F4 (G2) | `examples/g2_equivariant_orbit_probe.rs`, `examples/g2_wide_probe.rs`: block decomposition exact, 31–128× rank-time ratio vs dense symmetric (`RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md` §3) | steps 4–5 of the scope: inverse transform, `reduce_system` engine, agreement vs the gauge-fixed baseline | **WIP** (uncommitted worktree `ict-g2-equivariant-20261012`) | — |
| 5 | Shantz–Teske δ method | dimension is a parameter everywhere; nothing sweeps δ and reports the net | the sweep | see §5 | `sage/factor_base_dor.sage --delta-sweep`, `research/ic_techniques_20261012/a2_cyclotomic_and_delta/` |
| 6 | Galbraith–Gebregiyorgis special subspace, direct symmetric S₄ at m=3 | `koblitz_sym_s4` (quadratic symmetric S₄ with geometric `W₁,W₂,W₃`; agrees with chain and brute force) | the D\* measurement of the symmetric encoding against the chain, and the special-vs-generic subspace comparison | see §6 | `sage/factor_base_dor.sage --symmetric-s4`, `sage/factor_base_dor_results.json` |
| 7 | Higher-arity extraction (6-sum) | Boolean 6-sum systems (`wide_sixsum`, `koblitz_rotated_chain`, `native_signed_mitm`), the 4-sum pair-table kernel; feasibility note §4.7 calls the extraction kernel "the concrete open design" | the extraction kernel and its probe accounting | see §7 | `research/ic_techniques_20261012/sixsum_extraction/` |
| 8a | `Z[τ]` small-norm multiplier closure (A1) | `compact-orbit-scan:norm=N`, measured: 9.2× one-target total at n=41, m=3 (`research/ic_leads_20261010/RESULTS_ROUND2.md` §A1) | none — landed and measured | — | `RESULTS_ROUND2.md` |
| 8b | Cyclotomic-sparse bases (A2) | nothing (`cyclotomic_cosets` is a different object) | all of it | see §8 | `research/ic_techniques_20261012/a2_cyclotomic_and_delta/` |
| 9 | Filtering + sparse LA for a 129-bit modulus | `koblitz_sparse_la` (dedup, singleton, clique cascade, merge, block Wiedemann, block Lanczos) with `MAX_MODULUS_BITS = 63` | the wide-modulus path | see §9 | `research/ic_techniques_20261012/wide_sparse_la/` |

Sections 1–9 below are filled from the per-item evidence directories.

## 1. Double large prime on the Frobenius-folded pair table

**Built.** `ic_framework/large_prime.rs` (new) plus `ic_boundary::collect_and_solve_relations`,
`FrobeniusPairTable::build_filtered` / `orbit_key`, and `ic bench --large-primes {0,1,2}
--lp-domain-factor c`: partial relations `R − s₁ − s₂ = L` (one LP) and the two-LP
form are fed to `koblitz_large_prime::WideLargePrimeEliminator`, every complete row
is replayed in the counted group before it enters linear algebra, and LP=0 is
counter-identical to the classic loop on a fixed seed (tested). Tests: 3 new
(`ic_framework::large_prime`), `koblitz_large_prime` 4/4, `ic_framework` 44/44,
`ic_boundary` 35/35.

**Measured** (one-target cold totals, 51 cells / 260 targets, all verified, 0 rejected
partials): at n=23 and n=31 the rungs are saturated (`|F|³/6 ≥ r`, the plain table
decomposes nearly every trial) and every LP cell loses 2.2–750×, as lead-ladder §7.4
predicts above its break-even. At n=41 (r = 2³⁹, 10 columns, below break-even) single
LP cuts trials per relation 4.1–10.3× and double LP 8.7–138×; end to end the best
double-LP cell (|D|/|F| = 10) is 23.8× cheaper than the plain oracle at fixed columns
(2.21·10⁶ vs 5.25·10⁷ GAE) but only **2.1×** cheaper than the equal-envelope control
(plain oracle on a base the size of F ∪ D, 4.59·10⁶ GAE), with an 11-row instead of
101-row matrix. Single LP is 1.07–1.98× *worse* than that control.

**Verdict.** The double-LP graph turns §7.4's "a bounded D cannot cross over" into a
2× constant over spending the same table on more columns. Engineering, one rung, no
exponent; the best cell is still 17.7× the counted rho. n=41 at 5 and 20 columns was
not run (host contention; `run.sh 41` reproduces).

## 7. Six-sum extraction kernel

**Built.** `koblitz_sixsum_extract.rs` (new, framework oracle `sixsum-extract`) with a
materialised Frobenius pair table and five kernels: the 4-sum S₃ probe (baseline),
A = iterate one table (`R − s₁`, then the 4-sum scan, early abort), B = S₄ quartic per
ordered state pair in row order, B2 = the same in diagonal stride order, C = 5-sum
control (`R − P_e`, then the 4-sum scan); S₄ quartic and two root finders checked
against `binary_semaev_s4` and brute force over all of F₂¹⁷; every relation lifted and
re-added in the group with orbit-column coefficients (fail-closed; 0 rejects, 0 lift
failures across 47 cells). Example `koblitz_sixsum_extract_probe`. Tests 6/6.

**Measured** at n=41 (r ≈ 2³⁹), per verified relation, S₃ solves for 4sum/A/C and S₄
solves for B/B2 (wall is secondary; host contended):

| K | 4-sum S₃/rel | A S₃/rel | C S₃/rel | B2 S₄/rel | A outer probes/rel | feasibility §4.7 model |
|--:|--:|--:|--:|--:|--:|--:|
| 5 | 4.55·10⁷ | 1.37·10⁷ | 0 rel in 1.4·10⁸ | 1.12·10⁶ | 332 | 1.40·10⁴ |
| 10 | 6.19·10⁶ | 7.55·10⁶ | 7.84·10⁶ | 2.66·10⁵ | 46 | 875 |
| 20 | 1.01·10⁶ | 7.91·10⁵ | 8.55·10⁵ | 1.32·10⁵ | 1.97 | 54.7 |

**Verdict.** The feasibility note's `45·r/(n⁴K⁴)` was in *first-summand-probe* units
(kernel A's outer probes per relation equal the 4-sum's targets per relation), and
each such probe hides `K²n²` S₃ solves; in S₃-solve units A ≡ 4-sum ≡ C, because every
inner loop is the 4-sum scan of a shifted pseudo-random target and the 6-sum
abundance is wasted by the scan order. Only the stride-ordered S₄ kernel B2 converts
multiplicity into fewer solves (`r/(8A)` quartics, 7.6–40× fewer than the 4-sum), but
one S₄ root solve cost 20–60 S₃ solves here, so it is 1–8× *more* work; break-even
needs a quartic root solve under ~12 S₃. Classification: accounting correction to
§4.7 plus an unestablished engineering candidate (B2). n=23/31 are saturated (all 30
cells within 0.75–1.3× of the dense model). Several n=41 cells have < 32 relations;
the K=5 5-sum control found 0 relations where ~7 were expected (unresolved, ~0.1%
outcome).

## 5. Shantz–Teske δ method

**Built.** Sage `--delta-sweep` in `factor_base_dor.sage` (work-per-relation proxy
from Macaulay shape at D\*); Rust `delta_method.rs` + `examples/delta_method_sweep.rs`
on K₀ at n ∈ {17,23,31}.

**Measured.** Sage at n=8,9: δ=+1..+2 lowers D\* ~1 and the work-per-relation proxy
10–13×; δ=+3 kills coverage (`sage/factor_base_dor_results.json`). Rust at n=31, m=2:
measured fraction tracks the predicted envelope (δ=+1: 0.219 vs 0.223; δ=+3: 0.031 vs
0.061); F4 max deg stays 3; net ms/relation is lowest at δ=−1 among the completed
cells (213 ms) and rises as coverage collapses at δ=+3. m=3 cells are far more
expensive (net 5·10⁴–∞ ms/relation).

**Verdict.** A constant-factor tuning knob at small n; no exponent change. Full
README packaging for the Rust sweep is still thin (raw JSON + generated summaries
only).

## 6. Galbraith–Gebregiyorgis special subspace / direct symmetric S₄

**Built.** `sage/factor_base_dor.sage --symmetric-s4`: product spans W₁,W₂,W₃,
eleven-term S₄ encoding, PolyBoRi Boolean UNSAT, D\* via dense GF(2) Macaulay.

**Measured.** Underdetermined GG cells (6ℓ−3 > n) never UNSAT (spurious Boolean
solutions). Determined cells: n=9 ℓ=2 → D\*=4; n=15 ℓ=3 → D\*≥5; subfield base at
the small cell reaches D\*=3 with zero spurious. Unknown count drops vs the chained
S₃ system (15 vs 24 at n=9 ℓ=2) but the degree does not.

**Verdict.** Buys variables, not degree. Not a D\* drop for ECC2K-130 total work.

## 8. Cyclotomic-sparse bases (A2)

**Built.** `cyclotomic_factor_base.rs` (Frobenius-stable γ-weight sets, ONB-weight
matched controls), example `a2_cyclotomic_base`, plugin wiring.

**Measured** (census + partial solver cells at n=11,23; see
`research/ic_techniques_20261012/a2_cyclotomic_and_delta/`): bases are Frobenius-stable
at predicted size; K₁/n=23/w=1 has usable subgroup points, K₀/n=23/w=1 does not
(0 in subgroup). Several cyclotomic solver cells timed out or had null encodings in
the first n=23 pass; ONB-matched controls reached F4 max deg 3 with multi-second
walls. Yield table in the generated summary is still empty for the main n=23 file.

**Verdict.** Implementation + census landed; the falsifier ("cyclotomic solving
degree lower or yield/column higher than matched ONB") is **not yet established** —
treat as open engineering lead, not a closed negative or positive.

## 9. Wide-modulus sparse LA

**Built.** `koblitz_sparse_la` generic modular word: `u64` fast path unchanged;
`U128` / `U192` / `U256` cover odd moduli up to 255 bits (`MAX_MODULUS_BITS =
4·64−1`), including ECC2K-130's 129-bit r. Dispatch in `koblitz_index_calculus`
picks the sparse path for wide r.

**Measured.** Module tests for planted systems on the wide paths (agent commit);
dedicated `research/ic_techniques_20261012/wide_sparse_la/{results.json,README.md}`
benchmark pack is **not yet in this PR**.

**Verdict.** Engineering prerequisite for 130-bit sparse LA is in-tree; 10⁹-column
cost remains an extrapolation until the benchmark pack lands.

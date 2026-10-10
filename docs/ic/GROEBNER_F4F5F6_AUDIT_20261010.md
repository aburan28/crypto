# Audit of the F4, F5, F6-IC and Macaulay engines, 2026-10-10

Scope: every Gröbner-basis and Macaulay-matrix engine the index-calculus
pipeline can call, read for correctness first and for speed second, with
the fixes and measurements this audit landed. Paths are under
`src/cryptanalysis/`. Wall times are from one shared 4-core cloud host
(AVX-512) and are indicative only; counters are exact.

## 1. What runs where

| engine | file | ring | linear algebra | production use |
|:--|:--|:--|:--|:--|
| Buchberger / naive matrix-F4 | `groebner_f4.rs` | `F_p`, BigUint | dense, per pair | reference only (silent 5 000-step cap) |
| degree-bounded F4 | `f4_fp.rs` | `F_p`, p < 2³² | sparse-by-reducers, then dense RREF on free columns | coordinate thread, PKM pilots, GLV/JV callers |
| tower F4 | `f4_fp_tower.rs` | `F_p[y,u]/(towers)` | Faugère–Lachartre (`A⁻¹B`, `D − C·B'`) | PKM tower oracle |
| signature tower F4 | `sig_fp_tower.rs` | same | same | PKM tower oracle (variants) |
| Buchberger | `pq_groebner_f2.rs` | `F_2[v]/(v²−v)` | term lists | reference |
| boolean F4 | `pq_f4_f2.rs` | boolean, ≤ 64 vars | own M4RI kernels | PQ descent studies |
| Macaulay matrix-F4 | `koblitz_groebner.rs` (`matrix_f4_f2*`) | boolean, ≤ 64 vars | `gf2_elim` (≥ 128 × 256) or legacy | from-scratch reference engine |
| matrix-F5 criterion | `matrix_f5_f2.rs` | boolean | `gf2_elim` | `SolverEngine::MatrixF5` |
| inherited F4 | `inherited_f4.rs` | boolean | root echelon + per-node re-reduction | **default decomposition oracle** |
| F6-IC | `koblitz_index_calculus.rs` (`groebner_decompose_f6_ic`) | inherited F4 + exact geometry gate | as above + batched group additions | opt-in arm |
| wide solver | `wide_groebner.rs` | boolean, ≤ 128 vars | `gf2_elim` at every node | m ≥ 4 ladder |
| sparse Macaulay | `sparse_macaulay.rs` | boolean | structured sparse + dense finish | profiling |
| GF(2) kernel | `gf2_elim.rs` | — | Four Russians, 4 × 8-bit tables, AVX-512 | everything above |

The engines agree with each other and with brute force on every pinned
test (`cargo test --release --lib` on the modules above: 88 + 15 + 9 + 14
+ 13 + 12 + 21 + 6 passing). The audit found one soundness bug, two
contract defects and several latent holes, listed next.

## 2. Correctness findings

### 2.1 Fixed: `f4_fp::solve` certified solutions of an inconsistent system

A degree-bounded run (`pairs_above_bound > 0`) or a staircase stop returns
the basis of a **sub-ideal**, and `add` retires an input element whose
leading monomial becomes divisible and then drops its pair when the lcm
is above the bound (`f4_fp.rs`, `add`), so the retired element leaves
the basis unreduced. The substitution tree then enumerates a superset of
the variety and `solve` reported it as `Verdict::Solutions`.

Counterexample, now the test `truncated_run_must_not_certify_solutions`:
over `F_5` with `D = 2`, `{x³ + 1, xy − 1, x + y}` is inconsistent, but
the degree-2 step finds `x² + 1`, retires `x³ + 1`, drops the degree-3
pair, and `solve` returned `[2, 3]` and `[3, 2]`. Callers in
`jv_quintic.rs`, `jv_quartic.rs` and `glv_gaudry.rs` consumed
`Solutions` without re-evaluating the input.

Fix: `solve` evaluates every candidate on the input and drops the ones
that fail (`SolveReport::candidates_rejected`); an empty remainder is
`Inconsistent`. This is exact: every basis in the tree generates a
sub-ideal of its input, so the tree's candidates contain the variety,
and refutations need no check because `1` in a sub-ideal is `1` in the
ideal. Verdicts of runs that were never truncated are unchanged.

### 2.2 Open contract defects (documented, not fixed here)

- `pq_f4_f2.rs`: the docstring says a budget or size overrun returns a
  set generating the same ideal. It returns only the *active* elements;
  an element deactivated by a newer divisor is recoverable only through
  the unprocessed pair with it, so the timed-out set can generate a
  strictly smaller ideal. `timed_out` is set, so no verdict is wrong,
  but no caller may treat that set as generators. The echelon-timeout
  path also fails to put the drained pairs back, so `pairs_left`
  under-reports.
- `groebner_f4::buchberger`: a silent `step_cap = 5_000` returns a
  partial basis with no flag; LIFO pair order; product criterion only.
  Reference for tiny inputs only.
- `matrix_f5_f2::F5Criterion::new` is sound only when
  `occurring_vars(polys) ⊆ multiplier_mask`; the one production caller
  satisfies it, nothing asserts it, and a violation would silently drop
  rows (`pack(...)? → continue`).
- `matrix_f5_f2` certificate output forms (`CertifiedOriginalRows`,
  `SelectedColumnCertificate`, `SupportSeparatedOriginalRows`) return
  unreduced rows; the "leading column ≥ degree-≤1 boundary ⇒ linear
  consequence" invariant does not hold for them. Safe for rank
  consumers only.
- `inherited_f4::ReducedBasis::from_system` with `degree` below a
  generator's degree would skip completion for that generator after a
  degree drop (the ladder never does this; library callers could).

### 2.3 Verified sound

- Inherited specialisation (`inherited_f4.rs`): kept rows keep distinct
  pivots under `v ↦ c`, folding collisions are XOR-scattered and zero
  images detected, degree drops are completed, and the **linear tail**
  equals the from-scratch tail under the pinned policy. Under the
  production policy (`support_local`, `linear_elimination`) the row
  space is a documented subspace: roots are preserved, the tree differs.
- `gf2_elim.rs`: pivot search, table build and `pext` indexing are
  consistent (the tables are subset tables, not Gray-code tables, same
  cost); RREF is checked bit for bit against a naive reference at eight
  shapes and three densities, echelon at three shapes. Not tested: rows
  ≥ 512 against naive (where 8-bit tables arise), the AVX2 path (never
  selected without an env var), garbage padding bits.
- `sparse_macaulay.rs`: the tail-boundary argument is correct; rows are
  dropped only when exactly zero.
- Boolean F4 (`pq_f4_f2.rs`): field pairs, Gebauer–Möller `UPDATE`, the
  "new iff pivot has no active divisor" test and both M4RI kernels check
  out; the F5 criterion's `W_i(d)` construction and its rejection of the
  naive boolean rule are right.
- F6-IC's refutations are exact for the enumerated base: support checks
  fire only on fully assigned, undefined summands; one-fixed closure
  enumerates every lift and every compatible second point; residual
  closure walks every sign lift; every witness is group-verified.
  Pinned against brute force only at n = 9 (bases of 14 and 2^2…2^4
  points); nothing at |F| > 256 or n > 9.

## 3. What F6-IC is

F6-IC = inherited-F4 DPLL + a factor-base domain check + meet-in-the-middle
leaves. Concretely, at every node: (i) a fully assigned summand whose
abscissa is not a base point refutes; (ii) at m = 3 with one summand fixed
and |F| ≤ 256, every compatible second point is added and the third looked
up (Θ(|F|) batched additions a call); (iii) with m − 1 summands fixed the
last is looked up. The inherited specialisation is the genuinely new
engineering (it predates F6-IC); the closure is classical enumeration
plus hash lookup, which is why the ecbench ladder prices it as a
relabelling once geometry is charged
(`research/f6_ic_ecbench_ladder_20261005`: median log₂(W4/W6) 0.38–0.43
at |F| ≤ 139, 0.00 at |F| = 275, S up 1 292× at m = 19). Its exponent is
the DPLL tree's, about |F|^(m−2) live partial assignments
(`docs/ic/perf/OPTIMIZATION_PLAN.md` §7); geometry removes the last one or
two levels only.

## 4. Changes landed and measured

All verdicts, witnesses and the node-tree counters the tests pin are
unchanged; only work per node moves. Each change has a same-binary
control.

| change | where | control | effect |
|:--|:--|:--|:--|
| root echelon of the inherited engine and every other `echelon_f2_counted` caller now reaches `gf2_elim` at ≥ 128 × 256 (it never did: only `rref_f2_counted` routed there) | `koblitz_groebner.rs` | `KIC_F4_KERNEL=legacy` | kernel alone 1.8–3.2× on the oracle's own Macaulay matrices (`gf2_elim_bench --quick`, table below); reach ladder 1.0–1.1× (its matrices are mostly below the gate and build-dominated) |
| child nodes ask the node oracle **before** specialising the parent's bases | `koblitz_groebner.rs::solve_rec` | `KIC_F4_EARLY_ORACLE=0` | n = 17, dim 6, m = 3, 8 targets: F6-IC 197.7 → 185.4 ms with both changes; the same 463 reductions, 280 refutations, 58 341 additions |
| partial support pruning: a summand with some bits assigned refutes as soon as no base code agrees with them | `koblitz_index_calculus.rs::F6GeometricGate::decide` | `KIC_F6_PARTIAL_SUPPORT=0` | n = 17: reductions 508 → 463 (−9 %), 45 extra refutations; larger cells in §4.2 |
| `gf2_elim` selects AVX2 on hosts without AVX-512 (they ran scalar) | `gf2_elim.rs::simd_kind` | `KIC_GF2_SIMD=0` | §4.3; no change on AVX-512 hosts |
| `f4_fp::solve` candidate filter (§2.1) | `f4_fp.rs` | — | correctness |
| `examples/f6_ic_probe.rs`: same-binary F4/F6-IC timing probe with verdict agreement | new | — | tooling |

### 4.1 Kernel (L0), this host

| cell | rows × cols | legacy ms | gf2_elim ms | speedup |
|:--|--:|--:|--:|--:|
| K/2^23 m2 d3 | 529 × 1 464 | 1.6 | 0.9 | 1.75× |
| K/2^23 m2 d4 | 5 842 × 8 449 | 326.3 | 106.3 | 3.07× |
| K/2^5 m3 d4 | 860 × 2 466 | 4.8 | 1.5 | 3.22× |
| K/2^5 m3 d5 | 4 940 × 8 357 | 182.4 | 64.2 | 2.84× |
| rand 4096² ½ | 4 096 × 4 096 | 115.8 | 65.1 | 1.78× |

Reach ladder (`pdp_bench --tier reach --engines groebner`), s a target,
legacy → routed: K₁/2¹⁷ m2 0.0019 → 0.0018; K₀/2²³ m2 0.0137 → 0.0123;
K₁/2⁹ m3 0.0008 → 0.0009; K₁/2¹⁵ m3 0.0084 → 0.0080. Hits identical.

### 4.2 F6-IC probe

Standard-subspace bases, m = 3, `InheritedF4 { max_degree: 3 }`, node
budget 20 000, best of 3 runs a target (one run at n = 31), both arms in
one binary; "controls off" is F6-IC with `KIC_F4_EARLY_ORACLE=0
KIC_F6_PARTIAL_SUPPORT=0`. Verdicts agree across arms on every target.

| cell | \|F\| | arm | total ms | reductions | splits | gate refutations | group additions | word XORs |
|:--|--:|:--|--:|--:|--:|--:|--:|--:|
| K₁/2¹⁷ dim 6, 8 targets (1 hit) | 62 | inherited F4 | 255.6 | 1 201 | 504 | — | — | 15.02 M |
| | | F6-IC, controls off | 197.7 | 508 | 275 | 235 | 58 341 | 12.65 M |
| | | **F6-IC** | **185.4** | **463** | 275 | 280 | 58 341 | 12.64 M |
| K₁/2²³ dim 7, 4 targets (0 hits) | 107 | inherited F4 | 578.1 | 1 180 | 500 | — | — | 53.27 M |
| | | F6-IC, controls off | 531.0 | 542 | 288 | 216 | 91 592 | 51.97 M |
| | | **F6-IC** | **520.3** | **450** | 272 | 276 | 91 592 | 51.92 M |
| K₀/2²³ dim 8, 4 targets (0 hits), cap off | 275 | inherited F4 | 1 741.2 | 2 796 | 1 192 | — | — | 158.51 M |
| | | F6-IC, controls off | 1 789.6 | 2 796 | 1 192 | 0 | 0 | 158.51 M |
| | | **F6-IC** | **1 719.7** | **2 692** | 1 184 | 88 | 0 | 158.47 M |
| K₀/2³¹ dim 9, 1 target, budget 2 000 (exhausted) | 553 | inherited F4 | 2 253.3 | 1 410 | 602 | — | — | 281.58 M |
| | | F6-IC | 2 375.4 | 1 357 | 599 | 47 | 0 | 281.54 M |

Reading it: the early oracle and partial support pruning take 6 % (n = 17),
2 % (n = 23, dim 7) and 4 % (n = 23, dim 8) off F6-IC's wall and 9 %, 17 %
and 4 % off its reductions, with identical geometry. Above the 256-point
closure cap F6-IC was identical to inherited F4; the partial support test
is now the only thing it does there, and it is worth about 4 %. None of
this moves the exponent: at every cell the random targets are refuted by
walking the tree, and the n = 31 cell exhausts a 2 000-node budget in
2.3 s a target where meet in the middle answers in under a millisecond
(`OPTIMIZATION_PLAN.md` §6).

### 4.3 Kernel row update on this host

`gf2_elim_bench --quick`, ms per cell, two runs (3 and 5 repetitions):

| cell | scalar (`KIC_GF2_SIMD=0`) | AVX2 (`KIC_GF2_FORCE_AVX2=1`) | AVX-512 (default) |
|:--|--:|--:|--:|
| K/2^23 m2 d4 (5 842 × 8 449) | 79.8 / 72.2 | 75.4 / 79.0 | 110.8 / 99.7 |
| K/2^5 m3 d5 (4 940 × 8 357) | 48.2 / 55.9 | 56.7 / 44.7 | 58.3 / 46.8 |
| rand 4096² ½ | 57.8 / 56.1 | 33.9 / 34.0 | 70.9 / 36.6 |

On this cloud Xeon the AVX-512 path is never the fastest and is the
slowest on the largest Macaulay cell, by 1.3–1.4× against AVX2; on the
random matrix AVX2 and AVX-512 tie and beat scalar 1.6×. The default is
left as the repository measured it on its own runner
(`gf2-elim-reference-v1.json`); what this audit changed is that hosts
without AVX-512 now take the AVX2 path instead of scalar. Re-measuring
the preference on the CI runner is a one-line change in `simd_kind`
with these numbers as the reason.

## 5. The fastest version of this, in order of expected payoff

Engineering (keeps every verdict):

1. **Inherited engine, specialisation (92 % of build at n = 17).**
   Word-wise `rewrite_with` (`pext`/`pdep` compress and merge instead
   of bit scatter, 3–5× on that part), batch the displaced rows against
   the touched pivots instead of one `insert` at a time, AVX-512
   `xor_words`, fold forced-propagation chains into one layout step.
2. **Boolean F4 (`pq_f4_f2`)**: a Faugère–Lachartre split (reducers to
   RREF once, dense kernel on `D − C·A⁻¹B` over non-pivot columns only;
   3–5× fewer dense words at the n = 24, d = 4 shape), a divisor index
   for symbolic preprocessing (`reducer_among` is a linear scan per
   monomial), sparse rows until the dense block, and `gf2_elim`'s kernel
   in place of the two home-grown ones.
3. **`f4_fp`**: the default narrow RREF (`rref32_deferred`) is serial;
   row-blocked parallel deferred elimination as the tower engine does;
   precompute `A⁻¹B` instead of reducing every pair row through the
   reducer cascade; Gebauer–Möller chain criterion and sugar; packed
   monomial keys and a divisor index; per-task `FIELD_OPS` accumulation.
4. **`gf2_elim`**: select AVX2 by default where AVX-512 is absent (today
   those hosts run scalar table XORs), hoist the SIMD dispatch out of the
   per-row closure, 6–8 narrower tables per pass to halve matrix passes,
   test rows ≥ 512 against the naive reference.
5. **F6-IC**: a target-independent sign-folded pair index built once per
   base (then lift the 256-point cap, whose scan is Θ(|F|) a call), the
   `u128` gate so m ≥ 4 can use it, a geometry-aware split rule.

Algorithmic (changes the question; AGENTS.md §3):

- The measured obstacle is not the kernel. §7 of the plan and the n = 53
  table show every algebraic formulation enumerating about |F|^(m−2)
  partial decompositions before algebra finishes them, and prolonging to
  a higher Macaulay degree barely prunes above the leaves (1.2–1.7× fewer
  nodes for 20–100× more work). "Solving at a high degree of regularity"
  therefore does not buy the crossover; what would is a formulation whose
  refutation cost per random target at a balanced cell is below the
  enumeration reference on the same host, which none of the engines here
  is near. The F5 criterion as implemented prunes 4–8 % of rows at the
  production shapes (the trivial syzygies, exactly), so F5 is not the
  lever either.

## 6. Reproduce

```bash
cargo test --release --lib -- f4_fp:: koblitz_groebner:: inherited_f4:: f6_ic matrix_f5_f2:: gf2_elim
cargo run --release --example gf2_elim_bench -- --quick --reps 3
cargo run --release --example pdp_bench -- --tier reach --engines groebner --budget 60
KIC_F4_KERNEL=legacy cargo run --release --example pdp_bench -- --tier reach --engines groebner --budget 60
cargo run --release --example f6_ic_probe -- --n 17 --a 1 --dim 6 --m 3 --targets 8 --reps 2
KIC_F4_EARLY_ORACLE=0 KIC_F6_PARTIAL_SUPPORT=0 cargo run --release --example f6_ic_probe -- --n 17 --a 1 --dim 6 --m 3 --targets 8 --reps 2 --engines f6
```

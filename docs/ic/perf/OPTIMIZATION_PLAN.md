# Index-calculus performance plan

A standing plan for making the whole index-calculus attack faster, run as
a loop: **measure → change one thing → re-measure against the frozen
baseline → keep or revert → record**. Every step lands with its before
and after numbers, and every change that moves a pinned counter freezes a
new reference beside the old one rather than overwriting it
(`docs/ic/ci/README.md`, `AGENTS.md` §3 and §7).

Started 2026-09-24 on `main` at `9250d0cc` (after #694).

## 1. Where the time goes today

| stage | what exists | measured cost | measured in CI? |
|:--|:--|:--|:--|
| Whole pipeline, m = 3, pair table | `ic workflow`, three frozen rungs | k0n53: 1.65 s end to end (collection 0.95 s, descent 0.70 s) | **yes**, `ic-e2e-benchmark.yml`, reference v5 |
| Decomposition by summation polynomial + Gröbner, m = 2 | `koblitz_groebner` matrix-F4 + splitting | ≤ 0.2 s per stage ladder rung; Macaulay **build** ≈ 60× the **reduce** | no (`groebner_stage_bench`, compared by hand) |
| Same, m = 3 (chained S₃) | same | solving degree 6; 20,240 × 20,686 matrix, 27.4 s dense / 2.9 s sparse; 17-variable cell 82–494 s end to end; 27-variable cells estimated at days | no |
| m = 4 | enumeration arm only; Gröbner builder capped at 64 variables | no recorded Gröbner decomposition | no |
| m = 5, 6 | S₅, S₆ exist symbolically (`semaev_leading_form`, S₆ = 190,252 monomials) and are never Weil-descended | nothing recorded | no |
| GPU | `gpu/semaev`, `gpu/macaulay` | never run on a GPU; no GPU on the dev host | no |

The four blockers for m = 5 and 6, from the code survey:

1. **Monomials are `u64` masks** (`MAX_VARS = 64`) in every F₂ solver.
   Chained m = 4 already needs 66 variables at n = 21; m = 5, 6 cannot
   be expressed.
2. **Elimination** is M4RI with block width 4, scalar XOR loops, a
   serial pivot search, and no sparse–dense split in the solving path —
   although the diagnostic sparse path was already measured 9.6× faster
   on the m = 3 degree-6 matrix.
3. **Macaulay construction** (HashMap column index, full products
   `t·fᵢ` with no symbolic preprocessing) costs more than the
   reduction on every stage rung.
4. **No symmetrised binary S₅/S₆ system and no fast exhaustive or
   hybrid solver above degree 2**. FES/Monica is quadratic-only, and
   chaining pays `(m−2)·n` extra full-field unknowns to stay at degree 3.

## 2. The measurement suite (built first)

A change is kept only if it moves one of these, and it has to leave
correctness untouched.

| layer | harness | unit | gate |
|:--|:--|:--|:--|
| **L0 kernel** | `examples/gf2_elim_bench.rs` (new): reduces frozen Macaulay matrices dumped from real decomposition systems (m = 2, 3; degrees 3–6), plus random matrices of the same shape | ns per row-word; rank and a digest of the reduced basis must match | frozen JSON reference, compared by `scripts/pdp_bench.py` |
| **L1 oracle** | `examples/pdp_bench.rs` (new): a PDP ladder over m = 2…6 with **planted** targets (sum of m factor-base points, so a decomposition is known to exist) and random targets, per solver engine | seconds per decided target, success on planted targets, degree reached, largest matrix, word XORs | frozen reference `docs/ic/perf/pdp-reference-v1.json`; a cell that times out is recorded as *unreached*, never as a negative result |
| **L2 stage** | `groebner_stage_bench` (exists) | word XORs, build/reduce/read split | by hand against `research/groebner_stage_20260915` |
| **L3 end to end** | `ic workflow` + ρ, `ic-e2e-benchmark.yml` (exists) | pinned counters + same-host ρ/IC ratio | CI, reference v5 |

The PDP ladder has three tiers so the frontier stays visible:
**reach** (cells solved today, which pin regressions), **frontier**
(cells that are slow today: chained m = 3 at ℓ ≥ 5, m = 4 at small n)
and **aspiration** (m = 5, m = 6 at the smallest sizes that make sense).
A new reference is frozen whenever a cell moves tier.

## 3. Targets, in order

Each item says what it changes, which layer should move, and how the
change could turn out not to help.

**A. Elimination kernel (L0, then L1, L2).**
1. M4RI with block width 8 and multiple Gray-code tables; AVX-512 row
   XOR (512-bit lanes, `vpternlogq` for three-way XOR) with a portable
   fallback selected at runtime, as the scan kernel does.
2. A parallel pivot search: rows split across threads, and the table
   application parallelised from a lower threshold.
3. A Faugère–Lachartre sparse/dense split in the **solving** path:
   eliminate the sparse pivot block with sparse rows first, then run the
   dense kernel on the remaining Schur complement only. This is the
   measured 9.6× from the diagnostic path, moved into production.
   *Falsifier:* at small sizes the build dominates, so L1 may not move
   until B lands.

**B. Macaulay construction (L2, L1).**
1. Rank monomials with the combinatorial number system (graded
   reverse-lex rank is closed form) instead of hashing them.
2. Symbolic preprocessing: only the products `t·fᵢ` whose leading
   monomials are needed (true F4 selection; `pq_f4_f2` already has the
   skeleton), with the F5 criterion on by default where it proves
   useless rows.

**C. Wide monomials (unblocks m ≥ 4).** A width-generic monomial type
(`[u64; W]`, W = 1, 2, 4) for the F₂ solvers. W = 1 must stay exactly
as fast as today (checked at L0/L1), and W = 2 makes chained m = 4 at
n ≥ 21 and m = 5 representable.

**D. Algorithms for m = 5, 6 (L1: the aspiration tier).**
1. **Symmetrised binary S₅ and S₆.** Express them in e₁…e_m (FGHR
   symmetrisation) from the cached symbolic S₅/S₆ and Weil-descend them
   in the symmetric frame. This divides the solution count by m! and
   lowers the degree the solver has to reach.
2. **Hybrid guess-and-solve** (Bettale–Faugère–Perret): fix k
   variables, walking the 2^k guesses in Gray-code order so consecutive
   systems differ by one substitution (the `inherited_f4` idea, made
   systematic), and choose k by measured cost rather than by formula.
3. **Bitsliced exhaustive search for degree ≤ 4** (FES generalised
   beyond quadratic: Gray-code enumeration with k-th derivatives, 512
   candidates per AVX-512 register), as the leaf solver under the hybrid
   and as the baseline every Gröbner path must beat on the same cell.
4. Crossbred (exists, single-threaded) parallelised and put on the
   ladder as an engine.
   *Falsifier:* if the degree of regularity at m = 5 grows as the notes
   predict (Boolean degree m(m−1)), Gröbner cells may stay unreached,
   and the exhaustive and hybrid solvers set the frontier. That is a
   result, and it is recorded as one.

**E. End-to-end pipeline (L3).**
1. Descent is now 43 % of the k0n53 rung. Candidates: an AVX-512
   batched two-summand walk key, and a pair-table width priced for the
   descent as well as for collection.
2. Selection and table build at larger bases (`what_is_left` in
   `koblitz-collection-aim-20260922.json`: 23 % and 29 %).
3. When D reaches it, an **m = 4 rung** on the e2e gate, so the
   decomposition oracle's progress shows up end to end.

**F. GPU.** Port the L0 kernel and the D.3 exhaustive leaf to CUDA
(`gpu/macaulay` is the starting point). Everything is written to
compile in the existing `nvcc-compile` CI job. It **cannot be measured
here**, because the dev host has no GPU. It is benchmarked only once a
GPU host is available, and is never claimed from CPU numbers.

## 4. How each step is reported

Each step is a commit (or small PR) carrying a before/after table from
the layers it claims to move and a statement of the layers it did not
move. The step is classified by `AGENTS.md` §3: a faster kernel is
**engineering**; a lower degree of regularity or fewer unknowns for the
same decomposition question is an **algorithmic** change and says so. A
regression on any gated layer is reverted, not explained away.

## 5. Log

| date | step | layer moved | before → after | reference |
|:--|:--|:--|:--|:--|
| 2026-09-24 | #694 AVX-512 scan kernel | L3 | k0n53 collection 4.45 s → 3.55 s | v4 (counters identical) |
| 2026-09-24 | #694 windowed, aimed collection rungs | L3 | k0n53 IC 4.41 s → 1.65 s | v5 |
| 2026-09-24 | A.1 + A.2: `gf2_elim` Four Russians kernel (up to four adaptive-width Gray-code tables per pass, word-strip pivot search, BMI2 `pext` pattern gather, AVX-512 row update, rayon) replaces block-width-4 M4RI for the oracle's reduced row echelon form | L0; L1 dense solving | m = 3, degree 6 (20,240 × 20,686): elimination 4.12 s → 0.55 s (7.5×); dense `solving_profile` 3.84 s → 0.86 s, now ahead of the sparse path (1.98 s). Stage ladder unmoved: its matrices sit below the kernel's size gate, and the build dominates it (B) | `gf2-elim-reference-v1.json`; reduced form bit-identical on every cell |
| 2026-09-24 | B (part): decomposition systems instantiated from a per-thread `DecompositionTemplate` memo; Macaulay row cap read once per matrix | L2 | system build was 29 % of an m = 2 decomposition at n = 23 (callgrind); stage-ladder wall within its noise (tens of ms) | equations identical (new test) |
| 2026-09-24 | L1 frozen: `examples/pdp_bench.rs`, balanced ladder m = 2…6 | — | see §6 | `pdp-reference-v1.json` (a few frontier timings overlapped another job on the host; hits and statuses are unaffected) |
| 2026-09-24 | E.1: descent walk steps all 64 walks with `add_many_lazy(G, walks)` (the AVX-512 kernel) instead of `add_pairwise` | L3 | k0n53 descent 0.729 s → 0.622 s (−15 %, paired, 3 runs each); k0n41 unchanged within noise | v5, every counter identical |
| 2026-09-24 | E.1: the workflow's solve stage decides its targets in parallel (rayon, state written per target under a lock) | L3 wall only | solve-stage wall: k0n53 0.756 s → 0.183 s (4.1×), k0n41 0.212 s → 0.068 s (3.1×). The gated `descent_seconds_total` sums per-target walls, so the gate prices the same work as before — by design, since ρ runs on one core | v5, every counter identical |
| 2026-09-24 | D.2 probe: summand-first splitting (`SOLVER_SPLIT_RULE=lowest`) at n = 31, m = 3 | L1 | planted 0/2 → 1/2, at 63 s a target against meet in the middle's 0.1 ms | exploratory; default unchanged |
| 2026-09-24 | E.2: the windowed and aimed collection scan (`witnesses_fast_scan`) forms its rests with `add_many_lazy` — the AVX-512 kernel the full scan already used — completing ordinates only for keys the filter admits | L3 | k0n53 unit of 2048 probes on one thread 90 ms → 60 ms; precompute 0.95 s → 0.80 s, IC whole process ≈ 1.39 s | v5, every counter identical |
| 2026-09-24 | E.3: the folded table's orbit keys (`keys_of`, used by the full scan, the windowed/aimed scan and the descent) are computed in bulk by `FrobeniusCanon::canon_in_place`, a branch-free all-rotations minimum on AVX-512, sixteen keys per call | L3 | paired, 3 runs each: k0n53 IC whole process 1.43 s → 1.05 s (−26 %; descent 0.61 → 0.43 s, precompute 0.82 → 0.62 s); k0n41 0.30 s → 0.21 s (−29 %). In isolation the bulk key is slower than the scalar one (16.6 vs 13.4 ns); in the scan it is faster, for a reason not established | v5, every counter identical |

## 6. What the L1 ladder says (pdp-reference-v1)

Balanced cells (`m·ℓ ≈ n`), seconds per decided target, planted hits:

| cell | vars | Gröbner | SAT | meet in the middle | enumerate |
|:--|--:|:--|:--|:--|:--|
| n = 23, m = 2, ℓ = 11 | 22 | 10 ms, 8/8 | 1.5 s, 8/8 | 0.1 ms, 8/8 | 1.3 ms |
| n = 15, m = 3, ℓ = 5 | 30 | 35 ms, 8/8 | budget | < 0.1 ms | 0.9 ms |
| **n = 31, m = 3, ℓ = 10** | 61 | **3.7 s, 0/2 planted** (node budget) | budget | 0.1 ms, 2/2 | 0.9 s |
| n = 15, m = 4, ℓ = 4 | 46 | 69 ms, 4/4 | budget | < 0.1 ms | 2 ms |
| n = 15, m = 5, ℓ = 3 | 60 | 1.8 s, 1/1 | budget | 0.3 ms | 1.6 ms |
| n = 7, m = 6, ℓ = 1 | 34 | 14 ms, 4/4 | 2.3 ms | < 0.1 ms | < 0.1 ms |
| n = 31, m = 4 / 5 / 6 | 86 / 123 / 154 | not buildable (64-bit monomials) | — | 0.2 ms / 7 ms / 0.21 s | 55 ms / 0.75 s / — |
| n = 39, m = 3, ℓ = 13 | 78 | not buildable | — | 0.6 ms | — |

Reading it honestly:

1. On every cell the summation-polynomial + Gröbner oracle is **two to four
   orders of magnitude slower** than the exact meet-in-the-middle oracle
   the end-to-end pipeline already uses. The elimination kernel (A) moved
   a constant; it did not touch that gap.
2. At n = 31, m = 3 it is also **incomplete**: its splitting search
   exhausts 20,000 nodes (and 200,000, at 38 s a target) without finding
   decompositions that exist. The default split rule branches on the
   31-bit intermediate abscissa first. Guessing a summand instead and
   solving the m = 2 remainder is the obvious fix, and it would still cost
   about |F| × 5 ms ≈ 7 s a target, against 0.1 ms.
3. m ≥ 4 at n ≥ 31 cannot even be written down until the monomials are
   widened (C). Doing that measures the frontier; the numbers above give
   no reason to expect it to cross meet in the middle.

So the plan's order stands for the Gröbner track (C, then D.2 hybrid
guessing with summand-first splitting), with its success criterion stated
in advance: **a Gröbner cell counts as progress when it decides planted
targets that it previously missed, or when its seconds per target fall on
an unchanged cell; it counts as a crossover only when it beats the
meet-in-the-middle column on the same cell.** End-to-end work (E) proceeds
alongside it.

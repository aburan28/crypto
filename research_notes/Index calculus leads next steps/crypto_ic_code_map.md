# Index calculus code map for aburan28/crypto: what the code implements and where it falls short of the claims

Scope: read-only survey of `/home/user/crypto` on 2026-10-14 (src/, gpu/, sage/, scripts/, examples/, tests/, benches/, groebner_compare/; hdl/, isogeny_algos/ and quantum_circuit_secp256k1/ had nothing IC-specific). The sources are local files, cited as `path:line`. Sizes are `wc -l`. Test counts are the number of `#[test]` attributes in the file. Inline tests are often `#[ignore]`d or toy-scale, so a test count shows that tests exist, not that the code is validated at scale.

Scale: 231 files in src/, 296 in examples/, 225 in scripts/, 31 in gpu/, 12 in tests/ and 4 in sage/ match "index calculus | summation poly | semaev | factor base". `src/cryptanalysis/` holds 239 modules.

## Q1. Which modules implement each IC stage? (path, size, language, maturity)

### Takeaway
Every IC stage has real Rust code with unit tests, so very little exists only in documents. Most of it is toy-scale or binary/Koblitz-specific. Two pipelines run end to end: the binary Koblitz pipeline (`koblitz_index_calculus.rs`, 21k lines) and the pluggable `ic_framework` that ecbench calls. The strongest sparse linear algebra (block Wiedemann, block Lanczos, filtering) is wired only into the Koblitz pipeline. The framework and ecbench offer only `incremental-gauss*` and `structured-gauss`.

### Cited Findings

**Orchestration / pipelines**
- `src/cryptanalysis/koblitz_index_calculus.rs`: 21,321 lines, Rust. The end-to-end binary Koblitz IC. It imports `koblitz_sparse_la` for the log solve (`:164`, `:16785` sets `BlockWiedemannOptions`). Status: tested and benchmarked.
- `src/cryptanalysis/ic_framework/`: mod.rs 1,505 lines, plugins.rs 2,464, solvers.rs 1,431, linalg.rs, large_prime.rs, rootfind.rs, fiber.rs 1,277, shared_rank.rs, residual.rs, prime_tower.rs and stages.rs. This is the pluggable pipeline "factor base → targets → decomposition → relations → logarithm" with a `SystemSolver` plug for F4/F5/XL/SAT. Its doc says every phase counts toward the total `S` (`ic_framework/mod.rs:1-60`). Status: tested (`tests/ic_framework.rs`).
- `src/cryptanalysis/ic_boundary.rs`: 10,269 lines, 44 tests. The boundary ledger, which measures fixed variants end to end.
- CLIs: `src/bin/ic/*.rs` (experiment.rs 1,859, workflow.rs 2,480, price.rs 1,554, rho.rs 1,234, descent.rs, large_prime.rs), `src/bin/icprog*` (rule.rs 2,354, b5a.rs 2,382, rounds.rs 1,633), and `src/bin/icx.rs`. Integration tests: `tests/ic_fixed.rs`, `tests/ic_framework.rs`, `tests/ic_progress.rs`, `tests/icprog.rs`, `tests/icx.rs` and `tests/ecbench_large_prime.rs`.

**ecbench IC arms**
- Method IDs: `ic.pipeline`, `ic.pipeline_counted`, `ic.large_prime` and `ic.shared_rank`. `ic.pipeline` takes `factor_base`, `oracle`, `solver` and `linalg` (default `incremental-gauss`), plus `targets`, `max_trials` and `solver_budget_seconds` (`src/cryptanalysis/ecbench/methods.rs:297-410`).
- Prime factor bases: `prime-abscissa`, `glv-orbit` and `prime-kummer`. Prime oracles: `subtract`, `mitm`, `mitm-prime-batch`, `mitm-fiber`, `rootfind-mitm`, `rootfind-pairs` and `rootfind-subtract` (`methods.rs:1803-1827`).
- Binary factor bases: `binary-subspace`, `koblitz-orbit`, `compact-orbit-scan` and `cyclotomic`. Binary oracles: `subtract`, `mitm`, `mitm-frobenius(-counted|-fiber)`, `rootfind-*`, `descent-algebraic` and `descent-algebraic-fiber`, where the last two use a named system solver (`methods.rs:1915-2005`).
- Extra ecbench IC binaries: `src/bin/ecbench/pdp3_koblitz.rs` (1,085 lines) and `pdp_split_oracle.rs`. Tests: `tests/ecbench.rs` (2,080 lines).

**Summation-polynomial construction**
- Prime-field `S_3` and a univariate-in-`X_4` extraction of `S_4`: `src/cryptanalysis/ec_index_calculus.rs` (2,404 lines, 20 tests). Its own doc calls it "toy-scale … sub-40-bit" (`:1-15`).
- `S_3` through `S_6` over prime fields by Semaev's resultant recurrence: `src/cryptanalysis/semaev_higher.rs` (889 lines, 6 tests).
- Binary `S_4`: `binary_semaev.rs` (577 lines), `binary_semaev_s4.rs` (1,933 lines, 23 tests) and `semaev_leading_form.rs` (1,371 lines, 16 tests).
- Symmetrised (`S_m` / elementary symmetric) form: `src/cryptanalysis/symmetrized_semaev.rs` (563 lines, 6 tests; doc at `:1-20`) and `koblitz_sym_s4.rs` (832 lines; binary `S_4` in `e_1, e_2, e_3`).
- Koblitz 2-torsion invariant frame (`u = 1/(x+1)`, `w = u²+u`): `koblitz_symmetrised.rs` (3,468 lines, 17 tests).
- FGHR `D_m = (Z/2)^{m−1} ⋊ S_m`: `src/cryptanalysis/fghr_line.rs` (1,549 lines; doc `:1-8`), `fghr_full.rs` and `gaudry_quartic_d4.rs` (541 lines). These were benchmarked on 2026-10-13: the `C_4` cost fell 47.5–48.4× (`1.21·10¹²` → `2.54·10¹⁰` `F_p` multiplications) at `p = 269` (benchmarks/README.md row `20261013-gaudry-quartic-d4`).

**Decomposition solvers**
- Binary system solvers in the framework (`ic_framework/solvers.rs:15-25`): `f4-f2`, `buchberger-f2`, `matrix-f4`, `matrix-f5`, `inherited-f4`, `crossbred-f2` (Joux–Vitse crossbred), `xl-f2`, `sat-cdcl`, `fes-f2`, `fes-f2-wide` (AVX2/AVX-512) and `exhaustive`. The backing modules:
  - F4 and Gröbner: `koblitz_groebner.rs` (9,432 lines, 56 tests), `pq_f4_f2.rs` (2,610), `pq_groebner_f2.rs` (1,777), `inherited_f4.rs` (3,567), `matrix_f5_f2.rs` (1,696), `wide_groebner.rs` (2,124) and `sparse_macaulay.rs` (413 lines, 0 tests).
  - Other system solvers: `crossbred.rs` (1,217), `pq_xl.rs` (419), `sat.rs` (3,250), `semaev_sat.rs` (1,792), `binary_semaev_chain_sat.rs` (1,187), `wdsat_oracle.rs`, `mq_fes.rs` (1,228) and `mq_monica.rs`.
- Prime-field F4: `f4_fp.rs` (3,513 lines), `f4_fp_tower.rs` (3,227) and `groebner_f4.rs` (482).
  - `groebner_f4.rs` doc: "Out of scope. Real F4 (signature-based reduction, … resolution — all left as future work" (`src/cryptanalysis/groebner_f4.rs:36-39`).
  - External comparison: msolve on 1 thread is 12–24× faster than `f4_fp` on square generic systems (quad_n9: 246 ms against 5,693 ms). On Boolean `n = 17` the in-house decider takes 7–24 ms per target against msolve's 42.7–46.5 s for the full basis (benchmarks/README.md row `20261011-groebner-vs-msolve-singular-polybori`). That Python/Julia/C++ harness is in `groebner_compare/` (runner.py, m5gb_bridge.cpp, worker.jl).
- Root finding (linearized subspace polynomial `gcd(q, L_V mod q)`):
  - Binary `S_4`: `semaev_decomp.rs` (4,041 lines, 30 tests; doc says almost all IC time goes here, `:1-12`).
  - Prime 2-decomposition: `prime_fast.rs` (1,797 lines) and `ic_framework/rootfind.rs`.
  - Benchmarked: a whole IC solve with the root finder is 350× faster in wall time than with F4 at `n = 13`, but IC is still 2.9–25,116× strong rho in `S` (binary) and 29–141× rho (prime), with the gap widening (row `20261013-ic-rootfind-ecbench`).
- Resultants: implemented inside `semaev_higher.rs` (`S_4 = Res_X(S_3, S_3)` comment, `:300-316`), `koblitz_sym_s4.rs`, `gaudry_cubic.rs`, `fghr_line.rs`, `semaev_leading_form.rs` and `hyperelliptic_index_calculus.rs`. There is no standalone resultant-based decomposition solver plug-in.

**Factor-base constructions**
- Framework bases (`ic_framework/plugins.rs:12-18`): `prime-abscissa`, `binary-subspace`, `koblitz-orbit` (signed Frobenius orbit folding), `compact-orbit-scan`, `cyclotomic` (γ-weight abscissae), `glv-orbit` (automorphism orbits: j = 0 folds 6 points per column, j = 1728 folds 4) and `gls-line` (`⟨−1, ψ⟩`-orbits over `F_{p²}`).
- General endomorphism-invariant base: `glv_invariant_base.rs` (1,987 lines, 15 tests).
- Other base modules: `koblitz_factor_base_search.rs` (1,860), `cyclotomic_factor_base.rs` (1,135), `p256_dickson_factor_base.rs`, `p256_bitbox_factor_base.rs`, and the binaries `src/bin/p256_factor_base_symmetry.rs`, `p256_root_orbit_factor_base.rs` and `p256_fb_log_transport.rs`.
- Sage scripts: `sage/factor_base_dor.sage` (1,265 lines) and `sage/low_hw_factor_base.sage`.

**Relation collection and large primes**
- Framework collector `ic_framework/large_prime.rs` supports 0, 1 or 2 large primes (`:17`, `:82-83`: "large_primes must be 0, 1 or 2"). It is exposed as ecbench `ic.large_prime` (`methods.rs:350-390`).
- `koblitz_large_prime.rs` (1,397 lines, 4 tests) does exact elimination of 1LP/2LP relations over signed Frobenius orbits (`:1-8`).
- `ecbench_large_prime.rs` (966 lines) and `gaudry_cubic.rs` also contain double-large-prime variants.
- `large_prime_filter.rs` (489 lines, 5 tests) is a general 1/2/k-large-prime hypergraph peel plus Markowitz elimination with witness replay (`:1-14`).

**Sparse linear algebra**
- `koblitz_sparse_la.rs` (4,944 lines, 33 tests). Full NFS-style filtering: duplicates, singletons, the clique rule, structured Gaussian elimination under a fill bound (`filter_relations`, `:1417`). Block Wiedemann with matrix Berlekamp–Massey (`block_wiedemann_kernel`, `:2377`) and block Lanczos (`block_lanczos_solve`, `:2730`), generic over word widths, with a GPU/worker SpMV backend (`SpmvBackend::Worker`).
- `pq_wiedemann.rs` (1,053 lines, 12 tests) is scalar Wiedemann plus Berlekamp–Massey over `Z/N`. `gaudry_cubic.rs:3319,3439` has its own `wiedemann_u64` and `wiedemann_big`.
- `pq_sparse_la.rs` (831 lines) is Markowitz structured Gaussian elimination. `ic_framework/linalg.rs` provides `incremental-gauss`, `-full-rank`, `-full-rank-checked` and `structured-gauss` (`:9-15`, `:250-275`). `gf2_elim.rs` (1,428 lines) does `F_2` elimination.

**Descent and lifting**
- Weil descent:
  - `pq_descent.rs` (788 lines) and `pq_descent_symbolic.rs`: binary Semaev `S_3` to Boolean systems.
  - `ghs_descent.rs` (747 lines, 4 tests): GHS `E(F_{2^N}) → Jac(C)(F_{2^l})`.
  - `diem_descent.rs` (513 lines, 5 tests): "Diem-style … toy implementation" (`:1`).
  - `weil_charts.rs`.
- Cover and decomposition: `jv_cover.rs` (3,498 lines; JV12 genus-3 GHS cover with Nagao decompositions, `:1-14`), `jv_sieve.rs` (2,399 lines; the JV12 §3.2 sieving variant), `jv_quintic.rs`, `jv_quintic_edwards.rs`, `jv_quartic.rs` and `jv_isogeny_walk.rs` (2,734 lines).
- Gaudry `E(F_{p^k})`: `gaudry_cubic.rs` (4,679 lines), `gaudry_quartic.rs` (2,379) and `glv_gaudry.rs` (2,803; includes a "function-first (Nagao L(4O))" system, `:29`, `:1643`).
- Individual-log descent:
  - `recursive_descent.rs` (464 lines, 2 tests): a budgeted scheduler with Frobenius-orbit canonicalisation.
  - `ic_descent_degrees.rs` (1,068 lines).
  - `descent_expansion.rs`, `descent_lowgamma.rs`, `descent_treewidth.rs` and `descent_algebraic.rs`: rank-profile diagnostics.
- Lifting: `canonical_lift.rs` and `cm_canonical_lift.rs`. xedni appears only as a one-point "xedni-style" height computation (`examples/nist_cyclotomic_point_heights.rs:24,1164-1176`).

**GPU kernels** (CUDA C++, plus Metal for FES). Every README says no kernel has been run on a GPU. All are validated by host emulation against Python or CPU oracles.
- `gpu/semaev`: binary `S_4` pairs-and-solve sweep ("No kernel here has been run on a GPU", README). No Rust caller found.
- `gpu/rootfind`: prime 2-decomposition, with a separate cargo crate `rs/`. "the root crate has no feature/dep/build script for it" (README). `src/cryptanalysis/gpu_rootfind.rs` (258 lines) is only the CPU reference.
- `gpu/macaulay`: batched `F_p` row reduction for `gaudry_cubic`'s residual solves, which are 57% of the solve cost (README, citing RESIDUAL_WALKS §11.6). No Rust caller found.
- `gpu/fes`: FES Gray-code search, wired through `src/cryptanalysis/fes_gpu.rs` (462 lines) into `icx` as a subprocess. Its outputs are re-verified on the CPU.
- `gpu/spmv`: SpMV worker implementing the SPMV1 contract for `koblitz_sparse_la`. Wired; any mismatch falls back to the CPU (README).
- `gpu/ecc2k`: Koblitz rho and a pair table. This is rho, not IC.

### Inferences
- The repository's main gap is integration, not missing algorithms. Most stages have more than one implementation (three Wiedemann variants, four structured-Gauss variants, three large-prime collectors), but ecbench `ic.pipeline` can reach only the weakest linear algebra.
- Every throughput figure attached to GPU code is an estimate. No GPU measurement exists anywhere.

### Gaps
- I did not run any tests, so the "tested" label means only that tests exist. I did not check which tests are `#[ignore]`.
- I read only file headers for `ic_boundary.rs` and `koblitz_index_calculus.rs`, so their internal sub-stages are not mapped.

## Q2. Which algorithms from the literature are missing?

### Takeaway
Most of the named techniques are present in some form: block Wiedemann and block Lanczos, FGHR `(Z/2)^{m−1} ⋊ S_m` symmetrisation, Frobenius and GLV orbit factor-base reduction, the Joux–Vitse crossbred solver and cover-and-decomposition, double large primes, and Nagao and Diem decompositions. The real absences are these:
- xedni, beyond a single height computation.
- Joux–Vitse *variable splitting* / the `n−1` / "last fall degree" style splitting for the large-`m` regime. I found only the crossbred solver and covers.
- A signature-based F5 / msolve-class F4 for `F_p`.
- Diem decomposition at non-toy scale.
- A resultant or triangular-decomposition solver as a framework plug-in.
- A Petit-style smooth-model GHS descent (marked TODO).

### Cited Findings
- Block Wiedemann and block Lanczos are present: `src/cryptanalysis/koblitz_sparse_la.rs:2377` and `:2730`. A comment at `koblitz_sparse_la.rs:301` says modular division was "most of the block Wiedemann time" before a fix.
- Symmetrisation:
  - Plain `S_m` symmetrisation: `symmetrized_semaev.rs:1-20`.
  - FGHR `D_m = (Z/2)^{m−1} ⋊ S_m`: `fghr_line.rs:1-8`.
  - The binary 2-torsion frame: `koblitz_symmetrised.rs:1-12`.
- Orbit factor-base reduction:
  - Frobenius orbits: `koblitz-orbit`, `compact-orbit-scan` (`ic_framework/plugins.rs:14-15`) and `koblitz_large_prime.rs`.
  - GLV/automorphism orbits: `glv-orbit` and `gls-line` (`plugins.rs:17-18`), `glv_invariant_base.rs` and `ec_index_calculus_j0.rs` (S₃ / order-6 relation amplification, 1,152 lines). `ec_index_calculus_j0.rs:77` says it "exposes the S₃ equivariance so that future work could …", so the equivariant Gröbner side is not implemented.
- Joux–Vitse:
  - Crossbred solver: `crossbred-f2` in `ic_framework/solvers.rs:20`; `crossbred.rs` is 1,217 lines.
  - Cover and decomposition plus sieving: `jv_cover.rs` and `jv_sieve.rs`.
  - No module, solver name or comment matched "variable splitting".
- xedni: the only occurrences are `examples/nist_cyclotomic_point_heights.rs:24,1164,1176`, the "xedni-style lift of G" used for canonical-height statistics. There is no xedni lifting attack.
- Gaudry's double large prime: implemented as `large_primes = 2` (`ic_framework/large_prime.rs:17,82`), in `koblitz_large_prime.rs:1` ("one- and two-large-prime relations") and in the general k-large-prime filter `large_prime_filter.rs:1-14`.
- Nagao decomposition: `jv_cover.rs:12-14` (Nagao decompositions on a genus-3 Jacobian), `glv_gaudry.rs:1643` (Nagao `L(4O)` function-first system), and mentions in `jv_sieve.rs`.
- Diem decomposition: `diem_descent.rs:1` ("toy implementation", 513 lines). Its only caller is `examples/diem_descent_demo.rs`.
- F4/F5 quality: `groebner_f4.rs:36-39` leaves signature-based reduction "as future work". The Boolean side has `matrix-f5` (the F5 criterion rows-skip, `solvers.rs:18`), not a full F5.
- GHS smooth model: `ghs_descent.rs:63` and `:310-317` say that producing the smooth model "requires Newton-polygon / Puiseux desingularisation — left as a documented TODO".
- Hyperelliptic IC exists (`hyperelliptic_index_calculus.rs`, 2,801 lines, 22 tests; `hyperelliptic_ic_bench.rs`), as do the `bielliptic_quartic` and `jacobian_transfer` modules.

### Inferences
- Of the techniques the brief lists, the only ones absent are xedni (which is not useful for ECDLP anyway, so this is low priority) and Joux–Vitse splitting. Diem decomposition and GHS descent exist but only at toy or incomplete scale.
- The highest-leverage literature gap is solver quality on `F_p`: there is no msolve-class F4, and the benchmark shows a 12–24× deficit on square systems. The workaround already exists: `groebner_compare/` can call msolve as an external engine.

### Gaps
- I did not establish whether `fghr_full.rs` extends the `D_m` symmetry to `m ≥ 5` or to prime fields generally; I read only its existence and its caller.
- I did not check whether the "Summation polynomial with Frobenius action" constructions in the Koblitz path use the full `⟨π⟩ × D_m` group.

## Q3. TODO/FIXME/unimplemented markers near IC code

### Takeaway
The IC code carries almost no explicit TODO markers. There are no `todo!()` or `unimplemented!()` calls and essentially no FIXMEs. Gaps are recorded in prose ("future work", "not implemented", "placeholder") and through `#[allow(dead_code)]`.

### Cited Findings
- `todo!(` / `unimplemented!(`: zero hits in `src/cryptanalysis`, `src/bin/ic*`, `src/bin/ecbench*`, `gpu/`, `sage/`, `groebner_compare/`, `tests/` and `benches/`.
- `NotImplementedError("even m")`: `gpu/ecc2k/ecref2k.py:159`, in the Koblitz reference, which serves rho rather than IC.
- TODO: `src/cryptanalysis/ghs_descent.rs:63` and `:316` (smooth model of the GHS descent target).
- "future work": `groebner_f4.rs:36-39` (signature F4/F5) and `ec_index_calculus_j0.rs:77` (S₃ equivariance).
- "not implemented":
  - `src/bin/ic/params.rs:432` ("general binary-curve point counting is not implemented").
  - `:472` ("normal/abstract basis is described but not implemented by this inspector").
  - `:570` ("general prime-curve point counting is not implemented").
  - `descent_treewidth.rs:43`.
- Placeholder code: `semaev_higher.rs:313` (`let _ = (a_term, …); // placeholder`). The resultant-in-Y derivation was abandoned in favour of `s4_in_last_proper`.
- "unwired" or dead code:
  - `groebner_cache.rs:8-10` says "Fail-open cache infra, currently unwired: kept for upcoming use".
  - `#[allow(dead_code)]` appears in `koblitz_index_calculus.rs:695,783,5705,8677,9657`, `inherited_f4.rs:513` ("rewrite-plan cache reserved for a later path"), `:551,1021,1811`, `gaudry_cubic.rs:651,1110,2357`, `semaev_decomp.rs:1205`, `sat.rs:591` and `src/bin/ecbench/pdp3_koblitz.rs:538` ("future adapters").

### Inferences
- Because markers are rare, a TODO grep understates the unfinished work. The "No kernel has been run on a GPU" disclaimers and the orphan modules listed in Q5 are the more informative signal.

### Gaps
- I did not scan `examples/` (296 IC files) or `scripts/` (225) exhaustively for TODOs. Those are mostly research drivers.

## Q4. Performance bottlenecks according to code comments and benchmarks

### Takeaway
The code and benchmarks agree on three bottlenecks:
1. The decomposition oracle, especially on targets that do not decompose.
2. Macaulay elimination inside Gröbner solves.
3. Below that, the pair-table memory and probe cost.

Linear algebra is small at the measured sizes. IC remains 3× to more than 10⁴× slower than rho in `S` on every measured curve.

### Cited Findings
- `semaev_decomp.rs:1-12`: "Almost all of an attack's time goes here, and almost all of *that* on the ~5/6 of targets that do not decompose."
- `gpu/macaulay/README.md`: the `gaudry_cubic` residual solve is `0.88·10⁶` `F_p` multiplications, 57% of it in the Macaulay step (17% forward elimination, 40% normal forms). The matrix shape is 226×286 at rank 222 for `p = 2083`.
- `sparse_macaulay.rs:1-8`: the dense `F_2` reduction at `n = 5, m = 3` is 8,340×21,778 and takes about 8 minutes, and degree 6 at 27 unknowns needs 397,594 rows. That is the reason for structured elimination.
- `koblitz_sparse_la.rs:301`: modular division dominated block Wiedemann time before the fix.
- `ic_framework/linalg.rs:15-25`: Markowitz pivoting saves only a few percent of `row_ops` on IC relation matrices (538 against 582 on 63 rows × 121 columns). The real cost of a dense matrix is the column scan and the `rank×|F|` storage.
- `ic_framework/solvers.rs:837`: XL's work is dominated by the linearised elimination.
- Benchmarks (`benchmarks/README.md`):
  - `20261013-ic-rootfind-ecbench`: IC is 2.9–25,116× strong rho in `S` (binary) and 29–141× (prime), with the gap widening.
  - `20261013-ecbench-frontier`: the cold frontier is still signed-Frobenius rho. The compact orbit base is 9.985× matched rho at n17, and the retained-setup online frontier at n61 is 13.9× [9.8, 19.1] (rho/IC). 48 inherited-F4 solver-counter mismatches turned up in replay.
  - `20261013-prime-rootfinding-arity`: `m = 3` costs 1.43× `m = 2` in counted GAE.
  - `20261012-pdp-root-finder`: `rf-mitm` uses `5|F| − 3` multiplications per target, with "every slope at its generic floor".
  - `20261013-gaudry-quartic-solve-speed`: the `S_5` solve takes 164–203 s.

### Inferences
- Better linear algebra will not close the rho gap at current sizes. The levers that can change `S` are the decomposition probability and cost per call: symmetry folding (as in the FGHR `D_4` result) and faster Macaulay or F4.
- The 48 inherited-F4 counter mismatches in replay are a correctness or reproducibility issue for the main Gröbner oracle.

### Gaps
- There is no profile of the share of total `S` spent in linear algebra at a scale (`|F| ≥ 10⁵`) where block Wiedemann would matter. No benchmark row exercises `block_wiedemann_kernel` at scale.

## Q5. Unused or half-wired components that could be finished cheaply

### Takeaway
Six components could be wired in cheaply:
1. Expose `koblitz_sparse_la`'s filter plus block Wiedemann/Lanczos as an `ic_framework` linalg plug-in.
2. Wire the orphan modules `large_prime_filter` and `recursive_descent`.
3. Hook the `gpu/rootfind` crate and the `gpu/semaev` and `gpu/macaulay` kernels into Rust behind a feature flag.
4. Wire or delete `groebner_cache`.
5. Add tests to `sparse_macaulay` (it has none).
6. Finish the `semaev_higher` placeholder or remove it.

### Cited Findings
- ecbench and the framework can choose only `incremental-gauss`, `-full-rank`, `-full-rank-checked` and `structured-gauss` (`ic_framework/linalg.rs:250-275`; default in `ecbench/methods.rs:120,320`). `block_wiedemann_kernel` and `block_lanczos_solve` are called only from `koblitz_sparse_la.rs` itself, `koblitz_index_calculus.rs` (`:16785`), the `ic`/`icprog` CLIs and `examples/perfbench/relation.rs`. The framework's `RelationSolver` contract is incremental (`linalg.rs` "The contract both obey"), so a batch solver needs a deferred-solve adapter.
- `large_prime_filter.rs` (489 lines, 5 tests) has no reference anywhere in src/, examples/ or tests/ other than its `pub mod` line (`src/cryptanalysis/mod.rs:231`).
- `recursive_descent.rs` (464 lines, 2 tests) is in the same state (`mod.rs:296`).
- `groebner_cache.rs:8-10` is "currently unwired".
- `gpu/rootfind/README.md`: the root crate has no feature, dependency or build script for the `rs/` crate. `gpu/semaev` and `gpu/macaulay` have no Rust caller (a grep for `gpu/semaev` and `gpu/macaulay` in `.rs` found none), unlike `gpu/fes` (`fes_gpu.rs`) and `gpu/spmv` (`koblitz_sparse_la` `SpmvBackend::Worker`).
- `sparse_macaulay.rs` has 0 tests, though `koblitz_groebner.rs`, `koblitz_bench.rs` and `examples/perfbench/gf2_la.rs` use it.
- `diem_descent.rs` has only `examples/diem_descent_demo.rs` as a caller. It is not a framework oracle.
- `symmetrized_semaev.rs` is used by `groebner_f4.rs`, `coordinate_descent.rs` and the semaev root-vs-Gröbner examples, but it does not appear among the framework's prime oracle names (`methods.rs:1821-1827`).
- `bin/ic/params.rs:432,570`: point counting is not implemented for general binary or prime curves, which limits the inspector to curves with known orders.

### Inferences
- The cheapest high-value step is a `block-wiedemann`/`block-lanczos` linalg plug-in that buffers rows and solves on demand (reusing `filter_relations` plus `block_wiedemann_kernel`), together with an `ic.large_prime` path through `large_prime_filter`. Both are mostly glue around tested code.
- Wiring the GPU kernels is low value until someone has GPU access, since nothing has run on a device.
- Adding a `symmetrized`/`fghr` prime oracle to `ic.pipeline` would let ecbench measure the one symmetry result (FGHR `D_4`, 48×) that has moved counted costs.

### Gaps
- I could not check whether there is a reason, recorded in `docs/ic/FRAMEWORK.md` or PR history, that the block solvers were kept out of the framework. Neither was read.
- The wiring effort estimates here are judgments from file structure, not measured.

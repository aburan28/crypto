# Index calculus in aburan28/cryptanalysis (local clone /home/user/aburan28/cryptanalysis)

Scope note: I read local files only (clone at commit `9b45e8c`, shallow history of 1 commit, about 797 MB). Sources are given as repo-relative paths. `cryptanalysis:` means `/home/user/aburan28/cryptanalysis/` and `crypto:` means `/home/user/crypto/`. I ran no builds or computations. The repo is very large: about 1,046 files under `experiments/` match IC keywords, and `experiments/groebner-perf-20260924/` alone has 90 `round*` directories. This inventory reads each study's headline README/RESULT and conclusion lines. It did not read every round.

## Which IC methods are implemented, and how mature is each?

### Takeaway
IC in this repo is mature as research infrastructure: Rust modules in `suite/` total about 80k lines. It has several interchangeable point-decomposition (PDP) oracles, a fast Boolean F4/F5 kernel with CUDA offload, sparse linear algebra modulo r, and factor-base search. Every one of these stays a toy or stage-level tool for binary and Koblitz ECDLP. The only IC that beats rho is the C library's prime-field multiplicative-group IC, `(Z/pZ)^*`, which crosses rho at about 56 bits. The binary-curve IC pipeline has "no complete released DLP command yet" in the C CLI.

### Cited Findings
**C library (narrow, production-grade)**
- `ca_indexcalc.h` provides IC in `(Z/pZ)^*` only: linear sieve, Pohlig-Hellman for small factors, structured elimination plus Lanczos, Hensel lifting and verified logs. It handles `p < 2^63` at cost `L_p[1/2,1]` and solves a 56-bit p in 1.3 s — [cryptanalysis:README.md lines 44](README.md). The implementation is `src/indexcalc.c`, 835 lines, plus a 102-line header — `cryptanalysis:src/indexcalc.c`, `include/cryptanalysis/ca_indexcalc.h`.
- The README says: "`ic` currently solves DLPs in the multiplicative group of a prime field; the experimental binary-curve IC pipeline is separate and has no complete released DLP command yet" — [cryptanalysis:README.md ~line 146](README.md).
- `cuda/`, `fpga/`, `gpueff/` and `usecases/` contain no IC keyword hits (grep count 0). Their GPU and FPGA work is the ECC2K-130 rho walker. IC GPU code lives in `suite/cuda/` instead (`f4_gf2_device.cuh`, `f4_gf2_echelon.cuh`) — grep over the listed directories; [cryptanalysis:suite/docs/ic/README.md](suite/docs/ic/README.md).

**Rust suite (`suite/src/cryptanalysis/`, line counts from `wc -l`)**
- Koblitz/binary IC core: `koblitz_index_calculus.rs` 14,243 lines; `koblitz_groebner.rs` 7,443; `koblitz_symmetrised.rs` 3,455; `koblitz_sparse_la.rs` 2,601; `koblitz_factor_base_search.rs` 1,602; `koblitz_relation_solver.rs` 955; `koblitz_pdp_phase_a.rs` 1,140; `koblitz_isogeny_cost.rs` 1,701; `koblitz_fast*.rs` about 2,600 — `cryptanalysis:suite/src/cryptanalysis/`.
- Summation polynomials: `semaev_decomp.rs` 2,612; `binary_semaev_s4.rs` 1,919; `semaev_leading_form.rs` 1,792; `semaev_sat.rs` 1,525; `semaev_higher.rs` 889; `semaev_corpus.rs` 863; `binary_semaev.rs` 577; `symmetrized_semaev.rs` — same directory.
- Groebner and linear algebra: `f4_gf2.rs` 1,757 (fast Boolean F4); `f4_fp.rs` 1,493; `f4_fp_tower.rs` 2,553; `f4_gpu.rs` 2,228; `inherited_f4.rs` 2,005 (the default engine); `matrix_f5_f2.rs` 697; `f4_batch.rs` 573; `groebner_f4.rs` 482; `crossbred.rs` 1,214; `gf2_elim.rs` 637; `sparse_macaulay.rs`; `pq_wiedemann.rs`, `pq_sparse_la.rs` — same directory.
- Frameworks and accounting: `ic_boundary.rs` 7,484; `ec_index_calculus.rs` 1,440 plus `_j0` 884 and `_curves` 531; `prime_orbit_index_calculus.rs`; `hyperelliptic_index_calculus.rs` 1,825; `ic_oracle_pricing.rs` 1,198; `ic_descent_degrees.rs` 1,067; `ic_run.rs`, `ic_engine.rs`, `ic_progress.rs`, `index_calculus_budget.rs`; the pluggable `ic_framework/{linalg,plugins,solvers,stages}.rs` — same directory.
- Weil descent and GHS: `ghs_descent.rs` 703; `ghs_full_attack.rs` 479; `diem_descent.rs` 513; `weil_charts.rs` 654; `descent_{algebraic,expansion,lowgamma,treewidth}.rs`, `gaudry_cubic.rs`, `gaudry_quartic.rs` — same directory. No file name contains "xedni", and grep found no xedni implementation.
- The `ca-ic` CLI (`suite/src/cli/tools/ic/`: bench, boundary, budget, corpus, descent, experiment, fixed, params, prime, rho, workflow) offers six PDP solvers. The documented solvers are: `enumerate` (|F|^(m-1) group ops per target); `pair-table` (meet in the middle over |F|(|F|+1)/2 pair sums); `groebner` (the Weil-restricted Semaev system reduced by matrix-F4); `sat` (CDCL with parity rows); `wdsat` (external Trimoska WDSat); and `mq-fes` (libfes-style, m = 2 only). The CLI accepts degrees n odd from 3 to 63 for the Koblitz family, summands 2 to 4, and GGMP subfield curves with `--subfield k` up to 8 — [cryptanalysis:suite/docs/ic/README.md](suite/docs/ic/README.md).
- Relation columns are signed Frobenius orbits merged by cofactor projection. The relation matrix is kept in reduced echelon form over Z/rZ as relations arrive, and a run stops as soon as the scalar is pinned, then checks `[d]G = Q`. Factor-base materialization is capped at 4,096 abscissae — [cryptanalysis:suite/docs/ic/README.md](suite/docs/ic/README.md).
- `ca-ic search` ranks factor bases by expected trials, `(columns+1+extra)/coverage`. With `--solve-cost-targets`, it ranks by expected trials × measured F4 word XORs instead. At icv1-f2m15-t275-b7f03703 the two orderings disagree by 22.41× over 12 verified logs. With the legacy family, `ic run --degree 15 --curve-a 1` collects no relation in 20,000 trials, while `search` finds bases that cover all 210 targets — [cryptanalysis:suite/docs/ic/README.md](suite/docs/ic/README.md).
- Engine audit: binary/Koblitz and `F_{p^k}` (k ≥ 3) curves are "where this repository already has measured end-to-end results (largest solved: a 48-bit subgroup of `icv1-f2m61-t158598901-ab42b6c5`)". For prime-field curves "no sub-exponential IC is known". The `icx` engine labels every curve with one of three regimes: `attack`, `scaled` or `estimate` — [cryptanalysis:suite/docs/ic/ENGINE_AUDIT_AND_DESIGN.md §1](suite/docs/ic/ENGINE_AUDIT_AND_DESIGN.md).

**F4 and GPU maturity (measured)**
- The fast Boolean F4 (`f4_gf2`, combinatorial monomial ranking) has verdict "PASS_LOCAL as a PDP-stage improvement". With identical answers, the Groebner oracle runs 10.0–10.3× faster inside complete degree-23 `ca-ic run` DLPs. It is 2.1–12.2× faster on the stage ladder and 190–280× faster on degree sweeps. Relation collection at degree 23 on 4 cores takes 1.08 s, down from 41.9 s. The end-to-end S and the rho ratio remain "unknown" — [cryptanalysis:experiments/f4-gpu-20260925/RESULT.md](experiments/f4-gpu-20260925/RESULT.md).
- On an RTX 5090 (Runpod), large-matrix elimination of the icv1-f2m13-t181-515ee569 degree-6 matrix (168,363 × 188,203) takes 5.7 s against 426 s on the host Xeon, and `dreg_sweep` drops from 18.8 s to 8.2 s. The batched per-node decider is correct but only as fast as 13 host threads, and the host's inherited engine is 4.4× faster than the device path. The guidance is "leave `--f4-backend` off" — [cryptanalysis:experiments/f4-gpu-runpod-20261006/RESULT.md](experiments/f4-gpu-runpod-20261006/RESULT.md); [suite/docs/ic/README.md](suite/docs/ic/README.md).
- Early F4/F5 campaign wins: prime-field F4 symbolic preprocessing for Semaev fb=3 runs 2.34× faster (0.988 → 0.423 ms) and for fb=4 1.99× faster. Native Boolean F5B runs 1.36–1.49× faster and the F5/M4RI hybrid with 12 variables 1.42× faster. The campaign "establishes local progress, not a global speed record". 90 round directories exist — [cryptanalysis:experiments/groebner-perf-20260924/README.md](experiments/groebner-perf-20260924/README.md).
- `groebner_compare/` is "solver-stage infrastructure" for Boolean systems of 1–12 variables, with an msolve/M5GB bridge and an ECC2K-17 five-way comparison (`ecc2k17_five_way*.json`). It "does not implement an ECC attack ... or establish a speedup" — [cryptanalysis:groebner_compare/README.md](groebner_compare/README.md).
- A separator-conditioned (treewidth) Boolean query prototype exists — [cryptanalysis:experiments/groebner-f6-20261008/README.md](experiments/groebner-f6-20261008/README.md).
- Root-index and Bloom-filter engineering on N53:
  - A conjugate-x Bloom gate avoids 97.574% of orbit-key computations, 7,918,512 of 8,115,360 — [cryptanalysis:experiments/a1-allhit-conjugate-bloom-20261009/README.md](experiments/a1-allhit-conjugate-bloom-20261009/README.md).
  - A 4 MiB canonical-root Bloom filter removes 99.254% of exact S3 table probes, with 3,154,661 root keys — [cryptanalysis:experiments/koblitz-s3-root-bloom-20261009/RESULT.md](experiments/koblitz-s3-root-bloom-20261009/RESULT.md).
  - Sparse cyclic-zero-gap canonicalization replaces the 53-rotation canonicalizer — [cryptanalysis:experiments/koblitz-s3-sparse-canonical-20261009/RESULT.md](experiments/koblitz-s3-sparse-canonical-20261009/RESULT.md).
- Other implemented side mechanisms:
  - A bielliptic quartic `quartic-ic` diagnostic for `v^4 = x^3+ax+b`, run at toy size p = 53 — [cryptanalysis:experiments/bielliptic-quartic/README.md](experiments/bielliptic-quartic/README.md).
  - A hyperelliptic-cover library mirror pinned to crypto commit `3b8664a4`, which "intentionally adds no solver" — [cryptanalysis:experiments/hyperelliptic-cover-infrastructure/README.md](experiments/hyperelliptic-cover-infrastructure/README.md).

### Inferences
- Maturity ranking:
  - Prime-field `Z_p^*` IC in C: production, CLI-exposed and tested.
  - Koblitz/binary summation-polynomial IC in Rust: research-grade, with complete toy DLPs up to a 48-bit subgroup.
  - F4/F5 Boolean kernels: heavily optimized stage tools.
  - GPU F4: built and verified, but no end-to-end win.
  - GHS/Weil descent, hyperelliptic and quartic: diagnostics.
  - FPGA IC: none.
- "Root finding" for decomposition shows up as the S3/S4 root index and pair-table meet-in-the-middle (MITM) on N53 and N83, not as a general univariate root-finder module. Upstream crypto has `RESEARCH_GPU_ROOTFIND_DECOMPOSITION.md` and `RESEARCH_SEMAEV_ROOTFIND_VS_GROEBNER.md` at the root (`ls crypto:`), so that line of work probably lives upstream.

### Gaps
- I did not read the 90 `groebner-perf` rounds one by one, so per-round outcomes beyond the README summary are uncovered.
- I did not check test coverage or CI status for each IC module.

## Which experiments and results exist (sizes, timings, rates, comparisons with rho), and what do the papers say?

### Takeaway
Every complete binary/Koblitz IC DLP in this repo is slower than matched Pollard rho. On Koblitz curves the phase-aware IC is 1.3–300× rho for n ≤ 83. At N53 the observed rho/IC wall ratio is 0.256–0.319, so IC is about 3–4× slower. The extrapolations to ECC2K-130 leave gaps of 23.6 to 50+ bits. The only IC that crosses rho is the prime-field multiplicative group near 56 bits. Both papers are anonymized research drafts that have not been submitted; neither was withdrawn.

### Cited Findings
**Papers (`papers/README.md` status table)**
- `papers/ecdlp-constant/`, "The Constant in Front of √n: Measured Pollard Rho, ECC2K-130 Arithmetic, and Index Calculus over Prime Fields", has status "research draft, anonymized, not submitted" — [cryptanalysis:papers/README.md](papers/README.md). Its claims:
  - The abstract says that on prime-field curves and over `F_{p^3}` the IC method "sits hundreds of times above rho at every size that fits". The constant in front of √n "did not fall below one". The multiplicative-group linear sieve "does cross rho near 56 bits" — [cryptanalysis:papers/ecdlp-constant/abstract.txt](papers/ecdlp-constant/abstract.txt).
  - §5 table at 56 bits: fb 669, 6,270 unknowns, 7,794 relations, sieve 0.013 s, LA 0.793 s. Rho is about 1.0 s of multiply-equivalents. At 60 and 62 bits, rho is about 4.1 s and 8.3 s against IC's 1.63 s and the 62-bit figure — [cryptanalysis:papers/ecdlp-constant/sec-indexcalc.tex lines 36–70](papers/ecdlp-constant/sec-indexcalc.tex).
  - ECDLP IC values of S: prime curve with a 23.4-bit subgroup, MITM m=3, S = 56.8, which is 14.5× rho. Random binary curve with 24.4 bits, S = 206, 89× rho. Koblitz curve with a 39-bit subgroup, orbit-column MITM, S = 58.8, 300× rho. Matched rho has S = 3.93, 2.31 and 0.20 (signed-Frobenius). The fitted exponents, for example signed-orbit MITM at 0.54, are "all ... above one half. No cell is below matched rho" — [sec-indexcalc.tex lines 94–107](papers/ecdlp-constant/sec-indexcalc.tex).
  - `E(F_{p^3})` MITM triple oracle for p = 271..2083 (n = 2^24.2..2^33.1): S = 1,224 → 3,805, which is 1,057× → 2,240× rho, with fitted exponent 0.69 against a predicted 2/3. Engineering cut C_3 from 4.7–4.9M to 879,000 F_p-mults (5.5×), which moved the relation-only crossover from 2^111 to 2^96. In the genus-3/4 table, p = 271 dense gives S = 2,696 (2,328× rho), and p = 2,083 dense gives 898 (528× rho) — [sec-indexcalc.tex lines 110–168](papers/ecdlp-constant/sec-indexcalc.tex).
- `papers/frobenius-orbits/`, "Efficient Signed-Frobenius Orbit Canonicalization and Phase-Aware Representations for Elliptic-Curve Index Calculus", has status "research draft, not peer-reviewed, not submitted, not posted", category Implementation, about 25 pages — [cryptanalysis:papers/frobenius-orbits/README.md](papers/frobenius-orbits/README.md). Its claims:
  - Least-rotation canonicalization takes 4.78 rotations per 83-bit word against 82 for the plain scan.
  - The 18-byte root record with a 7-bit phase gives an index 2.5× smaller than an open-addressed table.
  - A Γ-orbit Groebner base solves a complete n = 19 ECDLP with 6.85× less CPU.
  - Normal-basis canonicalization makes the signed-Frobenius rho reference 4.1× cheaper, "which reverses an apparent index-calculus win".
  - Orbit folding divides pair-table entries by 2n and changes total operations by +1.2% at one width and −12% at another.
  - "Nothing here is faster than Pollard rho: at every degree priced, n ≤ 83, the phase-aware index calculus is 1.3–300× the matched rho" — [cryptanalysis:papers/frobenius-orbits/sec-abstract.tex](papers/frobenius-orbits/sec-abstract.tex).

**Toy end-to-end Koblitz IC**
- N53 weight-three root index: three replayed one-target DLPs (runs v3, v4, v5). IC took 3,202.6 / 4,572.4 / 3,147.0 ms against rho's 1,023.2 / 1,460.1 / 804.3 ms, for an observed rho/IC of 0.319 / 0.319 / 0.256. 721 of 721 ordinary queries gave verified relations. The "controlled single-target online speedup is unknown" because the Mac host had no isolation receipt. A prospective panel of 512 ordinary queries succeeded 512/512 — [cryptanalysis:experiments/hamming-ic-e2e-20260929/N53_W3_ROOT.md](experiments/hamming-ic-e2e-20260929/N53_W3_ROOT.md); [README.md](experiments/hamming-ic-e2e-20260929/README.md).
- N83 weight-three geometry gate: 89,474 usable points and a four-sum support ceiling of about 1.10 per million uniform subgroup targets. This is geometry only and solves no N83 DLP — [cryptanalysis:experiments/hamming-ic-e2e-20260929/README.md](experiments/hamming-ic-e2e-20260929/README.md).
- N83 costed screen on curve `EC1N83Ckb1h2bcb59d56ad6` (r = 2417851639230796216685689, cofactor 4): every online phase cost, yield, rank and rho reference is `null`, so "every online speedup is unknown" — [cryptanalysis:experiments/n83-costed-candidate-screen-20261005/RESULTS.md](experiments/n83-costed-candidate-screen-20261005/RESULTS.md).
- N83 pair-index screen: with 7 summands, B = 4,036 per slot and 1,128,870 fixed prefixes are needed for 50% coverage. With 8 summands, B = 2,018 and 17,901 prefixes are needed — [cryptanalysis:experiments/n83-residual-search-20261006/PAIR_INDEX_SCREEN.md](experiments/n83-residual-search-20261006/PAIR_INDEX_SCREEN.md).
- N13 complete candidate `IC1N13Ckb1fb6PDP4satRCsampleLAgaussTDdescentISO0h6adc50caf2cb` (four-sum selector CNF, 20 paired toy DLPs). The catalog is "1,000 design proposals, not 1,000 measured attacks" — [cryptanalysis:experiments/ic-candidate-catalog/README.md](experiments/ic-candidate-catalog/README.md), [COMPLETE_N13.md](experiments/ic-candidate-catalog/measurements/2026-09-25/COMPLETE_N13.md).
- `ic-bench`: a Python end-to-end harness for toy Koblitz `y^2+xy=x^3+1` with a Macaulay/XL PDP. It prices work in calibrated "rps" units against `rho_operations = round(sqrt(pi r/2))` × ec_add — [cryptanalysis:experiments/ic-bench/README.md](experiments/ic-bench/README.md). `fb-search` predicts end-to-end cost from the factor base at n = 19 and 23. `fb-archive` stores content-addressed bases, including n = 131 `geomtraceu` l = 12 — [fb-search/README.md](experiments/fb-search/README.md), [fb-archive/README.md](experiments/fb-archive/README.md).

**Scaling and extrapolation to ECC2K-130 (all rho-negative)**
- `pdp-scaling` takes rho on the <−1,τ> orbits as 2^60.81 iterations, about 6,417 core-years. The per-PDP budget needed to tie rho is 2^12.2 for m=4 (l=29), 2^32.8 for m=5 (l=28) and 2^36.8 for m=6 (l=24). For m=3, linear algebra alone exceeds rho. Measured decomposition costs extrapolate to 2^58 (m=4), 2^84 (m=5) and 2^72 (m=6) group additions per point, gaps of +46, +51 and +35 bits. IC "needs a decomposition algorithm about 2^35 to 2^50 times faster than meet-in-the-middle at l ~ 25–30". The solvers used were CryptoMiniSat 5.15, msolve 0.10.1 and WDSat — [cryptanalysis:experiments/pdp-scaling/README.md](experiments/pdp-scaling/README.md).
- `frobenius-quotient-m4`: no nontrivial Frobenius-stable subspace exists at n = 131, because 2 is primitive mod 131. The orbit quotient on a nonlinear weight base raises the per-decomposition budget to 2^31.1 (m=4) and 2^40.6–2^52.9 (m=5..12). Measured m=4 on that base still misses, since MITM is floored at 2^65.5 per attempt — [cryptanalysis:experiments/frobenius-quotient-m4/README.md](experiments/frobenius-quotient-m4/README.md).
- `linearized-half-decomposition` grows more slowly, but the total is floored at 2^(2n/3) = 2^87.3 at n = 131, "about 30 bits above rho". The README notes prior art in Courtois (ePrint 2016/003) and Karabina (2015/319) — [cryptanalysis:experiments/linearized-half-decomposition/README.md](experiments/linearized-half-decomposition/README.md).
- `orbit-quotient-summation`: the phase-free predicate "looks algebraically generic". At degree 131 there are 2^24.087 usable orbit keys and μ4 ≈ 2^-96.866, which leaves the approach short by about 59.6 bits — [cryptanalysis:experiments/orbit-quotient-summation/README.md](experiments/orbit-quotient-summation/README.md).
- `pdp-degree-heuristics`: across 656 factor bases (7 families, 9 fields), the linearization excess `e = n − dim V − dim V^(2)` sets the m = 2 solving degree. "Two summands cannot beat rho for any factor" base — [cryptanalysis:experiments/pdp-degree-heuristics/README.md](experiments/pdp-degree-heuristics/README.md).
- `shifted-base-geometry`: giving each summand a different Frobenius-shifted base raises coverage with no extra columns. This is the known GGMP construction (§3.1, Thm 3.2). Both chained-S3 SAT variants timed out on every smoke query — [cryptanalysis:experiments/shifted-base-geometry/README.md](experiments/shifted-base-geometry/README.md).

### Inferences
- The papers repackage the experiments above as measurement and engineering contributions. Neither claims an attack. In the README's naming, the "status" is draft/not submitted for both, and there is no withdrawn paper.
- The relation-rate data are reported as queries per verified relation, for example 721/721 at N53. I found no single "relations per second" headline figure.

### Gaps
- I did not extract the full N13 rho-paired ratios or every ic-bench `history.csv` row.
- No public relation rate (relations/s) at ECC2K-130 scale exists, because no relation collection was run at n = 131.

## What is the ECC2K-130 work here, and what did it conclude?

### Takeaway
There are two strands. The first is a production rho client (GPU in `ecc2k130/`, an FPGA core, and a live campaign shared with crypto). The second is an IC and isogeny-class hardness study whose answer is negative: no curve in the isogeny class is easier than E0, the best cost stays 2^60.81 rho iterations, and IC gives "0 bits".

### Cited Findings
- README status note dated 2026-10-01: "This repository contains no ECDLP speedup, and none is claimed". All volcano/descendant experiments (PRs #62, #65, #66, #70, #116) are "null or negative". E0 needs 0.996× the attempts of its descendants at n = 19 (95% CI 0.976–1.017). "Groebner degree stays the same on every curve tested" — [cryptanalysis:README.md lines 8–32](README.md).
- `experiments/factorbase-census` (n = 11, 13, 17, 19): curve-derived full-density subspaces come out "at or slightly below a random set of the same density", so they give "no structural gain" — [cryptanalysis:README.md](README.md), [experiments/factorbase-census/README.md](experiments/factorbase-census/README.md).
- Isogeny-class report:
  - E0 and all 262 floor curves cost 2^60.81 (2^60.91 with the Bailey et al. walk). The 263-isogeny transport costs ≤ 2^15.6 F_q-mults. Levels p and 263p sit in [60.81, 64.33].
  - IC (summation polynomials): "no method family beats rho, closest non-generic model 23.6 bits slower". Codex's E0 d_reg advantage was "a basis-presentation effect (0 bits after minimising over presentations)".
  - GHS magic number m = 1 (E0) or m = 131 (floor curves, genus 2^130), worth 0 bits. The break-even genus would be about 273–314.
  - Source: [cryptanalysis:experiments/ecc2k130-isogeny-class/REPORT.md §1, §3, §6](experiments/ecc2k130-isogeny-class/REPORT.md).
- The gap audits G2 (Weil descent/GHS), G3 (IC cost model; `margins.json`: rho_E0 = 60.809 log2, generic bound margin −0.326 bits) and G5 (per-curve summation polynomial at n = 131, including b-dependence by degree) — [cryptanalysis:experiments/ecc2k130-isogeny-class/gaps/](experiments/ecc2k130-isogeny-class/gaps/).
- Volcano IC at n = 19: E0 plus 456 floor descendants. Predicted attempts span 813–847, with E0 at 828. A trace-zero-hyperplane subspace cuts attempts from 824 to 602 (−27%) on every curve. Modelled Groebner CPU per ECDLP is 404–418 s — [cryptanalysis:experiments/volcano-ic/README.md](experiments/volcano-ic/README.md).
- The IC harness `ecc2k130/runner/codegen/indexcalc.py` states: "Index calculus does not beat rho at m=131 today and nothing here changes that". Its goal is a measured growth curve for a Frobenius-reduced base. The `nagao*` decomposition codegen sits alongside — [cryptanalysis:ecc2k130/runner/codegen/indexcalc.py](ecc2k130/runner/codegen/indexcalc.py).
- The rho GPU client reuses crypto's packed GF(2^131) arithmetic and walks the live campaign's `sigma^j(R)+R` with distinguished weight 32 — [cryptanalysis:README.md ~line 360](README.md).

### Inferences
- In this repo, ECC2K-130 IC is closed as a negative result. The open residual is the uninstantiated p and 263p levels. Reaching them needs a char-2 Kani pipeline, which is not an IC question.

### Gaps
- I did not open the typeset PDF `ecc2k130_isogeny_class_hardness.pdf` or the red-team scripts in detail.

## Which leads were declared negative or killed, and why?

### Takeaway
Nearly every binary-curve IC lead is recorded as negative or as "stage-only, no end-to-end speedup". The reasons repeat across leads: the linear-algebra or PDP budget exceeds rho, there are no stable subspaces at prime n = 131, structure turns out generic, or an apparent win is reversed by a better rho baseline.

### Cited Findings
- Curve choice in the volcano/isogeny class: null — 0.996× (CI 0.976–1.017) — [README.md](README.md).
- Curve-adapted factor bases: no structural gain — [experiments/factorbase-census](experiments/factorbase-census/README.md).
- Subspace IC at ECC2K-130 for m = 3..6: gaps of +35 to +51 bits — [pdp-scaling](experiments/pdp-scaling/README.md).
- Frobenius quotient with m = 4: still above budget — [frobenius-quotient-m4](experiments/frobenius-quotient-m4/README.md).
- Linearized half: about 30 bits above rho — [linearized-half-decomposition](experiments/linearized-half-decomposition/README.md).
- Orbit-quotient-native summation: generic, 59.6 bits short — [orbit-quotient-summation](experiments/orbit-quotient-summation/README.md).
- Two-summand IC: cannot beat rho for any base — [pdp-degree-heuristics](experiments/pdp-degree-heuristics/README.md).
- Codex d_reg advantage on E0: a presentation artefact. GHS: 0 bits — [ecc2k130-isogeny-class/REPORT.md](experiments/ecc2k130-isogeny-class/REPORT.md).
- An apparent IC win was reversed once rho got normal-basis canonicalization (4.1×) — [papers/frobenius-orbits/sec-abstract.tex](papers/frobenius-orbits/sec-abstract.tex).
- GPU batched F4 decider: no faster than 13 host threads, and 4.4× slower than the inherited host engine — [f4-gpu-runpod-20261006](experiments/f4-gpu-runpod-20261006/RESULT.md).
- Chained-S3 SAT on shifted bases: timed out on every query — [shifted-base-geometry](experiments/shifted-base-geometry/README.md).
- Sage storage candidate: put on HOLD at a geometric mean of 1.002× — [sage-binary-storage/RESULT.md](experiments/sage-binary-storage/RESULT.md).
- Homogeneous-fraction block solver: a PDP toy with `candidate_id` deliberately null — [homogeneous-fraction](experiments/homogeneous-fraction/README.md).

### Inferences
- The surviving positives are all engineering wins inside one stage:
  - F4 10× faster at degree 23.
  - Bloom-filter screens removing 97–99% of probes.
  - The 6.85× Γ-orbit Groebner base at n = 19.
  - The Sage batch arithmetic `PASS_LOCAL` results of 1.5–2.4×.

### Gaps
- `experiments/nonfrobenius-ic/` and `factor-base-yield*/` contain no README, only run data, so I could not determine their verdicts.

## Which TODOs, roadmaps, open docs issues and future-work items concern IC?

### Takeaway
The open work falls into four groups: `icx` engine sequencing (gf3m/fpk fields, CUDA/Metal FES backends), four documented but unfixed defects, factor-base heuristic follow-ups, IC in `F_{p^k}^*`, and wiring stage kernels into complete queries.

### Cited Findings
- The `icx` sequencing plan has 8 steps: catalog, progress/estimate, run with regimes, `gf3m` + `fpk` curves, field dispatch plus the aarch64 path, CUDA/Metal FES backends, bug fixes, and CI on 3 architectures. The defects documented but not yet fixed:
  - Oakley Group 3 order.
  - `fes-f2` returns 1 solution for an empty system.
  - `pack_point` key collision at n = 63.
  - `buchberger-f2` panics at n = 25..26.

  Source: [cryptanalysis:suite/docs/ic/ENGINE_AUDIT_AND_DESIGN.md §2.5, §5](suite/docs/ic/ENGINE_AUDIT_AND_DESIGN.md).
- Review remaining gaps: #3 binary and `F_{p^k}` curves in `crax challenge solve` ("reachable only through the `ic` tools"); #5 "Index calculus in F_{p^k}^*", to make MOV and the singular node case subexponential; #6 Block Lanczos and double large primes for factoring — [cryptanalysis:docs/REVIEW-2026-09-25.md](docs/REVIEW-2026-09-25.md).
- pdp-degree-heuristics next steps:
  - Widen paired workloads at n = 19 and 23.
  - Add m = 3 cells beyond n = 13.
  - Test the excess rule beyond N = 18 with msolve.
  - Add a trace-zero sibling of the `n131_poly_*` profiles.
  - Fit an m ≥ 3 rule on the symmetric model.

  Source: [experiments/pdp-degree-heuristics/README.md "Next steps"](experiments/pdp-degree-heuristics/README.md).
- Unintegrated kernels:
  - groebner-perf round67, round69 and round43: "not yet connected to the complete query path".
  - Metal batched point maps: "the IC runner does not yet call this".
  - N83 full W4: "does not yet make a viable PDP solver".
  - Catalog interfaces that implicit PDP solvers "do not yet" support.

  Sources: [groebner-perf-20260924/round67/RESULTS.md, round69/README.md, round43/INTEGRATION.md](experiments/groebner-perf-20260924/); [sage-ic-campaign/metal-batched-point-maps-20260924/RESULT.md](experiments/sage-ic-campaign/); [hamming-ic-e2e-20260929/N83_FULL_W4.md](experiments/hamming-ic-e2e-20260929/N83_FULL_W4.md); [ic-candidate-catalog/README.md](experiments/ic-candidate-catalog/README.md).
- The Sage IC campaign roadmap is an experiment queue. Item 1 (native Frobenius) reached `PASS_LOCAL` at 2.15×/2.31×, with PR #72 open. Item 2 reached 1.85×/1.79× and item 3 reached 1.337×/1.258× — [experiments/sage-ic-campaign/ROADMAP.md](experiments/sage-ic-campaign/ROADMAP.md).
- An open measurement at n = 53: windowed scan time splits into subtraction 10.9 ns, bulk key 9.5 ns and filter 3.1 ns, and resolving the 12.3% of keys the filter admits is unresolved — [suite/docs/ic/perf/OPTIMIZATION_PLAN.md line 148](suite/docs/ic/perf/OPTIMIZATION_PLAN.md).
- No GPU measurement on a device in this tree for the ECC2K-130 GPU kernel: "device run has not yet been measured" — [ecc2k130/README.md line 33](ecc2k130/README.md).
- The ECC2K-130 p-level route needs an "unimplemented char-2 higher-dimensional (Kani) isogeny pipeline" — [experiments/ecc2k130-isogeny-class/REPORT.md](experiments/ecc2k130-isogeny-class/REPORT.md).
- Paper checklists ("owed before submission") are in each paper README. Neither paper has been submitted, and both are anonymized — [papers/README.md](papers/README.md).

### Inferences
- Planned but not done:
  - `F_{p^k}^*` IC.
  - A released binary-curve IC CLI (`crax`/C) — the C `ic` command handles prime fields only.
  - Integration of the GPU/Metal kernels into the online query path.
  - An n = 131 relation-collection run.
  - The xedni calculus (no trace anywhere).
  - Isolated-host replays turning the N53 and N83 exploratory timings into controlled claims.
  - The m ≥ 3 degree heuristics.

### Gaps
- I found no GitHub-issue list locally. Open issues on GitHub were not checked, since the brief required local files only.

## Capabilities here but absent in /home/user/crypto, or the reverse (light check)

### Takeaway
cryptanalysis is a downstream port plus productionization of crypto's IC research. The Rust suite was ported from `crypto@54e977ab` (PR #121). It then added a native C prime-field IC, the fast Boolean F4 kernel with GPU echelon, the Python ic-bench/fb harnesses and both papers. crypto keeps growing as the upstream research repo: its `koblitz_index_calculus.rs` is about 50% larger, and it holds the boundary scoreboard, the autolab and the research notes.

### Cited Findings
- Port: "Everything cryptanalytic in crypto@54e977ab is now in `suite/`", with 216 files merged and 12 conflicting hunks in `koblitz_groebner.rs`. Not ported: cryptopals, zk, tls13, bls12_381 and the PQC CLIs — [cryptanalysis:docs/REVIEW-2026-09-25.md](docs/REVIEW-2026-09-25.md); [README.md "The attack suite"](README.md).
- "the per-stage scoreboard (`BOUNDARY_TARGETS.md`), the frozen run records and the autolab control plane ... stay in the upstream crypto research repository" — [cryptanalysis:suite/docs/ic/README.md](suite/docs/ic/README.md). Confirmed present in `crypto:docs/ic/BOUNDARY_TARGETS.md` (`ls`).
- File-presence check by name (`find`):
  - Line counts: `koblitz_index_calculus.rs` is 21,321 lines in crypto against 14,243 here; `ic_boundary.rs` 10,269 against 7,484; `koblitz_sparse_la.rs` 4,944 against 2,601; `semaev_decomp.rs` 4,041 against 2,612.
  - Absent in crypto by these names: `f4_gf2.rs`, `f4_gpu.rs`, `prime_orbit_index_calculus.rs`, `src/indexcalc.c`, and the experiment directories `ic-bench`, `fb-archive`, `factorbase-census`, `pdp-scaling`, `frobenius-quotient-m4`, `linearized-half-decomposition`, plus both paper directories.
  - Present in both: `groebner_compare`.
  - Upstream-only items (crypto root): `RESEARCH_ECC2K130_IC_FEASIBILITY.md`, `RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md`, `RESEARCH_SEMAEV_ROOTFIND_VS_GROEBNER.md`, `RESEARCH_GPU_ROOTFIND_DECOMPOSITION.md`, `research/notes/index-calculus/` and `research/KOBLITZ_SUBFIELD_INDEX_CALCULUS_QUEUE.md`.
  - Source: `ls`/`find` over `crypto:`.
- Some recent experiments here explicitly copy candidates from crypto, for example the sparse canonicalizer from `crypto/examples/koblitz_orbit_dlp_slice_ic_fastindex4.rs` at crypto rev `8ab924b9` — [experiments/koblitz-s3-sparse-canonical-20261009/RESULT.md](experiments/koblitz-s3-sparse-canonical-20261009/RESULT.md). The hyperelliptic cover is mirrored from crypto `3b8664a4` — [experiments/hyperelliptic-cover-infrastructure/README.md](experiments/hyperelliptic-cover-infrastructure/README.md).
- The C library and README cite the WDSat unification note and the factor-base solve-cost note as upstream `crypto` documents (`research/notes/ecc2k130/RESEARCH_WDSAT_IC_UNIFICATION.md`, `research/notes/index-calculus/RESEARCH_FACTOR_BASE_SOLVE_COST.md`) — [suite/docs/ic/README.md](suite/docs/ic/README.md).

### Inferences
- The relationship is a port that has since diverged in both directions:
  - cryptanalysis holds the hardened, CI-checked engine, the measured-negative ECC2K-130 studies, the fast F4/GPU kernels and the papers.
  - crypto holds the newer and larger Koblitz IC core, the autolab, the scoreboard and the research-lead notes.
- Missing files may only have been renamed. For example, crypto may hold F4 kernels under other names, because the check was by file name only. Kernel features may be out of sync between the two repos, so re-porting newer crypto Koblitz IC code (about +7k lines) is a probable pending task. I found no explicit plan for it.

### Gaps
- I did not diff the module contents or check whether crypto has an equivalent of `f4_gf2`, the fast Boolean F4, under another name.
- I did not check whether crypto has merged back any cryptanalysis-native code.

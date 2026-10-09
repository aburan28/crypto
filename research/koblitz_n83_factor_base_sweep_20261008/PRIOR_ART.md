# Source and literature inventory

Inspected source revision: `c70c32d486a3ac7531fe27f7193d9f09caa58344`. The isolated native baseline is recorded in RESULT.md. Every capability below describes source, not a completed degree-83 experiment.

`SOLVER_GATES.md` distinguishes modules present in that inspected main snapshot from modules actually present on the isolated study branch. In particular, the WDSat, FES and double-large-prime files named below are absent from this branch, and their pinned implementations have width or system-shape limits that bar direct N83 execution. The native S4 SAT source is present but its full retained-base model and capacity have not been validated.

| Source | Existing work | Consequence for this design |
| --- | --- | --- |
| `src/cryptanalysis/koblitz_index_calculus.rs` | Frobenius-divisor bases, explicit orbit domains, SAT decomposition, wide arithmetic, pinned degree-83 curve constructors and modular LA | Separate the invariant-subspace and compact-orbit families; bind exact representations and subgroup order |
| `src/cryptanalysis/koblitz_factor_base_search.rs` | Yield/coverage scoring and bounded candidate search | Avoid selecting a base by yield alone; include rank, full preprocessing and holdout costs |
| `src/cryptanalysis/koblitz_symmetrised.rs` | Rational-torsion coordinates and symmetrized algebraic systems | Symmetry is a domain/model change with exceptional-point and lifting obligations |
| `src/cryptanalysis/koblitz_groebner.rs` | F4, splitting heuristics and batch split bits | Freeze branch policies, degree/node/memory bounds and UNKNOWN outcomes |
| `src/cryptanalysis/wdsat_oracle.rs` | ANF export, capacity derivation and external WDSat model checks | Bind source ANF and binary/configuration; never silently truncate allocations |
| `src/cryptanalysis/mq_fes.rs` | Moebius, Monica and derivative Gray enumeration, packed-equation limits | Degree-83 equation representation and full verification are explicit adapter gates |
| `src/cryptanalysis/ecbench_large_prime.rs` | Exact zero/one/two-large-prime binary IC with single-word field/point/order types | The idea is implemented, but this adapter cannot be assumed to cover degree 83 |
| `src/cryptanalysis/ecbench/` and `docs/ecbench/README.md` | Native matched workloads, exclusive phase accounting, L0/L2 grades, replay and paired intervals | Required integration path for future complete comparisons; construction data is a stage diagnostic |
| `examples/koblitz_orbit_dlp_fast.rs` | Compact signed-Frobenius bases, S3 four-sum root index, parallel guided rank | Existing online-after-setup results do not identify the minimum cold runtime |
| `experiments/koblitz-single-target-n83-20261006/` | Frozen a=1 public fixture, K=600 base, three online IC/rho observations, separate preprocessing records | Different 53-bit subgroup; historical walls shared host load and excluded setup |
| `RESEARCH_ECC2K130_IC_FEASIBILITY.md` and `docs/ic/PLAN_IC_ACCOUNTING_FIXES_20261007.md` | Higher-arity proposals, implicit-index limits and corrections to total-cost extrapolations | Keep counts, empirical fits, proposed arities and end-to-end evidence separate |
| `docs/ic/boundary_targets.json` | Fail-closed evidence, distinct timing classes and promotion gates | No promotion from this construction panel or a compatibility count |

Primary literature inspected on 2026-10-08:

- Galbraith, Granger, Merz and Petit, [On Index Calculus Algorithms for Subfield Curves](https://sacworkshop.org/SAC20/files/preproceedings/18-IndexCalculus.pdf): invariant factor bases and Frobenius symmetry. Its domain conditions motivate the degree-83 invariant-module gate.
- Trimoska, Ionica and Dequen, [Parity (XOR) Reasoning for the Index Calculus Attack](https://arxiv.org/abs/2001.11229): ANF-aware WDSat and branching preprocessing. Its solver-stage results are not imported as complete-pipeline improvements.
- Bouillaguet, Cheng, Chou, Niederhagen and Yang, [Fast Exhaustive Search for Quadratic Systems in F2 on FPGAs](https://perso.lip6.fr/Charles.Bouillaguet/static/publis/SAC13b.pdf): derivative Gray-code evaluation for quadratic Boolean systems. Hardware and equation-family assumptions remain explicit.
- [Index calculus with double large prime variation for curves of small genus with cyclic class group](https://arxiv.org/abs/math/0606607): large-prime variation in curve class groups. Transfer to the exact elliptic point/orbit representation is an implementation and empirical obligation, not an automatic consequence.

The native order test derives the degree-83 cyclotomic-block fact directly. This study proposes no new polynomial identity, so the repository's prime-field `identity.certificate/v1` polynomial-evaluation certificate does not apply to that finite modular-order computation. No new isogeny edge or admitted speedup changes the canonical scoreboard, progress timeline, curve graphs or performance-gains graphs.

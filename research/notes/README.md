# Research notes

Every research note in this repository lives here, grouped by theme.  The
experiment directories that the notes cite (`research/<topic>_<date>/`,
`experiments/`, `docs/ic/runs/`) stay where they are; a note is the
reading, an experiment directory is the frozen evidence, and
`docs/index-calculus-scoreboard.html` is the drawn ledger.  The reporting
convention every note follows is in [`AGENTS.md`](../../AGENTS.md).

Filenames keep their historical `RESEARCH_` prefix so that commit
messages, pull-request threads and code comments that name a note by
its old name still identify it.

## ecdlp-general — surveys and generic-group levers

| note | subject |
|:--|:--|
| [`RESEARCH_ECDLP_STATE_OF_THE_ART.md`](ecdlp-general/RESEARCH_ECDLP_STATE_OF_THE_ART.md) | The ECDLP: state of the art, 2025–2026 |
| [`RESEARCH_TORSION_AUXILIARY_INPUTS.md`](ecdlp-general/RESEARCH_TORSION_AUXILIARY_INPUTS.md) | Auxiliary inputs (Cheon) and torsion points against rho and index calculus |
| [`RESEARCH_REPRESENTATION_STRUCTURE.md`](ecdlp-general/RESEARCH_REPRESENTATION_STRUCTURE.md) | Where exploitable structure can come from: the transfer pattern, an R1–R5 admissibility test for candidate handles, and why murmurations fail it |
| [`RESEARCH_BENCH_LOG.md`](ecdlp-general/RESEARCH_BENCH_LOG.md) | Cryptanalysis research bench: empirical log |

## index-calculus — Semaev decomposition, factor bases, Gröbner, first-fall degree

| note | subject |
|:--|:--|
| [`RESEARCH_RESIDUAL_WALKS.md`](index-calculus/RESEARCH_RESIDUAL_WALKS.md) | Residual walks over partial decompositions; the reference thread for the boundary-table-ratio rule |
| [`RESEARCH_IC_BOUNDARY_LEDGER.md`](index-calculus/RESEARCH_IC_BOUNDARY_LEDGER.md) | The boundary ledger: prime, random-binary and Koblitz index calculus end to end in one unit against the generic floor and a counted rho (`ic boundary`, frozen runs under `docs/ic/runs/`); §10 is the Round-2 engineering ledger (folded pair tables, walk targets, exact ceiling, balanced base) with every first-round row kept as its before mark |
| [`RESEARCH_SEMAEV_DECOMPOSITION.md`](index-calculus/RESEARCH_SEMAEV_DECOMPOSITION.md) | Fast factor-base decomposition for binary Semaev `S₄` |
| [`RESEARCH_SYMMETRIZED_SEMAEV.md`](index-calculus/RESEARCH_SYMMETRIZED_SEMAEV.md) | Symmetrised summation polynomials (FGHR) |
| [`RESEARCH_HIGHER_SEMAEV.md`](index-calculus/RESEARCH_HIGHER_SEMAEV.md) | Higher-order Semaev polynomials over prime fields |
| [`RESEARCH_SAT_SEMAEV.md`](index-calculus/RESEARCH_SAT_SEMAEV.md) | SAT-encoded binary Semaev systems |
| [`RESEARCH_INDEX_CALCULUS_FACTOR_BASE.md`](index-calculus/RESEARCH_INDEX_CALCULUS_FACTOR_BASE.md) | Structured factor bases over GF(p): four open fronts |
| [`RESEARCH_FACTOR_BASE_SOLVE_COST.md`](index-calculus/RESEARCH_FACTOR_BASE_SOLVE_COST.md) | Choosing the factor base for the solver, not for the yield |
| [`RESEARCH_GROEBNER_F4.md`](index-calculus/RESEARCH_GROEBNER_F4.md) | Minimal Buchberger / matrix-F4 solver |
| [`RESEARCH_DREG_MEASUREMENT.md`](index-calculus/RESEARCH_DREG_MEASUREMENT.md) | Solving degree vs first-fall degree on binary Semaev systems |
| [`RESEARCH_DESCENT_CROSSOVER.md`](index-calculus/RESEARCH_DESCENT_CROSSOVER.md) | The descent crossover: why a descended Semaev system's overdetermination and its decomposition yield are one parameter, and the scoping rule that follows for every refutation measurement |
| [`RESEARCH_FFD_MEASUREMENT.md`](index-calculus/RESEARCH_FFD_MEASUREMENT.md) | First-fall-degree measurement |
| [`RESEARCH_FFD_PROOF_COMPLEXITY.md`](index-calculus/RESEARCH_FFD_PROOF_COMPLEXITY.md) | A proof-complexity bridge for the first-fall-degree assumption |
| [`RESEARCH_FFD_WORKFLOW.md`](index-calculus/RESEARCH_FFD_WORKFLOW.md) | FFD falsification-driven experiment loop |
| [`RESEARCH_DEGREE_REDUCTION.md`](index-calculus/RESEARCH_DEGREE_REDUCTION.md) | Reducing the solving degree |
| [`RESEARCH_DIEM_DESCENT.md`](index-calculus/RESEARCH_DIEM_DESCENT.md) | Diem-style index calculus on `E/F_{p^k}` |
| [`RESEARCH_QUASI_SUBFIELD.md`](index-calculus/RESEARCH_QUASI_SUBFIELD.md) | Quasi-subfield polynomials over `F_{2^n}` |
| [`RESEARCH_HYPERELLIPTIC_IC_RHO.md`](index-calculus/RESEARCH_HYPERELLIPTIC_IC_RHO.md) | Index calculus vs rho on genus-2 and genus-3 Jacobians |
| [`RESEARCH_EXOTIC_COORDINATES.md`](index-calculus/RESEARCH_EXOTIC_COORDINATES.md) | Exotic coordinates for point decomposition |
| [`RESEARCH_AUTOLAB_LOG.md`](index-calculus/RESEARCH_AUTOLAB_LOG.md) | Research AutoLab log |
| [`RESEARCH_PKM_TOWER_ORACLE.md`](index-calculus/RESEARCH_PKM_TOWER_ORACLE.md) | Design and pre-registration of the prime-field algebraic oracle (Petit–Kosters–Messeng towers): construction, framework integration, the one-generator bound, and the solver-axis falsification test |

## ecc2k130 — the ECC2K-130 campaign, Koblitz curves, binary Weil descent

| note | subject |
|:--|:--|
| [`RESEARCH_ECC2K130_IC_LITERATURE.md`](ecc2k130/RESEARCH_ECC2K130_IC_LITERATURE.md) | What the literature has for index calculus on ECC2K-130 |
| [`RESEARCH_ECC2K130_ROUTES.md`](ecc2k130/RESEARCH_ECC2K130_ROUTES.md) / [`_ROUTE_TARGETS.md`](ecc2k130/RESEARCH_ECC2K130_ROUTE_TARGETS.md) | Routes and their experiments |
| [`RESEARCH_ECC2K130_DECOMPOSITION.md`](ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md) / [`_TARGETS.md`](ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION_TARGETS.md) / [`_RUNS.md`](ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION_RUNS.md) | Point decomposition: plan, targets, runs |
| [`RESEARCH_ECC2K130_RR_SOLVER_PANEL.md`](ecc2k130/RESEARCH_ECC2K130_RR_SOLVER_PANEL.md) | A different decomposition solver on the real curve |
| [`RESEARCH_ECC2K130_RELATION_SWEEPS.md`](ecc2k130/RESEARCH_ECC2K130_RELATION_SWEEPS.md) | Homogeneous relation sweeps and what rank is worth |
| [`RESEARCH_ECC2K130_EXTENSION.md`](ecc2k130/RESEARCH_ECC2K130_EXTENSION.md) | Raising ECC2K-130 to an extension field |
| [`RESEARCH_ECC2K130_HYPERELLIPTIC.md`](ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md) | Hyperelliptic covers of ECC2K-130 |
| [`RESEARCH_GROEBNER_STAGE.md`](ecc2k130/RESEARCH_GROEBNER_STAGE.md) | Optimising the Gröbner stage of the decomposition oracle |
| [`RESEARCH_ISOGENY_CLASS_SEARCH.md`](ecc2k130/RESEARCH_ISOGENY_CLASS_SEARCH.md) / [`RESEARCH_ISOGENY_DEGREE_SEARCH.md`](ecc2k130/RESEARCH_ISOGENY_DEGREE_SEARCH.md) | Searching the isogeny class for an easier Gröbner problem |
| [`RESEARCH_KOBLITZ_INDEX_CALCULUS.md`](ecc2k130/RESEARCH_KOBLITZ_INDEX_CALCULUS.md) | Frobenius-invariant factor bases on Koblitz curves |
| [`RESEARCH_KOBLITZ_SCALING_TARGET.md`](ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md) | Making the Koblitz decomposition oracle reach a useful `m` |
| [`RESEARCH_TRIMOSKA_BENCHMARKS.md`](ecc2k130/RESEARCH_TRIMOSKA_BENCHMARKS.md) | Bit-sliced Weil descent: EC-Index-Calculus-Benchmarks review |
| [`RESEARCH_WDSAT_IC_UNIFICATION.md`](ecc2k130/RESEARCH_WDSAT_IC_UNIFICATION.md) | Unifying Koblitz index calculus with the Trimoska WDSat solver |

## cm-isogeny — CM structure, isogeny graphs, covers, P-256 and secp256k1 audits

| note | subject |
|:--|:--|
| [`RESEARCH_SECP256K1_CM.md`](cm-isogeny/RESEARCH_SECP256K1_CM.md) | Structural audit of secp256k1's `j = 0` CM structure |
| [`RESEARCH_VOLCANO_FLOOR_RHO.md`](cm-isogeny/RESEARCH_VOLCANO_FLOOR_RHO.md) | Volcano-floor class-group-augmented rho |
| [`RESEARCH_MESTRE_HOWE.md`](cm-isogeny/RESEARCH_MESTRE_HOWE.md) | Mestre's algorithm for the explicit Howe-glued cover |
| [`PAPER_STRUCTURAL_COMPLETENESS.md`](cm-isogeny/PAPER_STRUCTURAL_COMPLETENESS.md) | Structural completeness of prime-field curves against isogeny-graph cryptanalysis |
| [`RESEARCH_P256.md`](cm-isogeny/RESEARCH_P256.md) / [`RESEARCH_P256_ISOGENY_COVER.md`](cm-isogeny/RESEARCH_P256_ISOGENY_COVER.md) | P-256 structural findings and the `(N, N)`-split-Jacobian cover search (PDF reports alongside) |
| [`RESEARCH_NIST_SOLINAS_STRUCTURE.md`](cm-isogeny/RESEARCH_NIST_SOLINAS_STRUCTURE.md) / [`_EXPERIMENTS.md`](cm-isogeny/RESEARCH_NIST_SOLINAS_EXPERIMENTS.md) | Cyclotomic structure in NIST Solinas primes |
| [`RESEARCH_PKM_CRITERION.md`](cm-isogeny/RESEARCH_PKM_CRITERION.md) | PKM-resistance audit across standardised curves |
| [`RESEARCH_EDS_RESIDUE.md`](cm-isogeny/RESEARCH_EDS_RESIDUE.md) | Elliptic divisibility sequences and elliptic nets |

## lattice-hnp — hidden-number problems, GLV leaks, LLL

| note | subject |
|:--|:--|
| [`RESEARCH_HNP_LANDSCAPE.md`](lattice-hnp/RESEARCH_HNP_LANDSCAPE.md) | Side-channel / HNP attack landscape for ECC |
| [`RESEARCH_GLV_HNP.md`](lattice-hnp/RESEARCH_GLV_HNP.md) / [`RESEARCH_GLV_HNP_PHASE2.md`](lattice-hnp/RESEARCH_GLV_HNP_PHASE2.md) | GLV-aware HNP on secp256k1 |
| [`RESEARCH_LLL_GS_ANALYSIS.md`](lattice-hnp/RESEARCH_LLL_GS_ANALYSIS.md) | Gram–Schmidt analysis of the secp256k1 LLL degeneracy |
| [`RESEARCH_CGA_HNC.md`](lattice-hnp/RESEARCH_CGA_HNC.md) | Class-group-amortised hidden-number cryptanalysis (formerly `docs/research/RESEARCH.md`) |

## code-based

| note | subject |
|:--|:--|
| [`RESEARCH_TII_MCELIECE.md`](code-based/RESEARCH_TII_MCELIECE.md) | TII McEliece key-recovery challenges |

## Adding a note

Put it in the theme it belongs to (or add a theme directory), add a row
here, and cite its frozen evidence directory by path.  Guides and primers
that are not research results go in `docs/guides/`; the library roadmap
is `docs/DEFERRED.md`.

A research note never lives at the repository root or loose in a code
directory: that is the one placement this index exists to prevent.  When a
note moves into its theme, update every inbound link in the same commit —
sibling notes, `docs/` pages, and the `docs/index-calculus-scoreboard.html`
citation — but leave the basename unchanged so the frozen experiment JSONs
that name the note by filename keep resolving.

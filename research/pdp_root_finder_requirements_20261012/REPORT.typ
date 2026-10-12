#set page(paper: "a4", margin: 1.7cm)
#set text(size: 9.5pt)
#set table(stroke: 0.4pt)
#show raw: set text(size: 8.5pt)
#show figure.caption: it => it.body

= Requirements for a dedicated root-finder for the point decomposition problem (prime and binary fields)

2026-10-12. Study directory `research/pdp_root_finder_requirements_20261012/` (the README is the full document; this report is its digest). Class: accounting. Nothing is measured here: every number is quoted from a frozen repository artefact or a cited paper, or is a calculation shown with its inputs in the README.

== The question

Scope the requirements for a new root-finding algorithm for the point decomposition problem (PDP) of Semaev-polynomial index calculus, specialised enough to be sped up, on prime fields and on binary fields; state how fast it must be against F4/F5 and general Gröbner-basis computation; and fix the baseline it is compared with. The request came with two pages of Petit–Kosters–Messeng (PKC 2016) §3: Algorithm 1 (the composed map $L = L_(n') compose dots compose L_1$ with $|L^(-1)(0)| approx p^(1\/m)$, the staircase System (4), and the sentence "we have not investigated the existence of any algorithm better than Gröbner basis algorithms to solve System (4)") and the "back door" remark on choosing standard primes through auxiliary curves of smooth order. The authors' own slides end on "dedicated algorithms to solve the systems: open problem". This document is the specification of that open problem as a measurable programme.

== What the repository already has, and what it measured

The PKM mechanism exists in the tree as the `prime-tower` factor base (Kummer, Dickson, 2-isogeny towers), the `pkm-tower` oracle (full, chain, tree presentations) and the engines `f4_fp` (level with msolve 0.6.5), `f4_fp_tower` and `sig_fp_tower`; the binary analogue is `descent-algebraic` with eleven `SystemSolver` engines, `semaev_decomp` (pairs-and-solve) and `mitm`. The measurements that decide what a dedicated root-finder must be:

#table(columns: 4,
  [track], [F4 solving degree (measured)], [growth of F4's per-target work], [generic oracles' growth (derived)],
  [binary, Koblitz $n = 11 dots 19$, unsplit $S_4$, $m = 3$ (split-dreg grid)], [7 at $ell = 4$ and $ell = 5$; Macaulay refutation degree 6 where the first-fall degree is 3], [$2^(8.6)$ – $2^(9.3)$ per rung of $ell$ (word XORs; two rungs, a ratio not a fit)], [enumeration $2^2$ per rung, pair table $2^1$],
  [prime, Kummer towers at $p_1$ (PKM rounds 2–7)], [5→6 at $m = 2$ ($N = 12 dots 22$); 6, 6, 7, 7 at $m = 3$ ($N = 9 dots 18$); 6, 7, 7 at $m = 4$ ($N = 8 dots 16$)], [1.6 / 2.1 / 2.4 bits per tower variable (multiply-adds)], [enumeration $(m-1)\/m$ = 0.50–0.75; pair table $floor(m\/2)\/m$ = 0.33–0.50],
)

F5/GVW signature criteria need rows 1–2 degrees higher than F4 on 18 of 18 prime systems and 1.7–28× the memory; F4 trace replay (Joux–Vitse) is bounded by the redundant rows: 4.3 % on the measured binary cell, at most 58 % of pair rows on generic $FF_p$ systems (Joux–Vitse's own figure: 2.1–3.4× their F4, up to 124× Magma); every SAT arm measured is enumeration with overhead (plain WDSat visits $2^(3 ell)\/3!$ nodes to within 0.2 %), and WDSat's published 60–380× is relative to Magma F4, itself above enumeration. On the exported $n = 17$ Boolean systems a full reduced basis costs msolve 42.7–46.5 s and PolyBoRi 40.4–52.3 s where the inherited F4 decides a target in 7–24 ms: a full basis is the wrong product.

#figure(image("figures/solving_degree.svg", width: 94%), caption: [Figure 2. Measured F4 solving degrees on both tracks (recorded values: split-dreg grid; PKM rounds) against the Macaulay degree the per-target budget admits at $m = 4$ (derived, README §6.4).])

== The requirements, condensed

*Functional (F1–F10).* Instances $(K, E, V, m, "presentation", "symmetry", x_R)$ on both tracks; a `SystemSolver` / `fp_solver_registry()` engine with a new registry id, registered in ecbench before any end-to-end row; verdicts `Decomposed(witnesses)` / `Refuted` / `BudgetExceeded`, witnesses group-verified, `BudgetExceeded` never a refutation; exactness against the exhaustive reference on the frozen corpora (digests equal, zero false refutations); seeded, deterministic counters, replayable; counted units (word XORs and field operations; $FF_p$ multiply-adds with the `ns_per_fp_muladd` calibration), `prepare` charged separately, wall time as a practicality note; Rust in-tree; ICV1 slugs; a frozen protocol before the first run; one design for both tracks with its field-specific parts named.

*Structural (S1–S8).* Membership consumed as structure (the tower or subspace as rewrite rules, or the group structure of $V$), never as one polynomial of degree $|V|$; the rational zeros of one function on the $2^N$-point grid $V^m$ are the target, not the ideal; no Macaulay-degree dependence in the exponent — any monomial linear algebra has a fixed degree within the admissible one, and the design names the mechanism that extracts zeros from the single $FF_p$ equation (one generator, Proposition 2 of the PKM note) or the $n$ descended equations; shape invariance exploited and charged (`prepare` once, lock-step batching, $(T_"prepare", T_"target")$ at $K in {1, 16, 256}$); symmetries stacked with claimed and measured factors; memory polynomial in $N$ or a declared, charged table; the track's obstruction named (no Frobenius-stable subspace at prime $n$; no subfield for $FF_p$); the inadmissible list of the PKM pre-registration.

== The gates

#table(columns: 3,
  [gate], [condition], [meaning],
  [G0 exactness], [digests equal the exhaustive reference; every witness group-verified; every `BudgetExceeded` listed], [a row may be read],
  [G1 vs enumeration floor], [work ÷ $T_"enum"(m, t) = 2^((m-1)t)\/(m-1)!$ below 1 and falling (≤ 0.6× per rung) over ≥ 4 rungs, uniform targets], [does algebra the enumerator does not],
  [G2 vs F4/F5 and general Gröbner (requested)], [(a) ≤ 0.1× the best Gröbner arm at the largest cell both finish; (b) lower fitted growth $beta$ with disjoint 95 % intervals over ≥ 4 sizes (to beat: 1.6 / 2.1 / 2.4 bits per tower variable; $2^(8.6)$–$2^(9.3)$ per rung); (c) not slower than msolve on exported systems], [faster in constant and exponent, not an artefact of our engines],
  [G3 vs the generic line], [$beta$ below the pair table's $floor(m\/2)\/m$ at poly($N$) memory, or $log_2(S T^2) - log_2(r\/|"Aut"|) <= 0$], [the only exponent gate: an advance],
  [G4 vs SAT (binary)], [fewer decisions than WDSat and the native CDCL+XOR on the frozen regression; decides $n = 59, ell = 9$], [the SAT line is not a competitor],
  [G5 amortised accounting], [$(S, T_"prepare", T_"target", K)$ and the $S T^2$ gap per rung; preprocessing-rho control at equal advice], [comparable with the ledger],
  [G6 end-to-end admissibility], [ECC2K-130: $T_"target" lt.tilde 2^35$ group ops at $m = 4$ (Macaulay degree 2); 256-bit prime: $T_"dec" lt.tilde p^(1\/4)\/32 approx 2^59$ at $m = 4$ (degree 3), $approx 2^76$ at $m = 5$; then ecbench claim / challenge / verdict], [bears on a deployed size],
  [G7 classification], [advance / engineering / relabelling / accounting], [only a G3 pass is a result],
)

#figure(image("figures/cost_growth.svg", width: 94%), caption: [Figure 3. Growth of F4's per-target work, indexed to the smallest cell of each line so that word XORs and multiply-adds share one dimensionless axis, against the enumeration and pair-table slopes (derived).])

== The baseline

Eight arms on five frozen corpora: R1 pairs-and-solve (PS-0: 337.7–452.3 ns/pair at $ell = 8 dots 12$), R2 `mitm` ($n = 59, ell = 9$: 80/80 in 3.86 core-s), R3 the repository's Boolean Gröbner engines (split-dreg grid), R4 WDSat / CryptoMiniSat / native CDCL+XOR (the 60-input regression), R5 msolve / Singular / PolyBoRi with Magma scripts written and not run, R6 `tower-exhaustive` / `mitm` / `subtract` on the `prime-tower` base, R7 `f4_fp_tower` (the PKM rows), R8 F4 trace replay. R8 and the prime R6 rows are to be built (about 3 and 1 person-days); a four-rung prime ladder per $m in {3, 4}$ is to be frozen. Procedure: AGENTS.md §8 and §10 (recorded baseline, A/A and ABAB ≥ 5 rounds on an isolated core, identical-output check), fits over ≥ 4 sizes with bootstrap intervals, the decision rule frozen before the candidate runs, one table with ratio columns to each boundary.

#figure(image("figures/gating_ladder.svg", width: 100%), caption: [Figure 1. The gate ladder with the two tracks' floors, reference arms, measured F4 behaviour and derived budgets.])

== Status and environment

Requirements and gates: complete as a document. Baseline: specified, not run. The algorithm itself: not attempted (out of scope). Measurements: none. Conductor's control plane was unreachable in this container (no Docker daemon or Postgres), so no task card is attached. `eprint.iacr.org` PDFs sit behind a Cloudflare challenge through this container's proxy and the PKC 2016 paper's host returned 503; the PKC 2016 slides, the CT-RSA 2011 Joux–Vitse paper and the arXiv WDSat text were read instead, and Amadori–Pintore–Sala is cited through the repository's PKM note. Graphs checked and left unchanged: the IC scoreboard, `progress-timeline.json`, the leaderboard and `docs/bounds/FRONTIER.md` (no measured $S$).

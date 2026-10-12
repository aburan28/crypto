#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9.5pt)
#set table(stroke: 0.4pt)

= lit2026 W5: lane-batched inherited-F4 replay, F4 cost-model check, gate-first F6-IC

2026-10-11. Branch `lit2026/w5-lane-batched-f4`, stacked on aburan28/crypto\#1681. Candidate commit `d22fb3b52`; program baseline `4b6e09b78`. Run reports: `reports/run-20261011T2225Z/report.pdf` and `reports/run-20261012T0008Z/report.pdf`. Host: shared Apple M4 Pro, aarch64, load 82–388; isolation *L0*.

== Question

Three ideas from the 2025–26 literature search and W4, applied to the production decomposition oracle `SolverEngine::InheritedF4 { max_degree: 3 }` (a splitting tree whose nodes read a support-local degree-3 Macaulay basis built at the root and specialised down the tree) and to F6-IC (the same tree with the exact geometry gate as node oracle):

+ *Lane-batched replay* (Demin, arXiv 2607.06372): replay one trace for k instances at once with packed coefficients. Over F#sub[2]: bit-slice k targets' matrices into lane words.
+ *F4-tracer cost model* (Kouba–Neiger–Safey El Din, arXiv 2603.16378): compare the repo's per-degree Macaulay rows and columns against the generic counts.
+ *Gate-first ordering* (W4 lead): ask the gate before the algebra a node would otherwise do for a refuted outcome.

== Answers

#table(columns: 3, align: (left, left, left),
  [task], [kind], [result],
  [1. lanes (`KIC_F4_LANES`)], [negative], [Exact (each lane is the unique RREF of its own root row space; tree identical to the RREF control), but 2.6–4.0× the control's stage unit on the three cells with a batch, 16–162× the sequential root's units in the root phase for k = 2…64, and slower in CPU time.],
  [2. cost model], [diagnostic], [Degree-3 root rows equal the generic Boolean count with occurring multipliers (0.67–0.70× with the production support-local ones); columns 0.28–0.49×; rank one below generic. Dense random systems of the same shape match generic exactly.],
  [3. gate-first (`KIC_F6_GATE_FIRST`)], [speedup], [F6-IC stage unit 1.57–1.78× lower than original on all four cells; process CPU 2.1–2.5× lower than previous (L0); verdicts, counters, oracle calls and the relation identical; F4 arm unchanged.],
)

== Why bit-slicing cannot pay over F#sub[2]

A packed row holds 64 coefficients per 64-bit word. A bit-sliced row over k lanes holds one column per lane word. A row operation on k packed rows of C columns touches $k C slash 64$ words; the bit-sliced one touches $C$ lane words of $k$ bits, the same number of bits. So at best it ties at $k = 64$ with every lane sharing every pivot, and it loses at any smaller k or any divergence. Demin's gain comes from packing small F#sub[p] coefficients that each occupy a word; over F#sub[2] that packing is already the representation. Measured: the lanes share their pivot row in only 3–14 % of pivot columns, and the lane root costs 16–162× the sequential root (lanes.svg in run 1). The splitting trees of different targets coincide in shape (94–97 % of F6-IC's stage unit is at nodes several targets reach), but their bases do not share pivots, so no tree lane engine was built.

#figure(image("reports/run-20261011T2225Z/arms.svg", width: 100%), caption: [run-20261011T2225Z, 5 interleaved rounds: stage unit original / arm, wall original / arm (L0), process CPU previous / arm (L0).])

== Gate-first

When a node's reduction forces variables, the default specialises its bases by each and asks the gate on the next pass. Gate-first assigns them (same order, same checks, same counters), asks the gate once, and specialises only if it continues. The gate reads only the assignment, so nothing it decides changes. Skipped per arm pass (3 repetitions): 447 / 432 / 840 / 405 nodes and 4 641 / 5 835 / 11 004 / 6 591 specialisations on K1/2#super[17], K1/2#super[23], K0/2#super[23], K0/2#super[31].

#table(columns: 5, align: (left, right, right, right, right),
  [cell (f6p)], [original = previous], [gate], [original / gate], [CPU previous / gate],
  [K1/2^17 dim 6], [12 635 914], [8 061 140], [1.5675], [2.22 [2.11–2.28]],
  [K1/2^23 dim 7], [51 918 433], [30 697 846], [1.6913], [2.49 [2.35–2.51]],
  [K0/2^23 dim 8], [124 949 801], [71 049 603], [1.7586], [2.24 [2.16–2.29]],
  [K0/2^31 dim 9], [278 237 634], [156 423 046], [1.7788], [2.13 [2.04–2.52]],
)

RUN2_SECTION

== Cost model

#figure(image("reports/run-20261011T2225Z/costmodel.svg", width: 100%), caption: [Degree-3 root Macaulay rows and columns, measured / generic Boolean count (rows $sum_j sum_(i <= D - d_j) binom(n, i)$, columns $sum_(i <= D) binom(n, i)$, rank from $(1+t)^n slash product_j (1 + t^(d_j))$ truncated).])

KNS state their counts over a field; the Boolean analogue is ours. The descended systems' sparse supports (mean 25–41 of 35–58 variables) halve the columns; support-local multipliers drop a further third of the rows; one syzygy appears at degree 3 in every cell. The tree's per-node step plans are smaller still (mean rank 134–1 036, columns 636–4 595).

== Identity

Per target, arm and round: hit, exhausted, reductions, splits, propagations, infeasible, geometric refutations and witnesses, partial refutations, group additions and pair-index counters identical across every arm and round of both runs; the relation of K1/2#super[17] target 3, base indices [21, 48, 57], identical across the candidate arms. Tests: `lane_roots_match_the_sequential_rref`, `w5_gate_first_and_lane_roots_keep_every_decomposition`.

== Limits

L0 only; wall time admits no claim (original/previous null control 0.24–4.60 per round); CPU time is L0 too. Four cells, m = 3, degree 3, aarch64 (lazy rewrite; the x86 eager rewrite is unmeasured). Only the roots are lane-batched. No IC-versus-rho statement; canonical graphs unchanged.

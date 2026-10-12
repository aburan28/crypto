#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#set table(stroke: 0.4pt)

= lit2026 W4: early abort at the refutation degree in GF(2) decomposition solves

2026-10-11 · branch `lit2026/w4-early-abort` · stacked on aburan28/crypto\#1681 · run report `reports/run-20261011T2025Z/report.pdf`

== Question

Horak and Semaev (arXiv 2603.24890) study how likely a sparse system over $FF_2$ is to be inconsistent. The program's reading (literature search, item 4) was: most index-calculus decomposition attempts have no solution, so the cost of proving $1 in I$ dominates, and an F4/F5 run should stop at the degree where $1$ appears. W4 asks where the production decomposition oracle decides its attempts, by degree, and whether an exact early abort at that degree saves work.

== What the production path is

The default decomposition oracle (`SolverEngine::InheritedF4 { max_degree: 3 }`, `koblitz_groebner.rs::solve_rec` and `reduce_inherited`) is not one Gröbner run per attempt. It is a splitting tree. Each node eliminates linear generators, then reads the linear tail of a support-local degree-3 Macaulay basis, and either refutes ($1$ in the row space), propagates forced variables, or splits. The basis is specialised from the parent's (inherited), or built from scratch after a linear elimination (fresh). On the descended systems of these cells ($m = 3$, chain with an auxiliary point, at least 24 variables) the degree ladder is degree 3 alone. A node whose basis holds $1$ is already closed at once: `reduce_inherited` stops its ladder on a decisive row, and an inherited child stops specialising at the row that makes it refuted. The degree-early-exit part of the idea was therefore already in the engine. W4 measures what is left.

== Method

- *Trace* (`RefutationDegreeTrace`, thread-local, off by default). At every reduction the node's own system is also reduced from scratch at degree 2, then 3. The node is filed under the first degree holding $1$ (refuted), else a forced variable (forced), else stalled. Refutations before any matrix are counted (constant, linear elimination, forced conflict, oracle). The stage-unit word operations each node spends are charged to "ends refuted" or "ends in a split". Root systems are tested at degrees 2–4. Probe: `examples/f4_refutation_degree_probe.rs`.
- *Policies* behind `KIC_F4_EARLY_ABORT`:
  - `1`: a fresh node first reduces its degree-2 space and is refuted if $1$ is there. Tree-preserving, because degree-2 support-local rows are rows of the degree-3 matrix.
  - `2`: a child that would specialise its parent's bases first reduces its own degree-2 space. Exact, because $1$ in a sub-ideal is $1$ in the ideal. The probe counts as the node's reduction.
  - `3`: both policies.
  - Pinned by `early_abort_tests` (roots, planted roots and tree counters identical across modes on 12 wide sparse systems; trace accounting).
- *Arms*, interleaved, 5 rounds, rotated order, `f6_ic_probe --reps 3`:
  - original: the `4b6e09b78` binary;
  - previous: `KIC_F4_EARLY_ABORT=0`;
  - candidates: `=1` and `=2`.
- *Gate*: `word_xors` (stage unit, deterministic). Wall is L0 only.

== Results

#figure(image("reports/run-20261011T2025Z/degrees.svg", width: 100%), caption: [Where the inherited-F4 tree decides its nodes (left) and where its stage unit goes (right), sums over attempts, `KIC_F4_EARLY_ABORT=0`, run-20261011T2025Z.])

#table(columns: 6,
  [cell (curve, subgroup of `KoblitzCurve::new`)], [F4 refutations: at deg 2 / deg 3 only], [share at deg 2], [forced (all at deg)], [stage unit in refuted nodes], [degree-2 probes / refutations, mode 1; mode 2],
  [K1/2^17 dim 6, `icv1-f2m17-tm101-00378d4e`, 8 att.], [203 / 305], [0.400], [189 (3)], [0.173], [254 / 10; 908 / 193],
  [K1/2^23 dim 7, `icv1-f2m23-tm5197-1f85e9e1`, 4 att.], [387 / 115], [0.771], [178 (3)], [0.038], [222 / 6; 924 / 381],
  [K0/2^23 dim 8, `icv1-f2m23-t5197-69e76b73`, 4 att.], [425 / 770], [0.356], [409 (3)], [0.203], [556 / 3; 2111 / 422],
  [K0/2^31 dim 9, `icv1-f2m31-tm90707-c95f16f5`, 1 att., budget 2000], [536 / 67], [0.889], [205 (3)], [0.0076], [279 / 2; 1061 / 534],
)

#figure(image("reports/run-20261011T2025Z/arms.svg", width: 100%), caption: [Arms, run-20261011T2025Z, 5 interleaved rounds each. Left: original / arm in the stage unit. Right: wall per round, L0, load 63.6–110.9.])

#table(columns: 5,
  [cell, engine], [original = previous `word_xors`], [cand1 (=1)], [cand2 (=2)], [ratio original/cand1; /cand2],
  [K1/2^17, f4], [15 015 618], [14 959 088], [15 058 684], [1.00378; 0.99714],
  [K0/2^23, f4], [158 512 775], [158 321 408], [158 846 162], [1.00121; 0.99790],
  [K0/2^31, f4], [281 583 026], [281 350 163], [281 787 203], [1.00083; 0.99928],
  [K1/2^17, f6p], [12 635 914], [12 636 832], [12 665 650], [0.99993; 0.99765],
  [K0/2^23, f6p], [124 949 801], [124 951 833], [125 135 729], [0.99998; 0.99851],
  [K0/2^31, f6p], [278 237 634], [278 238 729], [278 409 768], [1.00000; 0.99938],
)

*Identity.* Every per-target verdict (hit, exhausted) and every counter except `word_xors` (reductions, splits, propagations, infeasible) is identical across all four arms and all five rounds. The one relation found (K1/2^17, target 3, base indices [21, 48, 57]) is identical across modes 0, 1, 2. The original and previous arms agree to the word operation.

*Wall.* L0 medians and min–max are in the run report. The null control (original against previous, same tree and same stage unit) ranges 0.14–3.15 per round, wider than any arm difference, so no wall-time claim is made.

== Reading (interpretation)

- *The premise does not hold for this oracle.* Inconsistent nodes are most of the nodes (refuted reductions ≈ 42–43 % of reductions; every attempt here except one is a no-relation attempt). But they are cheap: 0.8–20 % of the stage unit at the four cells. The cost lies in nodes that do not refute (they propagate and split), and those are decided only at degree 3. No node was forced at degree 2, and no root system was refuted at degree 2 or 3. Capping the degree at 2 would lose all propagation; the attempt-level "refutation degree" is the depth of the splitting tree, not a Macaulay degree.
- *Mode 1* catches 2–10 fresh refutations per cell for 222–556 probes and saves 0.08–0.38 % of `word_xors` on the f4 arm. On f6p the change is within ±0.01 %. This is *no-change*: the switch stays off by default.
- *Mode 2* refutes 20–50 % of its probes, yet adds 0.06–0.29 %. An inherited child that holds $1$ already stops at the first row that refutes it, so the specialisation it skips costs less than the degree-2 build. This is *negative*: kept as a switch, off.
- *Lead, out of W4 scope.* In the F6-IC arm the F4 bases refute nothing: the geometry gate refutes the nodes. Yet 67–74 % of F6-IC's node word operations are spent in nodes the gate eventually refutes, that is, work done at a node before the gate's refuting call. The gate is already asked at node entry and after every propagation round, so these nodes are ones it can refute only once the algebra has forced more bits. This run does not say which part of that work (the entry specialisation, the reductions, or the specialisations by propagated variables) dominates. That is where an "abort before the algebra" could pay, for instance by a gate test on the forced bits before the bases are specialised by them. It would be a geometry-gate change, not a refutation-degree one.

== Limits

- L0 only.
- `word_xors` omits Macaulay row generation. Mode 1 also skips the row generation of the degree-3 builds it replaces (one or two builds per attempt), and that saving is not in the counter.
- Only $m = 3$, standard-subspace bases, `max_degree` 3, and the four cells shown. No F5 engine is on this path. W5's lane batching (inherited F4) is untouched: the change lives in `reduce_inherited` and `solve_rec` only.
- No IC-versus-rho statement.

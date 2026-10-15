#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#set heading(numbering: "1.")
#show table: set text(size: 8.5pt)

= Algebraic and SAT decomposition for $m >= 3$ against the meet-in-the-middle floor

*Stage diagnostic* · study `research/ic_m3_algebraic_vs_mitm_20261014/` · protocol frozen in commit `76e88621e` (2026-10-14, before the claim grid) · claim grid 2026-10-14 18:54–20:26 UTC · instrument `examples/ic_algebraic_floor/` · analysis `ANALYSIS.json` (`ic_algebraic_floor analyze`, 2,000 bootstrap replicates, seed 20261014).

== Question and answer

On binary Koblitz curves $K_a slash FF_(2^n)$ with factor base $F = {P : x(P) in V}$, $V = "span"(1, z, dots, z^(d-1))$, $d = ceil(n slash m)$, does any decomposition method that solves the summation-polynomial system on $V^m$ algebraically (Weil descent, then Gröbner or SAT) have a cost per *relation* that grows more slowly in $n$ than the pair-table (meet-in-the-middle) floor on the same cells?

*No, on the measured range.*
- *Grows faster than the floor*, with the 95% interval of slope(arm) − slope(floor) entirely above 0: F4 on the S#sub[3] chain and F6-IC, at both $m = 3$ and $m = 4$; and both unsplit-S#sub[4] hybrids at $m = 3$.
- *Undetermined* (fewer than 3 finished cells): SAT, at both $m$, and the S#sub[3]+S#sub[4] split hybrid. These censor too early for a fit; their lower bounds lie $10^6$–$10^7 times$ above the floor.
- *Amendment 1:* the wide `u128` chain solver reaches all seven $m = 4$ cells and also grows faster than the floor (§1.6).

The prediction made before the data held. The conditional end-to-end ecbench run (protocol §7) was therefore not triggered.

== Method (frozen protocol)

- *Cells.* $m = 3$: $n = 13, 17, 19, 23, 29, 31$. $m = 4$: $n = 11, 13, 17, 19, 23, 29, 31$. At each $n$, the Koblitz curve with the smaller cofactor.
- *Targets.* A pool of seeded uniform targets $[k]G$ per cell: 2,048 at $m = 3$, 1,024 at $m = 4$. The exact floor decides all of them, which fixes $p = P("decomposable")$ and the floor's cost.
- *Algebraic arms.* Each runs the same 6 decomposable and 6 non-decomposable targets. Every solve runs in its own process, on one thread, under a 120 s and 4 GiB budget. A solve that hits a budget is *censored*, and its cell becomes a *lower bound*.
- *Cost per relation:* $C = (p c_"yes" + (1-p) c_"no") slash p$, in wall-clock ns. The floor's cost includes the $|F|^2$ table build, spread over the $|F|$ relations a solve needs.
- *Slope.* OLS of $log_2 C$ on $n$. *Decision.* Δ = arm − floor, fitted on the same cells in the same bootstrap replicate.
- *Exactness.* Every finished verdict agreed with the floor (0 disagreements, 0 bad witnesses). The floor equals exhaustive enumeration on the 12 arm targets of every cell with $n <= 23$.

== Results

#figure(image("figures/cost_per_relation_m3.svg", width: 92%),
  caption: [$m = 3$: $log_2$ wall ns per relation against $n$, 95% bootstrap bars. Open markers with arrows are lower bounds (censored). Source: `cells.csv`, `ANALYSIS.json`.])

#table(columns: 8, align: right,
  [*n*], [*d*], [*p*], [*floor*], [*f4-chain*], [*f6-ic*], [*s4-hybrid-f12*], [*sat*],
  [13], [5], [0.494], [14.47], [23.56], [22.71], [30.12], [30.21],
  [17], [6], [0.274], [14.40], [26.12], [25.60], [32.77], [≥ 37.73],
  [19], [7], [0.404], [14.76], [27.28], [26.38], [34.24], [≥ 37.53],
  [23], [8], [0.247], [15.47], [29.55], [28.83], [≥ 38.64], [≥ 38.82],
  [29], [10], [0.279], [17.38], [32.66], [32.68], [≥ 38.81], [—],
  [31], [11], [0.479], [17.82], [33.08], [33.04], [—], [—],
)
$m = 3$, $log_2$ ns per relation. The cheapest arm divided by the floor is $3.0 dot 10^2$ at $n = 13$ and $3.8 dot 10^4$ at $n = 31$.

#figure(image("figures/slopes.svg", width: 88%),
  caption: [Slopes ($log_2$ ns per relation per unit $n$), 95% bootstrap, for arms with at least 3 measured cells.])

#table(columns: 6, align: (left, left, right, right, right, left),
  [*m*], [*arm*], [*slope*], [*floor, same cells*], [*Δ 95%*], [*verdict*],
  [3], [f4-chain (6 cells)], [0.534], [0.207], [[0.315, 0.338]], [grows faster],
  [3], [f6-ic (6 cells)], [0.583], [0.207], [[0.366, 0.384]], [grows faster],
  [3], [s4-hybrid-f12 (n ≤ 19)], [0.684], [0.040], [[0.608, 0.687]], [grows faster],
  [3], [s4-hybrid-f10 (n ≤ 19)], [1.320], [0.040], [[1.240, 1.321]], [grows faster],
  [3], [sat], [—], [—], [—], [undetermined],
  [4], [f4-chain (n ≤ 19)], [0.560], [0.156], [[0.342, 0.453]], [grows faster],
  [4], [f6-ic (n ≤ 19)], [0.573], [0.156], [[0.354, 0.469]], [grows faster],
  [4], [sat, s22-hybrid-g8], [—], [—], [—], [undetermined],
)
The floor's slope over all cells is 0.207 [0.204, 0.211] at $m = 3$ and 0.363 [0.342, 0.379] at $m = 4$.

*Post-hoc sensitivity (no verdict changes).*
- With the online floor (no table build), F4-chain's Δ is [0.233, 0.257] at $m = 3$ and [0.239, 0.296] at $m = 4$.
- Without the degenerate $m = 4$, $n = 11$ cell ($|F| = 3$), F4-chain's Δ is [0.41, 0.48].

#figure(image("figures/cost_per_relation_m4.svg", width: 92%),
  caption: [$m = 4$, main grid. For $n >= 23$ the chain system has $4d + 2n > 64$ unknowns, and the narrow engines cap Boolean systems at 64 (`MAX_VARS`), so only the floor ran there in this grid. Amendment 1 (§1.6) covers those cells with the `u128` solver.])

== Timeouts, censoring and limits

- *SAT.* At $m = 3$, refutations censor from $n = 17$ (3 of 6), then 6 of 6 at $n = 19$; the cell stopped at $n = 23$. At $m = 4$ it censors from $n = 13$ and stopped at $n = 19$.
- *Unsplit S#sub[4] hybrids.* Measured through $n = 19$; lower bounds at $n = 23$; stopped at $n = 29$.
- *Split hybrid.* Stopped at $n = 17$.
- *Not a cost.* A lower bound is never presented as a cost.
- *Scope.* Isolation level L0 (shared 4-vCPU container, the two grids concurrent), so wall time is descriptive. One seed and 12 targets per algebraic cell. Only the standard subspace base. No external solvers (msolve, CryptoMiniSat, Kissat) on this host. No symbolic binary S#sub[5], so no unsplit $m = 4$ system.
- *Range.* Slopes describe $n = 13 dots 31$ with these implementations; they are not asymptotic proofs.

== Interpretation (labelled)

- *F4 and F6-IC.* The repository's best algebraic decomposers grow at about $2^(0.53 n)$–$2^(0.58 n)$ per relation, against the floor's $2^(0.21 n)$ amortized, or $2^(0.29 n)$ online. The gap widens from $3 dot 10^2$ to $4 dot 10^4$ over the range.
- *What it would take.* Beating the floor's $|F| approx 2^(n slash 3)$ probes per target would need a solver whose per-target cost grows more slowly than that. None of these does.
- *Consistency.* This agrees with the earlier root-finder study (no slope below a floor) and with the SAT survey (SAT behaves as an enumerator). It adds the missing binary-field fit, per relation and on the same cells, for the Gröbner arms.

== Amendment 1: $m = 4$ with the wide `u128` solver (2026-10-15)

The owner asked for the 64-variable cap to be widened to 128 and the $m = 4$ grid rerun. The route chosen was the existing `u128` solver, `wide_groebner_decompose` (up to 128 unknowns), rather than a rewrite of the 64-bit engines.

It solves the same S#sub[3]-chain system. Its search is summand enumeration pruned by Macaulay linearisation, so it is a new arm, `wide-chain`. It ran on all seven $m = 4$ cells under the same rules (`AMENDMENT_1.md`, frozen in `882598e55`); the run is `reports/run-20261015T0024Z/`. All 70 finished verdicts agreed with the floor.

#table(columns: 6, align: right,
  [*n*], [*d*], [*unknowns*], [*floor*], [*wide-chain*], [*÷ floor*],
  [11], [3], [34], [18.09], [27.20 [26.45, 28.42]], [$5.5 dot 10^2$],
  [13], [4], [42], [14.99], [26.05 [25.81, 26.27]], [$2.1 dot 10^3$],
  [17], [5], [54], [16.96], [29.85 [29.62, 30.10]], [$7.6 dot 10^3$],
  [19], [5], [58], [18.37], [32.27 [31.99, 32.59]], [$1.5 dot 10^4$],
  [23], [6], [70], [21.06], [35.96 [35.45, 36.63]], [$3.1 dot 10^4$],
  [29], [8], [90], [21.68], [≥ 38.58 (censored)], [≥ $1.2 dot 10^5$],
  [31], [8], [94], [23.48], [≥ 40.52 (censored)], [≥ $1.3 dot 10^5$],
)

*Slope fits.*
- *Measured cells* ($n = 11$–$23$): slope 0.805 [0.717, 0.872], against the floor's 0.333 on the same cells. Δ is [0.433, 0.504]: *grows faster than the floor.*
- *Lower-bound fit* (all seven cells, valid): Δ [0.358, 0.390].
- *Without $n = 11$:* Δ [0.369, 0.410].

*Overlap with F4.* On $n <= 19$ the wide solver costs about what inherited F4 and F6-IC cost. It extends the measured $m = 4$ range to $n = 23$, and adds lower bounds at 29 and 31. The verdict is unchanged: no arm grows more slowly than the floor.

#figure(image("amendment-1/figures/cost_per_relation_m4.svg", width: 100%),
  caption: [$m = 4$ with the wide `u128` solver (orange), on the amendment run's floor. The original run's arms are faded. Open markers with arrows are lower bounds.])

#figure(image("amendment-1/figures/slopes.svg", width: 80%),
  caption: [Amendment 1 slopes, 95% bootstrap.])

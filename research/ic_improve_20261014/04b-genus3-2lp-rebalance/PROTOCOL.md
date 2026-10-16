# PROTOCOL: genus-3 IC-2LP rebalance (04b)

2026-10-15. Conductor T-1. Written and committed before any recorded run.
Changes after the first recorded run are dated amendments at the end of this
file. The text above them is never edited.

Tool: [`examples/ic_improve_04b_genus3_2lp.rs`](../../../examples/ic_improve_04b_genus3_2lp.rs)
(dependency-free Rust). Analysis: [`examples/ic_improve_04b_analyze.rs`](../../../examples/ic_improve_04b_analyze.rs).
Python only draws figures from the analysis output.

## 0. Provenance: what this study rebuilds, and what it cannot claim

The request is a follow-up to `research/ic_improve_20261014/04-genus3-2lp/`,
with its `PROTOCOL.md`, `research/ic_improve_20261014/RULES.md`,
`examples/ic_improve_04_genus3_2lp.rs` and `examples/ic_improve_04_analyze.rs`.
**None of these exist in this repository.** They were searched for on
2026-10-15 in:

- `main` and the working branch;
- every remote branch whose name carries 20261014 or 20261015, and `codex/ic-improvement-round3`;
- GitHub code search;
- pull-request search.

The user chose to have the baseline rebuilt (option 2 of the reply that
reported the gap).

Consequences, fixed here:

- The numbers in the request are the user's statement about the missing study.
  This study does not reproduce or check them: local exponents 1.25 to 1.46
  between 2^10 and 2^20, linear algebra at q^1.52, relation phase at q^1.30,
  row weight 7 to 25, F = ceil(2 q^{2/3}) tuned at 2^14. The rebuilt baseline is set
  beside them only qualitatively, labelled as such.
- "Same instances, seeds and counted unit" holds **within this study**: every
  arm runs on the same curves, targets, walk seeds and unit. It cannot hold
  with respect to the missing original.
- The user's rules are taken from the request:
  - this protocol is written before data;
  - every logarithm is verified;
  - memory cap 3 GB;
  - no job longer than 60 minutes;
  - `src/` is not edited;
  - work lives in this directory.

## 1. Question

Gaudry, Thomé, Thériault and Diem (GTTD) give a heuristic total cost of
Õ(q^{2−2/g}) for index calculus with double large primes on genus-g curves.
That is q^{4/3} for g = 3, with the factor base of size Θ(q^{1−1/g}) = Θ(q^{2/3}).

**Question.** In a rebuilt IC-2LP on genus-3 GHS covers, does the total cost
grow faster than q^{4/3} between q = 2^14 and 2^20? If it does, does the
departure go away when either of these changes:

- (i) the fixed factor-base constant κ = 2 is retuned at each size;
- (ii) the linear algebra after large-prime elimination is replaced by
  structured Gaussian elimination before Lanczos?

## 2. Instances

- **Sizes.** k = 10, …, 20, with q_k the largest prime below 2^k.
  - q_k: 1021, 2039, 4093, 8191, 16381, 32749, 65521, 131071, 262139, 524287, 1048573.
  - Factor-base sizes F = ceil(2 q^{2/3}): 203, 322, 512, 813, 1291, 2048, 3251, 5161, 8192, 13004, 20643.
- **Curves.** The hyperelliptic GHS cubic cover of type m = 3, r = 5, split arm.
  The construction is T-2463's, in `research/ghs_cubic_hyperelliptic_genus3_20261013/`.
  - Take K = F_{q³}, p ∈ K \ F_q, Q(x) = (x−a)(x−b) with a ≠ b in F_q.
  - The F_q-conic is Tr_{K/F_q}((σ²p − p) ξ²) = 0. Parametrize it by lines
    through a rational point. Then x = (r0 ξ² − r1 σ²(ξ)²)/(ξ² − σ²(ξ)²) and
    h(z) = (xn − a·xd)(xn − b·xd), with r_i = σ^i(p).
  - Move a rational Weierstrass point to infinity and scale. This gives the
    monic squarefree degree-7 model C : y² = f(x).
- **Acceptance.** A curve is accepted when all of these hold:
  - #Jac(C)(F_q) = N = h·ℓ with h ∈ {4, 8} and ℓ prime. This is T-2463's cofactor rule.
  - N is found by baby-step giant-step in q³ + 1 ± 2q^{3/2}. T-2463 verified
    the shape L(T) = 1 − tT³ + q³T⁶ (its P3) on its instances.
  - N is unique in that interval for the probe divisors.
  - N is confirmed by [N]D = 0 on four further random divisors.
  - N1 = #C(F_q) = q + 1 holds. The tool asserts this on every accepted curve.
- **Seed.** The base seed for recorded runs is **20261016**. Development runs used
  20261015 and are scratch: not committed, not reported as data.
- **Curve numbering and roles.** Curves are numbered 1, 2, … in acceptance
  order of the stream seeded by (seed, k). Curves 1–4 are **test**; curves 5–7
  are **training**.
- **Targets.** P = [h]D for a random divisor D, with P ≠ 0 and [ℓ]P = 0. The
  planted k is uniform in [1, ℓ−1], and Q = [k]P.
  - Test instances: targets 1, 2 and 3 on each test curve, so 12 per size.
  - Training instances: targets 1 and 2 on each training curve, so 6 per size and κ.
- **Walk seed.** The relation walk's seed depends only on (seed, k, curve, target).
  It is shared by every arm and every κ (common random numbers). Arms therefore
  differ only in what they change.
- **Transfer.** The transfer φ : E0(F_{q³}) → Jac(C)(F_q) is not executed per
  target.
  - T-2463 verified that φ is an injective group map on the ℓ-part, with image
    of index 2 or 4, both prime to ℓ. So Jac(C)[ℓ] = φ(E0(K)[ℓ]), and a target
    drawn in Jac(C)[ℓ] is a transferred target.
  - Curve construction, the group order and the transfer are inputs to the
    instance. They are not counted (§5).

## 3. Algorithm common to every arm

1. **Factor base.** Take the first F = ceil(κ q^{2/3}) x-coordinates
   α = 0, 1, 2, … for which f(α) is a nonzero square. There is one place per ±
   pair; the canonical y is min(y, q−y).
2. **Walk.**
   - Start at R0 = a0 P + b0 Q.
   - Use r = 32 multipliers M_i = α_i P + β_i Q.
   - Step R_{j+1} = R_j + M_{H(j)}, where H hashes the step counter j, not R_j.
     The walk is a fixed random sequence, so it cannot fall into a rho cycle.
3. **Decomposition.**
   - A Mumford u of degree ≤ 3 splits into distinct linear factors iff
     x^q ≡ x mod u, tested by square-and-multiply.
   - Roots come from Cantor–Zassenhaus.
   - A Weierstrass place (v(α) = 0) has order 2, prime to ℓ, so it is dropped.
   - The relation is *full* (no large prime), *partial* (1 or 2 large primes)
     or *discarded* (3).
4. **Large-prime graph.**
   - Vertices are the large primes plus a vertex 0; a partial is an edge (a
     1-LP partial joins vertex 0).
   - Collection stops when fulls + (E0 − V0) ≥ F + 10. Here E0 and V0 are the
     edge and large-prime-vertex counts of the component of vertex 0, and each
     non-tree edge of that component yields one relation.
5. **Logarithm.**
   - Find λ ≠ 0 over F_ℓ with Σ λ_j row_j = 0, then log = −Σ λ_j a_j / Σ λ_j b_j.
   - λ comes from Lanczos on the normal equations: A = M Mᵀ, solve A x = A y
     for random y, and take λ = x − y.
   - Mᵀλ = 0 is checked. If Lanczos fails, a new y is tried (at most 5).
6. **Verification of every logarithm.** The run checks [log]P = Q in Jac(C)
   and log = planted k. A run that fails either check counts as **failed**.
   Failed runs are listed and never dropped silently.

## 4. Arms

| arm | factor base | linear algebra | role |
|---|---|---|---|
| **BL** (rebuilt IC-2LP baseline) | κ = 2 at every size | **tree-cycle elimination + Lanczos**: a BFS spanning tree of the component of vertex 0 from vertex 0; one relation per non-tree edge (the edge plus tree paths, large primes cancelled with ±1 or ±2 scalars); matrix rows are fulls plus these relations, columns the factor base | reference for every ratio |
| **A1** (factor base retuned) | κ*_k ∈ {0.5, 0.7, 1, 1.4, 2}, tuned at each k = 14, …, 20 on training curves 5–7 | as BL | requested arm 1 |
| **A2** (cheaper linear algebra) | κ = 2; **the relation phase is BL's own record** (same run) | **structured Gaussian elimination + Lanczos**, details below | requested arm 2 |
| **A3** (both, secondary) | κ**_k tuned like A1 but on A2's pipeline | as A2 | secondary; outside the decision rule |
| **rho** (reference) | — | — | k = 10–14, test curves 1–4, target 1 |

**A2's elimination.**

1. Build the matrix from all fulls and all partials. Columns are the factor
   base plus every large prime.
2. Remove singletons repeatedly.
3. Prune excess: while rows − cols > 20, delete the heaviest row and redo
   singleton removal.
4. Merge in Markowitz order. Take the column with the least estimated fill,
   (w − 1)(L_p − 2) − L_p, where L_p is the lightest row in it, and pivot on
   that row. Keep merging while the Lanczos cost estimate
   cols · (2·nnz + 7·rows) falls; stop at the first merge that would raise it.
5. Run Lanczos on what remains.

**Tuning (A1, A3).**

- For each k from 14 to 20 and each κ in the grid, run the 6 training
  instances once. Each run does both linear-algebra routes on one relation phase.
- κ*_k is the arg min of the mean `total_mults` of BL's pipeline. κ**_k is the
  same for A2's pipeline. Ties go to the larger κ.
- The chosen values are written to `tuning/kappa_choice.json` by the analysis
  tool and committed before any test run of A1 or A3. They are not changed
  afterwards.
- If κ*_k = 2, A1 is BL at that size; the same deterministic run is reused, and
  the paired ratio is exactly 1.
- An optimum on the grid edge (0.5 or 2) is reported as such. The grid is the
  request's and is not widened.

**Rho.** An r-adding walk (r = 20) with Brent cycle detection, in the same
counted unit, on the same instances.

## 5. Counted unit

`total_mults` counts multiplications in F_q and multiplications mod ℓ:

- **Relation phase (F_q):**
  - setup: factor-base Legendre symbols, the start point and the 32 multipliers;
  - every walk addition;
  - every decomposition test and root finding.
- **Linear algebra (mod ℓ):**
  - large-prime elimination or merging;
  - Lanczos, including b = Ay, the inversions, the Mᵀλ check and the final division.
- **Verification (F_q):** the scalar multiplication behind [log]P = Q.

Counting rules:

- Squarings count as multiplications.
- Inversions are by Fermat, and their multiplications are counted.
- In elimination and merging, every scalar × entry product counts 1, including
  2 per row for the (a, b) pair, whatever the scalar's value.
- In Lanczos, every matrix-entry × vector-element product counts 1, whatever
  the entry.
- Additions and subtractions are not counted.
- Curve construction, group order and transfer are not counted (§2).

`S = total_mults / √ℓ`. Rho is measured in the same unit.

Wall-clock time is recorded and is never the metric. The host is a shared
cloud container with 4 vCPU and no isolation.

## 6. Metrics and analysis

For arm X, size k and phase φ ∈ {total, relation, la}:

- **Mean cost.** C_X(k) is the mean over the 12 test instances. The relation
  phase is `rel_mults`; la is elimination plus Lanczos.
- **Local exponents.**
  - Definition: e_X(k) = ln(C_X(k+1)/C_X(k)) / ln(q_{k+1}/q_k).
  - BL and A2: k = 10, …, 19. A1 and A3: k = 14, …, 19.
- **Trend.** T_X is the OLS slope of e_X(k) on k over the six local exponents
  k = 14, …, 19, per unit of log₂ q.
- **Fit.** E_X is the OLS slope of ln C_X(k) on ln q_k over k = 14, …, 20
  (7 sizes), compared with 4/3.
- **Paired ratio.** For each test instance and k, r = total_X / total_BL.
  Per size, report the geometric mean; also per phase.
- **Uncertainty.**
  - Nonparametric bootstrap with B = 10000 and seed 20261016.
  - Each replicate resamples, with replacement, the 12 test instances
    independently within each size.
  - The same resampled instances are used for every arm (paired), so
    differences between arms are paired.
  - Report percentile intervals at 95% and at 99%. The 99% interval is the
    decision interval.
- **Weight.** Report the row weight nnz/rows and nnz/cols of the matrix that
  enters Lanczos, per arm and size.

## 7. Decision rule (frozen)

**D0: does the rebuilt baseline depart?** BL departs if the 99% interval of
E_BL lies above 4/3, or the 99% interval of T_BL lies above 0. If BL does not
depart, the study reports that the rebuild does not show the departure. The
question then cannot be answered on this rebuild; the arms are still reported.

**D1: what each requested arm does** (A1 and A2 separately):

| verdict | condition |
|---|---|
| removes the departure | the 99% interval of E_X contains 4/3 or lies below it, **and** the 99% interval of T_X contains 0 or lies below it |
| reduces the exponent | the 99% interval of E_X − E_BL lies below 0, without meeting "removes" |
| constant only | the interval of E_X − E_BL contains 0, and the 99% interval of the paired ratio at k = 20 excludes 1 |
| no measured effect | both contain 0 or 1 |
| raises the exponent | the 99% interval of E_X − E_BL lies above 0 |

**Answer.** The departure is called an *artifact of choice (i) or (ii)* for
the arm that removes it while BL departs. It is *partly* an artifact if the
arm only reduces it, and *not* otherwise. Both arms are reported whichever
way this comes out.

**Multiple comparisons.** There are five decision tests: D0 and four D1
conditions over two arms. The 99% interval is Bonferroni over a family
α ≈ 0.05.

**Power check.** If the half-width of the 95% interval of E_BL exceeds 0.03,
the study cannot resolve the smallest departure of interest. That departure
is 1.37 − 4/3 ≈ 0.04, below the 1.46 of the request. In that case, verdicts
whose intervals straddle a boundary are reported as "inconclusive at this
resolution".

## 8. Controls

- **Negative control (pairing).** A2's relation phase must be BL's own record:
  the ratio of relation-phase costs is exactly 1 on every instance.
- **A/A control.** Fit E_BL on target 1 alone and on targets 2–3 alone. The
  95% interval of the difference should contain 0; report it either way.
- **Positive control (mechanism visible).** The baseline row weight is
  measured at every size. The study reports whether it grows (the mechanism
  named in the request) and whether A2 bounds it.
- **Accounting check.** `lanczos_mults` must lie within 2% of
  iters · (2·nnz + 7·rows) plus the b = Ay product, the Mᵀλ check and the
  inversions, on every run. Any violation is reported as an accounting failure.
- **Rho reference.** S_rho is measured at k = 10–14 and set beside each arm's S.

## 9. Execution limits

- **Jobs.** One job is one curve generation for one k, one solve
  (k, curve, target, κ, with both linear-algebra routes), or one rho run.
- **Limits.** Each job runs under `timeout 3600` and `prlimit --as=3 GiB`.
  Four jobs run in parallel.
- **Failures.** A job that hits a limit or fails a check is recorded with its
  reason in `results/failures.jsonl`. It is not rerun with other parameters,
  and the affected arm and size are marked incomplete.
- **Order.**
  1. Curves for k = 10, …, 20.
  2. Tuning runs for k = 14, …, 20.
  3. The κ choice, committed.
  4. Test runs: BL and A2 at κ = 2 for k = 10, …, 20; A1 at κ*; A3 at κ**.
  5. Rho.
  6. Analysis.

## 10. Pre-flight estimate (development seed 20261015, scratch, not data)

One instance per row, curve 1, target 1, single-threaded:

| k | κ | rel. mults | BL la mults | BL total | A2 total | BL row weight | wall (s) |
|---|---|---|---|---|---|---|---|
| 14 | 2 | 6.06e7 | 3.65e7 | 9.71e7 | 7.68e7 | 7.60 | 0.6 |
| 16 | 2 | 3.99e8 | 2.75e8 | 6.75e8 | 5.57e8 | 9.64 | 4.3 |
| 18 | 2 | 2.36e9 | 2.02e9 | 4.38e9 | 3.80e9 | 11.60 | 41 |
| 18 | 0.5 | 7.86e9 | 2.04e8 | 8.06e9 | 8.04e9 | 20.71 | 45 |
| 20 | 2 | 1.53e10 | 1.52e10 | 3.05e10 | 2.77e10 | 14.36 | 361 |
| 20 | 0.5 | 5.26e10 | 1.58e9 | 5.41e10 | 5.39e10 | 25.99 | 304 |

Peak memory was 119 MB, at k = 20 with κ = 2. Every job is projected below
the 60-minute and 3 GB limits.

Estimated CPU time:

| part | runs | CPU time |
|---|---|---|
| tuning | 210; 30 at k = 20 at about 5.5 min each | about 4 h |
| test | at most 36 relation phases at k = 20 | about 5 h |

On 4 parallel jobs that is about 2.5 to 3 hours of wall time.

## 11. What does not count

The following are new, labelled protocol versions, and both versions are
reported:

- changing, after seeing test data:
  - the κ grid, the training curves or the instance counts;
  - the stop margin F + 10, A2's parameters (excess 20 and the cost estimate)
    or the counted unit;
- choosing seeds;
- dropping failed runs;
- refitting over a different size range.

Rho, the A/A control and the arms' phase exponents are descriptive. They do
not enter the decision rule.

## Amendments

(none yet)

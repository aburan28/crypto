# Joint factor-base dimension / decomposition-length pilot

This round explores the higher-dimensional-factor-base idea from the supplied
textbook pages. It varies support dimension **ell** and number of summands **k**
independently, measures actual uniform-target coverage, and charges a complete
tiny ECDLP pipeline. It does **not** implement a general-k Nagao polynomial
solver or establish an asymptotic improvement.

The frozen [contract](contract.json), [targets](targets.json),
[comparison](comparison.json), and [complete results](RESULTS.md) are the
reproducible baseline for subsequent iterations. No old benchmark was replaced.

## Outcome

- 33 parameter cells: 1,584 generic stage calls, 504 matched algebraic controls,
  132 cold IC attempts, and eight matched rho references.
- All stage cells finished within their two-second budget and matched the exact
  domain oracle. Of 132 IC attempts, 86 verified the planted scalar; 46 retained
  failures comprise 24 zero-column cases, 12 insufficient-support cases and ten
  rank-deficient cases. All eight rho scalars verified.
- 774 exhaustive five-bit solver checks, four GF(8)-space closure checks, and
  4,195 independently checked full-run decomposition calls, including negatives.
- Four post-run audit tests also pass: exact DP multiplicities against direct
  signed enumeration, the subfield projection obstruction, timeout preservation,
  and actual-curve scalar verification for all 94 completed IC/rho runs.
- The previous full 240-stage / 48-DLP benchmark was replayed in `regression/`.
  Its comparator reproduced all 11 completed stage-baseline and six cold-baseline
  counter checks, and audited 384 oracle calls. Previous measured source files
  are unchanged. The compatible a=0 suite is used because the WDSat a=1 protocol
  is not the encoding or curve exercised by these implementations.

Classification: **accounting/coverage pilot, with component-wise engineering
candidates**. No calibrated S, rho ratio, floor ratio, runtime speedup or exponent
claim is established. Field and scalar arithmetic and Python control have no
measured common conversion here.

## What we learned

1. **Success probability is not monotone in k in this domain.** At binary5,
   ell=3, exact coverage is 28/43, 33/43 and 29/43 for k=2,3,4. Distinct nonzero
   abscissae and the target-abscissa exclusion matter. The paper intensity is not
   an exact probability, especially on tiny structured supports.
2. **More coverage can still cost more overall.** With pair-MITM at binary9,
   k=3, moving ell=4 to ell=5 raises exact coverage from 440/507 to 507/507,
   but projected columns grow from ten to 16. Across both verified seeds,
   field multiplications rise from 44,407 to 62,288 and scalar multiplications
   from 893 to 3,674. The reverse change passes the predeclared field-vector
   gate on each seed, with scalar nonincrease as well.
3. **Cold query and cold DLP rankings differ.** At binary9, ell=4,k=3, pair-MITM
   requires more field multiplications than enumerate-last on both twelve-target
   cold stage panels (58,672 vs 42,836 development; 59,068 vs 47,451 holdout).
   Across two full cold DLPs, where each run builds its table once and reuses it
   for its own relation/individual calls, field multiplications fall from
   94,192 to 44,407. Additions fall from 39,995 to 21,184; squarings remain
   13,708; scalar A/M/I remain 745/893/22. This is charged within-run reuse,
   not precomputation shared for free across target DLPs.
4. **A subfield factor base can be structurally useless.** For both GF(8)
   profiles, ell=1 yields zero nonidentity projected columns. For b=1,
   #E(GF(8))=4 and h=508/127=4; for b=22, #E(GF(8))=12 and h=468/13=36.
   These cofactors annihilate the entire subfield group. Rational abscissae
   over GF(8) cannot acquire new quadratic y-roots in the odd-degree extension
   GF(512)/GF(8). Moving to ell=2 yields 26 or six usable columns respectively;
   every k=2,3,4 / solver / seed full-DLP case then completes.
5. **MITM is not universally better.** At binary5, ell=4,k=3, its two-run
   field multiplications rise from 5,695 to 9,018 because setup dominates.
   At binary9, ell=4,k=4, they fall from 130,058 to 49,497, but scalar
   multiplications rise from 1,014 to 1,707 because the first relation selected
   can differ. All such regressions remain in the comparison.

The frozen neighbor gate passes 21 directed parameter changes on the measured
field vector; only 13 also have scalar nonincrease on each seed. These are
small-sample screening signals, not 21 discoveries or a fitted crossover.

## Domains and boundaries

All curves have equation `y^2 + x*y = x^3 + b` in GF(2^bits). The b=1 profiles
are the same a=0 Koblitz family as the prior pullback suite. In the normal
basis, one has all bits set. The fourth profile uses b=22, an embedded GF(8)
element outside GF(2), and is not labelled Koblitz.

| Profile | q | n | Ambient bits | b coordinates | #E | Prime subgroup N | ell |
|:--|--:|--:|--:|--:|--:|--:|:--|
| binary5 | 2 | 5 | 5 | 31 | 44 | 11 | 2,3,4 |
| binary9 | 2 | 9 | 9 | 511 | 508 | 127 | 2,3,4,5 |
| subfield8-k0 | 8 | 3 | 9 | 511 | 508 | 127 | 1,2 |
| subfield8-b22 | 8 | 3 | 9 | 22 | 468 | 13 | 1,2 |

For each dimension, k ranges independently over 2,3,4. GF(2) supports are spans
of the first ell normal-basis vectors. GF(8) supports start with the subfield
itself and extend by the first deterministic independent element. Their exact
sizes and scalar closure are checked. These are specific nested spaces, not a
search over all bases or subspaces.

Let M be the **actual** number of nonzero rational abscissae in the support.
Both point signs are allowed. We require k distinct abscissae, none equal to
the target abscissa. A counting ceiling on uniform affine-target coverage is

`min(1, 2^k * binomial(M,k) / (#E - 1))`.

The numerator bounds signed unordered decompositions before target exclusion;
collisions, sums at infinity and exclusions can only reduce coverage. The
paper estimate `q^(k*ell-n)/k!` is a sparse random-sum **intensity**, not a
rigorous probability or a statement about usable subgroup matrix rank.
`1-exp(-actual mean representations)` is also only a random-sum prediction.
An independent group-coordinate dynamic program gives exact representation
counts and covered target counts here; measured solvers never receive that map.

The predeclared engineering gate compares both fixed seeds separately for each
neighboring parameter change. Every field component must not increase, and at
least one must decrease. Scalar nonincrease is additional reported evidence.
No uncalibrated vector is turned into `S=operations/sqrt(N)`. The matched rho
reference includes the identical tiny group/order setup and the same target.
Ambient bit length is not subgroup bit length.

## Algorithms actually exercised

- **Enumerate-last:** enumerate k-1 signed distinct factors in canonical index
  order, then look up `R - partial_sum` in the single-point table. Worst-case
  query work O(M^(k-1)); factor-base/table storage O(M), for fixed k.
- **Pair-MITM:** retain every collision bucket of signed two-point sums.
  k=2 uses singles/singles; k=3 uses singles/pairs; k=4 uses pairs/pairs.
  Canonical index separation prevents repeated factors and permutation repeats.
  Pair precomputation/storage is O(M^2); query probes are O(M) for k=3 and
  O(M^2) for k=4, **plus collision-bucket filtering**, not a collision-free
  complexity guarantee. Every returned sum is checked on the curve.
- **Algebraic controls:** cached coefficient pullback, symmetric S4 and chained
  S3 through the existing adapters, only at q=2,k=3 on identical targets/support.
  Each control has 168 cells (87 found, 81 proven empty). Returned triples
  are checked with an independent exact pair oracle. Existing SAT diagnostics
  remain in raw adapter results; missing common operation units remain null.
- **Full DLP:** order/generator validation, cofactor projection, incremental
  dense modular row reduction, individual descent, final `[scalar]G=Q` check.
  This tiny matrix control is not sparse Wiedemann or a new F4 implementation.

For relation `R=s*A=sum P_i`, projection by h gives `s*G=sum h*P_i`, with
`G=h*A`. Individual `R=Q+s*A` gives projected logarithm `h*secret+s` modulo N;
the implementation checks `gcd(h,N)=1` before division. Zero projected columns
are rejected explicitly. Dependent rows are not counted as independent rank.

## Accounting, correctness and limits

Each cold full run charges group enumeration/order checks, generator and target
construction, support/table construction, failed searches, projection, exact
verification, modular matrix work, individual recovery and final verification.
The full-group generator algorithm requires a cyclic toy group and an invertible
cofactor for the selected prime; other groups are outside this pilot. Group
enumeration is intentionally expensive and not a scalable group-order method.

Field additions/multiplications/squarings are exclusive API counts. Binary-field
inversion work expands into these operations; `inversionCalls` is diagnostic.
Inherited curve-call counters are **not instrumented in the new generic curve
adapter**; their zero values do not mean no curve arithmetic. Scalar modular
A/M/I and search/table/bucket counters remain separate. Conversion, hashing,
integer order factoring and Python control are not priced in a common operation
unit; elapsed cold time includes them. Exclusive field phases and full-run wall
phases are preserved; interrupted-phase overhead has its own wall bucket.

Validation-only group-log maps, coverage DP and independent audits are outside
timed algorithm costs and are never solver inputs. Census setup cost is recorded
separately. The measured solvers themselves also verify every returned sum.
`compare.py` checks source hashes, exact cell-set completeness, same-target rho
matching, exclusive field totals, phase-time sums and paired completion status.

The 24 affine targets per profile were frozen uniformly without replacement,
with twelve development and twelve disjoint holdout points. No supported-target
conditioning or retuning after holdout results. Binary9 and subfield8-k0 use the
same curve and targets, so they are not independent curve samples. The cold DLP
uses two new fixed seeds, not a statistical performance distribution. One stage
repetition and these tiny subgroups do not justify confidence intervals or
extrapolation to cryptographic sizes. No generic-prime, double-large-prime,
general-k summation-polynomial or Diem-asymptotic claim is made.

## Rerun without overwriting evidence

Use Python 3.12 and `pycryptosat==5.14.7`. The existing vendored curve, SAT and
pullback adapters are loaded from this repository; there is no new dependency.
Run from the repository root, choosing **new** output paths each time:

```bash
python research/nagao_relations/pullback_incremental/run.py --targets research/nagao_relations/pullback_incremental/targets.json --output research/nagao_relations/parameter_sweep/regression_next
python research/nagao_relations/pullback_incremental/compare.py --input research/nagao_relations/parameter_sweep/regression_next --output research/nagao_relations/parameter_sweep/regression_next_comparison.json
python research/nagao_relations/parameter_sweep/run.py --targets research/nagao_relations/parameter_sweep/targets.json --validate-only --output research/nagao_relations/parameter_sweep/validation_next
python research/nagao_relations/parameter_sweep/run.py --targets research/nagao_relations/parameter_sweep/targets.json --output research/nagao_relations/parameter_sweep/results_next
python research/nagao_relations/parameter_sweep/compare.py --input research/nagao_relations/parameter_sweep/results_next --targets research/nagao_relations/parameter_sweep/targets.json --regression research/nagao_relations/parameter_sweep/regression_next_comparison.json --output research/nagao_relations/parameter_sweep/comparison_next.json --markdown research/nagao_relations/parameter_sweep/RESULTS_next.md
python research/nagao_relations/parameter_sweep/test_core.py
```

`--freeze new_targets.json` reproduces the target construction for the unchanged
contract and historical corpus; normal benchmark reruns should use committed
`targets.json`. A fresh future holdout requires a separately frozen contract,
not editing these inputs. Run timed campaigns serially on an otherwise idle
machine. All benchmark destinations refuse an existing directory. The initial
post-processing relative-path error is retained in `analysis_failure_01.txt`;
no measured source or run was changed to fix it.

## Next bounded experiment

Use the binary9 ell=4,k=3 case to test an adaptive **build-on-reuse** threshold
for pair tables, against both unconditional table construction and enumerate-last.
Freeze thresholds on development cases, retain new holdouts, and rerun the full
old suite plus this sweep. Before promoting the candidate as an end-to-end
speedup, calibrate field/scalar/control costs in one unit and repeat paired
timings. General-k/GF(8) polynomial adapters are a separate implementation task,
not something this generic search experiment establishes by proxy.

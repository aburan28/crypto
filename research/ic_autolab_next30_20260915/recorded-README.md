# AutoLab: another 30% from the verified winner

**Goal achieved: 47.72% fewer complete cold-solve instructions**, measured against the previous verified fast-lifting batch1/window8 winner. AutoLab confirmation and replay passed; the independent audit and canonical scoreboard are complete. The campaign finished and its container stopped after 23.5 minutes.

Confirmation cost ratio: **0.522831359**; paired 95% interval **[0.514432015, 0.532117819]**. This is a 1.913× instruction-cost speedup. Fresh-process replay ratio: **0.522829527**. Both beat the frozen 0.70 target; every cell improved.

The locked candidate was tested on 60 distinct fresh targets across five small Koblitz curve cells (degrees 13, 17, 19 and 23), with three repetitions, both IC arms and matched rho. Confirmation excludes all 289 earlier public targets and this round's earlier-stage targets. All **1,356 profile/native pairs** passed the independent group, relation, column-log, matrix-rank and final-scalar checks. The reporting audit also matched factor bases in all **912 IC profiles**.

## Confirmed comparison

| Variant | S (Ir / sqrt(r)) | Cost / previous winner | Cost / rho | Cost / floor | Verified pairs | Class |
|---|---:|---:|---:|---:|---:|---|
| Previous verified winner | 379671 | 1 | 9.52663 | 1.53287e+07 | 180/180 | reference |
| Fast scalar multiplication | 198504 | 0.522831 | 4.98082 | 8.01435e+06 | 180/180 | engineering experiment |
| Matched rho | 39853.6 | 0.104969 | 1 | unmeasured | 180/180 | reference |

Unit: **Valgrind 3.22 amd64 Ir**, covering the complete cold worker process, including setup, unsuccessful attempts, checks, linear algebra, final recovery and cleanup. Fixed support, m=3, trial cap and all accounting rules are unchanged. The fixed signed-base coverage boundary remains `binomial(B+m-1,m)`; the weak K-instruction floor applies only to this full-rank K-column collector.

This is implementation engineering. Rho still costs less: the candidate uses 4.98× its instructions on these matched targets. Curve diversity is limited; native timings are diagnostic. No arithmetic exponent, broad family result or asymptotic crossover is claimed. Prior rounds used different holdouts, so the successive percentage gains are not combined into a new measured cumulative result.

## Source change and validation

The sole runtime change routes `KoblitzCurve.mul` through the existing `FastCurve::mul` when the point exactly round-trips through its single-word representation. General arithmetic remains the fallback. Scalars retain arbitrary precision. Reduction-table setup stays inside the measured worker; every required check and certificate remains enabled.

Six targeted source tests passed, including comparison with unchanged general scalar arithmetic, random-point arithmetic, corrupted-log rejection, full log recovery, sparse/dense agreement and an even-degree subfield solve. Sixteen existing harness admission/pricing tests also passed. The tested patch and full candidate source are preserved; the library source is not automatically merged by this harness.

The initial eight-target development probe measured 47.17% lower cost with 72/72 pairs independently verified. Its baseline, candidate and rho rows remain in [development measurements](development-measurements.json) and the canonical scoreboard. One development candidate and one final campaign were needed.

## Research cost

Measured containers used **0.4662 CPU-hours** in this goal: 0.0631 for development and 0.4031 for the final campaign. Final peak RAM was 1.69 GiB; no OOM events occurred. No model API calls or new cloud resources were used. Local compute, preparation and this assistant session remain unpriced. See the [full cost assessment](COSTS.md) and [machine-readable ledger](cost-assessment.json).

## Reproduction and artifacts

- [Source patch](candidate-a.patch), [source tests](validation/source-tests.json), [harness tests](validation/harness-tests.json).
- [Audited tournament report](runs/autolab-20260915T224645Z/jobs/autolab-20260915T224645Z/task__4brp5nH/verifier/round-20260915t224645z/REPORT.md), [frozen measurement table](runs/autolab-20260915T224645Z/jobs/autolab-20260915T224645Z/task__4brp5nH/verifier/round-20260915t224645z/measurements.json), [decision](runs/autolab-20260915T224645Z/jobs/autolab-20260915T224645Z/task__4brp5nH/verifier/round-20260915t224645z/decision.json), [independent audit](runs/autolab-20260915T224645Z/jobs/autolab-20260915T224645Z/task__4brp5nH/verifier/round-20260915t224645z/audit.json).
- [Frozen contract](runs/autolab-20260915T224645Z/jobs/autolab-20260915T224645Z/task__4brp5nH/verifier/round-20260915t224645z/contract.json), [candidate source](runs/autolab-20260915T224645Z/jobs/autolab-20260915T224645Z/task__4brp5nH/verifier/round-20260915t224645z/source_candidates/autolab/source/), [raw receipts and profiles](runs/autolab-20260915T224645Z/jobs/autolab-20260915T224645Z/task__4brp5nH/verifier/round-20260915t224645z/runs/).
- [Preregistered proposal](proposal.json), [stricter promotion threshold](evaluator-threshold.patch), [fresh-target exclusions](previous-targets.json), [freshness adapter](fresh_targets.py).
- [Development evidence](development/), [development audit](development-audit.json), [complete confirmation phase costs](confirmation-phase-costs.json).
- [Canonical scoreboard](../../docs/index-calculus-scoreboard.html#ic-tournament-round-20260915t224645z), [completed supervisor](runs/autolab-20260915T224645Z/status.json).

AutoLab upstream commit: `4127da3dde8449be61a1cf9859473b9fbbd51751`; Harbor 0.3.0. Frozen verification used its nop agent with no code edits or model calls after candidate lock. Image: `sha256:68e32e2cf212ceda22065b836dd130e8c89c9172e17684f39da1046122352cb9`.

The frozen evaluator can re-audit the preserved round without rerunning measurements:

```bash
ROUND="/home/ubuntu/crypto/research/ic_autolab_next30_20260915/runs/autolab-20260915T224645Z/jobs/autolab-20260915T224645Z/task__4brp5nH/verifier/round-20260915t224645z"
python3 "$ROUND/evaluator/tournament.py" verify --round "$ROUND"
```

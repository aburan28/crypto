# Source improvement: factor-base point lifting

The same original batch16 source/configuration is the baseline throughout. The source-only change reduced complete instruction cost by 15.6%; the separately measured source plus batch1/window8 candidate reduced it by 28.5% on development inputs. The final campaign now confirms **26.9% lower complete instruction cost** (1.367×), with a 95% reduction interval of 24.0–29.4%; fresh-process replay agrees. The 20% goal is achieved on the five tested Koblitz cells.

## Final confirmation and replay

[Final research report](../runs/autolab-20260915T205319Z/jobs/autolab-20260915T205319Z/task__s4YLvUd/verifier/round-20260915t205319z/REPORT.md) · [Frozen measurements](../runs/autolab-20260915T205319Z/jobs/autolab-20260915T205319Z/task__s4YLvUd/verifier/round-20260915t205319z/measurements.json) · [Independent audit](../runs/autolab-20260915T205319Z/jobs/autolab-20260915T205319Z/task__s4YLvUd/verifier/round-20260915t205319z/audit.json).

All 1,356 profile/native pairs verified, including 60 fresh confirmation targets in five curve cells, three repetitions for each of the two IC arms and rho, and a fresh-process replay. All 912 complete IC profiles used matching factor-base support. Neither fresh target seeds nor target points overlapped the first campaign. The source was locked before these fixtures were created.

| Variant | S (Ir / sqrt(r)) | Cost / original baseline | Cost / rho | Cost / floor | Verified confirmation pairs | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 515445 | 1 | 13.0517 | 2.08105e+07 | 180/180 | reference |
| autolab | 376934 | 0.731278 | 9.54443 | 1.52182e+07 | 180/180 | engineering experiment |
| rho | 39492.5 | 0.0766183 | 1 | unmeasured | 180/180 | reference |

The confirmed ratio is 0.731278 (95% interval 0.706305–0.760035); replay is 0.731279. Every cell improved, with the least favorable cell at 0.781706. Rho still uses only 0.104773 times the candidate instructions: this is an engineering improvement, with no rho crossover or scaling claim. Native times remain diagnostics.

The [confirmation phase ledger](confirmation-phase-costs.json) records summed Ir over 180 profiles per arm. Factor-base/table setup fell 30.9%; verification/filtering/linear algebra fell 34.3%. These phase sums are descriptive; the headline ratio uses the frozen curve-weighted paired estimator.

The full gain requires both the [source patch](fastlift.patch) and [batch1/window8 configuration](candidate-config.json). The source is implemented and tested in the isolated candidate and frozen final snapshot; library defaults are not automatically changed. [Cost assessment](../COSTS.md).

## Mechanism

The private fast path uses the existing single-word field arithmetic for odd-degree factor-base lifts. It preserves the exact half-trace root and return order, checks the equivalent root equation, reuses field reduction tables, and keeps the general path for even/wide fields. It changes neither the public target constructor nor the rho implementation. No work is moved outside the measured process.

## Development evidence

All 144 profile/native pairs were independently reverified from raw artifacts. Four curve cells, two targets per cell and three repetitions per arm are exploratory evidence; they are not 24 independent curves.

| Experiment | Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified pairs | Class |
|---|---|---:|---:|---:|---:|---:|---|
| Source only | incumbent | 568508 | 1 | 13.0047 | 2.04519e+07 | 24/24 | reference |
| Source only | candidate | 479624 | 0.843653 | 10.9715 | 1.72543e+07 | 24/24 | engineering experiment |
| Source only | rho | 43715.6 | 0.0768952 | 1 | unmeasured | 24/24 | reference |
| Source + tuning | incumbent | 568507 | 1 | 13.0047 | 2.04519e+07 | 24/24 | reference |
| Source + tuning | candidate | 406674 | 0.715337 | 9.30272 | 1.463e+07 | 24/24 | engineering experiment |
| Source + tuning | rho | 43715.6 | 0.0768954 | 1 | unmeasured | 24/24 | reference |

The instruction unit covers complete worker startup through termination. Kernel/device, profiling implementation, builds and the independent audit are research overhead. The weak floor is K instructions for this full-rank K-column collector. Fixed signed-base coverage remains at most binomial(B+m-1,m); no coverage or exponent improvement is claimed. Rho remains substantially cheaper on these small fixtures.

## Correctness and cost

Five actual tests passed: ordered-root equivalence, factor-base invariance, pair-table/exhaustive-decomposition agreement, full logarithm recovery and even-degree subfield recovery. A mistyped test module initially selected zero tests; that invocation is retained and the corrected exact test passed. One development preparation failed because copied files were unreadable; it made no measurements and is retained.

Source development consumed 470.35 container CPU-seconds and peaked at 1.49 GiB RAM. It made no model API calls. Host test compilation and image preparation are separate overhead; no local CPU-dollar rate was supplied.

## Reproducible artifacts

- [Source patch](fastlift.patch) (applies cleanly to the current repository; not automatically merged).
- [Preregistered hypothesis](proposal.json), [development measurements](development-table.json), [independent development audit](development-audit.json).
- [Exact test results](source-tests.json), [raw development evidence](../../ic_autolab_evidence_20260915/README.md), [resource accounting](development-resources.json).
- [Completed final run](../runs/autolab-20260915T205319Z/status.json), [full report](../runs/autolab-20260915T205319Z/jobs/autolab-20260915T205319Z/task__s4YLvUd/verifier/round-20260915t205319z/REPORT.md), [matching configuration](candidate-config.json).

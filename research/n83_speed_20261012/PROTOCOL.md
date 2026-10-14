# Exact scalar arithmetic and expanded rho policy, frozen before measurement

This follow-up answers the request to make the degree-83 solver faster while
retaining a fair, optimized control. It is an engineering comparison on public
synthetic one-target workloads, not an exponent or deployed-curve claim.
The earlier report and measurements remain historical evidence.

## Identity and control

The curve is `icv1-f2m83-t6151469093347-9eba6cd5`: a=1, b=1,
polynomial z^83+z^7+z^4+z^2+1, subgroup order 8569786107849059,
cofactor 1128547018. The generator is
(6614244651398184860292015,8868260410029249391539529).
The signed Frobenius quotient has 166 phases and lambda=2998004823950457.
This legacy representation is distinct from the primary a=0 confidence gate.

The matched control is the previous selected Fold field backend, 64 lockstep
lanes, DP6, disabled fruitless escape, and reference u64 scalar arithmetic.
Both arms are monomorphizations in the same newly frozen executable. The
reference arithmetic is unchanged; candidate arithmetic uses exact Barrett
reduction and conditional modular addition. Both use the same random sampler.
Neither arm receives validation scalars. Larger lane/DP settings change walk
policy and may change trajectories; their total cost is measured accordingly.

## Fresh targets and fixed stages

Native C++ generates three training and five holdout points, with generator
seed 202610126201 and conditionally uniform nonzero subgroup scalar proposals.
It rejects the signed Frobenius orbits of the original target, all K600 points,
all eight prior optimization targets, and points accepted in this campaign.
Point-only files are solver inputs; a separate validation sidecar is checker-only.
The generator and independent verifier retain every proposal and rejection.

1. A/A noise: five pairs of identical reference runs on training_1 with seed
   202610126211, ten processes. Retain all repetitions, median, minimum and
   pair spread; these are excluded from selection and held-out estimates.
2. Arithmetic: reference versus Barrett, Fold64/DP6/off, three training targets
   times three seeds (202610126211 through 202610126219), 18 processes.
   Every pair must agree on answer and recorded walk counters and charges.
   This checks recorded semantics, not a complete step trace.
3. Policy: selected scalar kernel times lanes {64,128,256} and DP {6,8,10},
   escape fixed off, on the same nine training cases, 81 processes.
4. Holdout: reference Fold64/DP6/off versus the selected kernel and policy,
   five fresh points times three seeds (202610126241 through 202610126255),
   30 processes, alternating order within successive pairs.

Total: 139 rho processes. All matrices, seeds, orders, source hashes, binaries,
fixtures and this protocol are frozen before the first timed child. Arithmetic
and policy selections are immutable before holdout. No holdout tuning or
post hoc target/seed exclusion is permitted. No new IC run or new IC/rho ratio
is part of this follow-up; the earlier matched cross-method study is preserved.

## Selection and accounting

Select the arithmetic mean actual whole-process wall cost among configurations
that complete and verify all nine training cases. Prefer the reference kernel
or the selected-kernel 64/DP6 default if within 2% of the minimum. The arithmetic
stage also requires all nine recorded-semantic matches; an ineligible candidate
cannot win. Preserve all losing policies, failures, timeouts and rejected cases.
Charge failed/capped cases their 60-second cap in descriptive comparisons;
an incomplete arm cannot earn an unqualified win.

Whole-process wall includes startup, point parsing, normal-basis conversion,
table construction, walk, collision recovery, built-in scalar replay and output.
Independent C++ replay is outside both process intervals and separately reported.
Online times and deterministic counters are secondary diagnostics. A measured
whole-process gain is an engineering runtime gain on this host, not a reduction
in the generic counting boundary or a calibrated operation-equivalent S claim.
The generic reference scale sqrt(r/166) remains unchanged. No common calibrated
primitive-cost unit, L2 admission, or independent-host replay is asserted.

## Isolation and statistical estimand

One OS thread per child, CPU3, RLIMIT_AS 12 GiB, timeout 60 seconds, identical
clean environment and native 5 ms wait4 monitor. The native isolated-bench tool
reserves the core and observes contention for each stage. Builds, tests, plotting
and audits pause during timing. Retain per-process wall, CPU, RSS, exit status,
environment, commands and file hashes. VM host neighbors remain unobservable.
Contended reservations do not support wall claims and are never silently pooled.

The primary held-out estimand is the geometric mean of reference/selected
whole-process cost ratios after taking each target's three-seed arithmetic mean.
Native target-block bootstrap: 20,000 resamples of the five target blocks, seed
202610126299, percentile 95% interval. Seeds are nested in targets, not 15
independent target samples. Report each target and every seed, arithmetic mean
costs and ratio, online/CPU/RSS/counters, and A/A spread. Five targets limit
population generalization. An interval including one or a gain inside A/A noise
does not establish a faster runtime. All answers require independent native
scalar replay. Equality of scalar-kernel counters does not establish every
intermediate state or a full rank transcript.

## Deliverables and gates

Retain implementation, written arithmetic proof and formalization status,
meaningful exactness tests, native controller tests, full release library checks,
fixture and all-answer replay, immutable raw records and checksums. Produce a
source-linked visual report and PDF, prepend it to the existing comprehensive
PDF without altering its old pages, and copy visual records into benchmarks.
Update only dashboard fields this evidence can support; do not promote a ledger
row without missing admission and operation-accounting gates. Commit, push and
open a ready follow-up PR. Conductor was attempted but is unavailable here.

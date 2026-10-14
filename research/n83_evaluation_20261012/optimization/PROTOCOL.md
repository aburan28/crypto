# N83 fair-rho optimization follow-up — frozen before timing

User requirement: optimize the rho control too, so the IC comparison does not
select a weak reference. This follow-up tests a bounded, declared optimization
search and evaluates its single selected configuration on disjoint targets.
It does not claim the global optimum. Original five-pair records and reports
remain immutable historical evidence against the frozen default strong rho.

## Exact workload and implementation

Use the legacy degree-83 a=1 curve, b=1, modulus z^83+z^7+z^4+z^2+1,
prime subgroup r=8569786107849059, cofactor1128547018,
G=(6614244651398184860292015,8868260410029249391539529),
lambda2998004823950457, signed Frobenius orbit size166.
Identity: icv1-f2m83-t6151469093347-9eba6cd5. No substitution of the
registered field with tail [0,1,2,45] or the primary a=0 degree-83 curve.

Both rho backends execute the existing generic strong algorithm: signed
Frobenius quotient, normal-basis canonicalization, distinguished points,
G-only jumps and lockstep batch inversion. The default backend is WideGf2;
the candidate is a thin RhoField adapter to the existing Gf2_128 primitives
used by IC. It offers the existing carryless fold reduction and specialized
squaring. Tests must check field results and identical trajectories/charges
before any performance interpretation. No floating scalar approximation;
this subgroup exceeds the safe threshold. Both use the same standard release
compiler profile with runtime hardware kernel detection, without asymmetric
target-cpu, LTO, PGO or compiler flags.

The new point-only worker receives a public Q, walk seed and declared policy.
It never opens the known-answer sidecar. The parent and standalone C++ checker
may use the separate sidecar for validation. Source and executable identities,
all arguments, inputs and configuration hashes are recorded.

## Targets, seeds and candidate budget

Standalone native arithmetic generates three training and five held-out
public targets from uniform nonzero scalar proposals using unbiased rejection sampling
and SplitMix64 seed202610128301. It rejects shared signed-Frobenius orbits
with the original fixed target, retained K600 base or another generated target.
Thus the target law is conditional on those exclusions, without replacement;
the targets are not claimed to be unconditional independent draws.
The scalar sidecar is validation-only. Freeze files/hashes before timed runs;
hold-out costs cannot influence either selection stage.

1. Kernel stage: Wide and Fold, each with lanes32, DP4, escape disabled,
   step_cap_factor2000; three training targets by three seeds, 18 launches.
   Training walk seeds202610128311 through202610128319, three per target.
   Alternate kernel order between cases. Verify equal walk outcomes and charges
   across backends on every same-input case. Freeze one backend selection.
2. Policy stage: selected backend, lanes16/32/64 by DP4/6 by escape disabled/
   enabled, twelve candidates and the same nine training cases, 108 launches.
   Rotate candidate order between cases. Every candidate receives all nine
   cases; retain caps and failures. Freeze one configuration before holdout.
3. Held-out rho: original Wide lanes32/DP4/escape-disabled and the one selected
   configuration on five held-out targets by three fresh seeds per target,
   30 launches. Seeds202610128341 through202610128355. Alternate arm order.
4. Prepared IC: one fresh process of the immutable original IC binary on the
   same five held-out targets, retained K600 base, rank seed29, one thread.
   Its full preparation is measured and reported separately. Record each
   target's complete online interval, five phases, residual and verification.
   This batch does not replace five cold single-target IC measurements.

No fastest-seed selection, per-target policy selection, adaptive pruning or
post-holdout tuning. There are 126 training rho launches and30 hold-out rho
launches, plus one IC process. The finite maximum rho charge is9360s at60s
per launch. This is a declared optimization budget, not a runtime prediction.

## Selection rule and accounting

Kernel/policy eligibility requires all nine verified, uncapped completions;
any wrong answer rejects the candidate. Among eligible candidates choose the
lowest arithmetic mean whole-process wall cost, with the default candidate
preferred when within2% of the minimum. Save all costs and decisions. If no
candidate completes all nine, freeze that failure and block dependent stages.
Online mean is secondary; a cold-selected policy is not an online optimum.
Never use ns/step alone to select policy: parameters change collision progress.

Whole-process rho wall/CPU includes process launch, input load, curve/normal
basis, G-only jump preparation, Q checks, walk, collision, scalar replay and
receipt output. Record wait4 monitoring resolution and its possible exit-
detection/scheduling overhead. Producer online starts before Q validation and
ends after independent generic-library scalar replay; report excluded setup
separately. Dataset construction and training are explicit research costs,
not charged anew to a selected production solve. Memory and failures remain
part of the result. Cold scope still excludes original retained-base creation
and OS cache flushing, as in the original report.

Common rho envelope: CPU3, one OS thread, RLIMIT_AS12GiB,60s full-process
cap. IC envelope: same CPU/thread/memory,1800s process cap. Stock extraction
scans the finite index once per target and has no20-million-probe cap; the
third CLI argument is a rank seed. These caps are explicitly unequal; there is
no shared calibrated unit equating probes to rho steps. Cross-method cost
ratios require all five uncensored verified IC targets and all held-out rho
completions. Otherwise preserve failures and report conditional completed-
subset descriptions only, without an unqualified full-workload winner.

For held-out rho, average the three walk costs within each target; report
geometric mean default/selected ratios and20,000 percentile bootstrap draws
over the five target blocks (seed202610128399). Fresh prepared IC online
comparison uses each target's IC interval versus that target's mean of three
selected-rho online intervals, using bootstrap seed202610128400. Its target-block interval reflects five target
cases, with only one IC preparation and one IC measurement per target; it
does not establish preparation/host replication. Keep exact direction in
every ratio label. Report arithmetic means and individual cases too.

## Isolation, correctness and report gates

Native isolated_bench reserves/monitors the sequential campaign; all workers
pin to CPU3. Start no builds, tests or other heavy jobs during timed work.
Preserve reservation records and contention classifications. Independent C++
arithmetic must replay all completed scalars and all IC point decompositions.
Rank600/count5/basehash are semantic requirements. The full factor-log/rank
transcript remains unavailable, as in the earlier evaluation.

This exact legacy case lacks the applicable sealed ecbench wide IC adapter.
Local timings remain diagnostic: no admitted ledger promotion, common S,
global optimum, deployed-curve claim or independent-host claim. Report this
optimization evidence and the unchanged historical comparison separately;
refresh the comprehensive PDF with a leading source-linked addendum.

Conductor was rechecked for the expanded new-example and study-directory
scope and returned command-not-found. Metadata inspection found that the
original PR1719 had already merged. This authorized follow-up uses the new
feature branch codex/n83-rho-optimization-20261012 and a separate result PR;
old frozen evidence is not edited. These delivery metadata changes precede
the first timed run and do not alter the experiment or selection rules.
